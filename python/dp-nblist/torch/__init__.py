import torch
import torch.nn as nn

class TorchBaseNBL(nn.Module):

    def forward(
        self,
        positions
    ):
        return self.update(positions)
    
    def update(self, positions):
        pass

    def build(
        self,
        Z: torch.Tensor,
        positions: torch.Tensor,
        box: torch.Tensor,
        pbc: torch.Tensor,
        cutoff: float,
    ):
        """Override with specific neighbor list implementation"""
        raise NotImplementedError

class TorchNeighborList(TorchBaseNBL):
    """
    Environment provider making use of neighbor lists as implemented in TorchAni

    Supports cutoffs and PBCs and can be performed on either CPU or GPU.

    References:
        https://github.com/aiqm/torchani/blob/master/torchani/aev.py
    """
    def build(self, positions, box, cutoff, pbc):
        # Check if shifts are needed for periodic boundary conditions
        if torch.all(pbc == 0):
            shifts = torch.zeros(0, 3, device=box.device, dtype=torch.long)
        else:
            shifts = self._get_shifts(box, pbc, cutoff)
        idx_i, idx_j, offset = self._get_neighbor_pairs(positions, box, shifts, cutoff)

        # Create bidirectional id arrays, similar to what the ASE neighbor_list returns
        bi_idx_i = torch.cat((idx_i, idx_j), dim=0)
        bi_idx_j = torch.cat((idx_j, idx_i), dim=0)

        # Sort along first dimension (necessary for atom-wise pooling)
        sorted_idx = torch.argsort(bi_idx_i)
        idx_i = bi_idx_i[sorted_idx]
        idx_j = bi_idx_j[sorted_idx]

        bi_offset = torch.cat((-offset, offset), dim=0)
        offset = bi_offset[sorted_idx]
        offset = torch.mm(offset.to(box.dtype), box)

        return idx_i, idx_j, offset

    def _get_neighbor_pairs(self, positions, box, shifts, cutoff):
        """Compute pairs of atoms that are neighbors
        Copyright 2018- Xiang Gao and other ANI developers
        (https://github.com/aiqm/torchani/blob/master/torchani/aev.py)
        Arguments:
            positions (:class:`torch.Tensor`): tensor of shape
                (molecules, atoms, 3) for atom coordinates.
            box (:class:`torch.Tensor`): tensor of shape (3, 3) of the three vectors
                defining unit box: tensor([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]])
            shifts (:class:`torch.Tensor`): tensor of shape (?, 3) storing shifts
        """
        num_atoms = positions.shape[0]
        all_atoms = torch.arange(num_atoms, device=box.device)

        # 1 Central box
        pi_center, pj_center = torch.combinations(all_atoms).unbind(-1)
        shifts_center = shifts.new_zeros(pi_center.shape[0], 3)

        # 2 boxs with shifts
        # shape convention (shift index, molecule index, atom index, 3)
        num_shifts = shifts.shape[0]
        all_shifts = torch.arange(num_shifts, device=box.device)
        shift_index, pi, pj = torch.cartesian_prod(
            all_shifts, all_atoms, all_atoms
        ).unbind(-1)
        shifts_outside = shifts.index_select(0, shift_index)

        # 3 combine results for all boxs
        shifts_all = torch.cat([shifts_center, shifts_outside])
        pi_all = torch.cat([pi_center, pi])
        pj_all = torch.cat([pj_center, pj])

        # 4 Compute shifts and distance vectors
        shift_values = torch.mm(shifts_all.to(box.dtype), box)
        Rij_all = positions[pi_all] - positions[pj_all] + shift_values

        # 5 Compute distances, and find all pairs within cutoff
        distances = torch.norm(Rij_all, dim=1)
        in_cutoff = torch.nonzero(distances < cutoff, as_tuple=False)

        # 6 Reduce tensors to relevant components
        pair_index = in_cutoff.squeeze()
        atom_index_i = pi_all[pair_index]
        atom_index_j = pj_all[pair_index]
        offsets = shifts_all[pair_index]

        return atom_index_i, atom_index_j, offsets

    def _get_shifts(self, box, pbc, cutoff):
        """Compute the shifts of unit box along the given box vectors to make it
        large enough to contain all pairs of neighbor atoms with PBC under
        consideration.
        Copyright 2018- Xiang Gao and other ANI developers
        (https://github.com/aiqm/torchani/blob/master/torchani/aev.py)
        Arguments:
            box (:class:`torch.Tensor`): tensor of shape (3, 3) of the three
            vectors defining unit box: tensor([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]])
            pbc (:class:`torch.Tensor`): boolean vector of size 3 storing
                if pbc is enabled for that direction.
        Returns:
            :class:`torch.Tensor`: long tensor of shifts. the center box and
                symmetric boxs are not included.
        """
        reciprocal_box = box.inverse().t()
        inverse_lengths = torch.norm(reciprocal_box, dim=1)

        num_repeats = torch.ceil(cutoff * inverse_lengths).long()
        num_repeats = torch.where(
            pbc, num_repeats, torch.Tensor([0], device=box.device).long()
        )

        r1 = torch.arange(1, num_repeats[0] + 1, device=box.device)
        r2 = torch.arange(1, num_repeats[1] + 1, device=box.device)
        r3 = torch.arange(1, num_repeats[2] + 1, device=box.device)
        o = torch.zeros(1, dtype=torch.long, device=box.device)

        return torch.cat(
            [
                torch.cartesian_prod(r1, r2, r3),
                torch.cartesian_prod(r1, r2, o),
                torch.cartesian_prod(r1, r2, -r3),
                torch.cartesian_prod(r1, o, r3),
                torch.cartesian_prod(r1, o, o),
                torch.cartesian_prod(r1, o, -r3),
                torch.cartesian_prod(r1, -r2, r3),
                torch.cartesian_prod(r1, -r2, o),
                torch.cartesian_prod(r1, -r2, -r3),
                torch.cartesian_prod(o, r2, r3),
                torch.cartesian_prod(o, r2, o),
                torch.cartesian_prod(o, r2, -r3),
                torch.cartesian_prod(o, o, r3),
            ]
        )