#include "refCellList.h"

namespace dpnblist
{

CellList::CellList(Box *box, double r_cutoff) : _box(box), _r_cutoff(r_cutoff)
{
    auto cell_length = box->get_lengths();
    _cell_length = cell_length / r_cutoff;
}

size_t CellList::get_cell_index(Vec3<int> &cell_vector) const
{
    // c = cxLcyLcz + cyLcz + cz
    return cell_vector[0] * _cell_length[1] * _cell_length[2] + cell_vector[1] * _cell_length[2] + cell_vector[2];
}

Vec3<int> CellList::get_cell_vector(size_t cell_index) const
{
    // cx = c / (LcyLcz)
    // cy = (c / Lcz) - cxLcy or (c / Lcz) mod Lcy
    // cz = c mod Lcz
    Vec3<int> cell_vector;
    cell_vector[0] = cell_index / (_cell_length[1] * _cell_length[2]);
    cell_vector[1] = std::fmod((cell_index / _cell_length[2]), _cell_length[1]);
    cell_vector[2] = std::fmod(cell_index, _cell_length[2]);
    return cell_vector;
}

void CellList::build(std::vector<Vec3<double>> &xyz)
{

    reset();

    update(xyz);
}

void CellList::update(std::vector<Vec3<double>> &xyz)
{
    size_t n_atoms = xyz.size();
    _natoms += n_atoms;
    _lscl.resize(_natoms, EMPTY);
    Vec3<int> xyz_cell_index;
    xyz = _box->wrap(xyz);
    for (size_t i = 0; i < n_atoms; i++)
    {
        xyz_cell_index = Vec3<int>(xyz[i] / _r_cutoff);
        size_t cell_index = get_cell_index(xyz_cell_index);
        _lscl[i] = _head[cell_index];
        _head[cell_index] = i;
    }
}

void CellList::reset()
{
    _natoms = 0;
    _head.resize(_cell_length[0] * _cell_length[1] * _cell_length[2], EMPTY);
    _lscl.resize(_natoms, EMPTY);
}

size_t CellList::get_ncells()
{
    return _cell_length[0] * _cell_length[1] * _cell_length[2];
}

std::vector<size_t> CellList::get_atoms_in_cell(size_t cell_index) const
{
    std::vector<size_t> atoms_in_cell;
    size_t atom_index = _head[cell_index];
    while (atom_index != EMPTY)
    {
        atoms_in_cell.push_back(atom_index);
        atom_index = _lscl[atom_index];
    }
    return atoms_in_cell;
}

std::vector<size_t> CellList::get_neighbors(size_t cell_index) const
{
    std::vector<size_t> neighbors;
    Vec3<int> cell_vector = get_cell_vector(cell_index);
    int x1, y1, z1;
    for (int x = cell_vector[0] - 1; x <= cell_vector[0] + 1; x++)
    {
        for (int y = cell_vector[1] - 1; y <= cell_vector[1] + 1; y++)
        {
            for (int z = cell_vector[2] - 1; z <= cell_vector[2] + 1; z++)
            {
                x1 = x;
                y1 = y;
                z1 = z;
                // periodic
                if (x < 0)
                {
                    x1 = x + _cell_length[0];
                }
                else if (x >= _cell_length[0])
                {
                    x1 = x - _cell_length[0];
                }

                if (y < 0)
                {
                    y1 = y + _cell_length[1];
                }
                else if (y >= _cell_length[1])
                {
                    y1 = y - _cell_length[1];
                }

                if (z < 0)
                {
                    z1 = z + _cell_length[2];
                }
                else if (z >= _cell_length[2])
                {
                    z1 = z - _cell_length[2];
                }

                Vec3<int> neighbor_cell_vector = {x1, y1, z1};
                size_t neighbor_cell_index = get_cell_index(neighbor_cell_vector);
                if (neighbor_cell_index != cell_index)
                {
                    neighbors.push_back(neighbor_cell_index);
                }
            }
        }
    }

    return neighbors;
}

NeighborList::NeighborList(Box *box, double r_cutoff)
    : BaseNeighborList(box), _r_cutoff(r_cutoff), _cell_list(box, r_cutoff)
{
}

NeighborList::~NeighborList()
{
}

void NeighborList::reset()
{
    _neighborListArray.clear();
}

void NeighborList::build(std::vector<Vec3<double>> &xyz)
{
    size_t natoms = xyz.size();

    _neighborListArray.resize(natoms);

    _cell_list.build(xyz);

    update(xyz);
}

void NeighborList::update(std::vector<Vec3<double>> &xyz)
{
    size_t n_cells = _cell_list.get_ncells();
    for (size_t cell_index = 0; cell_index < n_cells; ++cell_index)
    {
        for (auto neighbor_cell : _cell_list.get_neighbors(cell_index))
        {
            for (size_t i : _cell_list.get_atoms_in_cell(cell_index))
            {
                for (size_t j : _cell_list.get_atoms_in_cell(neighbor_cell))
                {
                    if (i < j) // Avoid double counting
                    {
                        Vec3<double> pos_i = xyz[i];
                        Vec3<double> pos_j = xyz[j];
                        double r = _box->calc_distance(pos_i, pos_j);
                        if (r < _r_cutoff)
                        {
                            _neighborListArray[i].push_back(j);
                            _neighborListArray[j].push_back(i);
                        }
                    }
                }
            }
        }
    }
}

std::vector<std::vector<size_t>> NeighborList::get_listArray()
{
    return _neighborListArray;
}

} // namespace dpnblist
