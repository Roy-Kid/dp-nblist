#include "refCellList.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace dpnblist
{

CellList::CellList(Box *box, double r_cutoff) : _box(box), _r_cutoff(r_cutoff)
{
    auto cell_length = box->get_lengths();
    _cell_length = Vec3<int>(cell_length / r_cutoff);
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
        for (int j = 0; j < 3; j++){
            if (xyz_cell_index[j] == _cell_length[j])
                xyz_cell_index[j] = xyz_cell_index[j] - 1;
        }
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
    // the cell is described by the vector, the cell_vector is the vector of the central cell, the matrix is the offset vector, the (matrix + cell_vector) is the 26 neighbors cell vector.
    // the vector can be express as r = f * A, A is the box matrix, f is the fractional coordinate, can be express as f = inv_A * r
    // if periodic boundary condition is considered, wrapped_f = f - floor(f), wrapped_f is the wrapped fractional coordinate, so t = A * wrapped_f, t is the wrapped vector

    std::vector<size_t> neighbors;
    Vec3<int> cell_vector = get_cell_vector(cell_index);
    // define a offset matrix 26*3
    std::vector<Vec3<int>> matrix = {
        {-1, -1, -1},{-1, -1, 0},{-1, -1, 1},{-1, 0, -1},{-1, 0, 0},{-1, 0, 1},{-1, 1, -1},{-1, 1, 0},{-1, 1, 1},
        {0, -1, -1},{0, -1, 0},{0, -1, 1},{0, 0, -1},{0, 0, 1},{0, 1, -1},{0, 1, 0},{0, 1, 1},
        {1, -1, -1},{1, -1, 0},{1, -1, 1},{1, 0, -1},{1, 0, 0},{1, 0, 1},{1, 1, -1},{1, 1, 0},{1, 1, 1} 
    };
    Mat3<double> A = {_cell_length[0], 0, 0, 0, _cell_length[1], 0, 0, 0, _cell_length[2]};
    Mat3<double> inv_A = A.invert();
    Vec3<double> f;
    Vec3<double> wrapped_f;
    Vec3<int> warp_neb_vector;
    for (int i = 0; i < 26; i++){
        Vec3<int> neb_vector = cell_vector + matrix[i];
        f = inv_A * neb_vector;
        wrapped_f = f - floor(f);
        warp_neb_vector = A * wrapped_f;
        neighbors.push_back(get_cell_index(warp_neb_vector));
    }
    neighbors.push_back(cell_index);
    return neighbors;
}

NeighborList::NeighborList(Box *box, double r_cutoff)
    : BaseNBL(box), _r_cutoff(r_cutoff), _cell_list(box, r_cutoff)
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

namespace py = pybind11;
PYBIND11_MODULE(dpnblist, m) {
    py::class_<dpnblist::Vec3<double>>(m, "Vec3")
        .def(py::init<const double &, const double &, const double &>());
    
    py::class_<dpnblist::Box>(m, "Box")
        .def(py::init<dpnblist::Vec3<double>, dpnblist::Vec3<double>>(), py::arg("lengths"), py::arg("angles") = dpnblist::Vec3<double>(90, 90, 90));
    
    py::class_<dpnblist::NeighborList>(m, "NeighborList")
        .def(py::init<dpnblist::Box*, double>())
        .def("build", &dpnblist::NeighborList::build)
        .def("update", &dpnblist::NeighborList::update)
        .def("get_listArray", &dpnblist::NeighborList::get_listArray);
}