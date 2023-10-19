#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <omp.h>

namespace py = pybind11;

class GB_NBL_MP_cube {
public:
    GB_NBL_MP_cube(std::vector<double> cube_size, int num_particles, double cut_off_radius)
        : cube_size(cube_size), rc(cut_off_radius), num_particles(num_particles)
    {
        lc = {static_cast<int>(cube_size[0] / cut_off_radius),
              static_cast<int>(cube_size[1] / cut_off_radius),
              static_cast<int>(cube_size[2] / cut_off_radius)};

        grid_size = {cube_size[0] / lc[0],
                     cube_size[1] / lc[1],
                     cube_size[2] / lc[2]};

        num_cells = lc[0] * lc[1] * lc[2];

        particle_list.resize(num_particles);
        cell_list.resize(num_cells);
    }

    std::vector<int> xyz2ind(const std::vector<std::vector<double>>& xyz) {
        std::vector<int> cell_indices;
        cell_indices.reserve(xyz.size());
        for (const auto& point : xyz) {
            std::vector<int> cell_index(3);

            for (int i = 0; i < 3; i++) {
                double x = fmod(point[i], cube_size[i]);
                if (x < 0) x += cube_size[i];

                cell_index[i] = static_cast<int>(floor(x / grid_size[i]));
            }
            cell_indices.push_back(cell_index[0] * lc[1] * lc[2] + cell_index[1] * lc[2] + cell_index[2]);
        }
        return cell_indices;
    }

    std::vector<std::vector<double>> ind2xyz(const std::vector<int>& ind) {
        std::vector<std::vector<double>> xyz;
        xyz.reserve(ind.size());
        for (const auto& index : ind) {
            std::vector<double> point(3);
            point[0] = (index / (lc[1] * lc[2])) * grid_size[0] + grid_size[0] / 2;
            point[1] = (index / lc[2]) * grid_size[1] + grid_size[1] / 2;
            point[2] = (index % lc[2]) * grid_size[2] + grid_size[2] / 2;
            xyz.push_back(point);
        }
        return xyz;
    }

    std::vector<double> get_min_diff(const std::vector<double>& xyz1, const std::vector<double>& xyz2) {
        std::vector<double> difference(3);
        for (int i = 0; i < 3; i++) {
            double diff = xyz2[i] - xyz1[i];
            diff = std::fmod(diff + cube_size[i] / 2, cube_size[i]);
            if (diff < 0) diff += cube_size[i];
            diff -= cube_size[i] / 2;
            difference[i] = diff;
        }
        return difference;
    }

    std::unordered_set<int> get_neighbor_cells(int cell_ind) {
        std::vector<int> cell_ind_list(1, cell_ind);
        std::vector<std::vector<double>> xyz = ind2xyz(cell_ind_list);

        std::unordered_set<int> adjacent_cells;
        const int num_neighbors = 27; // Number of adjacent cells including the cell itself
        int neighbors[num_neighbors][3];
        int di, dj, dk, k;

        #pragma omp parallel for private(di, dj, dk, k)
        for (int ni = -1; ni <= 1; ni++) {
            for (int nj = -1; nj <= 1; nj++) {
                for (int nk = -1; nk <= 1; nk++) {
                    di = ni;
                    dj = nj;
                    dk = nk;
                    k = (di + 1) * 9 + (dj + 1) * 3 + (dk + 1);
                    neighbors[k][0] = di;
                    neighbors[k][1] = dj;
                    neighbors[k][2] = dk;
                }
            }
        }

        #pragma omp parallel for reduction(unordered_set_merge: adjacent_cells)
        for (int k = 0; k < num_neighbors; k++) {
            double ni = static_cast<int>(xyz[0][0]) + neighbors[k][0] * static_cast<int>(grid_size[0]);
            double nj = static_cast<int>(xyz[0][1]) + neighbors[k][1] * static_cast<int>(grid_size[1]);
            double nk = static_cast<int>(xyz[0][2]) + neighbors[k][2] * static_cast<int>(grid_size[2]);
            std::vector<int> ind_temp = xyz2ind({{ni, nj, nk}});
            adjacent_cells.insert(ind_temp[0]);
        }
        return adjacent_cells;
    }

    std::vector<double> get_neighbors(int particle_seq) {
        return particle_list[particle_seq];
    }

    // void constructor(const std::vector<std::vector<double>>& inputs) {
    //     // 把粒子添加到对应的cell之中
    //     std::vector<int> particle_inds = xyz2ind(inputs);

    //     #pragma omp parallel for
    //     for (int i = 0; i < num_particles; i++) {
    //         cell_list[particle_inds[i]].push_back(i);
    //     }

    //     // 建立链接关系，考虑rc
    //     #pragma omp parallel for
    //     for (int particle_seq = 0; particle_seq < num_particles; particle_seq++) {
    //         std::vector<double> particle_xyz = inputs[particle_seq];

    //         int particle_ind = particle_inds[particle_seq];
    //         std::unordered_set<int> adjacent_cells = get_neighbor_cells(particle_ind);

    //         for (const auto& cell_temp : adjacent_cells) {
    //             const std::vector<int>& neighbor_particles_temp = cell_list[cell_temp];
    //             if (!neighbor_particles_temp.empty()) {
    //                 for (const auto& neighbor_particle : neighbor_particles_temp) {
    //                     if (neighbor_particle == particle_seq) continue;
    //                     const std::vector<double>& neighbor_xyz = inputs[neighbor_particle];
    //                     std::vector<double> diff = get_min_diff(particle_xyz, neighbor_xyz);

    //                     double distance = 0;
    //                     for (int i = 0; i < 3; i++) {
    //                         distance += diff[i] * diff[i];
    //                     }
    //                     distance = sqrt(distance);

    //                     if (distance - rc < 1e-10) {
    //                         particle_list[particle_seq].push_back(neighbor_particle);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    void constructor(const std::vector<std::vector<double>>& inputs) {
    // 把粒子添加到对应的cell之中
    std::vector<int> particle_inds = xyz2ind(inputs);

    #pragma omp parallel for
    for (int i = 0; i < num_particles; i++) {
        cell_list[particle_inds[i]].push_back(i);
    }

    // 建立链接关系，考虑rc
    #pragma omp parallel for
    for (int particle_seq = 0; particle_seq < num_particles; particle_seq++) {
        std::vector<double> particle_xyz = inputs[particle_seq];
        int particle_ind = particle_inds[particle_seq];
        std::unordered_set<int> adjacent_cells = get_neighbor_cells(particle_ind);

        #pragma omp parallel for
        for (const auto& cell_temp : adjacent_cells) {
            const std::vector<int>& neighbor_particles_temp = cell_list[cell_temp];
            if (!neighbor_particles_temp.empty()) {
                for (const auto& neighbor_particle : neighbor_particles_temp) {
                    if (neighbor_particle == particle_seq) continue;
                    const std::vector<double>& neighbor_xyz = inputs[neighbor_particle];
                    std::vector<double> diff = get_min_diff(particle_xyz, neighbor_xyz);

                    double distance = 0;
                    for (int i = 0; i < 3; i++) {
                        distance += diff[i] * diff[i];
                    }
                    distance = sqrt(distance);

                    if (distance - rc < 1e-10) {
                        #pragma omp critical
                        {
                            particle_list[particle_seq].push_back(neighbor_particle);
                        }
                    }
                }
            }
        }
    }
}


    double rc;                                  // the cut-off radius
    std::vector<double> cube_size;              // cube domain的大小
    std::vector<int> lc;                        // 离散化后的cell尺寸
    std::vector<double> grid_size;              // cell的尺寸

    int num_particles;                          // 系统中的颗粒数量
    int num_cells;                              // 系统中的cell数量

    std::vector<std::vector<double>> particle_list;// 表示颗粒的邻接关系
    std::vector<std::vector<int>> cell_list;    // 表示cell的邻接关系
};

PYBIND11_MODULE(grid_based_nbl_mp, m) {
    py::class_<GB_NBL_MP_cube>(m, "GB_NBL_MP_cube")
        .def(py::init<std::vector<double>, int, double>())
        .def("constructor", &GB_NBL_MP_cube::constructor)
        .def("get_neighbors", &GB_NBL_MP_cube::get_neighbors);
}

