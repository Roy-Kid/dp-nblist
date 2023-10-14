#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class GB_NBL_cube {
public:
    GB_NBL_cube(std::vector<double> cube_size, int num_particles, double cut_off_radius)
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
        for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
                for (int dk = -1; dk <= 1; dk++) {
                    double ni = static_cast<int>(xyz[0][0]) + di * static_cast<int>(grid_size[0]);
                    double nj = static_cast<int>(xyz[0][1]) + dj * static_cast<int>(grid_size[1]);
                    double nk = static_cast<int>(xyz[0][2]) + dk * static_cast<int>(grid_size[2]);
                    std::vector<int> ind_temp = xyz2ind({{ni, nj, nk}});
                    adjacent_cells.insert(ind_temp[0]);
                }
            }
        }
        return adjacent_cells;
    }

    std::vector<double> get_neighbors(int particle_seq) {
        return particle_list[particle_seq];
    }

    void constructor(const std::vector<std::vector<double>>& inputs) {
        // 把粒子添加到对应的cell之中
        std::vector<int> particle_inds = xyz2ind(inputs);

        for (int i = 0; i < num_particles; i++) {
            cell_list[particle_inds[i]].push_back(i);
        }

        // 建立链接关系，考虑rc
        for (int particle_seq = 0; particle_seq < num_particles; particle_seq++) {
            std::vector<double> particle_xyz = inputs[particle_seq];

            int particle_ind = particle_inds[particle_seq];
            std::unordered_set<int> adjacent_cells = get_neighbor_cells(particle_ind);

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
                        // std::cout << distance << std::endl;
                        // std::cout << std::endl;
                        // std::cout << "particle seq: " << particle_seq << " "
                        // << "particle neighbor: " << neighbor_particle << " "
                        // << "distance: " << distance << std::endl;

                        if (distance - rc < 1e-10) {
                            particle_list[particle_seq].push_back(neighbor_particle);
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


PYBIND11_MODULE(grid_based_nbl_cpp, m) {
    py::class_<GB_NBL_cube>(m, "GB_NBL_cube")
        .def(py::init<std::vector<double>, int, double>())
        .def("constructor", &GB_NBL_cube::constructor)
        .def("get_neighbors", &GB_NBL_cube::get_neighbors);
}



// int main() {
//     // 测试 xyz - index转换
//     std::vector<double> cube_size = { 4.0, 4.0, 4.0 };
//     int num_particles = 6;
//     double cut_off_radius = 2.0;
//     GB_NBL_cube lc_cube(cube_size, num_particles, cut_off_radius);

//     std::vector<std::vector<double>> xyz = { {1.5, 2.2, 3.7}, {4.1, 5.9, 6.3}, {7.8, 8.4, 9.2},
//                                              {1.5, 2.2, 3.7}, {4.1, 5.9, 6.3}, {7.8, 8.4, 9.2} };

//     std::vector<int> cell_indices;
//     std::vector<std::vector<double>> xyz_restored;

//     cell_indices = lc_cube.xyz2ind(xyz);
//     std::cout << "Cell indices:" << std::endl;
//     for (const auto& index : cell_indices) 
//     {
//         std::cout << index << " ";
//     }
//     std::cout << std::endl;

//     std::vector<std::vector<double>> restored_xyz = lc_cube.ind2xyz(cell_indices);
//     std::cout << "Restored XYZ:" << std::endl;
//     for (const auto& points : restored_xyz) {
//         for (const auto& point : points) {
//             std::cout << point << " ";
//         }
//         std::cout << std::endl;
//     }

//     cell_indices = lc_cube.xyz2ind(restored_xyz);
//     std::cout << "Restored cell indices:" << std::endl;
//     for (const auto& index : cell_indices) 
//     {
//         std::cout << index << " ";
//     }
//     std::cout << std::endl;

//     // 测试get_neighbor_cells
//     int cell_ind = 0;  // 替换为要测试的 cell 索引
//     std::unordered_set<int> neighbors = lc_cube.get_neighbor_cells(cell_ind);

//     // 输出结果
//     std::cout << "Neighbor cells of cell " << cell_ind << ":" << std::endl;
//     for (const auto& cell : neighbors) {
//         std::cout << cell << " ";
//     }
//     std::cout << std::endl;


//     // 测试距离计算
//     std::vector<double> cube_size2 = { 20, 20, 20 };
//     int num_particles2 = 800;
//     double cut_off_radius2 = 2.0;
//     GB_NBL_cube lc_cube2(cube_size2, num_particles2, cut_off_radius2);

//     std::vector<double> p1 = {19.608,  7.38 ,  6.096}; 
//     std::vector<double> p2 = {0.   , 9.012, 6.944};
//     std::vector<double> res = lc_cube2.get_min_diff(p1, p2);
//     std::cout << "Real value: [" << 0.392 << " " << 1.632 << " " << 0.848 << "]" << std::endl;
//     std::cout << "Calculated value: [" << res[0] << " " << res[1] << "" << res[2] << "]" << std::endl;
    
//     // 测试constructor
//     lc_cube.constructor(xyz);

    
//     // 测试get_neighbors
//     for (int i = 0; i < 6; i++){
//         for (auto particle_a: lc_cube.get_neighbors(i)) {
//             std::cout << particle_a << " ";
//         }
//         std::cout << std::endl;
//     }

//     return 0;
// }

