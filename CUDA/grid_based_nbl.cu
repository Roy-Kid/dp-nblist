#include "tools.cuh" // 引入array.cuh头文件

#include <random>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <iostream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/unique.h>
#include <thrust/sort.h>
#include <chrono> 
// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>

// namespace py = pybind11;


// // GPU kernel to add particles to corresponding cells
// __global__ void add_particles_to_cells_gpu(const int* particle_inds, Dynamic2DArray<int>& cell_list, int num_particles) {
//     int tid = blockIdx.x * blockDim.x + threadIdx.x;
//     if (tid < num_particles) {
//         int particle_ind = particle_inds[tid];
//         cell_list.push_back(particle_ind, tid);
//         // atomicAdd(&cell_list[particle_ind], tid);
//     }
// }

// // GPU kernel to build neighbor relationships considering rc
// __global__ void build_neighbor_relationships_gpu(const double* inputs, const int* particle_inds,
//                                                  const int* cell_list, int* particle_list,
//                                                  int num_particles, double cube_size_x, double cube_size_y, double cube_size_z,
//                                                  double grid_size_x, double grid_size_y, double grid_size_z,
//                                                  double rc) {
//     int tid = blockIdx.x * blockDim.x + threadIdx.x;
//     if (tid < num_particles) {
//         double particle_xyz[3];
//         particle_xyz[0] = inputs[tid];
//         particle_xyz[1] = inputs[tid + num_particles];
//         particle_xyz[2] = inputs[tid + 2 * num_particles];

//         int particle_ind = particle_inds[tid];
//         // ... (继续核函数中的逻辑)
//     }
// }

class GB_NBL_GPU_cube {
public:
    GB_NBL_GPU_cube(const thrust::host_vector<double>& c_size, int num_particles, double cut_off_radius)
        : rc(cut_off_radius), num_particles(num_particles)
    {
        this->cube_size = c_size;

        lc = thrust::device_vector<int>(3);
        lc[0] = static_cast<int>(cube_size[0] / cut_off_radius);
        lc[1] = static_cast<int>(cube_size[1] / cut_off_radius);
        lc[2] = static_cast<int>(cube_size[2] / cut_off_radius);

        grid_size = thrust::device_vector<double>(3);
        grid_size[0] = cube_size[0] / lc[0];
        grid_size[1] = cube_size[1] / lc[1];
        grid_size[2] = cube_size[2] / lc[2];

        num_cells = lc[0] * lc[1] * lc[2];

        particle_list.resize(num_particles, 100);
        cell_list.resize(num_cells, 100);
    }

    // Convert xyz to cell indices using CUDA (GPU)
    thrust::device_vector<int> xyz2ind_gpu(Dynamic2DArray<double>& xyz) {
        thrust::device_vector<int> cell_indices(xyz.max_rows);
        for (int i = 0; i < xyz.max_rows; ++i) {
            thrust::device_vector<int> cell_index(3);

            for (int j = 0; j < 3; j++) {
                double x = fmod(xyz.getElement(i, j), cube_size[j]); // Access data on GPU

                if (x < 0) x += cube_size[j];

                cell_index[j] = static_cast<int>(floor(x / grid_size[j]));
            }
            cell_indices[i] = cell_index[0] * lc[1] * lc[2] + cell_index[1] * lc[2] + cell_index[2];
        }
        return cell_indices;
    }


    // Convert cell indices to xyz using CUDA (GPU)
    Dynamic2DArray<double> ind2xyz_gpu(const thrust::device_vector<int>& ind) {
        Dynamic2DArray<double> xyz(ind.size(), 3);

        for (int i = 0; i < ind.size(); ++i) {
            thrust::device_vector<double> point(3);
            int index = ind[i];

            point[0] = (index / (lc[1] * lc[2])) * grid_size[0] + grid_size[0] / 2;
            point[1] = (index / lc[2]) * grid_size[1] + grid_size[1] / 2;
            point[2] = (index % lc[2]) * grid_size[2] + grid_size[2] / 2;

            xyz.push_back(i, fmod(point[0], cube_size[0]));
            xyz.push_back(i, fmod(point[1], cube_size[1]));
            xyz.push_back(i, fmod(point[2], cube_size[2]));

            
            // std::cout << "___________" << std::endl;
            // std::cout << "0: " << xyz.helper[0] << " 1: " << xyz.helper[1] << " 2: " << xyz.helper[2] << std::endl;
            // xyz.print();
            // break;
        }
        return xyz;
    }

    // Get the minimum difference between two points using CUDA (GPU)
    thrust::device_vector<double> get_min_diff_gpu(const thrust::device_vector<double>& xyz1, const thrust::device_vector<double>& xyz2) {
        thrust::device_vector<double> difference(3);
        for (int i = 0; i < 3; i++) {
            double diff = xyz2[i] - xyz1[i];
            diff = std::fmod(diff + cube_size[i] / 2, cube_size[i]);
            if (diff < 0) diff += cube_size[i];
            diff -= cube_size[i] / 2;
            difference[i] = diff;
        }
        return difference;
    }

    // Get neighbor cells using CUDA (GPU)
    thrust::device_vector<int> get_neighbor_cells_gpu(int cell_ind) {
        thrust::device_vector<int> cell_ind_list(1);
        cell_ind_list[0] = cell_ind;
        Dynamic2DArray<double> xyz = ind2xyz_gpu(cell_ind_list);

        thrust::device_vector<int> adjacent_cells;
        for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
                for (int dk = -1; dk <= 1; dk++) {
                    double ni = static_cast<int>(xyz.getElement(0, 0)) + di * static_cast<int>(grid_size[0]);
                    double nj = static_cast<int>(xyz.getElement(0, 1)) + dj * static_cast<int>(grid_size[1]);
                    double nk = static_cast<int>(xyz.getElement(0, 2)) + dk * static_cast<int>(grid_size[2]);
                    Dynamic2DArray<double> point(1, 3);
                    point.push_back(0, ni);
                    point.push_back(0, nj);
                    point.push_back(0, nk);
                    thrust::device_vector<int> ind_temp = xyz2ind_gpu(point);

                    adjacent_cells.push_back(ind_temp[0]);
                }
            }
        }

        thrust::device_vector<int> result = vector2set_int(adjacent_cells);

        return result;
    }


    thrust::device_vector<int> get_neighbors(int particle_seq) {
        int neighbor_njm = particle_list.helper[particle_seq];
        thrust::device_vector<int> res(neighbor_njm);

        res = particle_list.getVector(particle_seq);   

        return res;
    }

    // Convert constructor to gpu with two different GPU kernel
    void constructor_gpu(Dynamic2DArray<double>& inputs) {
        // Step 1: Convert xyz to cell indices using CUDA (GPU)
        thrust::device_vector<int> particle_inds = xyz2ind_gpu(inputs);

        // // Step 2: Add particles to corresponding cells
        // for (int i = 0; i < num_particles; i++) {
        //     cell_list.push_back(particle_inds[i], i);
        //     // std::cout << "temp" << particle_inds[i] << " " << i << std::endl;
        // }

        // // Step 2: Add particles to corresponding cells
        // thrust::device_vector<int> device_cell_list(num_cells, -1);
        // add_particles_to_cells_gpu<<<num_blocks, block_size>>>(particle_inds.data().get(),
        //                                                        cell_list,
        //                                                        num_particles);

        // // Step 3: Build neighbor relationships considering rc
        // for (int particle_seq = 0; particle_seq < num_particles; particle_seq++) {
        //     thrust::device_vector<double> particle_xyz = inputs.getVector(particle_seq);

        //     // int particle_ind = particle_inds[particle_seq];
        //     // thrust::device_vector<int> adjacent_cells = get_neighbor_cells_gpu(particle_ind);

        //     // for (int i = 0; i < adjacent_cells.size(); i++) {
        //     //     thrust::device_vector<int> neighbor_particles_temp = cell_list.getVector(adjacent_cells[i]);
        //     //     if (neighbor_particles_temp.size()) {
        //     //         for (int j = 0; j < neighbor_particles_temp.size(); j++) {
        //     //             int neighbor_particle_seq = neighbor_particles_temp[j];
        //     //             if (neighbor_particle_seq == particle_seq) continue;
        //     //             thrust::device_vector<double> neighbor_xyz = inputs.getVector(neighbor_particle_seq);
        //     //             thrust::device_vector<double> diff = get_min_diff_gpu(particle_xyz, neighbor_xyz);

        //     //             double distance = 0;
        //     //             for (int k = 0; k < 3; k++) {
        //     //                 distance += diff[i] * diff[i];
        //     //             }

        //     //             if (distance - rc * rc < 1e-10) {
        //     //                 particle_list.push_back(particle_seq, neighbor_particle_seq);
        //     //             }
        //     //         }
        //     //     }
        //     // }
        // }
        // // Step 3: Build neighbor relationships considering rc
        // thrust::device_vector<int> device_particle_list(num_particles, -1);
        // build_neighbor_relationships_gpu<<<num_blocks, block_size>>>(inputs.data().get(),
        //                                                              particle_inds.data().get(),
        //                                                              device_cell_list.data().get(),
        //                                                              device_particle_list.data().get(),
        //                                                              num_particles,
        //                                                              cube_size[0], cube_size[1], cube_size[2],
        //                                                              grid_size[0], grid_size[1], grid_size[2],
        //                                                              rc);

        // // Copy the result back to the host
        // thrust::copy(device_cell_list.begin(), device_cell_list.end(), cell_list.begin());
        // thrust::copy(device_particle_list.begin(), device_particle_list.end(), particle_list.begin());
    }

    double rc;                                  // the cut-off radius
    int num_particles;                          // 系统中的颗粒数量
    int num_cells;                              // 系统中的cell数量

    thrust::device_vector<double> cube_size;     // cube domain的大小
    thrust::device_vector<int> lc;              // 离散化后的cell尺寸
    thrust::device_vector<double> grid_size;    // cell的尺寸

    Dynamic2DArray<int> particle_list; // 使用 Dynamic2DArray 表示 particle_list和cell_list
    Dynamic2DArray<int> cell_list; 

    // const int num_blocks = 12;
    // const int block_size = 1024;
};


double generate_random_double(double min_val, double max_val) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min_val, max_val);
    return dis(gen);
}


// 绑定 Dynamic2DArray 类模板到 Python
template <typename T>
void bind_Dynamic2DArray(py::module& m, const std::string& name) {
    // 定义中间函数来明确模板类型
    auto class_fn = [](py::module& m, const std::string& name) {
        return py::class_<Dynamic2DArray<T>>(m, name)
            .def(py::init<>())
            .def(py::init<int, int>())
            .def("push_back", &Dynamic2DArray<T>::push_back)
            .def("getElement", &Dynamic2DArray<T>::getElement)
            .def("getVector", &Dynamic2DArray<T>::getVector)
            .def("resize", &Dynamic2DArray<T>::resize)
            .def("print", (void (Dynamic2DArray<T>::*)()) &Dynamic2DArray<T>::print)
            .def("print", (void (Dynamic2DArray<T>::*)(int)) &Dynamic2DArray<T>::print);
    };

    // 调用中间函数
    class_fn(m, name);
}

PYBIND11_MODULE(grid_based_nbl_cpp, m) {
    // // Dynamic2DArray类模板的Pybind11绑定
    // bind_Dynamic2DArray<double>(m, "Dynamic2DArrayDouble");
    // bind_Dynamic2DArray<int>(m, "Dynamic2DArrayInt");

    // // 辅助函数vector2set_int的Pybind11绑定
    // m.def("vector2set_int", &vector2set_int, "Convert a vector to a set (remove duplicates) using Thrust.");

    // GB_NBL_GPU_cube类的Pybind11绑定
    py::class_<GB_NBL_GPU_cube>(m, "GB_NBL_GPU_cube")
        .def(py::init<const thrust::host_vector<double>&, int, double>())
        .def("constructor_gpu_wrapper", &GB_NBL_GPU_cube::constructor_gpu)
        .def("get_neighbors", &GB_NBL_GPU_cube::get_neighbors, py::arg("particle_seq"));
}




int main() {
    // 测试 xyz - index转换
    thrust::host_vector<double> c_size(3);
    for (int i = 0; i < 3; i++) {
        c_size[i] = 10.0; 
    }

    int num_particles = 3;
    double cut_off_radius = 1.0;

    // Call GB_NBL_GPU_cube constructor with thrust::device_vector
    GB_NBL_GPU_cube lc_cube(c_size, num_particles, cut_off_radius);

    // Convert xyz to thrust::device_vector of Dynamic2DArray
    std::vector<std::vector<double>> xyz_temp = { {1.5, 2.2, 3.7}, {4.1, 5.9, 6.3}, {7.8, 8.4, 9.2} };
    Dynamic2DArray<double> xyz(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            xyz.push_back(i, xyz_temp[i][j]);
        }
    }

    // Call the GPU wrapper functions to perform computations on the GPU
    thrust::device_vector<int> cell_indices = lc_cube.xyz2ind_gpu(xyz);
    thrust::host_vector<int> host_cell_indices = cell_indices;
    std::cout << "Cell indices:" << std::endl;
    for (const auto& index : host_cell_indices) {
        std::cout << index << " ";
    }
    std::cout << std::endl;


    // Call the GPU wrapper function for ind2xyz_gpu
    std::cout << "Restore: " << std::endl;
    Dynamic2DArray<double> device_restored_xyz = lc_cube.ind2xyz_gpu(cell_indices);
    device_restored_xyz.print();
    thrust::device_vector<double> temp = device_restored_xyz.getVector(0);
    for (int i = 0; i < temp.size(); i++) {
        std::cout << temp[i] << " ";
    }
    std::cout << "sighn" << std::endl;

    // Check the distance calculation
    std::vector<double> cube_size2 = { 20, 20, 20 };
    thrust::device_vector<double> cube_size2_dev(3);
    for(int i = 0; i < 3; i++) {
        cube_size2_dev[i] = cube_size2[i];
    }

    int num_particles2 = 800;
    double cut_off_radius2 = 2.0;
    GB_NBL_GPU_cube lc_cube2(cube_size2, num_particles2, cut_off_radius2);

    std::vector<double> p1 = {19.608,  7.38 ,  6.096}; 
    thrust::device_vector<double> p1_dev(3);
    for(int i = 0; i < 3; i++) {
        p1_dev[i] = p1[i];
    }

    std::vector<double> p2 = {0.   , 9.012, 6.944};
    thrust::device_vector<double> p2_dev(3);
    for(int i = 0; i < 3; i++) {
        p2_dev[i] = p2[i];
    }

    thrust::device_vector<double> res = lc_cube2.get_min_diff_gpu(p1, p2);
    std::cout << "Real value: [" << 0.392 << " " << 1.632 << " " << 0.848 << "]" << std::endl;
    std::cout << "Calculated value: [" << res[0] << " " << res[1] << "" << res[2] << "]" << std::endl;
    


    // Check get_neighbor_cells_gpu
    thrust::host_vector<double> c_size(3);
    for (int i = 0; i < 3; i++) {
        c_size[i] = 80.0; 
    }

    int num_particles = 5000;
    double cut_off_radius = 2.0;

    // Call GB_NBL_GPU_cube constructor with thrust::device_vector
    GB_NBL_GPU_cube lc_cube3(c_size, num_particles, cut_off_radius);

    // thrust::device_vector res = lc_cube3.get_neighbor_cells_gpu(1);
    // for(int i = 0; i < res.size(); i++) {
    //     std::cout << res[i] << " ";
    // }

    // Check the constructor
    std::vector<std::vector<double>> xyz_temp;
    for (int i = 0; i < num_particles; i++) {
        double x = generate_random_double(0.0, 80.0);
        double y = generate_random_double(0.0, 80.0);
        double z = generate_random_double(0.0, 80.0);
        xyz_temp.push_back({x, y, z});
    }

    // 将 xyz_temp 中的坐标导入到 Dynamic2DArray<double> xyz 中
    Dynamic2DArray<double> xyz(num_particles, 3);
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < 3; j++) {
            xyz.push_back(i, xyz_temp[i][j]);
        }
    }

    std::cout << "Data Preparation Ready!" << std::endl;

    // 统计时间：开始
    auto start_time = std::chrono::high_resolution_clock::now();
    // 执行命令
    lc_cube3.constructor_gpu(xyz);
    // 统计时间：结束
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    // 输出执行时间（以秒为单位）
    std::cout << "执行时间：" << duration.count() << " 秒" << std::endl;




    // 测试get_neighbors
    thrust::device_vector res_neighbor = lc_cube3.get_neighbors(0);
    for (int i = 0; i < res_neighbor.size(); i++){
        std::cout << res_neighbor[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Thrust version: " << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << std::endl;


    return 0;
}