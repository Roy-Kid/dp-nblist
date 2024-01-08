#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include "cuda_nblist.h"
#include "box.h"
#include "vec3.cpp"

// CUDA相关头文件
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace dpnblist
{
    /*
    * to get the neighbor cell index of each cell, and store them in nebcell_list
    * param: offset_vec: the offset vector of each neighbor cell
    * param: ncells: the number of cells
    * param: n_nebcells: the number of neighbor cells
    * param: L: the box length of each dimension
    */
    __global__ void get_neb(int* offset_vec, size_t *ncells, size_t *n_nebcells, int *L, int* nebcell) {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;

        if (tid < *ncells) {
            // int L[3] = {4, 4, 4};
            int cell_veci[3];
            int offset_veci[3];
            float neb_vec[3];
            float wrap_neb_vec[3];
            int round_wrap_neb_vec[3];
            float f[3];
            float wrapped_f[3];

            cell_veci[0] = tid / (L[1] * L[2]);
            cell_veci[1] = (tid - cell_veci[0] * L[1] * L[2]) / L[2];
            cell_veci[2] = tid - cell_veci[0] * L[1] * L[2] - cell_veci[1] * L[2];

            for (int i = 0; i < *n_nebcells; ++i) {
                for (int j = 0; j < 3; ++j) {
                    offset_veci[j] = offset_vec[i * 3 + j];
                    neb_vec[j] = cell_veci[j] + offset_veci[j];
                }
                for (int k = 0; k < 3; ++k) {
                    f[k] = neb_vec[k] / L[k];
                    wrapped_f[k] = f[k] - std::floor(f[k] + 0.000001);
                    wrap_neb_vec[k] = wrapped_f[k] * L[k];
                    round_wrap_neb_vec[k] = round(wrap_neb_vec[k]);
                }
                
                nebcell[*n_nebcells * tid + i] = round_wrap_neb_vec[0] * L[1] * L[2] + round_wrap_neb_vec[1] * L[2] + round_wrap_neb_vec[2];
            }
        }
    }
    /*
    * to build the linked list, and store the head, lscl, atom_cellindex, cell_count
    * param: xyz_dev: the coordinates of all atoms
    * param: head_dev: the head atom of each cell
    * param: lscl_dev: the next atom of each cell
    * param: atom_cellindex: the cell index of each atom
    * param: cell_count_dev: the number of atoms in each cell
    * param: r_cutoff: the cutoff radius
    * param: L: the box length of each dimension
    * param: natoms: the number of atoms
    */
    __global__ void build_linked_list_kernel(float *xyz_dev, int *head_dev, int *lscl_dev, int *atom_cellindex, int *cell_count_dev, float *_r_cutoff, int *L, size_t *natoms) {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < *natoms) {
            int cell_index_vec[3];
            for (int j = 0; j < 3; j++) {
                cell_index_vec[j] = xyz_dev[i * 3 + j] / *_r_cutoff;
                if (cell_index_vec[j] == L[j])
                    cell_index_vec[j] = cell_index_vec[j] - 1;
            }
            int cell_index = cell_index_vec[0] * L[1] * L[2] + cell_index_vec[1] * L[2] + cell_index_vec[2];
            atom_cellindex[i] = cell_index;
            lscl_dev[i] = atomicExch(&head_dev[cell_index], i);
            atomicAdd(&cell_count_dev[cell_index], 1);
        }
    }

    /*
    * to calculate the distance between two atoms, and return the square of the distance
    * param: pos_i: the coordinates of atom i
    * param: pos_j: the coordinates of atom j
    * param: length: the box length of each dimension
    */
    __device__ float calc_distance(float *pos_i, float *pos_j, float *length) {
        float difference[3] = {0.0, 0.0, 0.0};
        float diff = 0.0;

        for (int i = 0; i < 3; i++) {
            float dri = pos_i[i] - pos_j[i];
            diff = fmodf((dri + length[i] / 2), length[i]);
            if (diff < 0) diff += length[i];
            diff -= length[i] / 2;
            difference[i] = diff;
        }

        return (difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2]);
    }

    /*
    * to build the neighbor list array, and store the neighbor list of each atom to neighborListArray
    * param: natoms: the number of atoms
    * param: xyz: the coordinates of all atoms
    * param: head: the head atom of each cell
    * param: lscl: the next atom of each cell
    * param: atom_cellindex: the cell index of each atom
    * param: cell_count: the number of atoms in each cell
    * param: nebcell: the neighbor cell index of each cell
    * param: L: the box length of each dimension
    * param: box_length: the cell number of each dimension
    * param: r_cutoff2: the square of cutoff radius
    * param: neighborListArray: the neighbor list of each atom
    */
    __global__ void  buildListArray(size_t *natoms, float *xyz, int *head, int *lscl, int *atom_cellindex, int *cell_atoms_count, int *nebcell, int *L, float *box_length, float *r_cutoff2, int *neighborListArray){
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid < *natoms){
            size_t count = 0;                                                   //记录这是第几个邻居原子
            float pos_i[3] = {xyz[tid*3], xyz[tid*3+1], xyz[tid*3+2]};          //原子坐标
            size_t cell_index = atom_cellindex[tid];                            //原子属于哪个单元
            size_t neighbor_cell = 0;                                           //单元的邻居单元索引
            for (int i = cell_index * 27; i < cell_index * 27 + 27; ++i)        // Loop through the neighboring cells
            {                                                                   
                neighbor_cell = nebcell[i];                                     //27个邻居单元索引
                size_t num_atoms = cell_atoms_count[neighbor_cell];             //单元有几个原子
                size_t atom_id = head[neighbor_cell];                           //单元的第一个原子
                for (size_t j = 0; j < num_atoms; ++j)                          //循环单元内所有的原子
                {
                    float pos_j[3] = {xyz[atom_id*3], xyz[atom_id*3+1], xyz[atom_id*3+2]};
                    float r = calc_distance(pos_i, pos_j, box_length);
                    if (r < *r_cutoff2)
                    {
                        neighborListArray[tid * 100 + count] = atom_id;                   // Add j to the neighbor list of i
                        count++;
                    }
                    atom_id = lscl[atom_id];
                }
            }
        }
    }

    /*
    * the constructor of class CudaCellList, to initialize the parameters
    * param: box: the box of the system
    * param: r_cutoff: the cutoff radius
    * param: skin: the skin of the system
    * initialize the parameters of the system and allocate the memory of the device
    * param: r_cutoff: the cutoff radius
    * param: skin: the skin of the system
    * param: d_box_len: the box length of each dimension
    * param: d_cell_len: the cell number of each dimension
    * param: d_ncells: the number of cells
    */
    
    CudaCellList::CudaCellList(Box *box, float r_cutoff, float skin):_box(box), _r_cutoff(r_cutoff+skin),_skin(skin)
	{
	  Vec3<float> box_length = box->get_lengths();
        for (int i = 0; i < 3; ++i) {
            _box_len[i] = box_length[i];
        }
        //_r_cutoff = r_cutoff + skin;
        cudaMalloc((void**)&d_r_cutoff, sizeof(float));
        cudaMemcpy(d_r_cutoff, &_r_cutoff, sizeof(float), cudaMemcpyHostToDevice);

        //_skin = skin;
        cudaMalloc((void**)&d_skin, sizeof(float));
        cudaMemcpy(d_skin, &_skin, sizeof(float), cudaMemcpyHostToDevice);

        for (int i = 0; i < 3; ++i) {
            _cell_len[i] = _box_len[i]/_r_cutoff;
        }
        _ncells = _cell_len[0] * _cell_len[1] * _cell_len[2];

        cudaMalloc((void**)&d_box_len, 3 * sizeof(float));
        cudaMemcpy(d_box_len, _box_len, 3 * sizeof(float), cudaMemcpyHostToDevice);

        cudaMalloc((void**)&d_cell_len, 3 * sizeof(int));
        cudaMemcpy(d_cell_len, _cell_len, 3 * sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc((void**)&d_ncells, sizeof(size_t));
        cudaMemcpy(d_ncells, &_ncells, sizeof(size_t), cudaMemcpyHostToDevice);

    }

    /*
    * to build the neighbor list array, and initialize the parameters and allocate the memory of the device
    * variable: d_natoms: the number of atoms
    * variable: d_off_set_vec: the offset vector of each neighbor cell
    * variable: d_n_nebcells: the number of neighbor cells
    * variable: d_nebcell_list: the neighbor cell index of each cell
    * run kernel get_neb to get the neighbor cell index of each cell
    * run update(xyz) to build the neighbor list array
    */
    void CudaCellList::build(std::vector<std::vector<float>> &xyz) {
        auto start_time = std::chrono::high_resolution_clock::now();

        _natoms = xyz.size();
        cudaMalloc((void**)&d_natoms, sizeof(size_t));
        cudaMemcpy(d_natoms, &_natoms, sizeof(size_t), cudaMemcpyHostToDevice);

        
        std::vector<std::vector<int>> off_set_vec = {
            {-1, -1, -1},{-1, -1, 0},{-1, -1, 1},{-1, 0, -1},{-1, 0, 0},{-1, 0, 1},{-1, 1, -1},{-1, 1, 0},{-1, 1, 1},
            {0, -1, -1},{0, -1, 0},{0, -1, 1},{0, 0, -1},{0, 0, 0},{0, 0, 1},{0, 1, -1},{0, 1, 0},{0, 1, 1},
            {1, -1, -1},{1, -1, 0},{1, -1, 1},{1, 0, -1},{1, 0, 0},{1, 0, 1},{1, 1, -1},{1, 1, 0},{1, 1, 1}
        };
        size_t n_nebcells = off_set_vec.size();
        size_t *d_n_nebcells;
        cudaMalloc((void**)&d_n_nebcells, sizeof(size_t));
        cudaMemcpy(d_n_nebcells, &n_nebcells, sizeof(size_t), cudaMemcpyHostToDevice);

        int off_set_vec_1d[n_nebcells * 3];
        for (int i = 0; i < n_nebcells; ++i) {
            for (int j = 0; j< 3; ++j) {
                off_set_vec_1d[i * 3 + j] = off_set_vec[i][j];
            }
        }
        int *d_off_set_vec_1d;
        cudaMalloc((void**)&d_off_set_vec_1d, n_nebcells * 3 * sizeof(int));
        cudaMemcpy(d_off_set_vec_1d, off_set_vec_1d, n_nebcells * 3 * sizeof(int), cudaMemcpyHostToDevice);
        
        cudaMalloc((void**)&d_nebcell_list, _ncells * n_nebcells * sizeof(int));
        int threadsPerBlock = 256;
        int blocksPerGrid = (_ncells + threadsPerBlock - 1) / threadsPerBlock;
        //////////
        auto neb_start = std::chrono::high_resolution_clock::now();
        get_neb<<<blocksPerGrid, threadsPerBlock>>>(d_off_set_vec_1d, d_ncells, d_n_nebcells, d_cell_len, d_nebcell_list);
        auto neb_end = std::chrono::high_resolution_clock::now();
        auto neb_duration = std::chrono::duration_cast<std::chrono::microseconds>(neb_end - neb_start);
        std::cout << "Time taken by get nebcell: " << neb_duration.count() << " microseconds" << std::endl;
        //////////
        cudaFree(d_off_set_vec_1d);

        // int *d_head, *d_lscl, *d_atom_cellindex, *d_cell_atoms_count;
        cudaMalloc((void**)&d_head, _ncells * sizeof(int));                      // head[i] is the first atom in cell i
        cudaMalloc((void**)&d_lscl, _natoms * sizeof(int));                      // lscl is the atom linked list, lscl[i] is the next atom in the cell
        cudaMalloc((void**)&d_atom_cellindex, _natoms * sizeof(int));            // atom in which cell
        cudaMalloc((void**)&d_cell_atoms_count, _ncells * sizeof(int));          // cell_atoms_count[i] is the number of atoms in cell i
        
        _neighborListArray = new int[_natoms * 100];
        update(xyz);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        time_cost = duration.count();
        std::cout << "Time taken by cuda: " << time_cost << " milliseconds" << std::endl;
    }

    /*
    * to update the neighbor list array, and allocate the memory of the device
    * variable: d_xyz_1d: the coordinates of all atoms
    * variable: d_head: the head atom of each cell
    * variable: d_lscl: the next atom of each cell
    * variable: d_atom_cellindex: the cell index of each atom
    * variable: d_cell_atoms_count: the number of atoms in each cell
    * run kernel build_linked_list_kernel to build the linked list
    * run kernel buildListArray to build the neighbor list array
    */
    void CudaCellList::update(std::vector<std::vector<float>> &xyz) {

        std::vector<Vec3<float>> wrap_xyz;
        for (int i = 0; i < _natoms; ++i) {
            wrap_xyz.emplace_back(xyz[i][0], xyz[i][1], xyz[i][2]);
        }
	    wrap_xyz = _box->wrap(wrap_xyz);
        float *xyz_1d = new float[_natoms * 3];
        // float xyz_1d[_natoms*3];
        for (int i = 0; i < _natoms; ++i) {
            for (int j = 0; j< 3; ++j) {
                xyz_1d[i * 3 + j] = wrap_xyz[i][j];
            }
        }

        float *d_xyz_1d;
        cudaMalloc((void**)&d_xyz_1d, _natoms * 3 * sizeof(float));
        cudaMemcpy(d_xyz_1d, xyz_1d, _natoms * 3 * sizeof(float), cudaMemcpyHostToDevice);

        // int *d_head, *d_lscl, *d_atom_cellindex, *d_cell_atoms_count;
        _cell_atoms_count = new int[_ncells];
        for (int i = 0; i < _ncells; ++i) {
            _cell_atoms_count[i] = 0;
        }
        cudaMemcpy(d_cell_atoms_count, _cell_atoms_count, _ncells * sizeof(int), cudaMemcpyHostToDevice);
        // auto start_time = std::chrono::high_resolution_clock::now();
        int threadsPerBlock = 256;
        int blocksPerGrid = (_natoms + threadsPerBlock - 1) / threadsPerBlock;

        //////////
        auto linked_list_start = std::chrono::high_resolution_clock::now();
        build_linked_list_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_xyz_1d, d_head, d_lscl, d_atom_cellindex, d_cell_atoms_count, d_r_cutoff, d_cell_len, d_natoms);
        auto linked_list_end = std::chrono::high_resolution_clock::now();
        auto linked_list_duration = std::chrono::duration_cast<std::chrono::microseconds>(linked_list_end - linked_list_start);
        std::cout << "Time taken by build_linked_list: " << linked_list_duration.count() << " microseconds" << std::endl;
        //////////

        float r_cutoff2 = _r_cutoff * _r_cutoff;
        float *d_r_cutoff2;
        cudaMalloc((void**)&d_r_cutoff2, sizeof(float));
        cudaMemcpy(d_r_cutoff2, &r_cutoff2, sizeof(float), cudaMemcpyHostToDevice);

        int *temp_neighborListArray = new int[_natoms * 100]; // neighborListArray is the neighbor list of atoms, has a maximum 100 neighbors for each atom
        for (int i = 0; i < _natoms * 100; ++i) {
            temp_neighborListArray[i] = -1;
        }

        int *d_neighborListArray;
        cudaMalloc((void**)&d_neighborListArray, _natoms * 100 * sizeof(int));
        cudaMemcpy(d_neighborListArray, temp_neighborListArray, _natoms * 100 * sizeof(int), cudaMemcpyHostToDevice);

        blocksPerGrid = (_natoms + threadsPerBlock - 1) / threadsPerBlock;

        //////////
        auto build_start = std::chrono::high_resolution_clock::now();
        buildListArray<<<blocksPerGrid, threadsPerBlock>>>(d_natoms, d_xyz_1d, d_head, d_lscl, d_atom_cellindex, d_cell_atoms_count, d_nebcell_list, d_cell_len, d_box_len, d_r_cutoff2, d_neighborListArray);
        auto build_end = std::chrono::high_resolution_clock::now();
        auto build_duration = std::chrono::duration_cast<std::chrono::microseconds>(build_end - build_start);
        std::cout << "Time taken by build list array: " << build_duration.count() << " microseconds" << std::endl;
        //////////

        //////////
        auto copy_start = std::chrono::high_resolution_clock::now();
        cudaMemcpy(_neighborListArray, d_neighborListArray, _natoms * 100 * sizeof(int), cudaMemcpyDeviceToHost);
        auto copy_end = std::chrono::high_resolution_clock::now();
        auto copy_duration = std::chrono::duration_cast<std::chrono::milliseconds>(copy_end - copy_start);
        std::cout << "Time taken by Memcopy: " << copy_duration.count() << " milliseconds" << std::endl;
        //////////

        cudaFree(d_r_cutoff2);
        cudaFree(d_xyz_1d);
        cudaFree(d_neighborListArray);
        delete[] xyz_1d;
        delete[] temp_neighborListArray;
        // out();
    }
    // output the neighbor list array to a file
    void CudaCellList::out(){
        std::ofstream outfile;
        outfile.open("neb_list.txt");
        int neighborListArrayj = 0;
        for (int i = 0; i < _natoms; ++i) {
        // for (int i = 0; i < 10; ++i) {
            outfile << (i+1);
            for(int j = 0; j < 100; ++j) {
                neighborListArrayj = _neighborListArray[i * 100 + j];
                // outfile << "  " << (neighborListArrayj+1);
                if (neighborListArrayj >= 0 && neighborListArrayj != i) outfile << "\t" << (neighborListArrayj+1);
            }
            outfile << std::endl;
        }
        outfile.close();
    }

    // return the neighbor list array
    std::vector<std::vector<size_t>> CudaCellList::get_listArray() {
        std::vector<std::vector<size_t>> neighborListArray;
        for (int i = 0; i < _natoms; ++i) {
            std::vector<size_t> neighborListArrayi;
            for(int j = 0; j < 100; ++j) {
                int neighborListArrayj = _neighborListArray[i * 100 + j];
                if (neighborListArrayj >= 0 && neighborListArrayj != i) neighborListArrayi.push_back(neighborListArrayj);
            }
            neighborListArray.push_back(neighborListArrayi);
        }
        return neighborListArray;
    }

    // return the total time cost of the neighbor list array
    size_t CudaCellList::gettime() {
        return time_cost;
    }

    CudaCellList::~CudaCellList()
    {
        cudaFree(d_box_len);
        cudaFree(d_cell_len);
        cudaFree(d_ncells);
        cudaFree(d_natoms);
        cudaFree(d_nebcell_list);
        cudaFree(d_head);
        cudaFree(d_lscl);
        cudaFree(d_atom_cellindex);
        cudaFree(d_cell_atoms_count);
        delete[] _cell_atoms_count;
        delete[] _neighborListArray;
    }


    base_NBL* NeighborList::createNeighborList(std::string type, Box *box, float r_cutoff, float skin)
    {
        Vec3<float> vec3_box_length = box->get_lengths();
        std::vector<float> box_length = {vec3_box_length[0], vec3_box_length[1], vec3_box_length[2]};
        if (type == "celllist") {
            //return new CudaCellList(box_length, r_cutoff, skin);
            return new CudaCellList(box, r_cutoff, skin);
        }
        else {
            std::cout << "type error" << std::endl;
            return NULL;
        }
    }

    NeighborList::~NeighborList()
    {
    }

}//namespace dpnblist

