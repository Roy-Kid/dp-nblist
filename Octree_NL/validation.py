#%%
import numpy as np
# 导入编译后的C++模块
import grid_based_nbl_gpu

import asap3
import numpy as np
import ase
from ase.io import read, write, Trajectory
from pathlib import Path
import time


# glob all pdb in current folder
for pdb in Path('.').glob('*.pdb'):
    if pdb.stem == 'lj':
        continue
    print(pdb)
    # read pdb
    atoms = read(pdb)
    # get the name of pdb
    name = pdb.stem
    # get the positions of atoms
    xyz = atoms.arrays['positions']
    # set the size of cell
    atoms.set_cell(np.diag([20, 20, 20]))
    atoms.set_pbc([True, True, True])

    # get_neighborlist
    cutoff = 2.0  # cutoff radius      

    # 记录开始时间
    start_time = time.time()
    nl = asap3.FullNeighborList(cutoff, atoms)

    pair_i_idx = []
    pair_j_idx = []
    n_diff = []
    for i in range(len(atoms)):
        indices, diff, _ = nl.get_neighbors(i)
        pair_i_idx += [i] * len(indices)               # local index of pair i
        pair_j_idx.append(indices)   # local index of pair j
        n_diff.append(diff)
    # 计算代码执行所需的时间
    elapsed_time = time.time() - start_time
    print(f"asap3代码执行耗时: {elapsed_time} 秒")
    break

# %%
cube_size = (20, 20, 20)
num_particles = 800
cut_off_radius = 2.0
lc_cube = grid_based_nbl_gpu.GB_NBL_CUDA_cube(cube_size, num_particles, cut_off_radius)

# 生成测试数据
inputs = xyz

# 记录开始时间
start_time = time.time()

# 调用 constructor 方法进行初始化
lc_cube.constructor(inputs)

# 计算代码执行所需的时间
elapsed_time = time.time() - start_time
print(f"grid based nbl CUDA 代码执行耗时: {elapsed_time} 秒")


# # 测试 get_neighbors
# particle_seq = 2
# neighbors = lc_cube.get_neighbors(particle_seq)
# print(f"Neighbors of particle {particle_seq}: {neighbors}")

# # 逐个颗粒进行检测
# for i in range(len(atoms)):
#     idx1, diff, _ = nl.get_neighbors(i)
#     idx2 = lc_cube.get_neighbors(i)

#     if len(idx1) != len(idx2):
#         print(i)
#         print(f"idx1: {idx1}, idx2: {idx2}")
# %%


# %%
