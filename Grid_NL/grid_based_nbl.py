#%%
import numpy as np
import time 

class GB_NBL_cube:
    def __init__(self, cube_size, num_particles, cut_off_radius):
        self.cube_size = np.array(cube_size) # cube的大小

        # cube离散化后的cell数量
        self.lc = np.array(cube_size) / np.array(cut_off_radius) 
        self.rc = cut_off_radius # cut off radius
        self.num_particles = num_particles # 空间中颗粒的数量

        # grid_size = cube_size / (cube_size / rc)
        self.grid_size = np.array(cube_size) / (self.lc)
        
        # 计算总cell数量，此处不考虑ghost_cell
        self.num_cells = int(self.lc[0] * self.lc[1] * self.lc[2])

        # 表示颗粒的邻接关系
        self.particle_list = []   

        # 表示cell间的邻接关系
        self.cell_list = []  


    def _xyz2ind(self, xyz):
        """
        将坐标转换为对应的 cell 索引
        参数:
            - xyz: 三维坐标 (NumPy array 或单个值)
        返回:
            - cell_index: 对应的 cell 索引 (NumPy array 或单个值)
        """
        xyz = np.atleast_2d(xyz)

        # 将坐标限制在周期性边界内
        xyz = np.remainder(xyz, self.cube_size)

        cell_index = np.floor_divide(xyz, self.grid_size)
        cell_index = cell_index[:, 0] * self.lc[1] * self.lc[2] + cell_index[:, 1] * self.lc[2] + cell_index[:, 2]
        cell_index = cell_index.astype(int)
        if len(cell_index) == 1:
            return cell_index[0]
        return cell_index


    def _ind2xyz(self, ind):
        """
        将 cell 索引转换为对应的中心坐标
        参数:
            - ind: cell 索引 (NumPy array 或单个值)
        返回:
            - xyz: 对应的坐标 (NumPy array 或单个值)
        """
        ind = np.atleast_1d(ind)
        x = ind // (self.lc[1] * self.lc[2]) * self.grid_size[0] + self.grid_size[0] / 2
        y = (ind // self.lc[2]) * self.grid_size[1] + self.grid_size[1] / 2
        z = ind % self.lc[2] * self.grid_size[2] + self.grid_size[2] / 2
        xyz = np.column_stack((x % self.cube_size[0], y % self.cube_size[1], z % self.cube_size[2]))
        if len(xyz) == 1:
            return xyz[0]
        return xyz


    def _get_min_diff(self, xyz1, xyz2):
        """
        计算考虑周期边界的颗粒间最小距离
        参数:
            - xyz1: 第一个颗粒的坐标 (NumPy array 或单个值)
            - xyz2: 第二个颗粒的坐标 (NumPy array 或单个值)
        返回:
            - distance: 考虑周期边界的最小镜像距离 (NumPy array 或单个值)
        """
        difference = xyz2 - xyz1
        difference = np.remainder(difference + self.cube_size / 2, self.cube_size) - self.cube_size / 2
        return difference


    def constructor(self, inputs):
        # 初始化 self.particle_list， self.cell_list
        self.particle_list = [[] for _ in range(self.num_particles)]
        self.cell_list = [[] for _ in range(self.num_cells)]

        # 把所有粒子都添加到对应的 cell 中
        particle_inds = self._xyz2ind(inputs)
        for particle_seq, particle_ind in enumerate(particle_inds):
            self.cell_list[particle_ind].append(particle_seq)

        # 建立链接关系, 考虑rc
        for particle_seq, xyz in enumerate(inputs):
            particle_ind = particle_inds[particle_seq]
            adjacent_cells = self.get_neighbor_cells(particle_ind)

            for cell_temp in adjacent_cells:
                neighbor_particles_temp = self.cell_list[cell_temp]
                if neighbor_particles_temp:
                    for neighbor_particle in neighbor_particles_temp:
                        if neighbor_particle == particle_seq:
                            continue
                        neighbor_xyz = inputs[neighbor_particle]
                        distance = np.linalg.norm(self._get_min_diff(xyz, neighbor_xyz))
                        # print(f"particle seq: {particle_seq}, particle neighbor: {neighbor_particle}, diff: {self._get_min_diff(xyz, neighbor_xyz)}, distance: {distance}")
                        if distance - self.rc <= 1e-10:
                            self.particle_list[particle_seq].append(neighbor_particle)
                              
            
    def get_neighbor_cells(self, cell_ind):
        """
        给定一个 cell 索引，返回其 26 个相邻的 cell 索引，考虑周期性边界条件
        参数:
            - cell_ind: 给定的 cell 索引
            - lc: cell 的离散化数量 (NumPy array)
        返回:
            - adjacent_cells: 26 个相邻的 cell 索引列表 及 自己本身
        """
        i, j, k = self._ind2xyz(np.array(cell_ind))
        adjacent_cells = []
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                for dk in [-1, 0, 1]:
                    ni = (i + di * self.grid_size[0])
                    nj = (j + dj * self.grid_size[1])
                    nk = (k + dk * self.grid_size[2])

                    xyz_temp = np.column_stack((ni, nj, nk))
                    ind_temp = self._xyz2ind(xyz_temp)

                    adjacent_cells.append(ind_temp)
        return set(adjacent_cells)


    def get_neighbors(self, particle_seq):
        """
        获取给定颗粒序号的邻近颗粒列表
        参数:
            - particle_seq: 给定的颗粒序号
        返回:
            - neighbors: 邻近颗粒列表
        """
        return self.particle_list[particle_seq]

    # 为避免多次建表，只在颗粒变化时修改linked cell信息
    # def update_positions(self):
    #     # 更新粒子的位置和速度的方法
    #     # 遍历所有粒子，根据计算的力更新其位置和速度

    # def add_particle(self, particle):
    #     # 添加粒子的方法

    # def remove_particle(self, particle):
    #     # 移除粒子的方法



#%%

if __name__ == '__main__':
    # 简单测试
    
    ############################################
    # 测试坐标 / cell index转换
    # 创建 LinkedCell_cube 对象
    cube_size = (10, 10, 10)
    num_particles = 6
    cut_off_radius = 1.0
    lc_cube = GB_NBL_cube(cube_size, num_particles, cut_off_radius)

    # 生成测试数据
    xyz = np.array([[1.5, 2.2, 3.7], [4.1, 5.9, 6.3], [7.8, 8.4, 9.2], 
                    [1.5, 2.2, 3.7], [4.1, 5.9, 6.3], [7.8, 8.4, 9.2]])

    # 测试 xyz2ind 方法
    cell_indices = lc_cube._xyz2ind(xyz)
    print("Cell indices:")
    print(cell_indices)

    # 测试 ind2xyz 方法
    xyz_restored = lc_cube._ind2xyz(cell_indices)
    print("Restored XYZ:")
    print(xyz_restored)

    cell_indices_r = lc_cube._xyz2ind(xyz_restored)
    print("Cell index from restored XYZ:")
    print(cell_indices_r)

    print(lc_cube.get_neighbor_cells(0))

    lc_cube.constructor(xyz)
    for i in range(6):
        print(lc_cube.get_neighbors(i))

    # print(lc_cube._get_min_diff(xyz[1], xyz[0]))
    ################################################
    # 测试constructor和get_neighbor
    # 创建 LinkedCell_cube 对象
    cube_size = (10, 10, 10)
    num_particles = 800
    cut_off_radius = 1.0
    lc_cube = GB_NBL_cube(cube_size, num_particles, cut_off_radius)
    print(lc_cube.get_neighbor_cells(1))

    # # 生成测试数据
    # np.random.seed(1)
    # inputs = np.random.random((num_particles, 3)) * cube_size

    # # 记录开始时间
    # start_time = time.time()

    # # 调用 constructor 方法进行初始化
    # lc_cube.constructor(inputs)

    # # 计算代码执行所需的时间
    # elapsed_time = time.time() - start_time
    # print(f"代码执行耗时: {elapsed_time} 秒")

    # # 测试 constructor
    # print("Particle List:")
    # for i, particle_indices in enumerate(lc_cube.particle_list):
    #     print(f"Particle {i}: {particle_indices}")

    # print("Cell List:")
    # for i, cell_index in enumerate(lc_cube.cell_list):
    #     print(f"Cell {i}: Cell {cell_index}")

    # # 测试 get_neighbors()
    # particle_seq = 0
    # neighbors = lc_cube.get_neighbors(particle_seq)
    # print(f"Neighbors of particle {particle_seq}: {neighbors}")

    # # 测试get_neighbor_cells()
    # test_cell_ind = 0  # 测试用的 cell 索引
    # neighbor_cells = lc_cube.get_neighbor_cells(test_cell_ind)

    # # 打印结果
    # print("相邻单元的索引列表:")
    # print(neighbor_cells)

# %%
