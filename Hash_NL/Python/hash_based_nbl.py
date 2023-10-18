from hash_function import hash_function_float_2 as hash_function_float
import numpy as np
from collections import defaultdict
import time

class HB_NBL_cube:
    def __init__(self, cube_size, num_particles, cut_off_radius, coefficient=1, range_size=350, move_size=1):
        """
        初始化哈希表邻接格子对象

        参数:
            - cube_size (list or numpy.ndarray): 立方体的大小，包含三个值，用于定义模拟的三维空间的尺寸。

            - num_particles (int): 颗粒的数量，定义系统中存在的颗粒的总数。

            - cut_off_radius (float): 截断半径，用于确定颗粒之间的邻居关系的截断半径。
              颗粒之间的距离小于或等于截断半径的颗粒被认为是邻居。

            - coefficient (float, optional): 一个浮点数，是对于浮点数的放缩系数，用于使Z-order函数可以操作浮点数。
              默认为1，对于较小的系统，可能需要适当放大它。

            - range_size (int, optional): 一个整数，是建立哈希函数的桶大小。
              在几何空间中，它对应着一个自适应区域，与颗粒密度相关。
              桶越小遍历的粒子越少，效率也越高，反之亦然。
              默认为350。

            - move_size (int, optional): 一个整数，是颗粒映射时的边界范围。
              边界范围设置的越大，精确度越高，但随之带来更多的计算量，反之亦然。
              默认为1.
        """
        self.cube_size = np.array(cube_size)  # 立方体的大小

        self.rc = cut_off_radius  # 截断半径

        self.num_particles = num_particles

        # 表示颗粒的邻接关系
        self.particle_list = [[] for _ in range(num_particles)]

        # 哈希表: hashed_value : particle_seq，使用defaultdict构建
        self.hash_table = defaultdict(list)
        
        self.hash_list = [[] for _ in range(num_particles)]

        self.range_size = 50# 50 - rc = 2


    def _cal_hash_code(self, inputs):
        """
        计算颗粒的哈希值并构建哈希表
        参数:
            - inputs: 颗粒的坐标数据，一个大小为(num_particles, 3)的 NumPy array
        """
        def add_to_hash_table(input_temp, index):
            hash_code = hash_function_float(input_temp)
            hash_range_temp = hash_code // self.range_size * self.range_size

            if not self.hash_table[hash_range_temp]:
                self.hash_table[hash_range_temp] = []
            self.hash_table[hash_range_temp].append(index)   

            return hash_range_temp  # 返回集合



        # 构建邻接关系字典

        # # 非周期边界
        # for i in range(self.num_particles):
        #     hash_code = hash_function_float(inputs[i], self.coefficient)
        #     print(inputs[i])
        #     hash_range = hash_code // self.range_size
        #     print(hash_code)
        #     print(hash_range)
        #     if not self.hash_table[hash_range * self.range_size]:
        #         self.hash_table[hash_range * self.range_size] = []
        #     self.hash_table[hash_range * self.range_size].append(i)
        #     print(hash_range * self.range_size)
        #     print("___________")

        # def add_to_hash_table(input_temp, index):
        #     x, y, z = input_temp
        #     input_list = [[x + self.move_size, y, z],
        #                 [x - self.move_size, y, z],
        #                 [x, y + self.move_size, z],
        #                 [x, y - self.move_size, z],
        #                 [x, y, z + self.move_size],
        #                 [x, y, z - self.move_size]]
        #     # input_list = [[x, y, z]]
        #     hash_range_set = set()  # 使用集合来存储唯一的哈希范围值

        #     for value in input_list:
        #         hash_code = hash_function_float(value)
        #         hash_range_temp = hash_code // self.range_size
        #         hash_range_set.add(hash_range_temp * self.range_size)  # 使用add()方法来添加唯一值

        #     for hash_value in hash_range_set:  
        #         if not self.hash_table[hash_value]:
        #             self.hash_table[hash_value] = []
        #         self.hash_table[hash_value].append(index)   
        #     return hash_range_set  # 返回集合


        # # 周期边界
        for i in range(self.num_particles):
            hash_list_temp = []
            hash_temp = add_to_hash_table(inputs[i], i)
            # print(hash_temp)
            hash_list_temp.append(hash_temp)
            
            if any(np.less(inputs[i], self.rc)) or any(np.greater(inputs[i], self.cube_size - self.rc)):
                # 颗粒靠近边界，需要构建镜像

                # 构建并添加周期边界镜像
                for axis in range(3):
                    if inputs[i, axis] < self.rc:
                        mirrored_position = inputs[i].copy()
                        mirrored_position[axis] += self.cube_size[axis]
                        hash_temp = add_to_hash_table(mirrored_position, i)
                        hash_list_temp.append(hash_temp)

                        # hash_code = hash_function_float(mirrored_position)
                        # hash_range = hash_code // self.range_size
                        # add_to_hash_table(hash_range * self.range_size, i)
                        # hash_list_temp.append(hash_range * self.range_size)

                    if inputs[i, axis] > self.cube_size[axis] - self.rc:
                        mirrored_position = inputs[i].copy()
                        mirrored_position[axis] += self.rc
                        mirrored_position[axis] -= self.cube_size[axis]
                        hash_temp = add_to_hash_table(mirrored_position, i)
                        hash_list_temp.append(hash_temp)

            # print(hash_list_temp)
            self.hash_list[i] = set(hash_list_temp)
            # print(hash_list_temp)




    # # 计算所有hash values
    # def _cal_hashed_values(self, input_v):
    #     res = []
    #     # # 周期边界
    #     hash_code = hash_function_float(input_v, self.coefficient)
    #     hash_range = hash_code // self.range_size

    #     res.append(hash_range * self.range_size)

    #     if any(np.less(input_v, self.rc)) or any(np.greater(input_v, self.cube_size - self.rc)):
    #         # 颗粒靠近边界，需要构建镜像

    #         # 构建并添加周期边界镜像
    #         for axis in range(3):
    #             if input_v[axis] < self.rc:
    #                 mirrored_position = input_v.copy()
    #                 mirrored_position[axis] += self.cube_size[axis]
    #                 hash_code = hash_function_float(mirrored_position, self.coefficient)
    #                 hash_range = hash_code // self.range_size
    #                 res.append(hash_range * self.range_size)
    #                 print(mirrored_position)

    #             if input_v[axis] > self.cube_size[axis] - self.rc:
    #                 mirrored_position = input_v.copy()
    #                 # mirrored_position[axis] -= self.cube_size[axis]
    #                 mirrored_position[axis] += self.rc
    #                 mirrored_position[axis] -= self.cube_size[axis]

    #                 hash_code = hash_function_float(mirrored_position, self.coefficient)
    #                 hash_range = hash_code // self.range_size
    #                 res.append(hash_range * self.range_size)
    #                 print(mirrored_position)

    #     return res


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
        """
        构建哈希表和邻居列表，考虑截断半径
        参数:
            - inputs: 颗粒的坐标数据，一个大小为(num_particles, 3)的 NumPy array
        """

        # start_time = time.time()

        # 创建hash_code表
        self._cal_hash_code(inputs)

        # elapsed_time = time.time() - start_time
        # print(f"创建hash_table代码执行耗时: {elapsed_time} 秒")

        # start_time = time.time()
        # 建立链接关系, 考虑截断半径
        for particle_seq, particle_xyz in enumerate(inputs):
            particle_hash_values = self.hash_list[particle_seq]

            for particle_hash_value in particle_hash_values:
                # # particle_hash_range = particle_hash_value // self.range_size
                neighbor_particles = set()
                neighbor_particles.update(self.hash_table[particle_hash_value])
                neighbor_particles.update(self.hash_table[particle_hash_value + self.range_size])
                neighbor_particles.update(self.hash_table[particle_hash_value - self.range_size])

                # neighbor_particles.update(self.hash_table[particle_hash_range * self.range_size])
                # neighbor_particles.update(self.hash_table[(particle_hash_range + 1) * self.range_size])
                # neighbor_particles.update(self.hash_table[(particle_hash_range - 1) * self.range_size])
                # neighbor_particles.update(self.hash_table[(particle_hash_range + 2) * self.range_size])
                # neighbor_particles.update(self.hash_table[(particle_hash_range - 2) * self.range_size])

                # neighbor_particles = self.hash_table[particle_hash_value]


                for neighbor_particle_seq in neighbor_particles:
                    if neighbor_particle_seq <= particle_seq:
                        continue
                    
                    if neighbor_particle_seq in self.particle_list[particle_seq]:
                        continue

                    distance = np.linalg.norm(self._get_min_diff(particle_xyz, inputs[neighbor_particle_seq]))
                    # print(f"xyz: {particle_xyz}, neighbor_xyz: {inputs[neighbor_particle_seq]}")
                    # print(distance)
                    # if distance - self.rc <= 1e-10:
                    if distance <= self.rc:
                        self.particle_list[particle_seq].append(neighbor_particle_seq)
                        self.particle_list[neighbor_particle_seq].append(particle_seq)

                    # # 优化距离计算，不计算平方根
                    # squared_distance = np.sum((self._get_min_diff(particle_xyz, inputs[neighbor_particle_seq]))**2)
                    # if squared_distance <= self.rc**2:
                    #     self.particle_list[particle_seq].append(neighbor_particle_seq)
                    #     self.particle_list[neighbor_particle_seq].append(particle_seq)


        # elapsed_time = time.time() - start_time
        # print(f"建立临近表代码执行耗时: {elapsed_time} 秒")

    def get_neighbors(self, particle_seq):
        """
        获取给定颗粒序号的邻近颗粒列表
        参数:
            - particle_seq: 给定的颗粒序号
        返回:
            - neighbors: 邻近颗粒列表
        """
        return self.particle_list[particle_seq]


if __name__ == '__main__':
    # 简单测试
    
       ################################################
    # 测试constructor和get_neighbor
    # 创建 LinkedCell_cube 对象
    cube_size = (10, 10, 10)
    num_particles = 50000
    cut_off_radius = 1.0
    lc_cube = HB_NBL_cube(cube_size, num_particles, cut_off_radius)

    # 生成测试数据
    np.random.seed(1)
    inputs = np.random.random((num_particles, 3)) * cube_size
    print(inputs)

    # 记录开始时间
    start_time = time.time()

    # 调用 constructor 方法进行初始化
    lc_cube.constructor(inputs)

    # 计算代码执行所需的时间
    elapsed_time = time.time() - start_time
    print(f"代码执行耗时: {elapsed_time} 秒")

    # # 测试 constructor
    # print("Particle List:")
    # for i, particle_indices in enumerate(lc_cube.particle_list):
    #     print(f"Particle {i}: {particle_indices}")

    # # 测试 get_neighbors()
    # particle_seq = 0
    # neighbors = lc_cube.get_neighbors(particle_seq)
    # print(f"Neighbors of particle {particle_seq}: {neighbors}")

    # for key, value in lc_cube.hash_table.items():
    #     print(key)
    #     print(value)

