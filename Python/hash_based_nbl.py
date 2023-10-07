from hash_function import hash_function_float
import numpy as np
from collections import defaultdict
import time

class HB_NBL_cube:
    def __init__(self, cube_size, num_particles, cut_off_radius):
        """
        初始化哈希表邻接格子对象
        参数:
            - cube_size: 立方体的大小 (一个包含3个值的 NumPy array 或类似的列表)
            - num_particles: 颗粒的数量 (整数)
            - cut_off_radius: 截断半径，用于确定邻居关系 (浮点数)
        """
        self.cube_size = np.array(cube_size)  # 立方体的大小

        self.rc = cut_off_radius  # 截断半径

        self.num_particles = num_particles

        # 表示颗粒的邻接关系
        self.particle_list = [[] for _ in range(num_particles)]

        # 哈希表，使用defaultdict构建
        self.hash_table = defaultdict(list)

        self.coefficient = 15

        self.range_size = 200

    def _cal_hash_code(self, inputs):
        """
        计算颗粒的哈希值并构建哈希表
        参数:
            - inputs: 颗粒的坐标数据，一个大小为(num_particles, 3)的 NumPy array
        """
        def add_to_hash_table(hash_range, index):
            if not self.hash_table[hash_range]:
                self.hash_table[hash_range] = []
            self.hash_table[hash_range].append(index)

        # 构建邻接关系字典

        # 非周期边界
        for i in range(self.num_particles):
            hash_code = hash_function_float(inputs[i], self.coefficient)

            hash_range = hash_code // self.range_size
            if not self.hash_table[hash_range]:
                self.hash_table[hash_range] = []
            self.hash_table[hash_range].append(i)

        # # # 周期边界
        # for i in range(self.num_particles):
        #     hash_code = hash_function_float(inputs[i], self.coefficient)
        #     hash_range = hash_code // self.range_size

        #     if any(np.less(inputs[i], self.rc)) or any(np.greater(inputs[i], self.cube_size - self.rc)):
        #         # 颗粒靠近边界，需要构建镜像
        #         add_to_hash_table(hash_range, i)

        #         # 构建并添加周期边界镜像
        #         for axis in range(3):
        #             if inputs[i, axis] < self.rc:
        #                 mirrored_position = inputs[i].copy()
        #                 mirrored_position[axis] += self.cube_size[axis]
        #                 hash_code = hash_function_float(mirrored_position, self.coefficient)
        #                 hash_range = hash_code // self.range_size
        #                 add_to_hash_table(hash_range, i)

        #             if inputs[i, axis] > self.cube_size[axis] - self.rc:
        #                 mirrored_position = inputs[i].copy()
        #                 mirrored_position[axis] -= self.cube_size[axis]
        #                 hash_code = hash_function_float(mirrored_position, self.coefficient)
        #                 hash_range = hash_code // self.range_size
        #                 add_to_hash_table(hash_range, i)
        #     else:
        #         # 颗粒不靠近边界，直接添加到哈希表
        #         add_to_hash_table(hash_range, i)

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

        start_time = time.time()

        # 创建hash_code表
        self._cal_hash_code(inputs)

        elapsed_time = time.time() - start_time
        print(f"代码执行耗时: {elapsed_time} 秒")

        start_time = time.time()
        # 建立链接关系, 考虑截断半径
        for particle_seq, particle_xyz in enumerate(inputs):
            particle_hash_value = hash_function_float(particle_xyz, self.coefficient)
            particle_hash_range = particle_hash_value // self.range_size

            neighbor_particles = set()

            neighbor_particles.update(self.hash_table[particle_hash_range])
            neighbor_particles.update(self.hash_table[particle_hash_range + 1])
            neighbor_particles.update(self.hash_table[particle_hash_range - 1])

            # print(len(neighbor_particles))

            for neighbor_particle_seq, neighbor_particle_xyz in enumerate(neighbor_particles):
                if neighbor_particle_seq <= particle_seq:
                    continue

                distance = np.linalg.norm(self._get_min_diff(particle_xyz, neighbor_particle_xyz))

                if distance - self.rc <= 1e-10:
                    self.particle_list[particle_seq].append(neighbor_particle_seq)
                    self.particle_list[neighbor_particle_seq].append(particle_seq)
        
        elapsed_time = time.time() - start_time
        print(f"代码执行耗时: {elapsed_time} 秒")

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

