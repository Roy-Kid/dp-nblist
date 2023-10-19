from matplotlib.dviread import Box
import numpy as np
from collections import defaultdict

import numpy as np
from python.dp-nblist.baseNBL import BaseNBL

def hash_function_int(point):
    x, y, z = point
    # 将x、y、z分别转换为32位整数
    x = (x | (x << 16)) & 0x030000FF
    x = (x | (x << 8)) & 0x0300F00F
    x = (x | (x << 4)) & 0x030C30C3
    x = (x | (x << 2)) & 0x09249249

    y = (y | (y << 16)) & 0x030000FF
    y = (y | (y << 8)) & 0x0300F00F
    y = (y | (y << 4)) & 0x030C30C3
    y = (y | (y << 2)) & 0x09249249

    z = (z | (z << 16)) & 0x030000FF
    z = (z | (z << 8)) & 0x0300F00F
    z = (z | (z << 4)) & 0x030C30C3
    z = (z | (z << 2)) & 0x09249249

    # 将x、y、z按位交错合并为一个整数
    return x | (y << 1) | (z << 2)


def hash_function_float_1(point, coefficient):
    x, y, z = point
    # 将浮点数坐标乘以放大倍数，然后取整，得到离散化的整数坐标
    discrete_x = int(x * coefficient)
    discrete_y = int(y * coefficient)
    discrete_z = int(z * coefficient)

    # 对离散化的整数坐标进行位运算，生成哈希值
    x = (discrete_x | (discrete_x << 16)) & 0x030000FF
    x = (x | (x << 8)) & 0x0300F00F
    x = (x | (x << 4)) & 0x030C30C3
    x = (x | (x << 2)) & 0x09249249

    y = (discrete_y | (discrete_y << 16)) & 0x030000FF
    y = (y | (y << 8)) & 0x0300F00F
    y = (y | (y << 4)) & 0x030C30C3
    y = (y | (y << 2)) & 0x09249249

    z = (discrete_z | (discrete_z << 16)) & 0x030000FF
    z = (z | (z << 8)) & 0x0300F00F
    z = (z | (z << 4)) & 0x030C30C3
    z = (z | (z << 2)) & 0x09249249

    # 将离散化后的整数坐标按位交错合并为一个整数，作为最终的哈希值
    hash_value = x | (y << 1) | (z << 2)
    return hash_value

    # # 将离散化后的整数坐标按位交错合并为一个整数，作为最终的哈希值
    # hash_value = x | (y << 1) | (z << 2)

    # hash_bits = int(math.log10(hash_value) + 1) if hash_value != 0 else 1

    # return hash_bits


def hash_function_float_2(point):
    x, y, z = point
    # 将浮点坐标映射到整数坐标范围（例如，0到2^16 - 1）
    max_coord = 1  # 2^16 - 1
    x_int = int(x * max_coord)
    y_int = int(y * max_coord)
    z_int = int(z * max_coord)

    # 将整数坐标合并成一个一维哈希值，使用 Z-order curve 映射
    combined = 0
    for i in range(16):  # 假设每个坐标分量都是 16 位
        combined |= (x_int & 1) << (3 * i)
        combined |= (y_int & 1) << (3 * i + 1)
        combined |= (z_int & 1) << (3 * i + 2)
        x_int >>= 1
        y_int >>= 1
        z_int >>= 1

    return combined


def custom_hash(point, max_bits=32):
    x, y, z = point
    # Ensure x, y, and z are non-negative integers.
    x = int(x)
    y = int(y)
    z = int(z)

    # Calculate the maximum coordinate value based on the number of bits.
    max_coord = 2 ** max_bits - 1

    # Ensure the coordinates are within bounds.
    if x < 0 or x > max_coord or y < 0 or y > max_coord or z < 0 or z > max_coord:
        raise ValueError("Coordinates are out of bounds for the specified number of bits.")

    # Interleave the bits of x, y, and z to create a single integer.
    result = 0
    for i in range(max_bits):
        result |= ((x >> i) & 1) << (3 * i + 2)
        result |= ((y >> i) & 1) << (3 * i + 1)
        result |= ((z >> i) & 1) << (3 * i)

    return result


def hash_zhao(point):
    c1, c2, c3 = 1, 2**3, 4**3
    x, y, z = point

    return c1 * x + c2 * y + c3 * z



class HashNBL(BaseNBL):
    def __init__(self, box:Box, cut_off_radius, coefficient=1, range_size=350, move_size=1):
        """
        Initialize the hash table neighbor cube object.

        Args:
            - cube_size (list or numpy.ndarray): The size of the cube, containing three values, used to define the size of the three-dimensional space in the simulation.

            - num_particles (int): The number of particles, defining the total number of particles in the system.

            - cut_off_radius (float): The cut-off radius used to determine neighbor relationships between particles.
              Particles with a distance less than or equal to the cut-off radius are considered neighbors.

            - coefficient (float, optional): A floating-point scaling factor for handling floating-point numbers, used to make the Z-order function work with floating-point numbers.
              Defaults to 1, and may need to be increased for smaller systems.

            - range_size (int, optional): An integer representing the bucket size for constructing the hash function.
              In geometric space, it corresponds to an adaptive region related to particle density.
              A smaller bucket size results in traversing fewer particles, improving efficiency, and vice versa.
              Defaults to 350.

            - move_size (int, optional): An integer representing the boundary range for particle mapping.
              A larger boundary range provides higher accuracy but comes with more computational cost, and vice versa.
              Defaults to 1.
        """
        self.cube_size = np.array(cube_size)  # Size of the cube

        self.rc = cut_off_radius  # Cut-off radius

        self.num_particles = num_particles

        # Represents the adjacency relationships of particles
        self.particle_list = [[] for _ in range(num_particles)]

        # Hash table: hashed_value : particle_seq, constructed using defaultdict
        self.hash_table = defaultdict(list)

        self.hash_list = [[] for _ in range(num_particles)]

        self.range_size = 50  # 50 - rc = 2

    def _cal_hash_code(self, inputs):
        """
        Calculate the hash code of particles and construct the hash table.

        Args:
            - inputs: Particle coordinate data, a NumPy array of size (num_particles, 3).
        """
        def add_to_hash_table(input_temp, index):
            hash_code = hash_function_float(input_temp)
            hash_range_temp = hash_code // self.range_size * self.range_size

            if not self.hash_table[hash_range_temp]:
                self.hash_table[hash_range_temp] = []
            self.hash_table[hash_range_temp].append(index)

            return hash_range_temp  # Return the set

        # Periodic boundaries
        for i in range(self.num_particles):
            hash_list_temp = []
            hash_temp = add_to_hash_table(inputs[i], i)
            hash_list_temp.append(hash_temp)

            if any(np.less(inputs[i], self.rc)) or any(np.greater(inputs[i], self.cube_size - self.rc)):
                # Particles close to the boundaries need to build mirrors

                # Build and add periodic boundary mirrors
                for axis in range(3):
                    if inputs[i, axis] < self.rc:
                        mirrored_position = inputs[i].copy()
                        mirrored_position[axis] += self.cube_size[axis]
                        hash_temp = add_to_hash_table(mirrored_position, i)
                        hash_list_temp.append(hash_temp)

                    if inputs[i, axis] > self.cube_size[axis] - self.rc:
                        mirrored_position = inputs[i].copy()
                        mirrored_position[axis] += self.rc
                        mirrored_position[axis] -= self.cube_size[axis]
                        hash_temp = add_to_hash_table(mirrored_position, i)
                        hash_list_temp.append(hash_temp)

            self.hash_list[i] = set(hash_list_temp)

    def _get_min_diff(self, xyz1, xyz2):
        """
        Calculate the minimum distance between particles, considering periodic boundaries.

        Args:
            - xyz1: The coordinates of the first particle (NumPy array or a single value).
            - xyz2: The coordinates of the second particle (NumPy array or a single value).

        Returns:
            - distance: The minimum mirrored distance considering periodic boundaries (NumPy array or a single value).
        """
        difference = xyz2 - xyz1
        difference = np.remainder(difference + self.cube_size / 2, self.cube_size) - self.cube_size / 2
        return difference

    def constructor(self, inputs):
        """
        Build the hash table and neighbor list, considering the cut-off radius.

        Args:
            - inputs: Particle coordinate data, a NumPy array of size (num_particles, 3).
        """
        # Create the hash code table
        self._cal_hash_code(inputs)

        # Establish connection relationships, considering the cut-off radius
        for particle_seq, particle_xyz in enumerate(inputs):
            particle_hash_values = self.hash_list[particle_seq]

            for particle_hash_value in particle_hash_values:
                neighbor_particles = set()
                neighbor_particles.update(self.hash_table[particle_hash_value])
                neighbor_particles.update(self.hash_table[particle_hash_value + self.range_size])
                neighbor_particles.update(self.hash_table[particle_hash_value - self.range_size])

                for neighbor_particle_seq in neighbor_particles:
                    if neighbor_particle_seq <= particle_seq:
                        continue

                    if neighbor_particle_seq in self.particle_list[particle_seq]:
                        continue

                    distance = np.linalg.norm(self._get_min_diff(particle_xyz, inputs[neighbor_particle_seq]))
                    if distance <= self.rc:
                        self.particle_list[particle_seq].append(neighbor_particle_seq)
                        self.particle_list[neighbor_particle_seq].append(particle_seq)

    def get_neighbors(self, particle_seq):
        """
        Get the list of nearby particles for the given particle sequence.

        Args:
            - particle_seq: The given particle sequence.

        Returns:
            - neighbors: The list of nearby particles.
        """
        return self.particle_list[particle_seq]

