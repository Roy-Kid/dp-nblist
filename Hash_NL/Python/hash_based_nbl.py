from hash_function import hash_function_float_2 as hash_function_float
import numpy as np
from collections import defaultdict
import time

class HB_NBL_cube:
    def __init__(self, cube_size, num_particles, cut_off_radius, coefficient=1, range_size=350, move_size=1):
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



if __name__ == '__main__':  
       ################################################
    # Test functions: constructor and get_neighbor
    # Init LinkedCell_cube 
    cube_size = (10, 10, 10)
    num_particles = 50000
    cut_off_radius = 1.0
    lc_cube = HB_NBL_cube(cube_size, num_particles, cut_off_radius)

    # Generate testing dataset
    np.random.seed(1)
    inputs = np.random.random((num_particles, 3)) * cube_size
    print(inputs)

    # Start time 
    start_time = time.time()

    # Use the constructor function to do the initialization 
    lc_cube.constructor(inputs)

    # End time
    elapsed_time = time.time() - start_time
    print(f"代码执行耗时: {elapsed_time} 秒")

    # Test the constructor function
    # print("Particle List:")
    # for i, particle_indices in enumerate(lc_cube.particle_list):
    #     print(f"Particle {i}: {particle_indices}")

    # Test the get_neighbors function
    # particle_seq = 0
    # neighbors = lc_cube.get_neighbors(particle_seq)
    # print(f"Neighbors of particle {particle_seq}: {neighbors}")

    # for key, value in lc_cube.hash_table.items():
    #     print(key)
    #     print(value)

