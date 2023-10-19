import pytest

class TestHashFunc:

    pass

class TestHashNBL:

    def test_build(self):
        # we will use pytest-profiling to analysis code
        # for speed test of different size systems, using
        # pytert-benchmark and pytest.parametrize
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