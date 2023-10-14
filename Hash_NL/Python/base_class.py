class BaseHB_NBL:
    def __init__(self, cube_size, num_particles, cut_off_radius):
        self.cube_size = np.array(cube_size)
        self.rc = cut_off_radius
        self.num_particles = num_particles
        self.particle_list = [[] for _ in range(num_particles)]
        self.hash_table = defaultdict(list)
        self.coefficient = 15
        self.range_size = 200

    def _cal_hash_code(self, inputs):
        # 实现哈希值计算逻辑
        pass

    def _get_min_diff(self, xyz1, xyz2):
        # 实现最小距离计算逻辑
        pass

    def constructor(self, inputs):
        # 实现哈希表和邻居列表构建逻辑
        pass

    def get_neighbors(self, particle_seq):
        # 实现获取邻近颗粒的逻辑
        pass
