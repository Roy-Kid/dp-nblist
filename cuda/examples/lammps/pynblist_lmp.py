import cuda_nblist_py as nblist
import time

def get_nblist(length, xyz, cutoff, skin=0.0):
    """Get neighbor list from xyz coordinates and cutoff radius.

    Args:
        xyz (np.ndarray): Coordinates of atoms.
        cutoff (float): Cutoff radius.

    Returns:
        np.ndarray: Neighbor list.
    """
    # 记录开始时间
    start_time = time.time()
    # 调用 constructor 方法进行初始化
    box = nblist.Box(length)
    nb = nblist.createNeighborList("celllist",box, cutoff, skin)
    nb.build(xyz)
    # 计算代码执行所需的时间
    elapsed_time = time.time() - start_time
    print(f"build耗时: {elapsed_time} 秒")
    start_time = time.time()
    nb.update(xyz)
    elapsed_time = time.time() - start_time
    print(f"update耗时: {elapsed_time} 秒")
    nb.out()
    return nb.gettime()

def main(file, cutoff, skin=0.0):
    temp_time_list = []
    for f in file:
        lines = open(f,encoding='utf-8').readlines()
        #读取所有的原子坐标，找出盒子最大值
        xyz = []
        length = [0,0,0]
        # for i in range(1,len(lines)-1):#len(lines)-1
        natoms = int(lines[2].split()[0])
        i = 0
        spaceline = 0
        while True:
            i += 1
            linesp = lines[i].split()
            if len(linesp) == 0:
                spaceline += 1
                if spaceline == 2:
                    for j in range(3):
                        i += 1
                        length[j] = float(lines[i].split()[1]) - float(lines[i].split()[0])
                elif spaceline == 6:
                    for _ in range(natoms):
                        i += 1
                        linesp = lines[i].split()
                        xyz.append([float(linesp[2]),float(linesp[3]),float(linesp[4])])
                    break
        print(f"natoms: {natoms}")
        temp_time = get_nblist(length, xyz, cutoff, skin)

if __name__ == "__main__":
    time_list = []
    flie_norm = ["heterogeneous.lmp"]
    main(flie_norm, 3.0)
