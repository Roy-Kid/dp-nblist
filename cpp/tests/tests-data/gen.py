import chemfiles
from ase import Atoms
from asap3 import FullNeighborList
import numpy as np

def read(filepath):
    # 使用chemfiles打开PDB文件
    trajectory = chemfiles.Trajectory(filepath)
    # 获取第一帧（如果有多帧的话）
    frame = trajectory.read()
    # 获取原子数量
    num_atoms = len(frame.atoms)
    # 获取原子坐标
    pos = frame.positions
    print(f'原子数量:{len(pos)}')
    atoms = Atoms(positions = pos)
    length = np.max(pos,axis=0) - np.min(pos,axis=0)
    print(f'盒子大小:{length}')
    atoms.set_cell(length)
    atoms.set_pbc([True, True, True])
    print(f'周期边界:{atoms.get_pbc()}')
    # print(atoms.get_cell())
    r_cutoff = 2.6
    nb = FullNeighborList(r_cutoff, atoms)
    with open('out_nbl.txt', 'w') as fw:
        for i in range(len(nb)):
            fw.write(f'{i}')
            for inb in nb[i]:
                fw.write(f'\t{inb}')
            fw.write('\n')

if __name__=="__main__":
    # 指定要读取的PDB文件路径
    pdb_file_path = '50000.pdb'
    read(pdb_file_path)
