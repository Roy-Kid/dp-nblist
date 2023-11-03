#include "lammps.h"
#include "input.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
using namespace LAMMPS_NS;
int main(int argc, char **argv)
{
    const char *lmpargv[] {"liblammps", "-log", "none"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);
    MPI_Init(&argc, &argv);
    
    // create lammps obiect
    LAMMPS *lmp = new LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
    lmp->input->file("in.test");
    lmp->input->one("run 0 post no");

    // get neblist
    int nlist = lmp->neighbor->nlist;//生成了几组邻居列表，至少1组
    NeighList *list = lmp->neighbor->lists[0];
    int inum = list->inum; //原子个数
    std::cout << inum << std::endl;
    int *ilist = list->ilist; //局部原子索引
    int *numneigh = list->numneigh; //邻居个数
    int **firstneigh = list->firstneigh; //邻居原子的局部索引
    int *tags = lmp->atom->tag;//原子index映射
    
    // Open a file for writing
    std::ofstream outFile("output_full.txt");
    if (outFile.is_open()) {
        for (int ii = 0; ii < inum; ii++){
            int i = ilist[ii];
            outFile << tags[i] << '\t';//第i个原子的索引
            int jnum = numneigh[i];//第i个原子的邻居个数
            int *jlist = firstneigh[i];//第i个原子的邻居列表
            // std::cout << "atom "<< i << " neighbor:" << std::endl;
            for (int jj = 0; jj < jnum; jj++){
                int j = jlist[jj];
                // std::cout << j << std::endl;
                outFile << tags[j] << '\t';//原子index映射
            }
            outFile << '\n';
        }
        outFile.close();
        std::cout << "File writing complete." << std::endl;
    } else {
        std::cerr << "Unable to open the file for writing." << std::endl;
    }
    delete lmp;
    return 0;
}
