Here are several directories representing different types of models, including homogeneous, heterogeneous, interface, protein in water, and a series of scaling particles. The remaining files are used for invoking the LAMMPS dynamic library to generate neighbor lists. Below is an introduction on how to use them.
## Environment Setup
### 1. Compiling LAMMPS (Optional)
Extract LAMMPS, navigate to the LAMMPS directory, create a 'build-lib' directory, and refer to the [tutorial](https://lammpscn2.vercel.app/Tutorial/install/#step2a-%E4%BD%BF%E7%94%A8cmake%E7%BC%96%E8%AF%91lammps):
```
tar -zxvf lammps-stable-Aug2.tar.gz
cd lammps-stable-Aug2.tar.gz
mkdir build-so
cd build-so
```
Next, open `../cmake/presets/basic.cmake` with a text editor and add 'PYTHON' inside the 'set' section:
```
set(ALL_PACKAGES KSPACE MANYBODY MOLECULE RIGID PYTHON)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
```
Now, generate the makefile and compile as a dynamic library with additional flags:
```
cmake -C ../cmake/presets/basic.cmake -DLAMMPS_EXCEPTIONS=yes -DBUILD_LIB=yes -DBUILD_SHARED_LIBS=yes ../cmake
make -j8
``` 
### 2. Adding Environment Variables (Optional)
```
vim ~/.bashrc
export LIBRARY_PATH=/home/lammps/lammps-2Aug2023/build-so/:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/lammps/lammps-2Aug2023/build-so/:$LD_LIBRARY_PATH
source ~/.bashrc
```
### 3. Using the Dynamic Library to Retrieve Neighbor Lists
If you skip the first two steps, you can directly [download](https://neepueducn-my.sharepoint.com/:f:/g/personal/2016306040242_neepu_edu_cn/EvHuTvsjeo1GpHaDdD_KuHYBK03IKg9sXQao7LXXif6_ag?e=jqfL4H) the precompiled dynamic library. Place the library file in the same directory and compile as follows:
```
mpicxx -std=c++11 -o lmp_50000 lmp_neb.cpp -llammps -lmpi
```

## Explanation of Dynamic Library Invocation
First, create an 'lmp' object, then use the methods of the 'input' object to run the input file 'in.test'. The generated neighbor list can be obtained through the 'neighbor' object.
```
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
    int nlist = lmp->neighbor->nlist;// Number of neighbor lists generated, at least 1
    NeighList *list = lmp->neighbor->lists[0];
    int inum = list->inum; // Number of atoms
    std::cout << inum << std::endl;
    int *ilist = list->ilist; // Local atom indices
    int *numneigh = list->numneigh; // Number of neighbors
    int **firstneigh = list->firstneigh; // Local indices of neighboring atoms
    int *tags = lmp->atom->tag;// Atom index mapping
    
    // Open a file for writing
    std::ofstream outFile("output.txt");
    if (outFile.is_open()) {
        for (int ii = 0; ii < inum; ii++){
            int i = ilist[ii];
            outFile << i << '\t';// Index of the ith atom
            int jnum = numneigh[i];// Number of neighbors for the ith atom
            int *jlist = firstneigh[i];// List of neighbors for the ith atom
            // std::cout << "atom "<< i << " neighbor:" << std::endl;
            for (int jj = 0; jj < jnum; jj++){
                int j = jlist[jj];
                // std::cout << j << std::endl;
                outFile << tags[j] << '\t';// Atom index mapping
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
```
`int nlist = lmp->neighbor->nlist`: This return the number of generated neighbor lists, at least one.  
`NeighList *list = lmp->neighbor->lists[0]:`This return the first neighbor list.  
`int inum = list->inum`: This returns the number of atoms for which neighbor lists have been determined.  
`int *ilist = list->ilist`: This returns the local indices of atoms as 
`ilist[ii]`, where ii is the index of the same atom designated between 0 and 
(`inum-1`) according to its position in the neighbor list array.  
`int *numneigh = list->numneigh`: This returns the number of neighbors 
of each atom as `numneigh[i]`, where `i` is the local atom index.  
`int **firstneigh = list->firstneigh`: This returns the local index 
of the neighboring atom as `firstneigh[i][j]`, where `i` represents the central 
atom local index, and `j` represents the neighboring atom index. The `j` index is not 
the local atom index but is designated between 0 and (`numneigh[i]-1`), that 
is, according to its position in the neighbor list. The output value is the local atom 
index of the neighboring atom, `j`, of the central atom, `i`.  
```
for (int ii = 0; ii < inum; ii++){
    int i = ilist[ii];
    outFile << i << '\t';// Index of the ith atom
    int jnum = numneigh[i];// Number of neighbors for the ith atom
    int *jlist = firstneigh[i];// List of neighbors for the ith atom
    for (int jj = 0; jj < jnum; jj++){
        int j = jlist[jj];
        outFile << tags[j] << '\t';// Atom index mapping
    }
    outFile << '\n';
}
```
In the first loop over `inum` (with index `ii`), the local atom index (`i`) of each centrall atom is extracted (`i = ilist[ii]`), the number of neighbors for the atom (`jnum`) is retrieved (`jnum = numneigh[i]`), and an array of neighbors (`jlist`) is created (`jlist = firstneigh[i]`).  
A second loop (with index `jj`) located inside the first loop traverses over all neighbors `jnum` of atom `i`, and extracts the local atom index (`j`) of each neighbor (`j = jlist[jj]`). In effect, `jlist[jj]` returns the same output as `firstneigh[i][jj]`.

LAMMPS Input File:
```
# Simulation parameters
newton off
neighbor 0 bin
# Read the model
read_data 50000.lmp
# Set the potential function
pair_style zero 2.6 full
pair_coeff * *
```
`read_data 50000.lmp`means reading the model from the file`50000.lmp`.  
`pair_style zero 2.6 full`, where `2.6` represents the cutoff radius, specifies the pair style with a cutoff radius of 2.6.
