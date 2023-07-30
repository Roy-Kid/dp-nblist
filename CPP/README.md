# Introduction
This repository contains two implementations of neighbor list algorithms based on grid-based designs. The first implementation is an unoptimized version (grid_based_nbl_cpp.cpp), and the second one has undergone a simple OpenMP optimization (grid_based_nbl_mp.cpp). Both implementations have been wrapped with Pybind, enabling seamless execution within Python.

# Files
The repository consists of the following files:

__50000.pdb__: The dataset for testing of 50000 moleculars. 
__grid_based_nbl_cpp.cpp__: Unoptimized neighborlist implementation based on grid-based design.
__grid_based_nbl_cpp.so__: The compiled library of Unoptimized neighborlist for python
__grid_based_nbl_mp.cpp__: Neighborlist implementation with basic OpenMP optimization for improved performance.
__grid_based_nbl_mp.so__: The compiled library of OpenMP neighborlistfor python
__validation.py__: A test file comparing the performance of the implemented neighbor list algorithms with ASAP3 on a dataset of size 50000.


# How to Compile
These CPP files already have Pybind implementations and can be compiled as follows:
```
g++ -o grid_based_nbl_cpp.so -shared -std=c++11 -fPIC -I/root/anaconda3/include/python3.8 grid_based_nbl_cpp.cpp __your_python_path__
``` 
Please replace __your_python_path__ with the appropriate path to your Python installation.

Once compiled, the resulting .so files can be used directly in Python, as illustrated in validation.py.

Please note that it's essential to ensure compatibility and verify the performance of the implementations on different devices to avoid potential technical issues.

# Reference
1. An introduction to grid-based neighbor list: https://aiichironakano.github.io/cs596/01-1LinkedListCell.pdf
2. Behley, Jens, Volker Steinhage, and Armin B. Cremers. "Efficient radius neighbor search in three-dimensional point clouds." 2015 IEEE international conference on robotics and automation (ICRA). IEEE, 2015.[Link](https://jbehley.github.io/papers/behley2015icra.pdf)
3. Howard, Michael P., et al. "Efficient neighbor list calculation for molecular simulation of colloidal systems using graphics processing units." Computer Physics Communications 203 (2016): 45-52. [Link](https://www.sciencedirect.com/science/article/pii/S0010465516300182)
4. Donev, Aleksandar, Salvatore Torquato, and Frank H. Stillinger. "Neighbor list collision-driven molecular dynamics simulation for nonspherical hard particles. I. Algorithmic details." Journal of computational physics 202.2 (2005): 737-764. [Link](https://www.sciencedirect.com/science/article/pii/S0021999104003146)

Please refer to the above reference for more details on the algorithm design and principles.

