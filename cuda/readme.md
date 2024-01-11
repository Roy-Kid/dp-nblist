## Introduction  
This repository features the implementation of a grid-based neighbor list algorithm leveraging CUDA.
***
**Directory Structure**:

**src**: Contains the source code.  
**external**: Hosts external modules.  
**include**: Encompasses header files.  
**example**: Offers examples demonstrating the usage of the Python interface.  
**tests**: Holds unit tests for both C++ and CUDA code.  
***
## How to Compile
1. Download the code.
2. Navigate to the external directory:
```
cd external
tar -zxvf doctest-2.4.11.tar.gz
tar -zxvf pybind11-2.11.1.tar.gz
```
3. Return to the root directory:
```
cd ..
```
4. Create a build directory:
```
mkdir build
cd build
```
5. Run CMake with the desired CUDA architecture (modify as needed):
```
cmake -DCUDA_ARCHITECTURE=61 ..
```
6. Build the project:
```
make
```
Upon completion, two directories will be generated:  
 - **lib**: Contains a Python module wrapped by Pybind11, named "cuda_nblist_py".  
 - **bin**: Houses the executable file for testing purposes.
