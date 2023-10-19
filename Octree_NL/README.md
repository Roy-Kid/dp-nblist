# Introduction
This repository contains the implementation of a grid-based neighborlist algorithms based on CUDA with Thrust library. This method is tested with random number with C++.  
# Files
The repository consists of the following files:

__0.100.pdb__: The dataset for validation, consisting of 800 moleculars. 

__grid_based_nbl_gpu.cu__: The implementation of a grid-based neighborlist based on CUDA with Thrust library. 

__tools.cuh__: A CUDA head file for tools used in __grid_based_nbl_gpu.cu__.



# How to Compile
Those .cu files already have Pybind implementations. They can be compiled as following:
```
nvcc -arch=compute_89 -code=sm_89 your_file.cu -o executable
``` 
Please adjust the above code according to your specific hardware structure.

For the compilation of the .so file, use the following command:
```
nvcc -Xcompiler -fPIC -shared -std=c++11 -I/path/to/pybind11/include `python3 -m pybind11 --includes` -o grid_based_nbl_gpu.so grid_based_nbl.cu /root/anaconda3/envs/MD/lib/libpython3.8.so -Xcompiler -shared
```

The compiled .so files can be directly used in Python as demonstrated in __validation.py__. However, it's crucial to understand that the resulting .so file may not function as expected on various devices due to technical issues.

# Reference
1. Ferrando, Nestor, et al. "Octree-based, GPU implementation of a continuous cellular automaton for the simulation of complex, evolving surfaces." Computer Physics Communications 182.3 (2011): 628-640. [Link](https://pdf.sciencedirectassets.com/271575/1-s2.0-S0010465510X00145/1-s2.0-S0010465510004509/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPn%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQCDYj%2BLv57PuSkNQzs57zI2cfmHnHDDm84Sn2%2BapAWj2AIgRA1FA4iHUjAx%2BRhvGZW8UWHSsygjcvOf40cD6u80fQQqvAUIkv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDMa9LSlbVmcBFb4cPiqQBalnXD9%2B1T2Evd6JFSvE2WtGmT6DzcLdiQXcx%2FT58NjXSKLxNp7MKaUNX3kry0TUBn973mPi5OerDQ6w7i%2F1vX0Fel0wuVJnvR2kGoYw2TADNZWHY863rJYNvNfT2dwHAL1dVj8C5omN4fKU11cfRBbSi0DRrGCmZqNsCqjWAj%2Brg2J%2FBlojKl7b8OXUGVdoefwvUnWiPWx2KZ4iiZbGQZ9Tn2lrvYX4o1RxDm22LQGBZ9iH8RWGgXhBWUbz%2Br07uYS16HasIo5ElyWi5lbU%2BShyCx%2F29N6o%2BxSwmYaoYpj37QoLC7JSuRPqNaxmwphj8%2FLsja9%2BqP3zJ4qPo1dWgeGpyHIz6MGEkOT7PoG2GMzvhhWPpeAaz%2BvjIVq4vp0YcgjNY47yX6uCiVMco%2FUla7H97Di3d1suihnvzmBc6zNsekN0fk3sI50PVmT3o8Juk21hiDi3SKMh4YHsFsZxBTSQK80%2F1%2FF1onKBl0w8hkqcAJkr8n8eahlYVc%2FAFM%2FpIg705icrpadaXqRaQy3D0wknc4maq13ta0l%2BY8tDWiiaTYQ6gail9%2BSJ8qn1XXCnqXoS9OTXpu%2F13BlDtNg58p0WGY8fa8n5ZKWjJFtYNSZrCiU8B2zdIm42v6cn%2FpsnxvGX9ZTMtuvK55IlQpfn8npIcu%2F8KpI%2FBADcnKp9%2BR%2BFJ97EhOLlKGGB%2Fm0cDL2iyxEg87JfZtV4ftaZsrpxJD0wmI93vyXnMaq%2BOCnpM5g9w0E8DMvXcZCD%2FDeJs7itM2yf0l4WW0BYfpuvkT7A0UK32MEPGK8JsSyFgXSutahuX%2FmHkwcflCirnGbvuDyYEXauxtlOZBbG3o09b8Ps69EQkY337%2FLJfE9dDTvM9cGqMNWDlaYGOrEB3dUKkM0vTvNlkqVELzqmPe4A4GLMTymhQALZAAQgw3q%2FAPaMzn%2BiCD%2BqFolHz7iIjH%2BosaOwsaHQ8jJ6GDjZTeJCabKDnPkWTe0R7M%2BNhYTuz3NX1lJPYP18YaJa6xxxf%2FZbsRiS%2BC3GO8np1RkrqpzpRN5lP1PwawKDPJ%2BzGlzDACV8UAjV3BSN%2F0WDFANtTIEQR60j8MiA3puhO0d14nqnA3YpBF1VHG4nOEG7zvZI&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230729T171839Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYVSCFSRUF%2F20230729%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=d3e996546b9fb2fc6730ecec4ef32e70af9d274bdb2f6f9721e4d5363fafb841&hash=ab93a237b2b9fc8bda11f836462064bff6ce1c23ce72df9f7281933cf46899c4&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0010465510004509&tid=spdf-934037fa-0f3e-4fa5-9a9a-3dd8df09af23&sid=1c15774d53e2674b2989a4e4d32def85da1cgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=050c5101005951510304&rr=7ee705ba3d67043e&cc=cn)
2. Dice, Kevin, et al. "CUDA-Accelerated Simulation of Brownian Dynamics." Proceedings of the Practice and Experience on Advanced Research Computing. 2018. 1-3. [Link](https://dl.acm.org/doi/abs/10.1145/3219104.3229260)

Please refer to the above reference for more details on the algorithm design and principles.

