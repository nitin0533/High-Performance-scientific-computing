# High-performance-scientific-computing
## Resume points
 **Parallelization of 0th Order Generalised Mode Acceleration Method (GMAM) using C code**
1. Performed Model order reduction using Guyan Condensation of a 2-D cantilever subjected to harmonic nodal forces
2. Computed dynamic nodal responses by implementing GMAM using first 10 Eigen-modes found by QR decomposition
3. Acheived 5 times reduction in code runtime by parallelisation using OpenMP, MPI & CUDA C frameworks

**Parallelisation of 0th order Generalised Mode Acceleration Method (GMAM) using C code**
1. Performed Model order reduction using Guyan Condensation of a 2-D cantilever and determined the dynamic response using GMAM
2. Acheived 5 times reduction in code runtime using OpenMP, MPI CUDA C (GPGPU) frameworks in a supercomputer facility at IIT-K

**Parallelization of Generalized Mode Acceleration Method (GMAM)**
1. Developed C code to implement CUDA kernels on GPGPU, and MPI data-sharing frameworks to parallelize GMAM
2. Achieved 5 times reduction in code runtime with OpenMP, MPI and CUDA C (GPGPU) in a supercomputer facility

Reports and Presentation may be referred for verification of all the points. Snips are attached here for quick verification.

![image](https://user-images.githubusercontent.com/71177034/129440079-328fd2b6-98cc-42aa-9878-8a9009247e25.png)

Analysis done on 2-D cantilever, Model order reduction of 2-D cantilever using Guyan condensation, and used first 10 eigen modes for analysis. 

QR Decomposition was used
![image](https://user-images.githubusercontent.com/71177034/129440306-bab955f2-3f6a-4ae8-afab-7ce49a3644f4.png)

The timing study of the code runtime is tabulated below
![image](https://user-images.githubusercontent.com/71177034/129440250-6abf6760-5e42-4f72-9136-136ddbb47d66.png)
![image](https://user-images.githubusercontent.com/71177034/129440256-2f704119-3277-40ba-a95a-150056bb74e2.png)

CUDA kernels and MPI data sharing frameworks were used, kindly refer to these images
![image](https://user-images.githubusercontent.com/71177034/129440332-dcad5e02-0191-4aea-ae24-42800d4344af.png)
![image](https://user-images.githubusercontent.com/71177034/129440340-f0d85449-6f00-4a61-bff5-d7ec09e57ca0.png)
![image](https://user-images.githubusercontent.com/71177034/129440355-7b5bac10-b462-42b9-8a04-dfc2aaf97261.png)

