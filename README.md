# LBFS-IBM-ADE-SP
Lattice Bolzmann Flux Solver with Spring Particle RBC model

Dependencies:

* Boost
* CUDA
* TECIO - Tecplot Input/Output library https://github.com/su2code/SU2/tree/master/externals/tecio
          	
 
Set Up:
* Install Dependencies
* Modify Makefile provided to point to the correct librarys and sources
* Make

Input: 
* User Defined XML file
* OpenFoam Mesh directory
* Private message if required

Main Modules:
* The pre-processing takes place in the program.cpp file before calling the core GPU code. 
* The Solver_gpu.cu file contains the main code of the method. It is written in CUDA and enables the PLBFS to be run on GPUs. It calls all the GPU code from the host
* LBFS.cu contains all code required for LBFS including gradient calculations, flux calculations at face and time integration of cell volumes
* Immersed Boundary Method.cu contains IBM functions including force spreading, velocity integrations, force calc, and object node time integration

Disclaimer: This repository is very rough and ready in places outside the main CUDA kernels as I've been working on this as a solo developer. I will put an effort into refactoring some aspects of the code base if an interest is displayed.
