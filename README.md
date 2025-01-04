# The code will be made public after the paper is accepted.
# SIMP-HMGCG
A hybrid topology optimization method based on SIMP and MGCG: A robust topology optimization method for large-scale structures considering multiple loading conditions

# Single GPU
This is a method capable of achieving nearly 67 million elements and 201 million large-scale topology optimization iterations on a single personal computer equipped with a CPU with 32GB of memory and a GPU with 24GB of memory. 

It employs a hybrid CPU and GPU computing approach, with the parallel computing component implemented on the GPU, and commonly used function functionalities realized using MATLAB software. Scholars and designers can integrate this method with their own proposed new algorithms using MATLAB software to generate more intriguing 3D large-scale topology optimization methods. Consequently, many people will be able to utilize this low-cost large-scale topology optimization program to create interesting projects.

# Environmental configuration
Desktop computer :  GPU Nvidia 4090-24G  Matlab Version: R2023a OS: Windows 10  CUDA Version:DriverVersion: 12.0  MEX Compile Environment:MEX is configured to use 'Microsoft Visual C++ 2019 (C)' for C language compilation.

Laptop: GPU Nvidia 3060-6G  Matlab Version: R2022a OS: Windows 10   CUDA Version:DriverVersion: 12.1   MEX Compile Environment:MEX is configured to use 'Microsoft Visual C++ 2019 (C)' for C language compilation.

After the paper is accepted, I will be recording a video tutorial on how to download MATLAB software and CUDA onto a regular computer, and subsequently compiling a .cu file using MexCUDA to run large-scale topology optimization programs successfully.

# Programming Language
The majority of the code components are completed within the MATLAB software, the parallel computing tasks are performed on the GPU, and the results are returned to MATLAB software using the mexFunction function  for further CPU-based processing. 
