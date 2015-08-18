// fftGPU.cu - CUDA Code for FFTs on GPU using cuFFT
// Implemented by Alex Cerjanic on 8/17/2015
// Based on https://www.olcf.ornl.gov/tutorials/mixing-openacc-with-gpu-libraries/


template<typename T1>
void ifft2dGPU(T1 *d_data, int nx, int ny, void *stream);

template<typename T1>
void ifft3dGPU(T1 *d_data,  int nx, int ny, int nz, void *stream);

template<typename T1>
void fft3dGPU(T1 *d_data,  int nx, int ny, int nz, void *stream);

template<typename T1>
void fft2dGPU(T1 *d_data, int nx, int ny, void *stream);