// #include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "logbesselk.h"

#define CHUNKSIZE 16  // You can adjust this value

__device__ double modified_besselk(double nu, double x);

__global__ void dcmg_refined_matern_kernel(double *A, int m, int n,
        double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
        double localtheta0, double localtheta1, double localtheta2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        double sigma_square = localtheta0;
        double range = localtheta1;
        double smooth = localtheta2;

        double dx = l1_x_cuda[i] - l2_x_cuda[j];
        double dy = l1_y_cuda[i] - l2_y_cuda[j];
        double distance = sqrt(dx*dx + dy*dy);

        if (distance == 0) {
            A[i + j * m] = sigma_square;
        } else {
            double scaled_dist = distance / range;
            double coeff = sigma_square / (pow(2, smooth - 1) * tgamma(smooth));
            A[i + j * m] = coeff * pow(scaled_dist, smooth) * modified_besselk(smooth, scaled_dist);
        }
    }
}

extern "C" void cudacmg_refined_matern(double *A, int m, int n,
        double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
        double *localtheta, cudaStream_t stream)
{
    dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
    dim3 dimGrid((m + CHUNKSIZE - 1) / CHUNKSIZE, (n + CHUNKSIZE - 1) / CHUNKSIZE);

    dcmg_refined_matern_kernel<<<dimGrid, dimBlock, 0, stream>>>(
        A, m, n, l1_x_cuda, l1_y_cuda, l2_x_cuda, l2_y_cuda,
        localtheta[0], localtheta[1], localtheta[2]);

    cudaStreamSynchronize(stream);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err));
    }
}