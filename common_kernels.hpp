#ifndef COMMON_KERNELS
#define COMMON_KERNELS

#include "immersed_boundary_method.hpp"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

__device__ double3 operator+(const double3 &a, const double3 &b);
__device__ double3 operator-(const double3 &a, const double3 &b);
__device__ double3 operator/(const double3 &a, const double b);
__device__ double3 operator*(const double3 &a, const double b);
__device__ double3 set_equal(const double3 &a, const double b);

__device__ double dot_product(const double3 &a, const double3 &b);
__device__ double3 cross_product(const double3 &a, const double3 &b);
__device__ double myatomicAdd(double* address, double val);
__device__ double3 triangle_centroid(const double3 &a, const double3 &b, const double3 &c);
__device__ double norm(const double a, const double b, const double c);
__global__ void add_test(int n, double* delta_t, double3* area);
__global__ void clone_a_to_b(int n_cells, double4* a, double4* b);
__global__ void clone_a_to_b(int n_cells, double * a, double * b);
__global__ void fill_zero(int n_cells, double* a);

__global__ void fill_double(int n_cells, double* a, double val);
__global__ void fill_zero(int n_cells, double4* a);
__global__ void square(int n_cells, double* a);
__global__ void add(int n_cells, double* a, double *b);
__global__ void multiply(int n_cells, double* a, double b);
__global__ void divide (int n_cells, double* a, double * b);
__global__ void add_double(int n_cells, double* a, double val);

extern __constant__ double lattice_weight[15];


template <unsigned int blockSize>
__device__ void warpReduce(volatile double *sdata, unsigned int tid)
{
	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}

template <unsigned int blockSize>
__global__ void reduce6(double *g_idata, double *g_odata, unsigned int n) {
	extern __shared__ double sdata[2 * blockSize];
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockSize * 2) + tid;
	unsigned int gridSize = blockSize * 2 * gridDim.x;
	sdata[tid] = 0;
	while (i < n) {
		sdata[tid] += g_idata[i] + g_idata[i + blockSize];
		i += gridSize;
	}
	__syncthreads();
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
	if (tid < 32) warpReduce< blockSize>(sdata, tid);
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <unsigned int BLOCK_SIZE>
__global__ void total(double * input, double * output, int len) {
	//@@ Load a segment of the input vector into shared memory
	__shared__ double partialSum[2 * BLOCK_SIZE];
	unsigned int t = threadIdx.x, start = 2 * blockIdx.x * BLOCK_SIZE;
	if (start + t < len)
		partialSum[t] = input[start + t];
	else
		partialSum[t] = 0;
	if (start + BLOCK_SIZE + t < len)
		partialSum[BLOCK_SIZE + t] = input[start + BLOCK_SIZE + t];
	else
		partialSum[BLOCK_SIZE + t] = 0;
	//@@ Traverse the reduction tree
	for (unsigned int stride = BLOCK_SIZE; stride >= 1; stride >>= 1) {
		__syncthreads();
		if (t < stride)
			partialSum[t] += partialSum[t + stride];
	}
	//@@ Write the computed sum of the block to the output vector at the
	//@@ correct index
	if (t == 0)
		output[blockIdx.x] = partialSum[0];
}

template <unsigned int BLOCK_SIZE>
__global__ void sum_product(double * input, double * input2,  double * output, int len) {
	//@@ Load a segment of the input vector into shared memory
	__shared__ double partialSum[2 * BLOCK_SIZE];
	unsigned int t = threadIdx.x, start = 2 * blockIdx.x * BLOCK_SIZE;
	if (start + t < len)
		partialSum[t] = input[start + t]* input2[start +t];
	else
		partialSum[t] = 0;
	if (start + BLOCK_SIZE + t < len)
		partialSum[BLOCK_SIZE + t] = input[start + BLOCK_SIZE + t] * input2[start + BLOCK_SIZE + t];
	else
		partialSum[BLOCK_SIZE + t] = 0;
	//@@ Traverse the reduction tree
	for (unsigned int stride = BLOCK_SIZE; stride >= 1; stride >>= 1) {
		__syncthreads();
		if (t < stride)
			partialSum[t] += partialSum[t + stride];
	}
	//@@ Write the computed sum of the block to the output vector at the
	//@@ correct index
	if (t == 0)
		output[blockIdx.x] = partialSum[0];
}


template <unsigned int BLOCK_SIZE>
__global__ void sum_product_curvature(double * input, double * input2, double * output, int len) {
	//@@ Load a segment of the input vector into shared memory
	__shared__ double partialSum[2 * BLOCK_SIZE];
	unsigned int t = threadIdx.x, start = 2 * blockIdx.x * BLOCK_SIZE;
	if (start + t < len)
		partialSum[t] = input[start + t] * (input2[(start + t) *2 ] + input2[(start + t) * 2 + 1]) * 1/3;
	else
		partialSum[t] = 0;
	if (start + BLOCK_SIZE + t < len)
		partialSum[BLOCK_SIZE + t] = input[start + BLOCK_SIZE + t] * (input2[(start + BLOCK_SIZE + t) * 2] + input2[(start + BLOCK_SIZE + t) * 2 + 1]) * 1 / 3;
	else
		partialSum[BLOCK_SIZE + t] = 0;
	//@@ Traverse the reduction tree
	for (unsigned int stride = BLOCK_SIZE; stride >= 1; stride >>= 1) {
		__syncthreads();
		if (t < stride)
			partialSum[t] += partialSum[t + stride];
	}
	//@@ Write the computed sum of the block to the output vector at the
	//@@ correct index
	if (t == 0)
		output[blockIdx.x] = partialSum[0];
}

template <unsigned int BLOCK_SIZE>
__global__ void min(double * input, double * output, int len) {
	//@@ Load a segment of the input vector into shared memory
	__shared__ double partialMin[2 * BLOCK_SIZE];
	unsigned int t = threadIdx.x, start = 2 * blockIdx.x * BLOCK_SIZE;
	if (start + t < len)
		partialMin[t] = input[start + t];
	else
		partialMin[t] = 100000000;
	if (start + BLOCK_SIZE + t < len)
		partialMin[BLOCK_SIZE + t] = input[start + BLOCK_SIZE + t];
	else
		partialMin[BLOCK_SIZE + t] = 100000000;
	//@@ Traverse the reduction tree
	for (unsigned int stride = BLOCK_SIZE; stride >= 1; stride >>= 1) {
		__syncthreads();
		if (t < stride)
			partialMin[t] = min(partialMin[t + stride], partialMin[t]);
	}
	//@@ Write the computed sum of the block to the output vector at the
	//@@ correct index
	if (t == 0)
		output[blockIdx.x] = partialMin[0];
}


template <unsigned int BLOCK_SIZE>
__global__ void max(double * input, double * output, int len) {
	//@@ Load a segment of the input vector into shared memory
	__shared__ double partialMax[2 * BLOCK_SIZE];
	unsigned int t = threadIdx.x, start = 2 * blockIdx.x * BLOCK_SIZE;
	if (start + t < len)
		partialMax[t] = input[start + t];
	else
		partialMax[t] = -100000000;
	if (start + BLOCK_SIZE + t < len)
		partialMax[BLOCK_SIZE + t] = input[start + BLOCK_SIZE + t];
	else
		partialMax[BLOCK_SIZE + t] = -100000000;
	//@@ Traverse the reduction tree
	for (unsigned int stride = BLOCK_SIZE; stride >= 1; stride >>= 1) {
		__syncthreads();
		if (t < stride)
			partialMax[t] = max(partialMax[t + stride], partialMax[t]);
	}
	//@@ Write the computed sum of the block to the output vector at the
	//@@ correct index
	if (t == 0)
		output[blockIdx.x] = partialMax[0];
}


#endif