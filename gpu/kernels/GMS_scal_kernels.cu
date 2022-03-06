#include <cstdint>
// Modified by Bernard Gingold, contact: bgin@gmail.com
/*
    GLAS is licensed under the The MIT License:

Copyright (c) 2011 Hans Henrik Brandenborg Sørensen, DTU

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// Generic vector scale kernel.
//
// - Simple design for all capabilities
// - Handled by splitting matrix as
// --------------------------
// |                 |      |
// |                 |      |
// |                 |      |
// | left tiled part | rest | 
// |                 |      |
// |                 |      |
// |                 |      |
// --------------------------
// - Kernels are generic and work for all cases
//

// Version 1: as floats
// -----------------------------------------
// Threads per block    : BLOCKSIZE_n
// Work per thread      : WORKSIZE_n
// Shared mem per block : 0
//

__global__ void scal_kernel_r4(float alpha, 
                               float* __restrict__ x, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			x[0] = alpha * x[0]; x += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			x[i] = alpha * x[i];
		}
	}
}

__global__ void scal_kernel_c4(cuComplex alpha, 
                               cuComplex* __restrict__ x, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int32_t i = 0; i < WORKSIZE_n; i++ )
		{
			x[0] = alpha * x[0]; x += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int32_t i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			x[i] = alpha * x[i];
		}
	}
}

__global__ void scal_kernel_r8(double alpha, 
                               double* __restrict__ x, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			x[0] = alpha * x[0]; x += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			x[i] = alpha * x[i];
		}
	}
}


// Version 2: as doubles (requires WORKSIZE_n & 1 == 0)
// -----------------------------------------
// Threads per block    : BLOCKSIZE_n
// Work per thread      : WORKSIZE_n / 2
// Shared mem per block : 0
/*template <int BLOCKSIZE_n, int WORKSIZE_m, int WORKSIZE_n>
__launch_bounds__(BLOCKSIZE_n, MIN_BLOCKS_PER_MP)
__global__ void sscal_kernel2(float alpha, float* x, int n, int ntile)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx;
	double *xd = (double *) x + tid;
	__shared__ double dx_sh[BLOCKSIZE_n];
	double *dx_s = dx_sh + tid;
	float *x_s = (float *) dx_s;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < (WORKSIZE_n / 2); i++ )
		{
			dx_s[0] = xd[0];
			x_s[0] = alpha * x_s[0]; 
			x_s[1] = alpha * x_s[1]; 
			xd[0] = *dx_s; 
			xd += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < ((n - idx) >> 1) - tid; i += BLOCKSIZE_n)
		{
			dx_s[0] = xd[i]; 
			x_s[0] = alpha * x_s[0]; 
			x_s[1] = alpha * x_s[1]; 
			xd[i] = *dx_s;
		}
	}
	if (blockIdx.x == 0 && tid == 0 && n & 1) x[n - 1] = alpha * x[n - 1];
}*/

