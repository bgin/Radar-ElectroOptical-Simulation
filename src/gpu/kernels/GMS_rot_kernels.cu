

#include <cstdint>

// Modified by Bernard Gingold, contact: bgin@gmail.com
/*
    GLAS is licensed under the The MIT License:

Copyright (c) 2011 Hans Henrik Brandenborg Sørensen, DTU

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
// Generic vector scale and add kernel.
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

// Version 1: copy as floats
// -----------------------------------------
// Threads per block    : BLOCKSIZE_n
// Work per thread      : WORKSIZE_n
// Shared mem per block : 0
//


__global__ void rot_kernel_r4(float* __restrict__ x, 
                              float* __restrict__ y, 
                              const float c, 
                              const float s, 
                              const int32_t n, 
                              const int32_t ntile,
                              const int32_t BLOCKSIZE_n,
                              const int32_t WORKSIZE_m,
                              const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			float vx = x[0], vy = y[0]; 
			x[0] = c * vx + s * vy;
			y[0] = c * vy - s * vx;
			x += BLOCKSIZE_n; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			float vx = x[i], vy = y[i]; 
			x[i] = c * vx + s * vy;
			y[i] = c * vy - s * vx;
		}
	}
}

__global__ void rot_kernel_c4(cuComplex* __restrict__ x, 
                              cuComplex* __restrict__ y, 
                              const cuComplex c, 
                              const cuComplex s, 
                              const int32_t n, 
                              const int32_t ntile,
                              const int32_t BLOCKSIZE_n,
                              const int32_t WORKSIZE_m,
                              const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			cuComplex vx = x[0], vy = y[0]; 
			x[0] = c * vx + s * vy;
			y[0] = c * vy - s * vx;
			x += BLOCKSIZE_n; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			cuComplex vx = x[i], vy = y[i]; 
			x[i] = c * vx + s * vy;
			y[i] = c * vy - s * vx;
		}
	}
}

__global__ void rot_kernel_r8(double* __restrict__ x, 
                              double* __restrict__ y, 
                              const double c, 
                              const double s, 
                              const int32_t n, 
                              const int32_t ntile,
                              const int32_t BLOCKSIZE_n,
                              const int32_t WORKSIZE_m,
                              const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx + tid;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			double vx = x[0], vy = y[0]; 
			x[0] = c * vx + s * vy;
			y[0] = c * vy - s * vx;
			x += BLOCKSIZE_n; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			double vx = x[i], vy = y[i]; 
			x[i] = c * vx + s * vy;
			y[i] = c * vy - s * vx;
		}
	}
}

// Version 2: copy as doubles (requires WORKSIZE_n & 1 == 0)
// -----------------------------------------
// Threads per block    : BLOCKSIZE_n
// Work per thread      : WORKSIZE_n / 2
// Shared mem per block : 0
//
//template <int BLOCKSIZE_n, int WORKSIZE_m, int WORKSIZE_n>
//__launch_bounds__(BLOCKSIZE_n, MIN_BLOCKS_PER_MP)
//__global__ void srot_kernel2(float* x, float* y, const float c, const float s, int n, int ntile)
//{
/*
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	x += idx;
	y += idx;
	const double *xd = (double *) x + tid;
	double *yd = (double *) y + tid;
	__shared__ double dx_sh[BLOCKSIZE_n];
	__shared__ double dy_sh[BLOCKSIZE_n];
	double *dx_s = dx_sh + tid;
	double *dy_s = dy_sh + tid;
	float *x_s = (float *) dx_s;
	float *y_s = (float *) dy_s;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < (WORKSIZE_n / 2); i++ )
		{
			dx_s[0] = xd[0]; dy_s[0] = yd[0]; 
			y_s[0] += alpha * x_s[0]; y_s[1] += alpha * x_s[1]; 
			yd[0] = *dy_s; xd += BLOCKSIZE_n; yd += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < ((n - idx) >> 1) - tid; i += BLOCKSIZE_n)
		{
			dx_s[0] = xd[i]; 
			dy_s[0] = yd[i]; 
			y_s[0] += alpha * x_s[0]; 
			y_s[1] += alpha * x_s[1]; 
			yd[i] = *dy_s;
		}
	}
	if (blockIdx.x == 0 && tid == 0 && n & 1) y[n - 1] += alpha * x[n - 1];
*/
//}

