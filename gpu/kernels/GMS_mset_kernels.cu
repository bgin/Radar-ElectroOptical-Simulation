
#include <cstdint>
// Modified by Bernard Gingold, contact: bgin@gmail.com
/*
    GLAS is licensed under the The MIT License:

Copyright (c) 2011 Hans Henrik Brandenborg Sørensen, DTU

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
 

// Generic vector set kernel.
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


__global__ void mset_kernel_r4(float val, 
                               float* __restrict__ y, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			y[0] = val; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			y[i] = val;
		}
	}
}


__global__ void mset_kernel_c4(cuComplex val, 
                               cuComplex* __restrict__ y, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			y[0] = val; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			y[i] = val;
		}
	}
}

__global__ void mset_kernel_r8(double val, 
                               double* __restrict__ y, 
                               const int32_t n, 
                               const int32_t ntile,
                               const int32_t BLOCKSIZE_n,
                               const int32_t WORKSIZE_m,
                               const int32_t WORKSIZE_n)
{
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLOCKSIZE_n * WORKSIZE_n;
	y += idx + tid;
	if (idx != ntile)
	{
        #pragma unroll
		for (int i = 0; i < WORKSIZE_n; i++ )
		{
			y[0] = val; y += BLOCKSIZE_n;
		}
	}
	else
	{
        #pragma unroll 4
		for (int i = 0; i < n - idx - tid; i += BLOCKSIZE_n)
		{
			y[i] = val;
		}
	}
}



