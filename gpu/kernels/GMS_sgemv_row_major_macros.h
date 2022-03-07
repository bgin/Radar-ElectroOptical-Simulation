#define initialize_sums(WORKSIZE_m)								\
	float sum[WORKSIZE_m];										\
    if (WORKSIZE_m >= 1) sum[0] = 0.0;							\
	if (WORKSIZE_m >= 2) sum[1] = 0.0;							\
	if (WORKSIZE_m >= 3) sum[2] = 0.0;							\
	if (WORKSIZE_m >= 4) sum[3] = 0.0;							\
	if (WORKSIZE_m >= 5) sum[4] = 0.0;							\
	if (WORKSIZE_m >= 6) sum[5] = 0.0;							\
	if (WORKSIZE_m >= 7) sum[6] = 0.0;							\
	if (WORKSIZE_m >= 8) sum[7] = 0.0;

#define add_to_sums(lda, WORKSIZE_m)			\
	if (WORKSIZE_m >= 1) sum[0] += A[0] * x[0];							\
	if (WORKSIZE_m >= 2) sum[1] += A[lda] * x[0]; \
	if (WORKSIZE_m >= 3) sum[2] += A[2 * lda] * x[0]; \
	if (WORKSIZE_m >= 4) sum[3] += A[3 * lda] * x[0]; \
	if (WORKSIZE_m >= 5) sum[4] += A[4 * lda] * x[0]; \
	if (WORKSIZE_m >= 6) sum[5] += A[5 * lda] * x[0]; \
	if (WORKSIZE_m >= 7) sum[6] += A[6 * lda] * x[0]; \
	if (WORKSIZE_m >= 8) sum[7] += A[7 * lda] * x[0];

#define add_to_sums_check_bounds(row, mrest, lda, WORKSIZE_m)			\
	if (WORKSIZE_m >= 1) sum[0] += A[0] * x[0];							\
	if (WORKSIZE_m >= 2 && row + 1 < mrest) sum[1] += A[lda] * x[0]; \
	if (WORKSIZE_m >= 3 && row + 2 < mrest) sum[2] += A[2 * lda] * x[0]; \
	if (WORKSIZE_m >= 4 && row + 3 < mrest) sum[3] += A[3 * lda] * x[0]; \
	if (WORKSIZE_m >= 5 && row + 4 < mrest) sum[4] += A[4 * lda] * x[0]; \
	if (WORKSIZE_m >= 6 && row + 5 < mrest) sum[5] += A[5 * lda] * x[0]; \
	if (WORKSIZE_m >= 7 && row + 6 < mrest) sum[6] += A[6 * lda] * x[0]; \
	if (WORKSIZE_m >= 8 && row + 7 < mrest) sum[7] += A[7 * lda] * x[0];

#define atomicAdd_(arg, sum) atomicAdd(y + arg, sum)

#define reduce_and_atomicAdd_sums(BLOCKSIZE)							\
	if (WORKSIZE_m >= 1) { block_reduce<BLOCKSIZE>(sdata, smem, sum[0], tid); if (tid == 0) atomicAdd_(0, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 2) { block_reduce<BLOCKSIZE>(sdata, smem, sum[1], tid); if (tid == 0) atomicAdd_(1, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 3) { block_reduce<BLOCKSIZE>(sdata, smem, sum[2], tid); if (tid == 0) atomicAdd_(2, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 4) { block_reduce<BLOCKSIZE>(sdata, smem, sum[3], tid); if (tid == 0) atomicAdd_(3, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 5) { block_reduce<BLOCKSIZE>(sdata, smem, sum[4], tid); if (tid == 0) atomicAdd_(4, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 6) { block_reduce<BLOCKSIZE>(sdata, smem, sum[5], tid); if (tid == 0) atomicAdd_(5, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 7) { block_reduce<BLOCKSIZE>(sdata, smem, sum[6], tid); if (tid == 0) atomicAdd_(6, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 8) { block_reduce<BLOCKSIZE>(sdata, smem, sum[7], tid); if (tid == 0) atomicAdd_(7, alpha * sdata[0]); __syncthreads(); }

#define reduce_and_atomicAdd_sums_check_bounds(row, mrest, BLOCKSIZE)	\
	if (WORKSIZE_m >= 1) { block_reduce<BLOCKSIZE>(sdata, smem, sum[0], tid); if (tid == 0) atomicAdd_(0, alpha * sdata[0]); __syncthreads(); }	\
	if (WORKSIZE_m >= 2 && row + 1 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[1], tid); if (tid == 0) atomicAdd_(1, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 3 && row + 2 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[2], tid); if (tid == 0) atomicAdd_(2, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 4 && row + 3 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[3], tid); if (tid == 0) atomicAdd_(3, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 5 && row + 4 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[4], tid); if (tid == 0) atomicAdd_(4, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 6 && row + 5 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[5], tid); if (tid == 0) atomicAdd_(5, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 7 && row + 6 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[6], tid); if (tid == 0) atomicAdd_(6, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 8 && row + 7 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[7], tid); if (tid == 0) atomicAdd_(7, alpha * sdata[0]); __syncthreads(); }

#define reduce_2d_and_atomicAdd_sums(BLOCKSIZE)							\
	if (WORKSIZE_m >= 1) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[0], tid, tidx); if (tidx == 0) atomicAdd_(0, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 2) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[1], tid, tidx); if (tidx == 0) atomicAdd_(1, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 3) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[2], tid, tidx); if (tidx == 0) atomicAdd_(2, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 4) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[3], tid, tidx); if (tidx == 0) atomicAdd_(3, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 5) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[4], tid, tidx); if (tidx == 0) atomicAdd_(4, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 6) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[5], tid, tidx); if (tidx == 0) atomicAdd_(5, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 7) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[6], tid, tidx); if (tidx == 0) atomicAdd_(6, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 8) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[7], tid, tidx); if (tidx == 0) atomicAdd_(7, alpha * sdata[tid]); __syncthreads(); }

#define reduce_2d_and_atomicAdd_sums_check_bounds(row, mrest, BLOCKSIZE)	\
	if (WORKSIZE_m >= 1) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[0], tid, tidx); if (tidx == 0) atomicAdd_(0, alpha * sdata[tid]); __syncthreads(); }	\
	if (WORKSIZE_m >= 2 && row + 1 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[1], tid, tidx); if (tidx == 0) atomicAdd_(1, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 3 && row + 2 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[2], tid, tidx); if (tidx == 0) atomicAdd_(2, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 4 && row + 3 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[3], tid, tidx); if (tidx == 0) atomicAdd_(3, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 5 && row + 4 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[4], tid, tidx); if (tidx == 0) atomicAdd_(4, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 6 && row + 5 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[5], tid, tidx); if (tidx == 0) atomicAdd_(5, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 7 && row + 6 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[6], tid, tidx); if (tidx == 0) atomicAdd_(6, alpha * sdata[tid]); __syncthreads(); } \
	if (WORKSIZE_m >= 8 && row + 7 < mrest) { block_reduce_2d<BLOCKSIZE>(sdata, smem, sum[7], tid, tidx); if (tidx == 0) atomicAdd_(7, alpha * sdata[tid]); __syncthreads(); }

#define inner_loop(nrest)												\
	for (int k = 0; k < nrest; k += BLOCKSIZE_n)						\
	{																	\
		add_to_sums(lda, WORKSIZE_m);									\
		A += BLOCKSIZE_n; x += BLOCKSIZE_n;								\
	}

#define inner_loop_check_bounds(row, mrest, nrest)						\
	for (int k = 0; k < nrest; k += BLOCKSIZE_n)						\
	{																	\
		add_to_sums_check_bounds(row, mrest, lda, WORKSIZE_m);			\
		A += BLOCKSIZE_n; x += BLOCKSIZE_n;								\
	}

// Several rows per block special

#define atomicAdd_sr(arg, sum) y[arg] += sum

#define reduce_and_atomicAdd_sums_sr(BLOCKSIZE)							\
	if (WORKSIZE_m >= 1) { block_reduce<BLOCKSIZE>(sdata, smem, sum[0], tid); if (tid == 0) atomicAdd_sr(0, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 2) { block_reduce<BLOCKSIZE>(sdata, smem, sum[1], tid); if (tid == 0) atomicAdd_sr(1, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 3) { block_reduce<BLOCKSIZE>(sdata, smem, sum[2], tid); if (tid == 0) atomicAdd_sr(2, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 4) { block_reduce<BLOCKSIZE>(sdata, smem, sum[3], tid); if (tid == 0) atomicAdd_sr(3, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 5) { block_reduce<BLOCKSIZE>(sdata, smem, sum[4], tid); if (tid == 0) atomicAdd_sr(4, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 6) { block_reduce<BLOCKSIZE>(sdata, smem, sum[5], tid); if (tid == 0) atomicAdd_sr(5, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 7) { block_reduce<BLOCKSIZE>(sdata, smem, sum[6], tid); if (tid == 0) atomicAdd_sr(6, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 8) { block_reduce<BLOCKSIZE>(sdata, smem, sum[7], tid); if (tid == 0) atomicAdd_sr(7, alpha * sdata[0]); __syncthreads(); }

#define reduce_and_atomicAdd_sums_check_bounds_sr(row, mrest, BLOCKSIZE)	\
	if (WORKSIZE_m >= 1) { block_reduce<BLOCKSIZE>(sdata, smem, sum[0], tid); if (tid == 0) atomicAdd_sr(0, alpha * sdata[0]); __syncthreads(); }	\
	if (WORKSIZE_m >= 2 && row + 1 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[1], tid); if (tid == 0) atomicAdd_sr(1, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 3 && row + 2 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[2], tid); if (tid == 0) atomicAdd_sr(2, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 4 && row + 3 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[3], tid); if (tid == 0) atomicAdd_sr(3, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 5 && row + 4 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[4], tid); if (tid == 0) atomicAdd_sr(4, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 6 && row + 5 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[5], tid); if (tid == 0) atomicAdd_sr(5, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 7 && row + 6 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[6], tid); if (tid == 0) atomicAdd_sr(6, alpha * sdata[0]); __syncthreads(); } \
	if (WORKSIZE_m >= 8 && row + 7 < mrest) { block_reduce<BLOCKSIZE>(sdata, smem, sum[7], tid); if (tid == 0) atomicAdd_sr(7, alpha * sdata[0]); __syncthreads(); }

// Several rows per block shared mem special

#define copy_A_sh(lds, lda, WORKSIZE_m)			\
	if (WORKSIZE_m >= 1) A_s[i] = A[i];\
	if (WORKSIZE_m >= 2) A_s[i + lds] = A[i + lda];\
	if (WORKSIZE_m >= 3) A_s[i + 2 * lds] = A[i + 2 * lda];\
	if (WORKSIZE_m >= 4) A_s[i + 3 * lds] = A[i + 3 * lda];\
	if (WORKSIZE_m >= 5) A_s[i + 4 * lds] = A[i + 4 * lda];\
	if (WORKSIZE_m >= 6) A_s[i + 5 * lds] = A[i + 5 * lda];\
	if (WORKSIZE_m >= 7) A_s[i + 6 * lds] = A[i + 6 * lda];\
	if (WORKSIZE_m >= 8) A_s[i + 7 * lds] = A[i + 7 * lda];\

#define add_to_sums_sh(lda, WORKSIZE_m)			\
	if (WORKSIZE_m >= 1) sum[0] += A_s[0] * x_s[0];	  \
	if (WORKSIZE_m >= 2) sum[1] += A_s[lda] * x_s[0]; \
	if (WORKSIZE_m >= 3) sum[2] += A_s[2 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 4) sum[3] += A_s[3 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 5) sum[4] += A_s[4 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 6) sum[5] += A_s[5 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 7) sum[6] += A_s[6 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 8) sum[7] += A_s[7 * lda] * x_s[0];

#define add_to_sums_check_bounds_sh(row, mrest, lda, WORKSIZE_m)		\
	if (WORKSIZE_m >= 1) sum[0] += A_s[0] * x_s[0];						\
	if (WORKSIZE_m >= 2 && row + 1 < mrest) sum[1] += A_s[lda] * x_s[0];	\
	if (WORKSIZE_m >= 3 && row + 2 < mrest) sum[2] += A_s[2 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 4 && row + 3 < mrest) sum[3] += A_s[3 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 5 && row + 4 < mrest) sum[4] += A_s[4 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 6 && row + 5 < mrest) sum[5] += A_s[5 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 7 && row + 6 < mrest) sum[6] += A_s[6 * lda] * x_s[0]; \
	if (WORKSIZE_m >= 8 && row + 7 < mrest) sum[7] += A_s[7 * lda] * x_s[0];

#define inner_loop_sh(nrest)											\
	for (int k = 0; k < nrest; k += BLOCKSIZE_n)						\
	{																	\
		add_to_sums_sh(lds, WORKSIZE_m);								\
		A_s += BLOCKSIZE_n; x_s += BLOCKSIZE_n;							\
	}

#define inner_loop_check_bounds_sh(row, mrest, nrest)					\
	for (int k = 0; k < nrest; k += BLOCKSIZE_n)						\
	{																	\
		add_to_sums_check_bounds_sh(row, mrest, lds, WORKSIZE_m);		\
		A_s += BLOCKSIZE_n; x_s += BLOCKSIZE_n;							\
	}
