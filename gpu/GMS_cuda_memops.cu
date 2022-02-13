
#include <cuda.h>
#include <cuda_runtime.h>
#include "GMS_cuda_memops.cuh"


//
//	Implementation
//



void copy_int32_cpu_to_gpu(int32_t * __restrict d_ptr, 
			   const int32_t * __restrict h_ptr,
			   const size_t n,
			   int32_t * ierr ) {
	if(*ierr <= 0) *ierr = 0;
	if( NULL == h_ptr || n <= 0){ //  Host error handling
	   *ierr = -1;
	   return; 
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n * sizeof(int32_t),cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(int32_t),cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
Error:
	cudaFree(d_ptr);

	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
					status);
}






void copy_float_cpu_to_gpu(float * __restrict d_ptr,
			   const float * __restrict h_ptr,
			   const size_t n,
			   int32_t * ierr ) {
	if(*ierr <= 0) *ierr = 0;
	if(NULL == h_ptr || n <= 0){
	   *ierr = -1;
	   return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(float),
					cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(float),
					cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}






void copy_double_cpu_to_gpu(double * __restrict d_ptr,
			    const double * __restrict h_ptr,
			    const size_t n,
			    int32_t * ierr ) {
	if(*ierr <= 0) *ierr = 0;
	if(NULL == h_ptr ||
	   0 >= n     ) {
	    *ierr = -1;
	    return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(double),
					cudaMemcpyHostToDevice));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
	GMS_CUDA_CHECK(cudaMemcpy((void*)&d_ptr[0],(void*)&h_ptr[0],n*sizeof(double),
					cudaMemcpyHostToDevice));
#endif
	*ierr = 0;
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}



//
//	Allocate memory on GPU.
//

void alloc_int32_gpu(int32_t * __restrict d_ptr,
                     const size_t n,
		     int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(int32_t)));
#endif
	*ierr = 0;
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}







void alloc_float_gpu(float * __restrict d_ptr,
		     const size_t n,
		     int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n ) {

		*ierr = -1;
		return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(float)));
#endif
	*ierr = 0;
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}





void alloc_double_gpu(double * __restrict d_ptr,
		      const size_t n,
		      int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n ) {

		*ierr = -1;
		return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
#else
	GMS_CUDA_CHECK(cudaMalloc((void**)&d_ptr,n*sizeof(double)));
#endif
	*ierr = 0;
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}



//
// GPU to CPU memory copy routines
//


//
// Copy array  of int32_t from GPU to CPU.
//

void copy_int32_gpu_to_cpu(const int32_t * __restrict d_ptr,
			   int32_t * __restrict h_ptr,
			   const int32_t n,
			   int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(int32_t),
					cudaMemcpyDeviceToHost));
	*ierr = 0;
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(int32_t),
					cudaMemcpyDeviceToHost));
	*ierr = 0;
#endif
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}



void copy_float_gpu_to_cpu(const float * __restrict d_ptr,
			   float * __restrict h_ptr,
			   const int32_t n,
			   int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(float),
					cudaMemcpyDeviceToHost));
	*ierr = 0;
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(float),
					cudaMemcpyDeviceToHost));
	*ierr = 0;
#endif
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}



void copy_double_gpu_to_cpu(const double * __restrict d_ptr,
			    double * __restrict h_ptr,
			    const int32_t n,
			    int32_t * ierr ) {
	if(*ierr < 0) *ierr = 0;
	if(0 >= n) {

		*ierr = -1;
		return;
	}
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	GMS_CUDA_DEBUG_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(double),
					cudaMemcpyDeviceToHost));
	*ierr = 0;
#else
	GMS_CUDA_CHECK(cudaMemcpy(&h_ptr[0],&d_ptr[0],n*sizeof(double),
			   cudaMemcpyDeviceToHost));
	*ierr = 0;
#endif
Error:
	cudaFree((void*)&d_ptr[0]);
	*ierr = -2;
	fatal_gpu_error(__PRETTY_FUNCTION__,
		            status);
}

























