

#include <immintrin.h>
#include "GMS_init_gpus.cuh"
//
//	Implementation
//

#define TO_KiB(value) ((double)(value)) / (1024.0)

void devcaps_to_screen(const int dev_num, int * ierr, const bool all_info) {



#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of devcaps_to_screen");
#endif
	cudaDeviceProp devp;
	cudaError status;
	
	status = cudaGetDeviceProperties(&devp,dev_num);
	if(status != cudaSuccess) {
		fprintf(stderr,"cudaGetDeviceProperties failed with an error: %s\n",cudaGetErrorString(status));
		*ierr = -1;
		return;
	}
	// Check CUDA device no. x is installed
	if(devp.major == 9999 && devp.minor == 9999) {
		fprintf(stderr,"cudaGetDeviceProperties find out no CUDA device(s) ... terminating!!\n");
		exit_fatal("No CUDA device has been found",
			      "CuWRF/Run-time errors.txt"  );
	}
	printf("Begin dump of CUDA device no: %d\n", dev_num);
	printf("Device ASCII name:			   %s       \n",devp.name);
	printf("Device major compute caps:    %d       \n",devp.major);
	printf("Device minor compute caps:    %d       \n",devp.minor);
	printf("\n");
	printf("CUDA Device HW characteristics: \n");
	printf("Device clock rate:             %d Khz  .\n",devp.clockRate);
	printf("Device async engine count:     %d       \n",devp.asyncEngineCount);
	printf("Device L2 cache size:	        %d bytes.\n",devp.l2CacheSize);
	printf("Device multi-processor count:  %d       \n",devp.multiProcessorCount);
	printf("Device concurrent kernels:     %d       \n",devp.concurrentKernels);
	printf("Device compute mode:           %d       \n",devp.computeMode);
	printf("Device kernel timeout enabled: %d       \n",devp.kernelExecTimeoutEnabled);
	printf("Device PCI Device ID:          %d       \n",devp.pciDeviceID);
	printf("Device PCI Domain ID:		    %d       \n",devp.pciDomainID);
	printf("Device PCI Bus ID:             %d       \n",devp.pciBusID);
	printf("Device ECC enabled:		    %d       \n",devp.ECCEnabled);
	printf("Device host mem mapping:       %d       \n",devp.canMapHostMemory);
	printf("Device registers per block:    %d       \n",devp.regsPerBlock);
	printf("Device warp size:              %d       \n",devp.warpSize);
	printf("Device overlap:                %d       \n",devp.deviceOverlap);
	printf("\n");
	printf("CUDA Device memory characteristics: \n");
	printf("Device memory bus width:       %d       \n",devp.memoryBusWidth);
	printf("Device memory clock rate:      %d       \n",devp.memoryClockRate);
	printf("Device memory pitch:           %d bytes.\n",devp.memPitch);
	printf("Device shared memory/block:    %d bytes.\n",devp.sharedMemPerBlock);
	printf("Device total const memory:     %.9f KiB.\n",TO_KiB(devp.totalConstMem));
	printf("Device total global memory:    %.9f KiB.\n",TO_KiB(devp.totalGlobalMem));
	printf("Device unified addressing:     %d       \n",devp.unifiedAddressing);
	printf("\n");
	printf("CUDA Device compute characteristics: \n");
	printf("Device max threads/block:      %d    \n",devp.maxThreadsPerBlock);
	printf("Device max threads/dimension:        \n");
	printf("Dim 1 = %d threads.\n",devp.maxThreadsDim[0]);
	printf("Dim 2 = %d threads.\n",devp.maxThreadsDim[1]);
	printf("Dim 3 = %d threads.\n",devp.maxThreadsDim[2]);
	printf("Device max thread/multiproc:   %d    \n",devp.maxThreadsPerMultiProcessor);
	printf("Device max size of grid:             \n");
	printf("Grid dim 1 = %d \n",devp.maxGridSize[0]);
	printf("Grid dim 2 = %d \n",devp.maxGridSize[1]);
	printf("Grid dim 3 = %d \n",devp.maxGridSize[2]);

	if(all_info == true) {

	   printf("\n");
	   printf("CUDA Device texturing characteristics: \n");
	   printf("Device max surface 1D:              %d     \n",devp.maxSurface1D);
	   printf("Device max surface layered 1D:      %d     \n",devp.maxSurface1DLayered);
	   printf("Device max surface 2D:              %d     \n",devp.maxSurface2D);
	   printf("Device max surface layered 2D:      %d     \n",devp.maxSurface2DLayered);
	   printf("Device max surface 3D:              %d     \n",devp.maxSurface3D);
	   printf("Device max surface Cubemap:         %d     \n",devp.maxSurfaceCubemap);
	   printf("Device max surface Cubemap layered: %d     \n",devp.maxSurfaceCubemapLayered);
	   printf("Device max texture 1D:			   %d     \n",devp.maxTexture1D);
	   printf("Device max texture 1D layered:      %d     \n",devp.maxTexture1DLayered);
	   printf("Device max texture 1D linear:       %d     \n",devp.maxTexture1DLinear);
	   printf("Device max texture 1D mipmap:       %d     \n",devp.maxTexture1DMipmap);
	   printf("Device max texture 2D:              %d     \n",devp.maxTexture2D);
	   printf("Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	   printf("Device max texture 2D Linear:       %d     \n",devp.maxTexture2DLinear);
	   printf("Device max texture 2D Layered:      %d     \n",devp.maxTexture2DLayered);
	   printf("Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	   printf("Device max texture 2D Mipmap:       %d     \n",devp.maxTexture2DMipmap);
	   printf("Device max texture 3D:			   %d     \n",devp.maxTexture3D);
	   printf("Device max texture 3D Alt:          %d     \n",devp.maxTexture3DAlt);
	}
	printf("CUDA Device name: %s -- end of capabilities dump!! \n",devp.name);
	*ierr = 0;
}

void devcaps_to_file(const int dev_num, const char * fname, 
		     int * ierr,     const bool all_info) {



#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of devcaps_to_file");
#endif
	cudaDeviceProp devp;
	cudaError status;
#if (GMS_CUDA_DEBUG_ON) == 1
	 GMS_CUDA_DEBUG_CHECK(cudaGetDeviceProperties(&devp,dev_num));
#else
	 GMS_CUDA_CHECK(cudaGetDeviceProperties(&devp,dev_num));
#endif
	// Check CUDA device no. x is installed
	if(devp.major == 9999 && devp.minor == 9999) {
		fprintf(stderr,"cudaGetDeviceProperties find out no CUDA device(s) ... terminating!!\n");
		exit_fatal("No CUDA device has been found",
			      "CuWRF/Run-time errors.txt"  );
	}
	FILE *fp = NULL;
	fp = fopen(fname,"a+");
	if(NULL != fp) {
	   fprintf(fp,"Begin dump of CUDA device no: %d\n", dev_num);
	   fprintf(fp,"Device ASCII name:			   %s       \n",devp.name);
	   fprintf(fp,"Device major compute caps:      %d       \n",devp.major);
	   fprintf(fp,"Device minor compute caps:      %d       \n",devp.minor);
	   fprintf(fp,"\n");
	   fprintf(fp,"CUDA Device HW characteristics: \n");
	   fprintf(fp,"Device clock rate:             %d Khz  .\n",devp.clockRate);
	   fprintf(fp,"Device async engine count:     %d       \n",devp.asyncEngineCount);
	   fprintf(fp,"Device L2 cache size:	       %d bytes.\n",devp.l2CacheSize);
	   fprintf(fp,"Device multi-processor count:  %d       \n",devp.multiProcessorCount);
	   fprintf(fp,"Device concurrent kernels:     %d       \n",devp.concurrentKernels);
	   fprintf(fp,"Device compute mode:           %d       \n",devp.computeMode);
	   fprintf(fp,"Device kernel timeout enabled: %d       \n",devp.kernelExecTimeoutEnabled);
	   fprintf(fp,"Device PCI Device ID:          %d       \n",devp.pciDeviceID);
	   fprintf(fp,"Device PCI Domain ID:		  %d       \n",devp.pciDomainID);
	   fprintf(fp,"Device PCI Bus ID:             %d       \n",devp.pciBusID);
	   fprintf(fp,"Device ECC enabled:		       %d       \n",devp.ECCEnabled);
	   fprintf(fp,"Device host mem mapping:       %d       \n",devp.canMapHostMemory);
	   fprintf(fp,"Device registers per block:    %d       \n",devp.regsPerBlock);
	   fprintf(fp,"Device warp size:              %d       \n",devp.warpSize);
	   fprintf(fp,"Device overlap:                %d       \n",devp.deviceOverlap);
	   fprintf(fp,"\n");
	   fprintf(fp,"CUDA Device memory characteristics: \n");
	   fprintf(fp,"Device memory bus width:       %d       \n",devp.memoryBusWidth);
	   fprintf(fp,"Device memory clock rate:      %d       \n",devp.memoryClockRate);
	   fprintf(fp,"Device memory pitch:           %d bytes.\n",devp.memPitch);
	   fprintf(fp,"Device shared memory/block:    %d bytes.\n",devp.sharedMemPerBlock);
	   fprintf(fp,"Device total const memory:     %.9f KiB.\n",TO_KiB(devp.totalConstMem));
	   fprintf(fp,"Device total global memory:    %.9f KiB.\n",TO_KiB(devp.totalGlobalMem));
	   fprintf(fp,"Device unified addressing:     %d       \n",devp.unifiedAddressing);
	   fprintf(fp,"\n");
	   fprintf(fp,"CUDA Device compute characteristics: \n");
	   fprintf(fp,"Device max threads/block:      %d    \n",devp.maxThreadsPerBlock);
	   fprintf(fp,"Device max threads/dimension:        \n");
	   fprintf(fp,"Dim 1 = %d threads.\n",devp.maxThreadsDim[0]);
	   fprintf(fp,"Dim 2 = %d threads.\n",devp.maxThreadsDim[1]);
	   fprintf(fp,"Dim 3 = %d threads.\n",devp.maxThreadsDim[2]);
	   fprintf(fp,"Device max thread/multiproc:   %d    \n",devp.maxThreadsPerMultiProcessor);
	   fprintf(fp,"Device max size of grid:             \n");
	   fprintf(fp,"Grid dim 1 = %d \n",devp.maxGridSize[0]);
	   fprintf(fp,"Grid dim 2 = %d \n",devp.maxGridSize[1]);
	   fprintf(fp,"Grid dim 3 = %d \n",devp.maxGridSize[2]);

	   if(all_info == true) {

	      fprintf(fp,"\n");
	      fprintf(fp,"CUDA Device texturing characteristics: \n");
	      fprintf(fp,"Device max surface 1D:              %d     \n",devp.maxSurface1D);
	      fprintf(fp,"Device max surface layered 1D:      %d     \n",devp.maxSurface1DLayered);
	      fprintf(fp,"Device max surface 2D:              %d     \n",devp.maxSurface2D);
	      fprintf(fp,"Device max surface layered 2D:      %d     \n",devp.maxSurface2DLayered);
	      fprintf(fp,"Device max surface 3D:              %d     \n",devp.maxSurface3D);
	      fprintf(fp,"Device max surface Cubemap:         %d     \n",devp.maxSurfaceCubemap);
	      fprintf(fp,"Device max surface Cubemap layered: %d     \n",devp.maxSurfaceCubemapLayered);
	      fprintf(fp,"Device max texture 1D:			  %d     \n",devp.maxTexture1D);
	      fprintf(fp,"Device max texture 1D layered:      %d     \n",devp.maxTexture1DLayered);
	      fprintf(fp,"Device max texture 1D linear:       %d     \n",devp.maxTexture1DLinear);
	      fprintf(fp,"Device max texture 1D mipmap:       %d     \n",devp.maxTexture1DMipmap);
	      fprintf(fp,"Device max texture 2D:              %d     \n",devp.maxTexture2D);
	      fprintf(fp,"Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	      fprintf(fp,"Device max texture 2D Linear:       %d     \n",devp.maxTexture2DLinear);
	      fprintf(fp,"Device max texture 2D Layered:      %d     \n",devp.maxTexture2DLayered);
	      fprintf(fp,"Device max texture 2D Gather:       %d     \n",devp.maxTexture2DGather);
	      fprintf(fp,"Device max texture 2D Mipmap:       %d     \n",devp.maxTexture2DMipmap);
	      fprintf(fp,"Device max texture 3D:			  %d     \n",devp.maxTexture3D);
	      fprintf(fp,"Device max texture 3D Alt:          %d     \n",devp.maxTexture3DAlt);
	   }
	   fprintf(fp,"CUDA Device name: %s -- end of capabilities dump!! \n",devp.name);
	   *ierr = 0;
	}
	else {
		fprintf(stderr,"fopen returned invalid pointer: %p \n",fp);
		*ierr = -2;
		return;
	}
Error:
	 *ierr = -1;
	 return;
}

#if !defined (CHECK_FAILURE_GENERATOR)
#define CHECK_FAILURE_GENERATOR(ierr,msg) \
 do {                                  \
     if((ierr) < 0) {                     \
         fprintf(stderr,"%s %d\n",(msg),(ierr));      \
         return;                        \
	 }                                  \
 }while(0);
#endif

#if !defined (MALLOC_FAILED_GENERATOR)
#define MALLOC_FAILED_GENERATOR(ptr,msg) \
	do {                              \
		 if(NULL == ptr) {             \
		    fprintf(stderr, "%s %p\n",(msg),(ptr)); \
			return;								 \
		 }										 \
	}while(0);
#endif

#if !defined (INIT_VECI32_GENERATOR) 
#define INIT_VECI32_GENERATOR     \
	for(int32_t i = 0; i != len; ++i) {  \
        h_i32a[i] = i;                \
        h_i32b[i] = i;                \
		h_i32c[i] = 0;                \
	}                                  
#endif

#if !defined (INIT_VECR4_GENERATOR)
#define INIT_VECR4_GENEREATOR     \
	for(int32_t i = 0; i != len; ++i) { \
		h_r4a[i] = (float)i;       \
		h_r4b[i] = (float)i;       \
		h_r4c[i] = 0.F;            \
	}
#endif

#if !defined (INIT_VECR8_GENERATOR)
#define INIT_VECR8_GENERATOR      \
	for(int32_t i = 0; i != len; ++i) { \
		h_r8a[i] = (double)i;       \
		h_r8b[i] = (double)i;       \
		h_r8c[i] = 0.0;            \
	}
#endif

// Random value initialization

#if !defined (INIT_RAND_VECI4_GENERATOR)
#define INIT_RAND_VECI4_GENERATOR             \
	for(int32_t i = 0; i != len; ++i) {       \
	   uint32_t seed = 0;			           \
	   int32_t ret = _rdrand32_step(&seed);   \
	   if(ret == 1){                        \
	      h_i32a[i] = seed;                  \
		  h_i32b[i] = seed;                  \
	   }                                    \
	   else{                                \
	        h_i32a[i] = i;				    \
			h_i32b[i] = i;				    \
	   }								        \
	   h_i32c[i] = 0;                        \
    }
#endif
		
#if !defined (INIT_RAND_VECR4_GENERATOR)
#define INIT_RAND_VECR4_GENERATOR            \
	for(int32_t i = 0; i != len; ++i) {      \
		uint32_t seed = 0;                  \
		int32_t ret = _rdrand32_step(&seed);  \
		if(1 == ret) {                      \
		   h_r4a[i] = (float)seed;           \
		   h_r4b[i] = (float)seed;           \
		}                                  \
		else {                             \
		    h_r4a[i] = (float)i;            \
            h_r4b[i] = (float)i;            \
		}                                 \
		h_r4c[i] = 0.F;                    \
	}
#endif

#if !defined (INIT_RAND_VECR8_GENERATOR)
#define INIT_RAND_VECR8_GENERATOR            \
	for(int32_t i = 0; i != len; ++i) {      \
	   uint64_t seed = 0ULL;                \
	   int32_t ret = _rdrand64_step(&seed);  \
	   if(1 == ret) {                      \
          h_r8a[i] = (double)seed;          \
		  h_r8b[i] = (double)seed;          \
	   }									  \
	   else {                             \
	         h_r8a[i] = (double)i;         \
			 h_r8b[i] = (double)i;         \
	   }                                 \
	    h_r8c[i] = 0.0;                   \
	}
#endif

#if !defined (FREE_GPU_MEMORY) 
#define FREE_GPU_MEMORY               \
	 cudaFree(dev_r8c);			    \
	 cudaFree(dev_r8b);				\
	 cudaFree(dev_r8a);				\
	 cudaFree(dev_r4c);		        \
	 cudaFree(dev_r4b);		        \
	 cudaFree(dev_r4a);			    \
	 cudaFree(dev_i32c);              \
	 cudaFree(dev_i32b);              \
	 cudaFree(dev_i32a);              
#endif

#if !defined (FREE_CPU_MEMORY)
#define FREE_CPU_MEMORY               \
     _aligned_free(h_r8c);            \
	 _mm_free(h_r8b);			\
	 _mm_free(h_r8a);		    \
	 _mm_free(h_r4c);			\
	 _mm_free(h_r4b);			\
	 _mm_free(h_r4a);		    \
	 _mm_free(h_i32c);		    \
	 _mm_free(h_i32b);           \
	 _mm_free(h_i32a);			
#endif

#if !defined (PRINT_VECI4)
#define PRINT_VECI4	  \
	for(int32_t i = 0; i != len; ++i)  \
	    printf("GPU computed int32_t vector: %d\n",h_i32c[i]);
#endif

#if !defined (PRINT_VECR4)
#define PRINT_VECR4     \
	for(int32_t i = 0; i != len; ++i) \
		 printf("GPU computed REAL4 vector: %.6f\n",h_r4c[i]);
#endif

#if !defined (PRINT_VECR8)
#define PRINT_VECR8    \
	for(int32_t i = 0; i != len; ++i) \
		printf("GPU computed REAL8 vector: %.16f\n",h_r8c[i]);
#endif






__global__ void kvec_add_int32(int32_t * __restrict c,
						     const int32_t * __restrict b,
							 const int32_t * __restrict a) {
	int i = threadIdx.x;
	c[i] = b[i] + a[i];
}

__global__ void kvec_add_real4(REAL4 * __restrict c,
						    const REAL4 * __restrict b,
							const REAL4 * __restrict a) {
	int i = threadIdx.x;
	c[i] = b[i] + a[i];
}

__global__ void kvec_add_real8(REAL8 * __restrict c,
							const REAL8 * __restrict b,
							const REAL8 * __restrict a) {
	int i = threadIdx.x;
	c[i] = b[i] + a[i];
}

cudaError_t gpu_vec_add_tests(const int dev_num) {


#if (LOG_ACTIVITY) == 1
	GMS_CUDA_LOG("At prolog of gpu_vec_add_tests.\n");
#endif
	cudaError status;
	int32_t * __restrict dev_i32a = NULL;
	int32_t * __restrict dev_i32b = NULL;
	int32_t * __restrict dev_i32c = NULL;
	float   * __restrict dev_r4a  = NULL;
	float   * __restrict dev_r4b  = NULL;
	float   * __restrict dev_r4c  = NULL;
	double  * __restrict dev_r8a  = NULL;
	double  * __restrict dev_r8b  = NULL;
	double  * __restrict dev_r8c  = NULL;
	// Host arrays.
	int32_t * __restrict h_i32a = NULL;
	int32_t * __restrict h_i32b = NULL;
	int32_t * __restrict h_i32c = NULL;
	float *   __restrict h_r4a  = NULL;
	float *   __restrict h_r4b  = NULL;
	float *   __restrict h_r4c  = NULL;
	double *   __restrict h_r8a  = NULL;
	double *   __restrict h_r8b  = NULL;
	double *   __restrict h_r8c  = NULL;
	int32_t  ierr = 0;
	const int32_t len = 1 << 16;

	status = cudaSetDevice(dev_num);
	if(status != cudaSuccess) {
		fprintf(stderr,"CUDA Device: -- FATAL-ERROR -- : [%s]\n",cudaGetErrorString(status));
		return;
	}
	
	alloc1D_int32_gpu(&dev_i32a[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_int32_gpu failed allocating: dev_i32a")

	alloc1D_int32_gpu(&dev_i32b[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_int32_gpu failed allocating: dev_i32b")

	alloc1D_int32_gpu(&dev_i32c[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_int32_gpu failed allocating: dev_i32c")

	alloc1D_real4_gpu(&dev_r4a[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_real4_gpu failed allocating: dev_r4a")

	alloc1D_real4_gpu(&dev_r4b[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_real4_gpu failed allocating: dev_r4b")

	alloc1D_real4_gpu(&dev_r4c[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_real4_gpu failed allocating: dev_r4c")

	alloc1D_real8_gpu(&dev_r8a[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_real8_gpu failed allocating: dev_r8a")

	alloc1D_real8_gpu(&dev_r8b[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_real8_gpu failed allocating: dev_r8b")

	alloc1D_real8_gpu(&dev_r8c[0],len,&ierr);
	CHECK_FAILURE_GENERATOR(ierr,"alloc1D_real8_gpu failed allocating: dev_r8c")

    // Allocate host arrays.
	h_i32a = (int32_t*)_aligned_malloc(((size_t)len) * sizeof(int32_t),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_i32a,"_aligned_malloc failed to allocate: h_i32a")

	h_i32b = (int32_t*)_aligned_malloc(((size_t)len) * sizeof(int32_t),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_i32b,"_aligned_malloc failed to allocate: h_i32b")

	h_i32c = (int32_t*)_aligned_malloc(((size_t)len) * sizeof(int32_t),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_i32c,"_aligned_malloc failed to allocate: h_i32c")

	h_r4a = (REAL4*)_aligned_malloc(((size_t)len) * sizeof(REAL4),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r4a,"_aligned_malloc failed to allocate: h_r4a")

	h_r4b = (REAL4*)_aligned_malloc(((size_t)len) *sizeof(REAL4),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r4b,"_aligned_malloc failed to allocate: h_r4b")

	h_r4c = (REAL4*)_aligned_malloc(((size_t)len) * sizeof(REAL4),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r4c,"_aligned_malloc failed to allocate: h_r4c")

	h_r8a = (REAL8*)_aligned_malloc(((size_t)len) * sizeof(REAL8),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r8a,"_aligned_malloc failed to allocate: h_r8a")

	h_r8b = (REAL8*)_aligned_malloc(((size_t)len) * sizeof(REAL8),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r8b,"_aligned_malloc failed to allocate: h_r8b")

	h_r8c = (REAL8*)_aligned_malloc(((size_t)len) * sizeof(REAL8),HOST_ALIGN64);
	MALLOC_FAILED_GENERATOR(h_r8c,"_aligned_malloc failed to allocate: h_r8c")

	INIT_RAND_VECI4_GENERATOR
	INIT_RAND_VECR4_GENERATOR
	INIT_RAND_VECR8_GENERATOR

	// Copy host arrays to device arrays
	// Do not check for minor errors (non-critical)
	copy1D_int32_cpu_to_gpu(&dev_i32a[0],&h_i32a[0],len,&ierr);
	copy1D_int32_cpu_to_gpu(&dev_i32b[0],&h_i32b[0],len,&ierr);
	copy1D_int32_cpu_to_gpu(&dev_i32c[0],&h_i32c[0],len,&ierr);
	copy1D_real4_cpu_to_gpu(&dev_r4a[0],&h_r4a[0],len,&ierr);
	copy1D_real4_cpu_to_gpu(&dev_r4b[0],&h_r4b[0],len,&ierr);
	copy1D_real4_cpu_to_gpu(&dev_r4c[0],&h_r4c[0],len,&ierr);
	copy1D_real8_cpu_to_gpu(&dev_r8a[0],&h_r8a[0],len,&ierr);
	copy1D_real8_cpu_to_gpu(&dev_r8b[0],&h_r8b[0],len,&ierr);
	copy1D_real8_cpu_to_gpu(&dev_r8c[0],&h_r8c[0],len,&ierr);

	//
	// Launch simple kernels with single thread
	// per array element.
	//

	kvec_add_int32<<<1,len>>>(&dev_i32c[0],&dev_i32b[0],&dev_i32a[0]);
	kvec_add_real4<<<1,len>>>(&dev_r4c[0],&dev_r4b[0],&dev_r4a[0]);
	kvec_add_real8<<<1,len>>>(&dev_r8c[0],&dev_r8b[0],&dev_r8a[0]);
	
	CuWRF_DEBUG_CHECK(cudaGetLastError());
    CuWRF_DEBUG_CHECK(cudaDeviceSynchronize());

	//
	// Copy rresulting device arrays to host memory.
	//Do not check for minor errors (non-critical)
	//
	copy1D_int32_gpu_to_cpu(&dev_i32c[0],&h_i32c[0],len,&ierr);
	copy1D_real4_gpu_to_cpu(&dev_r4c[0],&h_r4c[0],len,&ierr);
	copy1D_real8_gpu_to_cpu(&dev_r8c[0],&h_r8c[0],len,&ierr);

	

	PRINT_VECI4
	PRINT_VECR4
	PRINT_VECR8

	CuWRF_DEBUG_CHECK(cudaDeviceReset());
	
	FREE_CPU_MEMORY
Error:
	
	FREE_GPU_MEMORY
	FREE_CPU_MEMORY
	return;
}
