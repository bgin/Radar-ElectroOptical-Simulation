
#include <malloc.h> // for _mm_malloc

#include "GMS_kmeans.h"
#include "GMS_kmeans_wrapper.h"
#if defined _WIN64
    #include "../GMS_config.h"
    #include "../GMS_error_macros.h"
    #if (PRINT_CALLSTACK_ON_ERROR) == 1
       #include "../System/StackWalker.h"
    #endif
    #include "../GMS_indices.h"
#elif defined __linux
    #include "GMS_config.h"
    #include "GMS_error_macros.h"
    #include "GMS_indices.h"
#endif
#include "GMS_fpexceptions.h"

#if !defined (GMS_KMEANS_FP_EXCEPT_PROLOG)
#define GMS_KMEANS_FP_EXCEPT_PROLOG(status)  \
	(status) = clear_fpexcepts();           \
    if (0 != (status))                     \
	    std::cerr << "clear_fpexcpets failed to clear FP environment!!!  Detection of floating-point exceptions will not be possible!!!\n"; 
#endif

#if !defined (GMS_KMEANS_FP_EXCEPT_EPILOG)
#define GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)          \
  do {		                                                    \
									                            \
		if (0 != (status)) {                                    \
		  (fpexcepts[0]) = test_feinvalid(FE_INVALID);		    \
		  (fpexcepts[1]) = test_feinexact(FE_INEXACT);           \
		  (fpexcepts[2]) = test_fedivbyzero(FE_DIVBYZERO);       \
		  (fpexcepts[3]) = test_feunormal(FE_DENORMAL);          \
		  (fpexcepts[4]) = test_feoverflow(FE_OVERFLOW);         \
		  (fpexcepts[5]) = test_feunderflow(FE_UNDERFLOW);       \
		}											             \
  } while (0); 
#endif			 		  


gms::math::stat
::KMEANS::KMEANS()
:
m_dim_num{},
m_samples_num{},
m_cluster_num{},
m_iter_max{},
m_iter_num{},
m_seed{},
m_cluster(NULL),
m_cluster_center(NULL),
m_cluster_population(NULL),
m_cluster_energy(NULL) {}

gms::math::stat
::KMEANS::KMEANS(
const int32_t dim_num,
const int32_t samples_num,
const int32_t cluster_num,
const int32_t iter_max,
const int32_t seed)
:
m_dim_num{ dim_num },
m_samples_num{ samples_num },
m_cluster_num{ cluster_num },
m_iter_max{ iter_max },
m_iter_num{},
m_seed{ seed },
m_cluster(reinterpret_cast<int32_t*>(_mm_malloc(m_samples_num*sizeof(int32_t), align64B))),
m_cluster_center(reinterpret_cast<double*>(_mm_malloc((m_dim_num*m_cluster_num)*sizeof(double), align64B))),
m_cluster_population(reinterpret_cast<int32_t*>(_mm_malloc(m_cluster_num*sizeof(int32_t), align64B))),
m_cluster_energy(reinterpret_cast<double*>(_mm_malloc(m_cluster_num*sizeof(double), align64B))) {
	if (NULL == m_cluster ||
	    NULL == m_cluster_center ||
            NULL == m_cluster_population ||
	    NULL == m_cluster_energy) {
#if defined _WIN64
    #if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
    #endif
#endif		  
			ABORT_ON_ERROR("KMEANS::KMEANS(int32_t,int32_t,int32_t,int32_t,int32_t) -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_samples_num; ++i)
		m_cluster[i] = 0;

	for (int32_t i = 0; i != m_dim_num; ++i)
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t j = 0; j != m_cluster_num; ++j)
		m_cluster_center[Ix2D(i, m_dim_num, j)] = 0.0;

#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_cluster_num; ++i) {
		m_cluster_population[i] = 0;
		m_cluster_energy[i] = 0.0;
	}
}

gms::math::stat
::KMEANS::KMEANS(const KMEANS &x)
:
m_dim_num{ x.m_dim_num },
m_samples_num{ x.m_samples_num },
m_cluster_num{ x.m_cluster_num },
m_iter_max{ x.m_iter_max },
m_iter_num{ x.m_iter_num },
m_seed{ x.m_seed },
m_cluster(reinterpret_cast<int32_t*>(_mm_malloc(m_samples_num*sizeof(int32_t), align64B))),
m_cluster_center(reinterpret_cast<double*>(_mm_malloc((m_dim_num*m_cluster_num)*sizeof(double), align64B))),
m_cluster_population(reinterpret_cast<int32_t*>(_mm_malloc(m_cluster_num*sizeof(int32_t), align64B))),
m_cluster_energy(reinterpret_cast<double*>(_mm_malloc(m_cluster_num*sizeof(double), align64B))) {
	if (NULL == m_cluster ||
	    NULL == m_cluster_center ||
	    NULL == m_cluster_population ||
	    NULL == m_cluster_energy) {
#if defined _WIN64
    #if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
    #endif
#endif		  
			ABORT_ON_ERROR("KMEANS::KMEANS(KMEANS &) -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_samples_num; ++i)
		m_cluster[i] = x.m_cluster[i];

	for (int32_t i = 0; i != m_dim_num; ++i)
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t j = 0; j != m_cluster_num; ++j)
		m_cluster_center[Ix2D(i, m_cluster_num, j)] = x.m_cluster_center[Ix2D(i, m_cluster_num, j)];

#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_cluster_num; ++i) {
		m_cluster_population[i] = x.m_cluster_population[i];
		m_cluster_energy[i] = x.m_cluster_energy[i];
	}
}

gms::math::stat
::KMEANS::KMEANS(KMEANS &&x)
:
m_dim_num{ x.m_dim_num },
m_samples_num{ x.m_samples_num },
m_cluster_num{ x.m_cluster_num },
m_iter_max{ x.m_iter_max },
m_iter_num{ x.m_iter_num },
m_seed{ x.m_seed },
m_cluster(NULL),
m_cluster_center(NULL),
m_cluster_population(NULL),
m_cluster_energy(NULL) {
	m_cluster = &x.m_cluster[0];
	m_cluster_center = &x.m_cluster_center[0];
	m_cluster_population = &x.m_cluster_population[0];
	m_cluster_energy = &x.m_cluster_energy[0];

	x.m_cluster_num = 0;
	x.m_dim_num = 0;
	x.m_samples_num = 0;
	x.m_iter_max = 0;
	x.m_iter_num = 0;
	x.m_seed = 0;
	x.m_cluster = NULL;
	x.m_cluster_center = NULL;
	x.m_cluster_population = NULL;
	x.m_cluster_energy = NULL;
}

gms::math::stat
::KMEANS::~KMEANS() {
	if (NULL != m_cluster_energy) {
		_mm_free(m_cluster_energy);
		m_cluster_energy = NULL;
	}
	if (NULL != m_cluster_population) {
		_mm_free(m_cluster_population);
		m_cluster_population = NULL;
	}
	if (NULL != m_cluster_center) {
		_mm_free(m_cluster_center);
		m_cluster_center = NULL;
	}
	if (NULL != m_cluster) {
		_mm_free(m_cluster);
		m_cluster = NULL;
	}
}

gms::math::stat::KMEANS &
gms::math::stat::KMEANS::operator=(const KMEANS &x) {
	if (this == &x) return (*this);
	m_dim_num = x.m_dim_num;
	m_samples_num = x.m_samples_num;
	m_cluster_num = x.m_cluster_num;
	m_iter_max = x.m_iter_max;
	m_iter_num = x.m_iter_num;
	m_seed = x.m_seed;
	_mm_free(m_cluster_energy);
	_mm_free(m_cluster_population);
	_mm_free(m_cluster_center);
	_mm_free(m_cluster);
	int32_t * __restrict cluster =
		reinterpret_cast<int32_t*>(_mm_malloc(m_samples_num*sizeof(int32_t), align64B));
	double * __restrict cluster_center =
		reinterpret_cast<double*>(_mm_malloc((m_dim_num*m_cluster_num)*sizeof(double), align64B));
	int32_t * __restrict cluster_population =
		reinterpret_cast<int32_t*>(_mm_malloc(m_cluster_num*sizeof(int32_t), align64B));
	double * __restrict cluster_energy =
		reinterpret_cast<double*>(_mm_malloc(m_cluster_num*sizeof(double), align64B));
	if (    NULL == cluster ||
		NULL == cluster_center ||
		NULL == cluster_population ||
		NULL == cluster_energy) {
#if defined _WIN64	  
    #if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
    #endif
#endif		  
			ABORT_ON_ERROR("HMEANS02::operator=(HMEANS &) -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != x.m_samples_num; ++i)
		cluster[i] = x.m_cluster[i];

	for (int32_t i = 0; i != x.m_dim_num; ++i)
#pragma simd vectorlengthfor(double)  vecremainder
	for (int32_t j = 0; j != x.m_cluster_num; ++j)
		cluster_center[Ix2D(i, x.m_dim_num, j)] = x.m_cluster_center[Ix2D(i, x.m_dim_num, j)];
#pragma simd vectorlengthfor(double)  vecremainder
	for (int32_t i = 0; i != x.m_cluster_num; ++i) {
		cluster_population[i] = x.m_cluster_population[i];
		cluster_energy[i] = x.m_cluster_energy[i];
	}
	m_cluster = cluster;
	m_cluster_center = cluster_center;
	m_cluster_population = cluster_population;
	m_cluster_energy = cluster_energy;
	return (*this);
}

gms::math::stat::KMEANS &
gms::math::stat::KMEANS::operator=(KMEANS &&x) {
	if (this == &x) return (*this);
	m_cluster_num = x.m_cluster_num;
	m_dim_num = x.m_dim_num;
	m_samples_num = x.m_samples_num;
	m_iter_max = x.m_iter_max;
	m_iter_num = x.m_iter_num;
	m_seed = x.m_seed;
	_mm_free(m_cluster_energy);
	_mm_free(m_cluster_population);
	_mm_free(m_cluster_center);
	_mm_free(m_cluster);
	m_cluster = &x.m_cluster[0];
	m_cluster_center = &x.m_cluster_center[0];
	m_cluster_population = &x.m_cluster_population[0];
	m_cluster_energy = &x.m_cluster_energy[0];

	x.m_cluster_num = 0;
	x.m_dim_num = 0;
	x.m_samples_num = 0;
	x.m_iter_max = 0;
	x.m_iter_num = 0;
	x.m_seed = 0;
	x.m_cluster = NULL;
	x.m_cluster_center = NULL;
	x.m_cluster_population = NULL;
	x.m_cluster_energy = NULL;
	return (*this);
}


void
gms::math::stat::
compute_hmeans01(KMEANS &hm01,
		 double * __restrict sig_samples,
		 double * __restrict cluster_variance,
		 std::int32_t fpexcepts[6] ) {

	int32_t status{};
	
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)
	cluster_initialize_5(&hm01.m_dim_num, &hm01.m_samples_num, &hm01.m_cluster_num,
		&sig_samples[0], &hm01.m_seed, &hm01.m_cluster_center[0]);

	hmeans_01(&hm01.m_dim_num, &hm01.m_samples_num, &hm01.m_cluster_num, &hm01.m_iter_max,
		&hm01.m_iter_num, &sig_samples[0], &hm01.m_cluster[0], &hm01.m_cluster_center[0],
		&hm01.m_cluster_population[0], &hm01.m_cluster_energy[0]);

	cluster_variance_compute(&hm01.m_dim_num, &hm01.m_samples_num, &hm01.m_cluster_num,
		&sig_samples[0], &hm01.m_cluster[0], &hm01.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&hm01.m_samples_num, &hm01.m_cluster_num, &hm01.m_cluster_population[0],
		&hm01.m_cluster_energy[0], &cluster_variance[0]);
	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::compute_hmeans02(KMEANS &hm02,
						         double * __restrict sig_samples,
							  double * __restrict cluster_variance,
						          std::int32_t fpexcept[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)
	cluster_initialize_5(&hm02.m_dim_num, &hm02.m_samples_num, &hm02.m_cluster_num,
		&sig_samples[0], &hm02.m_seed, &hm02.m_cluster_center[0]);

	hmeans_02(&hm02.m_dim_num, &hm02.m_samples_num, &hm02.m_cluster_num, &hm02.m_iter_max,
		&hm02.m_iter_num, &sig_samples[0], &hm02.m_cluster[0], &hm02.m_cluster_center[0],
		&hm02.m_cluster_population[0], &hm02.m_cluster_energy[0], &hm02.m_seed);

	cluster_variance_compute(&hm02.m_dim_num, &hm02.m_samples_num, &hm02.m_cluster_num,
		&sig_samples[0], &hm02.m_cluster[0], &hm02.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&hm02.m_samples_num, &hm02.m_cluster_num, &hm02.m_cluster_population[0],
		&hm02.m_cluster_energy[0], &cluster_variance[0]);
	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcept)
}

void
gms::math::stat::compute_kmeans01(       KMEANS &km01,
					 double * __restrict sig_samples,
					 double * __restrict cluster_variance,
				         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)
	cluster_initialize_5(&km01.m_dim_num,&km01.m_samples_num,&km01.m_cluster_num,
		&sig_samples[0], &km01.m_seed, &km01.m_cluster_center[0]);

	kmeans_01(&km01.m_dim_num,&km01.m_samples_num,&km01.m_cluster_num,&km01.m_iter_max,
		&km01.m_iter_num, &sig_samples[0], &km01.m_cluster[0], &km01.m_cluster_center[0],
		&km01.m_cluster_population[0], &km01.m_cluster_energy[0]);

	cluster_variance_compute(&km01.m_dim_num, &km01.m_cluster_num, &km01.m_cluster_num,
							 &sig_samples[0],&km01.m_cluster[0],&km01.m_cluster_center[0],
							 &cluster_variance[0] );

	cluster_print_summary(&km01.m_samples_num, &km01.m_cluster_num, &km01.m_cluster_population[0],
		&km01.m_cluster_energy[0], &cluster_variance[0]);
	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::compute_kmeans02(       KMEANS &km02,
				         double * __restrict sig_samples,
				         double * __restrict cluster_variance,
				         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)
	cluster_initialize_1(&km02.m_dim_num,&km02.m_samples_num,&km02.m_cluster_num,
						 &sig_samples[0], &km02.m_cluster_center[0]);

	kmeans_02(&km02.m_dim_num,&km02.m_samples_num,&km02.m_cluster_num,&km02.m_iter_max,
		&km02.m_iter_num, &sig_samples[0], &km02.m_cluster[0], &km02.m_cluster_center[0],
		&km02.m_cluster_population[0], &km02.m_cluster_energy[0]);

	cluster_variance_compute(&km02.m_dim_num,&km02.m_samples_num,&km02.m_cluster_num,
		&sig_samples[0], &km02.m_cluster[0], &km02.m_cluster_center[0],&cluster_variance[0]);

	cluster_print_summary(&km02.m_samples_num, &km02.m_cluster_num, &km02.m_cluster_population[0],
		&km02.m_cluster_energy[0], &cluster_variance[0]);
	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::compute_kmeans03(       KMEANS &km03,
				         double * __restrict sig_samples,
				         double * __restrict cluster_variance,
				         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)
	cluster_initialize_1(&km03.m_dim_num,&km03.m_samples_num,&km03.m_cluster_num,
							 &sig_samples[0],&km03.m_cluster_center[0]);
	
	kmeans_03(&km03.m_dim_num,&km03.m_samples_num,&km03.m_cluster_num,&km03.m_iter_max,
		&km03.m_iter_num, &sig_samples[0], &km03.m_cluster[0], &km03.m_cluster_center[0],
		&km03.m_cluster_population[0], &km03.m_cluster_energy[0]);

	cluster_variance_compute(&km03.m_dim_num,&km03.m_samples_num,&km03.m_cluster_num,
		&sig_samples[0], &km03.m_cluster[0], &km03.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km03.m_samples_num, &km03.m_cluster_num, &km03.m_cluster_population[0],
		&km03.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::
compute_hmeans01_and_kmeans01(           KMEANS &km,
				         double * __restrict sig_samples,
				         double * __restrict cluster_variance,
				         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)

	cluster_initialize_5(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	kmeans_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,&km.m_iter_num,
			  &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
			  &km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::
compute_hmeans01_and_kmeans02(  KMEANS &km,
			        double * __restrict sig_samples,
			        double * __restrict cluster_variance,
			        std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_01(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	kmeans_02(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::
compute_hmeans01_and_kmeans03(  KMEANS &km,
			        double * __restrict sig_samples,
			        double * __restrict cluster_variance,
			        std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_01(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	kmeans_03(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)

}

void
gms::math::stat::
compute_hmeans_w_01(     KMEANS &km,
		         double * __restrict sig_samples,
		         double * __restrict cluster_variance,
		         double * __restrict weight,
		         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)

	cluster_initialize_5(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_w_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &weight[0], &km.m_cluster[0],
		&km.m_cluster_center[0], &km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_dim_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::
compute_hmeans_w_02(     KMEANS &km,
		         double * __restrict sig_samples,
			 double * __restrict cluster_variance,
		         double * weight,
		         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_w_02(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &weight[0], &km.m_cluster[0],
		&km.m_cluster_center[0], &km.m_cluster_population[0], &km.m_cluster_energy[0],&km.m_seed);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_dim_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status,fpexcepts)
}

void
gms::math::stat::
compute_kmeans_w_01(     KMEANS &km,
		         double * __restrict sig_samples,
		         double * __restrict cluster_variance,
		         double * weight,
		         std::int32_t fpexcepts[6]) {

	std::int32_t status{};
	GMS_KMEANS_FP_EXCEPT_PROLOG(status)

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	kmeans_w_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &weight[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_dim_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	GMS_KMEANS_FP_EXCEPT_EPILOG(status, fpexcepts)
}


