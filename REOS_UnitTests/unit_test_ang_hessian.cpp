

#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle
#include "GMS_spher_grad_avx.h"

/*
   icpc -o unit_test_ang_hessian -fp-model fast=2 -ftz -ggdb -ipo  -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_spher_grad_avx.h GMS_spher_grad_avx.cpp unit_test_ang_hessian.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_spher_grad_avx.h GMS_spher_grad_avx.cpp unit_test_ang_hessian.cpp
*/

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_ang_hessian_ymm4r8();

void unit_test_ang_hessian_ymm4r8()
{
    constexpr int32_t nelems{72};
    __attribute__((aligned(32)))
    __m256d H[18];
    __m256d x;
    __m256d y;
    __m256d z;
    double * __restrict__ pH{nullptr};
    std::clock_t seed_x{0ull};
    std::clock_t seed_y{0ull};
    std::clock_t seed_z{0ull};
    std::clock_t seed_distro{0ull};
    int32_t which{-1};
    int32_t sysType{-1};
    seed_distro = std::clock();
    auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,3),
                         std::mt19937(seed_distro))};
    which = rand();

    switch (which)
    {
         case 0 : 
         {
             std::uniform_real_distribution<double> uniform_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(uniform_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             x = _mm256_setr_pd(rand_x(),rand_x(),rand_x(),rand_x());
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,6.283185307179586476925286766559),
                                   std::mt19937(seed_y))};
             y = _mm256_setr_pd(rand_y(),rand_y(),rand_y(),rand_y());
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<double>(0.0,50000),
                                   std::mt19937(seed_z))};
             z = _mm256_setr_pd(rand_z(),rand_z(),rand_z(),rand_z());
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_ang_hess_ymm4r8().\n");
             gms::math::spher_ang_hess_ymm4r8(H,x,y,z,sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             pH = (double* __restrict)&H[0];
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("H[%d]=%.16f\n",__i, pH[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_ang_hess_ymm4r8().\n");
         } 
         break;
         case 1 :
         {
             std::uniform_real_distribution<double> uniform_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(uniform_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             x = _mm256_setr_pd(rand_x(),rand_x(),rand_x(),rand_x());
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,6.283185307179586476925286766559),
                                   std::mt19937(seed_y))};
             y = _mm256_setr_pd(rand_y(),rand_y(),rand_y(),rand_y());
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<double>(0.0,50000),
                                   std::mt19937(seed_z))};
             z = _mm256_setr_pd(rand_z(),rand_z(),rand_z(),rand_z());
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_ang_hess_ymm4r8().\n");
             gms::math::spher_ang_hess_ymm4r8(H,x,y,z,sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             pH = (double* __restrict)&H[0];
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("H[%d]=%.16f\n",__i, pH[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_ang_hess_ymm4r8().\n");
         }
         break;
         case 2 : 
         {
             std::uniform_real_distribution<double> uniform_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(uniform_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             x = _mm256_setr_pd(rand_x(),rand_x(),rand_x(),rand_x());
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,6.283185307179586476925286766559),
                                   std::mt19937(seed_y))};
             y = _mm256_setr_pd(rand_y(),rand_y(),rand_y(),rand_y());
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<double>(0.0,50000),
                                   std::mt19937(seed_z))};
             z = _mm256_setr_pd(rand_z(),rand_z(),rand_z(),rand_z());
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_ang_hess_ymm4r8().\n");
             gms::math::spher_ang_hess_ymm4r8(H,x,y,z,sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             pH = (double* __restrict)&H[0];
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("H[%d]=%.16f\n",__i, pH[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_ang_hess_ymm4r8().\n");
         }
         break;
         case 3 :
         {
            std::uniform_real_distribution<double> uniform_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(uniform_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             x = _mm256_setr_pd(rand_x(),rand_x(),rand_x(),rand_x());
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,6.283185307179586476925286766559),
                                   std::mt19937(seed_y))};
             y = _mm256_setr_pd(rand_y(),rand_y(),rand_y(),rand_y());
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<double>(0.0,50000),
                                   std::mt19937(seed_z))};
             z = _mm256_setr_pd(rand_z(),rand_z(),rand_z(),rand_z());
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_ang_hess_ymm4r8().\n");
             gms::math::spher_ang_hess_ymm4r8(H,x,y,z,sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             pH = (double* __restrict)&H[0];
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("H[%d]=%.16f\n",__i, pH[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_ang_hess_ymm4r8().\n");
         }
         break;
         default : 
              return;
    }
}

int main()
{
    unit_test_ang_hessian_ymm4r8();
    return 0;
}