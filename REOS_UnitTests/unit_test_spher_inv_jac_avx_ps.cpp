

#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle
#include "GMS_spher_grad_avx.h"

/*
   icpc -o unit_test_spher_inv_jac_avx_ps -fp-model fast=2 -ftz -ggdb -ipo  -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_spher_grad_avx.h GMS_spher_grad_avx.cpp unit_test_spher_inv_jac_avx_ps.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_spher_grad_avx.h GMS_spher_grad_avx.cpp unit_test_spher_inv_jac_avx_ps.cpp
*/


__attribute__((hot))
__attribute__((noinline))
void unit_test_spher_inv_jac_ymm8r4_a();

void unit_test_spher_inv_jac_ymm8r4_a()
{
     __attribute__((aligned(32)))
     float x[8];
     __attribute__((aligned(32)))
     float y[8];
     __attribute__((aligned(32)))
     float z[8];
     __attribute__((aligned(32)))
     float J_out[72];
     float * __restrict__ J0{&J_out[0]};
     float * __restrict__ J1{&J_out[8]};
     float * __restrict__ J2{&J_out[16]};
     float * __restrict__ J3{&J_out[24]};
     float * __restrict__ J4{&J_out[32]};
     float * __restrict__ J5{&J_out[40]};
     float * __restrict__ J6{&J_out[48]};
     float * __restrict__ J7{&J_out[56]};
     float * __restrict__ J8{&J_out[64]};
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
             std::uniform_real_distribution<float> uniform_distro;
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
             auto rand_x{std::bind(std::uniform_real_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             
          }
          break;
          case 1 : 
          {
             std::normal_distribution<float> normal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(normal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
          }
          break;
          case 2 : 
          {
             std::lognormal_distribution<float> lognormal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(lognormal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognormal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
          }
          case 3 : 
          {
             std::cauchy_distribution<float> cauchy_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(cauchy_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cauchy_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::cauchy_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::cauchy_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::cauchy_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_ps(&x[0]),
                                               _mm256_load_ps(&y[0]),
                                               _mm256_load_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
          }
          break;
          default : 
             printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
             return;
     }

}

__attribute__((hot))
__attribute__((noinline))
void unit_test_spher_inv_jac_ymm8r4_u();

void unit_test_spher_inv_jac_ymm8r4_u()
{
    
     float x[4];
     float y[4];
     float z[4];
     float J_out[36];
     float * __restrict__ J0{&J_out[0]};
     float * __restrict__ J1{&J_out[8]};
     float * __restrict__ J2{&J_out[16]};
     float * __restrict__ J3{&J_out[24]};
     float * __restrict__ J4{&J_out[32]};
     float * __restrict__ J5{&J_out[40]};
     float * __restrict__ J6{&J_out[48]};
     float * __restrict__ J7{&J_out[56]};
     float * __restrict__ J8{&J_out[64]};
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
             std::uniform_real_distribution<float> uniform_distro;
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
             auto rand_x{std::bind(std::uniform_real_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             
          }
          break;
          case 1 : 
          {
             std::normal_distribution<float> normal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(normal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
          }
          break;
          case 2 : 
          {
             std::lognormal_distribution<float> lognormal_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(lognormal_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognormal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
          }
          case 3 : 
          {
             std::cauchy_distribution<float> cauchy_distro;
             char * distro_name{NULL};
             int32_t status{9999};
             distro_name = gms::common::demangle(typeid(cauchy_distro).name(),status);
             if(distro_name != NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cauchy_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             seed_x = std::clock();
             auto rand_x{std::bind(std::cauchy_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::cauchy_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::cauchy_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 8; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_a().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm8r4_u().\n");
             gms::math::spher_inv_jac_ymm8r4_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_ps(&x[0]),
                                               _mm256_loadu_ps(&y[0]),
                                               _mm256_loadu_ps(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 72; ++__i)
             {
                 printf("J_out[%d]=%.7f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm8r4_u().\n");
          }
          break;
          default : 
             printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
             return;
     }

}



int main()
{
    unit_test_spher_inv_jac_ymm8r4_a();
    unit_test_spher_inv_jac_ymm8r4_u();
    return 0;
}