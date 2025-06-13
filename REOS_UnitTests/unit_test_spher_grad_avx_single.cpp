
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle
#include "GMS_spher_grad_avx.h"

/*
   icpc -o unit_test_spher_grad_avx_single -fp-model fast=2 -ftz -ggdb -ipo  -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_spher_grad_avx.h GMS_spher_grad_avx.cpp unit_test_spher_grad_avx_single.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_spher_grad_avx.h GMS_spher_grad_avx.cpp unit_test_spher_grad_avx_single.cpp
*/


__attribute__((hot))
__attribute__((noinline))
void unit_test_spher_inv_jac_ymm4r8_a();

void unit_test_spher_inv_jac_ymm4r8_a()
{
     __attribute__((aligned(32)))
     double x[4];
     __attribute__((aligned(32)))
     double y[4];
     __attribute__((aligned(32)))
     double z[4];
     __attribute__((aligned(32)))
     double J_out[36];
     double * __restrict__ J0{&J_out[0]};
     double * __restrict__ J1{&J_out[4]};
     double * __restrict__ J2{&J_out[8]};
     double * __restrict__ J3{&J_out[12]};
     double * __restrict__ J4{&J_out[16]};
     double * __restrict__ J5{&J_out[20]};
     double * __restrict__ J6{&J_out[24]};
     double * __restrict__ J7{&J_out[28]};
     double * __restrict__ J8{&J_out[32]};
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_z{0ull};
     std::clock_t seed_distro{0ull};
     int32_t which{-1};
     int32_t sysType{-1};
     seed_distro = std::clock();
     auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,4),
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
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             
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
             auto rand_x{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
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
             auto rand_x{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
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
             auto rand_x{std::bind(std::cauchy_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::cauchy_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::cauchy_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_a(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_load_pd(&x[0]),
                                               _mm256_load_pd(&y[0]),
                                               _mm256_load_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
          }
          break;
          default : 
             printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
             std::terminate();
     }

}

__attribute__((hot))
__attribute__((noinline))
void unit_test_spher_inv_jac_ymm4r8_u();

void unit_test_spher_inv_jac_ymm4r8_u()
{
    
     double x[4];
     double y[4];
     double z[4];
     double J_out[36];
     double * __restrict__ J0{&J_out[0]};
     double * __restrict__ J1{&J_out[4]};
     double * __restrict__ J2{&J_out[8]};
     double * __restrict__ J3{&J_out[12]};
     double * __restrict__ J4{&J_out[16]};
     double * __restrict__ J5{&J_out[20]};
     double * __restrict__ J6{&J_out[24]};
     double * __restrict__ J7{&J_out[28]};
     double * __restrict__ J8{&J_out[32]};
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_z{0ull};
     std::clock_t seed_distro{0ull};
     int32_t which{-1};
     int32_t sysType{-1};
     seed_distro = std::clock();
     auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,4),
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
             auto rand_x{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::uniform_real_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             
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
             auto rand_x{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
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
             auto rand_x{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::normal_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_a().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
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
             auto rand_x{std::bind(std::cauchy_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_x))};
             seed_y = std::clock();
             auto rand_y{std::bind(std::cauchy_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_y))};
             seed_z = std::clock();
             auto rand_z{std::bind(std::cauchy_distribution<double>(0.0,1.57079632679489661923132169164),
                                   std::mt19937(seed_z))};
             for(int32_t __i{0}; __i != 4; ++__i)
             {
                 float rx{rand_x()};
                 x[__i]   = rx;
                 float ry{rand_y()};
                 y[__i]   = ry;
                 float rz{rand_z()};
                 z[__i]   = rz;
             }
             sysType = 0;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_a().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 1;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 2;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
             sysType = 3;
             printf("[UNIT-TEST]: -- START: spher_inv_jac_ymm4r8_u().\n");
             gms::math::spher_inv_jac_ymm4r8_u(J0,J1,J2,J3,
                                               J4,J5,J6,J7,
                                               J8,
                                               _mm256_loadu_pd(&x[0]),
                                               _mm256_loadu_pd(&y[0]),
                                               _mm256_loadu_pd(&z[0]),
                                               sysType);
             printf("[UNIT-TEST:] -- Dumping the results for sysType=%d\n",sysType);
             for(int32_t __i{0}; __i != 36; ++__i)
             {
                 printf("J_out[%d]=%.16f\n",__i, J_out[__i]);
             }
             printf("[UNIT-TEST]: -- END: spher_inv_jac_ymm4r8_u().\n");
          }
          break;
          default : 
             printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
             std::terminate();
     }

}



int main()
{
    unit_test_spher_inv_jac_ymm4r8_a();
    unit_test_spher_inv_jac_ymm4r8_u();
    return 0;
}