#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle
#include "GMS_complex_ymm8r4.h"


/*
   icpc -o unit_test_complex_ymm8r4 -fp-model fast=2 -ftz -ggdb -ipo  -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_complex_ymm8r4.h unit_test_complex_ymm8r4.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_simd_utils.h GMS_complex_ymm8r4.h unit_test_complex_ymm8r4.cpp
*/


__attribute__((hot))
__attribute__((noinline))
void unit_test_cadd_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_cadd_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_y[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: cadd_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      // scale by 10.0 
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cadd_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cadd_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cadd_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                 
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cadd_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cadd_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cadd_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      // scale by 10.0 
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cadd_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cadd_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cadd_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cadd_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cadd_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_csub_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_csub_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_y[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: csub_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                 
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = csub_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: csub_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: csub_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                   
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = csub_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: csub_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: csub_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = csub_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: csub_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: csub_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = csub_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: csub_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_cmul_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_cmul_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_y[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: cmul_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                 
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cmul_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cmul_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cmul_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                   
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cmul_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cmul_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cmul_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cmul_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cmul_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cmul_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cmul_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cmul_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_cdiv_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_cdiv_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_y[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: cdiv_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                 
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cdiv_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                   
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cdiv_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cdiv_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_cdiv_smith_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_cdiv_smith_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_y[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: cdiv_smith_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                 
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_smith_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_smith_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cdiv_smith_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                   
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_smith_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_smith_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cdiv_smith_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_smith_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_smith_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cdiv_smith_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   __m256 vry;
                   __m256 viy;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                      vry = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      viy = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_y[__j].re = vry;
                   buf_y[__j].im = viy;
                   buf_z[__j]    = cdiv_smith_ymm8c4(buf_x[__j],buf_y[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cdiv_smith_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_cabs_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_cabs_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: cabs_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                                                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = cabs_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cabs_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cabs_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                            
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = cabs_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cabs_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cabs_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                                    
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = cabs_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cabs_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: cabs_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                                      
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = cabs_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: cabs_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_carg_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_carg_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: carg_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                                                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = carg_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: carg_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: carg_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                            
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = carg_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: carg_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: carg_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                                    
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = carg_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: carg_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: carg_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                                      
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j].re = carg_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=%.7f\n",iter++,pre[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: carg_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_clog_ymm8c4();

__attribute__((hot))
__attribute__((noinline))
void unit_test_clog_ymm8c4()
{
     using namespace gms::math;
     constexpr int32_t buf_len{64};
     constexpr int32_t len_in_float{buf_len*8ull};
     __attribute__((aligned(32)))
     ymm8c4_t buf_x[buf_len];
     __attribute__((aligned(32)))
     ymm8c4_t alignas(32) buf_z[buf_len];
     std::clock_t seed_x{0ull};
     std::clock_t seed_y{0ull};
     std::clock_t seed_distro{};
     int32_t iter{-1};
     int32_t which{-1};
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
               printf("[UNIT-TEST]: -- START: clog_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                                                  
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j]    = clog_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: clog_ymm8c4().\n");
            }
            break;
            case 1 : 
            {
               std::normal_distribution<float> norm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(norm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(norm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::normal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: clog_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                     
                            
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j]    = clog_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: clog_ymm8c4().\n");
            }
            break;
            case 2 : 
            {
               std::lognormal_distribution<float> lognorm_distro;
               char * distro_name{NULL};
               int32_t status{9999};
               distro_name = gms::common::demangle(typeid(lognorm_distro).name(),status);
               if(distro_name != NULL && status == 0)
               {
                    printf("[UNIT-TEST]: Instantiated distribution of type: %s\n\n", distro_name);
                    gms::common::gms_mm_free(distro_name);
               }
               else
               {
                    printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lognorm_distro).name());
                    if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
               }
               seed_x = std::clock();
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: clog_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                                    
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j]    = clog_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: clog_ymm8c4().\n");
            }
            break;
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
               auto rand_x{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_x))};
               seed_y = std::clock();
               auto rand_y{std::bind(std::lognormal_distribution<float>(0.0f,1.57079632679489661923132169164f),
                                   std::mt19937(seed_y))};
               printf("[UNIT-TEST]: -- START: clog_ymm8c4().\n");
               for(int32_t __j{0}; __j != buf_len; ++__j)
               {
                   __m256 vrx;
                   __m256 vix;
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                      
                      vrx = _mm256_setr_ps(rand_x(),rand_x(),rand_x(),rand_x(),
                                           rand_x(),rand_x(),rand_x(),rand_x());
                      vix = _mm256_setr_ps(rand_y(),rand_y(),rand_y(),rand_y(),
                                           rand_y(),rand_y(),rand_y(),rand_y());
                      
                                      
                   }
                   buf_x[__j].re = vrx;
                   buf_x[__j].im = vix;
                   buf_z[__j]    = clog_ymm8c4(buf_x[__j]);
                   float * __restrict pre = reinterpret_cast<float *>(&buf_z[__j].re);
                   float * __restrict pim = reinterpret_cast<float *>(&buf_z[__j].im);
                   for(int32_t __i{0}; __i != 8; ++__i)
                   {
                       printf("[UNIT_TEST]: iter=%d,result=(%.7f,%.7f)\n",iter++,pre[__i],pim[__i]);
                   }
               }
               printf("[UNIT-TEST]: -- END: clog_ymm8c4().\n");
            }
            break;
            default : 
            {
                  printf("[UNIT-TEST]: Invalid switch variable=%d\n",which);
                  return;
            }
     }

}





int main()
{   
    unit_test_cadd_ymm8c4();
    unit_test_csub_ymm8c4();
    unit_test_cmul_ymm8c4();
    unit_test_cdiv_ymm8c4();
    unit_test_cdiv_smith_ymm8c4();
    unit_test_cabs_ymm8c4();
    unit_test_carg_ymm8c4();
    unit_test_clog_ymm8c4();

    return 0;
}
