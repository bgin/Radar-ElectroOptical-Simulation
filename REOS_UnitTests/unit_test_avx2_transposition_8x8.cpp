
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h" // for demangle 
#include "GMS_avx2_transposition_8x8.h"

/*
   icpc -o unit_test_avx2_transposition_8x8 -fp-model fast=2 -ftz -ggdb -ipo -xHost -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_avx2_transposition_8x8.h GMS_avx2_transposition_8x8.cpp unit_test_avx2_transposition_8x8.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -xHost -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_avx2_transposition_8x8.h GMS_avx2_transposition_8x8.cpp unit_test_avx2_transposition_8x8.cpp 
*/

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_transpose_ymm8r4_8x8_ip();

void unit_test_transpose_ymm8r4_8x8_ip()
{
    constexpr int32_t nelems{64};
    float __attribute__((aligned(32))) buf_init[nelems];
    float __attribute__((aligned(32))) buf_transposed[nelems];
    std::clock_t seed_b{0ull};
    std::clock_t seed_d{0ull};
    int32_t which{-1};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),
                std::mt19937(seed_d))};
    which = rand_d();
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
             __m256 x0;
             __m256 x1;
             __m256 x2;
             __m256 x3;
             __m256 x4;
             __m256 x5;
             __m256 x6;
             __m256 x7;
             seed_b = std::clock();
             auto rand_b{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),
                          std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                  float rd{rand_b()};
                  buf_init[__i] = rd;
             }
             x0 = _mm256_load_ps(&buf_init[0]);
             x1 = _mm256_load_ps(&buf_init[8]);
             x2 = _mm256_load_ps(&buf_init[16]);
             x3 = _mm256_load_ps(&buf_init[24]);
             x4 = _mm256_load_ps(&buf_init[32]);
             x5 = _mm256_load_ps(&buf_init[40]);
             x6 = _mm256_load_ps(&buf_init[48]);
             x7 = _mm256_load_ps(&buf_init[56]);
             printf("[UNIT-TEST]: -- START: transpose_ymm8r4_8x8_ip().\n");
             gms::math::transpose_ymm8r4_8x8_ip(x0,x1,x2,x3,
                                                x4,x5,x6,x7);
             _mm256_store_ps(&buf_transposed[0], x0);
             _mm256_store_ps(&buf_transposed[8], x1);
             _mm256_store_ps(&buf_transposed[16],x2);
             _mm256_store_ps(&buf_transposed[24],x3);
             _mm256_store_ps(&buf_transposed[32],x4);
             _mm256_store_ps(&buf_transposed[40],x5);
             _mm256_store_ps(&buf_transposed[48],x6);
             _mm256_store_ps(&buf_transposed[56],x7);
             printf("[UNIT-TEST:] -- Dumping results: Non-Transposed (initial state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_init[%d]=%.7f\n",__i,buf_init[__i]);
             }
             printf("[UNIT-TEST:] -- Dumping results: Transposed (final state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_transposed[%d]=%.7f\n",__i,buf_transposed[__i]);
             }
             printf("[UNIT-TEST]: -- END: transpose_ymm8r4_8x8_ip().\n");
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
             __m256 x0;
             __m256 x1;
             __m256 x2;
             __m256 x3;
             __m256 x4;
             __m256 x5;
             __m256 x6;
             __m256 x7;
             seed_b = std::clock();
             auto rand_b{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                          std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                  float rd{rand_b()};
                  buf_init[__i] = rd;
             }
             x0 = _mm256_load_ps(&buf_init[0]);
             x1 = _mm256_load_ps(&buf_init[8]);
             x2 = _mm256_load_ps(&buf_init[16]);
             x3 = _mm256_load_ps(&buf_init[24]);
             x4 = _mm256_load_ps(&buf_init[32]);
             x5 = _mm256_load_ps(&buf_init[40]);
             x6 = _mm256_load_ps(&buf_init[48]);
             x7 = _mm256_load_ps(&buf_init[56]);
             printf("[UNIT-TEST]: -- START: transpose_ymm8r4_8x8_ip().\n");
             gms::math::transpose_ymm8r4_8x8_ip(x0,x1,x2,x3,
                                                x4,x5,x6,x7);
             _mm256_store_ps(&buf_transposed[0], x0);
             _mm256_store_ps(&buf_transposed[8], x1);
             _mm256_store_ps(&buf_transposed[16],x2);
             _mm256_store_ps(&buf_transposed[24],x3);
             _mm256_store_ps(&buf_transposed[32],x4);
             _mm256_store_ps(&buf_transposed[40],x5);
             _mm256_store_ps(&buf_transposed[48],x6);
             _mm256_store_ps(&buf_transposed[56],x7);
             printf("[UNIT-TEST:] -- Dumping results: Non-Transposed (initial state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_init[%d]=%.7f\n",__i,buf_init[__i]);
             }
             printf("[UNIT-TEST:] -- Dumping results: Transposed (final state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_transposed[%d]=%.7f\n",__i,buf_transposed[__i]);
             }
             printf("[UNIT-TEST]: -- END: transpose_ymm8r4_8x8_ip().\n");
        }
        break;
        default : 
                 printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_transpose_u_ymm8r4_8x8_ip();

void unit_test_transpose_u_ymm8r4_8x8_ip()
{
    constexpr int32_t nelems{64};
    float buf_init[nelems];
    float buf_transposed[nelems];
    std::clock_t seed_b{0ull};
    std::clock_t seed_d{0ull};
    int32_t which{-1};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),
                std::mt19937(seed_d))};
    which = rand_d();
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
             float * __restrict__ x0{nullptr};
             float * __restrict__ x1{nullptr};
             float * __restrict__ x2{nullptr};
             float * __restrict__ x3{nullptr};
             float * __restrict__ x4{nullptr};
             float * __restrict__ x5{nullptr};
             float * __restrict__ x6{nullptr};
             float * __restrict__ x7{nullptr};
             seed_b = std::clock();
             auto rand_b{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),
                          std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                  float rd{rand_b()};
                  buf_init[__i] = rd;
             }
             x0 = &buf_init[0];
             x1 = &buf_init[8];
             x2 = &buf_init[16];
             x3 = &buf_init[24];
             x4 = &buf_init[32];
             x5 = &buf_init[40];
             x6 = &buf_init[48];
             x7 = &buf_init[56];
             printf("[UNIT-TEST]: -- START: transpose_u_ymm8r4_8x8_ip().\n");
             printf("[UNIT-TEST:] -- Dumping results: Non-Transposed (initial state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_init[%d]=%.7f\n",__i,buf_init[__i]);
             }
             gms::math::transpose_u_ymm8r4_8x8_ip(x0,x1,x2,x3,
                                                  x4,x5,x6,x7);
             _mm256_storeu_ps(&buf_transposed[0], _mm256_loadu_ps(x0));
             _mm256_storeu_ps(&buf_transposed[8], _mm256_loadu_ps(x1));
             _mm256_storeu_ps(&buf_transposed[16],_mm256_loadu_ps(x2));
             _mm256_storeu_ps(&buf_transposed[24],_mm256_loadu_ps(x3));
             _mm256_storeu_ps(&buf_transposed[32],_mm256_loadu_ps(x4));
             _mm256_storeu_ps(&buf_transposed[40],_mm256_loadu_ps(x5));
             _mm256_storeu_ps(&buf_transposed[48],_mm256_loadu_ps(x6));
             _mm256_storeu_ps(&buf_transposed[56],_mm256_loadu_ps(x7));
             
             printf("[UNIT-TEST:] -- Dumping results: Transposed (final state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_transposed[%d]=%.7f\n",__i,buf_transposed[__i]);
             }
             printf("[UNIT-TEST]: -- END: transpose_u_ymm8r4_8x8_ip().\n");
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
             float * __restrict__ x0{nullptr};
             float * __restrict__ x1{nullptr};
             float * __restrict__ x2{nullptr};
             float * __restrict__ x3{nullptr};
             float * __restrict__ x4{nullptr};
             float * __restrict__ x5{nullptr};
             float * __restrict__ x6{nullptr};
             float * __restrict__ x7{nullptr};
             seed_b = std::clock();
             auto rand_b{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                          std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                  float rd{rand_b()};
                  buf_init[__i] = rd;
             }
             x0 = &buf_init[0];
             x1 = &buf_init[8];
             x2 = &buf_init[16];
             x3 = &buf_init[24];
             x4 = &buf_init[32];
             x5 = &buf_init[40];
             x6 = &buf_init[48];
             x7 = &buf_init[56];
             printf("[UNIT-TEST]: -- START: transpose_u_ymm8r4_8x8_ip().\n");
             printf("[UNIT-TEST:] -- Dumping results: Non-Transposed (initial state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_init[%d]=%.7f\n",__i,buf_init[__i]);
             }
             printf("[UNIT-TEST]: -- START: transpose_ymm8r4_8x8_ip().\n");
             gms::math::transpose_u_ymm8r4_8x8_ip(x0,x1,x2,x3,
                                                  x4,x5,x6,x7);
             _mm256_storeu_ps(&buf_transposed[0], _mm256_loadu_ps(x0));
             _mm256_storeu_ps(&buf_transposed[8], _mm256_loadu_ps(x1));
             _mm256_storeu_ps(&buf_transposed[16],_mm256_loadu_ps(x2));
             _mm256_storeu_ps(&buf_transposed[24],_mm256_loadu_ps(x3));
             _mm256_storeu_ps(&buf_transposed[32],_mm256_loadu_ps(x4));
             _mm256_storeu_ps(&buf_transposed[40],_mm256_loadu_ps(x5));
             _mm256_storeu_ps(&buf_transposed[48],_mm256_loadu_ps(x6));
             _mm256_storeu_ps(&buf_transposed[56],_mm256_loadu_ps(x7));
             printf("[UNIT-TEST:] -- Dumping results: Transposed (final state)\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("buf_transposed[%d]=%.7f\n",__i,buf_transposed[__i]);
             }
             printf("[UNIT-TEST]: -- END: transpose_u_ymm8r4_8x8_ip().\n");
        }
        break;
        default : 
                 printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }
}



int main()
{
    unit_test_transpose_ymm8r4_8x8_ip();
    unit_test_transpose_u_ymm8r4_8x8_ip();
   
}