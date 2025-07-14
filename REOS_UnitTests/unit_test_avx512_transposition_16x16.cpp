#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h> // __rdtsc()
#include "GMS_malloc.h" // for demangle 
#include "GMS_avx512_transposition_16x16.h"

/*
   icpc -o unit_test_avx512_transposition_16x16 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_avx512_transposition_16x16.h GMS_avx512_transposition_16x16.cpp unit_test_avx512_transposition_16x16.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_avx512_transposition_16x16.h GMS_avx512_transposition_16x16.cpp unit_test_avx512_transposition_16x16.cpp 
*/
#define ZMM_LEN 16 

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_transpose_zmm16r4_16x16_ip();

void unit_test_transpose_zmm16r4_16x16_ip()
{
    constexpr int32_t nelems{256};
    float buf_init[nelems];
    float buf_transposed[nelems];
    unsigned long long seed_b{0ull};
    std::clock_t seed_d{0ull};
    int32_t which(-1);
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             __m512 x0;
             __m512 x1;
             __m512 x2;
             __m512 x3;
             __m512 x4;
             __m512 x5;
             __m512 x6;
             __m512 x7;
             __m512 x8;
             __m512 x9;
             __m512 x10;
             __m512 x11;
             __m512 x12;
             __m512 x13;
             __m512 x14;
             __m512 x15;
             seed_b = __rdtsc();
             auto rand_b{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd = rand_b();
                 buf_init[__i] = rd;
             }
              x0 = _mm512_load_ps(&buf_init[0*ZMM_LEN]);
             x1 = _mm512_load_ps(&buf_init[1*ZMM_LEN]);
             x2 = _mm512_load_ps(&buf_init[2*ZMM_LEN]);
             x3 = _mm512_load_ps(&buf_init[3*ZMM_LEN]);
             x4 = _mm512_load_ps(&buf_init[4*ZMM_LEN]);
             x5 = _mm512_load_ps(&buf_init[5*ZMM_LEN]);
             x6 = _mm512_load_ps(&buf_init[6*ZMM_LEN]);
             x7 = _mm512_load_ps(&buf_init[7*ZMM_LEN]);
             x8 = _mm512_load_ps(&buf_init[8*ZMM_LEN]);
             x9 = _mm512_load_ps(&buf_init[9*ZMM_LEN]);
             x10 = _mm512_load_ps(&buf_init[10*ZMM_LEN]);
             x11 = _mm512_load_ps(&buf_init[11*ZMM_LEN]);
             x12 = _mm512_load_ps(&buf_init[12*ZMM_LEN]);
             x13 = _mm512_load_ps(&buf_init[13*ZMM_LEN]);
             x14 = _mm512_load_ps(&buf_init[14*ZMM_LEN]);
             x15 = _mm512_load_ps(&buf_init[15*ZMM_LEN]);
             printf("[UNIT-TEST]: -- START: transpose_zmm16r4_16x16_ip().\n");
             gms::math::transpose_zmm16r4_16x16_ip(x0,x1,x2,x3,x4,x5,x6,x7,
                                                   x8,x9,x10,x11,x12,x13,x14,x15);
             _mm512_storeu_ps(&buf_transposed[0*ZMM_LEN],x0);
             _mm512_storeu_ps(&buf_transposed[1*ZMM_LEN],x1);
             _mm512_storeu_ps(&buf_transposed[2*ZMM_LEN],x2);
             _mm512_storeu_ps(&buf_transposed[3*ZMM_LEN],x3);
             _mm512_storeu_ps(&buf_transposed[4*ZMM_LEN],x4);
             _mm512_storeu_ps(&buf_transposed[5*ZMM_LEN],x5);
             _mm512_storeu_ps(&buf_transposed[6*ZMM_LEN],x6);
             _mm512_storeu_ps(&buf_transposed[7*ZMM_LEN],x7);
             _mm512_storeu_ps(&buf_transposed[8*ZMM_LEN],x8);
             _mm512_storeu_ps(&buf_transposed[9*ZMM_LEN],x9);
             _mm512_storeu_ps(&buf_transposed[10*ZMM_LEN],x10);
             _mm512_storeu_ps(&buf_transposed[11*ZMM_LEN],x11);
             _mm512_storeu_ps(&buf_transposed[12*ZMM_LEN],x12);
             _mm512_storeu_ps(&buf_transposed[13*ZMM_LEN],x13);
             _mm512_storeu_ps(&buf_transposed[14*ZMM_LEN],x14);
             _mm512_storeu_ps(&buf_transposed[15*ZMM_LEN],x15);
             printf("[UNIT-TEST:] -- Dumping results:\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("init[%d]=%.7f,transposed[%d]=%.7f\n",__i,buf_init[__i],__i,buf_transposed[__i]);
             }
             
             printf("[UNIT-TEST]: -- END: transpose_zmm16r4_16x16_ip().\n");
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             __m512 x0;
             __m512 x1;
             __m512 x2;
             __m512 x3;
             __m512 x4;
             __m512 x5;
             __m512 x6;
             __m512 x7;
             __m512 x8;
             __m512 x9;
             __m512 x10;
             __m512 x11;
             __m512 x12;
             __m512 x13;
             __m512 x14;
             __m512 x15;
             seed_b = __rdtsc();
             auto rand_b{std::bind(std::normal_distribution<float>(0.1f,2.5f),std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd = rand_b();
                 buf_init[__i] = rd;
             }
             x0 = _mm512_load_ps(&buf_init[0*ZMM_LEN]);
             x1 = _mm512_load_ps(&buf_init[1*ZMM_LEN]);
             x2 = _mm512_load_ps(&buf_init[2*ZMM_LEN]);
             x3 = _mm512_load_ps(&buf_init[3*ZMM_LEN]);
             x4 = _mm512_load_ps(&buf_init[4*ZMM_LEN]);
             x5 = _mm512_load_ps(&buf_init[5*ZMM_LEN]);
             x6 = _mm512_load_ps(&buf_init[6*ZMM_LEN]);
             x7 = _mm512_load_ps(&buf_init[7*ZMM_LEN]);
             x8 = _mm512_load_ps(&buf_init[8*ZMM_LEN]);
             x9 = _mm512_load_ps(&buf_init[9*ZMM_LEN]);
             x10 = _mm512_load_ps(&buf_init[10*ZMM_LEN]);
             x11 = _mm512_load_ps(&buf_init[11*ZMM_LEN]);
             x12 = _mm512_load_ps(&buf_init[12*ZMM_LEN]);
             x13 = _mm512_load_ps(&buf_init[13*ZMM_LEN]);
             x14 = _mm512_load_ps(&buf_init[14*ZMM_LEN]);
             x15 = _mm512_load_ps(&buf_init[15*ZMM_LEN]);
             printf("[UNIT-TEST]: -- START: transpose_zmm16r4_16x16_ip().\n");
             gms::math::transpose_zmm16r4_16x16_ip(x0,x1,x2,x3,x4,x5,x6,x7,
                                                   x8,x9,x10,x11,x12,x13,x14,x15);
              _mm512_storeu_ps(&buf_transposed[0*ZMM_LEN],x0);
             _mm512_storeu_ps(&buf_transposed[1*ZMM_LEN],x1);
             _mm512_storeu_ps(&buf_transposed[2*ZMM_LEN],x2);
             _mm512_storeu_ps(&buf_transposed[3*ZMM_LEN],x3);
             _mm512_storeu_ps(&buf_transposed[4*ZMM_LEN],x4);
             _mm512_storeu_ps(&buf_transposed[5*ZMM_LEN],x5);
             _mm512_storeu_ps(&buf_transposed[6*ZMM_LEN],x6);
             _mm512_storeu_ps(&buf_transposed[7*ZMM_LEN],x7);
             _mm512_storeu_ps(&buf_transposed[8*ZMM_LEN],x8);
             _mm512_storeu_ps(&buf_transposed[9*ZMM_LEN],x9);
             _mm512_storeu_ps(&buf_transposed[10*ZMM_LEN],x10);
             _mm512_storeu_ps(&buf_transposed[11*ZMM_LEN],x11);
             _mm512_storeu_ps(&buf_transposed[12*ZMM_LEN],x12);
             _mm512_storeu_ps(&buf_transposed[13*ZMM_LEN],x13);
             _mm512_storeu_ps(&buf_transposed[14*ZMM_LEN],x14);
             _mm512_storeu_ps(&buf_transposed[15*ZMM_LEN],x15);
             printf("[UNIT-TEST:] -- Dumping results:\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("init[%d]=%.7f,transposed[%d]=%.7f\n",__i,buf_init[__i],__i,buf_transposed[__i]);
             }
             
             printf("[UNIT-TEST]: -- END: transpose_zmm16r4_16x16_ip().\n");
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
void unit_test_transpose_u_zmm16r4_16x16();

void unit_test_transpose_u_zmm16r4_16x16()
{
    constexpr int32_t nelems{256};
    float buf_init[nelems];
    float buf_transposed[nelems];
    unsigned long long seed_b{0ull};
    std::clock_t seed_d{0ull};
    int32_t which(-1);
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             
             seed_b = __rdtsc();
             auto rand_b{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd = rand_b();
                 buf_init[__i] = rd;
             }
             
             printf("[UNIT-TEST]: -- START: transpose_u_zmm16r4_16x16().\n");
             gms::math::transpose_u_zmm16r4_16x16(&buf_init[0],&buf_transposed[0]);
             printf("[UNIT-TEST:] -- Dumping results:\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("init[%d]=%.7f,transposed[%d]=%.7f\n",__i,buf_init[__i],__i,buf_transposed[__i]);
             }
             
             printf("[UNIT-TEST]: -- END: transpose_u_zmm16r4_16x16().\n");
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
            
             seed_b = __rdtsc();
             auto rand_b{std::bind(std::normal_distribution<float>(0.1f,2.5f),std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd = rand_b();
                 buf_init[__i] = rd;
             }
             
             printf("[UNIT-TEST]: -- START: transpose_u_zmm16r4_16x16().\n");
             gms::math::transpose_u_zmm16r4_16x16(&buf_init[0],&buf_transposed[0]);
             printf("[UNIT-TEST:] -- Dumping results:\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("init[%d]=%.7f,transposed[%d]=%.7f\n",__i,buf_init[__i],__i,buf_transposed[__i]);
             }
             
             printf("[UNIT-TEST]: -- END: transpose_u_zmm16r4_16x16().\n");
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
void unit_test_transpose_a_zmm16r4_16x16();

void unit_test_transpose_a_zmm16r4_16x16()
{
    constexpr int32_t nelems{256};
    float __ATTR_ALIGN__(64) buf_init[nelems];
    float __ATTR_ALIGN__(64) buf_transposed[nelems];
    unsigned long long seed_b{0ull};
    std::clock_t seed_d{0ull};
    int32_t which(-1);
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
             
             seed_b = __rdtsc();
             auto rand_b{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd = rand_b();
                 buf_init[__i] = rd;
             }
             
             printf("[UNIT-TEST]: -- START: transpose_a_zmm16r4_16x16().\n");
             gms::math::transpose_a_zmm16r4_16x16(&buf_init[0],&buf_transposed[0]);
             printf("[UNIT-TEST:] -- Dumping results:\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("init[%d]=%.7f,transposed[%d]=%.7f\n",__i,buf_init[__i],__i,buf_transposed[__i]);
             }
             
             printf("[UNIT-TEST]: -- END: transpose_a_zmm16r4_16x16().\n");
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
             }
            
             seed_b = __rdtsc();
             auto rand_b{std::bind(std::normal_distribution<float>(0.1f,2.5f),std::mt19937(seed_b))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd = rand_b();
                 buf_init[__i] = rd;
             }
             
             printf("[UNIT-TEST]: -- START: transpose_a_zmm16r4_16x16().\n");
             gms::math::transpose_a_zmm16r4_16x16(&buf_init[0],&buf_transposed[0]);
             printf("[UNIT-TEST:] -- Dumping results:\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 printf("init[%d]=%.7f,transposed[%d]=%.7f\n",__i,buf_init[__i],__i,buf_transposed[__i]);
             }
             
             printf("[UNIT-TEST]: -- END: transpose_a_zmm16r4_16x16().\n");
        }
        break;
        default : 
        printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }
}





int main()
{
    unit_test_transpose_zmm16r4_16x16_ip();
    unit_test_transpose_u_zmm16r4_16x16();
    unit_test_transpose_a_zmm16r4_16x16();
}