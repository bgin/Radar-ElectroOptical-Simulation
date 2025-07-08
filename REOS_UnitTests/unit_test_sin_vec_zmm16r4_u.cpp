
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h>
#include "GMS_malloc.h" 
#include "GMS_sin_vec_zmm16r4.h"

/*
   icpc -o unit_test_sin_vec_zmm16r4_u -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sin_vec_zmm16r4.h GMS_sin_vec_zmm16r4.cpp unit_test_sin_vec_zmm16r4_u.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_sin_vec_zmm16r4.h GMS_sin_vec_zmm16r4.cpp unit_test_sin_vec_zmm16r4_u.cpp
*/



__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_sinv_zmm16r4_unroll_10x_u();

void unit_test_sinv_zmm16r4_unroll_10x_u()
{
    using namespace gms::common;
    __m512 a,b,c;
    float * __restrict__ p_x{NULL};
    float * __restrict__ p_y{NULL};
    FILE  *  __restrict__ fp{NULL};
   
    const char * __restrict__ fout{"UNIT-TEST_Output_sinv_zmm16r4_unroll_10x_u.csv"};
    constexpr std::size_t nelems{67456};
    constexpr std::size_t nbytes{4ull*nelems};
    std::clock_t seed_d{0ull};
    std::clock_t seed_x{0ull};
    seed_d = std::clock();
    int32_t which{-1};
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
            distro_name = demangle(typeid(uniform_distro).name(),status);
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
            a = _mm512_setr_ps(1.0f,2.0f,3.0f,4.0f,
                               5.0f,6.0f,7.0f,8.0f,
                               9.0f,10.0f,11.0f,12.0f,
                               13.0f,14.0f,15.0f,16.0f);
            b = _mm512_setr_ps(0.1f,0.2f,0.3f,0.4f,
                               0.5f,0.6f,0.7f,0.8f,
                               0.9f,0.11f,0.12f,0.13f,
                               0.14f,0.15f,0.16f,0.17f);
            p_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            p_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_x = std::clock();
            auto rand_x{std::bind(std::uniform_real_distribution<float>(0.0f,3.14159265358979323846264338328f),
                                  std::mt19937(seed_x))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float r0{rand_x()};
                 p_x[__i+0] = r0;
                 float r1{rand_x()};
                 p_x[__i+1] = r1;
                 float r2{rand_x()};
                 p_x[__i+2] = r2;
                 float r3{rand_x()};
                 p_x[__i+3] = r3;
                 float r4{rand_x()};
                 p_x[__i+4] = r4;
                 float r5{rand_x()};
                 p_x[__i+5] = r5;
                 float r6{rand_x()};
                 p_x[__i+6] = r6;
                 float r7{rand_x()};
                 p_x[__i+7] = r7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: sinv_zmm16r4_unroll_10x_u\n");
            gms::math::sinv_zmm16r4_unroll_10x_u(&p_x[0],&p_y[0],a,b,c,nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_x[%d]=%.7f,p_y[%d]=%.7f\n",__i,p_x[__i],__i,p_y[__i]);
            }
            if(p_x) gms_mm_free(p_x);
            if(p_y) gms_mm_free(p_y);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: sinv_zmm16r4_unroll_10x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
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
            a = _mm512_setr_ps(1.0f,2.0f,3.0f,4.0f,
                               5.0f,6.0f,7.0f,8.0f,
                               9.0f,10.0f,11.0f,12.0f,
                               13.0f,14.0f,15.0f,16.0f);
            b = _mm512_setr_ps(0.1f,0.2f,0.3f,0.4f,
                               0.5f,0.6f,0.7f,0.8f,
                               0.9f,0.11f,0.12f,0.13f,
                               0.14f,0.15f,0.16f,0.17f);
            p_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            p_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_x = std::clock();
            auto rand_x{std::bind(std::normal_distribution<float>(0.0f,2.5f),
                                  std::mt19937(seed_x))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float r0{rand_x()};
                 p_x[__i+0] = r0;
                 float r1{rand_x()};
                 p_x[__i+1] = r1;
                 float r2{rand_x()};
                 p_x[__i+2] = r2;
                 float r3{rand_x()};
                 p_x[__i+3] = r3;
                 float r4{rand_x()};
                 p_x[__i+4] = r4;
                 float r5{rand_x()};
                 p_x[__i+5] = r5;
                 float r6{rand_x()};
                 p_x[__i+6] = r6;
                 float r7{rand_x()};
                 p_x[__i+7] = r7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: sinv_zmm16r4_unroll_10x_u\n");
            gms::math::sinv_zmm16r4_unroll_10x_u(&p_x[0],&p_y[0],a,b,c,nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_x[%d]=%.7f,p_y[%d]=%.7f\n",__i,p_x[__i],__i,p_y[__i]);
            }
            if(p_x) gms_mm_free(p_x);
            if(p_y) gms_mm_free(p_y);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: sinv_zmm16r4_unroll_10x_u\n");
        }
    }
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_sinv_zmm16r4_unroll_6x_u();

void unit_test_sinv_zmm16r4_unroll_6x_u()
{
    using namespace gms::common;
    __m512 a,b,c;
    float * __restrict__ p_x{NULL};
    float * __restrict__ p_y{NULL};
    FILE  *  __restrict__ fp{NULL};
    const char * __restrict__ fout{"UNIT-TEST_Output_sinv_zmm16r4_unroll_6x_u.csv"};
    constexpr std::size_t nelems{67456};
    constexpr std::size_t nbytes{4ull*nelems};
    std::clock_t seed_d{0ull};
    std::clock_t seed_x{0ull};
    seed_d = std::clock();
    int32_t which{-1};
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
            distro_name = demangle(typeid(uniform_distro).name(),status);
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
            a = _mm512_setr_ps(1.0f,2.0f,3.0f,4.0f,
                               5.0f,6.0f,7.0f,8.0f,
                               9.0f,10.0f,11.0f,12.0f,
                               13.0f,14.0f,15.0f,16.0f);
            b = _mm512_setr_ps(0.1f,0.2f,0.3f,0.4f,
                               0.5f,0.6f,0.7f,0.8f,
                               0.9f,0.11f,0.12f,0.13f,
                               0.14f,0.15f,0.16f,0.17f);
            p_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            p_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_x = std::clock();
            auto rand_x{std::bind(std::uniform_real_distribution<float>(0.0f,3.14159265358979323846264338328f),
                                  std::mt19937(seed_x))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float r0{rand_x()};
                 p_x[__i+0] = r0;
                 float r1{rand_x()};
                 p_x[__i+1] = r1;
                 float r2{rand_x()};
                 p_x[__i+2] = r2;
                 float r3{rand_x()};
                 p_x[__i+3] = r3;
                 float r4{rand_x()};
                 p_x[__i+4] = r4;
                 float r5{rand_x()};
                 p_x[__i+5] = r5;
                 float r6{rand_x()};
                 p_x[__i+6] = r6;
                 float r7{rand_x()};
                 p_x[__i+7] = r7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: sinv_zmm16r4_unroll_6x_u\n");
            gms::math::sinv_zmm16r4_unroll_6x_u(&p_x[0],&p_y[0],a,b,c,nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_x[%d]=%.7f,p_y[%d]=%.7f\n",__i,p_x[__i],__i,p_y[__i]);
            }
            if(p_x) gms_mm_free(p_x);
            if(p_y) gms_mm_free(p_y);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: sinv_zmm16r4_unroll_6x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
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
            a = _mm512_setr_ps(1.0f,2.0f,3.0f,4.0f,
                               5.0f,6.0f,7.0f,8.0f,
                               9.0f,10.0f,11.0f,12.0f,
                               13.0f,14.0f,15.0f,16.0f);
            b = _mm512_setr_ps(0.1f,0.2f,0.3f,0.4f,
                               0.5f,0.6f,0.7f,0.8f,
                               0.9f,0.11f,0.12f,0.13f,
                               0.14f,0.15f,0.16f,0.17f);
            p_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            p_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_x = std::clock();
            auto rand_x{std::bind(std::normal_distribution<float>(0.0f,2.5f),
                                  std::mt19937(seed_x))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float r0{rand_x()};
                 p_x[__i+0] = r0;
                 float r1{rand_x()};
                 p_x[__i+1] = r1;
                 float r2{rand_x()};
                 p_x[__i+2] = r2;
                 float r3{rand_x()};
                 p_x[__i+3] = r3;
                 float r4{rand_x()};
                 p_x[__i+4] = r4;
                 float r5{rand_x()};
                 p_x[__i+5] = r5;
                 float r6{rand_x()};
                 p_x[__i+6] = r6;
                 float r7{rand_x()};
                 p_x[__i+7] = r7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: sinv_zmm16r4_unroll_6x_u\n");
            gms::math::sinv_zmm16r4_unroll_6x_u(&p_x[0],&p_y[0],a,b,c,nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_x[%d]=%.7f,p_y[%d]=%.7f\n",__i,p_x[__i],__i,p_y[__i]);
            }
            if(p_x) gms_mm_free(p_x);
            if(p_y) gms_mm_free(p_y);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: sinv_zmm16r4_unroll_6x_u\n");
        }
    }
}


int main()
{
    unit_test_sinv_zmm16r4_unroll_10x_u();
    unit_test_sinv_zmm16r4_unroll_6x_u();
    return 0;
}