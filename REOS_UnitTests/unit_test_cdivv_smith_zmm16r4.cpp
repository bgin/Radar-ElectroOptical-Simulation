#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h>
#include <omp.h>
#include <limits>
#include "GMS_malloc.h"
#include "GMS_cdiv_smith_vec_zmm16r4.h"

/*
   icpc -o unit_test_cdiv_smith_vec_zmm16r4_u -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_cdiv_smith_vec_zmm16r4.h GMS_cdiv_smith_vec_zmm16r4.cpp unit_test_cdiv_smith_vec_zmm16r4.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_cdiv_smith_vec_zmm16r4.h GMS_cdiv_smith_vec_zmm16r4.cpp unit_test_cdiv_smith_vec_zmm16r4.cpp
*/


inline static 
void zero_removal(float * __restrict__ data, int32_t n,
                  int32_t & zero_cnt)
{
    constexpr float eps{std::numeric_limits<float>::epsilon()};
    zero_cnt = 0;
    for(int32_t __i{0}; __i != n; ++__i)
    {
        if(data[__i] == 0.0f || data[__i] == -0.0f)
        {   zero_cnt++;
            data[__i] = eps;
        }
    }
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_cdivv_smith_zmm16r4_u10x_a();

void unit_test_cdivv_smith_zmm16r4_u10x_a()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    float * __restrict__ xre{NULL};
    float * __restrict__ xim{NULL};
    float * __restrict__ yre{NULL};
    float * __restrict__ yim{NULL};
    float * __restrict__ zre{NULL};
    float * __restrict__ zim{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_xre{0ull};
    unsigned long long seed_xim{0ull};
    unsigned long long seed_yre{0ull};
    unsigned long long seed_yim{0ull};
    float rxim,rxre;
    float ryim,ryre;
    int32_t zcnt1{0};
    int32_t zcnt2{0};
    int32_t zcnt3{0};
    int32_t zcnt4{0};
    int32_t which{-1};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{"UNIT-TEST_uni_distro_cdivv_smith_zmm16r4_u10x_a.csv"};
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
            xre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            xim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xre))};
            seed_yim = __rdtsc();
            auto rand_yim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_yim))};
            seed_yre = __rdtsc();
            auto rand_yre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_yre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                ryim = rand_yim();
                yim[__i+0] = ryim;
                ryre = rand_yre();
                yre[__i+0] = ryre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                ryim = rand_yim();
                yim[__i+1] = ryim;
                ryre = rand_yre();
                yre[__i+1] = ryre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                ryim = rand_yim();
                yim[__i+2] = ryim;
                ryre = rand_yre();
                yre[__i+2] = ryre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                ryim = rand_yim();
                yim[__i+3] = ryim;
                ryre = rand_yre();
                yre[__i+3] = ryre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                ryim = rand_yim();
                yim[__i+4] = ryim;
                ryre = rand_yre();
                yre[__i+4] = ryre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                ryim = rand_yim();
                yim[__i+5] = ryim;
                ryre = rand_yre();
                yre[__i+5] = ryre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                ryim = rand_yim();
                yim[__i+6] = ryim;
                ryre = rand_yre();
                yre[__i+6] = ryre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                ryim = rand_yim();
                yim[__i+7] = ryim;
                ryre = rand_yre();
                yre[__i+7] = ryre;
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                ryim = rand_yim();
                yim[__i] = ryim;
                ryre = rand_yre();
                yre[__i] = ryre;
            }
#pragma omp parallel sections 
{
   
        #pragma omp section
        {
            zero_removal(xre,nelems,zcnt1);
        }
        #pragma omp section
        {
            zero_removal(xim,nelems,zcnt2);
        }
        #pragma omp section 
        {
            zero_removal(yre,nelems,zcnt3);
        }
        #pragma omp section 
        {
            zero_removal(yim,nelems,zcnt4);
        }
   
}           
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("Number(s) of zeros found: xre=%d,xim%d,yre=%d,yim=%d\n",zcnt1,zcnt2,zcnt3,zcnt4);
            printf("[UNIT-TEST]: -- START: cdivv_smith_zmm16r4_u10x_a\n");
            gms::math::cdivv_smith_zmm16r4_u10x_a(&xre[0],&xim[0],&yre[0],&yim[0],
                                                  &zre[0],&zim[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,j%.7f]\n",zre[__i],zim[__i]);
            }
            if(zre) gms_mm_free(zre);
            if(zim) gms_mm_free(zim);
            if(yre) gms_mm_free(yre);
            if(yim) gms_mm_free(yim);
            if(xre) gms_mm_free(xre);
            if(xim) gms_mm_free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cdivv_smith_zmm16r4_u10x_a\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{"UNIT-TEST_norm_distro_cdivv_smith_zmm16r4_u10x_a.csv"};
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

            xre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            xim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_xre))};
            seed_yim = __rdtsc();
            auto rand_yim{std::bind(std::normal_distribution<float>(0.0f,5.0f),
                                     std::mt19937(seed_yim))};
            seed_yre = __rdtsc();
            auto rand_yre{std::bind(std::normal_distribution<float>(-2.0f,0.5f),
                                     std::mt19937(seed_yre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                ryim = rand_yim();
                yim[__i+0] = ryim;
                ryre = rand_yre();
                yre[__i+0] = ryre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                ryim = rand_yim();
                yim[__i+1] = ryim;
                ryre = rand_yre();
                yre[__i+1] = ryre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                ryim = rand_yim();
                yim[__i+2] = ryim;
                ryre = rand_yre();
                yre[__i+2] = ryre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                ryim = rand_yim();
                yim[__i+3] = ryim;
                ryre = rand_yre();
                yre[__i+3] = ryre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                ryim = rand_yim();
                yim[__i+4] = ryim;
                ryre = rand_yre();
                yre[__i+4] = ryre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                ryim = rand_yim();
                yim[__i+5] = ryim;
                ryre = rand_yre();
                yre[__i+5] = ryre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                ryim = rand_yim();
                yim[__i+6] = ryim;
                ryre = rand_yre();
                yre[__i+6] = ryre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                ryim = rand_yim();
                yim[__i+7] = ryim;
                ryre = rand_yre();
                yre[__i+7] = ryre;
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                ryim = rand_yim();
                yim[__i] = ryim;
                ryre = rand_yre();
                yre[__i] = ryre;
            }
#pragma omp parallel sections 
{
   
        #pragma omp section 
        {
            zero_removal(xre,nelems,zcnt1);
        }
        #pragma omp section 
        {
            zero_removal(xim,nelems,zcnt2);
        }
        #pragma omp section 
        {
            zero_removal(yre,nelems,zcnt3);
        }
        #pragma omp section 
        {
            zero_removal(yim,nelems,zcnt4);
        }
   
}
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("Number(s) of zeros found: xre=%d,xim%d,yre=%d,yim=%d\n",zcnt1,zcnt2,zcnt3,zcnt4);
            printf("[UNIT-TEST]: -- START: cdivv_smith_zmm16r4_u10x_a\n");
            gms::math::cdivv_smith_zmm16r4_u10x_a(&xre[0],&xim[0],&yre[0],&yim[0],
                                                  &zre[0],&zim[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,j%.7f]\n",zre[__i],zim[__i]);
            }
            if(zre) gms_mm_free(zre);
            if(zim) gms_mm_free(zim);
            if(yre) gms_mm_free(yre);
            if(yim) gms_mm_free(yim);
            if(xre) gms_mm_free(xre);
            if(xim) gms_mm_free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cdivv_smith_zmm16r4_u10x_a\n");
        }
        break;
        default : 
                 printf("Invalid switch variable=%d\n",which);
                 return;
    } 
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_cdivv_smith_zmm16r4_u8x_a();

void unit_test_cdivv_smith_zmm16r4_u8x_a()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    float * __restrict__ xre{NULL};
    float * __restrict__ xim{NULL};
    float * __restrict__ yre{NULL};
    float * __restrict__ yim{NULL};
    float * __restrict__ zre{NULL};
    float * __restrict__ zim{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_xre{0ull};
    unsigned long long seed_xim{0ull};
    unsigned long long seed_yre{0ull};
    unsigned long long seed_yim{0ull};
    float rxim,rxre;
    float ryim,ryre;
    int32_t zcnt1{0};
    int32_t zcnt2{0};
    int32_t zcnt3{0};
    int32_t zcnt4{0};
    int32_t which{-1};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{"UNIT-TEST_uni_distro_cdivv_smith_zmm16r4_u8x_a.csv"};
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
            xre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            xim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xre))};
            seed_yim = __rdtsc();
            auto rand_yim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_yim))};
            seed_yre = __rdtsc();
            auto rand_yre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_yre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                ryim = rand_yim();
                yim[__i+0] = ryim;
                ryre = rand_yre();
                yre[__i+0] = ryre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                ryim = rand_yim();
                yim[__i+1] = ryim;
                ryre = rand_yre();
                yre[__i+1] = ryre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                ryim = rand_yim();
                yim[__i+2] = ryim;
                ryre = rand_yre();
                yre[__i+2] = ryre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                ryim = rand_yim();
                yim[__i+3] = ryim;
                ryre = rand_yre();
                yre[__i+3] = ryre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                ryim = rand_yim();
                yim[__i+4] = ryim;
                ryre = rand_yre();
                yre[__i+4] = ryre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                ryim = rand_yim();
                yim[__i+5] = ryim;
                ryre = rand_yre();
                yre[__i+5] = ryre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                ryim = rand_yim();
                yim[__i+6] = ryim;
                ryre = rand_yre();
                yre[__i+6] = ryre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                ryim = rand_yim();
                yim[__i+7] = ryim;
                ryre = rand_yre();
                yre[__i+7] = ryre;
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                ryim = rand_yim();
                yim[__i] = ryim;
                ryre = rand_yre();
                yre[__i] = ryre;
            }
#pragma omp parallel sections 
{
   
        #pragma omp section
        {
            zero_removal(xre,nelems,zcnt1);
        }
        #pragma omp section
        {
            zero_removal(xim,nelems,zcnt2);
        }
        #pragma omp section 
        {
            zero_removal(yre,nelems,zcnt3);
        }
        #pragma omp section 
        {
            zero_removal(yim,nelems,zcnt4);
        }
   
}           
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("Number(s) of zeros found: xre=%d,xim%d,yre=%d,yim=%d\n",zcnt1,zcnt2,zcnt3,zcnt4);
            printf("[UNIT-TEST]: -- START: cdivv_smith_zmm16r4_u8x_a\n");
            gms::math::cdivv_smith_zmm16r4_u8x_a(&xre[0],&xim[0],&yre[0],&yim[0],
                                                  &zre[0],&zim[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,j%.7f]\n",zre[__i],zim[__i]);
            }
            if(zre) gms_mm_free(zre);
            if(zim) gms_mm_free(zim);
            if(yre) gms_mm_free(yre);
            if(yim) gms_mm_free(yim);
            if(xre) gms_mm_free(xre);
            if(xim) gms_mm_free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cdivv_smith_zmm16r4_u8x_a\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{"UNIT-TEST_norm_distro_cdivv_smith_zmm16r4_u8x_a.csv"};
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

            xre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            xim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            yim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zre = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            zim = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_xre))};
            seed_yim = __rdtsc();
            auto rand_yim{std::bind(std::normal_distribution<float>(0.0f,5.0f),
                                     std::mt19937(seed_yim))};
            seed_yre = __rdtsc();
            auto rand_yre{std::bind(std::normal_distribution<float>(-2.0f,0.5f),
                                     std::mt19937(seed_yre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                ryim = rand_yim();
                yim[__i+0] = ryim;
                ryre = rand_yre();
                yre[__i+0] = ryre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                ryim = rand_yim();
                yim[__i+1] = ryim;
                ryre = rand_yre();
                yre[__i+1] = ryre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                ryim = rand_yim();
                yim[__i+2] = ryim;
                ryre = rand_yre();
                yre[__i+2] = ryre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                ryim = rand_yim();
                yim[__i+3] = ryim;
                ryre = rand_yre();
                yre[__i+3] = ryre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                ryim = rand_yim();
                yim[__i+4] = ryim;
                ryre = rand_yre();
                yre[__i+4] = ryre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                ryim = rand_yim();
                yim[__i+5] = ryim;
                ryre = rand_yre();
                yre[__i+5] = ryre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                ryim = rand_yim();
                yim[__i+6] = ryim;
                ryre = rand_yre();
                yre[__i+6] = ryre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                ryim = rand_yim();
                yim[__i+7] = ryim;
                ryre = rand_yre();
                yre[__i+7] = ryre;
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                ryim = rand_yim();
                yim[__i] = ryim;
                ryre = rand_yre();
                yre[__i] = ryre;
            }
#pragma omp parallel sections 
{
   
        #pragma omp section 
        {
            zero_removal(xre,nelems,zcnt1);
        }
        #pragma omp section 
        {
            zero_removal(xim,nelems,zcnt2);
        }
        #pragma omp section 
        {
            zero_removal(yre,nelems,zcnt3);
        }
        #pragma omp section 
        {
            zero_removal(yim,nelems,zcnt4);
        }
   
}
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("Number(s) of zeros found: xre=%d,xim%d,yre=%d,yim=%d\n",zcnt1,zcnt2,zcnt3,zcnt4);
            printf("[UNIT-TEST]: -- START: cdivv_smith_zmm16r4_u8x_a\n");
            gms::math::cdivv_smith_zmm16r4_u8x_a(&xre[0],&xim[0],&yre[0],&yim[0],
                                                  &zre[0],&zim[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,j%.7f]\n",zre[__i],zim[__i]);
            }
            if(zre) gms_mm_free(zre);
            if(zim) gms_mm_free(zim);
            if(yre) gms_mm_free(yre);
            if(yim) gms_mm_free(yim);
            if(xre) gms_mm_free(xre);
            if(xim) gms_mm_free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cdivv_smith_zmm16r4_u8x_a\n");
        }
        break;
        default : 
                 printf("Invalid switch variable=%d\n",which);
                 return;
    } 
}


int main()
{
     unit_test_cdivv_smith_zmm16r4_u10x_a();
     unit_test_cdivv_smith_zmm16r4_u8x_a();
     return 0;
}