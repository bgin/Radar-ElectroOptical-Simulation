
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h>
#include <string>
#include "GMS_malloc.h"
#include "GMS_csinhv_svml_zmm16r4.h"

/*
   icpc -o unit_test_csinhv_svml_zmm16r4_u -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5 -fopenmp \
   GMS_config.h GMS_malloc.h GMS_csinhv_svml_zmm16r4.h GMS_csinhv_svml_zmm16r4.cpp unit_test_csinhv_svml_zmm16r4_u.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -fopenmp -falign-functions=32 GMS_config.h GMS_malloc.h GMS_csinhv_svml_zmm16r4.h GMS_csinhv_svml_zmm16r4.cpp unit_test_csinhv_svml_zmm16r4_u.cpp

*/

int32_t funi_distro_idx = 0;
int32_t fnorm_distr_idx = 0;

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void warm_cache_L1I();

void warm_cache_L1I()
{
    volatile float res_sinh[64];
    volatile float res_cosh[64];
    for(int32_t __i{0}; __i != 64; ++__i)
    {
        const float ri{static_cast<float>(__i)};
        res_sinh[__i] = std::sinh(ri);
        res_cosh[__i] = std::cosh(ri);
    }
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_csinhv_svml_zmm16r4_u10x_u();

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_csinhv_svml_zmm16r4_u10x_u()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    constexpr unsigned long long RDTSCP_LAT{42ull};
    std::string fname1 = "UNIT-TEST_uni_distro_csinhv_svml_zmm16r4_u10x_u_" + std::to_string(funi_distro_idx++)+".csv";
    std::string fname2 = "UNIT-TEST_norm_distro_csinhv_svml_zmm16r4_u10x_u_" + std::to_string(fnorm_distr_idx++)+".csv";
    float * __restrict__ xre{NULL};
    float * __restrict__ xim{NULL};
    float * __restrict__ zre{NULL};
    float * __restrict__ zim{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_xre{0ull};
    unsigned long long seed_xim{0ull};
    unsigned long long start{0ull};
    unsigned long long end{0ull};
    unsigned long long start_c{0ull};
    unsigned long long end_c{0ull};
    unsigned long long elapsed_tsc{0ull};
    uint32_t tsc_aux{9999};
    float rxim,rxre;
    int32_t which{-1};
    int32_t idx{};
    bool  xre_neq_0{};
    bool  xim_neq_0{};
    bool  zre_neq_0{};
    bool  zim_neq_0{};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{fname1.c_str()};
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
            xre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xre_neq_0 = (xre != NULL) && (nbytes > 0ull);
            xim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xim_neq_0 = (xim != NULL) && (nbytes > 0ull);
            zre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zre_neq_0 = (zre != NULL) && (nbytes > 0ull);
            zim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zim_neq_0 = (zim != NULL) && (nbytes > 0ull);
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
               
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: csinhv_svml_zmm16r4_u10x_u\n");
            start   = __rdtscp(&tsc_aux);
            gms::math::csinhv_svml_zmm16r4_u10x_u(&xre[0],&xim[0], &zre[0],&zim[0],nelems);
                                                 
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: csinhv_svml_zmm16r4_u10x_u -- TSC=%llu,TSC_AUX=%d\n",elapsed_tsc,tsc_aux);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(zre_neq_0) std::free(zre);
                 if(zim_neq_0) std::free(zim);
                 if(xre_neq_0) std::free(xre);
                 if(xim_neq_0) std::free(xim);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,%.7f]\n",zre[__i],zim[__i]);
            }
            if(zre_neq_0) std::free(zre);
            if(zim_neq_0) std::free(zim);
            if(xre_neq_0) std::free(xre);
            if(xim_neq_0) std::free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: csinhv_svml_zmm16r4_u10x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{fname2.c_str()};
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

            xre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_xre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
               
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: csinhv_svml_zmm16r4_u10x_u\n");
            start  = __rdtscp(&tsc_aux);
            gms::math::csinhv_svml_zmm16r4_u10x_u(&xre[0],&xim[0],
                                                  &zre[0],&zim[0],nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: csinhv_svml_zmm16r4_u10x_u() -- TSC=%llu,TSC_AUX=%d\n",elapsed_tsc,tsc_aux);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(zre) std::free(zre);
                 if(zim) std::free(zim);
                 if(xre) std::free(xre);
                 if(xim) std::free(xim);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,j%.7f]\n",zre[__i],zim[__i]);
            }
            if(zre) std::free(zre);
            if(zim) std::free(zim);
            if(xre) std::free(xre);
            if(xim) std::free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: csinhv_svml_zmm16r4_u10x_u\n");
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
void unit_test_csinhv_svml_zmm16r4_u8x_u();

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_csinhv_svml_zmm16r4_u8x_u()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    constexpr unsigned long long RDTSCP_LAT{42ull};
    std::string fname1 = "UNIT-TEST_uni_distro_csinhv_svml_zmm16r4_u8x_u_" + std::to_string(funi_distro_idx++)+".csv";
    std::string fname2 = "UNIT-TEST_norm_distro_csinhv_svml_zmm16r4_u8x_u_" + std::to_string(fnorm_distr_idx++)+".csv";
    float * __restrict__ xre{NULL};
    float * __restrict__ xim{NULL};
    float * __restrict__ zre{NULL};
    float * __restrict__ zim{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_xre{0ull};
    unsigned long long seed_xim{0ull};
    unsigned long long start{0ull};
    unsigned long long end{0ull};
    unsigned long long start_c{0ull};
    unsigned long long end_c{0ull};
    unsigned long long elapsed_tsc{0ull};
    uint32_t tsc_aux{9999};
    float rxim,rxre;
    int32_t which{-1};
    int32_t idx{};
    bool  xre_neq_0{};
    bool  xim_neq_0{};
    bool  zre_neq_0{};
    bool  zim_neq_0{};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{fname1.c_str()};
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
            xre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xre_neq_0 = (xre != NULL) && (nbytes > 0ull);
            xim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xim_neq_0 = (xim != NULL) && (nbytes > 0ull);
            zre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zre_neq_0 = (zre != NULL) && (nbytes > 0ull);
            zim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zim_neq_0 = (zim != NULL) && (nbytes > 0ull);
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
               
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: csinhv_svml_zmm16r4_u8x_u\n");
            start   = __rdtscp(&tsc_aux);
            gms::math::csinhv_svml_zmm16r4_u8x_u(&xre[0],&xim[0],
                                                  &zre[0],&zim[0],nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: csinhv_svml_zmm16r4_u8x_u -- TSC=%llu,TSC_AUX=%d\n",elapsed_tsc,tsc_aux);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(zre_neq_0) std::free(zre);
                 if(zim_neq_0) std::free(zim);
                 if(xre_neq_0) std::free(xre);
                 if(xim_neq_0) std::free(xim);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,%.7fj]\n",zre[__i],zim[__i]);
            }
            if(zre_neq_0) std::free(zre);
            if(zim_neq_0) std::free(zim);
            if(xre_neq_0) std::free(xre);
            if(xim_neq_0) std::free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: csinhv_svml_zmm16r4_u8x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{fname2.c_str()};
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

            xre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_xre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: csinhv_svml_zmm16r4_u8x_u\n");
            start  = __rdtscp(&tsc_aux);
            gms::math::csinhv_svml_zmm16r4_u8x_u(&xre[0],&xim[0],
                                                  &zre[0],&zim[0],nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: csinhv_svml_zmm16r4_u8x_u() -- TSC=%llu,TSC_AUX=%d\n",elapsed_tsc,tsc_aux);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(zre) std::free(zre);
                 if(zim) std::free(zim);
                 if(xre) std::free(xre);
                 if(xim) std::free(xim);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,%.7fj]\n",zre[__i],zim[__i]);
            }
            if(zre) std::free(zre);
            if(zim) std::free(zim);
            if(xre) std::free(xre);
            if(xim) std::free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: csinhv_svml_zmm16r4_u8x_u\n");
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
void unit_test_csinhv_svml_zmm16r4_u6x_u();

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_csinhv_svml_zmm16r4_u6x_u()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    constexpr unsigned long long RDTSCP_LAT{42ull};
    std::string fname1 = "UNIT-TEST_uni_distro_csinhv_svml_zmm16r4_u6x_u_" + std::to_string(funi_distro_idx++)+".csv";
    std::string fname2 = "UNIT-TEST_norm_distro_csinhv_svml_zmm16r4_u6x_u_" + std::to_string(fnorm_distr_idx++)+".csv";
    float * __restrict__ xre{NULL};
    float * __restrict__ xim{NULL};
    float * __restrict__ zre{NULL};
    float * __restrict__ zim{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_xre{0ull};
    unsigned long long seed_xim{0ull};
    unsigned long long start{0ull};
    unsigned long long end{0ull};
    unsigned long long start_c{0ull};
    unsigned long long end_c{0ull};
    unsigned long long elapsed_tsc{0ull};
    uint32_t tsc_aux{9999};
    float rxim,rxre;
    int32_t which{-1};
    int32_t idx{};
    bool  xre_neq_0{};
    bool  xim_neq_0{};
    bool  zre_neq_0{};
    bool  zim_neq_0{};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{fname1.c_str()};
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
            xre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xre_neq_0 = (xre != NULL) && (nbytes > 0ull);
            xim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xim_neq_0 = (xim != NULL) && (nbytes > 0ull);
            zre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zre_neq_0 = (zre != NULL) && (nbytes > 0ull);
            zim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zim_neq_0 = (zim != NULL) && (nbytes > 0ull);
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_xre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
                
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
               
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: csinhv_svml_zmm16r4_u6x_u\n");
            start   = __rdtscp(&tsc_aux);
            gms::math::csinhv_svml_zmm16r4_u6x_u(&xre[0],&xim[0],
                                                  &zre[0],&zim[0],nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: csinhv_svml_zmm16r4_u6x_u -- TSC=%llu,TSC_AUX=%d\n",elapsed_tsc,tsc_aux);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(zre_neq_0) std::free(zre);
                 if(zim_neq_0) std::free(zim);
                 if(xre_neq_0) std::free(xre);
                 if(xim_neq_0) std::free(xim);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,%.7fj]\n",zre[__i],zim[__i]);
            }
            if(zre_neq_0) std::free(zre);
            if(zim_neq_0) std::free(zim);
            if(xre_neq_0) std::free(xre);
            if(xim_neq_0) std::free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: csinhv_svml_zmm16r4_u6x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{fname2.c_str()};
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

            xre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            xim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zre = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            zim = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_xim = __rdtsc();
            auto rand_xim{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_xim))};
            seed_xre = __rdtsc();
            auto rand_xre{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_xre))};
            int32_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rxim = rand_xim();
                xre[__i+0] = rxim;
                rxre = rand_xre();
                xim[__i+0] = rxre;
                rxim = rand_xim();
                xre[__i+1] = rxim;
                rxre = rand_xre();
                xim[__i+1] = rxre;
                rxim = rand_xim();
                xre[__i+2] = rxim;
                rxre = rand_xre();
                xim[__i+2] = rxre;
                rxim = rand_xim();
                xre[__i+3] = rxim;
                rxre = rand_xre();
                xim[__i+3] = rxre;
                rxim = rand_xim();
                xre[__i+4] = rxim;
                rxre = rand_xre();
                xim[__i+4] = rxre;
                rxim = rand_xim();
                xre[__i+5] = rxim;
                rxre = rand_xre();
                xim[__i+5] = rxre;
                rxim = rand_xim();
                xre[__i+6] = rxim;
                rxre = rand_xre();
                xim[__i+6] = rxre;
                rxim = rand_xim();
                xre[__i+7] = rxim;
                rxre = rand_xre();
                xim[__i+7] = rxre;
            }
            for(;__i != nelems; ++__i)
            {
                rxim = rand_xim();
                xre[__i] = rxim;
                rxre = rand_xre();
                xim[__i] = rxre;
                
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: csinhv_svml_zmm16r4_u6x_u\n");
            start  = __rdtscp(&tsc_aux);
            gms::math::csinhv_svml_zmm16r4_u6x_u(&xre[0],&xim[0],
                                                  &zre[0],&zim[0],nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: csinhv_svml_zmm16r4_u6x_u() -- TSC=%llu,TSC_AUX=%d\n",elapsed_tsc,tsc_aux);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(zre) std::free(zre);
                 if(zim) std::free(zim);
                 if(xre) std::free(xre);
                 if(xim) std::free(xim);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"[%.7f,%.7fj]\n",zre[__i],zim[__i]);
            }
            if(zre) std::free(zre);
            if(zim) std::free(zim);
            if(xre) std::free(xre);
            if(xim) std::free(xim);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: csinhv_svml_zmm16r4_u6x_u\n");
        }
        break;
        default : 
                 printf("Invalid switch variable=%d\n",which);
                 return;
    } 
}

int main()
{  
   warm_cache_L1I();
   unit_test_csinhv_svml_zmm16r4_u10x_u();
   unit_test_csinhv_svml_zmm16r4_u8x_u();
   unit_test_csinhv_svml_zmm16r4_u6x_u();
   return 0;
}