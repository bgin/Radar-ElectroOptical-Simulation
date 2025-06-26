
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h>
#include "GMS_malloc.h"
#include "GMS_cabs_vec_zmm16r4.h"

/*
   icpc -o unit_test_cabs_vec_zmm16r4_u -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_cabs_vec_zmm16r4.h GMS_cabs_vec_zmm16r4.cpp unit_test_cabs_vec_zmm16r4_u.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_cabs_vec_zmm16r4.h GMS_cabs_vec_zmm16r4.cpp unit_test_cabs_vec_zmm16r4_u.cpp
*/

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_cabsv_zmm16r4_unroll_16x_u();

void unit_test_cabsv_zmm16r4_unroll_16x_u()
{
    using namespace gms::common;
    float * __restrict__ p_re{NULL};
    float * __restrict__ p_im{NULL};
    float * __restrict__ p_cabsv{NULL};
    FILE  *  __restrict__ fp{NULL};
    
    constexpr std::size_t nelems{67456};
    constexpr std::size_t nbytes{4ull*nelems};
    std::clock_t seed_d{0ull};
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
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
            const char * __restrict__ fout{"UNIT-TEST_uni_distro_cabsv_zmm16r4_unroll_16x_u.csv"};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(uniform_distro).name(),status);
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
            
            p_re    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_im    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_cabsv = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_re = __rdtsc();
            auto rand_re{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                  std::mt19937(seed_re))};
            seed_im = __rdtsc();
            auto rand_im{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                   std::mt19937(seed_im))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float re0{rand_re()};
                 p_re[__i+0] = re0;
                 float im0{rand_im()};
                 p_im[__i+0] = im0;
                 float re1{rand_re()};
                 p_re[__i+1] = re1;
                 float im1{rand_im()};
                 p_im[__i+1] = im1;
                 float re2{rand_re()};
                 p_re[__i+2] = re2;
                 float im2{rand_im()};
                 p_im[__i+2] = im2;
                 float re3{rand_re()};
                 p_re[__i+3] = re3;
                 float im3{rand_im()};
                 p_im[__i+3] = im3;
                 float re4{rand_re()};
                 p_re[__i+4] = re4;
                 float im4{rand_im()};
                 p_im[__i+4] = im4;
                 float re5{rand_re()};
                 p_re[__i+5] = re5;
                 float im5{rand_im()};
                 p_im[__i+5] = im5;
                 float re6{rand_re()};
                 p_re[__i+6] = re6;
                 float im6{rand_im()};
                 p_im[__i+6] = im6;
                 float re7{rand_re()};
                 p_re[__i+7] = re7;
                 float im7{rand_im()};
                 p_im[__i+7] = im7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cabsv_zmm16r4_unroll_16x_u\n");
            gms::math::cabsv_zmm16r4_unroll_16x_u(&p_re[0],&p_im[0],&p_cabsv[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_cabsv[%d]=%.7f\n",__i,p_cabsv[__i]);
            }
            if(p_re)    std::free(p_re);
            if(p_im)    std::free(p_im);
            if(p_cabsv) std::free(p_cabsv);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cabsv_zmm16r4_unroll_16x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{"UNIT-TEST_norm_distro_cabsv_zmm16r4_unroll_16x_u.csv"};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
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
            
            p_re    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_im    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_cabsv = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_re = __rdtsc();
            auto rand_re{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                  std::mt19937(seed_re))};
            seed_im = __rdtsc();
            auto rand_im{std::bind(std::normal_distribution<float>(0.0f,2.0f),
                                  std::mt19937(seed_im))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float re0{rand_re()};
                 p_re[__i+0] = re0;
                 float im0{rand_im()};
                 p_im[__i+0] = im0;
                 float re1{rand_re()};
                 p_re[__i+1] = re1;
                 float im1{rand_im()};
                 p_im[__i+1] = im1;
                 float re2{rand_re()};
                 p_re[__i+2] = re2;
                 float im2{rand_im()};
                 p_im[__i+2] = im2;
                 float re3{rand_re()};
                 p_re[__i+3] = re3;
                 float im3{rand_im()};
                 p_im[__i+3] = im3;
                 float re4{rand_re()};
                 p_re[__i+4] = re4;
                 float im4{rand_im()};
                 p_im[__i+4] = im4;
                 float re5{rand_re()};
                 p_re[__i+5] = re5;
                 float im5{rand_im()};
                 p_im[__i+5] = im5;
                 float re6{rand_re()};
                 p_re[__i+6] = re6;
                 float im6{rand_im()};
                 p_im[__i+6] = im6;
                 float re7{rand_re()};
                 p_re[__i+7] = re7;
                 float im7{rand_im()};
                 p_im[__i+7] = im7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cabsv_zmm16r4_unroll_16x_u\n");
            gms::math::cabsv_zmm16r4_unroll_16x_u(&p_re[0],&p_im[0],&p_cabsv[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_cabsv[%d]=%.7f\n",__i,p_cabsv[__i]);
            }
            if(p_re)     std::free(p_re);
            if(p_im)     std::free(p_im);
            if(p_cabsv)  std::free(p_cabsv);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cabsv_zmm16r4_unroll_16x_u\n");
        }
    }
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_cabsv_zmm16r4_unroll_10x_u();

void unit_test_cabsv_zmm16r4_unroll_10x_u()
{
    using namespace gms::common;
    float * __restrict__ p_re{NULL};
    float * __restrict__ p_im{NULL};
    float * __restrict__ p_cabsv{NULL};
    FILE  *  __restrict__ fp{NULL};
    
    constexpr std::size_t nelems{67456};
    constexpr std::size_t nbytes{4ull*nelems};
    std::clock_t seed_d{0ull};
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
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
            const char * __restrict__ fout{"UNIT-TEST_uni_distro_cabsv_zmm16r4_unroll_10x_u.csv"};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(uniform_distro).name(),status);
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
            
            p_re    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_im    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_cabsv = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_re = __rdtsc();
            auto rand_re{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                  std::mt19937(seed_re))};
            seed_im = __rdtsc();
            auto rand_im{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                   std::mt19937(seed_im))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float re0{rand_re()};
                 p_re[__i+0] = re0;
                 float im0{rand_im()};
                 p_im[__i+0] = im0;
                 float re1{rand_re()};
                 p_re[__i+1] = re1;
                 float im1{rand_im()};
                 p_im[__i+1] = im1;
                 float re2{rand_re()};
                 p_re[__i+2] = re2;
                 float im2{rand_im()};
                 p_im[__i+2] = im2;
                 float re3{rand_re()};
                 p_re[__i+3] = re3;
                 float im3{rand_im()};
                 p_im[__i+3] = im3;
                 float re4{rand_re()};
                 p_re[__i+4] = re4;
                 float im4{rand_im()};
                 p_im[__i+4] = im4;
                 float re5{rand_re()};
                 p_re[__i+5] = re5;
                 float im5{rand_im()};
                 p_im[__i+5] = im5;
                 float re6{rand_re()};
                 p_re[__i+6] = re6;
                 float im6{rand_im()};
                 p_im[__i+6] = im6;
                 float re7{rand_re()};
                 p_re[__i+7] = re7;
                 float im7{rand_im()};
                 p_im[__i+7] = im7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cabsv_zmm16r4_unroll_10x_u\n");
            gms::math::cabsv_zmm16r4_unroll_10x_u(&p_re[0],&p_im[0],&p_cabsv[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_cabsv[%d]=%.7f\n",__i,p_cabsv[__i]);
            }
            if(p_re)    std::free(p_re);
            if(p_im)    std::free(p_im);
            if(p_cabsv) std::free(p_cabsv);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cabsv_zmm16r4_unroll_10x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{"UNIT-TEST_norm_distro_cabsv_zmm16r4_unroll_10x_u.csv"};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
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
            
            p_re    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_im    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_cabsv = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_re = __rdtsc();
            auto rand_re{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                  std::mt19937(seed_re))};
            seed_im = __rdtsc();
            auto rand_im{std::bind(std::normal_distribution<float>(0.0f,2.0f),
                                  std::mt19937(seed_im))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float re0{rand_re()};
                 p_re[__i+0] = re0;
                 float im0{rand_im()};
                 p_im[__i+0] = im0;
                 float re1{rand_re()};
                 p_re[__i+1] = re1;
                 float im1{rand_im()};
                 p_im[__i+1] = im1;
                 float re2{rand_re()};
                 p_re[__i+2] = re2;
                 float im2{rand_im()};
                 p_im[__i+2] = im2;
                 float re3{rand_re()};
                 p_re[__i+3] = re3;
                 float im3{rand_im()};
                 p_im[__i+3] = im3;
                 float re4{rand_re()};
                 p_re[__i+4] = re4;
                 float im4{rand_im()};
                 p_im[__i+4] = im4;
                 float re5{rand_re()};
                 p_re[__i+5] = re5;
                 float im5{rand_im()};
                 p_im[__i+5] = im5;
                 float re6{rand_re()};
                 p_re[__i+6] = re6;
                 float im6{rand_im()};
                 p_im[__i+6] = im6;
                 float re7{rand_re()};
                 p_re[__i+7] = re7;
                 float im7{rand_im()};
                 p_im[__i+7] = im7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cabsv_zmm16r4_unroll_10x_u\n");
            gms::math::cabsv_zmm16r4_unroll_10x_u(&p_re[0],&p_im[0],&p_cabsv[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_cabsv[%d]=%.7f\n",__i,p_cabsv[__i]);
            }
            if(p_re)     std::free(p_re);
            if(p_im)     std::free(p_im);
            if(p_cabsv)  std::free(p_cabsv);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cabsv_zmm16r4_unroll_10x_u\n");
        }
    }
}

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_cabsv_zmm16r4_unroll_6x_u();

void unit_test_cabsv_zmm16r4_unroll_6x_u()
{
    using namespace gms::common;
    float * __restrict__ p_re{NULL};
    float * __restrict__ p_im{NULL};
    float * __restrict__ p_cabsv{NULL};
    FILE  *  __restrict__ fp{NULL};
    
    constexpr std::size_t nelems{67456};
    constexpr std::size_t nbytes{4ull*nelems};
    std::clock_t seed_d{0ull};
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
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
            const char * __restrict__ fout{"UNIT-TEST_uni_distro_cabsv_zmm16r4_unroll_6x_u.csv"};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(uniform_distro).name(),status);
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
            
            p_re    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_im    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_cabsv = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_re = __rdtsc();
            auto rand_re{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                  std::mt19937(seed_re))};
            seed_im = __rdtsc();
            auto rand_im{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                   std::mt19937(seed_im))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float re0{rand_re()};
                 p_re[__i+0] = re0;
                 float im0{rand_im()};
                 p_im[__i+0] = im0;
                 float re1{rand_re()};
                 p_re[__i+1] = re1;
                 float im1{rand_im()};
                 p_im[__i+1] = im1;
                 float re2{rand_re()};
                 p_re[__i+2] = re2;
                 float im2{rand_im()};
                 p_im[__i+2] = im2;
                 float re3{rand_re()};
                 p_re[__i+3] = re3;
                 float im3{rand_im()};
                 p_im[__i+3] = im3;
                 float re4{rand_re()};
                 p_re[__i+4] = re4;
                 float im4{rand_im()};
                 p_im[__i+4] = im4;
                 float re5{rand_re()};
                 p_re[__i+5] = re5;
                 float im5{rand_im()};
                 p_im[__i+5] = im5;
                 float re6{rand_re()};
                 p_re[__i+6] = re6;
                 float im6{rand_im()};
                 p_im[__i+6] = im6;
                 float re7{rand_re()};
                 p_re[__i+7] = re7;
                 float im7{rand_im()};
                 p_im[__i+7] = im7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cabsv_zmm16r4_unroll_6x_u\n");
            gms::math::cabsv_zmm16r4_unroll_6x_u(&p_re[0],&p_im[0],&p_cabsv[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_cabsv[%d]=%.7f\n",__i,p_cabsv[__i]);
            }
            if(p_re)    std::free(p_re);
            if(p_im)    std::free(p_im);
            if(p_cabsv) std::free(p_cabsv);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cabsv_zmm16r4_unroll_6x_u\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{"UNIT-TEST_norm_distro_cabsv_zmm16r4_unroll_6x_u.csv"};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
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
            
            p_re    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_im    = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            p_cabsv = reinterpret_cast<float * __restrict__>(std::malloc(nbytes));
            seed_re = __rdtsc();
            auto rand_re{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                  std::mt19937(seed_re))};
            seed_im = __rdtsc();
            auto rand_im{std::bind(std::normal_distribution<float>(0.0f,2.0f),
                                  std::mt19937(seed_im))};
            for(int32_t __i{0}; __i != nelems; __i += 8)
            {
                 float re0{rand_re()};
                 p_re[__i+0] = re0;
                 float im0{rand_im()};
                 p_im[__i+0] = im0;
                 float re1{rand_re()};
                 p_re[__i+1] = re1;
                 float im1{rand_im()};
                 p_im[__i+1] = im1;
                 float re2{rand_re()};
                 p_re[__i+2] = re2;
                 float im2{rand_im()};
                 p_im[__i+2] = im2;
                 float re3{rand_re()};
                 p_re[__i+3] = re3;
                 float im3{rand_im()};
                 p_im[__i+3] = im3;
                 float re4{rand_re()};
                 p_re[__i+4] = re4;
                 float im4{rand_im()};
                 p_im[__i+4] = im4;
                 float re5{rand_re()};
                 p_re[__i+5] = re5;
                 float im5{rand_im()};
                 p_im[__i+5] = im5;
                 float re6{rand_re()};
                 p_re[__i+6] = re6;
                 float im6{rand_im()};
                 p_im[__i+6] = im6;
                 float re7{rand_re()};
                 p_re[__i+7] = re7;
                 float im7{rand_im()};
                 p_im[__i+7] = im7;
            }
            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: cabsv_zmm16r4_unroll_6x_u\n");
            gms::math::cabsv_zmm16r4_unroll_6x_u(&p_re[0],&p_im[0],&p_cabsv[0],nelems);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"p_cabsv[%d]=%.7f\n",__i,p_cabsv[__i]);
            }
            if(p_re)     std::free(p_re);
            if(p_im)     std::free(p_im);
            if(p_cabsv)  std::free(p_cabsv);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: cabsv_zmm16r4_unroll_6x_u\n");
        }
    }
}


int main()
{
   unit_test_cabsv_zmm16r4_unroll_16x_u();
   unit_test_cabsv_zmm16r4_unroll_10x_u();
   unit_test_cabsv_zmm16r4_unroll_6x_u();
   return 0;

}
