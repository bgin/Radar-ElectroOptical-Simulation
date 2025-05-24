
#include <stdio.h>
#include <cstdlib>
#include <algorithm>
#include "GMS_sse_memset.h"

/*
   icpc -o unit_test_sse_memset -fp-model -fno-exceptions fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_sse_memset.h GMS_sse_memset.cpp unit_test_sse_memset.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 GMS_config.h  GMS_sse_memset.h GMS_sse_memset.cpp unit_test_sse_memset.cpp

*/
//#define ANSI_COLOR_RED          "\x1b[31m"
//#define ANSI_COLOR_GREEN        "\x1b[32m"
//#define ANSI_RESET_ALL          "\x1b[0m"

void unit_test_sse_memset_unroll8x_ps();

void unit_test_sse_memset_unroll8x_ps()
{
     using namespace gms::common;
     constexpr std::size_t size_1{1ULL};
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     
     float data_4000[size_4000] = {};
     float data_1024[size_1024] = {};
     float data_489[size_489]   = {};
     float data_145[size_145]   = {};
     float data_16[size_16]     = {};
     float data_5[size_5]       = {};
     float data_1[size_1]       = {};
     const float fill{3.14159265358979323846264338328F};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s() ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.7f\n",size_4000,fill);
     sse_memset_unroll8x_ps(&data_4000[0],fill,size_4000);
     for(std::size_t __i{0ULL}; __i != size_4000; ++__i) 
     {
        const float val{data_4000[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.7f\n",size_1024,fill);
     sse_memset_unroll8x_ps(&data_1024[0],fill,size_1024);
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i) 
     {
        const float val{data_1024[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.7f\n",size_489,fill);
     sse_memset_unroll8x_ps(&data_489[0],fill,size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i) 
     {
        const float val{data_489[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.7f\n",size_145,fill);
     sse_memset_unroll8x_ps(&data_145[0],fill,size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i) 
     {
        const float val{data_145[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.7f\n",size_16,fill);
     sse_memset_unroll8x_ps(&data_16[0],fill,size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i) 
     {
        const float val{data_16[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.7f\n",size_5,fill);
     sse_memset_unroll8x_ps(&data_5[0],fill,size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i) 
     {
        const float val{data_5[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %llu with value=%.7f\n",size_1,fill);
     sse_memset_unroll8x_ps(&data_1[0],fill,size_1);
     const float val{data_1[0ULL]};
     if(val==0.0f)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,0ULL);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() ---> ENDED\n", __PRETTY_FUNCTION__);


}


void unit_test_sse_memset_unroll8x_pd();

void unit_test_sse_memset_unroll8x_pd()
{
     using namespace gms::common;
     constexpr std::size_t size_1{1ULL};
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     
     double data_4000[size_4000] = {};
     double data_1024[size_1024] = {};
     double data_489[size_489]   = {};
     double data_145[size_145]   = {};
     double data_16[size_16]     = {};
     double data_5[size_5]       = {};
     double data_1[size_1]       = {};
     const  double fill{3.14159265358979323846264338328};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s() ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.16f\n",size_4000,fill);
     sse_memset_unroll8x_pd(&data_4000[0],fill,size_4000);
     for(std::size_t __i{0ULL}; __i != size_4000; ++__i) 
     {
        const double val{data_4000[__i]};
        if(val==0.0)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.16f\n",size_1024,fill);
     sse_memset_unroll8x_pd(&data_1024[0],fill,size_1024);
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i) 
     {
        const double val{data_1024[__i]};
        if(val==0.0)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.16f\n",size_489,fill);
     sse_memset_unroll8x_pd(&data_489[0],fill,size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i) 
     {
        const double val{data_489[__i]};
        if(val==0.0)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.16f\n",size_145,fill);
     sse_memset_unroll8x_pd(&data_145[0],fill,size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i) 
     {
        const double val{data_145[__i]};
        if(val==0.0)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.16f\n",size_16,fill);
     sse_memset_unroll8x_pd(&data_16[0],fill,size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i) 
     {
        const double val{data_16[__i]};
        if(val==0.0)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.16f\n",size_5,fill);
     sse_memset_unroll8x_pd(&data_5[0],fill,size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i) 
     {
        const double val{data_5[__i]};
        if(val==0.0)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %llu with value=%.16f\n",size_1,fill);
     sse_memset_unroll8x_pd(&data_1[0],fill,size_1);
     const double val{data_1[0ULL]};
     if(val==0.0)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.16f at pos: %llu" ANSI_RESET_ALL "\n",val,0ULL);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() ---> ENDED\n", __PRETTY_FUNCTION__);
}

void unit_test_sse_memset_unroll16x_ps();

void unit_test_sse_memset_unroll16x_ps()
{
      using namespace gms::common;
     constexpr std::size_t size_1{1ULL};
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     
     float data_4000[size_4000] = {};
     float data_1024[size_1024] = {};
     float data_489[size_489]   = {};
     float data_145[size_145]   = {};
     float data_16[size_16]     = {};
     float data_5[size_5]       = {};
     float data_1[size_1]       = {};
     const float fill{3.14159265358979323846264338328F};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s() ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.7f\n",size_4000,fill);
     sse_memset_unroll16x_ps(&data_4000[0],fill,size_4000);
     for(std::size_t __i{0ULL}; __i != size_4000; ++__i) 
     {
        const float val{data_4000[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.7f\n",size_1024,fill);
     sse_memset_unroll16x_ps(&data_1024[0],fill,size_1024);
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i) 
     {
        const float val{data_1024[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.7f\n",size_489,fill);
     sse_memset_unroll16x_ps(&data_489[0],fill,size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i) 
     {
        const float val{data_489[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.7f\n",size_145,fill);
     sse_memset_unroll16x_ps(&data_145[0],fill,size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i) 
     {
        const float val{data_145[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.7f\n",size_16,fill);
     sse_memset_unroll16x_ps(&data_16[0],fill,size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i) 
     {
        const float val{data_16[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.7f\n",size_5,fill);
     sse_memset_unroll16x_ps(&data_5[0],fill,size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i) 
     {
        const float val{data_5[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %llu with value=%.7f\n",size_1,fill);
     sse_memset_unroll16x_ps(&data_1[0],fill,size_1);
     const float val{data_1[0ULL]};
     if(val==0.0f)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %llu" ANSI_RESET_ALL "\n",val,0ULL);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() ---> ENDED\n", __PRETTY_FUNCTION__);
}



int main()
{
    unit_test_sse_memset_unroll8x_ps();
    unit_test_sse_memset_unroll8x_pd();
    unit_test_sse_memset_unroll16x_ps();
    return 0;
}