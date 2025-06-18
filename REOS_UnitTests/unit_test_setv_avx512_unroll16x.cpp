
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_setv_avx512_unroll16x.h"

/*
   icpc -o unit_test_setv_avx512_unroll16x -fp-model fast=2 -ftz -ggdb -ipo  -march=skylake-avx512 -mavx512f -qopt-zmm-usage=high -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_setv_avx512_unroll16x.h GMS_setv_avx512_unroll16x.cpp unit_test_setv_avx512_unroll16x.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -qopt-zmm-usage=high -falign-functions=32 GMS_config.h GMS_malloc.h GMS_setv_avx512_unroll16x.h GMS_setv_avx512_unroll16x.cpp unit_test_setv_avx512_unroll16x.cpp 
                        
*/

__attribute__((noinline))
__attribute__((aligned(32)))
__attribute__((hot))
void unit_test_ssetv_u_zmm16r4_unroll16x_incr_1();

void unit_test_ssetv_u_zmm16r4_unroll16x_incr_1()
{
     using namespace gms::math;
     constexpr int32_t size_1{1ULL};
     constexpr int32_t size_5{5ULL};
     constexpr int32_t size_16{16ULL};
     constexpr int32_t size_145{145ULL};
     constexpr int32_t size_489{489ULL};
     constexpr int32_t size_1024{1024ULL};
     constexpr int32_t size_4000{4000ULL};
          
     float data_4000[size_4000] = {};
     float data_1024[size_1024] = {};
     float data_489[size_489]   = {};
     float data_145[size_145]   = {};
     float data_16[size_16]     = {};
     float data_5[size_5]       = {};
     float data_1[size_1]       = {};
     float * __restrict__ pdata_4000{&data_4000[0]};
     float * __restrict__ pdata_1024{&data_1024[0]};
     float * __restrict__ pdata_489{&data_489[0]};
     float * __restrict__ pdata_145{&data_145[0]};
     float * __restrict__ pdata_16{&data_16[0]};
     float * __restrict__ pdata_5{&data_5[0]};
     float * __restrict__ pdata_1{&data_1[0]};

     const float fill{3.14159265358979323846264338328F};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s of increment: 1 ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %d with value=%.7f\n",size_4000,fill);
#if 0
     BREAK_INT3();
#endif
     ssetv_u_zmm16r4_unroll16x(size_4000,fill,&data_4000[0],1);
     for(int32_t __i{0}; __i != size_4000; ++__i) 
     {
        const float val{data_4000[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %d with value=%.7f\n",size_1024,fill);
     ssetv_u_zmm16r4_unroll16x(size_1024,fill,&data_1024[0],1);
     for(int32_t __i{0}; __i != size_1024; ++__i) 
     {
        const float val{data_1024[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %d with value=%.7f\n",size_489,fill);
     ssetv_u_zmm16r4_unroll16x(size_489,fill,&data_489[0],1);
     for(int32_t __i{0}; __i != size_489; ++__i) 
     {
        const float val{data_489[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %d with value=%.7f\n",size_145,fill);
     ssetv_u_zmm16r4_unroll16x(size_145,fill,&data_145[0],1);
     for(int32_t __i{0}; __i != size_145; ++__i) 
     {
        const float val{data_145[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %d with value=%.7f\n",size_16,fill);
     ssetv_u_zmm16r4_unroll16x(size_16,fill,&data_16[0],1);
     for(int32_t __i{0}; __i != size_16; ++__i) 
     {
        const float val{data_16[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %d with value=%.7f\n",size_5,fill);
     ssetv_u_zmm16r4_unroll16x(size_5,fill,&data_5[0],1);
     for(int32_t __i{0}; __i != size_5; ++__i) 
     {
        const float val{data_5[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %d with value=%.7f\n",size_1,fill);
     ssetv_u_zmm16r4_unroll16x(size_1,fill,&data_1[0],1);
     const float val{data_1[0ULL]};
     if(val==0.0f)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,0);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() of increment 1: ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}


__attribute__((noinline))
__attribute__((aligned(32)))
__attribute__((hot))
void unit_test_ssetv_u_zmm16r4_unroll16x_incr_3();

void unit_test_ssetv_u_zmm16r4_unroll16x_incr_3()
{
     using namespace gms::math;
     constexpr int32_t size_1{1};
     constexpr int32_t size_5{5};
     constexpr int32_t size_16{16};
     constexpr int32_t size_145{145};
     constexpr int32_t size_489{489};
     constexpr int32_t size_1024{1024};
     constexpr int32_t size_4000{4000};
          
     float data_4000[size_4000] = {};
     float data_1024[size_1024] = {};
     float data_489[size_489]   = {};
     float data_145[size_145]   = {};
     float data_16[size_16]     = {};
     float data_5[size_5]       = {};
     float data_1[size_1]       = {};
     float * __restrict__ pdata_4000{&data_4000[0]};
     float * __restrict__ pdata_1024{&data_1024[0]};
     float * __restrict__ pdata_489{&data_489[0]};
     float * __restrict__ pdata_145{&data_145[0]};
     float * __restrict__ pdata_16{&data_16[0]};
     float * __restrict__ pdata_5{&data_5[0]};
     float * __restrict__ pdata_1{&data_1[0]};

     const float fill{3.14159265358979323846264338328F};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s of increment: 3 ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST]: -- Warning: This test may be failing and this is expected state.\n");
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %d with value=%.7f\n",size_4000,fill);

     ssetv_u_zmm16r4_unroll16x(size_4000,fill,&data_4000[0],3);
     for(int32_t __i{0}; __i != size_4000; ++__i) 
     {
        const float val{data_4000[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %d with value=%.7f\n",size_1024,fill);
     ssetv_u_zmm16r4_unroll16x(size_1024,fill,&data_1024[0],3);
     for(int32_t __i{0}; __i != size_1024; ++__i) 
     {
        const float val{data_1024[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %d with value=%.7f\n",size_489,fill);
     ssetv_u_zmm16r4_unroll16x(size_489,fill,&data_489[0],3);
     for(int32_t __i{0}; __i != size_489; ++__i) 
     {
        const float val{data_489[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %d with value=%.7f\n",size_145,fill);
     ssetv_u_zmm16r4_unroll16x(size_145,fill,&data_145[0],3);
     for(int32_t __i{0}; __i != size_145; ++__i) 
     {
        const float val{data_145[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %d with value=%.7f\n",size_16,fill);
     ssetv_u_zmm16r4_unroll16x(size_16,fill,&data_16[0],3);
     for(int32_t __i{0}; __i != size_16; ++__i) 
     {
        const float val{data_16[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %d with value=%.7f\n",size_5,fill);
     ssetv_u_zmm16r4_unroll16x(size_5,fill,&data_5[0],3);
     for(int32_t __i{0}; __i != size_5; ++__i) 
     {
        const float val{data_5[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %d with value=%.7f\n",size_1,fill);
     ssetv_u_zmm16r4_unroll16x(size_1,fill,&data_1[0],3);
     const float val{data_1[0ULL]};
     if(val==0.0f)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,0);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() of increment 3: ---> ENDED\n", __PRETTY_FUNCTION__);
}

__attribute__((noinline))
__attribute__((aligned(32)))
__attribute__((hot))
void test_runner_ssetv_u_zmm16r4_unroll16x();

void test_runner_ssetv_u_zmm16r4_unroll16x()
{
     std::clock_t seed{0ULL};
     int32_t which{-1};
     seed = std::clock();
     auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,1),
                          std::mt19937(seed))};
     which = rand();
     switch (which)
     {
         case 0 : 
                   unit_test_ssetv_u_zmm16r4_unroll16x_incr_1();
         break;
         case 1 :
                   unit_test_ssetv_u_zmm16r4_unroll16x_incr_3();
         break;
         default : 
                   printf("[TEST-RUNNER]: %s -- invalid switch variable=%d\n",__PRETTY_FUNCTION__,which);
                   return;
     }
}

__attribute__((noinline))
__attribute__((aligned(32)))
__attribute__((hot))
void unit_test_ssetv_a_zmm16r4_unroll16x_incr_1();

void unit_test_ssetv_a_zmm16r4_unroll16x_incr_1()
{
     using namespace gms::math;
     constexpr int32_t size_1{1ULL};
     constexpr int32_t size_5{5ULL};
     constexpr int32_t size_16{16ULL};
     constexpr int32_t size_145{145ULL};
     constexpr int32_t size_489{489ULL};
     constexpr int32_t size_1024{1024ULL};
     constexpr int32_t size_4000{4000ULL};
          
     float __attribute__((aligned(64))) data_4000[size_4000] = {};
     float __attribute__((aligned(64))) data_1024[size_1024] = {};
     float __attribute__((aligned(64))) data_489[size_489]   = {};
     float __attribute__((aligned(64))) data_145[size_145]   = {};
     float __attribute__((aligned(64))) data_16[size_16]     = {};
     float __attribute__((aligned(64))) data_5[size_5]       = {};
     float __attribute__((aligned(64))) data_1[size_1]       = {};
     float * __restrict__ pdata_4000{&data_4000[0]};
     float * __restrict__ pdata_1024{&data_1024[0]};
     float * __restrict__ pdata_489{&data_489[0]};
     float * __restrict__ pdata_145{&data_145[0]};
     float * __restrict__ pdata_16{&data_16[0]};
     float * __restrict__ pdata_5{&data_5[0]};
     float * __restrict__ pdata_1{&data_1[0]};

     const float fill{3.14159265358979323846264338328F};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s of increment: 1 ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %d with value=%.7f\n",size_4000,fill);
#if 0
     BREAK_INT3();
#endif
     ssetv_a_zmm16r4_unroll16x(size_4000,fill,&data_4000[0],1);
     for(int32_t __i{0}; __i != size_4000; ++__i) 
     {
        const float val{data_4000[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %d with value=%.7f\n",size_1024,fill);
     ssetv_a_zmm16r4_unroll16x(size_1024,fill,&data_1024[0],1);
     for(int32_t __i{0}; __i != size_1024; ++__i) 
     {
        const float val{data_1024[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %d with value=%.7f\n",size_489,fill);
     ssetv_a_zmm16r4_unroll16x(size_489,fill,&data_489[0],1);
     for(int32_t __i{0}; __i != size_489; ++__i) 
     {
        const float val{data_489[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %d with value=%.7f\n",size_145,fill);
     ssetv_a_zmm16r4_unroll16x(size_145,fill,&data_145[0],1);
     for(int32_t __i{0}; __i != size_145; ++__i) 
     {
        const float val{data_145[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %d with value=%.7f\n",size_16,fill);
     ssetv_a_zmm16r4_unroll16x(size_16,fill,&data_16[0],1);
     for(int32_t __i{0}; __i != size_16; ++__i) 
     {
        const float val{data_16[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %d with value=%.7f\n",size_5,fill);
     ssetv_a_zmm16r4_unroll16x(size_5,fill,&data_5[0],1);
     for(int32_t __i{0}; __i != size_5; ++__i) 
     {
        const float val{data_5[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %d with value=%.7f\n",size_1,fill);
     ssetv_a_zmm16r4_unroll16x(size_1,fill,&data_1[0],1);
     const float val{data_1[0ULL]};
     if(val==0.0f)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,0);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() of increment 1: ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}


__attribute__((noinline))
__attribute__((aligned(32)))
__attribute__((hot))
void unit_test_ssetv_a_zmm16r4_unroll16x_incr_3();

void unit_test_ssetv_a_zmm16r4_unroll16x_incr_3()
{
     using namespace gms::math;
     constexpr int32_t size_1{1};
     constexpr int32_t size_5{5};
     constexpr int32_t size_16{16};
     constexpr int32_t size_145{145};
     constexpr int32_t size_489{489};
     constexpr int32_t size_1024{1024};
     constexpr int32_t size_4000{4000};
          
     float __attribute__((aligned(64))) data_4000[size_4000] = {};
     float __attribute__((aligned(64))) data_1024[size_1024] = {};
     float __attribute__((aligned(64))) data_489[size_489]   = {};
     float __attribute__((aligned(64))) data_145[size_145]   = {};
     float __attribute__((aligned(64))) data_16[size_16]     = {};
     float __attribute__((aligned(64))) data_5[size_5]       = {};
     float __attribute__((aligned(64))) data_1[size_1]       = {};
     float * __restrict__ pdata_4000{&data_4000[0]};
     float * __restrict__ pdata_1024{&data_1024[0]};
     float * __restrict__ pdata_489{&data_489[0]};
     float * __restrict__ pdata_145{&data_145[0]};
     float * __restrict__ pdata_16{&data_16[0]};
     float * __restrict__ pdata_5{&data_5[0]};
     float * __restrict__ pdata_1{&data_1[0]};

     const float fill{3.14159265358979323846264338328F};
     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     printf("[UNIT-TEST]: %s of increment: 3 ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST]: -- Warning: This test may be failing and this is expected state.\n");
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %d with value=%.7f\n",size_4000,fill);
#if 0
     BREAK_INT3();
#endif
     ssetv_a_zmm16r4_unroll16x(size_4000,fill,&data_4000[0],3);
     for(int32_t __i{0}; __i != size_4000; ++__i) 
     {
        const float val{data_4000[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b4000_fail = true;
            break;
        }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %d with value=%.7f\n",size_1024,fill);
     ssetv_a_zmm16r4_unroll16x(size_1024,fill,&data_1024[0],3);
     for(int32_t __i{0}; __i != size_1024; ++__i) 
     {
        const float val{data_1024[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b1024_fail = true;
            break;
        }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %d with value=%.7f\n",size_489,fill);
     ssetv_a_zmm16r4_unroll16x(size_489,fill,&data_489[0],3);
     for(int32_t __i{0}; __i != size_489; ++__i) 
     {
        const float val{data_489[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b489_fail = true;
            break;
        }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#3]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#4]: -- fill buffer of size: %d with value=%.7f\n",size_145,fill);
     ssetv_a_zmm16r4_unroll16x(size_145,fill,&data_145[0],3);
     for(int32_t __i{0}; __i != size_145; ++__i) 
     {
        const float val{data_145[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b145_fail = true;
            break;
        }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#5]: -- fill buffer of size: %d with value=%.7f\n",size_16,fill);
     ssetv_a_zmm16r4_unroll16x(size_16,fill,&data_16[0],3);
     for(int32_t __i{0}; __i != size_16; ++__i) 
     {
        const float val{data_16[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b16_fail = true;
            break;
        }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#6]: -- fill buffer of size: %d with value=%.7f\n",size_5,fill);
     ssetv_a_zmm16r4_unroll16x(size_5,fill,&data_5[0],3);
     for(int32_t __i{0}; __i != size_5; ++__i) 
     {
        const float val{data_5[__i]};
        if(val==0.0f)
        {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,__i);
            b5_fail = true;
            break;
        }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}

     printf("[UNIT-TEST--#7]: -- fill buffer of size: %d with value=%.7f\n",size_1,fill);
     ssetv_a_zmm16r4_unroll16x(size_1,fill,&data_1[0],3);
     const float val{data_1[0ULL]};
     if(val==0.0f)
     {
            printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %d" ANSI_RESET_ALL "\n",val,0);
            b1_fail = true;
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s() of increment 3: ---> ENDED\n", __PRETTY_FUNCTION__);
}

__attribute__((noinline))
__attribute__((aligned(32)))
__attribute__((hot))
void test_runner_ssetv_a_zmm16r4_unroll16x();

void test_runner_ssetv_a_zmm16r4_unroll16x()
{
     std::clock_t seed{0ULL};
     int32_t which{-1};
     seed = std::clock();
     auto rand{std::bind(std::uniform_int_distribution<int32_t>(0,1),
                          std::mt19937(seed))};
     which = rand();
     switch (which)
     {
         case 0 : 
                   unit_test_ssetv_a_zmm16r4_unroll16x_incr_1();
         break;
         case 1 :
                   unit_test_ssetv_a_zmm16r4_unroll16x_incr_3();
         break;
         default : 
                   printf("[TEST-RUNNER]: %s -- invalid switch variable=%d\n",__PRETTY_FUNCTION__,which);
                   return;
     }
}

int main()
{
    test_runner_ssetv_u_zmm16r4_unroll16x();
    test_runner_ssetv_a_zmm16r4_unroll16x();
    return 0;
}

