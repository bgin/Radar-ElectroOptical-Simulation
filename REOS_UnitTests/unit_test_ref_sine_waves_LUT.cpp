
#include <cstdio>
#include <cstdint>
#include "GMS_ref_sine_waves_LUT.h"

/*
   icpc -o unit_test_ref_sine_waves_LUT -fp-model fast=2 -std=c++17  -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_ref_sine_waves_LUT.h GMS_ref_sine_waves_LUT.cpp unit_test_ref_sine_waves_LUT.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_ref_sine_waves_LUT.h GMS_ref_sine_waves_LUT.cpp unit_test_ref_sine_waves_LUT.cpp

*/

void unit_test_ref_sine_waves_LUT();

void unit_test_ref_sine_waves_LUT()
{
     std::printf("[UNIT-TEST#1] -- Accessing a LUT: LUT_sin_1_1_32\n");
     for(int32_t __i{0}; __i != 32; ++__i)
     {
         std::printf("%.7f\n",gms::math::LUT_sin_1_1_32[__i]);
     }
}


int main()
{
    unit_test_ref_sine_waves_LUT();
    return 0;
}