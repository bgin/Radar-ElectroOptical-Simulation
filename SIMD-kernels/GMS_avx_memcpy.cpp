

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!*/


#include <immintrin.h>
#include "GMS_avx_memcpy.h"

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
void 
gms::common::
avx_memcpy_unroll8x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
{
     float * __restrict__ p_dst{dst};
     float * __restrict__ p_src{src};

     while(((uintptr_t)p_dst & 31) && sz)
     {
            const float t{*p_src};
            *p_dst = t;
             p_src++;
             p_dst++;
             sz--;
     }

     while(sz >= 64ull)
     {
           __m256 ymm0;
           __m256 ymm1;
           __m256 ymm2;
           __m256 ymm3;
           __m256 ymm4;
           __m256 ymm5;
           __m256 ymm6;
           __m256 ymm7;


#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
#endif
            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
#endif
            ymm4 = _mm256_load_ps(p_src+32ull);
            _mm256_store_ps(p_dst+32ull,ymm4);
            ymm5 = _mm256_load_ps(p_src+40ull);
            _mm256_store_ps(p_dst+40ull,ymm5);
            ymm6 = _mm256_load_ps(p_src+48ull);
            _mm256_store_ps(p_dst+48ull,ymm6);
            ymm7 = _mm256_load_ps(p_src+56ull);
            _mm256_store_ps(p_dst+56ull,ymm7);
#else 
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
#endif
            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
#endif
            ymm4 = _mm256_load_ps(p_src+32ull);
            ymm5 = _mm256_load_ps(p_src+40ull);
            ymm6 = _mm256_load_ps(p_src+48ull);
            ymm7 = _mm256_load_ps(p_src+56ull);

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
            _mm256_store_ps(p_dst+32ull,ymm4);
            _mm256_store_ps(p_dst+40ull,ymm5);
            _mm256_store_ps(p_dst+48ull,ymm6);
            _mm256_store_ps(p_dst+56ull,ymm7);
#endif 
            sz    -= 64ull;
            p_src += 64ull;
            p_dst += 64ull;
    
     }

     while(sz >= 48ull)
     {
           __m256 ymm0;
           __m256 ymm1;
           __m256 ymm2;
           __m256 ymm3;
           __m256 ymm4;
           __m256 ymm5;
           

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1


            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);
            ymm4 = _mm256_load_ps(p_src+32ull);
            _mm256_store_ps(p_dst+32ull,ymm4);
            ymm5 = _mm256_load_ps(p_src+40ull);
            _mm256_store_ps(p_dst+40ull,ymm5);
#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);
            ymm4 = _mm256_load_ps(p_src+32ull);
            ymm5 = _mm256_load_ps(p_src+40ull);
            
            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
            _mm256_store_ps(p_dst+32ull,ymm4);
            _mm256_store_ps(p_dst+40ull,ymm5);
           
#endif 
            sz    -= 48ull;
            p_src += 48ull;
            p_dst += 48ull;
     }

     while(sz >= 32ull)
     {
           __m256 ymm0;
           __m256 ymm1;
           __m256 ymm2;
           __m256 ymm3;
         
#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);

#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);
                       
            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
                       
#endif 
            sz    -= 32ull;
            p_src += 32ull;
            p_dst += 32ull;
     }

     while(sz >= 16ull)
     {
           __m256 ymm0;
           __m256 ymm1;
                   
#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);

#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
                                   
            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
                                   
#endif 
            sz    -= 16ull;
            p_src += 16ull;
            p_dst += 16ull;
     }

     while(sz >= 8ull)
     {
          __m256 ymm0;
          
          ymm0 = _mm256_load_ps(p_src+0ull);
          _mm256_store_ps(p_dst+0ull,ymm0);

          sz    -= 8ull;
          p_src += 8ull;
          p_dst += 8ull;
     }
     
     __asm__ __volatile__ ( "vzeroupper" : : : );

     while(sz)
     {
         const float t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }

}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
void 
gms::common::
avx_memcpy_unroll16x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
{
       float * __restrict__ p_dst{dst};
       float * __restrict__ p_src{src};

       while(((uintptr_t)p_dst & 31) && sz)
       {
             const float t{*p_src};
            *p_dst = t;
             p_src++;
             p_dst++;
             sz--;
       }

       while(sz >= 128ull)
       {

             __m256 ymm0;
             __m256 ymm1;
             __m256 ymm2;
             __m256 ymm3;
             __m256 ymm4;
             __m256 ymm5;
             __m256 ymm6;
             __m256 ymm7;
             __m256 ymm8;
             __m256 ymm9;
             __m256 ymm10;
             __m256 ymm11;
             __m256 ymm12;
             __m256 ymm13;
             __m256 ymm14;
             __m256 ymm15;

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
#endif
            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
#endif
            ymm4 = _mm256_load_ps(p_src+32ull);
            _mm256_store_ps(p_dst+32ull,ymm4);
            ymm5 = _mm256_load_ps(p_src+40ull);
            _mm256_store_ps(p_dst+40ull,ymm5);
            ymm6 = _mm256_load_ps(p_src+48ull);
            _mm256_store_ps(p_dst+48ull,ymm6);
            ymm7 = _mm256_load_ps(p_src+56ull);
            _mm256_store_ps(p_dst+56ull,ymm7);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x160,_MM_HINT_T0);
#endif
            ymm8 = _mm256_load_ps(p_src+64ull);
            _mm256_store_ps(p_dst+64ull,ymm8);
            ymm9 = _mm256_load_ps(p_src+72ull);
            _mm256_store_ps(p_dst+72ull,ymm9);
            ymm10= _mm256_load_ps(p_src+80ull);
            _mm256_store_ps(p_dst+80ull,ymm10);
            ymm11= _mm256_load_ps(p_src+88ull);
            _mm256_store_ps(p_dst+88ull,ymm11);   
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x240,_MM_HINT_T0);
#endif
            ymm12 = _mm256_load_ps(p_src+96ull);
            _mm256_store_ps(p_dst+96ull,ymm12);
            ymm13 = _mm256_load_ps(p_src+104ull);
            _mm256_store_ps(p_dst+104ull,ymm13);
            ymm14 = _mm256_load_ps(p_src+112ull);
            _mm256_store_ps(p_dst+112ull,ymm14);
            ymm15 = _mm256_load_ps(p_src+120ull);
            _mm256_store_ps(p_dst+120ull,ymm15);                    
#else 
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
#endif
            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
#endif
            ymm4 = _mm256_load_ps(p_src+32ull);
            ymm5 = _mm256_load_ps(p_src+40ull);
            ymm6 = _mm256_load_ps(p_src+48ull);
            ymm7 = _mm256_load_ps(p_src+56ull);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x160,_MM_HINT_T0);
#endif
            ymm8 = _mm256_load_ps(p_src+64ull);
            ymm9 = _mm256_load_ps(p_src+72ull);
            ymm10= _mm256_load_ps(p_src+80ull);
            ymm11= _mm256_load_ps(p_src+88ull);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x240,_MM_HINT_T0);
#endif
            ymm12= _mm256_load_ps(p_src+96ull);
            ymm13= _mm256_load_ps(p_src+104ull);
            ymm14= _mm256_load_ps(p_src+112ull);
            ymm15= _mm256_load_ps(p_src+120ull);

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
            _mm256_store_ps(p_dst+32ull,ymm4);
            _mm256_store_ps(p_dst+40ull,ymm5);
            _mm256_store_ps(p_dst+48ull,ymm6);
            _mm256_store_ps(p_dst+56ull,ymm7);
            _mm256_store_ps(p_dst+64ull,ymm8);
            _mm256_store_ps(p_dst+72ull,ymm9);
            _mm256_store_ps(p_dst+80ull,ymm10);
            _mm256_store_ps(p_dst+88ull,ymm11);
            _mm256_store_ps(p_dst+96ull,ymm12);
            _mm256_store_ps(p_dst+104ull,ymm13);
            _mm256_store_ps(p_dst+112ull,ymm14);
            _mm256_store_ps(p_dst+120ull,ymm15);
#endif 
            sz    -= 128ull;
            p_src += 128ull;
            p_dst += 128ull;
       }

       while(sz >= 112ull)
       {
             __m256 ymm0;
             __m256 ymm1;
             __m256 ymm2;
             __m256 ymm3;
             __m256 ymm4;
             __m256 ymm5;
             __m256 ymm6;
             __m256 ymm7;
             __m256 ymm8;
             __m256 ymm9;
             __m256 ymm10;
             __m256 ymm11;
             __m256 ymm12;
             __m256 ymm13;
             

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);

            ymm4 = _mm256_load_ps(p_src+32ull);
            _mm256_store_ps(p_dst+32ull,ymm4);
            ymm5 = _mm256_load_ps(p_src+40ull);
            _mm256_store_ps(p_dst+40ull,ymm5);
            ymm6 = _mm256_load_ps(p_src+48ull);
            _mm256_store_ps(p_dst+48ull,ymm6);
            ymm7 = _mm256_load_ps(p_src+56ull);
            _mm256_store_ps(p_dst+56ull,ymm7);

            ymm8 = _mm256_load_ps(p_src+64ull);
            _mm256_store_ps(p_dst+64ull,ymm8);
            ymm9 = _mm256_load_ps(p_src+72ull);
            _mm256_store_ps(p_dst+72ull,ymm9);
            ymm10= _mm256_load_ps(p_src+80ull);
            _mm256_store_ps(p_dst+80ull,ymm10);
            ymm11= _mm256_load_ps(p_src+88ull);
            _mm256_store_ps(p_dst+88ull,ymm11);   

            ymm12 = _mm256_load_ps(p_src+96ull);
            _mm256_store_ps(p_dst+96ull,ymm12);
            ymm13 = _mm256_load_ps(p_src+104ull);
            _mm256_store_ps(p_dst+104ull,ymm13);
                         
#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);

            ymm4 = _mm256_load_ps(p_src+32ull);
            ymm5 = _mm256_load_ps(p_src+40ull);
            ymm6 = _mm256_load_ps(p_src+48ull);
            ymm7 = _mm256_load_ps(p_src+56ull);

            ymm8 = _mm256_load_ps(p_src+64ull);
            ymm9 = _mm256_load_ps(p_src+72ull);
            ymm10= _mm256_load_ps(p_src+80ull);
            ymm11= _mm256_load_ps(p_src+88ull);

            ymm12= _mm256_load_ps(p_src+96ull);
            ymm13= _mm256_load_ps(p_src+104ull);
          

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
            _mm256_store_ps(p_dst+32ull,ymm4);
            _mm256_store_ps(p_dst+40ull,ymm5);
            _mm256_store_ps(p_dst+48ull,ymm6);
            _mm256_store_ps(p_dst+56ull,ymm7);
            _mm256_store_ps(p_dst+64ull,ymm8);
            _mm256_store_ps(p_dst+72ull,ymm9);
            _mm256_store_ps(p_dst+80ull,ymm10);
            _mm256_store_ps(p_dst+88ull,ymm11);
            _mm256_store_ps(p_dst+96ull,ymm12);
            _mm256_store_ps(p_dst+104ull,ymm13);
            
#endif 
            sz    -= 112ull;
            p_src += 112ull;
            p_dst += 112ull;
       }

       while(sz >= 80ull)
       {
             __m256 ymm0;
             __m256 ymm1;
             __m256 ymm2;
             __m256 ymm3;
             __m256 ymm4;
             __m256 ymm5;
             __m256 ymm6;
             __m256 ymm7;
             __m256 ymm8;
             __m256 ymm9;
                          

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);

            ymm4 = _mm256_load_ps(p_src+32ull);
            _mm256_store_ps(p_dst+32ull,ymm4);
            ymm5 = _mm256_load_ps(p_src+40ull);
            _mm256_store_ps(p_dst+40ull,ymm5);
            ymm6 = _mm256_load_ps(p_src+48ull);
            _mm256_store_ps(p_dst+48ull,ymm6);
            ymm7 = _mm256_load_ps(p_src+56ull);
            _mm256_store_ps(p_dst+56ull,ymm7);

            ymm8 = _mm256_load_ps(p_src+64ull);
            _mm256_store_ps(p_dst+64ull,ymm8);
            ymm9 = _mm256_load_ps(p_src+72ull);
            _mm256_store_ps(p_dst+72ull,ymm9);
                                    
#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);

            ymm4 = _mm256_load_ps(p_src+32ull);
            ymm5 = _mm256_load_ps(p_src+40ull);
            ymm6 = _mm256_load_ps(p_src+48ull);
            ymm7 = _mm256_load_ps(p_src+56ull);

            ymm8 = _mm256_load_ps(p_src+64ull);
            ymm9 = _mm256_load_ps(p_src+72ull);
                 

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
            _mm256_store_ps(p_dst+32ull,ymm4);
            _mm256_store_ps(p_dst+40ull,ymm5);
            _mm256_store_ps(p_dst+48ull,ymm6);
            _mm256_store_ps(p_dst+56ull,ymm7);
            _mm256_store_ps(p_dst+64ull,ymm8);
            _mm256_store_ps(p_dst+72ull,ymm9);
                        
#endif 
            sz    -= 80ull;
            p_src += 80ull;
            p_dst += 80ull;
       }

       while(sz >= 64ull)
       {
             __m256 ymm0;
             __m256 ymm1;
             __m256 ymm2;
             __m256 ymm3;
             __m256 ymm4;
             __m256 ymm5;
             __m256 ymm6;
             __m256 ymm7;
                                      

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);

            ymm4 = _mm256_load_ps(p_src+32ull);
            _mm256_store_ps(p_dst+32ull,ymm4);
            ymm5 = _mm256_load_ps(p_src+40ull);
            _mm256_store_ps(p_dst+40ull,ymm5);
            ymm6 = _mm256_load_ps(p_src+48ull);
            _mm256_store_ps(p_dst+48ull,ymm6);
            ymm7 = _mm256_load_ps(p_src+56ull);
            _mm256_store_ps(p_dst+56ull,ymm7);

                                               
#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);

            ymm4 = _mm256_load_ps(p_src+32ull);
            ymm5 = _mm256_load_ps(p_src+40ull);
            ymm6 = _mm256_load_ps(p_src+48ull);
            ymm7 = _mm256_load_ps(p_src+56ull);

                            

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
            _mm256_store_ps(p_dst+32ull,ymm4);
            _mm256_store_ps(p_dst+40ull,ymm5);
            _mm256_store_ps(p_dst+48ull,ymm6);
            _mm256_store_ps(p_dst+56ull,ymm7);
                                    
#endif 
            sz    -= 64ull;
            p_src += 64ull;
            p_dst += 64ull;
       }

       while(sz >= 32ull)
       {
             __m256 ymm0;
             __m256 ymm1;
             __m256 ymm2;
             __m256 ymm3;
                                                   

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
            ymm2 = _mm256_load_ps(p_src+16ull);
            _mm256_store_ps(p_dst+16ull,ymm2);
            ymm3 = _mm256_load_ps(p_src+24ull);
            _mm256_store_ps(p_dst+24ull,ymm3);

                                                          
#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
            ymm2 = _mm256_load_ps(p_src+16ull);
            ymm3 = _mm256_load_ps(p_src+24ull);

                                     

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
            _mm256_store_ps(p_dst+16ull,ymm2);
            _mm256_store_ps(p_dst+24ull,ymm3);
                                               
#endif 
            sz    -= 32ull;
            p_src += 32ull;
            p_dst += 32ull;
       }

       while(sz >= 16ull)
       {
             __m256 ymm0;
             __m256 ymm1;
                                                               

#if (AVX_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            ymm0 = _mm256_load_ps(p_src+0ull);
            _mm256_store_ps(p_dst+0ull,ymm0);
            ymm1 = _mm256_load_ps(p_src+8ull);
            _mm256_store_ps(p_dst+8ull,ymm1);
           
                                                          
#else 

            ymm0 = _mm256_load_ps(p_src+0ull);
            ymm1 = _mm256_load_ps(p_src+8ull);
                                               

            _mm256_store_ps(p_dst+0ull, ymm0);
            _mm256_store_ps(p_dst+8ull, ymm1);
                                                          
#endif 
            sz    -= 16ull;
            p_src += 16ull;
            p_dst += 16ull;
       }

      while(sz >= 8ull)
      {
          __m256 ymm0;
          
          ymm0 = _mm256_load_ps(p_src+0ull);
          _mm256_store_ps(p_dst+0ull,ymm0);

          sz    -= 8ull;
          p_src += 8ull;
          p_dst += 8ull;
     }
     
     __asm__ __volatile__ ( "vzeroupper" : : : );

     while(sz)
     {
         const float t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }
}