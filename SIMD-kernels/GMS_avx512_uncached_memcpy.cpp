

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
#include "GMS_avx512_uncached_memcpy.h"

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SKYLAKE-AVX512 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512f")
#endif
void 
gms::common::
avx512_uncached_memcpy_unroll8x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
{
     float * __restrict__ p_dst{dst};
     float * __restrict__ p_src{src};

     while(((uintptr_t)p_dst & 63) && sz)
     {
            const float t{*p_src};
            *p_dst = t;
             p_src++;
             p_dst++;
             sz--;
     }

     while(sz >= 128ull)
     {
           __m512 zmm0;
           __m512 zmm1;
           __m512 zmm2;
           __m512 zmm3;
           __m512 zmm4;
           __m512 zmm5;
           __m512 zmm6;
           __m512 zmm7;


#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);

#endif
            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
            
#endif
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+32,_MM_HINT_T0);
            
#endif
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+48,_MM_HINT_T0);
            
#endif
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0x64,_MM_HINT_T0);
#endif
            zmm4 = _mm512_load_ps(p_src+64ull);
            _mm512_stream_ps(p_dst+64ull,zmm4);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+80,_MM_HINT_T0);
            
#endif
            zmm5 = _mm512_load_ps(p_src+80ull);
            _mm512_stream_ps(p_dst+80ull,zmm5);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+96,_MM_HINT_T0);
            
#endif
            zmm6 = _mm512_load_ps(p_src+96ull);
            _mm512_stream_ps(p_dst+96ull,zmm6);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+112,_MM_HINT_T0);
            
#endif
            zmm7 = _mm512_load_ps(p_src+112ull);
            _mm512_stream_ps(p_dst+112ull,zmm7);
            
#else 
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+64ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+80ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+96ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+112ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);
            zmm4 = _mm512_load_ps(p_src+64ull);
            zmm5 = _mm512_load_ps(p_src+80ull);
            zmm6 = _mm512_load_ps(p_src+96ull);
            zmm7 = _mm512_load_ps(p_src+112ull);

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            _mm512_stream_ps(p_dst+112ull,zmm7);
#endif 
            sz    -= 128ull;
            p_src += 128ull;
            p_dst += 128ull;
    
     }
     
     _mm_sfence();

     while(sz >= 96ull)
     {
           __m512 zmm0;
           __m512 zmm1;
           __m512 zmm2;
           __m512 zmm3;
           __m512 zmm4;
           __m512 zmm5;
           

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1


            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            zmm4 = _mm512_load_ps(p_src+64ull);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            zmm5 = _mm512_load_ps(p_src+80ull);
            _mm512_stream_ps(p_dst+80ull,zmm5);
#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);
            zmm4 = _mm512_load_ps(p_src+64ull);
            zmm5 = _mm512_load_ps(p_src+80ull);
            
            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            _mm512_stream_ps(p_dst+80ull,zmm5);
           
#endif 
            sz    -= 96ull;
            p_src += 96ull;
            p_dst += 96ull;
     }

     _mm_sfence();

     while(sz >= 64ull)
     {
           __m512 zmm0;
           __m512 zmm1;
           __m512 zmm2;
           __m512 zmm3;
         
#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);

#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);
                       
            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
                       
#endif 
            sz    -= 64ull;
            p_src += 64ull;
            p_dst += 64ull;
     }

     _mm_sfence();

     while(sz >= 32ull)
     {
           __m512 zmm0;
           __m512 zmm1;
                   
#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);

#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
                                   
            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
                                   
#endif 
            sz    -= 32ull;
            p_src += 32ull;
            p_dst += 32ull;
     }

     _mm_sfence();

     while(sz >= 16ull)
     {
          __m512 zmm0;
          
          zmm0 = _mm512_load_ps(p_src+0ull);
          _mm512_stream_ps(p_dst+0ull,zmm0);

          sz    -= 16ull;
          p_src += 16ull;
          p_dst += 16ull;
     }

     _mm_sfence();
     
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
#pragma intel optimization_parameter target_arch=SKYLAKE-AVX512 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512f")
#endif
void 
gms::common::
avx512_uncached_memcpy_unroll16x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
{
       float * __restrict__ p_dst{dst};
       float * __restrict__ p_src{src};

       while(((uintptr_t)p_dst & 63) && sz)
       {
             const float t{*p_src};
            *p_dst = t;
             p_src++;
             p_dst++;
             sz--;
       }

       while(sz >= 256ull)
       {

             __m512 zmm0;
             __m512 zmm1;
             __m512 zmm2;
             __m512 zmm3;
             __m512 zmm4;
             __m512 zmm5;
             __m512 zmm6;
             __m512 zmm7;
             __m512 ymm8;
             __m512 ymm9;
             __m512 zmm10;
             __m512 zmm11;
             __m512 zmm12;
             __m512 zmm13;
             __m512 zmm14;
             __m512 zmm15;

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
#endif
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
#endif
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
#endif
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+64ull,_MM_HINT_T0);
#endif
            zmm4 = _mm512_load_ps(p_src+64ull);
            _mm512_stream_ps(p_dst+64ull,zmm4);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+80ull,_MM_HINT_T0);
#endif
            zmm5 = _mm512_load_ps(p_src+80ull);
            _mm512_stream_ps(p_dst+80ull,zmm5);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+96ull,_MM_HINT_T0);
#endif
            zmm6 = _mm512_load_ps(p_src+96ull);
            _mm512_stream_ps(p_dst+96ull,zmm6);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+112ull,_MM_HINT_T0);
#endif
            zmm7 = _mm512_load_ps(p_src+112ull);
            _mm512_stream_ps(p_dst+112ull,zmm7);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+128,_MM_HINT_T0);
#endif
            ymm8 = _mm512_load_ps(p_src+128ull);
            _mm512_stream_ps(p_dst+128ull,ymm8);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+144ull,_MM_HINT_T0);
#endif
            ymm9 = _mm512_load_ps(p_src+144ull);
            _mm512_stream_ps(p_dst+144ull,ymm9);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+160ull,_MM_HINT_T0);
#endif
            zmm10= _mm512_load_ps(p_src+160ull);
            _mm512_stream_ps(p_dst+160ull,zmm10);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+176ull,_MM_HINT_T0);
#endif
            zmm11= _mm512_load_ps(p_src+176ull);
            _mm512_stream_ps(p_dst+176ull,zmm11);   
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+192ull,_MM_HINT_T0);
#endif
            zmm12 = _mm512_load_ps(p_src+192ull);
            _mm512_stream_ps(p_dst+192ull,zmm12);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+208ull,_MM_HINT_T0);
#endif
            zmm13 = _mm512_load_ps(p_src+208ull);
            _mm512_stream_ps(p_dst+208ull,zmm13);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+224ull,_MM_HINT_T0);
#endif
            zmm14 = _mm512_load_ps(p_src+224ull);
            _mm512_stream_ps(p_dst+224ull,zmm14);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+240ull,_MM_HINT_T0);
#endif
            zmm15 = _mm512_load_ps(p_src+240ull);
            _mm512_stream_ps(p_dst+240ull,zmm15);                    
#else 
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+64ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+80ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+96ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+112ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+128ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+144ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+160ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+176ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+192ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+208ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+224ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+240ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);
            zmm4 = _mm512_load_ps(p_src+64ull);
            zmm5 = _mm512_load_ps(p_src+80ull);
            zmm6 = _mm512_load_ps(p_src+96ull);
            zmm7 = _mm512_load_ps(p_src+112ull);
            ymm8 = _mm512_load_ps(p_src+128ull);
            ymm9 = _mm512_load_ps(p_src+144ull);
            zmm10= _mm512_load_ps(p_src+160ull);
            zmm11= _mm512_load_ps(p_src+176ull);
            zmm12= _mm512_load_ps(p_src+192ull);
            zmm13= _mm512_load_ps(p_src+208ull);
            zmm14= _mm512_load_ps(p_src+224ull);
            zmm15= _mm512_load_ps(p_src+240ull);

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            _mm512_stream_ps(p_dst+112ull,zmm7);
            _mm512_stream_ps(p_dst+128ull,ymm8);
            _mm512_stream_ps(p_dst+144ull,ymm9);
            _mm512_stream_ps(p_dst+160ull,zmm10);
            _mm512_stream_ps(p_dst+176ull,zmm11);
            _mm512_stream_ps(p_dst+192ull,zmm12);
            _mm512_stream_ps(p_dst+208ull,zmm13);
            _mm512_stream_ps(p_dst+224ull,zmm14);
            _mm512_stream_ps(p_dst+240ull,zmm15);
#endif 
            sz    -= 256ull;
            p_src += 256ull;
            p_dst += 256ull;
       }

       _mm_sfence();

       while(sz >= 224ull)
       {
             __m512 zmm0;
             __m512 zmm1;
             __m512 zmm2;
             __m512 zmm3;
             __m512 zmm4;
             __m512 zmm5;
             __m512 zmm6;
             __m512 zmm7;
             __m512 ymm8;
             __m512 ymm9;
             __m512 zmm10;
             __m512 zmm11;
             __m512 zmm12;
             __m512 zmm13;
             

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);

            zmm4 = _mm512_load_ps(p_src+64ull);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            zmm5 = _mm512_load_ps(p_src+80ull);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            zmm6 = _mm512_load_ps(p_src+96ull);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            zmm7 = _mm512_load_ps(p_src+112ull);
            _mm512_stream_ps(p_dst+112ull,zmm7);

            ymm8 = _mm512_load_ps(p_src+128ull);
            _mm512_stream_ps(p_dst+128ull,ymm8);
            ymm9 = _mm512_load_ps(p_src+144ull);
            _mm512_stream_ps(p_dst+144ull,ymm9);
            zmm10= _mm512_load_ps(p_src+160ull);
            _mm512_stream_ps(p_dst+160ull,zmm10);
            zmm11= _mm512_load_ps(p_src+176ull);
            _mm512_stream_ps(p_dst+176ull,zmm11);   

            zmm12 = _mm512_load_ps(p_src+192ull);
            _mm512_stream_ps(p_dst+192ull,zmm12);
            zmm13 = _mm512_load_ps(p_src+208ull);
            _mm512_stream_ps(p_dst+208ull,zmm13);
                         
#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);

            zmm4 = _mm512_load_ps(p_src+64ull);
            zmm5 = _mm512_load_ps(p_src+80ull);
            zmm6 = _mm512_load_ps(p_src+96ull);
            zmm7 = _mm512_load_ps(p_src+112ull);

            ymm8 = _mm512_load_ps(p_src+128ull);
            ymm9 = _mm512_load_ps(p_src+144ull);
            zmm10= _mm512_load_ps(p_src+160ull);
            zmm11= _mm512_load_ps(p_src+176ull);

            zmm12= _mm512_load_ps(p_src+192ull);
            zmm13= _mm512_load_ps(p_src+208ull);
          

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            _mm512_stream_ps(p_dst+112ull,zmm7);
            _mm512_stream_ps(p_dst+128ull,ymm8);
            _mm512_stream_ps(p_dst+144ull,ymm9);
            _mm512_stream_ps(p_dst+160ull,zmm10);
            _mm512_stream_ps(p_dst+176ull,zmm11);
            _mm512_stream_ps(p_dst+192ull,zmm12);
            _mm512_stream_ps(p_dst+208ull,zmm13);
            
#endif 
            sz    -= 224ull;
            p_src += 224ull;
            p_dst += 224ull;
       }

       _mm_sfence();

       while(sz >= 160ull)
       {
             __m512 zmm0;
             __m512 zmm1;
             __m512 zmm2;
             __m512 zmm3;
             __m512 zmm4;
             __m512 zmm5;
             __m512 zmm6;
             __m512 zmm7;
             __m512 ymm8;
             __m512 ymm9;
                          

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

             zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);

            zmm4 = _mm512_load_ps(p_src+64ull);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            zmm5 = _mm512_load_ps(p_src+80ull);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            zmm6 = _mm512_load_ps(p_src+96ull);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            zmm7 = _mm512_load_ps(p_src+112ull);
            _mm512_stream_ps(p_dst+112ull,zmm7);

            ymm8 = _mm512_load_ps(p_src+128ull);
            _mm512_stream_ps(p_dst+128ull,ymm8);
            ymm9 = _mm512_load_ps(p_src+144ull);
            _mm512_stream_ps(p_dst+144ull,ymm9);
                                    
#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);
            zmm4 = _mm512_load_ps(p_src+64ull);
            zmm5 = _mm512_load_ps(p_src+80ull);
            zmm6 = _mm512_load_ps(p_src+96ull);
            zmm7 = _mm512_load_ps(p_src+112ull);
            ymm8 = _mm512_load_ps(p_src+128ull);
            ymm9 = _mm512_load_ps(p_src+144ull);

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            _mm512_stream_ps(p_dst+112ull,zmm7);
            _mm512_stream_ps(p_dst+128ull,ymm8);
            _mm512_stream_ps(p_dst+144ull,ymm9);
                        
#endif 
            sz    -= 160ull;
            p_src += 160ull;
            p_dst += 160ull;
       }

       _mm_sfence();

       while(sz >= 128ull)
       {
             __m512 zmm0;
             __m512 zmm1;
             __m512 zmm2;
             __m512 zmm3;
             __m512 zmm4;
             __m512 zmm5;
             __m512 zmm6;
             __m512 zmm7;
                                      

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);

            zmm4 = _mm512_load_ps(p_src+64ull);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            zmm5 = _mm512_load_ps(p_src+80ull);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            zmm6 = _mm512_load_ps(p_src+96ull);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            zmm7 = _mm512_load_ps(p_src+112ull);
            _mm512_stream_ps(p_dst+112ull,zmm7);

                                               
#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);
            zmm4 = _mm512_load_ps(p_src+64ull);
            zmm5 = _mm512_load_ps(p_src+80ull);
            zmm6 = _mm512_load_ps(p_src+96ull);
            zmm7 = _mm512_load_ps(p_src+112ull);

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
            _mm512_stream_ps(p_dst+64ull,zmm4);
            _mm512_stream_ps(p_dst+80ull,zmm5);
            _mm512_stream_ps(p_dst+96ull,zmm6);
            _mm512_stream_ps(p_dst+112ull,zmm7);                

                                              
#endif 
            sz    -= 128ull;
            p_src += 128ull;
            p_dst += 128ull;
       }

       _mm_sfence();

       while(sz >= 64ull)
       {
             __m512 zmm0;
             __m512 zmm1;
             __m512 zmm2;
             __m512 zmm3;
                                                   

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
            zmm2 = _mm512_load_ps(p_src+32ull);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            zmm3 = _mm512_load_ps(p_src+48ull);
            _mm512_stream_ps(p_dst+48ull,zmm3);

                                                          
#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
            zmm2 = _mm512_load_ps(p_src+32ull);
            zmm3 = _mm512_load_ps(p_src+48ull);

                                     

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
            _mm512_stream_ps(p_dst+32ull,zmm2);
            _mm512_stream_ps(p_dst+48ull,zmm3);
                                               
#endif 
            sz    -= 64ull;
            p_src += 64ull;
            p_dst += 64ull;
       }

       _mm_sfence();

       while(sz >= 32ull)
       {
             __m512 zmm0;
             __m512 zmm1;
                                                               

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_ps(p_src+0ull);
            _mm512_stream_ps(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_ps(p_src+16ull);
            _mm512_stream_ps(p_dst+16ull,zmm1);
           
                                                          
#else 

            zmm0 = _mm512_load_ps(p_src+0ull);
            zmm1 = _mm512_load_ps(p_src+16ull);
                                               

            _mm512_stream_ps(p_dst+0ull, zmm0);
            _mm512_stream_ps(p_dst+16ull, zmm1);
                                                          
#endif 
            sz    -= 32ull;
            p_src += 32ull;
            p_dst += 32ull;
       }

       _mm_sfence();

      while(sz >= 16ull)
      {
          __m512 zmm0;
          
          zmm0 = _mm512_load_ps(p_src+0ull);
          _mm512_stream_ps(p_dst+0ull,zmm0);

          sz    -= 16ull;
          p_src += 16ull;
          p_dst += 16ull;
     }

     _mm_sfence();
     
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
#pragma intel optimization_parameter target_arch=SKYLAKE-AVX512 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512f")
#endif
void 
gms::common::
avx512_uncached_memcpy_unroll8x_pd(double * __restrict__ dst,double * __restrict__ src,std::size_t sz)
{
     double * __restrict__ p_dst{dst};
     double * __restrict__ p_src{src};

     while(((uintptr_t)p_dst & 63) && sz)
     {
            const double t{*p_src};
            *p_dst = t;
             p_src++;
             p_dst++;
             sz--;
     }

     while(sz >= 64ull)
     {
           __m512d zmm0;
           __m512d zmm1;
           __m512d zmm2;
           __m512d zmm3;
           __m512d zmm4;
           __m512d zmm5;
           __m512d zmm6;
           __m512d zmm7;


#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
#endif
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
#endif
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
#endif
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
#endif
            zmm4 = _mm512_load_pd(p_src+32ull);
            _mm512_stream_pd(p_dst+32ull,zmm4);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+40ull,_MM_HINT_T0);
#endif
            zmm5 = _mm512_load_pd(p_src+40ull);
            _mm512_stream_pd(p_dst+40ull,zmm5);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
#endif
            zmm6 = _mm512_load_pd(p_src+48ull);
            _mm512_stream_pd(p_dst+48ull,zmm6);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+56ull,_MM_HINT_T0);
#endif
            zmm7 = _mm512_load_pd(p_src+56ull);
            _mm512_stream_pd(p_dst+56ull,zmm7);
#else 
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+40ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+56ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);
            zmm4 = _mm512_load_pd(p_src+32ull);
            zmm5 = _mm512_load_pd(p_src+40ull);
            zmm6 = _mm512_load_pd(p_src+48ull);
            zmm7 = _mm512_load_pd(p_src+56ull);

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            _mm512_stream_pd(p_dst+56ull,zmm7);
#endif 
            sz    -= 64ull;
            p_src += 64ull;
            p_dst += 64ull;
    
     }

     _mm_sfence();

     while(sz >= 48ull)
     {
           __m512d zmm0;
           __m512d zmm1;
           __m512d zmm2;
           __m512d zmm3;
           __m512d zmm4;
           __m512d zmm5;
           

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1


            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            zmm4 = _mm512_load_pd(p_src+32ull);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            zmm5 = _mm512_load_pd(p_src+40ull);
            _mm512_stream_pd(p_dst+40ull,zmm5);
#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);
            zmm4 = _mm512_load_pd(p_src+32ull);
            zmm5 = _mm512_load_pd(p_src+40ull);
            
            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            _mm512_stream_pd(p_dst+40ull,zmm5);
           
#endif 
            sz    -= 48ull;
            p_src += 48ull;
            p_dst += 48ull;
     }

     _mm_sfence();

     while(sz >= 32ull)
     {
           __m512d zmm0;
           __m512d zmm1;
           __m512d zmm2;
           __m512d zmm3;
         
#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);

#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);
                       
            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
                       
#endif 
            sz    -= 32ull;
            p_src += 32ull;
            p_dst += 32ull;
     }

     _mm_sfence();

     while(sz >= 16ull)
     {
           __m512d zmm0;
           __m512d zmm1;
                   
#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);

#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
                                   
            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
                                   
#endif 
            sz    -= 16ull;
            p_src += 16ull;
            p_dst += 16ull;
     }

     _mm_sfence();

     while(sz >= 8ull)
     {
          __m512d zmm0;
          
          zmm0 = _mm512_load_pd(p_src+0ull);
          _mm512_stream_pd(p_dst+0ull,zmm0);

          sz    -= 8ull;
          p_src += 8ull;
          p_dst += 8ull;
     }

     _mm_sfence();
     
     __asm__ __volatile__ ( "vzeroupper" : : : );

     while(sz)
     {
         const double t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }

}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SKYLAKE-AVX512 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx512f")
#endif
void 
gms::common::
avx512_uncached_memcpy_unroll16x_pd(double * __restrict__ dst,double * __restrict__ src,std::size_t sz)
{
       double * __restrict__ p_dst{dst};
       double * __restrict__ p_src{src};

       while(((uintptr_t)p_dst & 63) && sz)
       {
             const double t{*p_src};
            *p_dst = t;
             p_src++;
             p_dst++;
             sz--;
       }

       while(sz >= 128ull)
       {

             __m512d zmm0;
             __m512d zmm1;
             __m512d zmm2;
             __m512d zmm3;
             __m512d zmm4;
             __m512d zmm5;
             __m512d zmm6;
             __m512d zmm7;
             __m512d ymm8;
             __m512d ymm9;
             __m512d zmm10;
             __m512d zmm11;
             __m512d zmm12;
             __m512d zmm13;
             __m512d zmm14;
             __m512d zmm15;

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
#endif
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
#endif
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
#endif
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
#endif
            zmm4 = _mm512_load_pd(p_src+32ull);
            _mm512_stream_pd(p_dst+32ull,zmm4);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+40ull,_MM_HINT_T0);
#endif           
            zmm5 = _mm512_load_pd(p_src+40ull);
            _mm512_stream_pd(p_dst+40ull,zmm5);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
#endif
            zmm6 = _mm512_load_pd(p_src+48ull);
            _mm512_stream_pd(p_dst+48ull,zmm6);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+56ull,_MM_HINT_T0);
#endif
            zmm7 = _mm512_load_pd(p_src+56ull);
            _mm512_stream_pd(p_dst+56ull,zmm7);
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+64ull,_MM_HINT_T0);
#endif
            ymm8 = _mm512_load_pd(p_src+64ull);
            _mm512_stream_pd(p_dst+64ull,ymm8);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+72ull,_MM_HINT_T0);
#endif
            ymm9 = _mm512_load_pd(p_src+72ull);
            _mm512_stream_pd(p_dst+72ull,ymm9);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+80ull,_MM_HINT_T0);
#endif
            zmm10= _mm512_load_pd(p_src+80ull);
            _mm512_stream_pd(p_dst+80ull,zmm10);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+88ull,_MM_HINT_T0);
#endif
            zmm11= _mm512_load_pd(p_src+88ull);
            _mm512_stream_pd(p_dst+88ull,zmm11);   
#if (AVX_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+96ull,_MM_HINT_T0);
#endif
            zmm12 = _mm512_load_pd(p_src+96ull);
            _mm512_stream_pd(p_dst+96ull,zmm12);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+104ull,_MM_HINT_T0);
#endif
            zmm13 = _mm512_load_pd(p_src+104ull);
            _mm512_stream_pd(p_dst+104ull,zmm13);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+112ull,_MM_HINT_T0);
#endif
            zmm14 = _mm512_load_pd(p_src+112ull);
            _mm512_stream_pd(p_dst+112ull,zmm14);
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+120ull,_MM_HINT_T0);
#endif
            zmm15 = _mm512_load_pd(p_src+120ull);
            _mm512_stream_pd(p_dst+120ull,zmm15);                    
#else 
#if (AVX512_UNCACHED_MEMCPY_SOFT_PREFETCHING) == 1
            _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+40ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+56ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+64ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+72ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+80ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+88ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+96ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+104ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+112ull,_MM_HINT_T0);
            _mm_prefetch((char const*)p_src+120ull,_MM_HINT_T0);
#endif
            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);
            zmm4 = _mm512_load_pd(p_src+32ull);
            zmm5 = _mm512_load_pd(p_src+40ull);
            zmm6 = _mm512_load_pd(p_src+48ull);
            zmm7 = _mm512_load_pd(p_src+56ull);
            ymm8 = _mm512_load_pd(p_src+64ull);
            ymm9 = _mm512_load_pd(p_src+72ull);
            zmm10= _mm512_load_pd(p_src+80ull);
            zmm11= _mm512_load_pd(p_src+88ull);
            zmm12= _mm512_load_pd(p_src+96ull);
            zmm13= _mm512_load_pd(p_src+104ull);
            zmm14= _mm512_load_pd(p_src+112ull);
            zmm15= _mm512_load_pd(p_src+120ull);

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            _mm512_stream_pd(p_dst+56ull,zmm7);
            _mm512_stream_pd(p_dst+64ull,ymm8);
            _mm512_stream_pd(p_dst+72ull,ymm9);
            _mm512_stream_pd(p_dst+80ull,zmm10);
            _mm512_stream_pd(p_dst+88ull,zmm11);
            _mm512_stream_pd(p_dst+96ull,zmm12);
            _mm512_stream_pd(p_dst+104ull,zmm13);
            _mm512_stream_pd(p_dst+112ull,zmm14);
            _mm512_stream_pd(p_dst+120ull,zmm15);
#endif 
            sz    -= 128ull;
            p_src += 128ull;
            p_dst += 128ull;
       }

       _mm_sfence();

       while(sz >= 112ull)
       {
             __m512d zmm0;
             __m512d zmm1;
             __m512d zmm2;
             __m512d zmm3;
             __m512d zmm4;
             __m512d zmm5;
             __m512d zmm6;
             __m512d zmm7;
             __m512d ymm8;
             __m512d ymm9;
             __m512d zmm10;
             __m512d zmm11;
             __m512d zmm12;
             __m512d zmm13;
             

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);

            zmm4 = _mm512_load_pd(p_src+32ull);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            zmm5 = _mm512_load_pd(p_src+40ull);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            zmm6 = _mm512_load_pd(p_src+48ull);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            zmm7 = _mm512_load_pd(p_src+56ull);
            _mm512_stream_pd(p_dst+56ull,zmm7);

            ymm8 = _mm512_load_pd(p_src+64ull);
            _mm512_stream_pd(p_dst+64ull,ymm8);
            ymm9 = _mm512_load_pd(p_src+72ull);
            _mm512_stream_pd(p_dst+72ull,ymm9);
            zmm10= _mm512_load_pd(p_src+80ull);
            _mm512_stream_pd(p_dst+80ull,zmm10);
            zmm11= _mm512_load_pd(p_src+88ull);
            _mm512_stream_pd(p_dst+88ull,zmm11);   

            zmm12 = _mm512_load_pd(p_src+96ull);
            _mm512_stream_pd(p_dst+96ull,zmm12);
            zmm13 = _mm512_load_pd(p_src+104ull);
            _mm512_stream_pd(p_dst+104ull,zmm13);
                         
#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);

            zmm4 = _mm512_load_pd(p_src+32ull);
            zmm5 = _mm512_load_pd(p_src+40ull);
            zmm6 = _mm512_load_pd(p_src+48ull);
            zmm7 = _mm512_load_pd(p_src+56ull);

            ymm8 = _mm512_load_pd(p_src+64ull);
            ymm9 = _mm512_load_pd(p_src+72ull);
            zmm10= _mm512_load_pd(p_src+80ull);
            zmm11= _mm512_load_pd(p_src+88ull);

            zmm12= _mm512_load_pd(p_src+96ull);
            zmm13= _mm512_load_pd(p_src+104ull);
          

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            _mm512_stream_pd(p_dst+56ull,zmm7);
            _mm512_stream_pd(p_dst+64ull,ymm8);
            _mm512_stream_pd(p_dst+72ull,ymm9);
            _mm512_stream_pd(p_dst+80ull,zmm10);
            _mm512_stream_pd(p_dst+88ull,zmm11);
            _mm512_stream_pd(p_dst+96ull,zmm12);
            _mm512_stream_pd(p_dst+104ull,zmm13);
            
#endif 
            sz    -= 112ull;
            p_src += 112ull;
            p_dst += 112ull;
       }

       _mm_sfence();

       while(sz >= 80ull)
       {
             __m512d zmm0;
             __m512d zmm1;
             __m512d zmm2;
             __m512d zmm3;
             __m512d zmm4;
             __m512d zmm5;
             __m512d zmm6;
             __m512d zmm7;
             __m512d ymm8;
             __m512d ymm9;
                          

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);

            zmm4 = _mm512_load_pd(p_src+32ull);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            zmm5 = _mm512_load_pd(p_src+40ull);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            zmm6 = _mm512_load_pd(p_src+48ull);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            zmm7 = _mm512_load_pd(p_src+56ull);
            _mm512_stream_pd(p_dst+56ull,zmm7);

            ymm8 = _mm512_load_pd(p_src+64ull);
            _mm512_stream_pd(p_dst+64ull,ymm8);
            ymm9 = _mm512_load_pd(p_src+72ull);
            _mm512_stream_pd(p_dst+72ull,ymm9);
                                    
#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);

            zmm4 = _mm512_load_pd(p_src+32ull);
            zmm5 = _mm512_load_pd(p_src+40ull);
            zmm6 = _mm512_load_pd(p_src+48ull);
            zmm7 = _mm512_load_pd(p_src+56ull);

            ymm8 = _mm512_load_pd(p_src+64ull);
            ymm9 = _mm512_load_pd(p_src+72ull);
                 

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            _mm512_stream_pd(p_dst+56ull,zmm7);
            _mm512_stream_pd(p_dst+64ull,ymm8);
            _mm512_stream_pd(p_dst+72ull,ymm9);
                        
#endif 
            sz    -= 80ull;
            p_src += 80ull;
            p_dst += 80ull;
       }

       _mm_sfence();

       while(sz >= 64ull)
       {
             __m512d zmm0;
             __m512d zmm1;
             __m512d zmm2;
             __m512d zmm3;
             __m512d zmm4;
             __m512d zmm5;
             __m512d zmm6;
             __m512d zmm7;
                                      

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);

            zmm4 = _mm512_load_pd(p_src+32ull);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            zmm5 = _mm512_load_pd(p_src+40ull);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            zmm6 = _mm512_load_pd(p_src+48ull);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            zmm7 = _mm512_load_pd(p_src+56ull);
            _mm512_stream_pd(p_dst+56ull,zmm7);

                                               
#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);

            zmm4 = _mm512_load_pd(p_src+32ull);
            zmm5 = _mm512_load_pd(p_src+40ull);
            zmm6 = _mm512_load_pd(p_src+48ull);
            zmm7 = _mm512_load_pd(p_src+56ull);

                            

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
            _mm512_stream_pd(p_dst+32ull,zmm4);
            _mm512_stream_pd(p_dst+40ull,zmm5);
            _mm512_stream_pd(p_dst+48ull,zmm6);
            _mm512_stream_pd(p_dst+56ull,zmm7);
                                    
#endif 
            sz    -= 64ull;
            p_src += 64ull;
            p_dst += 64ull;
       }

       _mm_sfence();

       while(sz >= 32ull)
       {
             __m512d zmm0;
             __m512d zmm1;
             __m512d zmm2;
             __m512d zmm3;
                                                   

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
            zmm2 = _mm512_load_pd(p_src+16ull);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            zmm3 = _mm512_load_pd(p_src+24ull);
            _mm512_stream_pd(p_dst+24ull,zmm3);

                                                          
#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
            zmm2 = _mm512_load_pd(p_src+16ull);
            zmm3 = _mm512_load_pd(p_src+24ull);

                                     

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
            _mm512_stream_pd(p_dst+16ull,zmm2);
            _mm512_stream_pd(p_dst+24ull,zmm3);
                                               
#endif 
            sz    -= 32ull;
            p_src += 32ull;
            p_dst += 32ull;
       }

       _mm_sfence();

       while(sz >= 16ull)
       {
             __m512d zmm0;
             __m512d zmm1;
                                                               

#if (AVX512_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1

            zmm0 = _mm512_load_pd(p_src+0ull);
            _mm512_stream_pd(p_dst+0ull,zmm0);
            zmm1 = _mm512_load_pd(p_src+8ull);
            _mm512_stream_pd(p_dst+8ull,zmm1);
           
                                                          
#else 

            zmm0 = _mm512_load_pd(p_src+0ull);
            zmm1 = _mm512_load_pd(p_src+8ull);
                                               

            _mm512_stream_pd(p_dst+0ull, zmm0);
            _mm512_stream_pd(p_dst+8ull, zmm1);
                                                          
#endif 
            sz    -= 16ull;
            p_src += 16ull;
            p_dst += 16ull;
       }

       _mm_sfence();

      while(sz >= 8ull)
      {
          __m512d zmm0;
          
          zmm0 = _mm512_load_pd(p_src+0ull);
          _mm512_stream_pd(p_dst+0ull,zmm0);

          sz    -= 8ull;
          p_src += 8ull;
          p_dst += 8ull;
     }

     _mm_sfence();
     
     __asm__ __volatile__ ( "vzeroupper" : : : );

     while(sz)
     {
         const double t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }
}


