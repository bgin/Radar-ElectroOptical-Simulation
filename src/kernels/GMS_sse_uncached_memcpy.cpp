

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
#include "GMS_sse_uncached_memcpy.h"


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
void
gms::common::
sse_uncached_memcpy_unroll8x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
{
    float * __restrict__ p_dst{dst};
    float * __restrict__ p_src{src};

    while(((uintptr_t)&p_dst & 15) && sz)
    {
         float t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
    }

    while(sz >= 32ull)
    {
         __m128 t0,
                t1,
                t2,
                t3,
                t4,
                t5,
                t6,
                t7;

         _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+4ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+12ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+20ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+28ull,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
         t4  = _mm_load_ps(p_src+16ull);
         t5  = _mm_load_ps(p_src+20ull);
         t6  = _mm_load_ps(p_src+24ull);
         t7  = _mm_load_ps(p_src+28ull);

         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
         _mm_stream_ps(p_dst+16ull,t4);
         _mm_stream_ps(p_dst+20ull,t5);
         _mm_stream_ps(p_dst+24ull,t6);
         _mm_stream_ps(p_dst+28ull,t7);

         sz    -= 32ull;
         p_src += 32ull;
         p_dst += 32ull;
    }

    _mm_sfence();

    while(sz >= 24ull)
    {
         __m128 t0,
                t1,
                t2,
                t3,
                t4,
                t5;
               
         
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
         t4  = _mm_load_ps(p_src+16ull);
         t5  = _mm_load_ps(p_src+20ull);
         
         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
         _mm_stream_ps(p_dst+16ull,t4);
         _mm_stream_ps(p_dst+20ull,t5);
         
         sz    -= 24ull;
         p_src += 24ull;
         p_dst += 24ull;
    }

    _mm_sfence();

    while(sz >= 16ull)
    {
         __m128 t0,
                t1,
                t2,
                t3;
     
        
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
                  
         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
  
         sz    -= 16ull;
         p_src += 16ull;
         p_dst += 16ull;
    }

    _mm_sfence();

    while(sz >= 8ull)
    {
         __m128 t0,
                t1;
                
         
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
                          
         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
          
         sz    -= 8ull;
         p_src += 8ull;
         p_dst += 8ull;
    }

    _mm_sfence();

    while(sz >= 4ull)
    {
         __m128 t0;
               
         t0  = _mm_load_ps(p_src+0ull);
                        
         _mm_stream_ps(p_dst+0ull, t0);
       
         sz    -= 4ull;
         p_src += 4ull;
         p_dst += 4ull;
    }

    _mm_sfence();

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
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
void 
gms::common::
sse_uncached_memcpy_unroll16x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
{
     float * __restrict__ p_dst{dst};
     float * __restrict__ p_src{src};

     while(((uintptr_t)&p_dst & 15) && sz)
     {
         float t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }

     while(sz >= 64)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;
          __m128 t4;
          __m128 t5;
          __m128 t6;
          __m128 t7;

          __m128 t8;
          __m128 t9;
          __m128 t10;
          __m128 t11;
          __m128 t12;
          __m128 t13;
          __m128 t14;
          __m128 t15;
#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);
          _mm_prefetch((char const*)p_src+4ull,_MM_HINT_T0);
          t1 = _mm_load_ps(p_src+4ull);
          _mm_stream_ps(p_dst+4ull,t1);
          _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
          t2 = _mm_load_ps(p_src+8ull);
          _mm_stream_ps(p_dst+8ull,t2);
          _mm_prefetch((char const*)p_src+12ull,_MM_HINT_T0);
          t3 = _mm_load_ps(p_src+12ull);
          _mm_stream_ps(p_dst+12ull,t3);
          _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
          t4 = _mm_load_ps(p_src+16ull);
          _mm_stream_ps(p_dst+16ull,t4);
          _mm_prefetch((char const*)p_src+20ull,_MM_HINT_T0);
          t5 = _mm_load_ps(p_src+20ull);
          _mm_stream_ps(p_dst+20ull,t5);
          _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
          t6 = _mm_load_ps(p_src+24ull);
          _mm_stream_ps(p_dst+24ull,t6);
          _mm_prefetch((char const*)p_src+28ull,_MM_HINT_T0);
          t7 = _mm_load_ps(p_src+28ull);
          _mm_stream_ps(p_dst+28ull,t7);
          _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
          t8 = _mm_load_ps(p_src+32ull);
          _mm_stream_ps(p_dst+32ull,t8);
          _mm_prefetch((char const*)p_src+36ull,_MM_HINT_T0);
          t9 = _mm_load_ps(p_src+36ull);
          _mm_stream_ps(p_dst+36ull,t9);
          _mm_prefetch((char const*)p_src+40ull,_MM_HINT_T0);
          t10= _mm_load_ps(p_src+40ull);
          _mm_stream_ps(p_dst+40ull,t10);
          _mm_prefetch((char const*)p_src+44ull,_MM_HINT_T0);
          t11= _mm_load_ps(p_src+44ull);
          _mm_stream_ps(p_dst+44ull,t11);
          _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
          t12= _mm_load_ps(p_src+48ull);
          _mm_stream_ps(p_dst+48ull,t12);
          _mm_prefetch((char const*)p_src+52ull,_MM_HINT_T0);
          t13= _mm_load_ps(p_src+52ull);
          _mm_stream_ps(p_dst+52ull,t13);
          _mm_prefetch((char const*)p_src+56ull,_MM_HINT_T0);
          t14= _mm_load_ps(p_src+56ull);
          _mm_stream_ps(p_dst+56ull,t14);
          _mm_prefetch((char const*)p_src+60ull,_MM_HINT_T0);
          t15= _mm_load_ps(p_src+60ull);
          _mm_stream_ps(p_dst+60ull,t15);
#else 
         _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+4ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+12ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+20ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+28ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+32ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+36ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+40ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+44ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+48ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+52ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+56ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+60ull,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
         t4  = _mm_load_ps(p_src+16ull);
         t5  = _mm_load_ps(p_src+20ull);
         t6  = _mm_load_ps(p_src+24ull);
         t7  = _mm_load_ps(p_src+28ull);
         t8  = _mm_load_ps(p_src+32ull);
         t9  = _mm_load_ps(p_src+36ull);
         t10 = _mm_load_ps(p_src+40ull);
         t11 = _mm_load_ps(p_src+44ull);
         t12 = _mm_load_ps(p_src+48ull);
         t13 = _mm_load_ps(p_src+52ull);
         t14 = _mm_load_ps(p_src+56ull);
         t15 = _mm_load_ps(p_src+60ull);
         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
         _mm_stream_ps(p_dst+16ull,t4);
         _mm_stream_ps(p_dst+20ull,t5);
         _mm_stream_ps(p_dst+24ull,t6);
         _mm_stream_ps(p_dst+28ull,t7);
         _mm_stream_ps(p_dst+32ull,t8);
         _mm_stream_ps(p_dst+36ull,t9);
         _mm_stream_ps(p_dst+40ull,t10);
         _mm_stream_ps(p_dst+44ull,t11);
         _mm_stream_ps(p_dst+48ull,t12);
         _mm_stream_ps(p_dst+52ull,t13);
         _mm_stream_ps(p_dst+56ull,t14);
         _mm_stream_ps(p_dst+60ull,t15)
#endif 
         
         sz    -= 64ull;
         p_src += 64ull;
         p_dst += 64ull;
         
          
     }

     _mm_sfence();

     while(sz >= 48ull)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;
          __m128 t4;
          __m128 t5;
          __m128 t6;
          __m128 t7;

          __m128 t8;
          __m128 t9;
          __m128 t10;
          __m128 t11;
          __m128 t12;
#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);
          t1 = _mm_load_ps(p_src+4ull);
          _mm_stream_ps(p_dst+4ull,t1);
          t2 = _mm_load_ps(p_src+8ull);
          _mm_stream_ps(p_dst+8ull,t2);
          t3 = _mm_load_ps(p_src+12ull);
          _mm_stream_ps(p_dst+12ull,t3);
          t4 = _mm_load_ps(p_src+16ull);
          _mm_stream_ps(p_dst+16ull,t4);
          t5 = _mm_load_ps(p_src+20ull);
          _mm_stream_ps(p_dst+20ull,t5);
          t6 = _mm_load_ps(p_src+24ull);
          _mm_stream_ps(p_dst+24ull,t6);
          t7 = _mm_load_ps(p_src+28ull);
          _mm_stream_ps(p_dst+28ull,t7);
          t8 = _mm_load_ps(p_src+32ull);
          _mm_stream_ps(p_dst+32ull,t8);
          t9 = _mm_load_ps(p_src+36ull);
          _mm_stream_ps(p_dst+36ull,t9);
          t10= _mm_load_ps(p_src+40ull);
          _mm_stream_ps(p_dst+40ull,t10);
          t11= _mm_load_ps(p_src+44ull);
          _mm_stream_ps(p_dst+44ull,t11);
#else 
        
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
         t4  = _mm_load_ps(p_src+16ull);
         t5  = _mm_load_ps(p_src+20ull);
         t6  = _mm_load_ps(p_src+24ull);
         t7  = _mm_load_ps(p_src+28ull);
         t8  = _mm_load_ps(p_src+32ull);
         t9  = _mm_load_ps(p_src+36ull);
         t10 = _mm_load_ps(p_src+40ull);
         t11 = _mm_load_ps(p_src+44ull);

         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
         _mm_stream_ps(p_dst+16ull,t4);
         _mm_stream_ps(p_dst+20ull,t5);
         _mm_stream_ps(p_dst+24ull,t6);
         _mm_stream_ps(p_dst+28ull,t7);
         _mm_stream_ps(p_dst+32ull,t8);
         _mm_stream_ps(p_dst+36ull,t9);
         _mm_stream_ps(p_dst+40ull,t10);
         _mm_stream_ps(p_dst+44ull,t11);
#endif 

         sz    -= 48ull;
         p_src += 48ull;
         p_dst += 48ull;
         
     }

     _mm_sfence();

     while(sz >= 32ull)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;
          __m128 t4;
          __m128 t5;
          __m128 t6;
          __m128 t7;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);
          t1 = _mm_load_ps(p_src+4ull);
          _mm_stream_ps(p_dst+4ull,t1);
          t2 = _mm_load_ps(p_src+8ull);
          _mm_stream_ps(p_dst+8ull,t2);
          t3 = _mm_load_ps(p_src+12ull);
          _mm_stream_ps(p_dst+12ull,t3);
          t4 = _mm_load_ps(p_src+16ull);
          _mm_stream_ps(p_dst+16ull,t4);
          t5 = _mm_load_ps(p_src+20ull);
          _mm_stream_ps(p_dst+20ull,t5);
          t6 = _mm_load_ps(p_src+24ull);
          _mm_stream_ps(p_dst+24ull,t6);
          t7 = _mm_load_ps(p_src+28ull);
          _mm_stream_ps(p_dst+28ull,t7);
#else 
         
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
         t4  = _mm_load_ps(p_src+16ull);
         t5  = _mm_load_ps(p_src+20ull);
         t6  = _mm_load_ps(p_src+24ull);
         t7  = _mm_load_ps(p_src+28ull);
         
         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
         _mm_stream_ps(p_dst+16ull,t4);
         _mm_stream_ps(p_dst+20ull,t5);
         _mm_stream_ps(p_dst+24ull,t6);
         _mm_stream_ps(p_dst+28ull,t7);
#endif   

         sz    -= 32ull;
         p_src += 32ull; 
         p_dst += 32ull;       
     }

     _mm_sfence();

     while(sz >= 24ull)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;
          __m128 t4;
          __m128 t5;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);
          t1 = _mm_load_ps(p_src+4ull);
          _mm_stream_ps(p_dst+4ull,t1);
          t2 = _mm_load_ps(p_src+8ull);
          _mm_stream_ps(p_dst+8ull,t2);
          t3 = _mm_load_ps(p_src+12ull);
          _mm_stream_ps(p_dst+12ull,t3);
          t4 = _mm_load_ps(p_src+16ull);
          _mm_stream_ps(p_dst+16ull,t4);
          t5 = _mm_load_ps(p_src+20ull);
          _mm_stream_ps(p_dst+20ull,t5);
#else 
         
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
         t4  = _mm_load_ps(p_src+16ull);
         t5  = _mm_load_ps(p_src+20ull);

         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
         _mm_stream_ps(p_dst+16ull,t4);
         _mm_stream_ps(p_dst+20ull,t5);
#endif   
          
          sz    -= 24ull;
          p_src += 24ull;
          p_dst += 24ull;
     }

     _mm_sfence();

     while(sz >= 16ull)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);
          t1 = _mm_load_ps(p_src+4ull);
          _mm_stream_ps(p_dst+4ull,t1);
          t2 = _mm_load_ps(p_src+8ull);
          _mm_stream_ps(p_dst+8ull,t2);
          t3 = _mm_load_ps(p_src+12ull);
          _mm_stream_ps(p_dst+12ull,t3);
#else 
         
         t0  = _mm_load_ps(p_src+0ull);
         t1  = _mm_load_ps(p_src+4ull);
         t2  = _mm_load_ps(p_src+8ull);
         t3  = _mm_load_ps(p_src+12ull);
        
         _mm_stream_ps(p_dst+0ull, t0);
         _mm_stream_ps(p_dst+4ull, t1);
         _mm_stream_ps(p_dst+8ull, t2);
         _mm_stream_ps(p_dst+12ull,t3);
#endif        

         sz    -= 16ull;
         p_src += 16ull;
         p_dst += 16ull;
     }

     _mm_sfence();

     while(sz >= 8ull)
     {
          __m128 t0;
          __m128 t1;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);
          t1 = _mm_load_ps(p_src+4ull);
          _mm_stream_ps(p_dst+4ull,t1);
#else 
          t0  = _mm_load_ps(p_src+0ull);
          t1  = _mm_load_ps(p_src+4ull);
         
          _mm_stream_ps(p_dst+0ull, t0);
          _mm_stream_ps(p_dst+4ull, t1);
#endif    
         
         sz    -= 8ull;
         p_src += 8ull;
         p_dst += 8ull;
     }

     _mm_sfence();

     while(sz >= 4ull)
     {
           __m128 t0;
           t0 = _mm_load_ps(p_src+0ull);
          _mm_stream_ps(p_dst+0ull,t0);

          sz    -= 4ull;
          p_src += 4ull;
          p_dst += 4ull;
     }

     _mm_sfence();

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
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
void
gms::common::
sse_uncached_memcpy_unroll8x_pd(double * __restrict__ dst,double * __restrict__ src,std::size_t sz)
{
    double * __restrict__ p_dst{dst};
    double * __restrict__ p_src{src};

    while(((uintptr_t)&p_dst & 15) && sz)
    {
         double t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
    }

    while(sz >= 16ull)
    {
         __m128d t0,
                 t1,
                 t2,
                 t3,
                 t4,
                 t5,
                 t6,
                 t7;

         _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+2ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+4ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+6ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+10ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+12ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+14ull,_MM_HINT_T0);
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         t4  = _mm_load_pd(p_src+8ull);
         t5  = _mm_load_pd(p_src+10ull);
         t6  = _mm_load_pd(p_src+12ull);
         t7  = _mm_load_pd(p_src+14ull);

         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         _mm_stream_pd(p_dst+8ull,t4);
         _mm_stream_pd(p_dst+10ull,t5);
         _mm_stream_pd(p_dst+12ull,t6);
         _mm_stream_pd(p_dst+14ull,t7);

         sz    -= 16ull;
         p_src += 16ull;
         p_dst += 16ull;
    }

    _mm_sfence();

    while(sz >= 12ull)
    {
         __m128d t0,
                 t1,
                 t2,
                 t3,
                 t4,
                 t5;
               
       
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         t4  = _mm_load_pd(p_src+8ull);
         t5  = _mm_load_pd(p_src+10ull);
         
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         _mm_stream_pd(p_dst+8ull,t4);
         _mm_stream_pd(p_dst+10ull,t5);
         
         sz    -= 12ull;
         p_src += 12ull;
         p_dst += 12ull;
    }

    _mm_sfence();

    while(sz >= 8ull)
    {
         __m128d t0,
                t1,
                t2,
                t3;
     
         
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
                  
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
  
         sz    -= 8ull;
         p_src += 8ull;
         p_dst += 8ull;
    }

    _mm_sfence();

    while(sz >= 4ull)
    {
         __m128d t0,
                 t1;
                
         
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
                          
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
          
         sz    -= 4ull;
         p_src += 4ull;
         p_dst += 4ull;
    }

    _mm_sfence();

    while(sz >= 2ull)
    {
         __m128d t0;
               
         t0  = _mm_load_pd(p_src+0ull);
                        
         _mm_stream_pd(p_dst+0ull, t0);
       
         sz    -= 2ull;
         p_src += 2ull;
         p_dst += 2ull;
    }

    _mm_sfence();

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
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
void 
gms::common::
sse_uncached_memcpy_unroll16x_pd(double * __restrict__ dst,double * __restrict__ src,std::size_t sz)
{
     double * __restrict__ p_dst{dst};
     double * __restrict__ p_src{src};

     while(((uintptr_t)&p_dst & 15) && sz)
     {
         double t_src{*p_src};
         *p_dst = t_src;
          p_src++;
          p_dst++;
          sz--;
     }

     while(sz >= 32ull)
     {
          __m128d t0;
          __m128d t1;
          __m128d t2;
          __m128d t3;
          __m128d t4;
          __m128d t5;
          __m128d t6;
          __m128d t7;

          __m128d t8;
          __m128d t9;
          __m128d t10;
          __m128d t11;
          __m128d t12;
          __m128d t13;
          __m128d t14;
          __m128d t15;
#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
          t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);
          _mm_prefetch((char const*)p_src+2ull,_MM_HINT_T0);
          t1 = _mm_load_pd(p_src+2ull);
          _mm_stream_pd(p_dst+2ull,t1);
          _mm_prefetch((char const*)p_src+4ull,_MM_HINT_T0);
          t2 = _mm_load_pd(p_src+4ull);
          _mm_stream_pd(p_dst+4ull,t2);
          _mm_prefetch((char const*)p_src+6ull,_MM_HINT_T0);
          t3 = _mm_load_pd(p_src+6ull);
          _mm_stream_pd(p_dst+6ull,t3);
          _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
          t4 = _mm_load_pd(p_src+8ull);
          _mm_stream_pd(p_dst+8ull,t4);
          _mm_prefetch((char const*)p_src+10ull,_MM_HINT_T0);
          t5 = _mm_load_pd(p_src+10ull);
          _mm_stream_pd(p_dst+10ull,t5);
          _mm_prefetch((char const*)p_src+12ull,_MM_HINT_T0);
          t6 = _mm_load_pd(p_src+12ull);
          _mm_stream_pd(p_dst+12ull,t6);
          _mm_prefetch((char const*)p_src+14ull,_MM_HINT_T0);
          t7 = _mm_load_pd(p_src+14ull);
          _mm_stream_pd(p_dst+14ull,t7);
          _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
          t8 = _mm_load_pd(p_src+16ull);
          _mm_stream_pd(p_dst+16ull,t8);
          _mm_prefetch((char const*)p_src+18ull,_MM_HINT_T0);
          t9 = _mm_load_pd(p_src+18ull);
          _mm_stream_pd(p_dst+18ull,t9);
          _mm_prefetch((char const*)p_src+20ull,_MM_HINT_T0);
          t10= _mm_load_pd(p_src+20ull);
          _mm_stream_pd(p_dst+20ull,t10);
          _mm_prefetch((char const*)p_src+22ull,_MM_HINT_T0);
          t11= _mm_load_pd(p_src+22ull);
          _mm_stream_pd(p_dst+22ull,t11);
          _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
          t12= _mm_load_pd(p_src+24ull);
          _mm_stream_pd(p_dst+24ull,t12);
          _mm_prefetch((char const*)p_src+36ull,_MM_HINT_T0);
          t13= _mm_load_pd(p_src+26ull);
          _mm_stream_pd(p_dst+26ull,t13);
          _mm_prefetch((char const*)p_src+28ull,_MM_HINT_T0);
          t14= _mm_load_pd(p_src+28ull);
          _mm_stream_pd(p_dst+28ull,t14);
          _mm_prefetch((char const*)p_src+30ull,_MM_HINT_T0);
          t15= _mm_load_pd(p_src+30ull);
          _mm_stream_pd(p_dst+30ull,t15);
#else 
         _mm_prefetch((char const*)p_src+0ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+2ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+4ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+6ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+8ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+10ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+12ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+14ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+16ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+18ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+20ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+22ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+24ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+28ull,_MM_HINT_T0);
         _mm_prefetch((char const*)p_src+30ull,_MM_HINT_T0);
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         t4  = _mm_load_pd(p_src+8ull);
         t5  = _mm_load_pd(p_src+10ull);
         t6  = _mm_load_pd(p_src+12ull);
         t7  = _mm_load_pd(p_src+14ull);
         t8  = _mm_load_pd(p_src+16ull);
         t9  = _mm_load_pd(p_src+18ull);
         t10 = _mm_load_pd(p_src+20ull);
         t11 = _mm_load_pd(p_src+22ull);
         t12 = _mm_load_pd(p_src+24ull);
         t13 = _mm_load_pd(p_src+26ull);
         t14 = _mm_load_pd(p_src+28ull);
         t15 = _mm_load_pd(p_src+30ull);
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         _mm_stream_pd(p_dst+8ull,t4);
         _mm_stream_pd(p_dst+10ull,t5);
         _mm_stream_pd(p_dst+12ull,t6);
         _mm_stream_pd(p_dst+14ull,t7);
         _mm_stream_pd(p_dst+16ull,t8);
         _mm_stream_pd(p_dst+18ull,t9);
         _mm_stream_pd(p_dst+20ull,t10);
         _mm_stream_pd(p_dst+22ull,t11);
         _mm_stream_pd(p_dst+24ull,t12);
         _mm_stream_pd(p_dst+26ull,t13);
         _mm_stream_pd(p_dst+28ull,t14);
         _mm_stream_pd(p_dst+30ull,t15)
#endif 
         
         sz    -= 32ull;
         p_src += 32ull;
         p_dst += 32ull;
         
          
     }

     _mm_sfence();

     while(sz >= 24ull)
     {
          __m128d t0;
          __m128d t1;
          __m128d t2;
          __m128d t3;
          __m128d t4;
          __m128d t5;
          __m128d t6;
          __m128d t7;

          __m128d t8;
          __m128d t9;
          __m128d t10;
          __m128d t11;
          __m128d t12;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
      
          t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);
          t1 = _mm_load_pd(p_src+2ull);
          _mm_stream_pd(p_dst+2ull,t1);
          t2 = _mm_load_pd(p_src+4ull);
          _mm_stream_pd(p_dst+4ull,t2);
          t3 = _mm_load_pd(p_src+6ull);
          _mm_stream_pd(p_dst+6ull,t3);
          t4 = _mm_load_pd(p_src+8ull);
          _mm_stream_pd(p_dst+8ull,t4);
          t5 = _mm_load_pd(p_src+10ull);
          _mm_stream_pd(p_dst+10ull,t5);
          t6 = _mm_load_pd(p_src+12ull);
          _mm_stream_pd(p_dst+12ull,t6);
          t7 = _mm_load_pd(p_src+14ull);
          _mm_stream_pd(p_dst+14ull,t7);
          t8 = _mm_load_pd(p_src+16ull);
          _mm_stream_pd(p_dst+16ull,t8);
          t9 = _mm_load_pd(p_src+18ull);
          _mm_stream_pd(p_dst+18ull,t9);
          t10= _mm_load_pd(p_src+20ull);
          _mm_stream_pd(p_dst+20ull,t10);
          t11= _mm_load_pd(p_src+22ull);
          _mm_stream_pd(p_dst+22ull,t11);
#else 
         
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         t4  = _mm_load_pd(p_src+8ull);
         t5  = _mm_load_pd(p_src+10ull);
         t6  = _mm_load_pd(p_src+12ull);
         t7  = _mm_load_pd(p_src+14ull);
         t8  = _mm_load_pd(p_src+16ull);
         t9  = _mm_load_pd(p_src+18ull);
         t10 = _mm_load_pd(p_src+20ull);
         t11 = _mm_load_pd(p_src+22ull);
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         _mm_stream_pd(p_dst+8ull,t4);
         _mm_stream_pd(p_dst+10ull,t5);
         _mm_stream_pd(p_dst+12ull,t6);
         _mm_stream_pd(p_dst+14ull,t7);
         _mm_stream_pd(p_dst+16ull,t8);
         _mm_stream_pd(p_dst+18ull,t9);
         _mm_stream_pd(p_dst+20ull,t10);
         _mm_stream_pd(p_dst+22ull,t11);
#endif 

         sz    -= 24ull;
         p_src += 24ull;
         p_dst += 24ull;
         
     }

     _mm_sfence();

     while(sz >= 16ull)
     {
          __m128d t0;
          __m128d t1;
          __m128d t2;
          __m128d t3;
          __m128d t4;
          __m128d t5;
          __m128d t6;
          __m128d t7;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);
          t1 = _mm_load_pd(p_src+2ull);
          _mm_stream_pd(p_dst+2ull,t1);
          t2 = _mm_load_pd(p_src+4ull);
          _mm_stream_pd(p_dst+4ull,t2);
          t3 = _mm_load_pd(p_src+6ull);
          _mm_stream_pd(p_dst+6ull,t3);
          t4 = _mm_load_pd(p_src+8ull);
          _mm_stream_pd(p_dst+8ull,t4);
          t5 = _mm_load_pd(p_src+10ull);
          _mm_stream_pd(p_dst+10ull,t5);
          t6 = _mm_load_pd(p_src+12ull);
          _mm_stream_pd(p_dst+12ull,t6);
          t7 = _mm_load_pd(p_src+14ull);
          _mm_stream_pd(p_dst+14ull,t7);
#else 
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         t4  = _mm_load_pd(p_src+8ull);
         t5  = _mm_load_pd(p_src+10ull);
         t6  = _mm_load_pd(p_src+12ull);
         t7  = _mm_load_pd(p_src+14ull);
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         _mm_stream_pd(p_dst+8ull,t4);
         _mm_stream_pd(p_dst+10ull,t5);
         _mm_stream_pd(p_dst+12ull,t6);
         _mm_stream_pd(p_dst+14ull,t7);
#endif 


         sz    -= 16ull;
         p_src += 16ull; 
         p_dst += 16ull;       
     }

     _mm_sfence();

     while(sz >= 12ull)
     {
          __m128d t0;
          __m128d t1;
          __m128d t2;
          __m128d t3;
          __m128d t4;
          __m128d t5;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);
          t1 = _mm_load_pd(p_src+2ull);
          _mm_stream_pd(p_dst+2ull,t1);
          t2 = _mm_load_pd(p_src+4ull);
          _mm_stream_pd(p_dst+4ull,t2);
          t3 = _mm_load_pd(p_src+6ull);
          _mm_stream_pd(p_dst+6ull,t3);
          t4 = _mm_load_pd(p_src+8ull);
          _mm_stream_pd(p_dst+8ull,t4);
          t5 = _mm_load_pd(p_src+10ull);
          _mm_stream_pd(p_dst+10ull,t5);
#else 
        
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         t4  = _mm_load_pd(p_src+8ull);
         t5  = _mm_load_pd(p_src+10ull);
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         _mm_stream_pd(p_dst+8ull,t4);
         _mm_stream_pd(p_dst+10ull,t5);
#endif 
          
          sz    -= 12ull;
          p_src += 12ull;
          p_dst += 12ull;
     }

     _mm_sfence();

     while(sz >= 8ull)
     {
          __m128d t0;
          __m128d t1;
          __m128d t2;
          __m128d t3;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          
          t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);
          t1 = _mm_load_pd(p_src+2ull);
          _mm_stream_pd(p_dst+2ull,t1);
          t2 = _mm_load_pd(p_src+4ull);
          _mm_stream_pd(p_dst+4ull,t2);
          t3 = _mm_load_pd(p_src+6ull);
          _mm_stream_pd(p_dst+6ull,t3);
         
#else 
         
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
         t2  = _mm_load_pd(p_src+4ull);
         t3  = _mm_load_pd(p_src+6ull);
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         _mm_stream_pd(p_dst+4ull, t2);
         _mm_stream_pd(p_dst+6ull,t3);
         
#endif  

         sz    -= 8ull;
         p_src += 8ull;
         p_dst += 8ull;
     }

     _mm_sfence();

     while(sz >= 4ull)
     {
          __m128d t0;
          __m128d t1;

#if (SSE_UNCACHED_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
         
          t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);
          t1 = _mm_load_pd(p_src+2ull);
          _mm_stream_pd(p_dst+2ull,t1);
         
         
#else 
       
         t0  = _mm_load_pd(p_src+0ull);
         t1  = _mm_load_pd(p_src+2ull);
        
         _mm_stream_pd(p_dst+0ull, t0);
         _mm_stream_pd(p_dst+2ull, t1);
         
         
#endif  
         
         sz    -= 4ull;
         p_src += 4ull;
         p_dst += 4ull;
     }

     _mm_sfence();

     while(sz >= 2ull)
     {
           __m128d t0;
           t0 = _mm_load_pd(p_src+0ull);
          _mm_stream_pd(p_dst+0ull,t0);

          sz    -= 2ull;
          p_src += 2ull;
          p_dst += 2ull;
     }

     _mm_sfence();

     while(sz)
     {
         const double t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }
}


