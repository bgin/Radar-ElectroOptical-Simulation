

#include <immintrin.h>
#include "GMS_sse_memcpy.h"


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
void
gms::common::
sse_memcpy_unroll8x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
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

    while(sz >= 32ULL)
    {
         __m128 t0,
                t1,
                t2,
                t3,
                t4,
                t5,
                t6,
                t7;

         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
         t4  = _mm_load_ps(p_src+16ULL);
         t5  = _mm_load_ps(p_src+20ULL);
         t6  = _mm_load_ps(p_src+24ULL);
         t7  = _mm_load_ps(p_src+28ULL);

         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
         _mm_store_ps(p_dst+16ULL,t4);
         _mm_store_ps(p_dst+20ULL,t5);
         _mm_store_ps(p_dst+24ULL,t6);
         _mm_store_ps(p_dst+28ULL,t7);

         sz    -= 32ULL;
         p_src += 32ULL;
         p_dst += 32ULL;
    }

    while(sz >= 24ULL)
    {
         __m128 t0,
                t1,
                t2,
                t3,
                t4,
                t5;
               
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
         t4  = _mm_load_ps(p_src+16ULL);
         t5  = _mm_load_ps(p_src+20ULL);
         
         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
         _mm_store_ps(p_dst+16ULL,t4);
         _mm_store_ps(p_dst+20ULL,t5);
         
         sz    -= 24ULL;
         p_src += 24ULL;
         p_dst += 24ULL;
    }

    while(sz >= 16ULL)
    {
         __m128 t0,
                t1,
                t2,
                t3;
     
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
                  
         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
  
         sz    -= 16ULL;
         p_src += 16ULL;
         p_dst += 16ULL;
    }

    while(sz >= 8ULL)
    {
         __m128 t0,
                t1;
                
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
                          
         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
          
         sz    -= 8ULL;
         p_src += 8ULL;
         p_dst += 8ULL;
    }

    while(sz >= 4ULL)
    {
         __m128 t0;
               
         t0  = _mm_load_ps(p_src+0ULL);
                        
         _mm_store_ps(p_dst+0ULL, t0);
       
         sz    -= 4ULL;
         p_src += 4ULL;
         p_dst += 4ULL;
    }

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
sse_memcpy_unroll16x_ps(float * __restrict__ dst,float * __restrict__ src,std::size_t sz)
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
#if (SSE_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);
          t1 = _mm_load_ps(p_src+4ULL);
          _mm_store_ps(p_dst+4ULL,t1);
          t2 = _mm_load_ps(p_src+8ULL);
          _mm_store_ps(p_dst+8ULL,t2);
          t3 = _mm_load_ps(p_src+12ULL);
          _mm_store_ps(p_dst+12ULL,t3);
          t4 = _mm_load_ps(p_src+16ULL);
          _mm_store_ps(p_dst+16ULL,t4);
          t5 = _mm_load_ps(p_src+20ULL);
          _mm_store_ps(p_dst+20ULL,t5);
          t6 = _mm_load_ps(p_src+24ULL);
          _mm_store_ps(p_dst+24ULL,t6);
          t7 = _mm_load_ps(p_src+28ULL);
          _mm_store_ps(p_dst+28ULL,t7);
          _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
          t8 = _mm_load_ps(p_src+32ULL);
          _mm_store_ps(p_dst+32ULL,t8);
          t9 = _mm_load_ps(p_src+36ULL);
          _mm_store_ps(p_dst+36ULL,t9);
          t10= _mm_load_ps(p_src+40ULL);
          _mm_store_ps(p_dst+40ULL,t10);
          t11= _mm_load_ps(p_src+44ULL);
          _mm_store_ps(p_dst+44ULL,t11);
          t12= _mm_load_ps(p_src+48ULL);
          _mm_store_ps(p_dst+48ULL,t12);
          t13= _mm_load_ps(p_src+52ULL);
          _mm_store_ps(p_dst+52ULL,t13);
          t14= _mm_load_ps(p_src+56ULL);
          _mm_store_ps(p_dst+56ULL,t14);
          t15= _mm_load_ps(p_src+60ULL);
          _mm_store_ps(p_dst+60ULL,t15);
#else 
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
         t4  = _mm_load_ps(p_src+16ULL);
         t5  = _mm_load_ps(p_src+20ULL);
         t6  = _mm_load_ps(p_src+24ULL);
         t7  = _mm_load_ps(p_src+28ULL);
         _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
         t8  = _mm_load_ps(p_src+32ULL);
         t9  = _mm_load_ps(p_src+36ULL);
         t10 = _mm_load_ps(p_src+40ULL);
         t11 = _mm_load_ps(p_src+44ULL);
         t12 = _mm_load_ps(p_src+48ULL);
         t13 = _mm_load_ps(p_src+52ULL);
         t14 = _mm_load_ps(p_src+56ULL);
         t15 = _mm_load_ps(p_src+60ULL);
         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
         _mm_store_ps(p_dst+16ULL,t4);
         _mm_store_ps(p_dst+20ULL,t5);
         _mm_store_ps(p_dst+24ULL,t6);
         _mm_store_ps(p_dst+28ULL,t7);
         _mm_store_ps(p_dst+32ULL,t8);
         _mm_store_ps(p_dst+36ULL,t9);
         _mm_store_ps(p_dst+40ULL,t10);
         _mm_store_ps(p_dst+44ULL,t11);
         _mm_store_ps(p_dst+48ULL,t12);
         _mm_store_ps(p_dst+52ULL,t13);
         _mm_store_ps(p_dst+56ULL,t14);
         _mm_store_ps(p_dst+60ULL,t15)
#endif 
         
         sz    -= 64ULL;
         p_src += 64ULL;
         p_dst += 64ULL;
         
          
     }

     while(sz >= 48ULL)
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
#if (SSE_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);
          t1 = _mm_load_ps(p_src+4ULL);
          _mm_store_ps(p_dst+4ULL,t1);
          t2 = _mm_load_ps(p_src+8ULL);
          _mm_store_ps(p_dst+8ULL,t2);
          t3 = _mm_load_ps(p_src+12ULL);
          _mm_store_ps(p_dst+12ULL,t3);
          t4 = _mm_load_ps(p_src+16ULL);
          _mm_store_ps(p_dst+16ULL,t4);
          t5 = _mm_load_ps(p_src+20ULL);
          _mm_store_ps(p_dst+20ULL,t5);
          t6 = _mm_load_ps(p_src+24ULL);
          _mm_store_ps(p_dst+24ULL,t6);
          t7 = _mm_load_ps(p_src+28ULL);
          _mm_store_ps(p_dst+28ULL,t7);
          t8 = _mm_load_ps(p_src+32ULL);
          _mm_store_ps(p_dst+32ULL,t8);
          t9 = _mm_load_ps(p_src+36ULL);
          _mm_store_ps(p_dst+36ULL,t9);
          t10= _mm_load_ps(p_src+40ULL);
          _mm_store_ps(p_dst+40ULL,t10);
          t11= _mm_load_ps(p_src+44ULL);
          _mm_store_ps(p_dst+44ULL,t11);
#else 
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
         t4  = _mm_load_ps(p_src+16ULL);
         t5  = _mm_load_ps(p_src+20ULL);
         t6  = _mm_load_ps(p_src+24ULL);
         t7  = _mm_load_ps(p_src+28ULL);
         _mm_prefetch((char const*)p_src+0x80,_MM_HINT_T0);
         t8  = _mm_load_ps(p_src+32ULL);
         t9  = _mm_load_ps(p_src+36ULL);
         t10 = _mm_load_ps(p_src+40ULL);
         t11 = _mm_load_ps(p_src+44ULL);

         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
         _mm_store_ps(p_dst+16ULL,t4);
         _mm_store_ps(p_dst+20ULL,t5);
         _mm_store_ps(p_dst+24ULL,t6);
         _mm_store_ps(p_dst+28ULL,t7);
         _mm_store_ps(p_dst+32ULL,t8);
         _mm_store_ps(p_dst+36ULL,t9);
         _mm_store_ps(p_dst+40ULL,t10);
         _mm_store_ps(p_dst+44ULL,t11);
#endif 

         sz    -= 48ULL;
         p_src += 48ULL;
         p_dst += 48ULL;
         
     }

     while(sz >= 32)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;
          __m128 t4;
          __m128 t5;
          __m128 t6;
          __m128 t7;

#if (SSE_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);
          t1 = _mm_load_ps(p_src+4ULL);
          _mm_store_ps(p_dst+4ULL,t1);
          t2 = _mm_load_ps(p_src+8ULL);
          _mm_store_ps(p_dst+8ULL,t2);
          t3 = _mm_load_ps(p_src+12ULL);
          _mm_store_ps(p_dst+12ULL,t3);
          t4 = _mm_load_ps(p_src+16ULL);
          _mm_store_ps(p_dst+16ULL,t4);
          t5 = _mm_load_ps(p_src+20ULL);
          _mm_store_ps(p_dst+20ULL,t5);
          t6 = _mm_load_ps(p_src+24ULL);
          _mm_store_ps(p_dst+24ULL,t6);
          t7 = _mm_load_ps(p_src+28ULL);
          _mm_store_ps(p_dst+28ULL,t7);
#else 
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
         t4  = _mm_load_ps(p_src+16ULL);
         t5  = _mm_load_ps(p_src+20ULL);
         t6  = _mm_load_ps(p_src+24ULL);
         t7  = _mm_load_ps(p_src+28ULL);
         
         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
         _mm_store_ps(p_dst+16ULL,t4);
         _mm_store_ps(p_dst+20ULL,t5);
         _mm_store_ps(p_dst+24ULL,t6);
         _mm_store_ps(p_dst+28ULL,t7);
#endif   

         sz    -= 32ULL;
         p_src += 32ULL; 
         p_dst += 32ULL;       
     }

     while(sz >= 24ULL)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;
          __m128 t4;
          __m128 t5;

#if (SSE_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);
          t1 = _mm_load_ps(p_src+4ULL);
          _mm_store_ps(p_dst+4ULL,t1);
          t2 = _mm_load_ps(p_src+8ULL);
          _mm_store_ps(p_dst+8ULL,t2);
          t3 = _mm_load_ps(p_src+12ULL);
          _mm_store_ps(p_dst+12ULL,t3);
          t4 = _mm_load_ps(p_src+16ULL);
          _mm_store_ps(p_dst+16ULL,t4);
          t5 = _mm_load_ps(p_src+20ULL);
          _mm_store_ps(p_dst+20ULL,t5);
#else 
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
         t4  = _mm_load_ps(p_src+16ULL);
         t5  = _mm_load_ps(p_src+20ULL);

         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
         _mm_store_ps(p_dst+16ULL,t4);
         _mm_store_ps(p_dst+20ULL,t5);
#endif   
          
          sz    -= 24ULL;
          p_src += 24ULL;
          p_dst += 24ULL;
     }

     while(sz >= 16ULL)
     {
          __m128 t0;
          __m128 t1;
          __m128 t2;
          __m128 t3;

#if (SSE_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);
          t1 = _mm_load_ps(p_src+4ULL);
          _mm_store_ps(p_dst+4ULL,t1);
          t2 = _mm_load_ps(p_src+8ULL);
          _mm_store_ps(p_dst+8ULL,t2);
          t3 = _mm_load_ps(p_src+12ULL);
          _mm_store_ps(p_dst+12ULL,t3);
#else 
         _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
         t0  = _mm_load_ps(p_src+0ULL);
         t1  = _mm_load_ps(p_src+4ULL);
         t2  = _mm_load_ps(p_src+8ULL);
         t3  = _mm_load_ps(p_src+12ULL);
        
         _mm_store_ps(p_dst+0ULL, t0);
         _mm_store_ps(p_dst+4ULL, t1);
         _mm_store_ps(p_dst+8ULL, t2);
         _mm_store_ps(p_dst+12ULL,t3);
#endif        

         sz    -= 16ULL;
         p_src += 16ULL;
         p_dst += 16ULL;
     }

     while(sz >= 8ULL)
     {
          __m128 t0;
          __m128 t1;

#if (SSE_MEMCPY_INTERLEAVE_SIMD_OPS) == 1
          _mm_prefetch((char const*)p_src+0x0,_MM_HINT_T0);
          t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);
          t1 = _mm_load_ps(p_src+4ULL);
          _mm_store_ps(p_dst+4ULL,t1);
#else 
          t0  = _mm_load_ps(p_src+0ULL);
          t1  = _mm_load_ps(p_src+4ULL);
         
          _mm_store_ps(p_dst+0ULL, t0);
          _mm_store_ps(p_dst+4ULL, t1);
#endif    
         
         sz    -= 8ULL;
         p_src += 8ULL;
         p_dst += 8ULL;
     }

     while(sz >= 4ULL)
     {
           __m128 t0;
           t0 = _mm_load_ps(p_src+0ULL);
          _mm_store_ps(p_dst+0ULL,t0);

          sz    -= 4ULL;
          p_src += 4ULL;
          p_dst += 4ULL;
     }

     while(sz)
     {
         const float t_src{*p_src};
         *p_dst = t_src;
         p_src++;
         p_dst++;
         sz--;
     }
}