

#include <immintrin.h>
#include "GMS_sse_memset.h"

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
void 
gms::common::
sse_memset_unroll8x_ps(float * __restrict__ dst,const float filler,std::size_t size)
{
   
    const __m128 vfiller{_mm_setr_ps(filler,filler,filler,filler)};
    float * __restrict__ p_dst{dst};
    const float cfiller{filler};

    while(((uintptr_t)p_dst & 15) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 32ULL)
    {
         
        _mm_storeu_ps(p_dst+0, vfiller);
        _mm_storeu_ps(p_dst+4, vfiller);
        _mm_storeu_ps(p_dst+8, vfiller);
        _mm_storeu_ps(p_dst+12,vfiller);
        _mm_storeu_ps(p_dst+16,vfiller);
        _mm_storeu_ps(p_dst+20,vfiller);
        _mm_storeu_ps(p_dst+24,vfiller);
        _mm_storeu_ps(p_dst+28,vfiller);
         
        size  -= 32ULL;
        p_dst += 32ULL;
    }

    while(size >= 24ULL)
    {
        _mm_storeu_ps(p_dst+0, vfiller);
        _mm_storeu_ps(p_dst+4, vfiller);
        _mm_storeu_ps(p_dst+8, vfiller);
        _mm_storeu_ps(p_dst+12,vfiller);
        _mm_storeu_ps(p_dst+16,vfiller);
        _mm_storeu_ps(p_dst+20,vfiller);

        size  -= 24ULL;
        p_dst += 24ULL;
    }

    while(size >= 16ULL)
    {
        _mm_storeu_ps(p_dst+0, vfiller);
        _mm_storeu_ps(p_dst+4, vfiller);
        _mm_storeu_ps(p_dst+8, vfiller);
        _mm_storeu_ps(p_dst+12,vfiller);

        size  -= 16ULL;
        p_dst += 16ULL;
    }

    while(size >= 8ULL)
    {
        _mm_storeu_ps(p_dst+0, vfiller);
        _mm_storeu_ps(p_dst+4, vfiller);

        size  -= 8ULL;
        p_dst += 8ULL;
    }

    while(size >= 4ULL)
    {
        _mm_storeu_ps(p_dst+0, vfiller);

        size  -= 4ULL;
        p_dst += 4ULL;
    }

    while(size)
    {
        *p_dst = cfiller;
        size--;
        p_dst++;
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
sse_memset_unroll16x_ps(float * __restrict__ dst,const float filler,std::size_t size)
{
    const __m128 vfiller{_mm_setr_ps(filler,filler,filler,filler)};
    float * __restrict__ p_dst = dst;
    std::size_t __i;
    const float cfiller{filler};

   
    for(__i = 0ull; (__i+63ull) < size; __i += 64ull)
    {
          _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+4ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+12ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+16ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+20ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+24ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+28ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+32ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+36ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+40ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+44ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+48ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+52ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+56ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+60ULL], vfiller);
    }

    for(; (__i+47ull) < size;  __i += 48ull)
    {
          _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+4ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+12ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+16ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+20ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+24ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+28ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+32ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+36ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+40ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+44ULL], vfiller);
    }

    for(; (__i+31ull) < size; __i += 32ull)
    {
          _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+4ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+12ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+16ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+20ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+24ULL], vfiller);
          _mm_storeu_ps(&p_dst[__i+28ULL], vfiller);
    }

    for(; (__i+23ull) < size; __i += 24ull)
    {
          _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+4ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+12ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+16ull], vfiller);
          _mm_storeu_ps(&p_dst[__i+20ULL], vfiller);
    }

    for(; (__i+15ull) < size; __i += 16ull)
    {
          _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+4ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+8ull],  vfiller);
          _mm_storeu_ps(&p_dst[__i+12ull], vfiller);
    }

    for(; (__i+7ull) < size; __i += 8ull)
    {
           _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
           _mm_storeu_ps(&p_dst[__i+4ull],  vfiller);
    }

    for(; (__i+3ull) < size; __i += 4ull)
    {
           _mm_storeu_ps(&p_dst[__i+0ull],  vfiller);
    }

    for(; (__i+0ull) < size; __i += 1ull)
    {
          p_dst[__i] = cfiller;
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
sse_memset_unroll8x_pd(double * __restrict__ dst,const double filler,std::size_t size)
{
    const __m128d vfiller{_mm_setr_pd(filler,filler)};
    double * __restrict__ p_dst{dst};
    const double cfiller{filler};

    while(((uintptr_t)p_dst & 15) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 16ULL)
    {
         
        _mm_store_pd(p_dst+0ULL, vfiller);
        _mm_store_pd(p_dst+2ULL, vfiller);
        _mm_store_pd(p_dst+4ULL, vfiller);
        _mm_store_pd(p_dst+6ULL, vfiller);
        _mm_store_pd(p_dst+8ULL, vfiller);
        _mm_store_pd(p_dst+10ULL,vfiller);
        _mm_store_pd(p_dst+12ULL,vfiller);
        _mm_store_pd(p_dst+14ULL,vfiller);
         
        p_dst += 16ULL;
        size  -= 16ULL;
        
    }

    while(size >= 12ULL)
    {
        _mm_store_pd(p_dst+0ULL, vfiller);
        _mm_store_pd(p_dst+2ULL, vfiller);
        _mm_store_pd(p_dst+4ULL, vfiller);
        _mm_store_pd(p_dst+6ULL, vfiller);
        _mm_store_pd(p_dst+8ULL, vfiller);
        _mm_store_pd(p_dst+10ULL,vfiller);
        
        p_dst += 12ULL;
        size  -= 12ULL;
        
    }

    while(size >= 8ULL)
    {
        _mm_store_pd(p_dst+0ULL, vfiller);
        _mm_store_pd(p_dst+2ULL, vfiller);
        _mm_store_pd(p_dst+4ULL, vfiller);
        _mm_store_pd(p_dst+6ULL, vfiller);
        
        p_dst += 8ULL;
        size  -= 8ULL;
       
    }

    while(size >= 4ULL)
    {
        _mm_store_pd(p_dst+0ULL, vfiller);
        _mm_store_pd(p_dst+2ULL, vfiller);
        
        p_dst += 4ULL;
        size  -= 4ULL;
        
    }

    while(size >= 2ULL)
    {
        _mm_store_pd(p_dst+0ULL, vfiller);
        
        p_dst += 2ULL;
        size  -= 2ULL;
        
    }

    while(size)
    {
        *p_dst = cfiller;
         p_dst++;
         size--;
        
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
sse_memset_unroll16x_pd(double * __restrict__ dst,const double filler,std::size_t size)
{
     const __m128d vfiller{_mm_setr_pd(filler,filler)};
     double * __restrict__ p_dst{dst};
     double cfiller{filler};

     while(((uintptr_t)p_dst & 15) && size)
     {
          *p_dst++ = cfiller;
          --size;
     }

     while(size >= 32ULL)
     {
          _mm_store_pd(p_dst+0ULL, vfiller);
          _mm_store_pd(p_dst+2ULL, vfiller);
          _mm_store_pd(p_dst+4ULL, vfiller);
          _mm_store_pd(p_dst+6ULL, vfiller);
          _mm_store_pd(p_dst+8ULL, vfiller);
          _mm_store_pd(p_dst+10ULL,vfiller);
          _mm_store_pd(p_dst+12ULL,vfiller);
          _mm_store_pd(p_dst+14ULL,vfiller);
          _mm_store_pd(p_dst+16ULL,vfiller);
          _mm_store_pd(p_dst+18ULL,vfiller);
          _mm_store_pd(p_dst+20ULL,vfiller);
          _mm_store_pd(p_dst+22ULL,vfiller);
          _mm_store_pd(p_dst+24ULL,vfiller);
          _mm_store_pd(p_dst+26ULL,vfiller);
          _mm_store_pd(p_dst+28ULL,vfiller);
          _mm_store_pd(p_dst+30ULL,vfiller);

          size  -= 32ULL;
          p_dst += 32ULL;
      }

      while(size >= 24ULL)
      {
          _mm_store_pd(p_dst+0ULL, vfiller);
          _mm_store_pd(p_dst+2ULL, vfiller);
          _mm_store_pd(p_dst+4ULL, vfiller);
          _mm_store_pd(p_dst+6ULL, vfiller);
          _mm_store_pd(p_dst+8ULL, vfiller);
          _mm_store_pd(p_dst+10ULL,vfiller);
          _mm_store_pd(p_dst+12ULL,vfiller);
          _mm_store_pd(p_dst+14ULL,vfiller);
          _mm_store_pd(p_dst+16ULL,vfiller);
          _mm_store_pd(p_dst+18ULL,vfiller);
          _mm_store_pd(p_dst+20ULL,vfiller);
          _mm_store_pd(p_dst+22ULL,vfiller);

          size  -= 24ULL;
          p_dst += 24ULL;
      }

      while(size >= 16ULL)
      {
          _mm_store_pd(p_dst+0ULL, vfiller);
          _mm_store_pd(p_dst+2ULL, vfiller);
          _mm_store_pd(p_dst+4ULL, vfiller);
          _mm_store_pd(p_dst+6ULL, vfiller);
          _mm_store_pd(p_dst+8ULL, vfiller);
          _mm_store_pd(p_dst+10ULL,vfiller);
          _mm_store_pd(p_dst+12ULL,vfiller);
          _mm_store_pd(p_dst+14ULL,vfiller);

          size  -= 16ULL;
          p_dst += 16ULL;
      }

      while(size >= 8ULL)
      {
          _mm_store_pd(p_dst+0ULL, vfiller);
          _mm_store_pd(p_dst+2ULL, vfiller);
          _mm_store_pd(p_dst+4ULL, vfiller);
          _mm_store_pd(p_dst+6ULL, vfiller);

          size  -= 8ULL;
          p_dst += 8ULL;
      }

      while(size >= 4ULL)
      {
          _mm_store_pd(p_dst+0ULL, vfiller);
          _mm_store_pd(p_dst+2ULL, vfiller);

          size  -= 4ULL;
          p_dst += 4ULL;
      }

      while(size >= 2ULL)
      {
          _mm_store_pd(p_dst+0ULL, vfiller);

          size  -= 2ULL;
          p_dst += 2ULL;
      }

      while(size)
      {
          *p_dst = cfiller;
          size--;
          p_dst++;
      }
}