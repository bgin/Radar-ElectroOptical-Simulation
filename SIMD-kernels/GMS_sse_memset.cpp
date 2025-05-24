

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
sse_memset_unroll8x_ps(float * dst,const float filler,std::size_t size)
{
   
    const __m128 vfiller{_mm_setr_ps(filler,filler,filler,filler)};
    float * p_dst{dst};
    const float cfiller{filler};

    while(((uintptr_t)p_dst & 15) && size)
    {
        *p_dst++ = cfiller;
         --size;
    }

    while(size >= 32ULL)
    {
         
        _mm_store_ps(p_dst+0, vfiller);
        _mm_store_ps(p_dst+4, vfiller);
        _mm_store_ps(p_dst+8, vfiller);
        _mm_store_ps(p_dst+12,vfiller);
        _mm_store_ps(p_dst+16,vfiller);
        _mm_store_ps(p_dst+20,vfiller);
        _mm_store_ps(p_dst+24,vfiller);
        _mm_store_ps(p_dst+28,vfiller);
         
        size  -= 32ULL;
        p_dst += 32ULL;
    }

    while(size >= 24ULL)
    {
        _mm_store_ps(p_dst+0, vfiller);
        _mm_store_ps(p_dst+4, vfiller);
        _mm_store_ps(p_dst+8, vfiller);
        _mm_store_ps(p_dst+12,vfiller);
        _mm_store_ps(p_dst+16,vfiller);
        _mm_store_ps(p_dst+20,vfiller);

        size  -= 24ULL;
        p_dst += 24ULL;
    }

    while(size >= 16ULL)
    {
        _mm_store_ps(p_dst+0, vfiller);
        _mm_store_ps(p_dst+4, vfiller);
        _mm_store_ps(p_dst+8, vfiller);
        _mm_store_ps(p_dst+12,vfiller);

        size  -= 16ULL;
        p_dst += 16ULL;
    }

    while(size >= 8ULL)
    {
        _mm_store_ps(p_dst+0, vfiller);
        _mm_store_ps(p_dst+4, vfiller);

        size  -= 8ULL;
        p_dst += 8ULL;
    }

    while(size >= 4ULL)
    {
        _mm_store_ps(p_dst+0, vfiller);

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
sse_memset_unroll8x_pd(double * dst,const double filler,std::size_t size)
{
    const __m128d vfiller{_mm_setr_pd(filler,filler)};
    double * p_dst{dst};
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
         
        size  -= 16ULL;
        p_dst += 16ULL;
    }

    while(size >= 12ULL)
    {
        _mm_store_pd(p_dst+0ULL, vfiller);
        _mm_store_pd(p_dst+2ULL, vfiller);
        _mm_store_pd(p_dst+4ULL, vfiller);
        _mm_store_pd(p_dst+6ULL, vfiller);
        _mm_store_pd(p_dst+8ULL, vfiller);
        _mm_store_pd(p_dst+10ULL,vfiller);

        size  -= 12ULL;
        p_dst += 12ULL;
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