
#ifndef _ATMOSPHERE_LIB_DEFS_H_
#define _ATMOSPHERE_LIB_DEFS_H_



/*
    ---- Check for the Compiler support for C++11 features ----
*/
#if defined (__INTEL_COMPILER) && (__INTEL_COMPILER) < 1400
#error INTEL Compiler version 14.0 needed for partial support of C++11 and INTEL Compiler 15.0 needed for full support of C++11.
#elif defined (__INTEL_COMPILER) && (__INTEL_COMPILER) < 1000
#error Intel Compiler version 10.0 needed to support at least SSE4.
#elif defined (_MSC_VER) && _MSC_VER < 1500
#error MICROSOFT Visual Studio 2013 Compiler or later is required for MathLib compilation.
#endif

/*
    ---- Enable Support of  SIMD ISA ----
*/
#if defined (__INTEL_COMPILER)  && defined (__AVX__) || defined (__AVX2__)
#include <immintrin.h>
#elif defined (__SSE4_2__)
#include <nmmintrin.h>
#elif defined (__SSE4_1__)
#include <smmintrin.h>
#elif defined (__SSSE3__)
#include <tmmintrin.h>
#elif defined (__SSE3__)
#include <pmmintrin.h>
#elif defined (__SSE2__)
#include <emmintrin.h>
#elif defined (__SSE__)
#include <xmmintrin.h>
#elif defined (_MSC_VER)
#include <intrin.h>
#endif

#ifndef USE_OPENMP
#define USE_OPENMP 0x1
#include <omp.h>
#endif
#include <math.h>
#include <iostream>
#include <vector>
#include <memory>
#include <limits>


#endif /*_ATMOSPHERE_LIB_DEFS_H_*/