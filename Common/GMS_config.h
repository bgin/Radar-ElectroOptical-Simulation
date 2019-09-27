
#ifndef __GMS_CONFIG_H__
#define __GMS_CONFIG_H__

#if !defined(GMS_CONFIG_MAJOR)
#define GMS_CONFIG_MAJOR 1
#endif

#if !defined(GMS_CONFIG_MINOR)
#define GMS_CONFIG_MINOR 0
#endif

#if !defined(GMS_CONFIG_MICRO)
#define GMS_CONFIG_MICRO 0
#endif

#if !defined(GMS_CONFIG_FULLVER)
#define GMS_CONFIG_FULLVER 1000
#endif

#if !defined(GMS_CONFIG_CREATE_DATE)
#define GMS_CONFIG_CREATE_DATE "26-09-2018 21:06 +00200 (THR 26 SEP 2019 GMT+2)"
#endif
/*
Set this value to latest build date/time
*/
#if !defined(GMS_CONFIG_BUILD_DATE)
#define GMS_CONFIG_BUILD_DATE " "
#endif

#if !defined(GMS_CONFIG_AUTHOR)
#define GMS_CONFIG_AUTHOR  "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com"
#endif

#if !defined(GMS_CONFIG_DESCRIPT)
#define GMS_CONFIG_DESCRIPT "GMS configuration global settings."
#endif



/* Start of Compiler specific declarations.* /

/* Compiler supported CPP version
as reported by reading __cplusplus macro def.*/
#if defined (__cplusplus)
#define GMS_CXX_98 199711L
#define GMS_CXX_11 201103L
#define GMS_CXX_14 201420L
#endif

/* Deteremine current version supported by ICC.*/
#if defined (__cplusplus) && !defined (__INTEL_CXX11_MODE__)
#if GMS_CXX_98 < GMS_CXX_11
#define GMS_DEFAULT_CXX_VERSION 199711L
#else
#define GMS_DEFAULT_CXX_VERSION 201103L
#endif
#endif

// Is Intel Compiler choosen as default
// library compiler?

#if defined __INTEL_COMPILER
#define GMS_COMPILED_BY_ICC 1
#else
#define GMS_COMPILED_BY_ICC 0
#define GMS_COMPILED_BY_MSVC 1
#endif

/* Is 64bit mode current? */
#if defined (_Win64)
   #if (defined (_M_AMD64) || defined (_M_X64_) || defined (__amd64) ) \
        	&& !defined (__x86_64__)
    #define __x86_64__ 1
   #endif
#endif

/* Determine architectural support for full set
of GP registers*/
#if __x86_64__ == 1
#define GMS_HAS_FULL_GPR_SET 16
#elif __x86_64__ == 0
#define GMS_HAS_FULL_GPR_SET 8
#else
#error "COMPILE_TIME_ERROR: Cannot determine 32bit or 64bit mode!"
#endif

/* Determine architectural support for full set
of SIMD registers*/
#if __x86_64__ == 1
#define GMS_HAS_FULL_SIMD_REGSET 32
#elif __x86_64__ == 0
#define GMS_HAS_FULL_SIMD_REGSET 16
#else
#error "COMPILE_TIME_ERROR: Cannot determine 32bit or 64bit mode!"
#endif

/*
Compiler optimization settings.
*/
#if defined GMS_COMPILED_BY_ICC
#define GMS_NO_OPTIMIZATION 0
#define GMS_OPTIMIZATION_O1 1
#define GMS_OPTIMIZATION_O2 2
#define GMS_OPTIMIZATION_O3 3
#endif

/*
   CPU user-friemdly name:
   Current build machine CPU type.
*/
#if !defined (MACHINE_CPU_NAME)
#define MACHINE_CPU_NAME "Intel Core i7 4770 HQ"
#endif

/*
Using OpenMP.
*/
#if !defined(USE_OPENMP)
#define USE_OPENMP 1
#endif

#if USE_OPENMP == 1 && defined (GMS_COMPILED_BY_ICC)
#include <omp.h>
#elif defined (GMS_COMPILED_BY_MSVC)
#include <omp.h>
#else
#error "COMPILE_TIME_ERROR: Unsupported Compiler version!"
#endif

#if USE_OPENMP == 1 && __INTEL_COMPILER >= 1500
#define OMP_VER 40
#elif __INTEL_COMPILER < 1500
#define OMP_VER 10
#else
#error "COMPILE_TIME_ERROR: Unsupported Compiler version!"
#endif

/*
Intel MKL support.
Include all headers - master header file.
*/
#if !defined(USE_MKL)
#define USE_MKL 1
#endif

#if USE_MKL == 1 && defined (GMS_COMPILED_BY_ICC)
#include <mkl.h>
#else
#error "COMPILE_TIME_ERROR: Cannot include MKL: headers!"
#endif

// Use accurate floating-point library Intel lbfp754 library
#if (GMS_COMPILED_BY_ICC) == 1
#define USE_ACCURATE_IEEE754_2008_FP 1
#endif

// Debug mode.
#if defined (_DEBUG) && !defined (__linux)
#define GMS_DEBUG_ON 1
#include <crtdbg.h>
#elif defined (_DEBUG) && defined (__linux)
#define GMS_DEBUG_ON 1
#endif

#if defined (__linux) 
    #if (GMS_DEBUG_ON) == 1
        #if !defined (CALL_MALLOPT_M_CHECK_ACTION)
          #define CALL_MALLOPT_M_CHECK_ACTION 1
        #endif
    #endif
#endif

// Cache line size
#if !defined CACHE_LINE_SIZE
#define CACHE_LINE_SIZE   64
#endif

#if !defined (PAD_TO_CACHE_LINE)
#define PAD_TO_CACHE_LINE(var_size) \
	char pad[(CACHE_LINE_SIZE)-sizeof((var_size))];
#endif

#if !defined (PAD_TO)
#define PAD_TO(ordinal,size) \
  char pad##ordinal[(size)];
#endif

#if !defined (PAD_TO_ALIGNED) && !defined (__linux)
#define PAD_TO_ALIGNED(alignment,ordinal,size) \
	__declspec(align((alignment))) char pad##ordinal[(size)];
#elif !defined (PAD_TO_ALIGNED) && defined (__linux)
#define PAD_TO_ALIGNED(alignment,ordinal,size) #
        __attribute__(align((alignment))) char pad##ordinal[(size)];
#endif

#if !defined (ALIGN_AT) && !defined (__linux)
#define ALIGN_AT(alignment) \
	__declspec(align((alignment)))
#elif !defined (ALIGN_AT) && defined (__linux)
#define ALIGN_AT(alignment)  \
        __attribute__((align(alignment)))
#endif

#if !defined (USE_STRUCT_PADDING)
#define USE_STRUCT_PADDING 1
#endif

// Verbose debug
// Print every errorneous situation
#if GMS_DEBUG_ON == 1
#define GMS_VERBOSE_DEBUG_ON 1
#else
#define GMS_VERBOSE_DEBUG_ON 0
#endif

// Caching  stores
// default value is 1
// Set to '0' in order to bypass
// caching memory stores.
#ifndef GMS_CACHING_STORES
#define GMS_CACHE_MEM_STORES 1
#endif

#ifndef GMS_MAN_PREFETCH
#define GMS_MAN_PREFETCH 0
#endif

// Verbose mode on
#ifndef GMS_VERBOSE_ON
#define GMS_VERBOSE_ON 1
#endif

/*
Compiler software prefetching.
*/
#if GMS_COMPILED_BY_ICC == 1
#define ICC_PREFETCH_L1 1
#define ICC_PREFETCH_L2 2
#define ICC_PREFETCH_L3 3
#define ICC_PREFETCH_NTA 4
#endif

#if defined (ICC_PREFETCH_L1)
constexpr unsigned int L1_PREF_SHORT{ 4 };
constexpr unsigned int L2_PREF_LONG{ 8 };
#endif

// Highly CPU dependend
#if defined MACHINE_CPU_NAME  // Later should be encoded as ASCII integral type
#if defined (ICC_PREFETCH_L1)
constexpr unsigned int L1_MAX_FLOATS{ 8000 };
constexpr size_t L1_MAX_DOUBLES{ 4000 };
constexpr unsigned long long L1_MAX_CAPACITY_F32{ 8000ULL };
constexpr unsigned long long L1_MAX_CAPACITY_F64{4000ULL};
#endif
// Highly CPU dependend
#if defined ICC_PREFETCH_L2
constexpr unsigned int L2_MAX_FLOATS{ 8 * L1_MAX_FLOATS };
constexpr unsigned int L2_MAX_DOUBLES{ 8 * L1_MAX_DOUBLES };
#endif
// Highly CPU dependend
#if defined ICC_PREFETCH_L3
constexpr unsigned int L3_MAX_FLOATS{ 1572864 };
constexpr unsigned int L3_MAX_DOUBLES{786432};

#else
#error "COMPILE_TIME_ERROR: Please define your CPU type!!"
#endif
#endif

/*
	Performance measurement
*/
#if !defined(PRECISE_PERF_TIMING)
#define PRECISE_PERF_TIMING 1
#endif

#if PRECISE_PERF_TIMING == 0
#define CRUDE_PERF_TIMING 1
#endif

#if !defined (USE_PERF_PROFILER)
    #define USE_PERF_PROFILER 1
#endif

#if !defined (USE_MSR_TOOLS)
    #define USE_MSR_TOOLS 0
#endif

#if !defined (GMS_MANUAL_UNROLL)
#define GMS_MANUAL_UNROLL  1
#endif

#if GMS_MANUAL_UNROLL == 1
#define UNROLL_2X 2
#define UNROLL_4X 4
#define UNROLL_8X 8
#if __x86_64__ == 1
#define UNROLL_16X 16
#define UNROLL_32X 32
#endif
#endif

#if GMS_COMPILED_BY_ICC == 1
#ifndef _MM_MALLOC_ALIGNED_
#define _MM_MALLOC_ALIGNED_ 1
#endif
#endif

#if (_MM_MALLOC_ALIGNED_) == 1
constexpr unsigned long long align16B{ 16 };
constexpr unsigned long long align32B{ 32 };
constexpr unsigned long long align64B{ 64 };
#endif

#ifndef ADD_PADDING_64B_LOOP_PEEL
#define ADD_PADDING_64B_LOOP_PEEL 1
#endif

#if (USE_PERF_PROFILER) == 1
       #if !defined (PERF_PROFILE_FUNCTIONS)
           #define PERF_PROFILE_FUNCTIONS 1
       #endif
#endif

// Padding to add at the end of dynamically allocated array
// to enable high performing sequence for loop peeling
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
constexpr int padding64B{64};
#endif


#ifndef USE_AUTO_VEC
#define USE_AUTO_VEC 1
#endif

#if !defined (USE_TRIP_COUNT_HINT)
    #define USE_TRIP_COUNT_HINT 1
#endif

#if !defined (USE_VECTOR_REDUCE)
    #define USE_VECTOR_REDUCE 1
#endif



//
// Round down to 4 elements
//
#ifndef ROUND_TO_FOUR
#define ROUND_TO_FOUR(x,y) ((x) & ~((y)-1))
#endif
// Basically the same as above.
#ifndef ROUND_TO_EIGHT
#define ROUND_TO_EIGHT(x,y) ((x) & ~((y)-1))
#endif

#ifndef ROUND_TO_SIXTEEN
#define ROUND_TO_SIXTEEN(x,y) ((x) & ~((y)-1))
#endif



/* Fixing problem of defined in WinDef macros:
1) min.
2) max.
*/

#if defined _MSC_VER
#if defined (_WINDEF_) && defined (min) && defined (max)
#undef min
#undef max
#endif
#if !defined NOMINMAX
#define NOMINMAX
#endif
#endif

#if !defined (C_WRAPPER_ODEPACK_FPTR_WORKAROUND)
#define C_WRAPPER_ODEPACK_FPTR_WORKAROUND 1
#endif

// Set this to 1 in order to enable NT-stores.
#if !defined (USE_NT_STORES)
#define USE_NT_STORES 0
#endif





#define PRINT_MESSAGE(msg) std::cout << (msg) << "\n";

#if !defined (PRINT_CALLSTACK_ON_ERROR)
#define PRINT_CALLSTACK_ON_ERROR 1
#endif

	

#define PRINT_MESSAGE_VALUE(msg,val) std::cout << (msg) << (val) << "\n";

#define PRINT_MESSAGE_VALUE_2ARGS(msg,arg1,arg2) \
	std::cout << (msg) << std::dec <<  \
			"value 1: " << (arg1) << \
			"value 2: " << (arg2) << "\n";

#define PRINT_MESSAGE_VALUE_3ARGS(msg,arg1,arg2,arg3) \
	std::cout << (msg) << std::dec << \
			"value 1: " << (arg1)  << \
			"value 2: " << (arg2)  << \
			"value 3: " << (arg3)  << "\n";

#define PRINT_MESSAGE_VALUE_4ARGS(msg,arg1,arg2,arg3,arg4) \
	std::cout << (msg) << std::dec  << \
			"value 1: " << (arg1)	<< \
			"value 2: " << (arg2)   << \
			"value 3: " << (arg3)   << \
			"value 4: " << (arg4)   << "\n";

#define PRINT_MESSAGE_VALUE_5ARGS(msg,arg1,arg2,arg3,arg4,arg5) \
	std::cout << (msg) << std::dec << \
	"value 1: " << (arg1) << \
	"value 2: " << (arg2) << \
	"value 3: " << (arg3) << \
	"value 4: " << (arg4) << \
	"value 5: " << (arg5) << "\n";

#define PRINT_MESSAGE_HEX_1ARG(msg,arg1)  \
	std::cout << (msg) << std::hex << \
	"value 1: " << "0x" << (arg1) << "\n";

#define PRINT_MESSAGE_HEX_2ARGS(msg,arg1,arg2) \
	std::cout << (msg) << std::hex << \
	"value 1: " << "0x" << (arg1)  << \
	"value 2: " << "0x" << (arg2) << "\n";

#define PRINT_MESSAGE_HEX_3ARGS(msg,arg1,arg2,arg3) \
	std::cout << (msg) << std::hex << \
	"value 1: " << "0x" << (arg1) <<  \
	"value 2: " << "0x" << (arg2) <<  \
	"value 3: " << "0x" << (arg3) << "\n";

#define PRINT_MESSAGE_HEX_4ARGS(msg,arg1,arg2,arg3,arg4) \
	std::cout << (msg) << std::hex << \
	"value 1: " << "0x" << (arg1) << \
	"value 2: " << "0x" << (arg2) << \
	"value 3: " << "0x" << (arg3) << \
	"value 4: " << "0x" << (arg4) << "\n";

#define PRINT_MESSAGE_HEX_5ARGS(msg,arg1,arg2,arg3,arg4,arg5) \
	std::cout << (msg) << std::hex << \
	"value 1: " << "0x" << (arg1) << \
	"value 2: " << "0x" << (arg2) << \
	"value 3: " << "0x" << (arg3) << \
	"value 4: " << "0x" << (arg4) << \
	"value 5: " << "0x" << (arg5) << "\n";

#define PRINT_MESSAGE_FLOAT_1ARG(msg,precision,arg1) \
	std::cout << (msg) << std::setprecision((precision)) << \
	          << std::showpoint <<  \
			  "value 1: " << (arg1) << "\n";

#define PRINT_MESSAGE_FLOAT_2ARG(msg,precision,arg1,arg2) \
	std::cout << (msg) << std::setprecision((precision)) << \
			  << std::showpoint << \
			  "value 1: " << (arg1) << \
			  "value 2: " << (arg2) << "\n";

#define PRINT_MESSAGE_FLOAT_3ARG(msg,precision,arg1,arg2,arg3) \
	std::cout << (msg) << std::setprecision((precision)) << \
			  << std::showpoint << \
			  "value 1: " << (arg1) << \
			  "value 2: " << (arg2) << \
			  "value 3: " << (arg3) << "\n";

#define PRINT_MESSAGE_FLOAT_4ARG(msg,precision,arg1,arg2,arg3,arg4) \
	std::cout << (msg) << std::setprecision((precision)) << \
			  << std::showpoint << \
			  "value 1: " << (arg1) << \
			  "value 2: " << (arg2) << \
			  "value 3: " << (arg3) << \
			  "value 4: " << (arg4) << "\n";

#define PRINT_MESSAGE_FLOAT_5ARG(msg,arg1,arg2,arg3,arg4,arg5) \
	std::cout << (msg) << std::setprecision((precision)) << \
			  << std::showpoint << \
			  "value 1: " << (arg1) << \
			  "value 2: " << (arg2) << \
			  "value 3: " << (arg3) << \
			  "value 4: " << (arg4) << \
			  "value 5: " << (arg5) << "\n";
	


#if !defined (CHECK_FOR_NAN_GLOBALLY)
#define CHECK_FOR_NAN_GLOBALLY 1
#endif

#if !defined (CHECK_FP_EXCEPTIONS)
#define CHECK_FP_EXCEPTIONS 1
#if (GMS_COMPILED_BY_ICC) == 1
#include <fenv.h>
#else
#if defined (_WIN64) && defined (_WIN32)
#include <../../../Microsoft Visual Studio 12.0/VC/include/fenv.h>
#endif
#endif
#endif

#if !defined (SILENCE_COMPILER)
#define SILENCE_COMPILER 1
#endif

#if !defined (BE_VERBOSE)
#define BE_VERBOSE 1
#endif





#endif /*_GMS_CONFIG_H__*/
