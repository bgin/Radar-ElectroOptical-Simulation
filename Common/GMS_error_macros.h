
#ifndef __GMS_ERROR_MACROS_H__
#define __GMS_ERROR_MACROS_H__

#if !defined (GMS_ERROR_MACROS_MAJOR)
#define GMS_ERROR_MACROS_MAJOR 1
#endif

#if !defined (GMS_ERROR_MACROS_MINOR)
#define GMS_ERROR_MACROS_MINOR 0
#endif

#if !defined (GMS_ERROR_MACROS_MICRO)
#define GMS_ERROR_MACROS_MICRO 0
#endif

#if !defined (GMS_ERROR_MACROS_FULLVER)
#define GMS_ERROR_MACROS_FULLVER 1000
#endif

#if !defined (GMS_ERROR_MACROS_CREATE_DATE)
#define GMS_ERROR_MACROS_CREATE_DATE "28-09-2019 12:02 +00200 (SAT 28 SEP 2019 GMT+2)"
#endif
//
//  Set this value after successful build date/time
//
#if !defined (GMS_ERROR_MACROS_BUILD_DATE)
#define GMS_ERROR_MACROS_BUILD_DATE " "
#endif

#if !defined (GMS_ERROR_MACROS_AUTHOR)
#define GMS_ERROR_MACROS_AUTHOR "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com"
#endif

#if !defined (GMS_ERROR_MACROS_DESCRIPT)
#define GMS_ERROR_MACROS_DESCRIPT "Parametrized macros for displaying error messages."
#endif



#include <iostream>
#include <iomanip>
#if defined _WIN64
    #include "../LibAtmosModelCPP/System/StackWalker.h"
#endif

#if defined _WIN64
  #if !defined (DUMP_CALLSTACK_ON_ERROR)
      #define DUMP_CALLSTACK_ON_ERROR  \
	       StackWalker sw{};   \
	       sw.ShowCallstack();
  #endif
#endif

#if defined _WIN64
   #if !defined (PRINT_ERROR_INFO)
       #define PRINT_ERROR_INFO(msg) \
	       std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	       std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n";\
	       std::cerr << " Dumping callstack: \n"; \
	       StackWalker sw{};                      \
	       sw.ShowCallstack();
   #endif
#endif

#if defined _WIN64
   #if !defined (PRINT_ERROR_VALUE)
       #define PRINT_ERROR_VALUE(msg,value) \
	       std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << (value) << "\n"; \
	       std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	       std::cerr << " Dumping callstack: \n"; \
	       StackWalker sw{};                      \
	       sw.ShowCallstack();
    #endif
#endif

#if !defined (ABORT_ON_ERROR)
#define ABORT_ON_ERROR(msg,err) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	std::exit(err);
#endif	


#if !defined (THROW_ON_RUNTIME_ERROR)
#define THROW_ON_RUNTIME_ERROR(msg)  \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::runtime_error(msg);
#endif

#if !defined (THROW_ON_DOMAIN_ERROR)
#define THROW_ON_DOMAIN_ERROR(msg,value) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << " value: " << std::fixed << std::setprecision(15) << (value) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::domain_error(msg);
#endif

#if !defined (THROW_ON_OVERFLOW_ERROR)
#define THROW_ON_OVERFLOW_ERROR(msg,value) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << " value: " << std::fixed << std::setprecision << static_cast<double>((value)) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::overflow_error(msg);
#endif

#if !defined (THROW_ON_RANGE_ERROR)
#define THROW_ON_RANGE_ERROR(msg,value) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << " value: " << std::fixed << std::setprecision << static_cast<double>((value)) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::range_error(msg);
#endif

#if !defined (THROW_ON_INVALID_ARGUMENT_ERROR)
#define THROW_ON_INVALID_ARGUMENT_ERROR(msg) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::invalid_argument(msg);
#endif

#if !defined (THROW_ON_LOGIC_ERROR)
#define THROW_ON_LOGIC_ERROR(msg) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"; \
	throw std::logic_error(msg);
#endif

#if !defined (MALLOC_FAILED)
#define MALLOC_FAILED -1
#endif

#endif /*__GMS_ERROR_MACROS_H__*/
