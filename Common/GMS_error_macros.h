
#ifndef __GMS_ERROR_MACROS_H__
#define __GMS_ERROR_MACROS_H__ 280920191202







#include <iostream>
#include <iomanip>







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
