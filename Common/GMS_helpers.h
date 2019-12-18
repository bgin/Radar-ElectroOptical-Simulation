
#ifndef __GMS_HELPERS_H__
#define __GMS_HELPERS_H__





namespace file_info{


      const unsigned int gGMS_HELPERS_MAJOR = 1;
      const unsigned int gGMS_HELPERS_MINOR = 1;
      const unsigned int gGMS_HELPERS_MICRO = 0;
      const unsigned int gGMS_HELPERS_FULLVER =
      1000U*gGMS_HELPERS_MAJOR+100U*gGMS_HELPERS_MINOR+10U*gGMS_HELPERS_MICRO;
      const char * const pgGMS_HELPERS_CREATE_DATE = "Date: 08-04-2017 08 Apr 2017 , Time: 11:44 AM GMT+2 +00200";
      const char * const pgGMS_HELPERS_BUILD_DATE  = __DATE__ " " __TIME__;
      const char * const pgGMS_HELPERS_AUTHOR      = "Programmer: Bernard Gingold , contact: beniekg@gmail.com";
      const char * const pgGMS_HELPERS_SYNOPSIS    = "Various error condition checking and error raising helpers";
}

#include <type_traits>
#include "Config.h"



	
         /* 
		  *   Declaration of various error condition checking and error raising helpers.
		  *   
		 */


      /*
        @Description: Checks for C-style array's NULL pointer.
                      
        @Params:  source array of floats (32-bit, 24-bit precision)
        @Params:  none

        @Params:  destination array assumed to hold floats (32-bit, 24-bit precision)
        @Params:  none
        @Returns: bool (true on error)
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_2args( const float* __restrict,
                            float* __restrict,
			    const char * ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_2args( const float* __restrict,
                            float* __restrict,
			    const char * );
#endif		
     /*
        @Description: Checks for C-style array's NULL pointer.
                      
        @Params:  source array of doubles (64-bit, 53-bit precision)
        @Params:  none

        @Params:  destination array assumed to hold doubles (64-bit, 53-bit precision)
        @Params:  none
        @Returns: bool (true on error)
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER	
bool  check_null_ptr_2args( const double* __restrict,
                            double* __restrict,
			    const char *) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_2args( const double* __restrict,
                            double* __restrict,
			    const char *);

	 
     /*
        @Description: Checks for C-style array's NULL pointer.
                     
        @Params:  1st source array of floats (32-bit, 24-bit precision)
        @Params:  2nd source array of floats (32-bit, 24-bit precision)

        @Params:  destination array assumed to hold floats (32-bit, 24-bit precision)
        @Params:  executing function name, const char *
        @Returns: bool (true on error)
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_3args( const float* __restrict,
                            const float* __restrict, 
                            float* __restrict, const char *) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_3args( const float* __restrict,
                            const float* __restrict, 
                            float* __restrict, const char *);
#endif

	 
     /*
        @Description: Checks for C-style array's NULL pointer.
                      
        @Params:  1st source array of doubles (64-bit, 53-bit precision)
        @Params:  2nd source array of doubles (64-bit, 53-bit precision)

        @Params:  destination array assumed to hold doubles (64-bit, 53-bit precision)
        @Params:  none
        @Returns: bool (true on error)
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_3args( const double* __restrict,
                            const double* __restrict, 
                            double* __restrict, const char * ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_3args( const double* __restrict,
                            const double* __restrict, 
                            double* __restrict, const char * );
#endif

    /*
        @Description: Checks for C-style array's NULL pointer.
                     
        @Params:  source array of 32-bit signed integers
        @Params:  none

        @Params:  destination array assumed to hold 32-bit signed integers
        @Params:  const char * func_name , pointer to function name
        @Returns: bool (true on error)
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_2args( const int* __restrict,
                            int* __restrict,
			    const char *) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_2args( const int* __restrict,
                            int* __restrict,
			    const char *);
#endif

	
     /*
        @Description: Checks for C-style array's NULL pointer.
                       
        @Params:  source array of 64-bit signed integers
        @Params:  none

        @Params:  destination array assumed to hold 64-bit signed integers
        @Params:  const char * func_name , pointer to function name
        @Returns: bool (true on error)
        @Throws:  Nothing
        @Calls:   Nothing
*/
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_2args( const long long* __restrict,
                            long long* __restrict,
			    const char *)  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_2args( const long long* __restrict,
                            long long* __restrict,
			    const char *);
#endif
     /*
        @Description: Checks for C-style array's NULL pointer.
                     
        @Params:  1st source array of signed integers (32-bit)
        @Params:  2nd source array of signed integers (32-bit)

        @Params:  destination array assumed to hold signed integers (32-bit)
        @Params:  executing function name, const char *
        @Returns: bool (true on error)
        @Throws:  Nothing
      
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_3args(const int* __restrict,
                           const int* __restrict, 
                           int* __restrict, const char *)  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_3args(const int* __restrict,
                           const int* __restrict, 
                           int* __restrict, const char *);
#endif
     /*
        @Description: Checks for C-style array's NULL pointer.
                     
        @Params:  1st source array of signed integers (64-bit)
        @Params:  2nd source array of signed integers (64-bit)

        @Params:  destination array assumed to hold signed integers (64-bit)
        @Params:  executing function name, const char *
        @Returns: bool (true on error)
        @Throws:  Nothing
        
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
bool  check_null_ptr_3args(const long long* __restrict,
                           const long long* __restrict,
			   const long long* __restrict, const char *) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
bool  check_null_ptr_3args(const int* __restrict,
                           const int* __restrict, 
                           int* __restrict, const char *);
#endif



  


template<typename Ptr>  
std::enable_if<std::is_pointer<Ptr>::value,bool>
check_32alignment_2args( const Ptr* __restrict src, 
		         const Ptr* __restrict dst,
		         const char * func_name)   {

	 if (((reinterpret_cast<uintptr_t>(src)& 0x1F) != 0) ||
	     ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0)) {
#if (GMS_DEBUG_ON) == 1
		 std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Array(s) unaligned on 32-byte boundary in" << func_name << "\n"
			 << "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			 << "*****ERROR-DETAILS***** \n";
		   (reinterpret_cast<uintptr_t>(src)& 0x1F) == 0 ? std::cout << " src aligned on 32-byte boundary.\n" :
			 std::cout << " src unaligned on 32-byte boundary.\n";
		   (reinterpret_cast<uintptr_t>(dst)& 0x1F) == 0 ? std::cout << " dst aligned on 32-byte boundary.\n" :
			 std::cout << " dst unaligned on 32-byte boundary.\n";
		   
#endif
                 return (true);
	}
	return (false);
}



template<typename Ptr>
std::enable_if<std::is_pointer<Ptr>::value,bool> 
check_32alignment_3args(const Ptr* __restrict src1_ptr,
                        const Ptr* __restrict src2_ptr,
		        Ptr* __restrict dst_ptr,
			const char * func_name) {
	
	if (((reinterpret_cast<uintptr_t>(src1_ptr)& 0x1F) != 0) ||
		     ((reinterpret_cast<uintptr_t>(src2_ptr)& 0x1F) != 0) ||
		          ((reinterpret_cast<uintptr_t>(dst_ptr)& 0x1F) != 0)) {
#if (GMS_DEBUG_ON) == 1
		std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << "***NON-FATAL-ERROR***: Array(s) unaligned on 32-byte boundary in" << func_name << "\n"
			<< "at " << __FILE__ << ":" << __LINE__ << "(" << std::hex << "0x" << __FUNCTIONW__ << ")" << "\n"
			<< "*****ERROR-DETAILS***** \n";
		(reinterpret_cast<uintptr_t>(src1_ptr)& 0x1F) == 0 ? std::cout << " src1_ptr aligned on 32-byte boundary.\n" :
			std::cout << " src1_ptr unaligned on 32-byte boundary.\n";
		(reinterpret_cast<uintptr_t>(src2_ptr)& 0x1F) == 0 ? std::cout << " src2_ptr aligned on 32-byte boundary.\n" :
			std::cout << " src2_ptr unaligned on 32-byte boundary.\n";
		(reinterpret_cast<uintptr_t>(dst_ptr)& 0x1F) == 0 ? std::cout << " dst_ptr aligned on 32-byte boundary.\n" :
			std::cout << " dst_ptr unaligned on 32-byte boundary.\n";
#endif
                return (true);
	}
	return (false);
}

	
	 /*
		@Description: Dumps timing stats.

                      
        @Params:  TSC start value (wall clock)
        @Params:  TSC stop  value (wall clock)
        @Params:  Number of loop iterations
        @Params:  executing(calling) function name, const char *
        @Params:  none
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  dump_timing_stats(const unsigned __int64, const int,
		        const int, const int, const int,
                        const unsigned __int64, const int, const char * )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  dump_timing_stats(_In_ const unsigned __int64, _In_ const int,
		        _In_ const int, _In_ const int, _In_ const int,
                        _In_ const unsigned __int64, _In_ const int, _In_ const char * );
#endif

     /*
        @Description: Checks for input domain of x >= 0.0
					  Upon detecting error condition 'domain range'
                      indicator is set to true (1).
					  Manually vectorised version.
        @Params:  source array holding floats (32-bit)
        @Params:  lenght of source array
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Params:  none
        @Params:  none
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_lt_zero(const float* __restrict, 
			    const int, _Inout_ int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_lt_zero(_In_ const float* __restrict, 
			    _In_ const int, _Inout_ int* );
#endif

    /*
        @Description: Checks for input domain of x >= 0.0
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  2nd source array holding floats (32-bit)
        @Params:  length of both source arrays.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_lt_zero(const float* __restrict,
                            const float* __restrict,
			    const int,
			     int *    )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_lt_zero(_In_ const float* __restrict,
                            _In_ const float* __restrict,
			    _In_ const int,
			    _Inout_  int *    )	;
#endif


    /*
        @Description: Checks for input domain of x >= 0.0
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  source array holding floats (64-bit)
        @Params:  lenght of source array
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Params:  none
        @Params:  none
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_lt_zero(const double* __restrict,
	                    const int,
			    int *) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_lt_zero(_In_ const double* __restrict,
	                    _In_ const int,
			    _Inout_ int *);
#endif	 
     /*
        @Description: Checks for input domain of x >= 0.0
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding doubles (64-bit)
        @Params:  2nd source array holding doubles (64-bit)
        @Params:  length of both source arrays.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_lt_zero(const double* __restrict,
                            const double* __restrict,
			    const int,
			    int *   ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_lt_zero(_In_ const double* __restrict,
                            _In_ const double* __restrict,
			    _In_ const int,
			    _Inout_  int *   );
#endif
	  
     /*
        @Description: Checks for input domain of -1.0F <= x <= 1.0F
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  source array holding float (32-bit)
        @Params:  lenght of source array
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Params:  none
        @Params:  none
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
	 */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_bt_ones(const float* __restrict, 
			    const int,
			    int *) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_bt_ones(_In_ const float* __restrict, 
			    _In_ const int,
			    _Inout_ int *);
#endif
     /*
        @Description: Checks for input domain of -1.0F <= 0 <= 1.0F
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
       @Params:  1st source array holding doubles (32-bit)
       @Params:  2nd source array holding doubles (32-bit)
       @Params:  length of both source arrays.
       @Params:  none
       @Params:  pointer to variable holding 'Domain Range'error value.
       @Returns: Nothing
       @Throws:  Nothing
       @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_bt_ones(const float* __restrict,
                            const float* __restrict,
			    const int ,_Inout_  int*)  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_bt_ones(_In_ const float* __restrict,
                            _In_ const float* __restrict,
			    _In_ const int ,_Inout_  int*);
#endif

	 /*
        @Description: Checks for input domain of -1.0F <= x <= 1.0F
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  source array holding doubles (64-bit)
        @Params:  lenght of source array
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Params:  none
        @Params:  none
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_bt_ones(const double* __restrict, 
			    const int,
			    int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_bt_ones(_In_ const double* __restrict, 
			    _In_ const int,
			    _Inout_ int* );
#endif

    /*
        @Description: Checks for input domain of -1.0F <= 0 <= 1.0F
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding doubles (64-bit)
        @Params:  2nd source array holding doubles (64-bit)
        @Params:  length of both source arrays.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_bt_ones(const double* __restrict,
                            const double* __restrict,
			    const int, _Inout_  int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_bt_ones(_In_ const double* __restrict,
                            _In_ const double* __restrict,
			    _In_ const int, _Inout_  int* );
#endif

     /*
        @Description: Checks for input domain of  x != 0 && y != 0
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  2nd source array holding floats (32-bit)
        @Params:  length of both source arrays.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_ne_zero(const float* __restrict,
                            const float* __restrict,
			    const int ,
			    int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_ne_zero(_In_ const float* __restrict,
                            _In_ const float* __restrict,
			    _In_ const int ,
			    _Inout_  int* );
#endif
     /*
        @Description: Checks for input domain of  x != 0 && y != 0
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding doubles (64-bit)
        @Params:  2nd source array holding doubles (64-bit)
        @Params:  length of both source arrays.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
	 */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_ne_zero(const double* __restrict,
                            const double* __restrict,
			    const int,
			     int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined 
void  is_domain_f64_ne_zero(_In_ const double* __restrict,
                            _In_ const double* __restrict,
			    _In_ const int,
			    _Inout_  int* );
#endif


     /*
        @Description: Checks for input domain of  x != 0.F
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_ne_zero(const float* __restrict,
			    const int,
			    int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_ne_zero(_In_ const float* __restrict,
			    _In_ const int,
			    _Inout_ int* );
#endif

	 /*
        @Description: Checks for input domain of  x != 0.0L
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding doubles (64-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_ne_zero(const double* __restrict,
			    const int,
			    int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_ne_zero(_In_ const double* __restrict,
			    _In_ const int,
			    _Inout_ int* );
#endif
			    

     /*
        @Description: Checks for input domain of  x >= -1.F
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_ne_negone(const float* __restrict,
			      const int,
			      int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_ne_negone(_In_ const float* __restrict,
			      _In_ const int,
			      _Inout_ int* );
#endif
     /*
        @Description: Checks for input domain of  x >= -1.0
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_ne_negone(const double* __restrict,
			      const int,
			      int* ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_ne_negone(_In_ const double* __restrict,
			      _In_ const int,
			      _Inout_ int* );
#endif

	 
     /*
        @Description: Checks for existance of denormal values in
					  input range i.e. for some |x| < FLT_MIN
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_denormal_f32_present(const float* __restrict,
			      const int,
			      int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_denormal_f32_present(_In_ const float* __restrict,
			      _In_ const int,
			      _Inout_ int* );
#endif

     /*
        @Description: Checks for existance of denormal values in
                      input range i.e. for some |x| < DBL_MIN
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1).
                      Manually vectorised version.
        @Params:  1st source array holding doubles (64-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_denormal_f64_present(const double* __restrict,
			      const int,
			      int* ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_denormal_f64_present(_In_ const double* __restrict,
			      _In_ const int,
			      _Inout_ int* );
#endif
							  
	 /*
        @Description: Checks for existance of denormal values in
                      input range i.e. for some |x| < FLT_MIN
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1) and feraiseexcept(FE_DENORMAL)
					  is called.
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  check_denormf32_raise_except(const float* __restrict,
				   const int,
				   int* ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  check_denormf32_raise_except(_In_ const float* __restrict,
				   _In_ const int,
				   _Inout_ int* );
#endif

     /*
        @Description: Checks for existance of denormal values in
                      input range i.e. for some |x| < DBL_MIN
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1) and feraiseexcept(FE_DENORMAL)
                      is called.
                      Manually vectorised version.
        @Params:  1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  check_denormf64_raise_except(const double* __restrict,
				   const int,
				   int* )  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  check_denormf64_raise_except(_In_ const double* __restrict,
				   _In_ const int,
				   _Inout_ int* );
#endif

     /*
        @Description: Checks for existance of NaN and INF values
                      in input range i.e. for some |x| == NaN OR |x| == INF
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1) 
                      
                      Scalar version.
        @Params:      1st source array holding floats (32-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f32_nan_inf(_In_ const float* __restrict,
			    _In_ const int,
			    _Inout_ int*) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f32_nan_inf(_In_ const float* __restrict,
			    _In_ const int,
			    _Inout_ int*);
#endif

     /*
         @Description: Checks for existance of NaN and INF values
                      in input range i.e. for some |x| == NaN OR |x| == INF
                      Upon detecting error condition 'domain range'
                      indicator is set to true (1)

                     Scalar version
        @Params:      1st source array holding doubles (64-bit)
        @Params:  none
        @Params:  length of  source array.
        @Params:  none
        @Params:  pointer to variable holding 'Domain Range'error value.
        @Returns: Nothing
        @Throws:  Nothing
        @Calls:   Nothing
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
void  is_domain_f64_nan_inf(const double* __restrict,
			    const int,
			    int* ) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
void  is_domain_f64_nan_inf(_In_ const double* __restrict,
			    _In_ const int,
			    _Inout_ int* );
#endif
     /*
        @Description: Clears all floating-point state exceptions
                      Exception cleared:
					  FE_DENORMAL, FE_INVALID, FE_INEXACT, FE_UNDERFLOW
					  FE_OVERFLOW.
                      Scalar version
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
				  value from library 'feclearexcept' function
				  Non-zero value means error.
        @Throws:  Nothing
        @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  clear_all_feexcept() __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_all_feexcept();
#endif

     /*
        @Description: Clears only FE_DENORMAL exception
                      Exception cleared:
                      FE_DENORMAL
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
                  value from library 'feclearexcept' function
                  Non-zero value means error.
        @Throws:  Nothing
        @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  clear_fe_denormal() __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_fe_denormal();
#endif


	 /*
        @Description: Clears only FE_INEXACT exception
                      Exception cleared:
                      FE_INEXACT
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
                  value from library 'feclearexcept' function
                  Non-zero value means error.
        @Throws:  Nothing
        @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  clear_fe_inexact()  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_fe_inexact();
#endif

     /*
        @Description: Clears only FE_INVALID exception
                      Exception cleared:
                      FE_INVALID
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
                  value from library 'feclearexcept' function
                  Non-zero value means error.
        @Throws:  Nothing
        @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  clear_fe_invalid()  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_fe_invalid();
#endif
	 /*
        @Description: Clears only FE_DIVBYZERO exception
                      Exception cleared:
                      FE_DIVBYZERO
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
                  value from library 'feclearexcept' function
                  Non-zero value means error.
        @Throws:  Nothing
        @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  clear_fe_divbyezero() __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_fe_divbyzero();
#endif

     /*
        @Description: Clears only FE_OVERFLOW exception
                      Exception cleared:
                      FE_OVERFLOW
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
                  value from library 'feclearexcept' function
                  Non-zero value means error.
       @Throws:  Nothing
       @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER	
int  clear_fe_overflow() __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_fe_overflow();
#endif
	
     /*
        @Description: Clears only FE_UNDERFLOW exception
                      Exception cleared:
                      FE_UNDERFLOW
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'err' which indicates success or error as a return
                  value from library 'feclearexcept' function
                  Non-zero value means error.
        @Throws:  Nothing
        @Calls:   'feclearexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  clear_fe_underflow()  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  clear_fe_underflow();
#endif
     /*
        @Description: Tests all floating-point exceptions.
                      
        @Params:  All 7 floating-point exception types (exception values must be or'ed).
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function
                  
        @Throws:  Nothing
        @Calls:   'fetestexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  test_fe_excepts(const int) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  test_fe_excepts(const int);
#endif


     /*
        @Description: Tests for existance of FE_INVALID exception.

        @Params:  argument FE_INVALID macro.
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function

        @Throws:  Nothing
        @Calls:   'fetestexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  test_fe_invalid(const int) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  test_fe_invalid(const int);
#endif

     /*
        @Description: Tests for existance of FE_INEXACT exception.

        @Params:  argument FE_INEXACT macro.
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function

        @Throws:  Nothing
        @Calls:   'fetestexcept'
	 */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  test_fe_inexact(const int)  __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  test_fe_inexact(const int);
#endif

     /*
        @Description: Tests for existance of FE_DIVBYZERO exception.

        @Params:  argument FE_DIVBYZERO macro.
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function

        @Throws:  Nothing
        @Calls:   'fetestexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int test_fe_divbyzero(const int) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int test_fe_divbyzero(const int);
#endif

     /*
        @Description: Tests for existance of FE_UNNORMAL exception.

        @Params:  argument FE_UNNORMAL macro.
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function

        @Throws:  Nothing
        @Calls:   'fetestexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int test_fe_unnormal(const int) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int test_fe_unnormal(const int);
#endif

     /*
        @Description: Tests for existance of FE_UNDERFLOW exception.

        @Params:  argument FE_UNDERFLOW macro.
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function

        @Throws:  Nothing
        @Calls:   'fetestexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int test_fe_underflow(const int) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int test_fe_underflow(const int);
#endif

     /*
        @Description: Tests for existance of FE_OVERFLOW exception.

        @Params:  argument FE_OVERFLOW macro.
        @Params:  none
        @Params:  none
        @Params:  none
        @Params:  none
        @Returns: integer 'val' which indicates success or error as a return
                  value from library 'fetestexcept' function

        @Throws:  Nothing
        @Calls:   'fetestexcept'
     */
#if defined __GNUC__ || defined __INTEL_COMPILER
int  test_fe_overflow(const int) __attributes__((cold)) __attributes__((aligned(32)));
#elif defined _MSC_VER
int  test_fe_overflow(const int);
#endif

	


#endif /*__GMS_HELPERS_H__*/
