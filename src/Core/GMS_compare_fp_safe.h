
#ifndef __GMS_COMPARE_FP_SAFE_H__
#define __GMS_COMPARE_FP_SAFE_H__ 050920220913





namespace file_version 
{

    const unsigned int GMS_COMPARE_FP_SAFE_MAJOR = 1U;
    const unsigned int GMS_COMPARE_FP_SAFE_MINOR = 0U;
    const unsigned int GMS_COMPARE_FP_SAFE_MICRO = 0U;
    const unsigned int GMS_COMPARE_FP_SAFE_FULLVER =
      1000U*GMS_COMPARE_FP_SAFE_MAJOR+
      100U*GMS_COMPARE_FP_SAFE_MINOR+
      10U*GMS_COMPARE_FP_SAFE_MICRO;
    const char * const GMS_COMPARE_FP_SAFE_CREATION_DATE = "05-09-2022 09:13 AM +00200 (MON 05 SEP 2022 GMT+2)";
    const char * const GMS_COMPARE_FP_SAFE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPARE_FP_SAFE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPARE_FP_SAFE_DESCRIPTION   = "Based on book: Mathematical Theory of Electro-Optical Sensors (rus).";

}

#include <omp.h>
#include <cmath>
#include "GMS_config.h"

namespace gms
{


#pragma omp declare simd simdlen(16)
     bool approximatelyEqual(const float a,
		                     const float b,
					         const float epsilon) 
     {
			   const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
               return std::fabs(a - b) <=
                     ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
     }


#pragma omp declare simd simdlen(8)
	 bool approximatelyEqual(const double a,
		                     const double b,
					         const double epsilon) 
    {
			   const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
               return fabs(a - b) <=
                      ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
    }


#pragma omp declare simd simdlen(16)
     bool essentiallyEqual(const float a,
		                   const float b,
					       const float epsilon) 
    {
               const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return fabs(a - b) <=
			      ((fabsa > fabsb ? fabsb : fabsa) * epsilon);
    }


#pragma omp declare simd simdlen(8)
    bool essentiallyEqual(const double a,
		                  const double b,
					      const double epsilon) 
    {
               const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return fabs(a - b) <=
			       ((fabsa > fabsb ? fabsb : fabsa) * epsilon);
    }

#pragma omp declare simd simdlen(16)		   
    bool definitelyGreaterThan(const float a,
		                       const float b,
					           const float epsilon) 
    {
               const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return (a - b) >
			       ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
	}


#pragma omp declare simd simdlen(8)
	bool definitelyGreaterThan(const double a,
		                       const double b,
					           const double epsilon) 
    {
               const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return (a - b) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
	}


#pragma omp declare simd simdlen(16)
    bool definitelyLessThan(const float a,
		                    const float b,
					        const float epsilon) 
    {
               const float fabsa = std::fabs(a);
			   const float fabsb = std::fabs(b);
                           return (b - a) >
			    ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
    }


#pragma omp declare simd simdlen(8)
	bool definitelyLessThan(const double a,
		                    const double b,
					        const double epsilon) 
    {
               const double fabsa = std::fabs(a);
			   const double fabsb = std::fabs(b);
                           return (b - a) >
			   ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
    }

    

}















#endif /*__GMS_COMPARE_FP_SAFE_H__*/