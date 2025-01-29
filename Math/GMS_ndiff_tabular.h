

#ifndef __GMS_NDIFF_TABULAR_H__
#define __GMS_NDIFF_TABULAR_H__ 111020221056



namespace file_info {


        const unsigned int GMS_NDIFF_TABULAR_MAJOR = 1;

	const unsigned int GMS_NDIFF_TABULAR_MINOR = 1;

	const unsigned int GMS_NDIFF_TABULAR_MICRO = 0;

	const unsigned int GMS_NDIFF_TABULAR_FULLVER = 
		1000U*GMS_NDIFF_TABULAR_MAJOR+100U*GMS_NDIFF_TABULAR_MINOR+10U*GMS_NDIFF_TABULAR_MICRO;

	const char * const GMS_NDIFF_TABULAR_CREATE_DATE = "11-10-2022 09:56 +00200 (TUE 11 OCT 2022 GMT+2)";

	const char * const GMS_NDIFF_TABULAR_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_NDIFF_TABULAR_AUTHOR = "John Michael McNamee, York University Downsview, Ontario, Canada";

}


/*
       Purpose:
 !                    
 !                   NUMERTCAC DIFFERENTTATTON OF TARULAR FUNCTTON WITH AUTOMATIC
 !                   CHOICE OF OPTIMUM STEP SIZE. 
   Author: 
 !                           John Michael McNamee
 !
 !          Modified by:
 !                           Bernard Gingold
 !           
 !          References:
 !                            
 !                            NUMERICAL DIFFERENTIATION OF TABULATED
 !                            FUNCTIONS WITH AUTOMATIC CHOICE
 !                            OF STEP-SIZE
 !                            JOHN MICHAEL McNAMEE
 !                            York University Downsview, Ontario, Canada 
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
*/

#include <cstdint>
#include <cmath>
#include <algorithm>
#include "GMS_config.h"


namespace  gms {

        namespace math {
        

/*
     !n,            NUMBER OF POINTS AT WICH DERIVATIVE IS TO BE FOUND.
    !x0,           FIRST POINT AT WHICH DERIVATIVE IS TO BE FOUND
    !step,         INTERVAL BETWEEN POINTS AT WHICH DERIVATIVE IS TO BE found
    !h0,           INITIAL GUESS AT OPTIMUM STEP-SIZE(SEE ARTICLE 10
    !              OF' TEXT).
    !ne,           NUMRER OF POTNTS AT WHICH FUNCTTON IS TABULATED
    !xi,           FIRST POINT A WHICH FUNCTION IS TABULATED
    !ht,           INTERVAL AT WHICH FUNCTIONS IS TABULATED (MUST BE POSITIVE)
    !y,            (ARRAY OF DIMENSION AT LEAST NE). FUNCTION VALUES I.E.
    !              y(i) CONTAINS VALUE AT X=XI+(I-1)*HT
    !              (X0-XI),STEP,AND HO MUST BE MULTIPLES OF HT.
    !               X0 AND (X0+(N-1)*STEP) MUST LIE WITHIN RANGE OF TABLE
    !              Output
    !dy,           DERVATIVES AT REQUIRED POINTS
    !r1,           MAXIMUM ESTIMATED ERROR OVER THE WHOLE RANGE OF POINTS.
    !s1,           AVERAGE ESTIMATED ERROR OVER THE WHOLE RANGE OF POINTS.
    !npitl,        NUMBER OF POINTS AT WHICH INTERVAL OF TABULATION IS TOO
    !              LARGE.
    !nprts,        NUMBER OF POINTS AT WHICH RANGE OF TARLE IS TOO SMALL.
    !ierrfl,       A NON-ZERO VALUE INDICATES A FATAL ERROR AS FOLLOWS:
    !             -1 MEANS X DOES NOT COINCIDE WITH A TABULAR POINT
    !                (WITHIN RELATIVE TOLERANCE OF .00001).
    !             -2 MEANS X IS OUTSIDE RANGE OF TABLE.
    !             -3 MEANS HT IS NON-POSITIVE.
    !nerest,       NUMBER OF POINTS AT WHICH ERROR ESTIMATE COULD NOT BE
    !              MADE.
    !itmany,       NUMBER OF POINTS AT WHICH THERE WERE TOO MANY ITERATIONS
    !              (I.E., CHANGES OF STEP-SIZE).

       integer(kind=i4),               intent(in) :: n
       real(kind=sp),                  intent(in) :: x0
       real(kind=sp),                  intent(in) :: step
       real(kind=sp),                  intent(in) :: h0
       integer(kind=i4),               intent(in) :: ne
       real(kind=sp),                  intent(in) :: x1
       real(kind=sp),                  intent(in) :: ht
       real(kind=sp), dimension(1:ne), intent(in) :: y
       real(kind=sp), dimension(1:n),  intent(out):: dy
       real(kind=sp),                  intent(out):: r1
       real(kind=sp),                  intent(out):: s1
       integer(kind=i4),               intent(out):: npitl
       integer(kind=i4),               intent(out):: nprts
       integer(kind=i4),               intent(out):: ierrfl
       integer(kind=i4),               intent(out):: nerest
       integer(kind=i4),               intent(out):: itmany
*/
               
                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void diff(const int32_t n,
                                  const double x0,
                                  const double step,
                                  const double h0,
                                  const int32_t ne,
                                  const double x1,
                                  const double ht,
                                  const double * __restrict __ATTR_ALIGN__(64) y,
                                  double * __restrict __ATTR_ALIGN__(64) dy,
                                  double &r1,
                                  double &s1,
                                  int32_t &npitl,
                                  int32_t &nprts,
                                  int32_t &ierrfl,
                                  int32_t &nerest,
                                  int32_t &itmany); 
/*
   !z,  VALUE OF ABSCISSA.
   !x1, FIRST POINT AT WHICH FUNCTION IS TABULATED.
   !ht, INTERVAL OF TABULATION OF FUNCTION
   !y,  (ARRAY OF DIMENSION AT LEAST NE). TABLE OF VALUES OF FUNCTION.
*/
                        __ATTR_ALWAYS_INLINE__
                      
			static inline
                        double f(const double z,
                                 const double x1,
                                 const double ht,
                                 const double * __restrict y) {

                         double fz,t;
                         int32_t i;
                         fz = 0.0;
                         t  = (z-x1)/ht+1.0;
                         i  = std::abs(t)+0.5;
                         fz = y[i];
                         return (fz);
                   }  

/*
!x,  POINT AT WHICH DY/DX REQUIRED.
!h,  STEP-SIZE TO BE USED IN NUMERICAL DIFFERENTIATION FORMULAS.
!f,  TABLE-LOOKUP FUNCTION
!y,  ARRAY CONTAINING TABLE OY FUNCTION VALUES.
!x1, FIRST POINT AT WHICH FUNCTION IS TABULATED.
!ht, INTERVAL AT WHICH FUNCTION IS TABULATED.
!ne, NUMBER OF POINTS AT WHICH FUNCTION IS TABULATED.
!xe, LAST POINT AT WHICH FUNCTION IS TABULATED.
!d,  ESTIMATED VALUE OF DY/DX
!isf,FLAG TO INDICATE RANGE OF TABLE TOO SMALL.
!ief,FLAG TO INDICATE TABULAR INTERVAL LARGER THAN H.
!npitl, NUMBER OF POINTS AT WHICH TABULAR INTERVAT, TOO LARGE.
!nprts,NUMBER OF POINTS AT WHICH RANGE OF TABLE TOO SMALL.
*/
                        __ATTR_ALWAYS_INLINE__
                      
			static inline           
                        void deriv(const double x,
                                   const double h,
                                   double &d,
                                   int32_t &isf,
                                   int32_t &ief,
                                   const double * __restrict y,
                                   const double x1,
                                   const double ht,
                                   const int32_t ne,
                                   const double xe,
                                   int32_t &npitl,
                                   int32_t &nprts) {

                           double h,hc,t,t0,t1,t2,t3;
                           ief = 0;
                           hc  = ht*0.9;
                     l10:  if(h>hc) goto l20;
                           npitl=npitl+1;
                           h = ht;
                           ief = 1;
                           return;
                     l20:  t = x-4.0*h;
                           if(t<x1) goto l50;
                           t = x+4.0*h;
                           if(t>xe) goto l30;
                           t0 = 3.0*(f(x-4.0*h,x1,ht,y)-f(x+4.0*h,x1,ht,y));
                           t1 = 32.0*(f(x+3.0*h,x1,ht,y)-f(x-3.0*h,x1,ht,y));
                           t2 = 168.0*(f(x-2.0*h,x1,ht,y)-f(x+2.0*h,x1,ht,y));
                           t3 = 672.0*(f(x+h,x1,ht,y)-f(x-h,x1,ht,y));
                           d  = (t0+t1+t2+t3)/(840.0*h);  
                           return;
                     l30:  t = x-3.0*h;
                           if(t<x1) goto l40;
                           d = (11.0*f(x,x1,ht,y)-18.0*f(x-h,x1,ht,y)+ 
                               9.0*f(x-2.0*h,x1,ht,y)-2.0*f(x-3.0*h,x1,ht,y))/(6.0*h);
                           return;
                     l40:  nprts = nprts+1;
                           h = 0.5*h;
                           isf = 1;
                           goto l10;
                     l50:  t = x+3.0*h;
                           if(t>xe) goto l60;
                           d = (2.0*f(x+3.0*h,x1,ht,y)-9.0*f(x+2.0*h,x1,ht,y) 
                                + 18.0*f(x+h,x1,ht,y)-11.0*f(x,x1,ht,y))/(6.0*h);
                           return;
                     l60:  nprts = nprts+1;
                           h = 0.5*h;
                           isf = 1;
                           goto l10; 
                      
                      }  


              
                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void diff(const int32_t n,
                                  const float x0,
                                  const float step,
                                  const float h0,
                                  const int32_t ne,
                                  const float x1,
                                  const float ht,
                                  const float * __restrict __ATTR_ALIGN__(64) y,
                                  float * __restrict __ATTR_ALIGN__(64) dy,
                                  float &r1,
                                  float &s1,
                                  int32_t &npitl,
                                  int32_t &nprts,
                                  int32_t &ierrfl,
                                  int32_t &nerest,
                                  int32_t &itmany); 
                                  

                        __ATTR_ALWAYS_INLINE__
                      
			static inline
                        float f(const float z,
                                 const float x1,
                                 const float ht,
                                 const float * __restrict y) {

                         float fz,t;
                         int32_t i;
                         fz = 0.0f;
                         t  = (z-x1)/ht+1.0f;
                         i  = std::abs(t)+0.5f;
                         fz = y[i];
                         return (fz);
                   }  


                        __ATTR_ALWAYS_INLINE__
                       
			static inline           
                        void deriv(const float x,
                                   const float h,
                                   double &d,
                                   int32_t &isf,
                                   int32_t &ief,
                                   const float * __restrict y,
                                   const float x1,
                                   const float ht,
                                   const int32_t ne,
                                   const float xe,
                                   int32_t &npitl,
                                   int32_t &nprts) {

                           float h,hc,t,t0,t1,t2,t3;
                           ief = 0;
                           hc  = ht*0.9f;
                     l10:  if(h>hc) goto l20;
                           npitl=npitl+1;
                           h = ht;
                           ief = 1;
                           return;
                     l20:  t = x-4.0f*h;
                           if(t<x1) goto l50;
                           t = x+4.0f*h;
                           if(t>xe) goto l30;
                           t0 = 3.0f*(f(x-4.0f*h,x1,ht,y)-f(x+4.0f*h,x1,ht,y));
                           t1 = 32.0f*(f(x+3.0f*h,x1,ht,y)-f(x-3.0f*h,x1,ht,y));
                           t2 = 168.0f*(f(x-2.0f*h,x1,ht,y)-f(x+2.0f*h,x1,ht,y));
                           t3 = 672.0f*(f(x+h,x1,ht,y)-f(x-h,x1,ht,y));
                           d  = (t0+t1+t2+t3)/(840.0f*h);  
                           return;
                     l30:  t = x-3.0f*h;
                           if(t<x1) goto l40;
                           d = (11.0f*f(x,x1,ht,y)-18.0f*f(x-h,x1,ht,y)+ 
                               9.0f*f(x-2.0f*h,x1,ht,y)-2.0f*f(x-3.0f*h,x1,ht,y))/(6.0f*h);
                           return;
                     l40:  nprts = nprts+1;
                           h = 0.5f*h;
                           isf = 1;
                           goto l10;
                     l50:  t = x+3.0f*h;
                           if(t>xe) goto l60;
                           d = (2.0f*f(x+3.0f*h,x1,ht,y)-9.0f*f(x+2.0f*h,x1,ht,y) 
                                + 18.0f*f(x+h,x1,ht,y)-11.0f*f(x,x1,ht,y))/(6.0f*h);
                           return;
                     l60:  nprts = nprts+1;
                           h = 0.5f*h;
                           isf = 1;
                           goto l10; 
                      
                      }           

         

      } // math


} // gms





#endif /*__GMS_NDIFF_TABULAR_H__*/
