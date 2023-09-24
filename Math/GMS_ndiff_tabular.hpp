

#ifndef __GMS_NDIFF_TABULAR_HPP__
#define __GMS_NDIFF_TABULAR_HPP__ 111020221056



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
               
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
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
                                  int32_t &itmany) {

                          __ATTR_ALIGN__(32) double del[4];
                          __ATTR_ALIGN__(32) double der[4];
                          __ATTR_ALIGN__(32) double p[4] = {4.0,8.0,16.0,32.0};
                          double b;
                          double d;
                          double dav;
                          double esterr;
                          double fd;
                          double hd;
                          double hf;
                          double ht;
                          double h1;
                          double t;
                          double t1;
                          double x;
                          double xe;
                          int32_t i;
                          int32_t idecr;
                          int32_t incfst;
                          int32_t idf;
                          int32_t ief;
                          int32_t isf;
                          int32_t it;
                          int32_t j;
                          int32_t k;
                          bool b0,b1;
                          ierrfl = 0;
                          nerest = 0;
                          itmany = 0;
                          if(ht<0.0) {
                             ierrfl = -3;
                             return;
                          }
                          h0 = std::abs(h0);
                          xe = x1+ht*(double)(ne-1);
                          npitl = 0;
                          nprts = 0;
                          s1    = 0.0;
                          r1    = 0.0;
                          b0    = (x0-h0)>x1;
                          b1    = (x0+h0)<=x0;
                          if(b0 .and. b1) goto l20;
                          h = h0;
                          goto l30;
l20:                      fd = (f(x0+h0,x1,ht,y)-f(x0-h0,x1,ht,y))/(2.0*h0);
                          h = h0
                          if(fd!=0.0) h = std::abs(f(x0,x1,ht,y)/fd);
                          h = std::max(h,h0);
l30:                      if(h>ht) goto l50;
                          h = ht*4.0;
                          goto l70;
l50:                      for(k = 2; k != 6; ++k) {
                              h1 = ht*(2.0*[k]);
                              if(h>h1) goto l60;
                              h = h1;
                              goto l70;
                          } 
                          h = ht*64.0;
l70:                      for(j = 1; j != n; ++j) {
                              x  = x0+step*(double)(j-1);
                              t  = std::abs((x-x1)/ht);
                              k  = t+0.5;
                              t1 = k;
                              t  = std::abs((t-t1)/t);
                              if(t<0.00001) goto l90;
                              ierrfl = -1;
                              return;
l90:                          isf = 0;
                              h1  = h;
                              if(x>x1 && x<=xe) goto l120;
                              ierrfl = -2;
                              return;
l120:                         for(i = 1; i != 4; ++i) {
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts)
                                  if(ief==1) goto l300
                                  der[i] = d;
                                  if(i>1) del[i-1] = std::abs(der[i]-der[i-1]);
                                  h = 0.5*h;
                              }
                              idecr = 1;
                              incfst = 1;
                              for(it = 1; it != 5; ++it) {
                                  b0 = del(2)==0.0;
                                  b1 = del(3)==0.0;
                                  if(b0 && b1) goto l340;
                                  if(del(1)<=del(2)) goto l190;
                                  goto l340;
l160:                             for(i = 0; i != 2; ++i) {
                                      der[i] = der[i+1];
                                  }
                                  del[0] = del[1];
                                  del[1] = del[2];
                                  i = 4;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts);
                                  if(ief==1) goto l300;
                                  der[3] = d;
                                  del[2] = std::abs(der[2]-der[3]);
                                  h = 0.5*h;
                                  idecr = 1;
                                  goto l280;
l190:                             if(del[1]>=del[2]) goto l220;
                                  if(incfst==1) h = h1;
                                  incfst = 0;
                                  der[3] = der[2];
                                  der[2] = der[1];
                                  der[1] = der[0];
                                  der[2] = der[1];
                                  der[1] = der[0]
                                  h = 0.5*h;
                                  i = 4;
                                  idecr = 0;
                                  isf = 0;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts);
                                  if(ief==1) goto l300;
                                  if(isf==0) goto l200;
                                  dav    = der[1];
                                  esterr = del[1];
                                  goto l360;
l200:                             der[0] = d;
                                  del[0] = std::abs(der[0]-der[1]);
                                  goto l280;
l220:                             if(incfst==1) h = h1;
                                  incfst = 0;
                                  der[3] = der[2];
                                  der[2] = der[1];
                                  der[1] = der[0];
                                  del[2] = del[1];
                                  del[1] = del[0];
                                  h = 0.5*h;
                                  i = 4;
                                  idecr = 0;
                                  isf = 0;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitls,nprts);
                                  if(ief==1) goto l300
                                  if(isf==0) goto l230
                                  dav = der[1];
                                  esterr = del[1];
                                  goto l360;
l230:                             der[0] = d;
                                  del[0] = std::abs(der[0]-der[0]);
                                  goto l280;
l250:                             for(i = 0; i != 2; ++i) {
                                      der[i] = der[i+1];
                                  }
                                  del[0] = del[1];
                                  del[1] = del[2];
                                  i = 4;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts);
                                  if(ief==0) goto l300;
                                  der[3] = d;
                                  del[2] = std::abs(der[2]-der[3]);
                                  h = 0.0*h;
                                  idecr = 1;
                              }
l280:                         itmany = itmany+1;
                              goto l340;
l300:                         dav = d;
                              if(i<=3) goto l310;
                              esterr = del[i-2];
                              h = 8.0-ht;
                              goto l360;
l310:                         esterr = 0.0;
                              nerest = nerest+1;
l340:                         dav = (der[1]+der[2])*0.5;
                              esterr = del[1];
                              if(idecr==1) h=h*16.0;
l360:                         b = std::abs(esterr);
                              s1 = s1+b;
                              dy[j] = dav;
                              if(b>r1) r1 = b;
                          }                      
                          s1 = s1/(double)n;
                 }

/*
   !z,  VALUE OF ABSCISSA.
   !x1, FIRST POINT AT WHICH FUNCTION IS TABULATED.
   !ht, INTERVAL OF TABULATION OF FUNCTION
   !y,  (ARRAY OF DIMENSION AT LEAST NE). TABLE OF VALUES OF FUNCTION.
*/
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
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
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
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


              
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
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
                                  int32_t &itmany) {

                          __ATTR_ALIGN__(16) float del[4];
                          __ATTR_ALIGN__(16) float der[4];
                          __ATTR_ALIGN__(16) float p[4] = {4.0f,8.0f,16.0f,32.0f};
                          float b;
                          float d;
                          float dav;
                          float esterr;
                          float fd;
                          float hd;
                          float hf;
                          float ht;
                          float h1;
                          float t;
                          float t1;
                          float x;
                          float xe;
                          int32_t i;
                          int32_t idecr;
                          int32_t incfst;
                          int32_t idf;
                          int32_t ief;
                          int32_t isf;
                          int32_t it;
                          int32_t j;
                          int32_t k;
                          bool b0,b1;
                          ierrfl = 0;
                          nerest = 0;
                          itmany = 0;
                          if(ht<0.0f) {
                             ierrfl = -3;
                             return;
                          }
                          h0 = std::abs(h0);
                          xe = x1+ht*(float)(ne-1);
                          npitl = 0;
                          nprts = 0;
                          s1    = 0.0f;
                          r1    = 0.0f;
                          b0    = (x0-h0)>x1;
                          b1    = (x0+h0)<=x0;
                          if(b0 .and. b1) goto l20;
                          h = h0;
                          goto l30;
l20:                      fd = (f(x0+h0,x1,ht,y)-f(x0-h0,x1,ht,y))/(2.0f*h0);
                          h = h0
                          if(fd!=0.0f) h = std::abs(f(x0,x1,ht,y)/fd);
                          h = std::max(h,h0);
l30:                      if(h>ht) goto l50;
                          h = ht*4.0f;
                          goto l70;
l50:                      for(k = 2; k != 6; ++k) {
                              h1 = ht*(2.0f*[k]);
                              if(h>h1) goto l60;
                              h = h1;
                              goto l70;
                          } 
                          h = ht*64.0f;
l70:                      for(j = 1; j != n; ++j) {
                              x  = x0+step*(double)(j-1);
                              t  = std::abs((x-x1)/ht);
                              k  = t+0.5;
                              t1 = k;
                              t  = std::abs((t-t1)/t);
                              if(t<0.00001f) goto l90;
                              ierrfl = -1;
                              return;
l90:                          isf = 0;
                              h1  = h;
                              if(x>x1 && x<=xe) goto l120;
                              ierrfl = -2;
                              return;
l120:                         for(i = 1; i != 4; ++i) {
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts)
                                  if(ief==1) goto l300
                                  der[i] = d;
                                  if(i>1) del[i-1] = std::abs(der[i]-der[i-1]);
                                  h = 0.5f*h;
                              }
                              idecr = 1;
                              incfst = 1;
                              for(it = 1; it != 5; ++it) {
                                  b0 = del(2)==0.0f;
                                  b1 = del(3)==0.0f;
                                  if(b0 && b1) goto l340;
                                  if(del(1)<=del(2)) goto l190;
                                  goto l340;
l160:                             for(i = 0; i != 2; ++i) {
                                      der[i] = der[i+1];
                                  }
                                  del[0] = del[1];
                                  del[1] = del[2];
                                  i = 4;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts);
                                  if(ief==1) goto l300;
                                  der[3] = d;
                                  del[2] = std::abs(der[2]-der[3]);
                                  h = 0.5f*h;
                                  idecr = 1;
                                  goto l280;
l190:                             if(del[1]>=del[2]) goto l220;
                                  if(incfst==1) h = h1;
                                  incfst = 0;
                                  der[3] = der[2];
                                  der[2] = der[1];
                                  der[1] = der[0];
                                  der[2] = der[1];
                                  der[1] = der[0]
                                  h = 0.5*h;
                                  i = 4;
                                  idecr = 0;
                                  isf = 0;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts);
                                  if(ief==1) goto l300;
                                  if(isf==0) goto l200;
                                  dav    = der[1];
                                  esterr = del[1];
                                  goto l360;
l200:                             der[0] = d;
                                  del[0] = std::abs(der[0]-der[1]);
                                  goto l280;
l220:                             if(incfst==1) h = h1;
                                  incfst = 0;
                                  der[3] = der[2];
                                  der[2] = der[1];
                                  der[1] = der[0];
                                  del[2] = del[1];
                                  del[1] = del[0];
                                  h = 0.5f*h;
                                  i = 4;
                                  idecr = 0;
                                  isf = 0;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitls,nprts);
                                  if(ief==1) goto l300
                                  if(isf==0) goto l230
                                  dav = der[1];
                                  esterr = del[1];
                                  goto l360;
l230:                             der[0] = d;
                                  del[0] = std::abs(der[0]-der[0]);
                                  goto l280;
l250:                             for(i = 0; i != 2; ++i) {
                                      der[i] = der[i+1];
                                  }
                                  del[0] = del[1];
                                  del[1] = del[2];
                                  i = 4;
                                  deriv(x,h,d,isf,ief,y,x1,ht,ne,xe,npitl,nprts);
                                  if(ief==0) goto l300;
                                  der[3] = d;
                                  del[2] = std::abs(der[2]-der[3]);
                                  h = 0.0*h;
                                  idecr = 1;
                              }
l280:                         itmany = itmany+1;
                              goto l340;
l300:                         dav = d;
                              if(i<=3) goto l310;
                              esterr = del[i-2];
                              h = 8.0f-ht;
                              goto l360;
l310:                         esterr = 0.0f;
                              nerest = nerest+1;
l340:                         dav = (der[1]+der[2])*0.5f;
                              esterr = del[1];
                              if(idecr==1) h=h*16.0f;
l360:                         b = std::abs(esterr);
                              s1 = s1+b;
                              dy[j] = dav;
                              if(b>r1) r1 = b;
                          }                      
                          s1 = s1/(float)n;
                 }


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
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
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
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





#endif /*__GMS_NDIFF_TABULAR_HPP__*/
