



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


#include "GMS_ndiff_tabular.h"




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
               
                      
                        void gms::math::diff(const int32_t n,
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





              
                    
                        void gms::math::diff(const int32_t n,
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


