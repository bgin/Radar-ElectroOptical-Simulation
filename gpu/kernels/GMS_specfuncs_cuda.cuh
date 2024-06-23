
#ifndef __GMS_SPECFUNCS_CUDA_CUH__
#define __GMS_SPECFUNCS_CUDA_CUH__


#include <cuda/std/complex>

namespace file_info {

 const unsigned int GMS_SPECFUNCS_CUDA_MAJOR = 1U;
 const unsigned int GMS_SPECFUNCS_CUDA_MINOR = 0U;
 const unsigned int GMS_SPECFUNCS_CUDA_MICRO = 0U;
 const unsigned int GMS_SPECFUNCS_CUDA_FULLVER =
  1000U*GMS_SPECFUNCS_CUDA_MAJOR+100U*GMS_SPECFUNCS_CUDA_MINOR+10U*GMS_SPECFUNCS_CUDA_MICRO;
 const char * const GMS_SPECFUNCS_CUDA_CREATION_DATE = "07-04-2024 09:03AM +00200 (SUN 07 APR 2024 09:03 GMT+2)";
 const char * const GMS_SPECFUNCS_CUDA_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const GMS_SPECFUNCS_CUDA_AUTHOR        = "Shanjie Zhang, Jianming Jin, [Modified by Bernard Gingold, contact: beniekg@gmail.com]";
 const char * const GMS_SPECFUNCS_CUDA_SYNOPSIS      = "C-Cuda port of special functions library"


}

//using namespace cuda::std;

/*
   !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         specfuncs_cuda
 !          
 !          Purpose:
 !                        CUDA implementation of special functions library.
 !                        
 !                        
 !          History:
 !                        Date: 01-01-2024
 !                        Time: 09:46AM GMT+2
 !                        
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                     Shanjie Zhang, Jianming Jin
 ! 
 !          Modified:
                       Bernard Gingold 
 !          
 !                 
 !          References:
 !         
 !                     
 !                        Shanjie Zhang, Jianming Jin,
 !                        Computation of Special Functions,
 !                        Wiley, 1996,
 !                        ISBN: 0-471-11963-6,
 !                        LC: QA351.C45.
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
!==================================================================================85
*/


    

/*
    !*****************************************************************************80
!
!! AIRYA computes Airy functions and their derivatives.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!       
!  Parameters:
!
!    Input, float X, the argument of the Airy function.
!
!    Output, float AI, BI, AD, BD, the values of Ai(x), Bi(x),
!    Ai'(x), Bi'(x).
!
*/


      __device__ void airya(float x,
                            float  ai,
                            float  bi,
                            float  ad,
                            float  bd) {
            
           float vi1,vi2,vj1,vj2;
           float vk1,vk2,vy1,vy2;
           float xa,xq,z;
           constexpr float pir = 0.318309886183891f;
           constexpr float c1  = 0.355028053887817f;
           constexpr float c2  = 0.258819403792807f;
           constexpr float sr3 = 1.732050807568877f;
           xa = fabsf(x);
           z  = powf(xa,1.5f)*0.66666666666666666666667f;     
           xq = sqrtf(xa);
           ajyik(z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2);
           if(__builtin_expect((x==0.0f),0) {
              ai = c1;
              bi = sr3*c1;
              ad = -c2;
              bd = sr3*c2;
           }
           else if(0.0f<x) {
              ai = pir*xq/sr3*vk1;
              bi = xq*(__fmaf_ru(pir,vk1,2.0f/sr3*vi1)); //! pir * vk1 + 2.0 * invsr3 * vii
              ad = -xa/sr3*pir*vk2;
              bd = xa*(__fmaf_ru(pir,vk2,2.0f/sr3*vi2));
           }  
           else {
              ai = 0.5f*xq*(vj1-vy1/sr3);
              bi = -0.5f*xq*(vj1/sr3+vy1);
              ad = 0.5f*xa*(vj2+vy2/sr3);
              bd = 0.5f*xa*(vj2/sr3-vy2);
           }     
     }
     
     
/*
    !*****************************************************************************80
!
!! AIRYB computes Airy functions and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 June 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, float X, argument of Airy function.
!
!    Output, float AI, Ai(x).
!
!    Output, float BI, Bi(x).
!
!    Output, float AD, Ai'(x).
!
!    Output, float BD, Bi'(x).
!
*/


       __device__ void airyb(const float x,
                             float  ai,
                             float  bi,
                             float  ad,
                             float  bd) {
                             
             float ck[41];
             float dk[41];
             float df,dg,fx,gx;
             float r,rp,sad,sai;
             float sbd,sbi,sda,sdb;
             float ssa,ssb,xa,xar;
             float xcs,xe,xf,xm;
             float xp1,xq,xr1,xr2;
             float xss;
             int k,km;
             constexpr float eps = 1.0e-15f;
             constexpr float pi  = 3.14159265358979323846264f;
             constexpr float c1  = 0.355028053887817f;
             constexpr float c2  = 0.258819403792807f;
             constexpr float sr3 = 1.732050807568877f;
             xa = fabsf(x);
             if(__builtin_expect(x<=0.0f,0))
                xm = 8.0f;
             else
                xm = 0.5f;
             xq = sqrtf(xa);
             if(xa<=xm) {
                fx = 1.0f;
                r  = 1.0f;
                #pragma unroll
                for(k=1;k!=40;++k) {
                    float tk = (float)tk;
                    r  = r*x/(3.0f*tk)*x/(3.0f*tk-1.0f)*x;
                    fx = fx+r;
                    if ( fabsf(r)<fabsf(fx)*eps) break;
                }
                gx = x;
                r  = x;
                #pragma unroll
                for(k=1;k!=40;++k) {
                    float tk = (float)k;
                    r = r*x/(3.0f*tk)*x/(__fmaf_ru(3.0f,tk,1.0f))*x;
                    gx = gx+r;
                    if(fabsf(r)<fabsf(gx)*eps) break;
                }
                 ai = c1 * fx - c2 * gx;
                 bi = sr3 * ( __fmaf_ru(c1,fx ,c2 * gx) );
                 df = 0.5f * x * x;
                 r = df;
                 #pragma unroll
                for(k=1;k!=40;++k) {
                    float tk = (float)k;
                    r = r * x / ( 3.0f * tk ) * x / ( __fmaf_ru(3.0f,tk,2.0f) ) * x;
                    df = df + r; 
                    if (fabsf( r ) < fabsf( df ) * eps ) break; 
                }
                dg = 1.0f;
                r  = 1.0f;
                #pragma unroll
                for(k=1;k!=40;++k) {
                    float tk = (float)k;
                    r = r * x / ( 3.0f * tk ) * x / ( 3.0f * tk - 2.0f ) * x;
                    dg = dg + r;
                    if ( fabsf( r ) < fabsf( dg ) * eps ) break;
                }
                 ad = c1 * df - c2 * dg;
                 bd = sr3 * ( __fmaf_ru(c1,df,c2 * dg) );
             }
             else {
                 xe = xa * xq / 1.5f;
                 xr1 = 1.0f / xe;
                 xar = 1.0f / xq;
                 xf = sqrtf( xar );
                 rp = 0.5641895835477563f;
                 r = 1.0f;
                 #pragma unroll
                 for(k=1;k!=40;++k) {
                     float tk = (float)k;
                     r = r * ( 6.0f * tk - 1.0f ) * 
                     0.00462962962962962962963f * ( 6.0f * tk - 3.0f ) /
                     k * ( 6.0f * tk - 5.0f ) / ( 2.0f * tk - 1.0f );
                     ck[k] = r;
                     dk[k] = - ( 6.0f * tk + 1.0f ) / ( 6.0f * tk - 1.0f ) * ck[k];
                 }
                 km = (int)(24.5f-xa);
                 if(xa<6.0f)  km = 14;
                 if(15.0f<xa) km = 10;
                 if(0.0f<x) {
                    sai = 1.0f;
                    sad = 1.0f;
                    r   = 1.0f;
                    #pragma unroll
                    for(k=1;k!=km;++k) {
                        r = -r*xr1;
                        sai = sai + ck[k] * r;
                        sad = sad + dk[k] * r;
                    }
                    sbi = 1.0f;
                    sbd = 1.0f;
                    r   = 1.0f;
                    #pragma unroll
                    for(k=1;k!=km;++k) {
                         r = r * xr1;
                         sbi = sbi + ck[k] * r;
                         sbd = sbd + dk[k] * r;
                    }
                    xp1 = expf( - xe );
                    ai = 0.5f * rp * xf * xp1 * sai;
                    bi = rp * xf / xp1 * sbi;
                    ad = -0.5f * rp / xf * xp1 * sad;
                    bd = rp / xf / xp1 * sbd;
                  } 
                  else {
                    xcs = cosf( xe + pi * 0.25f );
                    xss = sinf( xe + pi * 0.25f );
                    ssa = 1.0f;
                    sda = 1.0f;
                    r = 1.0f;
                    xr2 = 1.0f / ( xe * xe );
                    #pragma unroll
                    for(k=1;k!=km;++k) {
                        r = - r * xr2;
                        ssa = ssa + ck[2*k] * r;
                        sda = sda + dk[2*k] * r;
                    }
                    ssb = ck[1] * xr1;
                    sdb = dk[1] * xr1;
                    r = xr1;
                    #pragma unroll
                    for(k=1;k!=km;++k) {
                        r = -r * xr2;
                        ssb = ssb + ck[2*k+1] * r;
                        sdb = sdb + dk[2*k+1] * r;
                    }
                    ai = rp * xf * ( xss * ssa - xcs * ssb );
                    bi = rp * xf * ( __fmaf_ru(xcs,ssa,xss * ssb) );
                    ad = -rp / xf * ( __fmaf_ru(xcs,sda,xss * sdb) );
                    bd =  rp / xf * ( xss * sda - xcs * sdb );
                  }
             }
      }
      
      
/*
    
      !*****************************************************************************80
!
!! AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
!
!  Discussion: 
!
!    Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions 
!    Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, float X, the argument.  X should not be zero.
!
!    Output, float VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2,
!    the values of J1/3(x), J2/3(x), Y1/3(x), Y2/3(x), I1/3(x), I2/3(x),
!    K1/3(x), K2/3(x).
!

*/     


       __device__ void ajyik(const float x,
                             float  vj1,
                             float  vj2,
                             float  vy1,
                             float  vy2,
                             float  vi1,
                             float  vi2,
                             float  vk1,
                             float  vk2) {
          
             float a0;
             float b0;
             float c0;
             float ck;
             float gn;
             float gn1;
             float gn2;
             float gp1;
             float gp2;
             int   k;
             int   k0;
             int   l;
             float pi;
             float pv1;
             float pv2;
             float px;
             float qx;
             float r;
             float rp;
             float rp2;
             float rq;
             float sk;
             float sum;
             float uj1;
             float uj2;
             float uu0;
             float vl;
             float vsl;
             float vv;
             float vv0;
             float vy1;
             float vy2;
             float x2;
             float xk; 
             
             pi = 3.141592653589793f;
             rp2 = 0.63661977236758f;
             gp1 = 0.892979511569249f;
             gp2 = 0.902745292950934f;
             gn1 = 1.3541179394264f;
             gn2 = 2.678938534707747f;
             vv0 = 0.444444444444444f;
             uu0 = 1.1547005383793f;
             x2 = x * x;

             if ( x < 35.0f )
                k0 = 12;
             else if ( x < 50.0f )
                k0 = 10;
             else
                k0 = 8;
             
             if ( x <= 12.0f ) {
             
                 for(l=1;l!=2;++l) {
                      vl = l * 0.333333333333333333333333333333333f;
                      vjl = 1.0f;
                      r = 1.0f;
                      #pragma unroll
                      for(k=1;k!=40;++k) {
                          const float kidx = (float)k;
                          r = -0.25f * r * x2 / ( kidx * ( kidx + vl ) );
                          vjl = vjl + r;
                          if ( fabsf( r ) < 1.0e-15f ) break;
                      }
                      a0 = powf(0.5 * x, vl);
                      if ( l == 1 )
                          vj1 = a0 / gp1 * vjl;
                      else
                          vj2 = a0 / gp2 * vjl;
                      }
             }
             else {
             
                  for(l=1;l!=2;++l) {
                       vv = vv0 * l * l;
                       px = 1.0f;
                       rp = 1.0f;
                       #pragma unroll
                       for(k=1;k!=k0;++k) {
                           const float kidx = (float)k;
                           float t0 = 4.0f * kidx - 3.0f;
                           float t1 = 4.0f * kidx - 1.0f;
                           rp = -0.78125e-02f * rp *
                           ( vv - t0*t0 ) * (vv - t1*t1) / 
                           ( kidx * ( 2.0f * kidx - 1.0f ) * x2 );
                           px += rp;
                       }
                       qx = 1.0f;
                       rq = 1.0f;
                       #pragma unroll
                       for(k=1;k!=k0;++k) {
                           const float kidx = (float)k;
                           float t0 = 4.0f * kidx - 1.0f;
                           float t1 = __fmaf_ru(4.0f,kidx,1.0f);
                           rq = -0.78125e-2f * rq * ( vv - t0*t0 ) *
                           ( vv - t1*t1 ) / ( k * ( __fmaf_ru(2.0f,kidx,1.0f) ) * x2 );
                           qx += rq;
                       }
                       qx = 0.125f * ( vv - 1.0f ) * qx / x;
                       xk = x - ( __fmaf_ru(0.5f,(l * 0.3333333333333333f),0.25f ) * pi;
                       a0 = sqrtf( rp2 / x );
                       ck = cosf( xk );
                       sk = sinf( xk );
                       if(l == 1) {
                           vj1 = a0 * ( px * ck - qx * sk );
                           vy1 = a0 * ( __fmaf_ru(px,sk,qx * ck) );
                       }
                       else {
                           vj2 = a0 * ( px * ck - qx * sk );
                           vy2 = a0 * ( __fmaf_ru(px,sk,qx * ck) )
                       }
                  }
             }
             
             if(x<=12.0f) {
                
                for(l=1;l!=2;++l) {
                    const float lidx = (float)l;
                    vl = lidx * 0.33333333333333333333333333333f;
                    vjl = 1.0f;
                    r = 1.0f;
                    #pragma unroll
                    for(k=1;k!=40;++k) {
                        const float kidx = (float)k;
                        r = -0.25f * r * x2 / ( kidx * ( kidx - vl ) );
                        vjl += r;
                        if(fabsf( r ) < 1.0e-15f ) break;
                    }
                    b0 = powf(2.0f/x,vl);
                    if( l == 1 )
                       uj1 = b0 * vjl / gn1;
                    else
                       uj2 = b0 * vjl / gn2;
     
                }
                   pv1 = pi * 0.33333333333333333333333f;
                   pv2 = pi * 0.66666666666666666666667f
                   vy1 = uu0 * ( vj1 * cosf( pv1 ) - uj1 );
                   vy2 = uu0 * ( vj2 * cosf( pv2 ) - uj2 );
             }
                
             if(x<=18.0f) {
                
                 for(l=1;l!=2;++l) {
                    const float lidx = (float)l;
                    vl = lidx * 0.33333333333333333333333333333f;
                    vil = 1.0f;
                    r = 1.0f;
                    #pragma unroll
                    for(k=1;k!=40;++k) {
                        const float kidx = (float)k;
                        r = 0.25f * r * x2 / ( kidx * ( kidx + vl ) );
                        vil += r;
                        if(fabsf( r ) < 1.0e-15f ) break;
                    }
                    a0 = powf(0.5f*x,vl);
                    if( l == 1 )
                       vi1 = a0 / gp1 * vil;
                    else
                       vi2 = a0 / gp2 * vil;
     
                }
                   
             }
             else {
                
                 c0 = expf( x ) / sqrtf( 2.0f * pi * x );
                 for(l=1;l!=2;++l) {
                     const float lidx = (float)l;
                     vv = vv0 * lidx * lidx;
                     vsl = 1.0f;
                     r = 1.0f;
                     #pragma unroll
                     for(k=1;k!=k0;++k) {
                         const float kidx = (float)k;
                         float t0 = 2.0f * k - 1.0f;
                         r = -0.125f * r * (vv - (t0*t0)) / (kidx*x);
                         vsl += r;
                     }
                     if(l == 1 )
                        vi1 = c0 * vsl;
                     else
                        vi2 = c0 * vsl;
                }
             }
             
            if(x<=9.0f) {
               
                for(l=1;l!=2;++l) {
                    const float lidx = (float)l;
                    vl = lidx * 0.3333333333333333333f;
                    if( l == 1 )
                       gn = gn1;
                    else
                       gn = gn2;
                    a0 = powf(2.0f / x,vl / gn);
                    sum = 1.0f;
                    r = 1.0f;
                    #pragma unroll
                    for(k=1;k!=60;++k) {
                        const float kidx = (float)k;
                        r = 0.25f * r * x2 / ( kidx * ( kidx - vl ) );
                        sum += r;
                        if(fabsf( r ) < 1.0e-15f) break; 
                    }
                   if( l == 1 ) 
                       vk1 = 0.5f * uu0 * pi * ( sum * a0 - vi1 );
                   else
                       vk2 = 0.5f * uu0 * pi * ( sum * a0 - vi2 );
    
               }
            }  
            else {
              
                 c0 = expf( -x ) * sqrtf(0.5f * pi / x );
                 for(l=1;l!=2;++l) {
                     const float lidx = (float)l;
                     vv = vv0 * lidx * lidx;
                     sum = 1.0f;
                     r = 1.0f;
                     #pragma unroll
                     for(k=1;k!=k0;++k) {
                         const float kidx = (float)k;
                         float t0 = 2.0f * kidx - 1.0f;
                         r = 0.125f * r * ( vv - t0*t0 ) / ( kidx * x );
                         sum += r;
                     }
                    if( l == 1 ) 
                       vk1 = c0 * sum;
                    else
                       vk2 = c0 * sum; 
                 }
            }
      }
      
      
/*
     !*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real(kind=sp) ::  GA, the value of the Gamma function.
!
*/

#include <cuda/std/limits>

        __device__ float gamma(const float x) {
                
                const float g[26] = {
                      1.0f, 
                      0.5772156649015329f,
                     -0.6558780715202538f,
                     -0.420026350340952e-01f
                      0.1665386113822915f, 
                     -0.421977345555443e-01f, 
                     -0.96219715278770e-02f, 
                      0.72189432466630e-02f, 
                     -0.11651675918591e-02f, 
                     -0.2152416741149e-03f, 
                      0.1280502823882e-03f,  
                     -0.201348547807e-04f, 
                     -0.12504934821e-05f, 
                      0.11330272320e-05f, 
                     -0.2056338417e-06f,  
                      0.61160950e-08f, 
                      0.50020075e-08f, 
                     -0.11812746e-08f, 
                      0.1043427e-09f,  
                      0.77823e-11f, 
                     -0.36968e-11f, 
                      0.51e-12f, 
                     -0.206e-13f, 
                     -0.54e-14f, 
                      0.14e-14f, 
                      0.1e-15f  
                };
                constexpr float pi = 3.14159265358979323846264f;
                float   gr;
                float   r;
                float   z;
                int32_t k;
                int32_t m;
                int32_t m1;
                if(x==truncf(x)) {
                
                   if(0.0f<x) {
                      ga = 1.0f;
                      m1 = (int32_t)(x-1.0f);
                      #pragma unroll
                      for(k=2;k!=m1;++k) {
                          const float kidx = (float)kidx;
                          ga *= kidx;
                      }
                   }
                   else {
                       ga = ::cuda::std::numeric_limits<float>::max();
                   }
                }
                else {
                
                    if(1.0f<fabsf(x)) {
                       z = fabsf(x);
                       m = (int32_t)z;
                       r = 1.0f;
                       #pragma unroll
                       for(k=1;k!=m;++k) {
                           const float kidx = (float)k;
                           r = r*z-kidx;
                       }
                       z = z - (float)m;
                    }
                    else {
                       z = x;
                    }
                    
                    gr = g[25];
                    #pragma unroll
                    for(k=25;k!=1;--k) {
                        float gk = g[k];
                        gr = __fmaf_ru(gr,z,gk);
                    }
                    ga = 1.0f/(gr*z);
                    if(1.0f<fabsf(x)) {
                       ga *= r;
                       if(x<0.0f) ga = -pi/(x*ga*sinf(pi*x));
                    }
                }
                
                return (ga);
        }


/*
       !*****************************************************************************80
!
!! BETA computes the Beta function B(p,q).
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, float P, Q, the parameters.
!    0 < P, 0 < Q.
!
!    Output, float BT, the value of B(P,Q).
!
*/


          __device__ float beta(const float p,
                                const float q) {
                   
                 float bt;
                 float gp;
                 float gpq;
                 float gq;
                 float ppq;
                 
                 gp = gamma(p);
                 gq = gamma(q);
                 ppq= p+q;
                 gpq= gamma(ppq);
                 bt = gp*gq/gpq;     
                 return (bt); 
        }
        
        
/*
      !*****************************************************************************80
!
!! CERF computes the error function and derivative for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex(kind=sp) :: z , the argument.
!
!    Output, complex(kind=sp) ::  CER, CDER, the values of erf(z) and erf'(z).
!    
*/


         __device__ void cerf(const cuda::std::complex<float> z,
                              cuda::std::complex<float> & cer,
                              cuda::std::complex<float> & cder) {
            
               cuda::std::complex<float> c0;
               cuda::std::complex<float> cs
               constexpr float pi   =  3.14159265358979323846264f;
               constexpr float eps  = 1.0e-12f;  
               float                     ei1;
               float                     ei2;
               float                     er;
               float                     er0;
               float                     er1;
               float                     er2;
               float                     eri;
               float                     err;
               float                     r;
               float                     ss;
               float                     w;
               float                     w1;
               float                     w2;
               float                     x;
               float                     x2;
               float                     y;
               int32_t                   k,n;
                
               x = real(z);
               x2= x*x;
               y = imag(z);
               if(x<=3.5f) {
                  
                  er = 1.0f;
                  r  = 1.0f;
                  for(k=1;k!=100;++k) {
                      float tk = (float)k;
                      r  = r*x2/tk+0.5f;
                      er += r;
                      if(fabsf(er-w)<=eps*fabsf(er)) break;
                      w = er;
                  }
                  c0 = 2.0f/sqrtf(pi)*x*expf(-x2);
                  er0 = real(c0*er);
               }
               else {
                  
                  er = 1.0f;
                  r  = 1.0f;
                  for(k=1;k!=12;++k) {
                      float tk = (float)k;
                      r  = -r*(tk-0.5f)/x2;
                      er += r;
                  }
                  c0 = expf(-x2)/(x*sqrtf(pi));
                  er0= 1.0f-real(c0*er);
               }
               
               if(y==0.0f) {
                  err = er0;
                  eri = 0.0f;
               }
               else {
                  
                  float t2    = 2.0f*pi*x;
                  float t0    = 2.0f*x*y;
                  cs          = cosf(t0);
                  float t1    = expf(-x2);
                  ss          = sinf(t0);
                  er1         = t1*(1.0f-cs)/t2;
                  ei1         = t1*ss/t2;
                  er2         = 0.0f;
                  for(n=1;n!=100;++n) {
                      float tn = (float)n;
                      float ny = tn*y;
                      float nn = tn*tn;
                      float t0 = expf(-0.25f*tn*tn);
                      float t1 = nn+4.0f*x2;
                      float t2 = 2.0f*x-2.0f*x*coshf(ny)*cs;
                      float t3 = tn*sinhf(ny)*ss;
                      float t4 = t2+t3;
                      er2      = er+t0/t1*t4;
                      if(fabsf((er2-w1)/er2)<eps) break;
                      w1 = er2;
                  }
                  
                  c0 = 2.0f*expf(-x2)/pi;
                  err= er0+er1+real(c0*er2);
                  ei2= 0.0f;
                  for(n=1;n!=100;++n) {
                      float tn = (float)n;
                      float ny = tn*y;
                      float nn = tn*tn;
                      float t0 = expf(-0.25f*tn*tn);
                      float t1 = nn+4.0f*x2;
                      float t2 = 2.0f*x*coshf(ny)*ss;
                      float t3 = tn*sinhf(ny)*cs;
                      float t4 = t2+t3;
                      ei2      = ei2+t0/t1*t4;
                      if(fabsf((ei2-w2)/ei2)<eps) break;
                      w2 = ei2;
                  }
                  
                  eri = ei1+real(c0*ei2);
               }
               
               cer = {err,eri};
               cder= 2.0f/sqrtf(pi)*exp(-z*z);     
       }


/*
      !*****************************************************************************80
!
!! CERROR computes the error function for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex(kind=sp) ::  Z, the argument.
!
!    Output, complex(kind=sp) ::  CER, the function value.
!
*/


        __device__ cuda::std::complex<float> 
                   cerror(const cuda::std::complex<float> z) {
           
              cuda::std::complex<float> c0,cl;
              cuda::std::complex<float> cr,cs;
              cuda::std::complex<float> z1;
              float                     a0;
              int32_t                   k;
              constexpr float           pi = 3.14159265358979323846264f;
              a0 = fabsf(z);
              c0 = exp(-z*z);
              z1 = z;
              
              if(real(z)<0.0f) z1 = -z;
              if(a0<=5.8f) {
                 
                 cs = z1;
                 cr = z1;
                 for(k=1;k!=120;++k) {
                     float tk = (float)k;
                     cr = cr*z1*z1/(tk+0.5f);
                     cs += cr;
                     if(abs(cr/cs)<1.0e-15f) break;
                 }
                 
                 cer = 2.0f*c0*cs/1.77245385090551602729817f;
              }
              else {
                 
                 cl = 1.0f/z1;
                 cr = cl;
                 for(k=1;k!=13;++k) {
                     float tk = (float)k;
                     cr = -cr*(tk-0.5f)/(z1*z1);
                     cl += cr;
                     if(abs(cr/cl)<1.0e-15f) break;
                 }
                 
                 cer = 1.0f-c0*cl/1.77245385090551602729817f;
              }  
              
              if(real(z)<0.0f) cer = -cer;
              return (cer);
                       
       }
       
       
/*
     !*****************************************************************************80
!
!! CFC computes the complex Fresnel integral C(z) and C'(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, complex(kind=sp) ::  Z, the argument.
!
!    Output, complex(kind=sp) ::  ZF, ZD, the values of C(z) and C'(z).
!

*/ 

      
        __device__ void cfc(const cuda::std::complex<float> z,
                            cuda::std::complex<float> & zf,
                            cuda::std::complex<float> & zd)  {
            
              cuda::std::complex<float> c,cf,cf0,cf1;
              cuda::std::complex<float> cg,cr,z0,zp;
              cuda::std::complex<float> zp2;
              float                     w0,wa,wa0;
              int                       k,m;
              constexpr float eps = 1.0e-14f;
              constexpr float pi  = 3.14159265358979323846264f;
              w0                  = abs(z);
              zp                  = 0.5f*pi*z*z;
              z0                  = {0.0f,0.0f};
              zp2                 = zp*zp;
              if(z==z0) {
                 c = z0;
              }
              else if(w0<=2.5f) {
                 cr = z;
                 c  = cr;
                 #pragma unroll
                 for(k=1; k!=80; ++k) {
                     float tk = (float)k;
                     float t0 = 4.0f*tk-3.0f;
                     float t1 = 2.0f*tk-1.0f;
                     float t2 = __fmaf_ru(4.0f,tk,1.0f);
                     float t3 = t0/tk/t1/t2;
                     cr       = -0.5f*t3*zp2;
                     c        += cr;
                     wa       = abs(c);
                     if(abs((wa-w0)/wa) < eps && 10 < k) break;
                     wa0 = wa;
                 }
              } 
              else if(2.5f<w0 && w0<4.5f) {
                     m   = 85;
                     c   = z0;
                     cf0 = {1.0e-30f,0.0f};
                     cf1 = z0;
                     #pragma unroll
                     for(k=m; k!=0; --k) {
                         float tk = (float)k;
                         float t0 = __fmaf_ru(2.0f,tk,3.0f);
                         cf       = t0*cf0/zp-cf1;
                         if(k==(int)((k/2)*2)) c += cf;
                         cf1 = cf0;
                         cf0 = c0; 
                     }
                     c = sqrt(2.0f/(pi*zp))*sin(zp)/cf*c;
              }   
              else {
                     cr = {1.0f,0.0f};
                     cf = {1.0f,0.0f};
                     cuda::std::complex<float> izp2 = 1.0f/zp2;
                     #pragma unroll
                     for(k=1; k!=20; ++k) {
                         float tk = (float)k;
                         float t0 = 4.0f*tk-1.0f;
                         float t1 = 4.0f*tk-3.0f;
                         cr = -0.25f*cr*t0*t1*izp2;
                         cf += cr;
                     }
                     cr = 1.0f/(pi*z*z);
                     cg = cr;
                     #pragma unroll
                     for(k=1; k!=12; ++k) {
                         float tk = (float)k;
                         float t0 = __fmaf_ru(4.0f,tk+1.0f);
                         float t1 = 4.0f*tk-1.0f;
                         cr = -0.25f*cr*t0*t1*izp2;
                         cg += cr;
                     }
                     c = 0.5f+(cf*sin(zp)-cg*cos(zp))/(pi*z);
              }   
              zf = c;
              zd = cos(0.5f*pi*z*z);
       }
       
       
/*
    !*****************************************************************************80
!
!! CFS computes the complex Fresnel integral S(z) and S'(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, complex(kind=sp) ::  Z, the argument.
!
!    Output, complex(kind=sp) ::  ZF, ZD, the values of S(z) and S'(z).
!
*/   


         __device__ void cfc(const cuda::std::complex<float> z,
                            cuda::std::complex<float> & zf,
                            cuda::std::complex<float> & zd)  {
                  
                 cuda::std::complex<float> cf,cf0,cf1,cg;
                 cuda::std::complex<float> cr,s,z0,zp,zp2;
                 float                     w0,wb,wb0;
                 int                       k,m;
                 constexpr float eps = 1.0e-14f;
                 constexpr float pi  = 3.14159265358979323846264f;  
                 w0                  = abs(z);
                 z0                  = {0.0f,0.0f};
                 zp                  = 0.5f*pi*z*z;
                 zp2                 = zp*zp;
                 if(z==z0) {
                    s = z0;
                 }
                 else if(w0<=2.5f) {
                    s  = z*zp*0.33333333333333333333333333333f;
                    cr = s;
                    for(k=1; k!=80; ++k) {
                        float tk = (float)k;
                        float t0 = 4.0f*tk-1.0f;
                        float t1 = __fmaf_ru(2.0f,tk,1.0f);
                        float t2 = __fmaf_ru(4.0f,tk,3.0f);
                        float t3 = t0/tk/t1/t2;
                        cr       = -0.5f*cr*t3*zp2;
                        s        += cr;
                        wb       = abs(s);
                        if(abs(wb-wb0)<eps && 10<k) break;
                        wb0 = wb;
                    } 
                 } 
                 else if(2.5f<w0 && w0<4.5f) {
                    m   = 85;
                    s   = z0;
                    cf1 = z0;
                    cf0 = {1.0e-30f,0.0f};
                    for(k=m; k!=0; --k) {
                        float tk = (float)k;
                        float t0 = __fmaf_ru(2.0f,tk,3.0f);
                        cf       = t0*cf0/zp-cf1;
                        if(k!=(int)((k/2)*2)) s += cf;
                        cf1 = cf0;
                        cf0 = cf;
                    } 
                    s = sqrt(2.0f/(pi*zp))*sin(zp)/cf*s;
                 }  
                 else {
                    cr = {1.0f,0.0f};
                    ci = {1.0f,0.0f};
                    cuda::std::complex<float> izp2 = 1.0f/zp2;
                    for(k=1; k!=20; ++k) {
                        float tk = (float)k;
                        float t0 = __fmaf_ru(4.0f,tk,1.0f);
                        float t1 = 4.0f*tk-3.0f;
                        cr       = -0.25f*cr*t0*t1*izp2;
                        cf       += cr;
                    }
                    cr = 1.0f/(pi*z*z);
                    cg = cr;
                    for(k=1; k!=12; ++k) {
                        float tk = (float)k;
                        float t0 = __fmaf_ru(4.0f,tk,1.0f);
                        float t1 = 4.0f*tk-1.0f;
                        cr       = -0.25f*cr*t0*t1*izp2;
                        cg       += cr;
                    }
                    s = 0.5f-(__fmaf_ru(cf,cos(zp),cg*sin(zp)))/(pi*z);
                 }  
                 zf = s;
                 zd = sin(0.5f*pi*z*z);              
        } 
        
        
/*
    
    !*****************************************************************************80
!
!! CGAMA computes the Gamma function for complex argument.
!
!  Discussion:
!
!    This procedcure computes the gamma function \E2(z) or ln[\E2(z)]
!    for a complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, Y, the real and imaginary parts of 
!    the argument Z.
!
!    Input, integer(kind=i4) ::  KF, the function code.
!    0 for ln[\E2(z)]
!    1 for \E2(z)
!
!    Output, real(kind=sp) ::  GR, GI, the real and imaginary parts of
!    the selected function.
!
*/ 


        __device__ void cgamma(float x,
                               float y,
                               const int   kf,
                               float &     gr,
                               float &     gi) {
            
            const float a[10] = {
                 8.333333333333333e-02f,-2.777777777777778e-03f, 
                 7.936507936507937e-04f, -5.952380952380952e-04f, 
                 8.417508417508418e-04f, -1.917526917526918e-03f,
                 6.410256410256410e-03f, -2.955065359477124e-02f,
                 1.796443723688307e-01f, -1.39243221690590f}; 
            float g0,gi1,gr1,si;
            float sr,t,th,th1;
            float th2,x0,x1,y1;
            float y2,z1,z2;
            int   j,k,na,idx;
            constexpr float pi = 3.14159265358979323846264f;
            
            if(x<0.0f) {
               x1 = x;
               y1 = y;
               x  = -x;
               y  = -y;
            }
            x0 = x;
            
            if(x<=7.0f) {
               na = (int)(7.0f-x);
               x0 = x+na;
            }
            
            z1 = sqrtf(__fmaf_ru(x0,x0,y*y));
            th = atanf(y/x0);
            gr = (x-0.5f)*logf(z1)-th*y-x0+
                 0.5f*logf(2.0f*pi);
            gi = __fmaf_ru(th,(x0-0.5f),y*log(z1)-y);
            
            idx = 0;
            for(k=1, k!=10; ++k) {
                ++idx;
                t = powf(z1,1-2*k);
                float tidx = (float)idx;
                float t0   = 2.0f*tidx-1.0f;
                gr += a[k]*t*cos(t0*th);
                gi -= a[k]*t*sin(t0*th);
            }
            
            if(x<=7.0f) {
               gr1 = 0.0f;
               gi1 = 0.0f;
               for(j=0; j!=na-1; ++j) {
                   float tj = (float)j;
                   float t0 = (x+tj)*(x+tj);
                   gr1 += 0.5f*logf(t0+y*y);
                   gi1 += atanf(y/(x+tj));
               }
               gr -= gr1;
               gi -= gi1;
            }
            
            if(x1<0.0f) {
               float t0 = pi*x;
               float t1 = pi*y;
               z1 = sqrtf(__fmaf_ru(x,x,y*y));
               th1= atanf(y/x);
               sr = -sinf(t0)*coshf(t1);
               si = -cosf(t0)*sinhf(t1);
               z2 = sqrtf(__fmaf_ru(sr,sr,si*si));
               if(sr<0.0f) th2 += pi;
               gr = logf(pi/(z1*z2))-gr;
               gi = -th1-th2-gi;
               x  = x1;
               y  = y1;
            }
            
            if(kf==1) {
               g0 = expf(gr);
               gr = g0*cosf(gi);
               gi = g0*sinf(gi);
            }
      }
      
      
/*
   
   !*****************************************************************************80
!
!! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes the modified Bessel functions I0(z), I1(z), 
!    K0(z), K1(z), and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex(kind=sp) Z, the argument.
!
!    Output, complex(kind=sp) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1, 
!    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z), 
!    and K1'(z).
!  
  
*/    

   
        __device__ void cik01(const cuda::std::complex<float> z,
                              cuda::std::complex<float> & cbi0,
                              cuda::std::complex<float> & cdi0,
                              cuda::std::complex<float> & cbi1,
                              cuda::std::complex<float> & cdi1,
                              cuda::std::complex<float> & cbk0,
                              cuda::std::complex<float> & cdk0,
                              cuda::std::complex<float> & cbk1) {
           
              const float a[12] = {
                  0.125f,               7.03125e-02f,
                  7.32421875e-02,       1.1215209960938e-01f,
                  2.2710800170898e-01f, 5.7250142097473e-01f,
                  1.7277275025845f,     6.0740420012735f,
                  2.4380529699556e+01f, 1.1001714026925e+02f,
                  5.5133589612202e+02f, 3.0380905109224e+03f};
             const float a1[10] = {
                   0.125f,              0.2109375f,
                   1.0986328125f,       1.1775970458984e+01f, 
                   2.1461706161499f,    5.9511522710323e+03f, 
                   2.3347645606175e+05f,1.2312234987631e+07f, 
                   8.401390346421e+08f, 7.2031420482627e+10f};
             const float b[12]   = {
                  -0.375f,              -1.171875e-01f, 
                  -1.025390625e-01f,    -1.4419555664063e-01f,
                  -2.7757644653320e-01f,-6.7659258842468e-01f,
                  -1.9935317337513f,    -6.8839142681099f,
                  -2.7248827311269e+01f, -1.2159789187654e+02f,
                  -6.0384407670507e+02f, -3.3022722944809e+03f};
              cuda::std::complex<float> ca,cb,ci,cr;
              cuda::std::complex<float> cs,cr,ct,cw;
              cuda::std::complex<float> z1,z2,zr,zr2;
              float                     w0;
              int                       k,k0;
              constexpr float pi = 3.14159265358979323846264f;
              a0                 = abs(z);
              ci                 = {0.0f,1.0f};
              z1                 = z;
              z2                 = z*z;
              
              if(real(z)<0.0f) z1 = -z;
              if(a0<=18.0f) {
                 cbi0 = {1.0f,0.0f};
                 cr   = {1.0f,0.0f};
                 #pragma unroll
                 for(k=1; k!=50; ++k) {
                     float tk = (float)k;
                     cr       =  0.25f*cr*z2/(tk*tk);
                     cbi0     += cr;
                     if(abs(cr/cbi0)<1.0e-15f) break;
                 }
                 cbi1 = {1.0f,0.0f};
                 cr   = {1.0f,0.0f};
                 for(k=1; k!=50; ++k) {
                     float tk = (float)k;
                     cr       =  0.25f*cr*z2/(tk*(tk+1.0f));
                     cbi1     += cr;
                     if(abs(cr/cbi1)<1.0e-15f) break;
                 }
                 cbi1 = 0.5f*z1*cbi1;
              }
              else {
              
                 if(a0<35.0f) {
                    k0 = 12;
                 }
                 else if(a0<50.0f) {
                    k0 = 9;
                 }
                 else {
                    k0 = 7;
                 }
                 
                 ca   = exp(z1)/sqrt(2.0f*pi*z1);
                 cbi0 = {1.0f,0.0f};
                 zr   = 1.0f/z1;
                 int idx = -1;
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     ++idx;
                     cbi0     = cbi0+a[idx]*pow(zr,tk);
                 }
                 cbi0 = ca+cbi0;
                 cbi1 = {1.0f,0.0f};
                 idx = -1;
                 for(k=1; k!=k0; ++k) {
                     ++idx;
                     float tk = (float)k;
                     cbi1     = __fmaf_ru(cbi1,b[idx],pow(zr,tk));
                 }
                 cbi1 = ca*cbi1;
              }
              
              if(a0<=9.0f) {
                 cs = {0.0f,0.0f};
                 ct = -log(0.5f*z1)-0.5772156649015329f;
                 w0 = 0.0f;
                 cr = {1.0f,0.0f};
                 for(k=1; k!=50; ++k) {
                     float tk = (float)k;
                     w0       = w0+1.0f/tk;
                     cr       = 0.25f*cr/(tk*tk)*z2;
                     cs       = cs+cr*(w0+ct);
                     if(abs((cs-cw)/cs)<1.0e-15f) break;
                     cw = cs;
                 }
                 cbk0 = ct+cs;
              }
              else {
                 cb  = 0.5f/z1;
                 zr2 = 1.0f/z2;
                 cbk0= {1.0f,0.0f};
                 int idx = -1;
                 for(k=1; k!=10; ++k) {
                     float tk = (float)k;
                     ++idx;
                     cbk0     = __fmaf_ru(cbk0,a1[idx],pow(zr2,tk));
                 }
                 cbk0 = cb*cbk0/cbi0;
              }
              
              cbk1 = (1.0f/z1-cbi1*cbk0)/cbi0;
              if(real(z)<0.0f) {
                 if(imag(z)<0.0f) {
                    cbk0 =  cbk0+ci*pi*cbi0;
                    cbk1 = -cbk1+ci*pi*cbi1;
                 }
                 else {
                    cbk0 = cbk0-ci*pi*cbi0;
                    cbk1 = -cbk1-ci*pi*cbi1;
                 }
                 cbi1 = -cbi1;
              }
              
              cdi0 = cbi1;
              cdi1 = cbi0-1.0f/z*cbi1;
              cdk0 = -cbk1;
              cdk1 = -cbk0-1.0f/z*cbk1;
              
      }
      
      
/*
    !*****************************************************************************80
!
!! CJK: asymptotic expansion coefficients for Bessel functions of large order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: KM, the maximum value of K.
!
!    Output, real(kind=sp) ::  A(L), the value of Cj(k) where j and k are 
!    related to L by L = j+1+[k*(k+1)]/2; j,k = 0,1,...,Km.
!
*/
      
        
        __device__ void cjk(const int km,
                            float * __restrict__ a) {
             
              float f,f0,g,g0;
              int   j,k,l1,l2;
              int   l3,l4;
              
              a[0] = 1.0f;
              f0   = 1.0f;
              g0   = 1.0f;
              
              for(k=0; k!=km; ++k) {
                  int   t0 = (k+1)*(k+2);
                  float tk = (float)k;
                  l1       = t0/2+1;
                  l2       = t0/2+k+2;
                  f        = (0.5f*tk+0.125f/(tk+1.0f))*f0;
                  a[l1]    = f;
                  f0       = f;
                  g        = -(1.5f*tk+0.625f/(3.0f*(tk+1.0f)))*g0;
                  a[l2]    = g;
                  g0       = g;
              } 
              
              for(k=1; k!=km-1; ++k) {
                   float tk = (float)k;
                  for(j=1; j!=k; ++j) {
                      l3       = k*(k+1)/2+j+1;
                      float tj = (float)j;
                      l4       = (k+1)*(k+2)/2+j+1;
                      float t0 = (__fmaf_ru(2.0f,tj,tk+1.0f));
                      float t1 = (tj+0.5f*tk+0.125f/t0)*a[l3];
                      float t2 = (tj+0.5f*tk-1.0f+0.625f/t0)*a[l3-1];
                      a[l4]    = t1-t2;
                  }
              }                   
       }
        
/*
   
      !*****************************************************************************80
!
!! CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
!
!  Discussion:
!
!    This procedure computes modified Bessel functions Iv(z) and
!    Kv(z) and their derivatives with a complex argument and a large order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) V, the order of Iv(z) and Kv(z).
!
!    Input, complex(kind=sp) Z, the argument.
!
!    Output, real(kind=sp) CBIV, CDIV, CBKV, CDKV, the values of
!    Iv(z), Iv'(z), Kv(z), Kv'(z).
!

*/


        __device__ void ciklv(const cuda::std::complex<float> z,
                              cuda::std::complex<float> &     cbiv,
                              cuda::std::complex<float> &     cdiv,
                              cuda::std::complex<float> &     cbkv,
                              cuda::std::complex<float> &     cdkv) {
           
                float a[91];
                cuda::std::complex<float> cf[12];
                int p1k[12] = {-1,1,-1,1,-1,1,-1,1,-1,1,-1,1};
                cuda::std::complex<float> ceta,cfi,cfk,csi;
                cuda::std::complex<float> csk,ct,ct2,cws;
                float                     v,v0,vr;
                int                       i,k,km;
                int                       l,l0,lf,idx,idx2;
                constexpr float pi = 3.14159265358979323846264f;
                km                 = 12;
                cjk(km,a);
                idx = -1;
                for(l=1; l!=0; --l) {
                    float tl = (float)l;
                    v0       = v-tl;
                    cuda::std::complex<float> zv0 = z/v0;
                    cws      = sqrt(1.0f+zv0*zv0);
                    ct       = 1.0f/cws;
                    ct2      = ct*ct;
                    for(k=1; k!=km; ++k) {
                        float tk = (float)k;
                        ++idx;
                        l0      = k*(k+1)/2+1;
                        lf      = l0+k;
                        cf[idx] = a[lf];
                        for(i=lf-1; i!=l0; --i) cf[idx] = __fmaf_ru(cf[idx],ct,a[i]);
                        cf[idx] = cf[idx]*pow(ct,tk);  
                      }
                    vr = 1.0f/v0;
                    csi= {1.0f,0.0f};
                    idx2 = -1;
                    for(k=1; k!=km; ++k) {
                        ++idx2;
                        csi = csi+cf[idx]*pow(vr,(float)k);
                    }
                    cbiv = sqrt(ct/(2.0f*pi*v0))*exp(v0*ceta)*csi;
                    if(l==1) cfi = cbiv;
                    csk = {1.0f,0.0f};
                    idx2 = -1;
                    for(k=1; k!=km; ++k) {
                        ++idx2;
                        csk = csk+p1k[idx]*cf[idx]*pow(vr,(float)k);
                    }
                    cbkv = sqrt(pi*ct/(2.0f*v0))*exp(-v0*ceta)*csk;
                    if(l==1) cfk = cbkv;
                }     
                
                cdiv = cfi-v/z*cbiv;
                cdkv = -cfk-v/z*cbkv;                 
      }
      
      
/*
     
    !*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, ?
!
!    Input, real(kind=sp) ::  X, ?
!
!    Output, real(kind=sp) ::  ENVJ, ?
! 

*/


        __device__ float envj(const int   n,
                              const float x) {
                              
                 float fn,l101,l102;
                 float result;
                 fn     = (float)n;
                 l101   = log10f(6.28f*fn);
                 l102   = log10f(1.36*x/fn);
                 result = 0.5*l101-fn*f102;
                 return (result);                       
       }
      
      
/*
   
    !*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward  
!    recurrence such that the magnitude of    
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Input, integer(kind=i4) :: MP, the negative logarithm of the 
!    desired magnitude.
!
!    Output, integer(kind=i4) :: MSTA1, the starting point.
! 
  
*/

        __device__ int msta1(const float x,
                             const int   mp) {
            
               float a0,f,f0,f1,fmp;
               int   it,n0,n1,nn;
               int result;
               fmp= (float)mp;
               a0 = fabsf(x);
               n0 = (int)(1.1e+00f*a0)+1;
               f0 = envj(n0,a0)-fmp; 
               n1 = n0+5;
               f1 = envj(n1,a0)-fmp;
               for(it=1; it!=20; ++it) {
                   nn = n1-(n1-n0)/(1.0f-f0/f1);
                   f  = envj(nn,a0)-fmp;
                   if(abs(nn-n1)<1) break;
                   n0 = n1;
                   f0 = f1;
                   n1 = nn;
                   f1 = f;
               }          
               result = nn;
               return (result);          
       }
       
       
/*

     !*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument of Jn(x).
!
!    Input, integer(kind=i4) :: N, the order of Jn(x).
!
!    Input, integer(kind=i4) :: MP, the number of significant digits.
!
!    Output, integer(kind=i4) :: MSTA2, the starting point.
!

*/


        __device__ int msta2(const float x,
                             const int   n,
                             const int   mp) {
            
             float a0,ejn,f,f0;
             float f1,hmp,obj,fmp;
             fmp = (float)mp;
             a0  = fabsf(x);
             hmp = 0.5f*fmp;
             ejn = envj(n,a0);
             if(ejn<=hmp) {
                obj = fmp;
                n0  = (int)(1.1f*a0);
             }  
             else {
                obj = hmp+ejn;
                n0  = n;
             }          
             f0 = envj(n0,a0)-obj;
             n1 = n0+5;
             f1 = envj(n1,a0)-obj;
             for(it=1; it!=20; ++it) {
                 nn = n1-(n1-n0)/(1.0f-f0/f1);
                 f  = envj(nn,a0)-obj;
                 if(abs(nn-n1)<1) break;
                 n0 = n1;
                 f0 = f1;
                 n1 = nn;
                 f1 = f;
             }   
             result = nn+10;
             return (result);       
      }
      
      
/*
     !*****************************************************************************80
!
!! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
!
!  Discussion:
!
!    This procedure also evaluates the derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    30 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order of In(z) and Kn(z).
!
!    Input, complex(kind=sp) ::  Z, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, complex(kind=sp) ::  CB((0:N), CDI(0:N), CBK(0:N), CDK(0:N), 
!    the values of In(z), In'(z), Kn(z), Kn'(z).
!
*/


        __device__ void ciknb(const int                       n,
                              const cuda::std::complex<float> z,
                              cuda::std::complex<float> * __restrict__ cbi,
                              cuda::std::complex<float> * __restrict__ cdi,
                              cuda::std::complex<float> * __restrict__ cbk,
                              cuda::std::complex<float> * __restrict__ cdk) {
              
             cuda::std::complex<float> c,cf,cf0,cf1;
             cuda::std::complex<float> cg,cg0,cg1;
             cuda::std::complex<float> ci,cr,cs0,csk0;
             cuda::std::complex<float> z1;
             float                     a0,fac,vt;
             int                       k,k0,l,m,nm;
             constexpr float pi = 3.14159265358979323846264f;
             constexpr float el = 0.57721566490153f;
             ci                 = {0.0f,1.0f};
             nm                 = n;
             a0                 abs(z);
             
             if(real(z)<0.0f) 
                z1 = -z;
             else
                z1 = z;
             if(n==0) nm=1;
             
             m = msta1(a0,200);
             if(m<nm)
                nm = m;
             else
                m  = msta2(a0,nm,15);
             
             cbs   = 0.0f;
             csk0  = 0.0f;
             cf0   = 0.0f;
             cf1   = 1.0e+30f;
             for(k=m; k!=0; --k) {
                 float tk = (float)k;
                 cf       = 2.0f*(tk+1.0f)*cf1/z1+cf0;
                 if(k<=nm) cbi[k] = cf;
                 if(k!=0 && k==2*(int)(k/2)) csk0 = csk0+4.0f*cf/tk;
                 cbs = cbs+2.0f*cf;
                 cf0 = cf1;
                 cf1 = cf;
             }    
             cs0 = exp(z1)/(cbs-cf);
             for(k=0; k!=nm; ++k) cbi[k] *= cs0;
             
             if(a0<=9.0f) {
                cbk[0] = __fmaf_ru(-(log(0.5f*z1)+el),cbi[0],cs0*csk0);
                cbk[1] =  (1.0f/z1-cbi[1]*cbk0[0])/cbi[0];
             }  
             else {
                ca0 = sqrt(pi/(2.0f*z1))*exp(-z1);
                if(a0<25.0f) 
                   k0 = 16;
                else if(a0<80.0f)
                   k0 = 10;
                else if(a0<200.0f)
                   k0 = 8;
                else
                   k0 = 6;
                
                for(l=0; l!=2; ++l) {
                    float tl = (float)l;
                    cbk1 = 1.0f;
                    vt   = 4.0f*tl;
                    cr   = {1.0f,0.0f};
                    for(k=1; k!=k0; ++k) {
                        float tk = (float)k;
                        float t0 = 2.0f*tk-1.0f;
                        cr       = 0.25f*cr*(vt-t0*t0)/(tk*z1);
                        cbkl     +=cr;
                    }
                    cbkl[l] = ca0*cbkl;
                }
             } 
             
             cg0 = cbk[0];
             cg1 = cbk[1];
             for(k=2; k!=nm; ++k) {
                 float tk = (float)k;
                 cg       = 2.0f*(tk-1.0f)/z1*cg1+cg0;
                 cbk[k]   = cg;
                 cg0      = cg1;
                 cg1      = cg; 
             }  
             if(real(z)<0.0f) {
                fac = 1.0f;
                if(imag(z)<0.0f)
                   for(k=0; k!=nm; ++k) {
                       cbk[k] = __fmaf_ru(fac,cbk[k],ci*pi*cbi[k]);
                       cbi[k] *=fac;
                       fac    = -fac;
                   }
                else
                   for(k=0; k!=nm; ++k) {
                       cbk[k] = fac*cbk[k]-ci*pi*cbi[k];
                       cbi[k] *=fac;
                       fac    = -fac;
                   }
             }   
             
             cdi[0] = cbi[1];
             cdk[0] = -cbk[1];
             for(k=1; k!=nm; ++k) {
                 float tk = (float)k;
                 cdi[k] = cbi[k-1]-tk/z*cbi[k];
                 cdk[k] = -cbk[k-1]-tk/z*cbk[k];
             }           
      }
      
      
/*
    !*****************************************************************************80
!
!! CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
!
!  Discussion:
!
!    Compute the modified Bessel functions Iv(z), Kv(z)
!    and their derivatives for an arbitrary order and
!    complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:       
!
!    Input, real ( kind = 8 ) V, the order of the functions.
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, real ( kind = 8 ) VM, the highest order computed.
!
!    Output, real ( kind = 8 ) CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
!    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
!
*/     


        __device__ void cikva(const float                     v,
                              const cuda::std::complex<float> z,
                              float                         & vm,
                              cuda::std::complex<float> * __restrict__ cbi,
                              cuda::std::complex<float> * __restrict__ cdi,
                              cuda::std::complex<float> * __restrict__ cbk,
                              cuda::std::complex<float> * __restrict__ cdk) {
              
              cuda::std::complex<float> ca,ca1,ca2,cb;
              cuda::std::complex<float> cbi0,cbk0,cbk1,cf;
              cuda::std::complex<float> cf1,cf2,cg0,cg1;
              cuda::std::complex<float> cgk,ci,ci0,cp;
              cuda::std::complex<float> cr,cr1,cr2,cs;
              cuda::std::complex<float> z1,z2;
              cuda::std::complex<float> csu,ct,cvk;
              float                     a0,gan,gap,piv,v;
              float                     v0,v0n,v0p,vm;
              float                     vt,w0,ws,ws0;
              int                       k,k0,m,n;
              constexpr float pi = 3.14159265358979323846264f;
              ci                 = {0.0f,1.0f};
              a0                 = fabsf(z);
              z1                 = z;
              z2                 = z*z;
              n                  = (int)v;
              v0                 = v-n;
              piv                = pi*v0;
              vt                 = 4.0f*v0*v0;
              
              if(n==0) n=1;
              if(a0<35.0f)
                 k0 = 14;
              else if(a0<50.0f)
                 k0 = 10;
              else
                 k0 = 8;
              if(real(z)<0.0f) z1 = -z;
              
              if(a0<18.0f) {
              
                 if(v0==0.0)
                    ca1 = {1.0f,0.0f};
                 else {
                    v0p = 1.0f+v0;
                    gap = gamma(v0p);
                    ca1 = cuda::std::pow(0.5f*z1,v0/gap);
                 }
                 
                 ci0 = {1.0f,0.0f};
                 cr  = {1.0f,0.0f};
                 #pragma unroll
                 for(k=1; k!=50; ++k) {
                     float tk = (float)k;
                     cr = 0.25f*cr*z2/(tk*(tk+v0));
                     ci0 += cr;
                     if(abs(cr)<abs(ci0)*1.0e-15) break;
                 }
                 
                 cbi0 = ci0*ca1;
              } 
              else {
                 
                 ca = exp(z1)/sqrt(2.0f*pi*z1);
                 cs = {1.0f,0.0f};
                 cr = {1.0f,0.0f};
                 #pragma unroll
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     float t0 = 2.0f*tk-1.0f;
                     cr = -0.125f*cr*(vt-t0*t0)/(tk*z1);
                     cs += cr;
                 }
                 cbi0 = ca*cs;
              }
              
              m = msta1(a0,200);
              
              if(m<n)
                 n = m;
              else
                 m = msta2(a0,n,15);
              
              cf2 = {0.0f,0.0f};
              cf1 = {(float)1.0e-30,0.0f};
              for(k=m; k!=0; --k) {
                  float tk = (float)k;
                  cf = 2.0f*(v0+tk+1.0f)/z1*cf1+cf2;
                  if(k<=n) cbi[k] = cf;
                  cf2 = cf1;
                  cf1 = cf;
              }
              cs = cbi0/cf;
              #pragma unroll
              for(k=0; k!=n; ++k) cbi[k] *= cs;
              
              if(a0<=9.0f) {
                 
                 if(v0==0.0) {
                    ct = -log(0.5f*z1)-0.5772156649015329f;
                    cs = {0.0f,0.0f};
                    w0 = 0.0f;
                    cr = {1.0f,0.0f};
                    #pragma unroll
                    for(k=1; k!=50; ++k) {
                        float tk = (float)k;
                        w0       = w0+1.0f/tk;
                        cr       = 0.25f*cr/(tk*tk)*z2;
                        cp       = cr*(w0+ct);
                        cs       = cs+cp;
                        if(10<=k && abs(cp/cs)<1.0e-15f) break;
                    }
                    cbk0 = ct+cs;
                 }
                 else {
                    
                    v0n = 1.0f-v0;
                    gan = gamma(v0n);
                    ca2 = 1.0f/(gan*pow(0.5f*z1,v0));
                    ca1 = pow(0.5f*z1,v0/gap);
                    csu = ca2-ca1;
                    cr1 = {1.0f,0.0f};
                    cr2 = {1.0f,0.0f};
                    #pragma unroll
                    for(k=1; k!=50; ++k) {
                        float tk = (float)k;
                        cr1      = 0.25f*cr1*z2/(tk*(tk-v0));
                        cr2      = 0.25f*cr2*z2/(tk*(tk+v0));
                        csu      = csu+ca2*cr1-ca1*cr2;
                        ws       = abs(csu);
                        if(10<=k && abs(ws-ws0)/ws<1.0e-15f) break;
                        ws0 = ws;
                    }
                    cbk0 = 0.5f*pi*csu/sin(piv);
                 }
              }
              else {
                 
                 cb = exp(-z1)*sqrt(0.5f*pi/z1);
                 cs = {1.0f,0.0f};
                 cr = {1.0f,0.0f};
                 #pragma unroll
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     float t0 = 2.0f*tk-1.0f;
                     cr = 0.125f*cr*(vt-t0*t0)/(tk*z1);
                     cs += cr;
                 }
                 cbk0 = cb*cs;
              }
              
              cbk1   = (1.0f/z1-cbi[1]*cbk0)/cbi[0];
              cbk[0] = cbk0;
              cbk[1] = cbk1;
              cg0    = cbk0;
              cg1    = cbk1;
              #pragma unroll
              for(k=2; k!=n; ++k) {
                  float tk = (float)k;
                  cgk      = 2.0f*(v0+tk-1.0f)/z1*cg1+cg0;
                  cbk[k]   = cgk;
                  cg0      = cg1;
                  cg1      = cgk
              }
              
              if(real(z)<0.0f) {
                 if(imag(z)<0.0f) {
                    #pragma unroll
                    for(k=0; k!=n; ++k) {
                        float tk = (float)k;
                        cvk      = exp((tk+v0)*pi*ci);
                        cbk[k]   = __fmaf_ru(cvk,cbk[k],pi*ci*cbi[k]);
                        cbi[k]   = cbi[k]/cvk;
                    }
                 }
                 else if(0.0f<imag(z)) {
                     #pragma unroll
                     for(k=0; k!=n; ++k) {
                         float tk = (float)k;
                         cvk      = exp((tk+v0)*pi*ci);
                         cbk[k]   = cbk[k]/cvk-pi*ci*cbi[k];
                         cbi[k]   = cvk*cbi[k];
                     } 
                 }
              }
              
              cdi[0] = v0/z*cbi[0]+cbi[1];
              cdk[0] = v0/z*cbk[0]-cbk[1];
              #pragma unroll
              for(k=1; k!=n; ++k) {
                  float tk = (float)k;
                  cdi[k] = -(tk+v0)/z*cbi[k]+cbi[k-1];
                  cdk[k] = -(tk+v0)/z*cbk[k]-cbk[k-1];
              }
              vm = n+v0;
       }
       
       
/*
   
     !*****************************************************************************80
!
!! IKNB compute Bessel function In(x) and Kn(x).
!
!  Discussion:
!
!    Compute modified Bessel functions In(x) and Kn(x),
!    and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order of In(x) and Kn(x).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  BI(0:N), DI(0:N), BK(0:N), DK(0:N),
!    the values of In(x), In'(x), Kn(x), Kn'(x).
!

*/
       
        __device__ void iknb(const int   n,
                             const float x,
                             int        &nm,
                             float * __restrict__ bi,
                             float * __restrict__ di,
                             float * __restrict__ bk,
                             float * __restrict__ dk) {
                             
                float a0,bkl,bs,f;
                float f0,f1,g,g0;
                float r,s0,sk0,vt;
                float invx;
                int   k,k0,l,m;
                constexpr float pi = 3.14159265358979323846264f; 
                constexpr float el = 0.5772156649015329f;
                nm                 = n;
                
                if(n==0) nm = 1;
                m = msta1(x,200);
                if(m<nm)
                   nm = m;
                else
                   nm = msta2(x,nm,15);
                bs  = 0.0f;
                sk0 = 0.0f;
                f0  = 0.0f;
                f1  = 1.1754943508e-38f;
                invx= 1.0f/x;
                for(k=m; k>=0; --k) {
                    float tk = (float)k;
                    f        = 2.0f*(tk+1.0f)*invx*f1+f0;
                    if(k<=nm) bi[k] = f;
                    if(k!=0 && k==2*(int)(k/2))
                       sk0 = sk0+4.0f*f/tk;
                    bs = __fmaf_ru(f,2.0f,bs);
                    f0 = f1;
                    f1 = f;
                }    
                  s0 = expf(x)/(bs-f);
                  for(k=0; k!=nm; ++k) s0 *= bi[k];
                  
                if(x<=8.0f) {
                   bk[0] = __fmaf_ru(-(logf(0.5f*x)+el),bi[0],s0*sk0);
                   bk[1] = (invx-bi[1]*bk[0])/bi[0];
                } 
                else {
                   a0 = sqrtf(pi/(2.0f*x))*expf(-x); 
                } 
                  if(x<25.0f)
                     k0 = 16;
                  else if(x<80.0f)
                     k0 = 10;
                  else if(x<200.0f) 
                     k0 = 8;
                  else
                     k0 = 6;
                  for(l=0; l<=1; ++l) {
                      float tl = (float)l;
                      bkl      = 1.0f;
                      vt       = 4.0f*tl;
                      r        = 1.0f;
                      for(k=1; k<=k0; ++k) {
                          float tk = (float)k;
                          float t0 = 2.0f*tk-1.0f;
                          r        = 0.125f*r*(vt-t0*t0)/(tk*x);
                          bkl      +=r;
                      }
                     bkl[l] *= a0;
                  }
               }
                g0 = bk[0];
                g1 = bk[1];
                for(k=2; k<=nm; ++k) {
                    float tk = (float)k;
                    g        = 2.0f*(tk-1.0f)*__fmaf_ru(invx,g1,g0);
                    bk[k]    = g;
                    g0       = g1;
                    g1       = g;
                }
                di[0] = bi[1];
                dk[0] = -bk[1];
                for(k=1; k<=nm; ++k) {
                    float tk = (float)k;
                    di[k] =  bi[k-1]-tk/x*bi[k];
                    dk[k] = -bk[k-1]-tk/x*bk[k];
                }
       }
       
/*     
!*****************************************************************************80
!
!! IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
!
!  Discussion:
!
!    This procedure computes modified Bessel functions I0(x), I1(x),
!    K0(x) and K1(x), and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
!    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
!
*/


         __device__ void ik01a(const float x,
                               float      &bi0,
                               float      &di0,
                               float      &bi1,
                               float      &di1,
                               float      &bk0,
                               float      &dk0,
                               float      &bk1,
                               float      &dk1) {
                               
                const float a[12] = {
                     0.125f, 7.03125e-02f, 
                     7.32421875e-02f, 1.1215209960938e-01f, 
                     2.2710800170898e-01, 5.7250142097473e-01f, 
                     1.7277275025845f, 6.0740420012735f, 
                     2.4380529699556e+01f, 1.1001714026925e+02f,
                     5.5133589612202e+02f, 3.0380905109224e+03f};    
               const float a1[8] =  {
                     0.125f, 0.2109375f, 
                     1.0986328125f, 1.1775970458984e+01f, 
                     2.1461706161499e+02f, 5.9511522710323e+03f,
                     2.3347645606175e+05f, 1.2312234987631e+07f};  
               const float b[12] =  {
                     -0.375f, -1.171875e-01f, 
                     -1.025390625e-01f, -1.4419555664063e-01f, 
                     -2.7757644653320e-01f, -6.7659258842468e-01f,
                     -1.9935317337513f, -6.8839142681099f,
                     -2.7248827311269e+01f, -1.2159789187654e+02f,
                     -6.0384407670507e+02f, -3.3022722944809e+03f};   
               float ca,cb,ct,r;
               float w0,ww,x2,xr;
               float xr2;
               int   idx,k,k0;
               constexpr float pi = 3.14159265358979323846264f;
               constexpr float el = 0.5772156649015329f;
               x2                 = x*x;
               if(x<=18.0f) {
                  bi0 = 1.0f;
                  r   = 1.0f;
                  for(k=1; k<=50; ++k) {
                      float tk = (float)k;
                      r        =  0.25f*r*x2/(tk*tk);
                      bi0      += r;
                      if(fabsf(r/bi0)<1.0e-15f) break;
                  }
                  bi1 = 1.0f;
                  r   = 1.0f;
                  for(k=1; k<=50; ++k) {
                      float tk = (float)k;
                      r        =  0.25f*r*x2/(tk*(tk+1.0f));
                      bi1      += r;
                      if(fabsf(r/bi1)<1.0e-15f) break;
                  }
                    bi1 = 0.5f*x*bi1;
               } 
               else {
                   ca   = expf(x)/sqrtf(2.0f*pi*x);
                   bi0  = 1.0f;
                   xr   = 1.0f/x;
                   if(x<35.0f) 
                     k0 = 12;
                   else if(x<50.0f)
                      k0 = 9;
                   else
                     k0 = 7;
                   idx = -1;
                   for(k=1; k<=k0; ++k) {
                       ++idx;
                       bi0 = __fmaf_ru(powf(xr,k),a[idx],bi0);
                   }
                    bi0 = ca*bi0;
                    bi1 = 1.0f;
                    idx = -1;
                    for(k=1; k<=k0; ++k) {
                       ++idx;
                       bi1 = __fmaf_ru(powf(xr,k),b[idx],bi1);
                   }
                    bi1 = ca*bi1;
               }     
              if(x<=9.0f) {
                 ct  = -(logf(x*0.5f)+el);
                 bk0 = 0.0f;
                 w0  = 0.0f;
                 r   = 1.0f;
                 for(k=1; k<=50; ++k) {
                     float tk = (float)k;
                     w0       = w0+1.0f/tk;
                     r        = 0.25f*r/(tk*tk)*x2;
                     bk0      = __fmaf_ru(w0+ct,r,bk0);
                     if(fabsf((bk0-ww)/bk0)<1.0e-15f) break;
                     ww       = bk0;
                 }
                 bk0 += ct;
              } 
              else {
                 cb  = 0.5f/x;
                 xr2 = 1.0f/x2;
                 bk0 = 1.0f;
                 idx = -1;
                 for(k=1; k<=8; ++k) {
                     ++idx;
                     bk0 = __fmaf_ru(powf(xr2,k),a1[idx],bk0);
                 }
                 bk0 = cb*bk0/bi0;
              }       
               bk1 = (1.0f/x-bi1*bk0)/bi0;
               di0 = bi1;
               di1 = bi0-bi1/x;
               dk0 = -bk1;
               dk1 = -bk0-bk1/x;
       }
       
       
/*

    !*****************************************************************************80
!
!! IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
!    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
!

*/


        __device__ void ik01b( const float x,
                               float      &bi0,
                               float      &di0,
                               float      &bi1,
                               float      &di1,
                               float      &bk0,
                               float      &dk0,
                               float      &bk1,
                               float      &dk1) {
            
            float t0;
            if(x==0.0f) return;
            if(x<=3.75f) {
               t   = x/3.75f;
               t2  = t*t;
               bi0 = __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(0.0045813f,t2,0.0360768f),
                                          t2,0.2659732f),
                                          t2,1.2067492f),
                                          t2,3.0899424f),
                                          t2,3.5156229f),
                                          t2,1.0f);
               bi1 = x*__fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(0.00032411f,t2,0.00301532f),
                                             t2,0.02658733f),
                                             t2,0.15084934f),
                                             t2,0.51498869f),
                                             t2,0.87890594f),
                                             t2,0.5f);
                                           
          }  
          else {
                t0  = expf(x)/sqrtf(x);
                t   = 3.75f/x;
                bi0 = (((((((( 
                      0.00392377f    * t 
                    - 0.01647633f)   * t 
                    + 0.02635537f)   * t 
                    - 0.02057706f)   * t 
                    + 0.916281e-02f) * t 
                    - 0.157565e-02f) * t 
                    + 0.225319e-02f) * t 
                    + 0.01328592f)   * t 
                    + 0.39894228f)   * t0;
               bi1  = (((((((( 
                    - 0.420059e-02f * t 
                    + 0.01787654f ) * t 
                    - 0.02895312f ) * t 
                    + 0.02282967f ) * t 
                    - 0.01031555f ) * t 
                    + 0.163801e-02f) * t 
                    - 0.00362018f ) * t 
                    - 0.03988024f ) * t 
                    + 0.39894228f ) * t0;
          }   
          
           if(x<=2.0f) {
               t   = x*0.5f;
               t2  = t*t;
               t0  = logf(t);
               bk0 = 
                     __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(
                     __fmaf_ru(0.0000074f,t2,0.0001075f),
                                          t2,0.00262698f),
                                          t2,0.0348859f),
                                          t2,0.23069756f),
                                          t2,0.4227842f)*t2-
                                          0.57721566f-bi0*t0;
               bk1 = 
                     (((((( 
                   - 0.00004686f  * t2 
                   - 0.00110404f) * t2 
                   - 0.01919402f) * t2 
                   - 0.18156897f) * t2 
                   - 0.67278579f) * t2 
                   + 0.15443144f) * t2
                   + 1.0f)/x+bi1*t0;
           }    
           else {
                t0  = expf(x)/sqrtf(x);
                t2  = t*t;
                bk0 = (((((( 
                      0.00053208f  * t 
                    - 0.0025154f)  * t 
                    + 0.00587872f) * t 
                    - 0.01062446f) * t 
                    + 0.02189568f) * t 
                    - 0.07832358f) * t 
                    + 1.25331414f) * t0;
                bk1 = (((((( 
                    - 0.00068245f   * t 
                    + 0.00325614f ) * t 
                    - 0.00780353f)  * t 
                    + 0.01504268f)  * t 
                    - 0.0365562f)   * t  
                    + 0.23498619f)  * t 
                    + 1.25331414f)  * t0;
           }   
           di0 = bi1;
           di1 = bi0-bi1/x;
           dk0 = -bk1;
           dk1 = -bk0-bk1/x;            
     }
     
     
/*

     !*****************************************************************************80
!
!! IKNA compute Bessel function In(x) and Kn(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order of In(x) and Kn(x).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  BI(0:N), DI(0:N), BK(0:N), DK(0:N),
!    the values of In(x), In'(x), Kn(x), Kn'(x).
!

*/


        __device__ void ikna(const int32_t n,
                             const float   x,
                             int32_t     &nm,
                             float * __restrict__ bi,
                             float * __restrict__ di,
                             float * __restrict__ bk,
                             float * __restrict__ dk) {
                             
             float bi0,bi1,bk0,bk1;
             float di1,di0,dk1,dk0;
             float f,f0,f1;
             float g,g0,g1;
             float h,h0,h1;
             float s0,invx;
             int32_t k,m,nm;
             
             if(n<=1) return;
             nm = n;
             ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1);
             bi[0] = bi0;
             bi[1] = bi1;
             bk[0] = bk0;
             bk[1] = bk1;
             di[0] = di0;
             di[1] = di1;
             dk[0] = dk0;
             dk[1] = dk1;  
             if(40.0f<x && n<(int32_t)(0.25f*x)) {
                 invx = 1.0f/x;
                 h0   = bi0;
                 h1   = bi1;
                 #pragma unroll
                 for(k=2; k<=n; ++k) {
                     float tk = (float)k;
                     h        = __fmaf_ru(-2.0f*(tjk-1.0f)*invx,h1,h0);
                     bi[k]    = h;
                     h0       = h1;
                     h1       = h;
                 }
             } 
             else {
                  invx = 1.0f/x;
                  m = msta1(x,200);
                  if(m<n)
                     nm = n;
                  else
                     m =  msta2(x,n,15);
                  f0   =  0.0f;
                  f1   =  1.1754943508e-38f;
                  for(k=m; k>=0; --k) {
                      float tk = (float)k;
                      f        = __fmaf_ru(2.0f*(tk+1.0f)*f1,invx,f0);
                      if(k<=nm) bi[k] = f;
                      f0       = f1;
                      f1       = f; 
                  }
                  s0 = bi0/f;
                  for(k=0; k<=nm; ++k) bi[k] *= s0;
                  g0 = bk0;
                  g1 = bk1;
                  #pragma unroll
                  for(k=2; k<=nm; ++k) {
                      float tk = (float)k;
                      g        = __fmaf_ru(2.0f*(tk-1.0f)*invx,g1,g0);
                      bk[k]    = g;
                      g0       = g1;
                      g1       = g;
                  }
                  #pragma unroll
                  for(k=2; k<=nm; ++k) {
                      float tk =  (float)k;
                      float t0 =  tk*invx;
                      di[k]    =  bi[k-1]-t0*bi[k];
                      dk[k]    = -bk[k-1]-t0*bk[k];
                  }
             }             
      }
      
      
/*

      !*****************************************************************************80
!
!! FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: KS, the sign code.
!    0, to calculate F+(x) and K+(x);
!    1, to calculate F_(x) and K_(x).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  FR, FI, FM, FA, the values of
!    Re[F+/-(x)], Im[F+/-(x)], |F+/-(x)|, Arg[F+/-(x)]  (Degs.).
!
!    Output, real(kind=sp) ::  GR, GI, GM, GA, the values of
!    Re[K+/-(x)], Im[K+/-(x)], |K+/-(x)|, Arg[K+/-(x)]  (Degs.).
!        

*/


        __device__ void ffk(const int32_t ks,
                            const float   x,
                            float         &fr,
                            float         &fi,
                            float         &fm,
                            float         &fa,
                            float         &gr,
                            float         &gi,
                            float         &gm,
                            float         &ga) {
                            
              float c1,cs,fi0;
              float s1,ss,invx;
              float x2,x4,xa,xc;
              float xf,xf0,xf1,xg;
              float xp,xq,xq2,xr;
              float sx2,cx2;
              float xs,xsu,xw;
              int32_t k,m;
              
              constexpr float srd = 57.29577951308233f;
              constexpr float eps = 1.0e-15f;
              constexpr float pi  = 3.141592653589793f;
              constexpr float pp2 = 1.2533141373155f;
              constexpr float p2p = 0.7978845608028654f;
              xa = fabsf(x);
              x2 = x*x;
              x4 = x2*x2;

              if(x==0.0f) {
                 fr = 0.62665706865775012560394f;
                 fi = powf(-1.0f,ks*fr);
                 fm = 1,25331413731550025120788f;
                 fa = powf(-1.0,ks*45.0f);
                 gr = 0.5f;
                 gi = 0.0f;
                 gm = 0.5f;
                 ga = 0.0f;
             }
             else {
                 if(xa<2.5f) {
                    xr = p2p*xa;
                    c1 = xr;
                    for(k=1; k<=50; ++k) {
                        float tk = (float)k;
                        float t0 = 4.0f*tk-3.0f;
                        float t1 = 2.0f*tk-1.0f;
                        float t2 = __fmaf_ru(4.0f,tk,1.0f);
                        xr       = -0.5f*xr*t0/tk/t1/t2*x4;
                        c1       += xr;
                        if(fabsf(xr/c1)<eps) break
                    }
                     s1 = p2p*xa*xa*xa*0.3333333333333333333333f;
                     xr = s1;
                     for(k=1; k<=50; ++k) {
                         float tk = (float)k;
                         float t0 = 4.0f*tk-1.0f;
                         float t1 = __fmaf_ru(2.0f,tk,1.0f);
                         float t2 = __fmaf_ru(4.0f,tk,3.0f);
                         xr       = -0.5f*xr*t0/tk/t1/t2*x4;
                         s1       += xr;
                         if(fabsf(xr/s1)<eps) break;
                     }
                 }
                 else if(x<5.5f) {
                     m   = (int32_t)(__fmaf_ru(x2,1.75f,42.0f);
                     invx= 1.0f/x2;
                     xsu = 0.0f;
                     xc  = 0.0f;
                     xs  = 0.0f;
                     xf1 = 0.0f;
                     xf0 = 1.0e-15f;
                     for(k=m; k>=0; --k) {
                         float tk = (float)k;
                         xf       = __fmaf_ru(2.0f,tk,3.0f)*xf0*invx-xf1;
                         if(k==2*(int32_t)(k/2))
                            xc    += xf;
                         else
                            xc    += xs;
                         xsu      = __fmaf_ru(xf,xf*__fmaf_ru(2.0f,tk,1.0f),xsu);
                         xf1      = xf0;
                         xf0      = xf;
                     }
                       xq = sqrtf(xsu);
                       xw = p2p*xa/xq;
                       c1 = xc*xw;
                       s1 = xs*xw;
                     }
                     else {
                        xr   = 1.0f;
                        xf   = 1.0f;
                        invx = 1.0f/x4;
                        #pragma unroll
                        for(k=1; k<=12; ++k) {
                            float tk = (float)k;
                            float t1 = 4.0f*tk-1.0f;
                            float t2 = t1-2.0f;
                            xr       = -0.25f*xr*t1*t2*invx;
                            xf       += xr;
                        }
                        xr = 1.0f/(2.0f*xa*xa);
                        xg = xr;
                        cx2= cosf(x2);
                        for(k=1; k<=12; ++k) {
                            float tk = (float)k;
                            float t1 = __fmaf_ru(4.0f,tk,1.0f);
                            float t2 = 4.0f*tk-1.0f;
                            xr       = -0.25f*xr*t1*t2*invx;
                            xg       += xr;
                        }
                        sx2 = sinf(x2);
                        c1  = 0.5f+xf*sx2-xg*cx2*0.39894228040143267793995f/(xa);
                        s1  = 0.5f-__fmaf_ru(xf,cx2,xg*sx2)*0.39894228040143267793995f/(xa);
                     }
                      fr  = pp2*(0.5f-c1);
                      fi0 = pp2*(0.5f-s1);
                      fi  = powf(-1.0f,ks)*fi0;
                      fm  = sqrtf(__fmaf_ru(fr,fr,fi*fi));
                      if(0.0<=fr) 
                         fa = srd*atanf(fi/fr);
                      else if(0.0f<fi)
                         fa = srd*(atanf(fi/fr)+pi);
                      else if(fi<0.0f) 
                         fa = srd*(atan(fi/fr)-pi);
                      xp  = __fmaf_ru(x,x,0.78539816339744830961566f);
                      cs  = cosf(xp);
                      ss  = sinf(xp);
                      xq2 = 0.56418958354775628694808f;
                      gr  = xq2*__fmaf_ru(fr,cs,fi0 * ss);
                      gi  = powf(-1.0f,ks)*xq2*(fi0*cs-fr*ss);
                      gm  = sqrtf( _fmaf_ru(gr,gr,gi*gi);
                      if(0.0f<=gr)
                         ga = srd*atanf(gi/gr);
                      else if(0.0f<gi)
                         ga = srd*(atanf(gi/gr)+pi);
                      else if(gi<0.0f) 
                         ga = srd*(atanf(gi/gr)-pi);
   
                      if(x<0.0f) {
                         fr = pp2-fr;
                         fi = powf(-1.0f,ks)*pp2-fi;
                         fm = sqrtf(__fmaf_ru(fr,fr,fi*fi));
                         fa = srd*atanf(fi/fr);
                         gr = cosf(x*x)-gr;
                         gi = -powf(-1.0f,ks)*sinf(x*x)-gi;
                         gm = sqrtf( __fmaf_ru(gr,gr,gi * gi));
                         ga = srd*atanf(gi/gr);
                      }
                     
               }  
                       
     } 
                            
       
       
/*
  
     !*****************************************************************************80
!
!! CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, 
!    CDY1, the values of J0(z), J0'(z), J1(z), J1'(z), Y0(z), Y0'(z), 
!    Y1(z), Y1'(z).

*/


        __device__ void cjy01(const cuda::std::complex<float> z,
                              cuda::std::complex<float> &     cbj0,
                              cuda::std::complex<float> &     cdj0,
                              cuda::std::complex<float> &     cbj1,
                              cuda::std::complex<float> &     cdj1,
                              cuda::std::complex<float> &     cby0,
                              cuda::std::complex<float> &     cdy0,
                              cuda::std::complex<float> &     cby1,
                              cuda::std::complex<float> &     cdy1) {
                              
            
              float a[13] = {0.0f,
                -0.703125e-01f,0.112152099609375e+00f, 
                -0.5725014209747314e+00f,0.6074042001273483e+01f, 
                -0.1100171402692467e+03f,0.3038090510922384e+04f, 
                -0.1188384262567832e+06f,0.6252951493434797e+07f, 
                -0.4259392165047669e+09f,0.3646840080706556e+11f, 
                -0.3833534661393944e+13f,0.4854014686852901e+15f}; 
              float a1[13] = {0.0f,
                 0.1171875e+00f,-0.144195556640625e+00f, 
                 0.6765925884246826e+00f,-0.6883914268109947e+01f, 
                 0.1215978918765359e+03f,-0.3302272294480852e+04f, 
                 0.1276412726461746e+06f,-0.6656367718817688e+07f, 
                 0.4502786003050393e+09f,-0.3833857520742790e+11f, 
                 0.4011838599133198e+13f,-0.5060568503314727e+15f};
              float b[13] = {0.0f,
                 0.732421875e-01f,-0.2271080017089844e+00f, 
                 0.1727727502584457e+01f,-0.2438052969955606e+02f, 
                 0.5513358961220206e+03f,-0.1825775547429318e+05f, 
                 0.8328593040162893e+06f,-0.5006958953198893e+08f, 
                 0.3836255180230433e+10f,-0.3649010818849833e+12f, 
                 0.4218971570284096e+14f,-0.5827244631566907e+16f};
              float b1[13] = {0.0f,
                 -0.1025390625e+00f,0.2775764465332031e+00f, 
                 -0.1993531733751297e+01f,0.2724882731126854e+02f, 
                 -0.6038440767050702e+03f,0.1971837591223663e+05f, 
                 -0.8902978767070678e+06f,0.5310411010968522e+08f, 
                 -0.4043620325107754e+10f,0.3827011346598605e+12f, 
                 -0.4406481417852278e+14f,0.6065091351222699e+16};
              cuda::std::complex<float> ci,cp,cp0,cp1;
              cuda::std::complex<float> cq0,cq1,cr,cs;
              cuda::std::complex<float> ct1,ct2,cu,z1,z2;
              float a0,w0,w1;
              int   k,k0;
              
              constexpr float pi = 3.141592653589793f;
              constexpr float el = 0.5772156649015329f;
              constexpr float rp2= 0.63661977236758134307554f;
              ci                 = {0.0f,1.0f};
              a0                 = abs(z);
              z2                 = z*z;
              z1                 = z;
              
              if(real(z)<0.0f) z1 = -z;
              if(a0<=12.0f) {
                 
                 cbj0 = {1.0f,0.0f};
                 cr   = {1.0f,0.0f};
                 for(k=1; k!=40; ++k) {
                     float tk = (float)k;
                     cr       = -0.25f*cr*z2/(tk*tk);
                     cbj0     += cr;
                     if(abs(cr)<abs(cbj0)*1.0e-15f) break;
                 }
                 
                 cbj1 = {1.0f,0.0f};
                 cr   = {1.0f,0.0f};
                 for(k=1; k!=40; ++k) {
                     float tk = (float)k;
                     cr       = -0.25f*cr*z2/(tk*(tk+1.0f));
                     cbj1     += cr;
                     if(abs(cr)<abscbj1)*1.0e-15f) break;
                 }
                 
                 cbj1 = 0.5f*z1*cbj1;
                 w0   = 0.0f;
                 cr   = {1.0f,0.0f};
                 cs   = {0.0f,0.0f};
                 for(k=1; k!=40; ++k) {
                     float tk = (float)k;
                     w0       = w0+1.0f/tk;
                     cr       = -0.25f*cr/(tk*tk)*z2;
                     cp       = cr*w0;
                     cs       = cs+cp;
                     if(abs(cp)<abs(cs)*1.0e-15f) break;
                 }
                 
                 cby0 = __fmaf_ru(rp2,(log(z1/2.0f),el)*cbj0-rp2*cs;
                 w1   = 0.0f;
                 cr   = {1.0f,0.0f};
                 cs   = {1.0f,0.0f};
                 for(k=1; k!=40; ++k) {
                     float tk = (float)k;
                     w1       = w1+1.0f/tk;
                     cr       = -0.25f*cr/(tk*(tk+1.0f))*z2;
                     cp       = cr*(2.0f*w1+1.0f/(tk+1.0f));
                     cs       = cs+cp;
                     if(abs(cp)<abs(cs)*1.0e-15f) break; 
                 }
                 
                 cby1 = rp2*((log(z1*0.5f)+el)*cbj1 -
                        1.0f/z1-0.25f*z1*cs);

              }
              else {
                 
                  if(a0<35.0f) 
                     k0 = 12;
                  else if(a0<50.0f) 
                     k0 = 10;
                  else
                     k0 = 8;
                 ct1 = z1-0.25f*pi;
                 cp0 = {1.0f,0.0f};
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     cp0      = __fmaf_ru(cp0,a[k],pow(z1,-2.0f*tk));
                 }
                 cq0 = -0.125f/z1;
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     cq0      = __fmaf_ru(cq0,b[k],pow(z1,-2.0f*t-1.0f));
                 }
                 cu   = sqrt(rp2/z1);
                 cbj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
                 cby0 = cu*(__fmaf_ru(cp0,sin(ct1),cq0*cos(ct1)));
                 ct2  = z1-0.75f*pi;
                 cp1  = {1.0f, 0.0f};
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     cp1      = __fmaf_ru(cp1,a1[k],pow(z1,-2.0f*tk));
                 }
                 cq1 = 0.375f/z1;
                 for(k=1; k!=k0; ++k) {
                     float tk = (float)k;
                     cq1      = __fmaf_ru(cq1,b1[k],pow(z1,-2.0f*tk-1.0f));
                 }
                 cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
                 cby1 = cu*(__fmaf_ru(cp1,sin(ct2),cq1*cos(ct2)));
              }
              
             if(real(z)<0.0f){
                if(imag(z)<0.0f) {
                   cby0 = cby0-2.0f*ci*cbj0;
                   cby1 = -(cby1-2.0f*ci*cbj1);
                }
                else {
                   cby0 = cby0+2.0f*ci*cbj0;
                   cby1 = -(cby1+2.0f*ci*cbj1);
                }
                cbj1 = -cbj1;
             }

             cdj0 = -cbj1;
             cdj1 = cbj0-1.0f/z*cbj1;
             cdy0 = -cby1;
             cdy1 = cby0-1.0f/z*cby1;
      }
       
       
/*
   
    !*****************************************************************************80
!
!! CJYLV: Bessel functions Jv(z), Yv(z) of complex argument and large order v.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V, the order of Jv(z) and Yv(z).
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBJV, CDJV, CBYV, CDYV, the values of Jv(z), 
!    Jv'(z), Yv(z), Yv'(z).
!

*/


        __device__ void cjylv(const float v,
                              const cuda::std::complex<float> z,
                              cuda::std::complex<float> &     cbjv,
                              cuda::std::complex<float> &     cdjv,
                              cuda::std::complex<float> &     cbyv,
                              cuda::std::complex<float> &     cdyv) {
                              
              float a[91];
              cuda::std::complex<float> cf[13];
              int p1k[12] = {-1,1,-1,1,-1,1,-1,1,-1,1,-1,1};
              cuda::std::complex<float> ceta,cfj,cfy,csj;
              cuda::std::complex<float> csy,ct,ct2,cws;
              float                     v0,vr;
              int                       i,k,km,l;
              int                       l0,lf;
              constexpr float pi = 3.14159265358979323846264f;
              km                 = 12;
              cjk(km,a);
              
              for(l=1; l!=0; --l) {
                  
                  v0   = v-l;
                  cws  = sqrt(1.0f-(z/v0)*(z/v0));
                  ceta = cws+log(z/v0/(1.0f+cws));
                  ct   = 1.0f/cws;
                  ct2  = ct*ct;
                  for(k=1; k!=km; ++k) {
                      l0    = k*(k+1)/2+1;
                      lf    = l0+k;
                      cf[k] = a[lf];
                      for(i=lf-1; i!=l0; --i)  
                          cf[k] = cf[k]*ct2+a[i];
                      cf[k] = cf[k]*pow(ct,(float)k);   
                  }
                  
                  vr = 1.0f/v0;
                  csj= {1.0f,0.0f};
                  for(k=1; k!=km; ++k) 
                      csj = csj+cf[k]*pow(vr,(float)k);
                  cbjv = sqrt(ct/(2.0f*pi*v0))*exp(v0*ceta)*csj;
                  if(l==1) cfj = cbjv;
                  csy = {1.0f,0.0f};
                  for(k=1; k!=km; ++k)
                      csy = csy+p1k[k]*cf[k]*pow(vr,(float)k);
                  cbyv = -sqrt(2.0f*ct/(pi*v0))*exp(-v0*ceta)*csy;
                  if(l==1) cfy = cbyv;
              }
              
              cdjv = -v/z*cbjv+cfj;
              cdyv = -v/z*cbyv+cfy;
                         
      }
      
      
/*
    !*****************************************************************************80
!
!! CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of Jn(z) and Yn(z).
!
!    Input, complex ( kind = 8 ) Z, the argument of Jn(z) and Yn(z).
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, complex ( kind = 8 ), CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
!    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
!
*/
    
    
        __device__ void cjyna(const int n,
                              const cuda::std::complex<float> z,
                              int                           & nm,
                              cuda::std::complex<float> * __restrict__ cbj,
                              cuda::std::complex<float> * __restrict__ cdj,
                              cuda::std::complex<float> * __restrict__ cby,
                              cuda::std::complex<float> * __restrict__ cdy) {
                              
               
              cuda::std::complex<float> cbj0,cbj1,cby0,cby1;
              cuda::std::complex<float> cdj0,cdj1,cdy0,cdy1;
              cuda::std::complex<float> cf,cf1,cf2,cg0;
              cuda::std::complex<float> cg1,ch0,ch1,ch2;
              cuda::std::complex<float> cj0,cj1,cjk,cp11;
              cuda::std::complex<float> cp12,cp21,cp22;
              cuda::std::complex<float> cs,cyk,cyl1,cyl2;
              cuda::std::complex<float> cylk;
              float                     a0,wa,ya0,ya1,yak;
              int                       k,lb,lb0;
              
              if(n<1) return;
              constexpr float pi = 3.14159265358979323846264;
              a0                 = abs(z);
              nm                 = n;
              cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1);  
              cbj[0] = cbj0;
              cbj[1] = cbj1;
              cby[0] = cby0;
              cby[1] = cby1;
              cdj[0] = cdj0;
              cdj[1] = cdj1;
              cdy[0] = cdy0;
              cdy[1] = cdy1;
               
              if(n<(int)(0.25f*a0)) {
                 cj0 = cbj0;
                 cj1 = cbj1;
                 for(k=2; k!=n; ++k) {
                     float tk = (float)k;
                     cjk      = 2.0f*(tk-1.0f)/z*cj1-cj0;
                     cbj[k]   = cjk;
                     cj0      = cj1;
                     cj1      = cjk;
                 }
              }  
              else {
                 
                 m = msta1(a0,200);
                 if(m<n)
                     nm = m;
                 else
                     m  = msta2(a0,n,15);
                 cf2 = {0.0f,0.0f}; 
                 cf1 = {(float)1.0e-30f, 0.0f};
                 for(k=m; k!=0; --k) {
                     float tk = (float)k;
                     cf = __fmaf_ru(2.0f,tk,1.0f)/z*cf1-cf2;
                     if(k<=nm) 
                        cbj[k] = cf;
                     cf2 = cf1;
                     cf1 = cf;
                 }
                 
                 if(abs(cbj1)<abs(cbj0))
                    cs = cbj0/cf;
                 else
                    cs = cbj1/cf2;
                 
                 for(k=0; k!=nm; ++k) cbj[k] *= cs;
             } 
             
             for(k=2; k!=nm; ++k) {
                 float tk = (float)k;
                 cdj[k] = cbj[k-1]-k/z*cbj[k];
             }
             ya0 = abs(cby0);
             lb = 0;
             cg0 = cby0;
             cg1 = cby1;
             for(k=2; k!=nm; ++k) {
                 float tk = (float)k;
                 cyk      = 2.0f*(tk-1.0f)/z*cg1-cg0;
                 if(abs(cyk)<=3.4028234664e+38f) {
                    yak = abs(cyk);
                    ya1 = abs(cg0);
                    if(yak<ya0 && yak<ya1) lb = k;
                    cby[k] = cyk;
                    cg0 = cg1;
                    cg1 = cyk;
                }
            }
             
           if(4<lb && imag(z)!=0.0f) {
              
              while(true) {
                 
                 if(lb==lb0) exit
                 
                 ch2 = {1.0f,0.0f};
                 ch1 = {0.0f,0.0f};
                 lb0 = lb;
                 for(k=lb; k!=1; ++k) {
                     float tk = (float)k;
                     ch0      = 2.0f*tk/z*ch1-ch2;
                     ch2      = ch1;
                     ch1      = ch0;
                 }
                 cp12 = ch0;
                 cp22 = ch2;
                 ch2  = {0.0f,0.0f};
                 ch1  = {1.0f,0.0f};
                 for(k=lb; k!=1; ++k) {
                     float tk = (float)k;
                     ch0      = 2.0f*tk/z*ch1-ch2;
                     ch2      = ch1;
                     ch1      = ch0;
                 }
                 cp11 = ch0;
                 cp21 = ch2;
                 if(lb==nm) {
                    cbj[lb+1] = 2.0f*lb/z*cbj[lb]-cbj[lb-1];
                 }
                 if(abs(cbj[1])<abs(cbj[0])) {
                    cby[lb+1] = (cbj[lb+1]*cby0-2.0f*cp11/(pi*z))/cbj[0];
                    cby[lb]   = (__fmaf_ru(cbj[lb],cby0,2.0f)*cp12/(pi*z))/cbj[0[;
                 }
                 else {
                    cby[lb+1] = (cbj[lb+1]*cby1-2.0f*cp21/(pi*z))/cbj[1];
                    cby[lb]   = (__fmaf_ru(cbj[lb],cby1,2.0f)*cp22/(pi*z))/cbj[1];
                 }
                cyl2 = cby[lb+1];
                cyl1 = cby[lb];
                for(k=lb-1; k!=0; --k) {
                    float tk = (float)k;
                    cylk     = __fmaf_ru(2.0f,tk,1.0f)/z*cyl1-cyl2;
                    cby[k]   = cylk;
                    cyl2     = cyl1;
                    cyl1     = cylk;
                }
                cyl1 = cby[lb];
                cyl2 = cby[lb+1];
                for(k=lb+1; k!=nm-1; ++k) {
                    float tk = (float)k;
                    cylk     = 2.0f*tk/z*cyl2-cyl1;
                    cby[k+1] = cylk;
                    cyl1     = cyl2;
                    cyl2     = cylk;
                }
                for(k=2; k!=nm; ++k) {
                    wa = abs(cby[k]);
                    if(wa<abs(cby[k-1])) lb = k;
                }
              }
           }
           
           for(k=2; k!=nm; ++k) {
             float tk = (float)k;
             cdy[k]   = cby[k-1]-tk/z*cby[k];
           }
     }
     
     
/*

   !*****************************************************************************80
!
!! CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer KF, the function choice.
!    0 for ZF = Y0(z) and ZD = Y0'(z);
!    1 for ZF = Y1(z) and ZD = Y1'(z);
!    2 for ZF = Y1'(z) and ZD = Y1''(z).
!
!    Input, complex(kind=sp) ::  Z, the argument.
!
!    Output, complex(kind=sp) ::  ZF, ZD, the values of the requested function 
!    and derivative.   


*/

        __device__ void cy01(const int kf,
                             const cuda::std::complex<float> z,
                             cuda::std::complex<float>      &zf,
                             cuda::std::complex<float>      &zd) {
                             
                const float a[12] = {
                     -0.703125e-01f,           0.112152099609375f,
                     -0.5725014209747314f,     0.6074042001273483e+01f,
                     -0.1100171402692467e+03f, 0.3038090510922384e+04f, 
                     -0.1188384262567832e+06f, 0.6252951493434797e+07f, 
                     -0.4259392165047669e+09f, 0.3646840080706556e+11f, 
                     -0.3833534661393944e+13f, 0.4854014686852901e+15f};
                const float a1[12] = {
                      0.1171875f, -0.144195556640625f,
                      0.6765925884246826f,       -0.6883914268109947e+01f, 
                      0.1215978918765359e+03f, -0.3302272294480852e+04f, 
                      0.1276412726461746e+06f, -0.6656367718817688e+07f, 
                      0.4502786003050393e+09f, -0.3833857520742790e+11f, 
                      0.4011838599133198e+13f, -0.5060568503314727e+15f};     
                const float b[12]  = {
                      0.732421875e-01f, -0.2271080017089844f,
                      0.1727727502584457e+01f, -0.2438052969955606e+02f, 
                      0.5513358961220206e+03f, -0.1825775547429318e+05f, 
                      0.8328593040162893e+06f, -0.5006958953198893e+08f, 
                      0.3836255180230433e+10f, -0.3649010818849833e+12f, 
                      0.4218971570284096e+14f, -0.5827244631566907e+16f};
                const float b1[12]  = {
                     -0.1025390625f,           0.2775764465332031f,
                     -0.1993531733751297e+01f, 0.2724882731126854e+02f,
                     -0.6038440767050702e+03f, 0.1971837591223663e+05f,
                     -0.8902978767070678e+06f, 0.5310411010968522e+08f,
                     -0.4043620325107754e+10f, 0.3827011346598605e+12f,
                     -0.4406481417852278e+14f, 0.6065091351222699e+16f};
                cuda::std::complex<float>  cbj0,cbj1;
                cuda::std::complex<float>  cby0,cby1;
                cuda::std::complex<float>  cdy0,cdy1;
                cuda::std::complex<float>  ci,cp;
                cuda::std::complex<float>  cp0,cp1;
                cuda::std::complex<float>  cq0,cq1;
                cuda::std::complex<float>  cr,cs;
                cuda::std::complex<float>  ct1,ct2;
                cuda::std::complex<float>  cu,z1,z2;
                cuda::std::complex<float>  ctx,stx;
                float                      w0,w1,a0;
                int32_t                    k,k0,idx;
                
                constexpr float pi =  3.14159265358979323846264f;
                constexpr float el =  0.5772156649015329f;
                constexpr float rp2=  0.63661977236758134307554f;
                ci                 =  {0.0f,1.0f};
                a0                 =  abs(z);
                z2                 =  z*z;
                z1                 =  z;
                
                if(a0==0.0f){ 
                   return;
                }
                else {
                   if(real(z)<0.0f) z1 = -z; 
                   if(a0<=12.0f) {
                      cbj0 = {1.0f,0.0f};
                      cr   = {1.0f,0.0};
                      for(k=1; k<=40; ++k) {
                          float tk = (float)k;
                          float t0 = tk*tk;
                          cr       = -0.25f*cr*z2/t0;
                          cbj0     += cr;
                          if(abs(cr)<abs(cbj0)*1.0e-15f) break;
                      }
                      cbj1 = {1.0f,0.0f};
                      cr   = {1.0f,0.0f};
                      for(k=1; k<=40; ++k) {
                          float tk = (float)k;
                          float t0 = tk*(tk+1.0f);
                          cr       = -0.25f*cr*z2/t0;
                          cbj1     += cr;
                          if(abs(cr)<abs(cbj1)*1.0e-15f) break;
                      }
                      cbj1 = 0.5f*z1*cbj1;
                      w0   = 0.0f;
                      cr   = {1.0f,0.0f};
                      cs   = {0.0f,0.0f};
                      for(k=1; k<=40; ++k) {
                          float tk = (float)k;
                          float t0 = tk*tk;
                          cr       = -0.25f*cr/t0*z2;
                          cp       = cr*w0;
                          cs       += cp;
                          if(abs(cp)<abs(cs)*1.0e-15f) break;
                      }
                       cby0 = rp2*(log(z1*0.5f)+el)*cbj0-rp2*cs;
                       w1   = 0.0f;
                       cr   = {1.0f,0.0f};
                       cs   = {1.0f,0.0f};
                       for(k=1; k<=40; ++k) {
                           float tk = (float)k;
                           float t0 = tk*(tk+1.0f);
                           float t1 = tk+1.0f;
                           w1       = w1+1.0f/tk;
                           cr       = -0.25f*cr/t0*z2;
                           cp       = cr*(2.0f*w1+1.0f/t1);
                           cs       = cs+cp;
                           if(abs(cp)<abs(cs)*1.0e-15f) break;
                       }
                       cby1 = rp2*((logf(z1*0.5f)+el)*cbj1-
                              1.0f/z1-0.25f*z1*cs);
                   }
                   else {
                       if(a0<35.0f)
                          k0 = 12;
                       else if(a0<50.0f)
                          k0 = 10;
                       else
                          k0 = 8;
                        ct1 = z1-0.78539816339744830961566f;
                        ctx = cos(ct1);
                        cp0 = {1.0f,0.0f};
                        cq0 = -0.125f/z1;
                        idx = -1;
                        for(k=1; k<=k0; ++k) {
                            ++idx;
                            cp0 = cp0+a[idx]*pow(z1,-2*k);
                            cq0 = cq0+b[idx]*pow(z1,-2*k-1);
                        }
                        stx  = sin(ct1);
                        cu   = sqrt(rp2/z1);
                        cbj0 = cu*cp0*ctx-cq0*stx;
                        cby0 = cu*cp0*stx+cq0*ctx;
                        ct2  = z1-2.35619449019234492884698f;
                        ctx  = cos(ct2);
                        cq1  = 0.375f/z1;
                        cp1  = {1.0f,0.0f};
                        idx  = -1;
                        for(k=1; k<=k0; ++k) {
                            ++idx;
                            cp1 = cp1+a1[idx]*pow(z1,-2*k);
                            cq1 = cq1+b1[idx]*pow(z1,-2*k-1);
                        }
                        stx  = sin(ct2);
                        cbj1 = cu*cp1*ctx-cq1*stx;
                        cby1 = cu*cp1*stx+cq1*ctx;
                   }
                   if(real(z)<0.0f) {
                      if(imag(z)<0.0f) {
                          cby0 = cby0-2.0f*ci*cbj0;
                          cby1 = -(cby1-2.0f*ci*cbj1);
                      }
                      else {
                          cby0 = cby0+2.0f*ci*cbj0;
                          cby1 = -(cby1+2.0f*ci*cbj1);
                      }
                       cbj1 = -cbj1;
                   }
                    cdy0 = -cby1;
                    cdy1 = cby0-1.0f/z*cby1;
                   
                }      
                 
                 if(kf==0) {
                    zf = cby0;
                    zd = cdy0;
                 }
                 else if(kf==1) {
                    zf = cby1;
                    zd = cdy1;
                 }
                 else if(kf==2) {
                    zf = cdy1
                    zd = -cdy1/z-(1.0f-1.0f/(z*z))*cby1;
                 }
       }  
       
       
/*

    !*****************************************************************************80
!
!! EIX computes the exponential integral Ei(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  EI, the function value.
!  

*/


        __device__ float eix(const float x) {
        
             float ga,r,invx,ei;
             int32_t k;
             
             if(x==0.0f) {
                ei = -3.4028234664e+38f;
             } 
             else if(x<=40.0f) {
                ei = 1.0f;
                r  = 1.0f;
                for(k=1; k<=100; ++k) {
                    float tk = (float)k;
                    float t0 = (tk+1.0f);
                    r        = r*tk*x/(t0*t0);
                    ei       +=r;
                    if(fabsf(r/ei)<=1.0e-15f) break;
                }
                ga = 0.5772156649015328f;
                ei = ga+logf(x)+x*ei;
             }
             else {
                invx = 1.0f/x;
                ei   = 1.0f;
                r    = 1.0f;
                r    = r*invx;
                ei   +=r;
                r    = r*2.0f*invx;
                ei   +=r;
                r    = r*3.0f*invx;
                ei   +=r;
                r    = r*4.0f*invx;
                ei   +=r;
                r    = r*5.0f*invx;
                ei   +=r;
                r    = r*6.0f*invx;
                ei   +=r;
                r    = r*7.0f*invx;
                ei   +=r;
                r    = r*8.0f*invx;
                ei   +=r;
                r    = r*9.0f*invx;
                ei   +=r;
                r    = r*10.0f*invx;
                ei   +=r;
                r    = r*11.0f*invx;
                ei   +=r;
                r    = r*12.0f*invx;
                ei   +=r;
                r    = r*13.0f*invx;
                ei   +=r;
                r    = r*14.0f*invx;
                ei   +=r;
                r    = r*15.0f*invx;
                ei   +=r;
                r    = r*16.0f*invx;
                ei   +=r;
                r    = r*17.0f*invx;
                ei   +=r;
                r    = r*18.0f*invx;
                ei   +=r;
                r    = r*19.0f*invx;
                ei   +=r;
                r    = r*20.0f*invx;
                ei   +=r;
                ei   = expf(x)*invx*ei;
             }
             return (ei);
        }
     
     
/*
    
     !*****************************************************************************80
!
!! COMELP computes complete elliptic integrals K(k) and E(k).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HK, the modulus.  0 <= HK <= 1.
!
!    Output, real ( kind = 8 ) CK, CE, the values of K(HK) and E(HK).
!
*/

       __device__ void comelp(const float hk,
                              float     & ck,
                              float     & ce) {
                              
             float ae,ak,be,bk;
             float ce,ck,pk;
             if(hk==1.0f) return;
             pk = 1.0f-hk*hk;  
             ak = 
                  __fmaf_ru(__fmaf_ru(__fmaf_ru(__fmaf_ru(0.01451196212e+00f, pk, 
                                                 ,0.03742563713e+00f ),pk, 
                                          0.03590092383e+00f ), pk, 
                                0.09666344259e+00f ),pk, 
                             1.38629436112e+00f);

            bk = 
                  __fmaf_ru(__fmaf_ru(__fmaf_ru(__fmaf_ru(0.00441787012e+00f,pk,
                                                 0.03328355346e+00f ),pk, 
                                          0.06880248576e+00f ),pk, 
                                0.12498593597e+00f ),pk, 
                              0.5e+00);

           ck = ak-bk*log(pk);

           ae = 
                __fmaf_ru(__fmaf_ru(__fmaf_ru(__fmaf_ru(0.01736506451e+00f,pk, 
                                              0.04757383546e+00f ),pk, 
                                        0.0626060122e+00f  ),pk, 
                               0.44325141463e+00f ),pk,
                         1.0e+00f);

           be = 
              __fmaf_ru(__fmaf_ru(__fmaf_ru(__fmaf_ru(0.00526449639e+00f,pk, 
                                             0.04069697526e+00f ),pk, 
                                     0.09200180037e+00f ),pk, 
                            0.2499836831e+00f  ),pk);

           ce = ae-be*log(pk);                    
     }
     
     
/*
   
     !****************************************************************************80
!
!! CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer N, the order.
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CDN, the function value.
!

*/

        __device__ cuda::std::complex<float> 
                                cpdla(const int n,
                                      const cuda::std::complex<float> cdn) {
            
               cuda::std::complex<float> cb0,cr,cdn;
               float                     fn;
               int                       k;
               
               fn  = (float)n;
               cb0 = pow(z,fn)*exp(-0.25f*z*z);
               cr  = {1.0f,0.0f};
               cdn = {1.0f,0.0f};
               #pragma unroll
               for(k=1; k!=16; ++k) {
                   float tk = (float)k;
                   float t0 = 2.0f*tk-fn;
                   cr       = -0.5f*cr*t0-1.0f*t0/(tk*z*z);
                   cdn      += cr;
                   if(abs(cr)<abs(cdn)*1.0e-12f) break;
               }   
               cdn = cb0*cdn;
               return (cdn);                         
       }
       
       
/*
     
*/
      
/*
     !*****************************************************************************80
!
!! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  SJ(0:N), the values of jn(x).
!
!    Output, real(kind=sp) ::  DJ(0:N), the values of jn'(x).
!
*/


        __device__ void sphj(const int n,
                             const float x,
                             int &       nm,
                             float * __restrict__ sj,
                             float * __restrict__ dj) {
                             
             
             float cs,f,f0,f1;
             float sa,sb,invx;
             int   k,m;
             
             invx  = 1.0f/x;
             nm    = n;
             sj[0] = sinf(x)*invx;
             sj[1] = (sj[0]-cosf(x))*invx;
             if(2<=n) {
                sa = sj[0];
                sb = sj[1];
                m  = msta1(x,200);
                if(m<n)
                   nm = n;
                else
                   m  = msta2(x,n,15);
                f0 = 0.0f;
                f1 = 1.0e-30f;
                for(k=m; k!=0; --k) {
                    float tk = (float)k;
                    f        = __fmaf_ru(2.0f,tk,3.0f)*f1*invx-f0;
                    if(k<=nm) sj[k] = f;
                    f0 = f1;
                    f1 = f;
                }
                
                if(fabsf(sa)<=fabsf(sb)
                   cs = sb/f0;
                else
                   cs = sa/f;
                for(k=0; k!=nm; ++k) sj[k] *= cs;
             }   
             
             dj[0] = (cosf(x)-sinf(x)*invx)*invx;
             for(k=1; k!=nm; ++k) {
                 float tk = (float)k;
                 dj[k]    = sj[k-1]-(tk+1.0f)*sj[k]*invx;
             }                     
     } 
     
/*

       !*****************************************************************************80
!
!! SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order of In(X).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  SI(0:N), DI(0:N), the values and derivatives
!    of the function of orders 0 through N.
! 
  
*/


        __device__ void sphi(const int   n,
                             const float x,
                             int &       nm,
                             float * __restrict__ si,
                             float * __restrict__ di) {
             
             float cs,f,f0,f1;
             float si0,invx;
             int   k,m,nm;
             
             invx   = 1.0f/x;
             nm     = n;
             si[0]  = sinhf(x)*invx;
             si[1]  = -(sinhf(x)*invx-coshf(x))*invx;
             si0    = si[0];
             if(2<=n) {
                
                m = msta1(x,200);
                if(m<n)
                   nm = m;
                else
                   m  = msta2(x,n,15);
                f0 = 0.0f;
                f1 = 1.0e-30f;
                for(k=m; k!=0; --k) {
                    float tk = (float)k;
                    f        = __fmaf_ru(__fmaf_ru(2.0f,tk,3.0f),(f1*invx),f0);
                    if(k<=nm) si[k] = f;
                    f0 = f1;
                    f1 = f;
                }
                cs = si0/f;
                for(k=0; k!=nm; ++k) si[k] *= cs;
             } 
             
             di[0] = si[1];
             for(k=1; k!=nm; ++k) {
                 float tk = (float)k;
                 di[k]    = si[k-1]-(tk+1.0f)*invx*si[k];
             }                     
      }
      
      
/*
    !*****************************************************************************80
!
!! SPHK computes modified spherical Bessel functions kn(x) and derivatives.
!
!  Discussion:
!
!    This procedure computes modified spherical Bessel functions
!    of the second kind, kn(x) and kn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  SK(0:N), DK(0:N), the values of kn(x) and kn'(x).
!
*/


        __device__ void sphk(const int   k,
                             const float x,
                             int &       nm
                             float * __restrict__ sk,
                             float * __restrict__ dk) {
             
              float f,f0,f1,invx;
              int   k;
              constexpr float pi = 3.14159265358979323846264f;
              invx               = 1.0f/x;
              nm                 = n;
              sk[0]              = 0.5f*pi*invx*expf(-x);
              sk[1]              = sk[0]*(1.0f+invx);
              f0                 = sk[0];
              f1                 = sk[1];
              #pragma unroll
              for(k=2; k!=n; ++k) {
                  float tk = (float)k;
                  f        = __fmaf_ru((2.0f*tk-1.0f),(f1*invx),f0);
                  sk[k]    = f;
                  if(3.4028234664e+38f<fabsf(f)) break;
                  f0 = f1;
                  f1 = f;
              }          
              nm    = k-1;
              dk[0] = -sk[1];
              #pragma unroll
              for(k=1; k!=nm; ++k) {
                  float tk = (float)k;
                  dk[k]    = -sk[k-1]-(tk+1.0f)*invx*sk[k]; 
              }       
       }
       
       
/*
   
   !*****************************************************************************80
!
!! SPHY computes spherical Bessel functions yn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
!   

*/

 
        __device__ void sphy(const int   n,
                             const float x,
                             int &       nm,
                             float * __restrict__ sy,
                             float * __restrict__ dy) {
                             
              float f,f0,f1,invx;
              int   k;
              
              nm    = n;
              invx  = 1.0f/x;
              sy[0] = -cosf(x)*invx;
              sy[1] = (sy[0]-sinf(x))*invx;
              f0    = sy[0];
              f1    = sy[1];
              #pragma unroll
              for(k=2; k!=n; ++k) {
                  float tk = (float)k;
                  f        = (2.0f*tk-1.0f)*f1*invx-f0;
                  sy[k]    = f;
                  if(3.4028234664e+38f <= fabsf(f)) break;
                  f0 = f1;
                  f1 = f;
              } 
              nm    = k-1;
              dy[0] = (sinf(x)+cosf(x)*invx)*invx;
              #pragma unroll
              for(k=1; k!=nm; ++k) {
                  float tk = (float)k;
                  dy[k]    = sy[k-1]-(tk+1.0f)*sy[k]*invx;
              }                     
       }
       
       
/*
   
     !*****************************************************************************80
!
!! STVH0 computes the Struve function H0(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  SH0, the value of H0(x).
! 

*/


        __device__ float stvh0(const float x) {
        
              float a0,by0,p0,q0;
              float r,s,t,t2;
              float ta0,invx,sh0;
              int   k,km;
              constexpr float pi  = 3.14159265358979323846264f;
              s                   = 1.0f;
              r                   = 1.0f;
              if(x<=20.0f) {
              
                 a0 = 2.0f*x/pi;
                 for(k=1; k!=60; ++k) {
                     float tk = (float)k;
                     float t0 = __fmaf_ru(2.0f,tk,1.0f);
                     r        = -r*x/t0*x/t0;
                     s        += r;
                     if(fabsf(r) < fabsf(s)*1.0e-12f) break;
                 }
                 sh0 = a0*s;
              }
              else {
              
                 invx = 1.0f/x;
                 if(x < 50.0f) 
                    km = (int)(0.5f*(x+1.0f));
                 else
                    km = 25;
                 
                 for(k=1; k!=km; ++k) {
                      float tk = (float)k;
                      float t0 = (2.0f*tk-1.0f)*invx;
                      float t1 = t0*t0;
                      r        = -r*t1;
                      s        += r;
                      if(fabsf(r) < fabsf(s)*1.0e-12f) break;
                  }
                  t  = 4.0f/x;
                  t2 = t*t;
                  p0 = (((( 
                       - 0.37043e-05f     * t2 
                       + 0.173565e-04f )  * t2 
                       - 0.487613e-04f )  * t2 
                       + 0.17343e-03f )   * t2 
                       - 0.1753062e-02f ) * t2 
                       + 0.3989422793f;
                  
                  q0 = t * ((((( 
                         0.32312e-05f     * t2 
                       - 0.142078e-04f )  * t2 
                       + 0.342468e-04f )  * t2 
                       - 0.869791e-04f )  * t2 
                       + 0.4564324e-03f ) * t2 
                       - 0.0124669441f );
                  ta0 = x-0.25f*pi;
                  by0 = __fmaf_ru(2.0f/sqrt(x),(p0*sinf(ta0),q0*cosf(ta0)));
                  sh0 = __fmaf_ru(2.0f/(pi*x),s,by0);
   
              }
              
              return (sh0);
        }
        
        
/*

     !*****************************************************************************80
!
!! STVH1 computes the Struve function H1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  SH1, the value of H1(x).
!   

*/


        __device__ float stvh1(const float x) {
          
              float a0,by1,p1,q1;
              float r,s,t,t2;
              float ta1,sh1;
              int   k,km; 
              constexpr float pi = 3.14159265358979323846264f;
              r                  = 1.0f;
              if(x <= 20.0f) {
                 
                 s0 = 0.0f;
                 a0 = -0.63661977236758134307554f;
                 for(k=1; k!=60; ++k) {
                     float tk = (float)k;
                     float t0 = (4.0f*tk*tk-1.0f);
                     r        = -r*x*x/t0;
                     s        += r;
                     if(fabsf(r) < fabsf(s)*1.0e-12f) break;
                 }
                 sh1 = a0*s;
             }
             else {
                 
                 s = 1.0f;
                 if(x <= 50.0f)
                    km = (int)(0.5f*x);
                 else
                    km = 25;
                 
                 for(k=1; k!=km; ++k) {
                     float tk = (float)k;
                     float t0 = (4.0f*tk*tk-1.0f);
                     r        = -r*t0/(x*x);
                     s        += r;
                     if(fabsf(r) < fabsf(s)*1.0e-12f) break;
                 }
                 t  = 4.0f/x;
                 t2 = t*t;
                 p1 = (((( 
                      0.42414e-05f      * t2 
                    - 0.20092e-04f  )   * t2 
                    + 0.580759e-04f )   * t2 
                    - 0.223203e-03f )   * t2 
                    + 0.29218256e-02f ) * t2 
                    + 0.3989422819;

                q1 = t * ((((( 
                   - 0.36594e-05f     * t2 
                   + 0.1622e-04f )    * t2 
                   - 0.398708e-04f )  * t2 
                   + 0.1064741e-03f ) * t2 
                   - 0.63904e-03f )   * t2 
                   + 0.0374008364f );
                ta1 = x-0.75f*pi;
                by1 = __fmaf_ru(2.0f/sqrtf(x),p1*sinf(ta1),q1*cosf(ta1));
                sh1 = __fmaf_ru(0.63661977236758134307554f,(1.0f+s/(x*x)),by1);
             }
             
             return (sh1);
        }
        
        
        
/*
   
      !*****************************************************************************80
!
!! STVL0 computes the modified Struve function L0(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  SL0, the function value.
!

*/


        __device__ float stvl0(const float x) {
        
              float a0,a1,bi0,r;
              float s,invx;
              float sl0;
              int   k,km;
              constexpr float pi =  3.14159265358979323846264f;
              s                  = 1.0f;
              r                  = 1.0f;
              if(x <= 20.0f) {
                 a0 = 2.0f*x/pi;
                 for(k=1; k!=60; ++k) {
                     float tk = (float)k;
                     float t0 = __fmaf_ru(2.0f,tk,1.0f);
                     float t1 = x/(t0*t0);
                     r        = r*t1*t1;
                     s        += r;
                     if(fabsf(r/s) < 1.0e-12f) break;
                 }
                 sl0 = a0*s;
              }  
              else {
                 
                 if(x < 50.0f) 
                    km = (int)(0.5f*(x+1.0f));
                 else
                    km = 25;   
                 for(k=1; k!=km; ++k) {
                     float tk = (float)k;
                     float t0 = (2.0f*tk-1.0f);
                     float t1 = x/(t0*t0);
                     r        = r*t1*t1;
                     s        += r;
                     if(fabsf(r/s) < 1.0e-12f) break;
                 }
                 
                 a1 = expf(x)/sqrtf(2.0f*pi*x);
                 r  = 1.0f;
                 bi0= 1.0f;
                 for(k=1; k!=16; ++k) {
                     float tk = (float)k;
                     float t0 = (2.0f*tk-1.0f);
                     r        = 0.125f*r*t0*t0/(tk*x);
                     bi0      += r;
                     if(fabsf(r/bi0) < 1.0e-12f) break;
                 }
                  bi0 = a1*bi0;
                  sl0 = -2.0f/__fmaf_ru((pi*x),s,bi0);
              }
              
              return (sl0);
        }
        
        
/*
     
        !*****************************************************************************80
!
!! STVL1 computes the modified Struve function L1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    05 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  SL1, the function value.
!

*/


        __device__ float stvl1(const float x) {
        
              float a1,bi1,r,s;
              float invxx,sl1;
              float x4;
              int   k,km;
              constexpr float pi = 3.14159265358979323846264f;
              r                  = 1.0f;
              if(x <= 20.0f) {
              
                 s = 0.0f;
                 for(k=1; k!=60; ++k) {
                     float tk = (float)k;
                     float t0 = (4.0f*tk*tk-1.0f);
                     r        = r*x*x/t0;
                     s        += r;
                     if(fabsf(r/s) < 1.0e-12f) break;
                 }
                 sl1 = 2.0f/pi*s;
              } 
              else {
                 
                 s  = 1.0f;
                 km = (int)(0.5f*x);
                 km = min(km,25);
                 invxx = 1.0f/(x*x);
                 for(k=1; k!=km; ++k) {
                     float tk = (float)k;
                     float t0 = __fmaf_ru(2.0f,tk,3.0f);
                     float t1 = __fmaf_ru(2.0f,tk,1.0f);
                     r        = r*t0*t1/invxx;
                     s        += r;
                     if(fabsf(r/s) < 1.0e-12f) break;
                 }
                 x4  = x*x*x*x;
                 sl1 = __fmaf_ru(2.0f/pi,-1.0f+invxx,3.0f*s/x4);
                 a1  = expf(x)/sqrtf(2.0f*pi*x);
                 r   = 1.0f;
                 bi1 = 1.0f;
                 for(k=1; k!=16; ++k) {
                     float tk = (float)k;
                     float t0 = 2.0f*tk-1.0f;
                     r        = -0.125f*r*t0*t0/(tk*x);
                     if(fabsf(r/bi1) < 1.0e-12f) break;
                 }
                 sl1 = sl1+a1*b1;
              }
              
              return (sl1);
        }
        
        
/*
    
      !*****************************************************************************80
!
!! STVHV computes the Struve function Hv(x) with arbitrary order v.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  V, the order of the function.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  HV, the value of Hv(x).
!
    
*/        
        
        __device__ float stvhv(const float v,
                               const float x) {
               
               float bf,bf0,bf1;
               float by0,by1,byv;
               float ga,gb,pu0,pu1;
               float qu0,qu1,r1,r2;
               float s,s0,sa,sr;
               float t0,t1,u,u0;
               float v,va,vb,vt;
               int   k,l,n;
               constexpr float pi = 3.141592653589793f;
               
               if(x<=20.0f) {
                  
                  v0 = v+1.5f;
                  ga = gamma(v0);
                  s  = 2.0f/(1.77245385090551602729817f*ga);
                  r1 = 1.0f;
                  for(k=1; k!=100; ++k) {
                      float tk = (float)k;
                      va       = tk+1.5f;
                      ga       = gamma(va);
                      vb       = v+tk+1.5f;
                      gb       = gamma(vb);
                      r1       = -r1*(0.5f*x)*(0.5f*x);
                      r2       = r1/(ga*gb);
                      s        += r2;
                      if(fabsf(r2)<fabsf(s)*1.0e-12f) break;
                  }
                  hv = pow((0.5f*x),(v+1.0f))*s;
               }
               else {
                  
                    sa = pow(0.5f*x,v-1.0f)/pi;
                    v0 = v+0.5f;
                    ga = gamma(v0);
                    s  = 1.77245385090551602729817f/ga;
                    r1 = 1.0f;
                    for(k=1; k!=12; ++k) {
                        float tk = (float)k;
                        va       = tk+0.5f;
                        ga       = gamma(va);
                        vb       = -tk+v+0.5f;
                        gb       = gamma(vb);
                        r1       = r1/(0.5f*x)*(0.5f*x);
                        s        = s+r1*ga/gb;
                    }
                    
                    s0 = sa*s;
                    u  = fabsf(v);
                    n  = (int)u;
                    u0 = u-n;
                    
                    for(l=0; l!=1; ++l) {
                        float tl = (float)l;
                        vt       = 4.0f*(u0+l)*(u0+l);
                        r1       = 1.0f;
                        pu1      = 1.0f;
                        for(k=1; k!=12; ++k) {
                            float tk   = (float)k;
                            float tmp0 = (4.0f*tk-3.0f);
                            float tmp1 = (4.0f*tk-1.0f);
                            float tmp2 = (2.0f*tk-1.0f);
                            r1         = -0.0078125f*r1*tmp0*tmp1/(tmp2*tk*x*x);
                            pu1        += r1;
                        }
                        qu1 = 1.0f;
                        r2  = 1.0f;
                        for(k=1; k!=12; ++k) {
                            float tk   = (float)k;
                            float tmp0 = (4.0f*tk-1.0f);
                            float tmp1 = __fmaf_ru(4.0f,tk,1.0f);
                            float tmp2 = __fmaf_ru(2.0f,tk,1.0f);
                            r2         = -0.0078125f*r2*tmp0*tmp1/(tmp2*tk*x*x);
                            pu1        += r1;
                        }
                        qu1 = 0.125f*(vt-1.0f)/x*qu1;
                        if(l==0) {
                           pu0 = pu1;
                           qu0 = qu1;
                        }
                    }
                    
                    t0  = x- __fmaf_ru(0.5f,u0,0.25f)*pi;
                    t1  = x- __fmaf_ru(0.5f,u0,0.75f)*pi;
                    sr  = sqrtf(2.0f/(pi*x));
                    by0 = sr*__fmaf_ru(pu0,sinf(t0),qu0*cosf(t0));
                    by1 = sr*__fmaf_ru(pu1,sinf(t1),qu1*cosf(t1));
                    bf0 = by0;
                    bf1 = by1;
                    if(n==0)
                       byv = by0;
                    else if(n==1)
                       byv = by1;
                    else
                       byv = bf;
                    hv = byv + s0;
               }
                 
                return (hv);                     
       }
       
       
/*
     !*****************************************************************************80
!
!! STVLV computes the modified Struve function Lv(x) with arbitary order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  V, the order of Lv(x).
!
!    Input, real(kind=sp) ::  X, the argument of Lv(x).
!
!    Output, real(kind=sp) ::  SLV, the value of Lv(x).
!
*/


        __device__ float stvlv(const float v,
                               const float x) {
                               
              
              float bf,bf0,bf1,biv;
              float biv0,ga,gb,r;
              float r1,r2,s,s0;
              float sa,u,u0,v;
              float v0,va,vb,vt;
              float slv;
              int k,l,n;
              constexpr float pi =  3.14159265358979323846264f;
              
              if(x==0.0f) {
              
                  if(-1.0f<v || (int)v-v == 0.5f)
                      slv = 0.0f;
                  else if(v < -1.0f)
                      slv = ::cuda::std::numeric_limits<float>::max();
                  else if(v == -1.0f)
                      slv = 0.63661977236758134307554f;
                      
             }    
              else if(x <= 40.0f) {
                  
                  v0 = v+1.5f;
                  ga = gamma(v0);
                  s  = 2.0f/(1.77245385090551602729817f*ga);
                  r1 = 1.0f;
                  for(k=1; k!=100; ++k) {
                      float tk = (float)k;
                      va       = tk+1.5f;
                      ga       = gamma(va);
                      vb       = v+t4+1.5f;
                      gb       = gamma(vb);
                      r1       = r1*(0.5f*x)*(0.5f*x);
                      r2       = r1/(ga*gb);
                      s        += r2;
                      if(fabsf(r2/s)<1.0e-12f) break;
                  }
                  slv = pow(0.5f*x,v+1.0f)*s;
             }   
             else {
                  
                  sa = -1.0f/pi*pow(0.5f*x,v-1.0f);
                  v0 = v+0.5f;
                  ga = gamma(v0);
                  s  = -1.77245385090551602729817f/ga;
                  r1 = -1.0f;
                  for(k=1; k!=12; ++k) {
                      float tk = (float)k;
                      va       = tk+0.5f;
                      ga       = gamma(va);
                      vb       = -tk+v+0.5f;
                      gb       = gamma(vb);
                      r1       = r1/(0.5f*x)*(0.5f*x);
                      s        = s+r1*ga/gb;
                  }
                    s0 = sa*s;
                    u  = fabsf(v);
                    n  = (int)u;
                    u0 = u-n;
                    
                    for(l=0; l!=1; ++l) {
                    
                        float tl = (float)l;
                        vt       =  u0+tl;
                        r        = 1.0f;
                        biv      = 1.0f;
                        for(k=1; k!=16; ++k) {
                            float tk = (float)k;
                            float t0 = 2.0f*tk-1.0f;
                            float t1 = 4.0f*vt*vt;
                            r        = -0.125f*r*t1*t0/(tk*x);
                            biv      += r;
                            if(fabsf(r/biv)<1.0e-12f) break;
                        }
                        if(l==0) biv0 = biv;
                    }
                    
                      bf0 = biv0;
                      bf1 = biv;
                      for(k=2; k!=n; ++k) {
                          float tk = (float)k;
                          float t0 = tk-1.0f+u0;
                          bf       = -2.0f*t0/x*bf1+bf0;
                          bf0      = bf1;
                          bf1      = bf;
                      }
                      
                      if(n==0)
                         biv = biv0;
                      else if(1<n)
                         biv = bf;
                      slv = expf(x)/sqrtf(2.0f*pi*x)*biv+s0;
             }    
             
             return (slv);          
      }
      
      
/*
     !*****************************************************************************80
!
!! CALCI0 computes various I0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the first kind
!    and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  If JINT = 1, then
!    the argument must be less than XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = I0(x);
!    2, RESULT = exp(-x) * I0(x);
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    1, I0(x);
!    2, exp(-x) * I0(x);
!
*/
       
        __device__ float calci0(const int jint,
                                const float arg) {
               
              const float p[15] = {
                     -5.2487866627945699800e-18f,-1.5982226675653184646e-14f,
                     -2.6843448573468483278e-11f,-3.0517226450451067446e-08f,
                     -2.5172644670688975051e-05f,-1.5453977791786851041e-02f, 
                     -7.0935347449210549190e+00f,-2.4125195876041896775e+03f, 
                     -5.9545626019847898221e+05f,-1.0313066708737980747e+08f,
                     -1.1912746104985237192e+10f,-8.4925101247114157499e+11f, 
                     -3.2940087627407749166e+13f,-5.5050369673018427753e+14f, 
                     -2.2335582639474375249e+15f};
             const float q[5] =   {
                     -3.7277560179962773046e+03f, 6.5158506418655165707e+06f,
                     -6.5626560740833869295e+09f, 3.7604188704092954661e+12f,
                     -9.7087946179594019126e+14f};
             const float pp[8] =  {
                     -3.9843750000000000000e-01f, 2.9205384596336793945e+00f,
                     -2.4708469169133954315e+00f, 4.7914889422856814203e-01f,
                     -3.7384991926068969150e-03f,-2.6801520353328635310e-03,
                      9.9168777670983678974e-05f,-2.1877128189032726730e-06f};
             const float qq[7] =  {
                     -3.1446690275135491500e+01f, 8.5539563258012929600e+01f,
                     -6.0228002066743340583e+01f, 1.3982595353892851542e+01f,
                     -1.1151759188741312645e+00f, 3.2547697594819615062e-02f,
                     -5.5194330231005480228e-04f};
             constexpr float one   = 1.0f;
             constexpr float one5  = 15.0f;
             constexpr float exp40 = 2.353852668370199854e+17f;
             constexpr float forty = 40.0f;
             constexpr float rec15 = 6.6666666666666666666e-2f;
             constexpr float two25 = 225.0f;
             constexpr float xsmall= 5.55e-17f;
             constexpr float xinf  = 3.4028235e+38f;
             constexpr float xmax  = 713.986e+00f;
             float           a,b,sump,sumq;
             float           xx,result;
             
             
             x = fabsf(arg);
             if(x<small) {
                result = one;
             }
             else if(x<one5) {
                xx   = x*x;
                sump = p[0];
                sump = __fmaf_ru(sump,xx,p[1]);
                sump = __fmaf_ru(sump,xx,p[2]);
                sump = __fmaf_ru(sump,xx,p[3]);
                sump = __fmaf_ru(sump,xx,p[4]);
                sump = __fmaf_ru(sump,xx,p[5]);
                sump = __fmaf_ru(sump,xx,p[6]);
                sump = __fmaf_ru(sump,xx,p[7]);
                sump = __fmaf_ru(sump,xx,p[8]);
                sump = __fmaf_ru(sump,xx,p[9]);
                sump = __fmaf_ru(sump,xx,p[10]);
                sump = __fmaf_ru(sump,xx,p[11]);
                sump = __fmaf_ru(sump,xx,p[12]);
                sump = __fmaf_ru(sump,xx,p[13]);
                sump = __fmaf_ru(sump,xx,p[14]); 
                xx   -= two25;
                sumq = 
                 __fmaf_ru(
                 __fmaf_ru(
                 __fmaf_ru(
                 __fmaf_ru(xx + q[0],xx,q[1]),
                                    ,xx,q[2]),
                                    ,xx,q[3]),   
                                    ,xx,q[4]);     
                result = sump/sumq;
                if(jint==2) result *= expf(-x);  
             }
             else if(one5<=x) {
                 if(jint==1 && xmax<x) 
                    result = xinf;
                 else {
                    xx   = one/x-rec15;
                    sump = 
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru( 
                        __fmaf_ru( 
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(pp[0],xx,pp[1]),
                                        xx,pp[2]),
                                        xx,pp[3]),
                                        xx,pp[4]),
                                        xx,pp[5]), 
                                        xx,pp[6]),  
                                        xx,pp[7]);   
                    sumq = 
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(xx+qq[0],xx,qq[1]),    
                                          xx,qq[2]),    
                                          xx,qq[3]),
                                          xx,qq[4]),
                                          xx,qq[5]),
                                          xx,qq[6]);
                    result = sump/sumq;
                    if(jint==2) 
                       result = (result-pp[0])/sqrtf(x);  
                    else {
                       if(x<=(xmax-one5)) {
                          a = expf(x);
                          b = one;
                       }
                       else {
                          a = expf(x-forty);
                          b = exp40;
                       }
                       result = ((result*a-pp[0]*a)/sqrtf(x))*b;
                    }                 
               }
                 
         }
         
          return (result);
     }       
     
     
/*
    
       !*****************************************************************************80
!
!! CALCI1 computes various I1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functioons of the first kind
!    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  If JINT = 1, then
!    the argument must be less than XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = I1(x);
!    2, RESULT = exp(-x) * I1(x);
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    1, I1(x);
!    2, exp(-x) * I1(x);
!

*/  


        __device__ float calci1(const float arg,
                                const int   jint) {
                                
             const float p[15] = {
                      1.9705291802535139930e-19f,-6.5245515583151902910e-16f, 
                     -1.1928788903603238754e-12f,-1.4831904935994647675e-09f, 
                     -1.3466829827635152875e-06f,-9.1746443287817501309e-04f, 
                     -4.7207090827310162436e-01f,-1.8225946631657315931e+02f, 
                     -5.1894091982308017540e+04f,-1.0588550724769347106e+07f, 
                     -1.4828267606612366099e+09f,-1.3357437682275493024e+11f, 
                     -6.9876779648010090070e+12f,-1.7732037840791591320e+14f, 
                     -1.4577180278143463643e+15f};
             const float q[5]  = {
                     -4.0076864679904189921e+03f, 7.4810580356655069138e+06f,
                     -8.0059518998619764991e+09f, 4.8544714258273622913e+12f,
                     -1.3218168307321442305e+15f};   
             const float pp[8] = {
                     -6.0437159056137600000e-02f, 4.5748122901933459000e-01f, 
                     -4.2843766903304806403e-01f, 9.7356000150886612134e-02f, 
                     -3.2457723974465568321e-03f,-3.6395264712121795296e-04f, 
                      1.6258661867440836395e-05f,-3.6347578404608223492e-07f};  
             const float qq[6] = {
                     -3.8806586721556593450e+00f, 3.2593714889036996297e+00f, 
                     -8.5017476463217924408e-01f, 7.4212010813186530069e-02f, 
                     -2.2835624489492512649e-03f, 3.7510433111922824643e-05f};   
             constexpr float one   = 1.0f;
             constexpr float one5  = 15.0f;
             constexpr float exp40 = 2.353852668370199854e+17f;
             constexpr float forty = 40.0f;
             constexpr float rec15 = 6.6666666666666666666e-2f;
             constexpr float two25 = 225.0f;
             constexpr float half  = 0.5f;
             constexpr float zero  = 0.f;
             constexpr float xsmall= 5.55e-17f;
             constexpr float xinf  = 3.4028235e+38f;
             constexpr float xmax  = 713.987f;
             constexpr float pbar  = 3.98437500e-01f;
             float           a,b,xx;
             float           result;
             
             x = fabsf(arg);
             if(x<xsmall) {
                result = half*x;
             }
             else if(x<one5) {
                xx   = x*x;
                sump = p[0];
                sump = __fmaf_ru(sump,xx,p[1]);
                sump = __fmaf_ru(sump,xx,p[2]);
                sump = __fmaf_ru(sump,xx,p[3]);
                sump = __fmaf_ru(sump,xx,p[4]);
                sump = __fmaf_ru(sump,xx,p[5]);
                sump = __fmaf_ru(sump,xx,p[6]);
                sump = __fmaf_ru(sump,xx,p[7]);
                sump = __fmaf_ru(sump,xx,p[8]);
                sump = __fmaf_ru(sump,xx,p[9]);
                sump = __fmaf_ru(sump,xx,p[10]);
                sump = __fmaf_ru(sump,xx,p[11]);
                sump = __fmaf_ru(sump,xx,p[12]);
                sump = __fmaf_ru(sump,xx,p[13]);
                sump = __fmaf_ru(sump,xx,p[14]); 
                xx   = xx-two25;
                sumq = __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(xx + q[0],xx,q[1]),
                                          ,xx,q[2]),
                                          ,xx,q[3]),   
                                          ,xx,q[4]);  
                result = (sump/sumq)*x;
                if(jint==2) result *= exp(-x);         
          }
          else if(jint==1 && xmax<x) {
                
                result = xinfl
          }
          else {
                
                xx   = one/x-rec15;
                sump = 
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru( 
                        __fmaf_ru( 
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(pp[0],xx,pp[1]),
                                        xx,pp[2]),
                                        xx,pp[3]),
                                        xx,pp[4]),
                                        xx,pp[5]), 
                                        xx,pp[6]),  
                                        xx,pp[7]); 
                sumq =    
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(
                       __fmaf_ru(xx+qq[0],xx,qq[1]),    
                                          xx,qq[2]),    
                                          xx,qq[3]),
                                          xx,qq[4]),
                                          xx,qq[5]);
                result = sump/sumq;
                if(jint!=1) 
                   result = (result+pbar)/sqrtf(x);
                else {
                     if(xmax-one5<x) {
                        a = expf(x-forty);
                        b = exp40;
                     }
                     else {
                        a = expf(x);
                        b = one;
                     }
                     result = ((__fmaf_ru(result,a,pbar*a)/sqrtf(x))*b;

                }
          }
          
          if(arg<zero) result = -result;
          return (result);
      }    
      
      
/*
     !*****************************************************************************80
!
!! CALCK0 computes various K0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K0(x);
!    2, RESULT = exp(x) * K0(x);
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
!    
*/ 


        __device__ float calck0(const float arg,
                                const int   jint) {
                                
               const float p[6] = { 
                            5.8599221412826100000e-04f, 1.3166052564989571850e-01f,
                            1.1999463724910714109e+01f, 4.6850901201934832188e+02f,
                            5.9169059852270512312e+03f, 2.4708152720399552679e+03};
               const float q[2] = {
                            -2.4994418972832303646e+02f, 2.1312714303849120380e+04f};   
               const float f[4] = {
                            -1.6414452837299064100e+00f,-2.9601657892958843866e+02f,
                            -1.7733784684952985886e+04f,-4.0320340761145482298e+05f}; 
               const float g[3] = {
                            -2.5064972445877992730e+02f, 2.9865713163054025489e+04f,
                                                        -1.6128136304458193998e+06f};     
               const float pp[10] = {
                             1.1394980557384778174e+02f, 3.6832589957340267940e+03f,
                             3.1075408980684392399e+04f, 1.0577068948034021957e+05f,
                             1.7398867902565686251e+05f, 1.5097646353289914539e+05f,
                             7.1557062783764037541e+04f, 1.8321525870183537725e+04f,
                             2.3444738764199315021e+03f, 1.1600249425076035558e+02f};
               const float qq[10] = {
                             2.0013443064949242491e+02f, 4.4329628889746408858e+03f, 
                             3.1474655750295278825e+04f, 9.7418829762268075784e+04f, 
                             1.5144644673520157801e+05f, 1.2689839587977598727e+05f, 
                             5.8824616785857027752e+04f, 1.4847228371802360957e+04f, 
                             1.8821890840982713696e+03f, 9.2556599177304839811e+01f};
               constexpr float one    = 1.0f;
               constexpr float zero   = 0.0f;
               constexpr float xsmall = 5.96e-08f;
               constexpr float xinf   = 3.4028235e+38f;
               constexpr float xmax   = 705.342f;
               float  sumf,sumg,sump,sumq;
               float  temp,xx;
               
               x = arg;
               if(zero<x) {
                  if(x<=one) {
                     temp = logf(x);
                     if(x<xsmll)
                        result = p[5]/q[1]-temp;
                     else {
                        xx = x*x;
                        sump = 
                               __fmaf_ru(
                               __fmaf_ru(
                               __fmaf_ru(
                               __fmaf_ru(
                               __fmaf_ru(p[0],xx,p[1]),
                                              xx,p[2]),
                                              xx,p[3]),
                                              xx,p[4]),
                                              xx,p[5]);
                        sumq = __fmaf_ru(xx+q[0],xx,q[1]);
                        sumf = __fmaf_ru(
                               __fmaf_ru(
                               __fmaf_ru(f[0],xx,f[1]),
                                              xx,f[2]),
                                              xx,f[3]);
                        sumg   = __fmaf_ru(__fmaf_ru(xx+g[0],xx,g[1]),xx,g[2]);
                        result = sump/sumq-xx*sumf*temp/sumg-temp;
                        if(jint==2) result *= exp(x); 
                     }
                     else if(jint==1 && xmax<x)
                           result = zero
                     else {
                           
                           xx   = one/x;
                           sump = pp[0];
                           sump = __fmaf_ru(sump,xx,pp[1]);
                           sump = __fmaf_ru(sump,xx,pp[2]);
                           sump = __fmaf_ru(sump,xx,pp[3]);
                           sump = __fmaf_ru(sump,xx,pp[4]);
                           sump = __fmaf_ru(sump,xx,pp[5]);
                           sump = __fmaf_ru(sump,xx,pp[6]);
                           sump = __fmaf_ru(sump,xx,pp[7]);
                           sump = __fmaf_ru(sump,xx,pp[8]);
                           sump = __fmaf_ru(sump,xx,pp[9]);
                           sumq = xx;
                           sumq = (sumq+qq[0])*xx;
                           sumq = (sumq+qq[1])*xx;
                           sumq = (sumq+qq[2])*xx;
                           sumq = (sumq+qq[3])*xx;
                           sumq = (sumq+qq[4])*xx;
                           sumq = (sumq+qq[5])*xx;
                           sumq = (sumq+qq[6])*xx;
                           sumq = (sumq+qq[7])*xx;
                           sumq = (sumq+qq[8])*xx;
                           sumq = sumq+qq[9];
                           result= sump/qumq/sqrtf(x);
                           if(jint==1) result *= expf(-x);
                     }
                  }
                  
                   
               }
               else {
                  result = xinf;
               }
               
               return (result);
       }
       
       
/*
    
     !*****************************************************************************80
!
!! CALCK1 computes various K1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K1(x);
!    2, RESULT = exp(x) * K1(x);
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);
!

*/ 


            __device__ float calck1(const float arg,
                                    const int   jint) {
                                    
                      const float p[5] = {
                            4.8127070456878442310e-1f, 9.9991373567429309922e+1f,
                            7.1885382604084798576e+3f, 1.7733324035147015630e+5f, 
                            7.1938920065420586101e+5f};
                      const float q[3] = {
                            -2.8143915754538725829e+2f, 3.7264298672067697862e+4f,
                            -2.2149374878243304548e+6f};
                      const float f[5] = {
                             -2.2795590826955002390e-1f,-5.3103913335180275253e+1f, 
                             -4.5051623763436087023e+3f,-1.4758069205414222471e+5f,
                             -1.3531161492785421328e+6f};
                      const float g[3] = {
                             -3.0507151578787595807e+2f, 4.3117653211351080007e+4f,
                             -2.7062322985570842656e+6f};
                      const float pp[11] = {
                              6.4257745859173138767e-2f, 7.5584584631176030810e+0f, 
                              1.3182609918569941308e+2f, 8.1094256146537402173e+2f, 
                              2.3123742209168871550e+3f, 3.4540675585544584407e+3f, 
                              2.8590657697910288226e+3f, 1.3319486433183221990e+3f, 
                              3.4122953486801312910e+2f, 4.4137176114230414036e+1f, 
                              2.2196792496874548962e+0f};
                      const float qq[9]  = {
                              3.6001069306861518855e+1f, 3.3031020088765390854e+2f,
                              1.2082692316002348638e+3f, 2.1181000487171943810e+3f, 
                              1.9448440788918006154e+3f, 9.6929165726802648634e+2f, 
                              2.5951223655579051357e+2f, 3.4552228452758912848e+1f, 
                              1.7710478032601086579e+0f};
                      constexpr float one = 1.0f;
                      constexpr float zero= 0.0f;
                      constexpr float xmax= 705.343f;
                      float           sumf,sumg,sumq,sump;
                      float           xx,result;
                      
                      x = arg;
                      if(x<=one) {
                         
                         xx    = x*x;
                         sump  = 
                                 __fmaf_ru(
                                 __fmaf_ru(
                                 __fmaf_ru(
                                 __fmaf_ru(
                                 __fmaf_ru(p[0],xx,p[1]),
                                                xx,p[2]),
                                                xx,p[3]),
                                                xx,p[4]),
                                                xx,q[2]);
                         sumq  = 
                                 __fmaf_ru(
                                 __fmaf_ru(xx+q[0],xx,q[1]),
                                                   xx,q[2]);
                         sumf  = 
                                 __fmaf_ru(
                                 __fmaf_ru(
                                 __fmaf_ru(
                                 __fmaf_ru(f[0],xx,f[1]),  
                                                xx,f[2]),  
                                                xx,f[3]), 
                                                xx,f[4]);
                         sumg  = 
                                 __fmaf_ru(
                                 __fmaf_ru(xx+g[0],xx,g[1]),
                                                   xx,g[2]);
                         result = (xx*logf(x)*sumf/sumg+sump/sumq)/x;
                         if(jint==2) result *= expf(x);
                  }
                  else if(jint==1 && xmax<x) {
                         
                         result = zero;
                  }
                  else {
                         
                         xx   = one/x;
                         sump = pp[0];
                         sump = __fmaf_ru(sump,xx,pp[1]);
                         sump = __fmaf_ru(sump,xx,pp[2]);
                         sump = __fmaf_ru(sump,xx,pp[3]);
                         sump = __fmaf_ru(sump,xx,pp[4]);
                         sump = __fmaf_ru(sump,xx,pp[5]);
                         sump = __fmaf_ru(sump,xx,pp[6]);
                         sump = __fmaf_ru(sump,xx,pp[7]);
                         sump = __fmaf_ru(sump,xx,pp[8]);
                         sump = __fmaf_ru(sump,xx,pp[9]);
                         sump = __fmaf_ru(sump,xx,pp[10]);
                         sumq = xx;
                         sumq = (sumq+qq[0])*xx;
                         sumq = (sumq+qq[1])*xx;
                         sumq = (sumq+qq[2])*xx;
                         sumq = (sumq+qq[3])*xx;
                         sumq = (sumq+qq[4])*xx;
                         sumq = (sumq+qq[5])*xx;
                         sumq = (sumq+qq[6])*xx;
                         sumq = (sumq+qq[7])*xx;
                         sumq += qq[8];
                         result= sumq/sumq/sqrtf(x);
                         if(jint==1) result *= expf(-x);
                  }
                  
                  return (result);
          }  
          
          
/*
     
      !*****************************************************************************80
!
!! ERROR evaluates the error function.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  ERR, the function value.
!
  
*/


        __device__ float error(const float x) {
        
               float c0,er,r,x2;
               float invx,err;
               int32_t k;
               constexpr float pi = 3.14159265358979323846264f;
               constexpr float eps= 1.0e-15f;
               
               x2 = x*x;
               if(fabsf(x)<3.5f) {
                  er = 1.0f;
                  r  = 1.0f;
                  for(k=1; k<=50; ++k) {
                      float tk = (float)k;
                      float t0 = tk+0.5f;
                      r        = r*x2/t0;
                      er       +=r;
                      if(fabsf(r)<=fabsf(er)*eps) break;
                  }
                  c0 = 1.12837916709551257389616f*x*expf(-x2);
                  err= c0*er;
               } 
               else {
                  invx =1.0f/x2;
                  er   = 1.0f;
                  r    = 1.0f;
                  r    = -r*0.5f*invx;
                  er   += r;
                  r    = -r*1.5f*invx;
                  er   += r;
                  r    = -r*2.5f*invx;
                  er   += r;
                  r    = -r*3.5f*invx;
                  er   += r;
                  r    = -r*4.5f*invx;
                  er   += r;
                  r    = -r*5.5f*invx;
                  er   += r;
                  r    = -r*6.5f*invx;
                  er   += r;
                  r    = -r*7.5f*invx;
                  er   += r;
                  r    = -r*8.5f*invx;
                  er   += r;
                  r    = -r*9.5f*invx;
                  er   += r;
                  r    = -r*10.5f*invx;
                  er   += r;
                  r    = -r*11.5f*invx;
                  er   += r;
                  c0   = expf(-x2)/(fabsf(x)*1.77245385090551602729817f);
                  err  = 1.0f-c0*er;
                  if(x<0.0f) err = -err;
                  return (err);
               }
       }
       
       
/*

   !*****************************************************************************80
!
!! E1XB computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  E1, the function value.
!
*/
       
         __device__ float e1xb(const float x) {
         
                float ga,r,t0,e1;
                int32_t k,m;
                
                if(x==0.0) {
                   e1 = 3.4028235e+38f);
                }
                else if(x<=1.0f) {
                   e1 = 1.0f;
                   r  = 1.0f;
                   for(k=1; k<=25; ++k) {
                       float tk = (float)k;
                       float t1 = tk+1.0f;
                       r        = -r*tk*x/(t1*t1);
                       e1       += r;
                       if(fabsf(r)<=fabsf(e1)*1.0e-15f) break;
                   }
                   ga = 0.5772156649015328f;
                   e1 = -ga-logf(x)+x*e1;
                }
                else {
                   m  = 20+(int32_t)(80.0f/x);
                   t0 = 0.0f;
                   for(k=m; k>=1; --k) {
                       float tk = (float)k;
                       t0       = tk/(1.0f+tk/(x+t0));
                   }
                   t  = 1.0f/(x+t0);
                   e1 = expf(-x)*t;
              }
              return (e1);
         }   
         
         
/*

      
!*****************************************************************************80
!
!! E1XA computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  E1, the function value.
!

*/


        __device__ float e1xa(const float x) {
        
                float es1,es2;
                
                if(x==0.0f) {
                    e1 = 3.4028235e+38f
                }
                else if(x<=1.0f) {
                    e1 = -logf(x)+((((
                          1.07857e-03f  * x
                         -9.76004e-03f) * x
                         + 5.519968e-02f ) * x
                         - 0.24991055f) * x
                         + 0.99999193f) * x 
                         - 0.57721566f;
                }
                else {
                    es1 = 
                     ((( x 
                         + 8.5733287401f) * x
                         +18.059016973f)  * x
                         + 8.6347608925f) * x 
                         + 0.2677737343f;

                   es2 = ((( x
                         +  9.5733223454f) * x
                         + 25.6329561486f) * x 
                         + 21.0996530827f) * x 
                         + 3.9584969228f;

                   e1 = expf(-x)/x*es1/es2;
                }
                return (e1);
        }
          
/*
       !*****************************************************************************80
!
!! CALERF computes various forms of the error function.
!
!  Discussion:
!
!    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
!    for a real argument x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Function,
!    Mathematics of Computation,
!    Volume 23, Number 107, July 1969, pages 631-638.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  If JINT is 1, the
!    argument must be less than XBIG.  If JINT is 2, the argument
!    must lie between XNEG and XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = erf(x);
!    1, RESULT = erfc(x) = 1 - erf(x);
!    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    0, erf(x);
!    1, erfc(x);
!    2, exp(x*x)*erfc(x).
!
*/


        __device__ float calerf(const float arg,
                                const int   jint) {
                                
               const float a[5] = {
                         3.16112374387056560e+00f,1.13864154151050156e+02f,
                         3.77485237685302021e+02f,3.20937758913846947e+03f, 
                         1.85777706184603153e-1f};                      
               const float b[4] = {
                         2.36012909523441209e+01f,2.44024637934444173e+02f,
                         1.28261652607737228e+03f,2.84423683343917062e+03f};
               const float c[9] = {
                         5.64188496988670089e-1f, 8.88314979438837594e+00f,
                         6.61191906371416295e+01f,2.98635138197400131e+02f, 
                         8.81952221241769090e+02f,1.71204761263407058e+03f, 
                         2.05107837782607147e+03f,1.23033935479799725e+03f,
                         2.15311535474403846e-8f};
               const float d[8] = {
                         1.57449261107098347e+01f,1.17693950891312499e+02f,
                         5.37181101862009858e+02f,1.62138957456669019e+03f, 
                         3.29079923573345963e+03f,4.36261909014324716e+03f,
                         3.43936767414372164e+03f,1.23033935480374942e+03f};
               const float p[6] = {
                         3.05326634961232344e-1f,3.60344899949804439e-1f,
                         1.25781726111229246e-1f,1.60837851487422766e-2f, 
                         6.58749161529837803e-4f,1.63153871373020978e-2f};
               const float q[5] = {
                         2.56852019228982242e+00f,1.87295284992346047e+00f,
                         5.27905102951428412e-1f, 6.05183413124413191e-2f,
                         2.33520497626869185e-3f};
               constexpr float four  = 4.0f;
               constexpr float one   = 1.0f;
               constexpr float half  = 0.5f;
               constexpr float two   = 2.0f;
               constexpr float zero  = 0.0f;
               constexpr float sqrpi = 5.6418958354775628695e-1f;
               constexpr float thresh= 0.46875f;
               constexpr float sixth = 16.0f;
               constexpr float xinf  = 3.4028235e+38f;
               constexpr float xneg  = -26.628f;
               constexpr float xsmall= 1.11e-16f;
               constexpr float xbig  = 26.543f;
               constexpr float xhuge = 6.71e+7f;
               float           del,sixten,y,ysq;
               float           result;
               
               x = arg;
               y = fabsf(x);
               if(y<=thresh) {
                  
                  ysq = zero;
                  if(xsmall<y) ysq = y*y;
                  xnum   = a[4]*ysq;
                  xden   = ysq;
                  xnum   = (xnum+a[0])*ysq;
                  xden   = (xden+b[0])*ysq;
                  xnum   = (xnum+a[1])*ysq;
                  xden   = (xden+b[1])*ysq;
                  xnum   = (xnum+a[2])*ysq;
                  xden   = (xden+b[2])*ysq;
                  result = x*((xnum+a[3])/(xden+b[3]));
                  if(jint!=0) result -= one;
                  if(jint==2) result *= expf(ysq);
                  return (result); 
               }
               else if(y<=four) {
                  
                  xnum = c[8]*y;
                  xden = y;
                  xnum = (xnum+c[0])*y;
                  xden = (xden+d[0])*y;
                  xnum = (xnum+c[1])*y;
                  xden = (xden+d[1])*y;
                  xnum = (xnum+c[2])*y;
                  xden = (xden+d[2])*y;
                  xnum = (xnum+c[3])*y;
                  xden = (xden+d[3])*y;
                  xnum = (xnum+c[4])*y;
                  xden = (xden+d[4])*y;
                  xnum = (xnum+c[5])*y;
                  xden = (xden+d[5])*y;
                  xnum = (xnum+c[6])*y;
                  xden = (xden+d[6])*y;
                  result = ((xnum+c[7])/(xden+d[7));
                  if(jint!=2) {
                     ysq    = (int)(y*sixten)/sixten;
                     del    = (y-ysq)*(y+ysq);
                     result = expf(-ysq*ysq)*expf(-del)*result; 
                  }
               }
               else {
                   
                   result = zero;
                   if(xbig<=y) {
                      if(jint!=2) goto L300;
                      if(xhuge<=y) {
                         result = sqrpi/y;
                         goto L300;
                      }
                   }
                   
                   ysq = one/(y*y);
                   xnum= p[5]*ysq;
                   xden= ysq;
                   xnum= (xnum+p[0])*ysq;
                   xden= (xden+q[0])*ysq;
                   xnum= (xnum+p[1])*ysq;
                   xden= (xden+q[1])*ysq;
                   xnum= (xnum+p[2])*ysq;
                   xden= (xden+q[2])*ysq;
                   xnum= (xnum+p[3])*ysq;
                   xden= (xden+q[3])*ysq;
                   result = ysq*(xnum+p[4])/(xden+q[4]);
                   result = (sqrpi-result)/y;
                   if(jint!=2) {
                     ysq    = (int)(y*sixten)/sixten;
                     del    = (y-ysq)*(y+ysq);
                     result = expf(-ysq*ysq)*expf(-del)*result; 
                   }
               }
L300:
                   if(jint==0) {
                      result = (half-result)+half;
                      if(x<zero) result = -result;
                   }  
                   else if(jint==1) {
                      if(x<zero) result = two-result;
                   }   
                   else {
                      if(x<zero) {
                         if(x<xneg) 
                            result = xinf;
                         else {
                            ysq    = (int)(y*sixten)/sixten;
                            del    = (y-ysq)*(y+ysq);
                            result = expf(-ysq*ysq)*expf(-del)*result; 
                            result = (y+y)-result;
                         }
                      }
                   }   
                   
                   return (result);     
       }
       
       
       
/*
    !*****************************************************************************80
!
!! CALJY0 computes various J0 and Y0 Bessel functions.
!
!  Discussion:
!
!    This routine computes zero-order Bessel functions of the first and
!    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!    for Y0, and |X| <= XMAX for J0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J0(x);
!    1, RESULT = Y0(x);
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    0, J0(x);
!    1, Y0(x);
*/  
     
       
        __device__ float caljy0(const float arg,
                                const int   jint) {
                                
               
               const float plg[4] = {
                      -2.4562334077563243311e+01f,2.3642701335621505212e+02f,
                      -5.4989956895857911039e+02f,3.5687548468071500413e+02f};
               const float qlg[4] = {
                       -3.5553900764052419184e+01f,1.9400230218539473193e+02f, 
                       -3.3442903192607538956e+02f,1.7843774234035750207e+02f};
               const float pj0[7] = {
                        6.6302997904833794242e+06f,-6.2140700423540120665e+08f,
                        2.7282507878605942706e+10f,-4.1298668500990866786e+11f,
                       -1.2117036164593528341e-01f, 1.0344222815443188943e+02f,
                       -3.6629814655107086448e+04f};
               const float qj0[5] = {
                        4.5612696224219938200e+05f, 1.3985097372263433271e+08f,
                        2.6328198300859648632e+10f, 2.3883787996332290397e+12f,
                        9.3614022392337710626e+02f};                      
               const float pj1[8] = {
                        4.4176707025325087628e+03f, 1.1725046279757103576e+04f,
                        1.0341910641583726701e+04f,-7.2879702464464618998e+03f,
                       -1.2254078161378989535e+04f,-1.8319397969392084011e+03f, 
                        4.8591703355916499363e+01f, 7.4321196680624245801e+02f};
               const float qj1[7] = {
                        3.3307310774649071172e+02f,-2.9458766545509337327e+03f,
                        1.8680990008359188352e+04f,-8.4055062591169562211e+04f,
                        2.4599102262586308984e+05f,-3.5783478026152301072e+05f, 
                       -2.5258076240801555057e+01f};
               const float py0[6] = {
                        1.0102532948020907590e+04f,-2.1287548474401797963e+06f,
                        2.0422274357376619816e+08f,-8.3716255451260504098e+09f,
                        1.0723538782003176831e+11f,-1.8402381979244993524e+01f};
               const float qy0[6] = {
                        6.6475986689240190091e+02f, 2.3889393209447253406e+05f,
                        5.5662956624278251596e+07f, 8.1617187777290363573e+09f,
                        5.8873865738997033405e+11f};
               const float py1[7] = {
                       -1.4566865832663635920e+04f, 4.6905288611678631510e+06f, 
                       -6.9590439394619619534e+08f, 4.3600098638603061642e+10f, 
                       -5.5107435206722644429e+11f,-2.2213976967566192242e+13f, 
                        1.7427031242901594547e+01f};
               const float qy1[6] = {
                        8.3030857612070288823e+02f, 4.0669982352539552018e+05f, 
                        1.3960202770986831075e+08f, 3.4015103849971240096e+10f, 
                        5.4266824419412347550e+12f, 4.3386146580707264428e+14f};
               const float py2[8] = {
                        2.1363534169313901632e+04f,-1.0085539923498211426e+07f,
                        2.1958827170518100757e+09f,-1.9363051266772083678e+11f,
                       -1.2829912364088687306e+11f, 6.7016641869173237784e+14f,
                       -8.0728726905150210443e+15f,-1.7439661319197499338e+01f};
               const float qy2[7] = {
                         8.7903362168128450017e+02f, 5.3924739209768057030e+05f, 
                         2.4727219475672302327e+08f, 8.6926121104209825246e+10f, 
                         2.2598377924042897629e+13f, 3.9272425569640309819e+15f, 
                         3.4563724628846457519e+17f};
               const float p0[6]  = {
                         3.4806486443249270347e+03f, 2.1170523380864944322e+04f, 
                         4.1345386639580765797e+04f, 2.2779090197304684302e+04f, 
                         8.8961548424210455236e-01f, 1.5376201909008354296e+02f};
               const float q0[5]  = {
                         3.5028735138235608207e+03f, 2.1215350561880115730e+04f,
                         4.1370412495510416640e+04f, 2.2779090197304684318e+04f, 
                         1.5711159858080893649e+02f};
               const float p1[6]  = {
                        -2.2300261666214198472e+01f,-1.1183429920482737611e+02f,
                        -1.8591953644342993800e+02f,-8.9226600200800094098e+01f,
                        -8.8033303048680751817e-03f,-1.2441026745835638459e+00f};
               const float q1[5]  = {
                         1.4887231232283756582e+03f, 7.2642780169211018836e+03f,
                         1.1951131543434613647e+04f, 5.7105024128512061905e+03f,
                         9.0593769594993125859e+01f}; 
              constexpr float zero   = 0.0f;
              constexpr float one    = 1.0f;
              constexpr float three  = 3.0f;
              constexpr float four   = 4.0f;
              constexpr float eight  = 8.0f;
              constexpr float five5  = 5.5f;
              constexpr float sixty4 = 64.0f
              constexpr float oneov8 = 0.125f;
              constexpr float p17    = 1.716e-1f;
              constexpr float two56  = 256.0f;
              constexpr float cons   = -1.1593151565841244881e-1f;
              constexpr float pi2    = 6.3661977236758134308e-1f;
              constexpr float twopi  = 6.2831853071795864769f;
              constexpr float twopi1 = 6.28125f;
              constexpr float twopi2 = 1.9353071795864769253e-3f;   
              constexpr float xmax   = 1.07e+09f;
              constexpr float xsmall = 9.31e-10f;
              constexpr float xinf   = 1.7e+38f;
//!
//!  Zeroes of Bessel functions
//!
              constexpr float xj0    = 2.4048255576957727686e+0f;
              constexpr float xj1    = 5.5200781102863106496e+0f;
              constexpr float xy0    = 8.9357696627916752158e-1f;
              constexpr float xy1    = 3.9576784193148578684e+0f;
              constexpr float xy2    = 7.0860510603017726976e+0f;
              constexpr float xj01   = 616.0e+0f;
              constexpr float xj02   = -1.4244423042272313784e-03f;
              constexpr float xj11   = 1413.0e+0f;
              constexpr float xj12   = 5.4686028631064959660e-04f;
              constexpr float xy01   = 228.0e+0f;
              constexpr float xy02   = 2.9519662791675215849e-03f;
              constexpr float xy11   = 1013.0e+0f;
              constexpr float xy12   = 6.4716931485786837568e-04f;
              constexpr float xy21   = 1814.0e+0f;
              constexpr float xy22   = 1.1356030177269762362e-04f;
              float           arg,ax,down,prod;
              float           resj,r0,r1,up;
              float           w,wsq,xden,z,zsq;
              float           result;
              
              ax  = fabsf(arg);
              if(jint==1 && arg<=zero)
                 return (-xinf);
              else if(xmax<ax)
                 return (zero);
              if(eight<ax) goto L800;
              
              if(ax<=xsmall) {
                 if(jint==0)
                    return (zero);
                 else
                    return (pi2*logf(ax)+cons);
             }
             
              //  Calculate J0 for appropriate interval, preserving
              //!  accuracy near the zero of J0.  
              zsq = ax*ax;
              if(ax<=four) {
                 xnum = __fmaf_ru(__fmaf_ru(pj0[4],zsq,pj0[5]),zsq,pj0[6]);
                 xden = zsq+qj0[4];
                 xnum = __fmaf_ru(xnum,zsq,pj0[0]);
                 xden = __fmaf_ru(xden,zsq,qj0[0]);
                 prod = ((ax-xj01/two56)-xj02)*(ax+xj0);
              }
              else {
                 wsq  = one-zsq/sixty4;
                 xnum = __fmaf_ru(pj1[6],wsq,pj1[7]);
                 xden = wsq+qj1[6];
                 xnum = __fmaf_ru(xnum,wsq,pj1[0]);
                 xden = __fmaf_ru(xden,wsq,qj1[0]);
                 xnum = __fmaf_ru(xnum,wsq,pj1[1]);
                 xden = __fmaf_ru(xden,wsq,qj1[1]);
                 xnum = __fmaf_ru(xnum,wsq,pj1[2]);
                 xden = __fmaf_ru(xden,wsq,qj1[2]);
                 xnum = __fmaf_ru(xnum,wsq,pj1[3]);
                 xden = __fmaf_ru(xden,wsq,qj1[3]);
                 xnum = __fmaf_ru(xnum,wsq,pj1[4]);
                 xden = __fmaf_ru(xden,wsq,qj1[4]);
                 xnum = __fmaf_ru(xnum,wsq,pj1[5]);
                 xden = __fmaf_ru(xden,wsq,qj1[5]);
                 prod = (ax+xj1)*((ax-xj11/two56)-xj12);
              }
              result  = prod*xnum/xden;
              if(jint==0) return (result);
              
              /*
                     Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
                     !  where xn is a zero of Y0.
              */
              
              if(ax<=three) {
                 up = (ax-xy01/two56)-xy02;
                 xy = xy0;
              }
              else if(ax<=five5) {
                 up = (ax-xy11/two56)-xy12;
                 xy = xy1;
              }
              else {
                 up = (ax-xy21/two56)-xy22;
                 xy = xy2;
              }
              down = ax+xy;
              
              if(fabsf(up)<p17*down) {
                 w    = up/down;
                 wsq  = w*w;
                 xnum = plg[0];
                 xden = wsq+qlg[0];
                 xnum = __fmaf_ru(xnum,wsq,plg[1]);
                 xden = __fmaf_ru(xden,wsq,qlg[1]);
                 xnum = __fmaf_ru(xnum,wsq,plg[2]);
                 xden = __fmaf_ru(xden,wsq,qlg[2]);
                 xnum = __fmaf_ru(xnum,wsq,plg[3]);
                 xden = __fmaf_ru(xden,wsq,qlg[3]);
                 resj = pi2*result*w*xnum/xden;
              }
              else {
                 resj = pi2*result*logf(ax/xy);
              }
              /*
                   Now calculate Y0 for appropriate interval, preserving
                   !  accuracy near the zero of Y0.
              */
              
              if(ax<=three) {
                 xnum = __fmaf_ru(py0[5],zsq,py0[0]);
                 xden = zsq+qy0[0];
                 xnum = __fmaf_ru(xnum,zsq,py0[1]);
                 xden = __fmaf_ru(xden,zsq,qy0[1]);
                 xnum = __fmaf_ru(xnum,zsq,py0[2]);
                 xden = __fmaf_ru(xden,zsq,qy0[2]);
                 xnum = __fmaf_ru(xnum,zsq,py0[3]);
                 xden = __fmaf_ru(xden,zsq,qy0[3]);
                 xnum = __fmaf_ru(xnum,zsq,py0[4]);
                 xden = __fmaf_ru(xden,zsq,qy0[4]);
            }
             else if(ax<=five5) {
                  xnum = __fmaf_ru(py1[6],zsq,py1[0]);
                  xden = zsq+qy1[0];
                  xnum = __fmaf_ru(xnum,zsq,py1[1]);
                  xden = __fmaf_ru(xden,zsq,qy1[1]);
                  xnum = __fmaf_ru(xnum,zsq,py1[2]);
                  xden = __fmaf_ru(xden,zsq,qy1[2]);
                  xnum = __fmaf_ru(xnum,zsq,py1[3]);
                  xden = __fmaf_ru(xden,zsq,qy1[3]);
                  xnum = __fmaf_ru(xnum,zsq,py1[4]);
                  xden = __fmaf_ru(xden,zsq,qy1[4]);
                  xnum = __fmaf_ru(xnum,zsq,py1[5]);
                  xden = __fmaf_ru(xden,zsq,qy1[5]);
            }
            else {
                  xnum = __fmaf_ru(py2[7],zsq,py2[0]);
                  xden = zsq + qy2[0];
                  xnum = __fmaf_ru(xnum,zsq,py2[1]);
                  xden = __fmaf_ru(xden,zsq,qy2[1]);
                  xnum = __fmaf_ru(xnum,zsq,py2[2]);
                  xden = __fmaf_ru(xden,zsq,qy2[2]);
                  xnum = __fmaf_ru(xnum,zsq,py2[3]);
                  xden = __fmaf_ru(xden,zsq,qy2[3]);
                  xnum = __fmaf_ru(xnum,zsq,py2[4]);
                  xden = __fmaf_ru(xden,zsq,qy2[4]);
                  xnum = __fmaf_ru(xnum,zsq,py2[5]);
                  xden = __fmaf_ru(xden,zsq,qy2[5]);
                  xnum = __fmaf_ru(xnum,zsq,py2[6]);
                  xden = __fmaf_ru(xden,zsq,qy2[6]);
            }
            result = resj+up*down*xnum/xden;
            return (result);
L800:
            z = eight/ax;
            w = ax/twopi;
            w = (int)w+oneov8;
            w = (ax-w*twopi1)-w*twopi2;
            zsq = z*z;
            xnum = __fmaf_ru(p0[4],zsq,p0[5]);
            xden = zsq+q0[4];
            up   = __fmaf_ru(p1[4],zsq,p1[5]);
            down = zsq+q1[4];
            xnum = __fmaf_ru(xnum,zsq,p0[0]);
            xden = __fmaf_ru(xden,zsq,q0[0]);
            xnum = __fmaf_ru(xnum,zsq,p0[1]);
            xden = __fmaf_ru(xden,zsq,q0[1]);
            xnum = __fmaf_ru(xnum,zsq,p0[2]);
            xden = __fmaf_ru(xden,zsq,q0[2]);
            xnum = __fmaf_ru(xnum,zsq,p0[3]);
            xden = __fmaf_ru(xden,zsq,q0[3]);
            r0   = xnum/xden;
            r1   = up/down;

            if(jint==0 ){
               result = sqrtf(pi2/ax) *
                          (r0* cosf(w)-z*r1*sinf(w));
               return (result);
           }
            else {
               result = sqrtf(pi2/ax) *
                          (r0*sinf(w)+z*r1*cosf(w));
               return (result);
         }
 
      }
      
      
/*
      
     !*****************************************************************************80
!
!! CALJY1 computes various J1 and Y1 Bessel functions.
!
!  Discussion:
!
!    This routine computes first-order Bessel functions of the first and
!    second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!    for Y1, and |X| <= XMAX for J1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real(kind=sp) ::  RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J1(x);
!    1, RESULT = Y1(x);
!
!    Input, integer(kind=i4) :: JINT, chooses the function to be computed.
!    0, J1(x);
!    1, Y1(x);
! 

*/


        __device__ float caljy1(const float arg,
                                const int   jint) {
                                
 
                const float plg[4] = {
                           -2.4562334077563243311e+01f,2.3642701335621505212e+02f,
                           -5.4989956895857911039e+02f,3.5687548468071500413e+02f};
                const float qlg[4] = {
                           -3.5553900764052419184e+01f,1.9400230218539473193e+02f, 
                           -3.3442903192607538956e+02f,1.7843774234035750207e+02f};
                const float pj0[7] = {
                            9.8062904098958257677e+05f,-1.1548696764841276794e+08f,
                            6.6781041261492395835e+09f,-1.4258509801366645672e+11f, 
                           -4.4615792982775076130e+03f, 1.0650724020080236441e+01f,
                           -1.0767857011487300348e-02f};  
                const float qj0[5] = {
                            5.9117614494174794095e+05f, 2.0228375140097033958e+08f,
                            4.2091902282580133541e+10f, 4.1868604460820175290e+12f,
                            1.0742272239517380498e+03f};
                const float pj1[8] = {
                            4.6179191852758252280e+00f,-7.1329006872560947377e+03f, 
                            4.5039658105749078904e+06f,-1.4437717718363239107e+09f, 
                            2.3569285397217157313e+11f,-1.6324168293282543629e+13f, 
                            1.1357022719979468624e+14f, 1.0051899717115285432e+15f};
                const float qj1[7] = {
                            1.1267125065029138050e+06f, 6.4872502899596389593e+08f, 
                            2.7622777286244082666e+11f, 8.4899346165481429307e+13f, 
                            1.7128800897135812012e+16f, 1.7253905888447681194e+18f,
                            1.3886978985861357615e+03f};
                const float py0[7] = {
                            2.2157953222280260820e+05f,-5.9157479997408395984e+07f,
                            7.2144548214502560419e+09f,-3.7595974497819597599e+11f, 
                            5.4708611716525426053e+12f, 4.0535726612579544093e+13f, 
                           -3.1714424660046133456e+02f};
                const float qy0[6] = {
                            8.2079908168393867438e+02f, 3.8136470753052572164e+05f,
                            1.2250435122182963220e+08f, 2.7800352738690585613e+10f,
                            4.1272286200406461981e+12f, 3.0737873921079286084e+14f};
                const float py1[9] = {
                            1.9153806858264202986e+06f,-1.1957961912070617006e+09f, 
                            3.7453673962438488783e+11f,-5.9530713129741981618e+13f, 
                            4.0686275289804744814e+15f,-2.3638408497043134724e+16f, 
                           -5.6808094574724204577e+18f, 1.1514276357909013326e+19f, 
                           -1.2337180442012953128e+03f};
                const float qy1[8] = {
                            1.2855164849321609336e+03f, 1.0453748201934079734e+06f, 
                            6.3550318087088919566e+08f, 3.0221766852960403645e+11f, 
                            1.1187010065856971027e+14f, 3.0837179548112881950e+16f, 
                            5.6968198822857178911e+18f, 5.3321844313316185697e+20f};
                const float p0[6]  = {
                           -1.0982405543459346727e+05f,-1.5235293511811373833e+06f, 
                           -6.6033732483649391093e+06f,-9.9422465050776411957e+06f, 
                           -4.4357578167941278571e+06f,-1.6116166443246101165e+03f};
                const float q0[6]  = {
                           -1.0726385991103820119e+05f,-1.5118095066341608816e+06f,
                           -6.5853394797230870728e+06f,-9.9341243899345856590e+06f, 
                           -4.4357578167941278568e+06f,-1.4550094401904961825e+03f};
                const float p1[6]  = {
                            1.7063754290207680021e+03f, 1.8494262873223866797e+04f, 
                            6.6178836581270835179e+04f, 8.5145160675335701966e+04f, 
                            3.3220913409857223519e+04f, 3.5265133846636032186e+01f};
                const float q1[6]  = {
                            3.7890229745772202641e+04f, 4.0029443582266975117e+05f,
                            1.4194606696037208929e+06f, 1.8194580422439972989e+06f,
                            7.0871281941028743574e+05f, 8.6383677696049909675e+02f};
                constexpr float eight = 8.0f;
                constexpr float four  = 4.0f;
                constexpr float half  = 0.5f;
                constexpr float throv8= 0.375f;
                constexpr float pi2   = 6.3661977236758134308e-1f;
                constexpr float p17   = 1.716e-1f;
                constexpr float twopi = 6.2831853071795864769e+0f;
                constexpr float zero  = 0.0f;
                constexpr float twopi1= 6.28125f;
                constexpr float twopi2= 1.9353071795864769253e-03f;
                constexpr float two56 = 256.0e+0f;
                constexpr float rtpi2 = 7.9788456080286535588e-1f;
                constexpr float xmax  = 1.07e+09f;
                constexpr float xsmall= 9.31e-10f;
                constexpr float xinf  = 1.7e+38f;
                constexpr float xj0   = 3.8317059702075123156e+0f;
                constexpr float xj1   = 7.0155866698156187535e+0f;
                constexpr float xy0   = 2.1971413260310170351e+0f;
                constexpr float xy1   = 5.4296810407941351328e+0f;
                constexpr float xj01  = 981.0e+0f;
                constexpr float xj02  = -3.2527979248768438556e-04f;
                constexpr float xj11  = 1796.0e+0f;
                constexpr float xj12  = -3.8330184381246462950e-05f;
                constexpr float xy01  = 562.0e+0f;
                constexpr float xy02  = 1.8288260310170351490e-03f;
                constexpr float xy11  = 1390.0e+0f;
                constexpr float xy12  = -6.4592058648672279948e-06f;
                float           ax,down,prod,resj;
                float           r0,r1,up,w;
                float           wsq,xden,xnum,xy;
                float           z,zsq,result;
                bool            b;
                
                ax = fabsf(arg);
                b  = (arg<half) && (ax*xinf<pi2);
                if(jint==1 && arg <= zero || b)
                   return (-xinf);
                else if(xmax<ax) 
                   return (zero);
                if(eight<ax) 
                   goto L800;
                else if(ax<=xsmall) {
                   if(jint==0)
                      result = arg*half;
                   else
                      result = -pi2/ax;
                   return (result);
                }
                
                zsq = ax*ax;
                if(ax<=four) {
                   xnum = __fmaf_ru(__fmaf_ru(pj0[6],zsq,pj0[5]),zsq,pj0[4]);
                   xden = zsq+qj0[4];
                   xnum = __fmaf_ru(xnum,zsq,pj0[0]);
                   xden = __fmaf_ru(xden,zsq,qj0[0]);
                   xnum = __fmaf_ru(xnum,zsq,pj0[1]);
                   xden = __fmaf_ru(xden,zsq,qj0[1]);
                   xnum = __fmaf_ru(xnum,zsq,pj0[2]);
                   xden = __fmaf_ru(xden,zsq,qj0[2]);
                   xnum = __fmaf_ru(xnum,zsq,pj0[3]);
                   xden = __fmaf_ru(xden,zsq,qj0[3]);
                   prod = arg*((ax-xj01/two56)-xj02)*(ax+xj0);
                }
                else {
                   xnum = pj1[0];
                   xden = __fmaf_ru(zsq+qj1[6],zsq,qj1[0]);
                   xnum = __fmaf_ru(xnum,zsq,pj1[1]);
                   xden = __fmaf_ru(xden,zsq,qj1[1]);
                   xnum = __fmaf_ru(xnum,zsq,pj1[2]);
                   xden = __fmaf_ru(xden,zsq,qj1[2]);
                   xnum = __fmaf_ru(xnum,zsq,pj1[3]);
                   xden = __fmaf_ru(xden,zsq,qj1[3]);
                   xnum = __fmaf_ru(xnum,zsq,pj1[4]);
                   xden = __fmaf_ru(xden,zsq,qj1[4]);
                   xnum = __fmaf_ru(xnum,zsq,pj1[5]);
                   xden = __fmaf_ru(xden,zsq,qj1[5]);
                   xnum = xnum*(ax-eight)*(ax+eight)+pj1[6];
                   xnum = xnum*(ax-four)*(ax+four)+pj1[7];
                   prod = arg*((ax-xj11/two56)-xj12)*(ax+xj1);
                }
                result = prod*(xnum/xden);
                if(jint==0) return (result);
                
                if(ax<=four) {
                   up = (ax-xy01/two56)-xy02;
                   xy = xy0;
                }
                else {
                   up = (ax-xy11/two56)-xy12;
                   xy = xy1;
               }
               down = ax+xy;
               if(fabsf(up)<p17*down) {
                  w   = up/down;
                  wsq = w*w;
                  xnum = plg[0];
                  xden = wsq+qlg[0];
                  xnum = __fmaf_ru(xnum,wsq,plg[1]);
                  xden = __fmaf_ru(xden,wsq,qlg[1]);
                  xnum = __fmaf_ru(xnum,wsq,plg[2]);
                  xden = __fmaf_ru(xden,wsq,qlg[2]);
                  xnum = __fmaf_ru(xnum,wsq,plg[3]);
                  xden = __fmaf_ru(xden,wsq,qlg[3]);
                  resj = pi2*result*w*xnum/xden;
               }
               else {
                  resj = pi2*result*log(ax/xy);
               }
               
               if(ax<=four) {
                  xnum = __fmaf_ru(py0[6],zsq,py0[0]);
                  xden = zsq+qy0[0];
                  xnum = __fmaf_ru(xnum,zsq,py0[1]);
                  xden = __fmaf_ru(xden,zsq,qy0[1]);
                  xnum = __fmaf_ru(xnum,zsq,py0[2]);
                  xden = __fmaf_ru(xden,zsq,qy0[2]);
                  xnum = __fmaf_ru(xnum,zsq,py0[3]);
                  xden = __fmaf_ru(xden,zsq,qy0[3]);
                  xnum = __fmaf_ru(xnum,zsq,py0[4]);
                  xden = __fmaf_ru(xden,zsq,qy0[4]);
                  xnum = __fmaf_ru(xnum,zsq,py0[5]);
                  xden = __fmaf_ru(xden,zsq,qy0[5]);
               }
               else {
                  xnum = __fmaf_ru(py1[8],zsq,py1[0]);
                  xden = zsq+qy1[0];
                  xnum = __fmaf_ru(xnum,zsq,py1[1]);
                  xden = __fmaf_ru(xden,zsq,qy1[1]);
                  xnum = __fmaf_ru(xnum,zsq,py1[2]);
                  xden = __fmaf_ru(xden,zsq,qy1[2]);
                  xnum = __fmaf_ru(xnum,zsq,py1[3]);
                  xden = __fmaf_ru(xden,zsq,qy1[3]);
                  xnum = __fmaf_ru(xnum,zsq,py1[4]);
                  xden = __fmaf_ru(xden,zsq,qy1[4]);
                  xnum = __fmaf_ru(xnum,zsq,py1[5]);
                  xden = __fmaf_ru(xden,zsq,qy1[5]);
                  xnum = __fmaf_ru(xnum,zsq,py1[6]);
                  xden = __fmaf_ru(xden,zsq,qy1[6]);
                  xnum = __fmaf_ru(xnum,zsq,py1[7]);
                  xden = __fmaf_ru(xden,zsq,qy1[7]);
               }
               result = resj+(up*down/ax)*xnum/xden;
               return (result);
L800:
               z = eight/ax;
               w = (int)(ax/twopi)+throv8;
               w = (ax-w*twopi1)-w*twopi2;
               zsq = z*z;
               xnum = p0[5];
               xden = zsq+q0[5];
               up = p1[5];
               down = zsq+q1[5];
               r0   = xnum/xden;
               r1   = up/down;
               if(jint==0) {
                  result = (rtpi2/sqrtf(ax)) *
                           (r0*cosf(w)-z*r1*sinf(w));
                  return (result);
               }
               else {
                  result = (rtpi2/sqrtf(ax)) *
                           (r0*sinf(w)+z*r1*cosf(w));
                  return (result);
               }
               if(jint==0 && arg<zero) {
                 result = -result;
                 return (result);
              }

       }
       
       
/*
    
     !*****************************************************************************80
!
!! DLGAMA evaluates log ( Gamma ( X ) ) for a real argument.
!
!  Discussion:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2.  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0 are similar to those in
!    reference 1, but are unpublished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument of the function.
!
!    Output, real(kind=sp) ::  DLGAMA, the value of the function.
!

*/


        __device__ float dlgama(const float x) {
        
               const float p1[8] = {
                         4.945235359296727046734888f,2.018112620856775083915565e+2f, 
                         2.290838373831346393026739e+3f,1.131967205903380828685045e+4f, 
                         2.855724635671635335736389e+4f,3.848496228443793359990269e+4f, 
                         2.637748787624195437963534e+4f,7.225813979700288197698961e+3f};
               const float q1[8] = {
                         6.748212550303777196073036e+1f,1.113332393857199323513008e+3f, 
                         7.738757056935398733233834e+3f,2.763987074403340708898585e+4f, 
                         5.499310206226157329794414e+4f,6.161122180066002127833352e+4f, 
                         3.635127591501940507276287e+4f,8.785536302431013170870835e+3f};
               const float p2[8] = {
                         4.974607845568932035012064f,5.424138599891070494101986e+2f, 
                         1.550693864978364947665077e+4f,1.847932904445632425417223e+5f, 
                         1.088204769468828767498470e+6f,3.338152967987029735917223e+6f,
                         5.106661678927352456275255e+6f,3.074109054850539556250927e+6f};
               const float q2[8] = {
                         1.830328399370592604055942e+2f,7.765049321445005871323047e+3f, 
                         1.331903827966074194402448e+5f,1.136705821321969608938755e+6f, 
                         5.267964117437946917577538e+6f,1.346701454311101692290052e+7f, 
                         1.782736530353274213975932e+7f,9.533095591844353613395747e+6f};
               const float p4[8] = {
                         1.474502166059939948905062e+4f,2.426813369486704502836312e+6f, 
                         1.214755574045093227939592e+8f,2.663432449630976949898078e+9f, 
                         2.940378956634553899906876e+10f,1.702665737765398868392998e+11f, 
                         4.926125793377430887588120e+11f,5.606251856223951465078242e+11f};
               const float q4[8] = {
                         2.690530175870899333379843e+3f,6.393885654300092398984238e+5f, 
                         4.135599930241388052042842e+7f,1.120872109616147941376570e+9f, 
                         1.488613728678813811542398e+10f,1.016803586272438228077304e+11f, 
                         3.417476345507377132798597e+11f,4.463158187419713286462081e+11f};
               const float c[7]  = {
                        -1.910444077728e-03f,8.4171387781295e-04f, 
                        -5.952379913043012e-04f,7.93650793500350248e-04f,
                        -2.777777777777681622553e-03f,8.333333333333333331554247e-02f,
                         5.7083835261e-03f};
               constexpr float one    = 1.0
               constexpr float half   = 0.5
               constexpr float twelve = 12.0
               constexpr float zero   = 0.0
               constexpr float four   = 4.0
               constexpr float thrhal = 1.5
               constexpr float two    = 2.0
               constexpr float pnt68  = 0.6796875f
               constexpr float sqrtpi = 0.9189385332046727417803297f
               constexpr float eps    = 1.19e-07f
               constexpr float frtbig = 3.4028234664e+38f
               constexpr float d1     = -5.772156649015328605195174e-1f
               constexpr float d2     = 4.227843350984671393993777e-1f
               constexpr float d4     = 1.791759469228055000094023f
               float     corr,xden,xnum;
               float     xm1,xm2,xm4;
               float     y,ysq,res;
               
               y  = x;
               if(zero<y) {
                  if(y<=eps) {
                     res = -logf(y);
                  }
                  else if(y<=thrhal) {
                     if(y<pnt68) {
                        corr = -logf(y);
                        xm1  = y;
                     }
                     else {
                        corr = zero;
                        xm1  = (y-half)-half;
                     }
                     
                     if(y<=half || pnt68<=y) {
                        xden = one;
                        xnum = zero;
                        xnum = __fmaf_ru(xnum,xm1,p1[0]);
                        xden = __fmaf_ru(xden,xm1,q1[0]);
                        xnum = __fmaf_ru(xnum,xm1,p1[1]);
                        xden = __fmaf_ru(xden,xm1,q1[1]);
                        xnum = __fmaf_ru(xnum,xm1,p1[2]);
                        xden = __fmaf_ru(xden,xm1,q1[2]);
                        xnum = __fmaf_ru(xnum,xm1,p1[3]);
                        xden = __fmaf_ru(xden,xm1,q1[3]);
                        xnum = __fmaf_ru(xnum,xm1,p1[4]);
                        xden = __fmaf_ru(xden,xm1,q1[4]);
                        xnum = __fmaf_ru(xnum,xm1,p1[5]);
                        xden = __fmaf_ru(xden,xm1,q1[5]);
                        xnum = __fmaf_ru(xnum,xm1,p1[6]);
                        xden = __fmaf_ru(xden,xm1,q1[6]);
                        xnum = __fmaf_ru(xnum,xm1,p1[7]);
                        xden = __fmaf_ru(xden,xm1,q1[7]);
                        res = corr+(xm1*(d1+xm1*(xnum/xden)));
                     }
                     else {
                        xm2  = (y-half)-half;
                        xden = one;
                        xnum = zero;
                        xnum = __fmaf_ru(xnum,xm2,p2[0]);
                        xden = __fmaf_ru(xden,xm2,q2[0]);
                        xnum = __fmaf_ru(xnum,xm2,p2[1]);
                        xden = __fmaf_ru(xden,xm2,q2[1]);
                        xnum = __fmaf_ru(xnum,xm2,p2[2]);
                        xden = __fmaf_ru(xden,xm2,q2[2]);
                        xnum = __fmaf_ru(xnum,xm2,p2[3]);
                        xden = __fmaf_ru(xden,xm2,q2[3]);
                        xnum = __fmaf_ru(xnum,xm2,p2[4]);
                        xden = __fmaf_ru(xden,xm2,q2[4]);
                        xnum = __fmaf_ru(xnum,xm2,p2[5]);
                        xden = __fmaf_ru(xden,xm2,q2[5]);
                        xnum = __fmaf_ru(xnum,xm2,p2[6]);
                        xden = __fmaf_ru(xden,xm2,q2[6]);
                        xnum = __fmaf_ru(xnum,xm2,p2[7]);
                        xden = __fmaf_ru(xden,xm2,q2[7]);
                        res = corr+xm2*(d2+xm2*(xnum/xden));
                     }
                  }
                  else if(y<=four) {
                   
                        xm2  = y-two;
                        xden = one;
                        xnum = zero;
                        xnum = __fmaf_ru(xnum,xm2,p2[0]);
                        xden = __fmaf_ru(xden,xm2,q2[0]);
                        xnum = __fmaf_ru(xnum,xm2,p2[1]);
                        xden = __fmaf_ru(xden,xm2,q2[1]);
                        xnum = __fmaf_ru(xnum,xm2,p2[2]);
                        xden = __fmaf_ru(xden,xm2,q2[2]);
                        xnum = __fmaf_ru(xnum,xm2,p2[3]);
                        xden = __fmaf_ru(xden,xm2,q2[3]);
                        xnum = __fmaf_ru(xnum,xm2,p2[4]);
                        xden = __fmaf_ru(xden,xm2,q2[4]);
                        xnum = __fmaf_ru(xnum,xm2,p2[5]);
                        xden = __fmaf_ru(xden,xm2,q2[5]);
                        xnum = __fmaf_ru(xnum,xm2,p2[6]);
                        xden = __fmaf_ru(xden,xm2,q2[6]);
                        xnum = __fmaf_ru(xnum,xm2,p2[7]);
                        xden = __fmaf_ru(xden,xm2,q2[7]);
                        res  = xm2*(d2+xm2*(xnum/xden));
               }   
               else if(y<=twelve) {
                        
                        xm4 = y-four;
                        xden = -one;
                        xnum = zero;
                        xnum = __fmaf_ru(xnum,xm4,p2[0]);
                        xden = __fmaf_ru(xden,xm4,q2[0]);
                        xnum = __fmaf_ru(xnum,xm4,p2[1]);
                        xden = __fmaf_ru(xden,xm4,q2[1]);
                        xnum = __fmaf_ru(xnum,xm4,p2[2]);
                        xden = __fmaf_ru(xden,xm4,q2[2]);
                        xnum = __fmaf_ru(xnum,xm4,p2[3]);
                        xden = __fmaf_ru(xden,xm4,q2[3]);
                        xnum = __fmaf_ru(xnum,xm4,p2[4]);
                        xden = __fmaf_ru(xden,xm4,q2[4]);
                        xnum = __fmaf_ru(xnum,xm4,p2[5]);
                        xden = __fmaf_ru(xden,xm4,q2[5]);
                        xnum = __fmaf_ru(xnum,xm4,p2[6]);
                        xden = __fmaf_ru(xden,xm4,q2[6]);
                        xnum = __fmaf_ru(xnum,xm4,p2[7]);
                        xden = __fmaf_ru(xden,xm4,q2[7]);
                        res  = d4+xm4*(xnum/xden);
               }  
               else {
                        res = zero;
                        if(y<=frtbig) {
                           res = c[6];
                           ysq = y*y;
                           res = res/ysq+c[0];
                           res = res/ysq+c[1];
                           res = res/ysq+c[2];
                           res = res/ysq+c[3];
                           res = res/ysq+c[4];
                           res = res/ysq+c[5];
                        }
                        res /= y;
                        corr = logf(y);
                        res  = res+sqrtpi-half*corr;
                        res  = res+y*(corr-one);
               }  
               
            }
            else {
                   
                   res = x;
            }
            
            return (res);
        }
        
        
/*

      !*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This function was originally named DGAMMA.
!
!    However, a number of Fortran compilers now include a library
!    function of this name.  To avoid conflicts, this function was
!    renamed R8_GAMMA.
!
!    This routine calculates the GAMMA function for a real argument X.
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument of the function.
!
!    Output, real(kind=sp) ::  R8_GAMMA, the value of the function.

*/
          
        __device__ float r4_gamma(const float x) {
        
                const float c[7] = {
                     -1.910444077728e-03f,
                      8.4171387781295e-04f,
                     -5.952379913043012e-04f,
                      7.93650793500350248e-04f,
                     -2.777777777777681622553e-03f,
                      8.333333333333333331554247e-02f,
                      5.7083835261e-03f};
                const float p[8] = {
                      -1.71618513886549492533811e+0f,   2.47656508055759199108314e+01f,
                      -3.79804256470945635097577e+02f,  6.29331155312818442661052e+02f,
                       8.66966202790413211295064e+02f, -3.14512729688483675254357e+04f,
                      -3.61444134186911729807069e+04f,  6.64561438202405440627855e+04f};
                const float q[8] = {
                      -3.08402300119738975254353e+01f,  3.15350626979604161529144e+02f, 
                      -1.01515636749021914166146e+03f, -3.10777167157231109440444e+03f, 
                       2.25381184209801510330112e+04f,  4.75584627752788110767815e+03f, 
                      -1.34659959864969306392456e+05f, -1.15132259675553483497211e+05f};
                 constexpr float one    = 1.0e+00f;
                 constexpr float half   = 0.5e+00f;
                 constexpr float twelve = 12.0e+00f;
                 constexpr float two    = 2.0e+00f;
                 constexpr float zero   = 0.0e+00f;
                 constexpr float sqrtpi = 0.9189385332046727417803297e+00f;
                 constexpr float pi     = 3.1415926535897932384626434e+00f;
                 constexpr float xbig   = 171.624e+00f;
                 constexpr float xminin = 1.1754943508e-38f;
                 constexpr float eps    = 1.19e-07f;
                 float           fact,res,sum;
                 float           xden,xnum,y,y1;
                 float           ysq,z;
                 int             n;
                 bool            parity;
                 
                 parity = false;
                 fact = one;
                 n = 0;
                 y = x;
                 
                 if(y<=zero) {
                    y   = -x;
                    y1  = truncf(y);
                    res = y-y1;
                    if(res!=zero) {
                       if(y1!=truncf(y1*half)*two) parity = true;
                       fact = -pi/sinf(pi*res);
                       y    =+ one;
                    }
                    else {
                       res = xminin;
                       return (res);
                    }
                 }
                 
                 if(y<eps) {
                    if(xminin<=y)
                       res = one/y;
                    else {
                       res = xminin;
                       return (res);
                    }
                 }
                 else if(y<twelve) {
                    y1 = y;
                    if(y<one) {
                       z = y;
                       y += one;
                    }
                    else {
                       n = (int)y-1;
                       y = y-(float)n;
                       z = y-one;
                    }
                    xnum = zero;
                    xden = one;
                    xnum = (xnum+p[0])*z;
                    xden = __fmaf_ru(xden,z,q[0]);
                    xnum = (xnum+p[1])*z;
                    xden = __fmaf_ru(xden,z,q[1]);
                    xnum = (xnum+p[2])*z;
                    xden = __fmaf_ru(xden,z,q[2]);
                    xnum = (xnum+p[3])*z;
                    xden = __fmaf_ru(xden,z,q[3]);
                    xnum = (xnum+p[4])*z;
                    xden = __fmaf_ru(xden,z,q[4]);
                    xnum = (xnum+p[5])*z;
                    xden = __fmaf_ru(xden,z,q[5]);
                    xnum = (xnum+p[6])*z;
                    xden = __fmaf_ru(xden,z,q[6]);
                    xnum = (xnum+p[7])*z;
                    xden = __fmaf_ru(xden,z,q[7]);
                    res  = xnum/xden+one;
                    if(y1<y) {
                       res /= y1;
                    }
                    else if(y<y1) {
                       for(int i=1; i!=n; ++i) {
                           res *= y;
                           y   += one;
                       }
                    }
                    else {
                       if(x<=xbig) {
                           ysq = y*y;
                           sum = c[6];
                           sum = sum/ysq+c[0];
                           sum = sum/ysq+c[1];
                           sum = sum/ysq+c[2];
                           sum = sum/ysq+c[3];
                           sum = sum/ysq+c[4];
                           sum = sum/ysq+c[5];
                           sum = sum/y-y+sqrtpi;
                           sum = sum+(y-half)*logf(y);
                           res = expf(sum);
                       }
                       else {
                           res = 3.4028234664e+38f;
                           return (res);
                       }
                    }
                 }
                 
                if(parity)    res = -res;
                if(fact!=one) res /= fact;
                return (res);
        }       
        
        

/*

      !*****************************************************************************80
!
!! R8_PSI evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
!             = d/dX LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  XX, the argument of the function.
!
!    Output, real(kind=sp) ::  R8_PSI, the value of the function.
!

*/ 


        __device__ float r4_psi(const float xx) {
        
                const float p1[9] = {
                      4.5104681245762934160e-03f,
                      5.4932855833000385356f, 
                      3.7646693175929276856e+02f, 
                      7.9525490849151998065e+03f, 
                      7.1451595818951933210e+04f, 
                      3.0655976301987365674e+05f, 
                      6.3606997788964458797e+05f, 
                      5.8041312783537569993e+05f, 
                      1.6585695029761022321e+05f};
                const float p2[7]  = {
                     -2.7103228277757834192f,
                     -1.5166271776896121383e+01f, 
                     -1.9784554148719218667e+01f,
                     -8.8100958828312219821f, 
                     -1.4479614616899842986f,
                     -7.3689600332394549911e-02f,
                     -6.5135387732718171306e-21f};
                const float q1[8]  = {
                      9.6141654774222358525e+01f, 
                      2.6287715790581193330e+03f, 
                      2.9862497022250277920e+04f, 
                      1.6206566091533671639e+05f,
                      4.3487880712768329037e+05f,
                      5.4256384537269993733e+05f,
                      2.4242185002017985252e+05f, 
                      6.4155223783576225996e-08f};
                const float q2[6]  = {
                      4.4992760373789365846e+01f,
                      2.0240955312679931159e+02f,
                      2.4736979003315290057e+02f,
                      1.0742543875702278326e+02f,
                      1.7463965060678569906e+01f, 
                      8.8427520398873480342e-01f};
                constexpr float  four  = 4.0f;
                constexpr float  fourth= 0.25f;
                constexpr float  half  = 0.5f;
                constexpr float  one   = 1.0f;
                constexpr float  piov4 = 0.78539816339744830962f;
                constexpr float  three = 3.0f;
                constexpr float  zero  = 0.0f;
                constexpr float  x01   = 187.0f;
                constexpr float  x01d  = 128.0f;
                constexpr float  x02   = 6.9464496836234126266e-04f;
                constexpr float  xlarge= 2.04e+15f;
                constexpr float  xmax1 = 3.60e+16f;
                constexpr float  xsmall= 2.05e-09f;
                float            sgn,upper,w,x;
                float            z,aug,den,result;
                int              n,nq,it;
                
                x  = xx;
                w  = fabsf(x);
                aug= zero;
                if(zero<x) return (xx);
                
                if(x<half) {
                
                   if(w<=xsmall) {
                      aug = -one/w;
                   }
                   else {
                      if(x<zero)
                         sgn = piov4;
                      else
                         sgn = -piov4;
                      it = (int)w;
                      w  = w-(float)it;
                      nq = (int)(w*four);
                      w  = four*(w-((float)nq)*fourth);
                      n  = nq/2;
                      if(n+n!=nq) w = one-w;
                      z  = piov4*w;
                      if((n%2)!=0) sgn = -sgn;
                      n  = (nq+1)/2;
                      if((n%2)==0) 
                         aug = sgn*(four/tanf(z));
                      else
                         aug = sgn*four*tanf(z);
                   }
                }
                x = one-x;
                
                if(x<=three) {
                   den   = x;
                   upper = p1[0]*x;
                   den   = (den+q1[0])*x;
                   upper = (upper+p1[1])*x;
                   den   = (den+q1[1])*x;
                   upper = (upper+p1[2])*x;
                   den   = (den+q1[2])*x;
                   upper = (upper+p1[3])*x;
                   den   = (den+q1[4])*x;
                   upper = (upper+p1[5])*x;
                   den   = (den+q1[5])*x;
                   upper = (upper+p1[6])*x;
                   den   = (den+q1[6])*x;
                   upper = (upper+p1[7])*x;
                   den   = (upper+p1[8])/(den+q1[7]);
                   x     = (x-x01/x01d)-x02;
                   result= __fmaf_ru(den,x,aug);
                   return (result);
                }
                
                if(x<xlarge) {
                   w     = one/(x*x);
                   den   = w;
                   upper = p2[0]*w;
                   den   = (den+q2[0])*w;
                   upper = (upper+p2[1])*w;
                   den   = (den+q2[1])*w;
                   upper = (upper+p2[2])*w;
                   den   = (den+q2[2])*w;
                   upper = (upper+p2[3])*w;
                   den   = (den+q2[3])*w;
                   upper = (upper+p2[4])*w;
                   den   = (den+q2[4])*w;
                   upper = (upper+p2[5])*w;
                   aug   = (upper+p2[6])/(den+q2[5])-half/x+aug; 
                }
                
                result = aug+logf(x);
                return (result);
        }
        
        
/*
    
     !*****************************************************************************80
!
!! VVSA computes parabolic cylinder function V(nu,x) for small arguments.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Input, real(kind=sp) ::  VA, the order nu.
!
!    Output, real(kind=sp) ::  PV, the value of V(nu,x).
!

*/


        __device__ float vvsa(const float x,
                              const float va) {
                              
               float a0,ep,fac,g1;
               float ga0,gm,gw,pv;
               float r,r1,sq2,sv;
               float sv0,v1,va0,vb0;
               float vm,pv;
               int   m;
               
               constexpr float eps = 1.0e-15f;
               constexpr float pi  = 3.14159265358979323846264f;
               ep                  = expf(-0.25f*x*x);
               va0                 = 1.0f+0.5f*va;
               if(x==0.0f) {
                  if((va0<=0.0f && va0==truncf(va0)) || va==0.0f) {
                     pv = 0.0f;
                  }
                  else {
                     vb0 = -0.5f*va;
                     sv0 = sinf(va0*pi);
                     ga0 = gamma(va0);
                     pv  = powf(2.0f,vb0)*sv0/ga0;
                  }
               }
               else {
                   for(m=1; m!=250; ++m) {
                       float tm = (float)m;
                       vm       = 0.5f*(tm-va);
                       gm       = gamma(vm);
                       r        = r*sq2*x/tm;
                       fac      = -fac;
                       gw       = fac*sv+1.0f;
                       r1       = gw*r*gm;
                       pv       = pv+r1;
                       if(fabsf(r1/pv)<eps && gw!=0.0f) break;
                   }
                   pv = a0*pv;
               }
               
               return (pv);
      }
      
      
/*

      !*****************************************************************************80
!
!! VVLA computes parabolic cylinder function Vv(x) for large arguments.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Input, real(kind=sp) ::  VA, the order nu.
!
!    Output, real(kind=sp) ::  PV, the value of V(nu,x).
!
 
*/


        __device__ float vvla(const float x,
                              const float va) {
                              
               float a0,dsl,gl,pdl;
               float qe,r,x1,pv;
               int   k;
               constexpr float pi   = 3.14159265358979323846264f;
               constexpr float eps  = 1.0e-12f;
               if(x<0.0f) return (x);
               qe                   = expf(0.25f*x*x);
               a0                   = powf(fabsf(x),-va-1.0f)*2.50662827463100050241577f*qe;
               r                    = 1.0f;
               pv                   = 1.0f;
               for(k=1; k!=18; ++k) {
                   float tk = (float)k;
                   float t0 = 2.0f*tk+va;
                   r        = 0.5f*r*t0-1.0f*t0/(tk*x*x);
                   pv       += r;
                   if(fabsf(r/pv)<eps) break;
               }  
               pv = a0*pv;
               return (pv);               
      }
      
      
/*
   
    !*****************************************************************************80
!
!! RCTY computes Riccati-Bessel function of the second kind, and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order of yn(x).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  RY(0:N), the values of x yn(x).
!
!    Output, real(kind=sp) ::  DY(0:N), the values of [x yn(x)]'.
!

*/


        __device__ void rcty(const int n,
                             const float x,
                             int        &nm,
                             float * __restrict__ ry,
                             float * __restrict__ dy) {
                             
              
              float  rf0,rf1,rf2;
              int    k;
              
              nm    = n;
              ry[0] = -cosf(x);
              ry[1] = ry[0]/x-sinf(x);
              rf0   = ry[0];
              rf1   = ry[1];
              for(k=2; k!=n; ++k) {
                  float tk = (float)k;
                  float t0 = 2.0f*tk-1.0f;
                  rf2      = tk*rf1/x-rf0;
                  if(3.4028234664e+38f<fabsf(rf2)) break;
                  ry[k]    = rf2;
                  rf0      = rf1;
                  rf1      = rf2;
              }     
              nm    = k-1;
              dy[0] = sinf(x);
              for(k=1; k!=nm; ++k) {
                  float tk = (float)k;
                  dy[k]    = -tk*ry[k]/x+ry[k-1];
              }           
      }
      
      
/*

    !*****************************************************************************80
!
!! RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order of jn(x).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  RJ(0:N), the values of x jn(x).
!
!    Output, real(kind=sp) ::  DJ(0:N), the values of [x jn(x)]'.
 
*/


        __device__ void rctj(const int n,
                             const float x,
                             int        &nm,
                             float * __restrict__ rj,
                             float * __restrict__ dj) {
                             
              float cs,f,f0,f1;
              float rj0,rj1;
              int   k,m;
              
              nm    = n;
              rj[0] = sinf(x);
              rj[1] = rj[0]/x-cosf(x);
              rj0 = rj[0];
              rj1 = rj[1];
              
              if(2<=n) {
                 m = msta1(x,200);
                 if(m<n)
                    nm = m;
                 else
                    m  = msta2(x,n,15);
                 f0 = 0.0f;
                 f1 = 1.1754943508e-38f;
                 for(k=m; k!=0; --k) {
                     float tk = (float)k;
                     float t0 = __fmaf_ru(2.0f,tk,3.0f)*f1/x-f0;
                     if(k<=nm) rj[k] = f;
                     f0 = f1;
                     f1 = f;
                 }
                 
                 if(fabsf(rj1)<fabsf(rj0))
                    cs = rj0/f;
                 else
                    cs = rj1/f0;
                 for(k=0; k!=nm; ++k) 
                     rj[k] *= cs;
             }
             
             dj[0] = cosf(x);
             for(k=1; k!=nm; ++k) {
                 float tk = (float)k;
                 dj[k]    = -tk*rj[k]/x+rj[k-1];
             }                    
      }
      
      
/*

    !*****************************************************************************80
!
!! OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
!
!  Discussion:
!
!    This procedure computes orthogonal polynomials: Tn(x) or Un(x),
!    or Ln(x) or Hn(x), and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: KT, the function code:
!    1 for Chebyshev polynomial Tn(x)
!    2 for Chebyshev polynomial Un(x)
!    3 for Laguerre polynomial Ln(x)
!    4 for Hermite polynomial Hn(x)
!
!    Input, integer(kind=i4) :: N, the order.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  PL(0:N), DPL(0:N), the value and derivative of
!    the polynomials of order 0 through N at X. 
  
*/


        __device__ void othpol(const int            n,
                               const float          x,
                               float * __restrict__ pl,
                               float * __restrict__ dpl) {
                               
               float a,b,c;
               float dy0,dy1,dyn;
               float y0,y1,yn;
               int   k,kf;
               
               a      = 2.0f
               b      = 0.0f;
               c      = 1.0f;
               y0     = 1.0f;
               y1     = 2.0f*x;
               dy0    = 0.0f;
               dy1    = 2.0f;
               pl[0]  = 1.0f;
               pl[1]  = 2.0f*x;
               dpl[0] = 0.0f;
               dpl[1] = 2.0f;
               
               if(kf==1) {
                  y1    = x;
                  dy1   = 1.0f;
                  pl[1] = x;
                  dpl[1]= 1.0f; 
               }  
               else if(kf==3) {
                  y1     = 1.0f-x;
                  dy1    = -1.0f;
                  pl[1]  = 1.0f-x;
                  dpl[1] = -1.0f;
               } 
               
               if(kf==3) {
                  for(k=2; k!=n; ++k) {
                      float tk = (float)k;
                      a        = -1.0f/tk;
                      b        = 2.0f+a;
                      c        = 1.0f+a;
                      yn       = __fmaf_ru(a,x,b)*y1-c*y0;
                      dyn      = __fmaf_ru(a,y1,__fmaf_ru(a,x,b))*dy1-c*dy0;
                      pl[k]    = yn;
                      dpl[k]   = dyn;
                      y0       = y1;
                      y1       = yn;
                      dy0      = dy1;
                      dy1      = dyn;
                  }
               } 
               else if(kf==4) {
                      for(k=2; k!=n; ++k) {
                      float tk = (float)k;
                      c        = 2.0f*(tk-1.0f);
                      yn       = __fmaf_ru(a,x,b)*y1-c*y0;
                      dyn      = __fmaf_ru(a,y1,__fmaf_ru(a,x,b))*dy1-c*dy0;
                      pl[k]    = yn;
                      dpl[k]   = dyn;
                      y0       = y1;
                      y1       = yn;
                      dy0      = dy1;
                      dy1      = dyn;
                  } 
               }                
      }
      
      
/*
     
     !*****************************************************************************80
!
!! LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the maximum degree.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  PN(0:N), PD(0:N), the values and derivatives
!    of the polyomials of degrees 0 to N at X.
!

*/


        __device__ void lpn(const int   n,
                            const float x,
                            float * __restrict__ pn,
                            float * __restrict__ pd) {
                            
                float p0,p1,pf;
                int   k;
                
                pn[0] = 1.0f;
                pn[1] = x;
                pd[0] = 0.0f;
                pd[1] = 1.0f;
                p0    = 1.0f;
                p1    = x;
                
                if(fabsf(x)==1.0f) {
                   for(k=2; k!=n; ++k) {
                       float tk = (float)k;
                       float t0 = 2.0f*tk-1.0f;
                       float t1 = tk-1.0f;
                       float t2 = tk+1.0f;
                       pf       = t0/(tk*x*p1)-t1/tk*p0;
                       pd[k]    = 0.5f*powf(x,t2)*tk*t2;
                       p0       = p1;
                       p1       = pf;
                   }
                } 
                else {
                    for(k=2; k!=n; ++k) {
                       float tk = (float)k;
                       float t0 = 2.0f*tk-1.0f;
                       float t1 = tk-1.0f;
                       pf       = t0/(tk*x*p1)-t1/tk*p0;
                       pd[k]    = tk*(p1-x*pf)/(1.0f*x*x);
                       p0       = p1;
                       p1       = pf;
                   }  
                }                   
      }
      
      
/*
    
       !*****************************************************************************80
!
!! LGAMA computes the gamma function or its logarithm.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: KF, the argument code.
!    1, for gamma(x);
!    2, for ln(gamma(x)).
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  GL, the function value.
!

*/


        __device__ float lgama(const int   kf,
                               const float x) {
                               
              const float a[10] = {
                      8.333333333333333e-02f, 
                     -2.777777777777778e-03f, 
                      7.936507936507937e-04f, 
                     -5.952380952380952e-04f, 
                      8.417508417508418e-04f,
                     -1.917526917526918e-03f,
                      6.410256410256410e-03f,
                     -2.955065359477124e-02f, 
                      1.796443723688307e-01f, 
                     -1.39243221690590f}; 
             float gl,gl0,x0;
             float x2,xp,t0,t1;
             int   k,n;
             
             x0  = x;
             if(x==1.0f || x==2.0f) {
                gl = 0.0f;
                if(kf==1) gl = 1.0f;
                return (gl);
             }
             else if(x<=7.0f) {
                n   = (int)(7.0f-x);
                x0  = x+(float)n;
             }
             
             x2 = 1.0f/(x0*x0);
             xp = 6.283185307179586477f;
             gl0= a[9];
             gl0= __fmaf_ru(gl0,x2,a[8]);
             gl0= __fmaf_ru(gl0,x2,a[7]);
             gl0= __fmaf_ru(gl0,x2,a[6]);
             gl0= __fmaf_ru(gl0,x2,a[5]);
             gl0= __fmaf_ru(gl0,x2,a[4]);
             gl0= __fmaf_ru(gl0,x2,a[3]);
             gl0= __fmaf_ru(gl0,x2,a[2]);
             gl0= __fmaf_ru(gl0,x2,a[1]);
             gl0= __fmaf_ru(gl0,x2,a[0]);
             t0 = (x0+0.5f)*logf(x0)-x0;
             gl = gl0/x0+0.5f*logf(xp)+t0;
             if(x<=7.0f) {
                for(k=1; k!=n; ++k) {
                    gl =  gl-logf(x-1.0f);
                    x0 -= 1.0f;
                }
             }
             if(kf==1) gl = expf(gl);
             return (gl);
             
      }
      
      
/*
    
       
!*****************************************************************************80
!
!! KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
!    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
!

*/


        __device__ void klvnb(const float x,
                              float      &ber,
                              float      &bei,
                              float      &ger,
                              float      &gei,
                              float      &der,
                              float      &dei,
                              float      &her,
                              float      &hei) {
                              
                float csn,csp,fxi,fxr;
                float pni,pnr,ppi,ppr;
                float ssn,ssp,t,t2;
                float tni,tnr,tpi,tpr;
                float u,v,yc1,yc2;
                float yci,ye1,ye2,yei,yd;
                int   l;
                constexpr float pi = 3.141592653589793f;
                
                if(x==0.0f) return;
                if(x<8.0f) {
                   t = x/8.0f;
                   t2 = t*t;
                   u = t2*t2;
                   ber = ((((((
                       - 0.901e-05f*u 
                       + 0.122552e-02f)*u 
                       - 0.08349609f)*u 
                       + 2.64191397f)*u 
                       - 32.36345652f)*u 
                       + 113.77777774f)*u 
                       - 64.0f)*u 
                       + 1.0_sp;

                   bei = t*t*((((((
                         0.11346e-03f*u 
                       - 0.01103667f)*u
                       + 0.52185615f)*u
                       - 10.56765779f)*u
                       + 72.81777742f)*u
                       - 113.77777774f)*u
                       + 16.0f);

                   ger = (((((( 
                       - 0.2458e-04f*u
                       + 0.309699e-02f)*u
                       - 0.19636347f)*u
                       + 5.65539121f)*u
                       - 60.60977451f)*u
                       + 171.36272133f)*u
                       - 59.05819744f)*u
                       - 0.57721566f;
                   ger = __fmaf_ru((ger-logf(0.5f*x))*ber,(0.25f*pi*bei));
                   
                   gei = t2*((((((
                         0.29532e-03f*u
                       - 0.02695875f)*u
                       + 1.17509064f)*u
                       - 21.30060904f)*u
                       + 124.2356965f)*u
                       - 142.91827687f)*u
                       + 6.76454936f);
                   gei = gei-logf(0.5f*x)*bei-0.25f*pi*ber;
                   
                   der = x*t2*((((((
                       - 0.394e-05f*u
                       + 0.45957e-03f)*u
                       - 0.02609253f)*u
                       + 0.66047849f)*u
                       - 6.0681481f)*u
                       + 14.22222222f)*u
                       - 4.0f);

                   dei = x * ((((((
                         0.4609e-04f*u
                       - 0.379386e-02f)*u
                       + 0.14677204f)*u
                       - 2.31167514f)*u
                       + 11.37777772f)*u 
                       - 10.66666666f)*u
                       + 0.5f); 

                   her = x*t2*((((((
                       - 0.1075e-04f*u
                       + 0.116137e-02f)*u
                       - 0.06136358f)*u
                       + 1.4138478f)*u
                       - 11.36433272f)*u 
                       + 21.42034017f)*u 
                       - 3.69113734f);
                   her = her-logf(0.5f*x)*der-ber/x +
                        0.25f*pi*dei;
                        
                   hei = x*((((((
                         0.11997e-03f*u
                       - 0.926707e-02f)*u
                       + 0.33049424f)*u
                       - 4.65950823f)*u
                       + 19.41182758f)*u
                       - 13.39858846f)*u
                       + 0.21139217f);

                   hei = hei-logf(0.5f*x)*dei-bei/x - 
                         0.25f*pi*der;
                }
                else {
                   
                   t = 8.0f/x;
                   for(l=1; l!=2; ++l) {
                       float tl = (float)l;
                       v        = powf(-1.0f,tl)*t;
                       tpr      = (((( 
                                  0.6e-06f*v
                                - 0.34e-05f)*v 
                                - 0.252e-04f)*v
                                - 0.906e-04f)*v*v
                                + 0.0110486f)*v

                       tpi      = ((((
                                  0.19e-05f*v
                                + 0.51e-05f)*v*v
                                - 0.901e-04f)*v
                                - 0.9765e-03f)*v
                                - 0.0110485f)*v
                                - 0.3926991f;
                       if(l==1) {
                          tnr = tpr;
                          tni = tpi;
                       }
                   }
                   
                   yd  = x/1.41421356237309504880169f;
                   ye1 = expf(yd+tpr);
                   ye2 = expf(-yd+tnr);
                   yc1 = 1.0f/sqrtf(2.0f*pi*x);
                   yc2 = sqrtf(pi/(2.0f*x));
                   csp = cosf(yd+tpi);
                   ssp = sinf(yd+tpi);
                   csn = cosf(-yd+tni);
                   ssn = sinf(-yd+tni);
                   ger = yc2*ye2*csn;
                   gei = yc2*ye2*ssn;
                   fxr = yc1*ye1*csp;
                   fxi = yc1*ye1*ssp;
                   ber = fxr-gei*0.31830988618379067153777f;
                   bei = fxi+ger*0.31830988618379067153777f;
                   
                   for(l=1; l!=2; ++l) {
                       float tl  = (float)l;
                       v         = powf(-1.0f,tl)*t;
                       ppr       = (((((
                                   0.16e-05f*v
                                 + 0.117e-04f)*v
                                 + 0.346e-04f)*v 
                                 + 0.5e-06f)*v
                                 - 0.13813e-02f)*v
                                 - 0.0625001f)*v
                                 + 0.7071068f;

                       ppi       = (((((
                                 - 0.32e-05f*v
                                 - 0.24e-05f)*v
                                 + 0.338e-04f)*v
                                 + 0.2452e-03f)*v
                                 + 0.13811e-02f)*v
                                 - 0.1e-06f)*v
                                 + 0.7071068f;
                        if(l==1) {
                          pnr = ppr;
                          pni = ppi;
                       }
                   }
                   
                   her = gei*pni-ger*pnr;
                   hei = -__fmaf_ru(gei,pnr,ger*pni);
                   der = fxr*ppr-fxi*ppi-hei*0.31830988618379067153777f;
                   dei = __fmaf_ru(__fmaf_ru(fxi,ppr,fxr),ppi,her*0.31830988618379067153777f);
                   
                }
                   
                
                                       
      }
      
      
/*

      !*****************************************************************************80
!
!! ITSL0 integrates the Struve function L0(t) from 0 to x.
!
!  Discussion:
!
!    This procedure evaluates the integral of modified Struve function
!    L0(t) with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  TL0, the integral of L0(t) from 0 to x.
! 

*/


       __device__ float itsl0(const float x) {
       
               float a[19];
               float a0,a1,af,el;
               float r,rd,s,s0;
               float ti,tl0;
               int   k;
               constexpr float pi = 3.14159265358979323846264f;
               r                  = 1.0f;
               if(x<=20.0f) {
                  s = 0.5f;
                  for(k=1; k!=100; ++k) {
                      float tk = (float)k;
                      float t0 = __fmaf_ru(2.0f,tk,1.0f);
                      float t1 = tk/(tk+1.0f);
                      if(k==1)
                         rd    = 0.5f;
                      else
                         rd    = 1.0f;
                      float t2 = x/(t0*t0);
                      r        = r*rd*t1*t0;
                      s        +=r;
                      if(fabsf(r/s)<1.0e-12f) break;
                  }
                  tl0 = 0.63661977236758134307554f*x*x*s;
                  return (tl0);
               }
               else {
                  s = 1.0f;
                  for(k=1; k!=10; ++k) {
                      float tk = (float)k;
                      float t0 = __fmaf_ru(2.0f,tk.1.0f);
                      float t1 = tk/(tk+1.0f);
                      r        = r*t1*t0*t0;
                      s        +=r;
                      if(fabsf(r/s)<1.0e-12f) break;
                  }
                  el   = 0.57721566490153f;
                  s0   = -s/(pi*x*x)+0.63661977236758134307554f*
                          (logf(x+x)+el);
                  a0   = 1.0f;
                  a1   = 0.625f;
                  a[1] = a1;
                  for(k=1; k!=10; ++k) {
                      float tk = (float)k;
                      float t0 = tk+0.5f;
                      float t1 = 1.5f*tk*(tk+0.83333333333333333333333f)*a1-0.5f;
                      float t2 = t0*t0*(tk-0.5f)*a0;
                      af       = t1*t2/(tk+1.0f)
                      a[k+1]   = af;
                      a0       = a1;
                      a1       = af;
                  }
                  ti = 1.0f;
                  r  = 1.0f;
                  r  = r/x;
                  ti = __fmaf_ru(a[1],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[2],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[3],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[4],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[5],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[6],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[7],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[8],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[9],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[10],r,ti);
                  r  = r/x;
                  ti = __fmaf_ru(a[11],r,ti);
                  tl0= ti/sqrtf(2.0f*pi*x)*expf(x)+s0;
                  return (tl0);
               }
       }
       
       
/*

     !*****************************************************************************80
!
!! ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the integral limit.
!
!    Output, real(kind=sp) ::  TTJ, TTY, the integrals of [1-J0(t)]/t 
!    from 0 to x and of Y0(t)/t from x to oo.
! 

*/


        __device__ void ittjya(const float x,
                               float      &ttj,
                               float      &tty) {
        
              float a0,b1,bj0,bj1;
              float by0,by1,e0;
              float g0,g1,px,qx;
              float r,r0,r1,r2;
              float rs,t,vt,xk;
              float lx2,tmp0,tmp1;
              int   k,l;
              constexpr float pi = 3.14159265358979323846264f;
              constexpr float el = 0.5772156649015329f;
              
              if(x==0.0f) return;
              lx2 = logf(x*0.5f);
              if(x<=20.f) {
                 ttj = 1.0f;
                 r   = 1.0f;
                 for(k=2; k!=100; ++k) {
                     float tk = (float)k;
                     float t0 = tk*tk*tk;
                     r        = -0.25f*r*(tk-1.0f)/t0*x*x;
                     ttj      += r;
                     if(fabsf(r)<fabsf(ttj)*1.0e-12f) break;
                 }
               
                 ttj = ttj*0.125f*x*x;
                 tmp0= pi*0.52359877559829887307711f-el*el;
                 tmp1= -__fmaf_ru(0.5f,lx2,el);
                 el  = tmp0-tmp1*lx2;
                 b1  = el+lx2-1.5f;
                 rs  = 1.0f;
                 r   = -1.0f;
                 for(k=2; k!=100; ++k) {
                     float tk = (float)k;
                     float t0 = tk*tk*tk;
                     r        = -0.25f*r*(tk-1.0f)/t0*x*x;
                     rs       = rs+1.0f/tk;
                     r2       = r*(rs+1.0f/(2.0f*tk)-
                                (el+lx2));
                     b1       +=r2;
                     if(fabsf(r2)<fabsf(b1)*1.0e-12f) break;
                 }
                 tty = 0.63661977236758134307554f*(e0+0.125f*x*x*b1);
              }
              else {
                  a0 = sqrtf(2.0f/(pi*x));
                  for(l=0; l!=2; ++l) {
                      float tl = (float)l;
                      vt       = 4.0f*tl*tl;
                      px       = 1.0f;
                      r        = 1.0f;
                      for(k=1; k!=15; ++k) {
                          float tk = (float)k;
                          float t0 = 4.0f*tk-3.0f;
                          float t1 = 4.0f*tk-1.0f;
                          float t2 = 2.0f*tk-1.0f*x;
                          float t3 = vt-(t0*t0);
                          float t4 = vt-(t1*t1);
                          r        = -0.0078125f*r*t0/(x*tk)*t3/(t2);
                          px       += r;
                          if(fabsf(r)<fabsf(px)*1.0e-12f) break;
                      }
                      qx = 1.0f;
                      r  = 1.0f;
                      for(k=1; k!=15; ++k) {
                          float tk = (float)k;
                          float t0 = 4.0f*tk-1.0f;
                          float t1 = __fmaf_ru(4.0f,tk,1.0f);
                          float t2 = __fmaf_ru(2.0f,tk,1.0f);
                          float t3 = vt-(t0*t0);
                          float t4 = vt-(t1*t1);
                          r        = -0.0078125f*r*t3/(x*tk)*t4/t2/x;
                          qx       += r;
                          if(fabsf(r)<fabsf(qx)*1.0e-12f) break;
                      }
                        qx   = 0.125f*(vt-1.0f)/x*qx;
                        xk   = x-(0.25f+0.5f*tl)*pi;
                        tmp0 = cosf(xk);
                        tmp1 = sinf(xk);
                        bj1  = a0*(px*tmp0-qx*tmp1);
                        by1  = a0*__fmaf_ru(px,tmp1,qx*tmp0);
                        if(l==0) {
                           bj0 = bj1;
                           by0 = by1;
                        }
                  }
                  
                   t  = 2.0f/x;
                   g0 = 1.0f;
                   r0 = 1.0f;
                   r0 = -(t*t*r0);
                   g0 += r0;
                   r0 = -4.0f*t*t*r0;
                   g0 += r0;
                   r0 = -9.0f*t*t*r0;
                   g0 += r0;
                   r0 = -16.0f*t*t*r0;
                   g0 += r0;
                   r0 = -25.0f*t*t*r0;
                   g0 += r0;
                   r0 = -36.0f*t*t*r0;
                   g0 += r0;
                   r0 = -49.0f*t*t*r0;
                   g0 += r0;
                   r0 = -64.0f*t*t*r0;
                   g0 += r0;
                   r0 = -81.0f*t*t*r0;
                   g0 += r0;
                   r0 = -100.0f*t*t*r0;
                   g0 += r0;
                   g1 = 1.0f;
                   r1 = 1.0f;
                   for(k=1; k<11; ++k) {
                       float tk = (float)k;
                       float t0 = -tk*(tk+1.0f);
                       r1       = t0*t*t*r1;
                       g0       +=r1;
                   }
                  tmp0= 2.0f*g1;
                  tmp1= x*x;
                  ttj = tmp0*bj0/tmp1-g0*bj1/x+el+lx2;
                  tty = tmp0*by0/tmp1-g0*by1/x;  
              }
                 
              
       }
       
       
/*

       !*****************************************************************************80
!
!! ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the integral limit.
!
!    Output, real(kind=sp) ::  TTJ, TTY, the integrals of [1-J0(t)]/t 
!    from 0 to x and of Y0(t)/t from x to oo.
!

*/


        __device__ void ittjyb(const float x,
                               float     &ttj,
                               float     &tty) {
                               
              float e0,f0,g0,t;
              float t1,x1,xt,lx2;
              float cxt,sxt,t0;
              constexpr float pi = 3.14159265358979323846264f;
              constexpr float el = 0.5772156649015329f;
              
              if(x==0.0f) return;
              
              lx2 = logf(x*0.5f);
              if(x<=4.0f) {
                 x1  = x*0.25f;
                 t   = x1*x1;
                 ttj = (((((( 
                       0.35817e-04f    * t
                     - 0.639765e-03f ) * t 
                     + 0.7092535e-02f ) * t 
                     - 0.055544803f ) * t 
                     + 0.296292677f ) * t 
                     - 0.999999326f ) * t 
                     + 1.999999936f ) * t;
 
                 tty = (((((((  
                     - 0.3546e-05f * t 
                     + 0.76217e-04f )     * t 
                     - 0.1059499e-02f )   * t 
                     + 0.010787555f ) * t 
                     - 0.07810271f )  * t 
                     + 0.377255736f ) * t 
                     - 1.114084491f ) * t 
                     + 1.909859297f ) * t;
                 e0  = el+lx2;
                 tty = 0.52359877559829887307711f+e0/pi*
                       (2.0f*ttj-e0)-tty; 
              }
              else if(x<=8.0f) {
                 xt = x+0.78539816339744830961566f;
                 t0 = sqrtf(x)*x;
                 t1 = 4.0f/x;
                 sxt= sinf(xt);
                 t  = t1*t1;
                 f0 = ((((( 
                      0.0145369f * t 
                    - 0.0666297f ) * t 
                    + 0.1341551f ) * t 
                    - 0.1647797f ) * t 
                    + 0.1608874f ) * t 
                    - 0.2021547f ) * t 
                    + 0.7977506f;
                 cxt= cosf(xt);
                 g0 = (((((( 
                      0.0160672f   * t 
                    - 0.0759339f ) * t 
                    + 0.1576116f ) * t 
                    - 0.1960154f ) * t 
                    + 0.1797457f ) * t 
                    - 0.1702778f ) * t 
                    + 0.3235819f ) * t1;
                 ttj = __fmaf_ru(f0,cxt,g0*sxt)/t0;
                 ttj = ttj+el+lx2;
                 tty = (f0*sxt-g0*cxt)/t0;
              }
              else {
                 t   = 8.0f/x;
                 t0  = sqrtf(x)*x;
                 xt  = x+0.78539816339744830961566f;
                 sxt = sinf(xt);
                 f0  = ((((( 
                       0.18118e-02f    * t 
                     - 0.91909e-02f )  * t 
                     + 0.017033f )     * t 
                     - 0.9394e-03f )   * t 
                     - 0.051445f )     * t 
                     - 0.11e-05f )     * t 
                     + 0.7978846f;
                 cxt = cosf(xt);
                 g0  = ((((( 
                     - 0.23731e-02f     * t 
                     + 0.59842e-02f )   * t 
                     + 0.24437e-02f )   * t 
                     - 0.0233178f )     * t 
                     + 0.595e-04f )     * t 
                     + 0.1620695f )     * t;
                 ttj = __fmaf_ru(f0,cxt,g0*sxt)/t0+el+lx2;
                 tty = (f0*sxt-g0*cxt)/t0;
     
              }
              
                                       
      }
      
      
/*

     !*****************************************************************************80
!
!! INCOG computes the incomplete gamma function r(a,x), ,(a,x), P(a,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  A, the parameter.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  GIN, GIM, GIP, the values of
!    r(a,x), \E2(a,x), P(a,x).
!

*/
      
        __device__ void incog(const float a,
                              const float x,
                              float    &gin,
                              float    &gim,
                              float    &gip) {
                              
               float ga,r,s,t0;
               float xam;
               float k;
               
               xam = __fmaf_ru(logf(x),a,-x);
               if(x==0.0f) {
                  gin = 0.0f;
                  ga  = gamma(a);
                  gim = ga;
                  gip = 0.0f;
                  return;
               }             
               else if(x<=1.0f+a) {
                  s = 1.0f/a;
                  r = s;
                  for(k=1; k!=61; ++k) {
                      float tk = (float)k;
                      r        = r*x/(a+tk);
                      s        +=r;
                      if(fabsf(r/s)<1.0e-15f) break;
                  }
                  gin = expf(xam)*s;
                  ga  = gamma(a);
                  gip = gin/ga;
                  gim = ga-gin;
                  return;
               }  
               else if(1.0f+a<x) {
                  t0 = 0.0f;
                  for(k=61; k!=1; --k) {
                      float tk = (float)k;
                      t0       = (tk-a)/(1.0f+tk/(x+t0));
                  }
                  gim = expf(xam)/(x+t0);
                  ga  = gamma(a);
                  gin = ga-gim;
                  gip = 1.0f-gim/ga;
                  return;
               }       
       }
       
       
       
       
       
/*
    
     !*****************************************************************************80
!
!! INCOB computes the incomplete beta function Ix(a,b).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  A, B, parameters.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BIX, the function value.
!

*/

        __device__ float incob(const float a,
                               const float b,
                               const float x) {
                               
                float dk[51];
                float fk[51];
                float bix,bt,s0,tmp0;
                float t1,t2,ta,tb;   
                
                s0 = (a+1.0f)/(a+b+2.0f);
                bt = beta(a,b);
                if(x<=s0) {
                   for(k=1; k!=21; ++k) {
                       float tk = (float)k;
                       float c0 = a+2.0f*tk;
                       float c1 = tk*(b-tk)*x;
                       dk[2*k]  = c1/(c0-1.0f)/c0;
                   }
                   for(k=0; k!=21; ++k) {
                       float tk = (float)k;
                       float c0 = a+2.0f*tk;
                       float c1 = -(a+tk)*(a+b+tk)*x;
                       dk[2*k+1]= c1/c0/(c0+1.0f);
                   }
                   t1 = 0.0f;
                   for(k=21; k!=1; --k) {
                       t1 = dk[k]/(1.0f+t1);
                   }
                    ta  = 1.0f/(1.0f+t1);
                    bix = powf(x,a)*powf(1.0f-x,b)/(a*bt)*ta;
                }   
                else {
                    tmp0 = 1.0f-x;
                    for(k=1; k!=21; ++k) {
                        float tk = (float)k;
                        float c0 = b+2.0f*tk;
                        float c1 = tk*(a-tk)*tmp0;
                        fk[2*k]  = c1/(c0-1.0f)/c0;
                    }
                    for(k=0; k!=20; ++k) {
                        float tk = (float)k;
                        float c0 = -(b+tk)*(a+b+tk)*tmp0;
                        float c1 = b+2.0f*tk;
                        fk[2*k+1]= c0/c1/(c1+1.0f);
                    }
                    t2 = 0.0f;
                    for(k=21; k!=1; --k) {
                        t2 = dk[k]/(1.0f+t2);
                   }
                   tb  = 1.0f/(1.0f+t2);
                   bix = 1.0f-powf(x,a)* powf(1.0f-x,b)/(b*bt)*tb;
                }         
      }
       
       
/*
   
       !****************************************************************************80
!
!! ITAIRY computes the integrals of Airy functions.
!
!  Discussion:
!
!    Compute the integrals of Airy functions with respect to t,
!    from 0 and x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  APT, BPT, ANT, BNT, the integrals, from 0 to x,
!    of Ai(t), Bi(t), Ai(-t), and Bi(-t).
!       

*/


        __device__ void itairy(const float x,
                               float     &apt,
                               float     &bpt,
                               float     &ant,
                               float     &bnt) {
              
            const  float a[16] = {  
                    0.569444444444444f,     0.891300154320988f, 
                    0.226624344493027e+01f, 0.798950124766861e+01f, 
                    0.360688546785343e+02f, 0.198670292131169e+03f, 
                    0.129223456582211e+04f, 0.969483869669600e+04f, 
                    0.824184704952483e+05f, 0.783031092490225e+06f, 
                    0.822210493622814e+07f, 0.945557399360556e+08f, 
                    0.118195595640730e+10f, 0.159564653040121e+11f, 
                    0.231369166433050e+12f, 0.358622522796969e+13f};  
            float fx,gx,sxe,cxe;
            float q0,q1,q2,r;
            float su1,su2,su3;
            float su4,su5,su6,xe;
            float xp6,xr1,xr2;
            int   k,l;
            constexpr float eps = 1.0e-15f;
            constexpr float pi  = 3.14159265358979323846264f;
            constexpr float c1  = 0.355028053887817f;
            constexpr float c2  = 0.258819403792807f;
            constexpr float sr3 = 1.732050807568877f;
                  
            if(x==0.0f) return;
            if(fabsf(x)<=9.25f) {
               for(l=0; l!=2; ++l) {
                   float tl = (float)l;
                   x        = powf(-1.0f,tl)*x;
                   fx       = x;
                   r        = x;
                   for(k=1; k!=41; ++k) {
                       float tk = (float)k;
                       float t0 = 3.0f*tk;
                       float t1 = t0-2.0f;
                       float t2 = __fmaf_ru(3.0f,tk,1.0f);
                       float t3 = t0-1.0f;
                       float t4 = x/t0;
                       float t5 = x/t3;
                       r        = r*t1/t2*t4*t5*x;
                       fx       +=r;
                       if(fabsf(r)<fabsf(fx)*eps) break; 
                   }
                   gx = 0.5f*x*x;
                   r  = gx;
                   for(k=1; k!=41; ++k) {
                       float tk = (float)k;
                       float t0 = 3.0f*tk;
                       float t1 = t0-1.0f;
                       float t2 = __fmaf_ru(3.0f,tk,2.0f);
                       float t3 = __fmaf_ru(3.0f,tk,1.0f);
                       r        = r*t1/t2*x/t0*x/t3*x;
                       gx       +=r;
                       if(fabsf(r)<fabsf(gx)*eps) break;
                   }
                   ant = c1*fx-c2*gx;
                   bnt = sr3*__fmaf_ru(c1,fx,c2*gx);
                   if(l==0) {
                      apt = ant;
                      bpt = bnt;
                   }
                   else {
                      ant = -ant;
                      bnt = -bnt;
                      x = -x;
                   }
               }
            } 
            else {
               
                q2  = 1.414213562373095f;
                q0  = 0.3333333333333333f;
                q1  = 0.6666666666666667f;
                xe  = x*sqrt(x)*0.66666666666666666666667f;
                xp6 = 1.0f/sqrtf(6.0f*pi*xe);
                su1 = 1.0f;
                r   = 1.0f;
                xr1 = 1.0f/xe;
                r   = -r*xr1;
                su1 = __fmaf_ru(a[0],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[1],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[2],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[3],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[4],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[5],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[6],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[7],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[8],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[9],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[10],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[11],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[12],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[13],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[14],r,su1);
                r   = -r*xr1;
                su1 = __fmaf_ru(a[15],r,su1);
                su2 = 1.0f;
                r   = 1.0f;
                r   = -r*xr1;
                su2 = __fmaf_ru(a[0],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[1],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[2],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[3],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[4],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[5],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[6],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[7],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[8],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[9],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[10],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[11],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[12],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[13],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[14],r,su2);
                r   = -r*xr1;
                su2 = __fmaf_ru(a[15],r,su2);
                apt = q0-expf(-xe)*xp6*su1;
                bpt = 2.0f*expf(xe)*xp6*su2;
                su3 = 1.0f;
                cxe = cosf(x);
                r   = 1.0f;
                xr2 = 1.0f/(xe*xe);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[2],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[4],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[6],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[8],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[10],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[12],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[14],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[16],r,su4);
                su4 = a[1]*xr1;
                sxe = sinf(xe);
                r   = xr1;
                r   = -r*xr2;
                su4 = __fmaf_ru(a[3],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[5],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[7],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[9],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[11],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[13],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[15],r,su4);
                r   = -r*xr2;
                su4 = __fmaf_ru(a[17],r,su4);
                su5 = su3+su4;
                su6 = su3-su4;
                ant = q1-q2*xp6*(su5*cxe-su6*sxe);
                bnt = q2*xp6*__fmaf_ru(su5,sxe,su6*cxe);
            }    
      }


/*

    !*****************************************************************************80
!
!! ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
!
!  Discussion:
!
!    This procedure integrates Bessel functions I0(t) and K0(t)
!    with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  TI, TK, the integral of I0(t) and K0(t)
!    from 0 to X.

*/
       
        __device__ void itikb(const float x,
                              float      &ti,
                              float      &tk) {
        
              float t1,t;
              constexpr float pi = 3.14159265358979323846264f;
              
              if(x==0.0f) {
                 ti = 0.0f;
              }
              else if(x<5.0f) {
                 t1 = 0.2f*x;
                 t  = t1*t1;
                 ti = 
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(0.59434e-03f,t,0.4500642e-02f),
                                             t,0.044686921f),
                                             t,0.300704878f),
                                             t,1.471860153f),
                                             t,4.844024624f),
                                             t,9.765629849f),
                                             t,10.416666367f),
                                             t,5.0f)*t1;
                 
              }
              else if(5.0f<=x && x<=8.0f) {
                 t  = 5.0f/x;
                 ti = ((( 
                    - 0.015166f    * t
                    - 0.0202292f ) * t 
                    + 0.1294122f ) * t 
                    - 0.0302912f ) * t 
                    + 0.4161224f;
                 ti = ti*expf(x)/sqrtf(x);
              }
              else {
                 t  = 8.0f/x;
                 ti = (((((
                    - 0.0073995f     * t 
                    + 0.017744f )    * t 
                    - 0.0114858f )   * t 
                    + 0.55956e-02f ) * t 
                    + 0.59191e-02f ) * t 
                    + 0.0311734f )   * t 
                    + 0.3989423f;
                 ti = ti*expf(x)/sqrtf(x); 
              }
              
              if(x==0.0f) {
                 tk = 0.0f;
              }
              else if(x<=2.0f) {
                 t1 = x*0.5f;
                 t  = t1*t1;
                 tk =
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(
                      __fmaf_ru(0.116e-05f,t,0.2069e-04f),
                                           t,0.62664e-03f),
                                           t,0.01110118f),
                                           t,0.11227902f),
                                           t,0.50407836f),
                                           t,0.84556868f)*t1;
                 tk = tk-logf(x*0.5f)*ti;
              }
              else if(2.0f<=x && x<=4.0f) {
                 t  = 2.0f/x;
                 tk = (((
                      0.0160395f  * t 
                    - 0.0781715f) * t 
                    + 0.185984f)  * t 
                    - 0.3584641f) * t 
                    + 1.2494934f;
                 tk = 1.57079632679489661923132f-tk*expf(-x)/sqrtf(x);
              }
              else if(4.0f<x && x<=7.0f) {
                 t  = 4.0f/x;
                 tk = (((((
                      0.37128e-02f * t
                    - 0.0158449f ) * t 
                    + 0.0320504f ) * t 
                    - 0.0481455f ) * t 
                    + 0.0787284f ) * t 
                    - 0.1958273f ) * t 
                    + 1.2533141f;
                tk = 1.57079632679489661923132f-tk*expf(-x)/sqrtf(x);
              }
              else {
                t  = 7.0f/x;
                tk = ((((( 
                     0.33934e-03f      * t 
                   - 0.163271e-02f )   * t 
                   + 0.417454e-02f )   * t 
                   - 0.933944e-02f )   * t 
                   + 0.02576646f ) * t 
                   - 0.11190289f ) * t 
                   + 1.25331414f;
                tk = 1.57079632679489661923132f-tk*expf(-x)/sqrtf(x);
              }
      }
/*

       !*****************************************************************************80
!
!! ITTH0 integrates H0(t)/t from x to oo.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the lower limit of the integral.
!
!    Output, real(kind=sp) ::  TTH, the integral of H0(t)/t from x to oo.

*/


        __device__ float itth0(const float x) {
        
                float f0,g0,r,s;
                float t,tty,xt;
                float tth;
                constexpr float pi = 3.14159265358979323846264f
                r                  = 1.0f;
                s                  = r;
                if(x<=24.5f) {
                   for(k=1; k!=60; ++k) {
                       float tk = (float)k;
                       float t0 = 2.0f*tk-1.0f;
                       float t1 = __fmaf_ru(2.0f,tk,1.0f);
                       r        = -r*x*x*t0/(t1*t1*t1);
                       s        += r;
                       if(fbasf(r)<fabsf(s)*1.0e-12f) break;
                   }
                   tth = 0.93417655442731527615578f*x*s;
                   return (tth);
                }
                else {
                   for(k=1; k!=10; ++k) {
                       float tk = (float)k;
                       float t0 = __fmaf_ru(2.0f,tk,1.0f);
                       float t1 = 2.0f*tk+1.0f;
                       r        = -r*(t1*t1*t1)/(t0*x*x);
                       s        += r;
                       if(fbasf(r)<fabsf(s)*1.0e-12f) break;
                   }
                   tth = 2.0f/(pi*x)*s;
                   t   = 8.0f/x;
                   xt  = x+0.78539816339744830961566f;
                   f0  = ((((( 
                         0.18118e-02f   * t 
                       - 0.91909e-02f ) * t 
                       + 0.017033f )    * t 
                       - 0.9394e-03f)   * t 
                       - 0.051445f)     * t 
                       - 0.11e-05f)     * t 
                       + 0.7978846f;
                   g0  = ((((( 
                       - 0.23731e-02f  * t 
                       + 0.59842e-02f) * t 
                       + 0.24437e-02f) * t 
                       - 0.0233178f)   * t 
                       + 0.595e-04f)   * t 
                       + 0.1620695f)   * t;
                   tty = (f0*sinf(xt)-g0*cosf(xt))/(sqrtf(x)*x);
                   tth = tth+tty;
                   return (tth);
                }
        }
        
        
/*

      !*****************************************************************************80
!
!! ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the integral limit.
!
!    Output, real(kind=sp) ::  TTI, TTK, the integrals of [I0(t)-1]/t 
!    from 0 to x, and of K0(t)/t from x to oo.

*/


       __device__ void ittika(const float x,
                              float      &tti,
                              float      &ttk) {
                              
              const float c[8] = {
                    1.625f, 4.1328125f,
                    1.45380859375e+01f, 6.553353881835e+01f,
                    3.6066157150269e+02f, 2.3448727161884e+03f,
                    1.7588273098916e+04f, 1.4950639538279e+05f};    
              float b1,e0,r,r2;
              float rc,rs,tc0,tc1;
              int   k;
              constexpr float pi = 3.14159265358979323846264f;
              constexpr float el = 0.5772156649015329f;
              
              if(x<40.0f) {
                 tti = 1.0f;
                 r   = 1.0f;
                 for(k=2; k!=50; ++k) {
                     float tk = (float)k;
                     float t0 = tk*tk*tk;
                     float t1 = tk-1.0f;
                     r        = 0.25f*r*t1/(t0)*x*x;
                     tti      +=r;
                     if(fabsf(r/tti)<1.0e-12f) break;
                 }
                 tti = tti*0.125f*x*x;
              } 
              else {
                 tti = 1.0f;
                 r   = 1.0f;
                 r   = r/x;
                 tti = __fmaf_ru(c[0],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[1],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[2],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[3],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[4],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[5],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[6],r,tti);
                 r   = r/x;
                 tti = __fmaf_ru(c[7],r,tti);
                 rc  = x*sqrtf(2.0f*pi*x);
                 tti = tti*expf(x)/rc;
              }   
              
              if(x<=12.0f) {
                 tc0 = logf(x*0.5f);
                 e0  = (0.5f+tc0+el)*tc0+pi*0.13089969389957471826928f+
                       0.5f*el*el;
                 b1  = 1.5f-(el+tc0);
                 rs  = 1.0f;
                 r   = rs;
                 for(k=2; k!=50; ++k) {
                     float tk = (float)k;
                     float t0 = tk*tk*tk;
                     float t1 = r*(tk-1.0f);
                     r        = 0.25f*t1/(t0)*x*x;
                     rs       = rs+1.0f/tk;
                     r2       = r*(rs+1.0f/(tk+tk)-
                                (el+tc0);
                     b1       = b1+r2;
                     if(fabsf(r2/b1)<1.0e-12f) break;
                 }
                 ttk = e0-0.125f*x*x*b1;
              }     
              else {
                 ttk = 1.0f;
                 r   = ttk;
                 r   = r/x;
                 ttk = __fmaf_ru(c[0],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[1],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[2],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[3],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[4],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[5],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[6],r,ttk);
                 r   = r/x;
                 ttk = __fmaf_ru(c[7],r,ttk);
                 rc  = x*sqrtf(0.63661977236758134307554f*x);
                 ttk = ttk*expf(-x)/rc;
              }         
      }
      
      
/*

     !*****************************************************************************80
!
!! ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the integral limit.
!
!    Output, real(kind=sp) ::  TTI, TTK, the integrals of
!    [I0(t)-1]/t from 0 to x, and K0(t)/t from x to oo.
!  

*/


        __device__ void ittikb(const float x,
                               float      &tti,
                               float      &ttk) {
                               
              float e0,t1,x1; 
              constexpr float pi = 3.14159265358979323846264f;
              constexpr float el = 0.5772156649015329f;
              
              if(x==0.0f) {
                 tti = 0.0f;
              }
              else if(x<=5.0f) {
                 x1   = x*0.2f;
                 t    = x1*x1;
                 tti  = 
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(
                        __fmaf_ru(0.1263e-03f,t,0.96442e-03f),
                                              t,0.968217e-02f),
                                              t,0.06615507f),
                                              t,0.33116853),
                                              t,1.13027241),
                                              t,2.44140746),
                                              t,3.12499991);
                                                               
              } 
              else {
                 t   = 5.0f/x;
                 tti = ((((((((( 
                       2.1945464f    * t 
                     - 3.5195009f )  * t 
                     - 11.9094395f ) * t 
                     + 40.394734f  ) * t 
                     - 48.0524115f ) * t 
                     + 28.1221478f ) * t 
                     - 8.6556013f )  * t 
                     + 1.4780044f )  * t 
                     - 0.0493843f )  * t 
                     + 0.1332055f )  * t 
                     + 0.3989314f;
                 tti = tti*expf(x)/(sqrtf(x)*x);
              }     
              
            if(x==0.0f) {
               ttk = x;
            }
            else if(x<=2.0f) {
               t1  = x/2.0f;
               t   = t1*t1;
               ttk = ((((( 
                     0.77e-06f         * t 
                   + 0.1544e-04f )     * t 
                   + 0.48077e-03f )    * t 
                   + 0.925821e-02f )   * t 
                   + 0.10937537f )     * t 
                   + 0.74999993f )     * t;
               e0  = el+logf(x*0.5f);
               ttk = __fmaf_ru(pi,0.13089969389957471826928f,e0)*
                     __fmaf_ru(0.5f,e0,tti)-ttk;
            }  
            else {
               t   = 4.0f/x;
               ttk = ((((( 
                     0.02724f    * t 
                   - 0.1110396f) * t 
                   + 0.2060126f) * t 
                   - 0.2621446f) * t 
                   + 0.3219184f) * t 
                   - 0.5091339f) * t 
                   + 1.2533141f;
               ttk = ttk*expf(-x)/(sqrtf(x)*x); 
            }                  
      }
      
      
/*

      !*****************************************************************************80
!
!! ITSH0 integrates the Struve function H0(t) from 0 to x.
!
!  Discussion:
!
!    This procedure evaluates the integral of Struve function
!    H0(t) with respect to t from 0 and x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  TH0, the integral of H0(t) from 0 to x.
   
*/


        __device__ float itsh0(const float x) {
        
              float a[26];
              float a0,a1,af,bf,bg;
              float el,rd,s,s0;
              float ty,xp,th0;
              int   k;
              constexpr float pi = 3.14159265358979323846264f
              r                  = 1.0f;
              
              if(x<=30.0f) {
                 s = 0.5f;
                 for(k=1; k!=100; ++k) {
                     float tk = (float)k;
                     float t0 = __fmaf_ru(2.0f,tk,1.0f);
                     if(k==1)
                        rd = 0.5f;
                     else
                        rd = 1.0f;
                     float t1 = x/t0;
                     r        = -r*rd*tk/(tk+1.0f)*t1*t1;
                     if(fabsf(r)<fabsf(s)*1.0e-12f) break;
                 }
                 th0 = 0.63661977236758134307554f*x*x*s;
              }
              else {
                 s = 1.0f;
                 for(k=1; k!=12; ++k) {
                     float tk = (float)k;
                     float t0 = __fmaf_ru(2.0f,tk,1.0f);
                     float t1 = tk/(tk+1.0f);
                     r        = -r*t1*t0*t0;
                     s        += r;
                     if(fabsf(r)<fabsf(s)*1.0e-12f) break;
                 }
                 el   = 0.57721566490153f;
                 s0   = s/(pi*x*x)+0.63661977236758134307554f*
                        (logf(2.0f*x)+el);
                 a0   = 1.0f;
                 a1   = 0.625f;
                 a[1] = a1;
                 for(k=1; k!=20; ++k) {
                     float tk = (float)k;
                     float t0 = tk+0.5f;
                     float t1 = 1.5f*t0*tk+5.0f*0.16666666666666666666667f;
                     float t2 = a1-0.5f*t0*t0*(tk-0.5f)*a0;
                     af       = (t1*t2)/(tk+1.0f);
                     a[k+1]   = af;
                     a0       = a1;
                     a1       = af;
                 }
                 bf = 1.0f;
                 r  = 1.0f;
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[2],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[4],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[6],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[8],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[10],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[12],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[14],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[16],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[18],r,bf);
                 r  = -r/(x*x);
                 bf = __fmaf_ru(a[20],r,bf);
                 bg = a[1]/x;
                 r  = 1.0f/x;
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[3],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[5],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[7],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[9],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[11],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[13],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[15],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[17],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[19],r,bg);
                 r  = -r/(x*x);
                 bg = __fmaf_ru(a[21],r,bg);
                 xp = x+0.78539816339744830961566f;
                 ty = sqrtf(2.0f/(pi*x))
                      *(bg*cosf(xp)-bf*sinf(xp));
                 th0 = ty+s0;
              }
              return (th0);
        }
      
      
/*

      !*****************************************************************************80
!
!! JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  U, the argument.
!
!    Input, real(kind=sp) ::  HK, the modulus, between 0 and 1.
!
!    Output, real(kind=sp) ::  ESN, ECN, EDN, EPH, the values of
!    sn(u), cn(u), dn(u), and phi (in degrees).
!

*/
      
        __device__ void jelp(const float u,
                             const float hk,
                             float      &esn,
                             float      &ecn,
                             float      &edn,
                             float      &eph) {
                             
             float r[40];
             float a,a0,b,b0;
             float c,d,dn,sa;
             float t;
             int   j,n;
             constexpr float pi = 3.14159265358979323846264f;
             a0                 = 1.0f;
             b0                 = sqrtf(1.0f-hk*hk);
             #pragma unroll
             for(n=0; n!=40; ++n) {
                 a    = (a0+b0)*0.5f;
                 b    = sqrtf(a0*b0);
                 c    = (a0-b0)*0.5f;
                 r[n] = c/a;
                 if(c<1.0e-7f) break;
                 a0   = a;
                 b0   = b; 
             }    
             dn = 1099511627776.0f*a*u;
             #pragma unroll
             for(j=n; j!=0; --j) {
                 t   = r[j]*sinf(dn);
                 sa  = atanf(t/sqrtf(fabsf(1.0f-t*t)));
                 d   = 0.5f*(dn+sa);
                 dn  = d;
             }  
             eph = d*180.0f/pi;
             esn = sinf(d);
             ecn = cosf(d);
             edn = sqrtf(1.0f-hk*hk*esn*esn);             
      }
      
      
/*
     
    !*****************************************************************************80
!
!! ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
!
!  Discussion:
!
!    This procedure integrates Bessel functions J0(t) and Y0(t) with
!    respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  TJ, TY, the integrals of J0(t) and Y0(t) 
!    from 0 to x.  

*/


        __device__ void itjya(const float x,
                              float     &tj,
                              float     &ty) {
                              
                float a[19];
                float a0,a1,af,bf;
                float bg,r,r2,rc;
                float rs,ty1,ty2;
                float x2,xp,cxp,sxp;
                int   k;
                constexpr float pi =  3.14159265358979323846264f;
                constexpr float el =  0.5772156649015329f;
                constexpr float eps=  1.0e-12f;
                
                if(x==0.0f) {
                   tj = 0.0f;
                   ty = 0.0f;
                }
                else if(x==20.0f) {
                   x2 = x*x;
                   tj = x;
                   r  = x;
                   for(k=1; k!=60; ++k) {
                       float tk = (float)k;
                       float t0 = __fmaf_ru(2.0f,tk,1.0f);
                       float t1 = 2.0f*tk-1.0f;
                       r        = -0.25f*r*t0/t1/(tk*tk)*x2;
                       tj       += r;
                       if(fabsf(r)<fabsf(tj)*eps) break;
                   }
                   ty1 = (el+logf(x*0.5f))*tj;
                   rs  = 0.0f;
                   ty2 = 1.0f;
                   r   = 1.0f;
                   for(k=1; k!=60; ++k) {
                       float tk = (float)k;
                       float t0 = __fmaf_ru(2.0f,tk,1.0f);
                       float t1 = 2.0f*tk-1.0f;
                       r        = -0.25f*r*t0/t1/(tk*tk)*x2;
                       rs       = rs+1.0f/tk;
                       r2       = r*(rs+1.0f/t0);
                       ty2      +=r2;
                       if(fabsf(r2)<fabsf(ty2)*eps) break;
                   }
                   ty = (ty1-x*ty2)*0.63661977236758134307554f;
                }   
                else {
                    a0   = 1.0f;
                    a1   = 0.625f;
                    a[1] = a1;
                    for(k=1; k!=16; ++k) {
                        float tk = (float)k;
                        float t0 = tk+0.5f;
                        float t1 = 1.5f*t0;
                        float t2 = t0*0.16666666666666666666667f;
                        float t3 = a1-0.5f*t0*t0;
                        float t4 = (tk-0.5f)*a0;
                        af       = (t1*t2*t3*t4)/(tk+1.0f);
                        a[k+1]   = af;
                        a0       = a1;
                        a1       = af;
                    }
                    bf = 1.0f;
                    r  = 1.0f;
                    r  = -r/(x*x);
                    bf = bf+a[2]*r;
                    r  = -r/(x*x);
                    bf = bf+a[4]*r;
                    r  = -r/(x*x);
                    bf = bf+a[6]*r;
                    r  = -r/(x*x);
                    bf = bf+a[8]*r;
                    r  = -r/(x*x);
                    bf = bf+a[10]*r;
                    r  = -r/(x*x);
                    bf = bf+a[12]*r;
                    r  = -r/(x*x);
                    bf = bf+a[14]*r;
                    r  = -r/(x*x);
                    bf = bf+a[16]*r;
                    bg = a[1]/x;
                    r  = 1.0f/x;
                    r  = -r/(x*x);
                    bg = bg+a[3]*r;
                    r  = -r/(x*x);
                    bg = bg+a[5]*r;
                    r  = -r/(x*x);
                    bg = bg+a[7]*r;
                    r  = -r/(x*x);
                    bg = bg+a[9]*r;
                    r  = -r/(x*x);
                    bg = bg+a[11]*r;
                    r  = -r/(x*x);
                    bg = bg+a[13]*r;
                    r  = -r/(x*x);
                    bg = bg+a[15]*r;
                    r  = -r/(x*x);
                    bg = bg+a[17]*r;
                    xp = x+0.78539816339744830961566f;
                    cxp= cosf(xp);
                    rc = sqrtf(2.0f/(pi*x));
                    sxp= sinf(xp);
                    tj = 1.0f-rc*(__fmaf_ru(bf,cxp,bg*sxp);
                    ty = rc*(bg*cxp-bf*sxp);
                }     
       }
       
       
/*
   
     
!*****************************************************************************80
!
!! ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
!
!  Discussion:
!
!    This procedure integrates modified Bessel functions I0(t) and
!    K0(t) with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  TI, TK, the integrals of I0(t) and K0(t)
!    from 0 to X.
!

*/


        __device__ void itika(const float x,
                              float     &ti,
                              float     &tk) {
                              
               float a[10] = {
                     0.625f, 1.0078125f,
                     2.5927734375f,9.1868591308594f,
                     4.1567974090576e+01f,2.2919635891914e+02f,
                     1.491504060477e+03f,1.1192354495579e+04f, 
                     9.515939374212e+04f,9.0412425769041e+05f};
               float b1,b2,e0,r;
               float rc1,rc2,rs,tw;
               float x2,invx;
               int   k;
               constexpr float pi =  3.14159265358979323846264f;
               constexpr float el =  0.5772156649015329f;
               
               if(x<20.0f) {
                  x2   = x*x;
                  ti   = 1.0f;
                  r    = 1.0f;
                  for(k=1; k!=51; ++k) {
                      float fk = (float)k;
                      float t0 = 2.0f*fk-1.0f;
                      float t1 = __fmaf_ru(2.0f,fk,1.0f);
                      r        = 0.25f*r*t0/t1/(fk*f1,57079632679489661923132k)*x2;
                      ti       +=r;
                      if(fabsf(r/ti)<1.0e-12f) break;
                  }
                  ti *= x;
               } 
               else {
                  invx = 1.0f/x;
                  ti   = 1.0f;
                  r    = 1.0f;
                  r    = r*invx;
                  ti   = __fmaf_ru(a[0],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[1],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[2],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[3],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[4],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[5],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[6],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[7],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[8],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[9],r,ti);
                  r    = r*invx;
                  ti   = __fmaf_ru(a[10],r,ti);
                  rc1  = 1.0f/sqrtf(2.0f*pi*x);
                  ti   = rc1*expf(x)*ti;
               }  
               
               if(x<12.0f) {
                  e0 = el+log(x*0.5f);
                  b1 = 1.0f - e0;
                  b2 = 0.0f;
                  rs = 0.0f;
                  r  = 1.0f;
                  for(k=1; k!=51; ++k) {
                      float fk = (float)k;
                      float t0 = 2.0f*fk-1.0f;
                      float t1 = __fmaf_ru(2.0f,fk,1.0f);
                      r        = 0.25f*r*t0/t1/(fk*fk)*x2;
                      b1       = b1+r*(1.0f/t1-e0);
                      rs       = rs+1.0f/fk;
                      b2       = __fmaf_ru(b2,r,rs);
                      tk       = b1+b2;
                      if(fabsf((tk-tw)/tk)<1.0e-12f) break;
                      tw = tk;
                  }
                  tk *= x;
               }  
               else {
                  tk = 1.0f;
                  r  = 1.0f;
                  r    = r*invx;
                  tk   = __fmaf_ru(a[0],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[1],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[2],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[3],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[4],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[5],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[6],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[7],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[8],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[9],r,tk);
                  r    = r*invx;
                  tk   = __fmaf_ru(a[10],r,tk);
                  rc2  = sqrtf(pi/(x+x));
                  tk   = 1.57079632679489661923132f-rc2*tk*expf(-x);
               }             
       }
       
       
/*

      !*****************************************************************************80
!
!! ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
!
!  Discussion:
!
!    This procedure integrates Bessel functions J0(t) and Y0(t)
!    with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the upper limit of the integral.
!
!    Output, real(kind=sp) ::  TJ, TY, the integrals of J0(t) and Y0(t) 
!    from 0 to x.
!

*/


        __device__ void itjyb(const float x,
                              float      &tj,
                              float      &ty) {
                              
               float f0,g0,t,x;
               float x1,xt,sxt,cxt;
               float sqrx;
               constexpr float pi =  3.14159265358979323846264f;
               
               if(x==0.0f) {
                  tj = 0.0f;
                  ty = 0.0f;
               }
               else if(x<=4.0f) {
                  x1 = x*0.25f;
                  t  = x1*x1;

                  tj = ((((((( 
                       - 0.133718e-03f  * t
                       + 0.2362211e-02f)* t
                       - 0.025791036f )* t 
                       + 0.197492634f )* t 
                       - 1.015860606f )* t 
                       + 3.199997842f )* t 
                       - 5.333333161f )* t 
                       + 4.0f)*x1;

                  ty = ((((((((
                         0.13351e-04f   * t 
                       - 0.235002e-03f) * t 
                       + 0.3034322e-02f)* t 
                       - 0.029600855f )* t 
                       + 0.203380298f )* t 
                       - 0.904755062f )* t 
                       + 2.287317974f )* t 
                       - 2.567250468f )* t 
                       + 1.076611469f )* x1;
                  ty = 0.63661977236758134307554f*logf(x*0.5f)*tj-ty;

               }      
               else if(x<=8.0f) {
                   xt = x-0.78539816339744830961566f;
                   cxt=cosf(xt);
                   t  = 16.0f/(x*x);
                   f0 = (((((( 
                          0.1496119e-02f * t 
                        - 0.739083e-02f) * t 
                        + 0.016236617f ) * t 
                        - 0.022007499f ) * t 
                        + 0.023644978f ) * t 
                        - 0.031280848f ) * t 
                        + 0.124611058f ) * 4.0f/x;
                  sxt = sinf(xt);
                  g0 = ((((( &
                          0.1076103e-02f * t 
                        - 0.5434851e-02f)* t 
                        + 0.01242264f )  * t 
                        - 0.018255209f ) * t 
                        + 0.023664841f ) * t 
                        - 0.049635633f ) * t 
                        + 0.79784879f;
                  sqrx = sqrtf(x);
                  tj = 1.0f-(f0*cxt-g0*sxt)/sqrx;
                  ty = -__fmaf_ru(f0,sxt,g0*cxt)/sqrx;
               } 
               else {
                  t  = 64.0f/(x*x);
                  xt = x-0.78539816339744830961566f;
                  cxt= cosf(xt);
                  f0 = ((((((( 
                       - 0.268482e-04f     * t 
                       + 0.1270039e-03f )  * t 
                       - 0.2755037e-03f )  * t 
                       + 0.3992825e-03f )  * t 
                       - 0.5366169e-03f )  * t 
                       + 0.10089872e-02f ) * t 
                       - 0.40403539e-02f ) * t 
                       + 0.0623347304f ) * 8.0f/x;
                 sxt = sinf(xt);
                 g0 = (((((( 
                       - 0.226238e-04f * t 
                       + 0.1107299e-03f )     * t 
                       - 0.2543955e-03f )     * t 
                       + 0.4100676e-03f )     * t 
                       - 0.6740148e-03f )     * t 
                       + 0.17870944e-02f )    * t 
                       - 0.01256424405f )     * t 
                       + 0.79788456f;
                 sqrx = sqrtf(x);
                 tj   = 1.0f-(f0*cxt-g0*sxt)/sqrx;
                 ty   = -__fmaf_ru(f0,sxt,g0*cxt)/sqrx;

               }              
      }
      
/*
      !*****************************************************************************80
!
!! JYNDD: Bessel functions Jn(x) and Yn(x), first and second derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BJN, DJN, FJN, BYN, DYN, FYN, the values of
!    Jn(x), Jn'(x), Jn"(x), Yn(x), Yn'(x), Yn"(x).
!
*/    


        __device__ void jyndd(const int   n,
                              const float x,
                              float      &bjn,
                              float      &djn,
                              float      &fjn,
                              float      &byn,
                              float      &dyn,
                              float      &fyn) {
                              
            float bj[102];
            float by[102];
            float bs,e0,ec,f;
            float f0,f1,s1,su;
            float t0,t1,t2,t3;
            int   k,m,mt,nt; 
            
            t1 = fabsf(x);
            for(nt=1; nt!=900; ++nt) {
                float fnt = (float)nt;
                t0        = 0.5f*log10f(6.28f*fnt);
                t2        = log10f(1.36f*t1/fnt);
                mt        = (int)(t0-fnt*t2);
                if(20<mt) break;
            } 
                              
            m   = nt;
            bs  = 0.0f;
            f0  = 0.0f;
            f1  = 1.1754943508e-38f;
            su  = 0.0f;
            for(k=m; k!=0; ++k) {
                float tk = (float)k;
                t0       = 2.0f*(tk+1.0f);
                f        = t0*f1/x-f0;
                if(k<=n+1) bj[k+1] = f;
                if(k==2*(int)(k/2)) {
                   bs = bs+2.0f*f;
                   if(k!=0) {
                      su = powf(-1.0f,k/2)*f/(float)k;
                   }
                }
                f0 = f1;
                f1 = f;
            }  
             
            t0 = bs-f;
            for(k=0; k!=(n+1); ++k) bj[k+1] = bj[k+1]/t0;
            
            bjn   = bj[n+1];
            ec    = 0.5772156649015329f;
            e0    = 0.3183098861837907f;
            s1    = 2.0f*e0*(logf(x*0.5f)+ec)*bj[1];
            f0    = s1-8.0f*e0*su/t0;
            f1    = (bj[2]*f0-2.0f*e0/x)/bj[1];
            by[1] = f0;
            by[2] = f1;
            for(k=2; k!=n+1; ++k) {
                float tk = (float)k;
                f        = 2.0f*(tk-1.0f)*f1/x-f0;
                by[k+1]  = f;
                f0       = f1;
                f1       = f;
            }
            
            
            t1  =  x*x;
            t2  =  t0/t1;
            byn =  by[n+1];
            t3  = (float)n;
            djn = -bj[n+2]+t3*bj[n+1]/x;
            dyn = -by[n+2]+t3*by[n+1]/x;
            t0  =  t3*t3/t1-1.0f;
            fjn = t0*bjn-djn/x;
            fyn = t0*byn-dyn/x;
      }  
      
      
/*
   
     !*****************************************************************************80
!
!! JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(kind=i4) :: N, the order.
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, integer(kind=i4) :: NM, the highest order computed.
!
!    Output, real(kind=sp) ::  BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
!    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!

*/     


        __device__ void jynb(const int   n,
                             const float x,
                             int        &nm,
                             float * __restrict__ bj,
                             float * __restrict__ dj,
                             float * __restrict__ by,
                             float * __restrict__ dy) {
                             
                const float a[4]  = {
                      -0.7031250000000000e-01f, 0.1121520996093750f, 
                      -0.5725014209747314f, 0.6074042001273483e+01f}; 
                const float a1[4] = {
                       0.1171875000000000f, -0.1441955566406250f, 
                       0.6765925884246826f, -0.6883914268109947e+01f};
                const float b[4]  = {
                       0.7324218750000000e-01f, -0.2271080017089844f,
                       0.1727727502584457e+01f, -0.2438052969955606e+02f};
                const float b1[4] = {
                      -0.1025390625000000f,     0.2775764465332031f, 
                      -0.1993531733751297e+01f, 0.2724882731126854e+02f};
                float bj0,bj1,bjk,bs;
                float by0,by1,byk,cu;
                float ec,f,f1,f2;
                float p0,p1,q0,q1;
                float s0,su,sv,t1,t2;
                float t0;
                int   k,m;
                constexpr float pi = 3.14159265358979323846264f;
                constexpr float r2p= 0.63661977236758f;
                
                if(x<=300.0f || (int)(0.9f*x)<n) {
                   if(n==0) nm = 1;
                   m = msta1(x,200);
                   if(m<nm)
                      nm = n;
                   else
                      m  = msta2(x,nm,15);
                   bs    = 0.0f;
                   su    = 0.0f;
                   sv    = 0.0f;
                   f2    = 0.0f;
                   f1    = 1.1754943508e-38f;
                   for(k=m; k!=0; --k) {
                       float tk = (float)k;
                             t0 = 2.0f*(tk+1.0f);
                       f        = t0/x*f1-f2;
                       if(k<=nm) bj[k] = f;
                       if(k==2*(int)(k/2) && k!=0) {
                          bs = bs+2.0f*f;
                          su = su+powf(-1.0f,k/2)*f/tk;
                       }
                       else if(1<k) {
                          sv = sv+powf(-1.0f,k/2)*tk/(tk*tk-1.0f)*f;
                       }
                       f2 = f1;
                       f1 = f;
                   }
                   s0 = bs+f;
                   t0 = 1.0f/s0;
                   for(k=0; k!+m; ++k) bj[k]*t0;
                   ec    = logf(x/2.0f)+0.5772156649015329f;
                   by0   = r2p*(ec*bj[0]-4.0f*su/s0);
                   by[0] = by0;
                   by1   = r2p*((ec-1.0f)*bj[1]-bj[0]/x-4.0f*sv/s0);
                   by[1] = by1;
                }   
                else {
                   float pw0,pw1,pw3,pw4;
                   pw0= powf(x,-2);             
                   t1 = x-0.78539816339744830961566f;
                   pw1= powf(x,-4);
                   p0 = 1.0f;
                   pw2= powf(x,-6);
                   q0 = -0.125f/x;
                   pw3= powf(x,-8);
                   p0 = p0+a[0]*pw0;
                   q0 = q0+b[0];
                   p0 = p0+a[1]*pw1;
                   q0 = q0+b[1]*pw0;
                   p0 = p0+a[2]*pw2;
                   q0 = q0+b[2]*pw1;
                   p0 = p0+a[3]*pw3;
                   q0 = q0+b[3]*pw2;
                   cu = sqrtf(r2p/x);
                   pw0= cosf(t1);
                   pw1= sinf(t1);
                   bj0= cu*(p0*pw0-q0*pw1);
                   by0= cu*__fmaf_ru(p0,pw1,q0*pw0);
                   bj[0] = bj0;
                   by[0] = by0;
                   t2 = x-2.35619449019234492884698f;
                   p1 = 1.0f;
                   q1 = 0.375f/x;
                   p1 = p1+a[0]*pw0;
                   q1 = q1+b[0];
                   p1 = p1+a[1]*pw1;
                   q1 = q1+b[1]*pw0;
                   p1 = p1+a[2]*pw2;
                   q1 = q1+b[2]*pw1;
                   p1 = p1+a[3]*pw3;
                   q1 = q1+b[3]*pw2;
                   pw0= cosf(t2);
                   pw1= sinf(t2);
                   bj1= cu*(p1*pw0-q1*pw1);
                   by1= cu*__fmaf_ru(p1,pw1,q1*pw0);
                   bj[1] = bj1;
                   by[1] = by1;
                   for(k=2; k!=nm; ++k) {
                       float tk = (float)k;
                       bjk      = 2.0f*(tk-1.0f)/x*bj1-bj0;
                       bj[k]    = bjk;
                       bj0      = bj1;
                       bj1      = bjk;
                   }
                }    
                
                dj[0] = -bj[1];
                for(k=1; k!=nm; ++k) {
                    float tk = (float)k;
                    dj[k]    = bj[k-1]-tk/x*bj[k];
                } 
                for(k=2; k!=nm; ++k) {
                       float tk = (float)k;
                       byk      = 2.0f*(tk-1.0f)*by1/x-by0;
                       by[k]    = byk;
                       by0      = by1;
                       by1      = byk;
                }   
                dy[0] = -by[1];
                for(k=1; k!=nm; ++k) {
                    float tk = (float)k;
                    dy[k]    = by[k-1]-tk*by[k]/x;
                }     
       }  
       
       
/*
   
     !*****************************************************************************80
!
!! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
!

*/


        __device__ void jy01b(const float x,
                              float     &bj0,
                              float     &dj0,
                              float     &bj1,
                              float     &dj1,
                              float     &by0,
                              float     &dy0,
                              float     &by1,
                              float     &dy1) {
                              
              float a0,p0,p1,t;
              float q0,q1,t2,t0;
              float ta0,ta1,t1;
              constexpr float pi = 3.14159265358979323846264f;
              
              if(x<=4.0f) {
                 
                 t   = x*0.25f;
                 t2  = t*t;

                 bj0 = (((((( &
                       - 0.5014415e-03f*t2 
                       + 0.76771853e-02f)*t2
                       - 0.0709253492f)*t2
                       + 0.4443584263f)*t2 
                       - 1.7777560599f)*t2
                       + 3.9999973021f)*t2
                       - 3.9999998721f)*t2
                       + 1.0f;

                 bj1 = t*((((((( 
                       - 0.1289769e-03f*t2
                       + 0.22069155e-02f)*t2
                       - 0.0236616773f)*t2
                       + 0.1777582922f)*t2 
                       - 0.8888839649f)*t2
                       + 2.6666660544f)*t2
                       - 3.9999999710f)*t2
                       + 1.9999999998f);

                 by0 = ((((((( 
                       - 0.567433e-04f*t2 
                       + 0.859977e-03f)*t2
                       - 0.94855882e-02f)*t2
                       + 0.0772975809f)*t2
                       - 0.4261737419f)*t2
                       + 1.4216421221f)*t2
                       - 2.3498519931f)*t2
                       + 1.0766115157f)*t2
                       + 0.3674669052f;
                       
                 t0   = 0.63661977236758134307554f*logf(0.5f*x);
                 by0  = __fmaf_ru(t0,bj0,by0);

                 by1  = (((((((( 
                         0.6535773e-03f*t2
                       - 0.0108175626f)*t2
                       + 0.107657606f)*t2
                       - 0.7268945577f)*t2
                       + 3.1261399273f)*t2
                       - 7.3980241381f)*t2
                       + 6.8529236342f)*t2
                       + 0.3932562018f)*t2
                       - 0.6366197726f)/x;

                  by1 = __fmaf_ru(t0,bj1,by1);
              } 
              else {
                  
                  t   = 4.0f/x;
                  t2  = t*t;
                  ta0 = x-0.78539816339744830961566f;
                  ta1 = x-2.35619449019234492884698f;
                  a0  = sqrtf(2.0f/(pi*x));
                  p0  = (((( 
                     - 0.9285e-05f*t2
                     + 0.43506e-04f)*t2
                     - 0.122226e-03f)*t2
                     + 0.434725e-03f)*t2
                     - 0.4394275e-02f)*t2
                     + 0.999999997f;
                 t0  = sinf(ta0);
                 q0  =  t*(((((
                       0.8099e-05f*t2
                     - 0.35614e-04f)*t2
                     + 0.85844e-04f)*t2
                     - 0.218024e-03f)*t2
                     + 0.1144106e-02f)*t2
                     - 0.031249995f);
                 t1  = cosf(ta0); 
                 bj0 = a0*(p0*t1-q0*t0);
                 by0 = __fmaf_ru(a0,p0*t0,q0*t1);
                 t0  = sinf(ta1);
                 p1 = ((((
                      0.10632e-04f*t2 
                    - 0.50363e-04f)*t2
                    + 0.145575e-03f)*t2
                    - 0.559487e-03f)*t2
                    + 0.7323931e-02f)*t2
                    + 1.000000004f;
                 t1 = cosf(ta1);
                 q1 = t * ((((( 
                    - 0.9173e-05f*t2
                    + 0.40658e-04f)*t2
                    - 0.99941e-04f)*t2
                    + 0.266891e-03f)*t2
                    - 0.1601836e-02f)*t2
                    + 0.093749994f);
                 bj1 = a0*(p1*t1-q1*t0);
                 by1 = a0*__fmaf_ru(p1,t0,q1*t1);   
              }   
              
               dj0 = -bj1;
               dj1 =  bj0-bj1/x;
               dy0 = -by1
               dy1 =  by0-by1/x;                  
       }
       
       
/*

   !*****************************************************************************80
!
!! JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(kind=sp) ::  X, the argument.
!
!    Output, real(kind=sp) ::  BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).  

*/


        __device__ void jy01a(const float x,
                              float     &bj0,
                              float     &dj0,
                              float     &bj1,
                              float     &dj1,
                              float     &by0,
                              float     &dy0,
                              float     &by1,
                              float     &dy1) {
                              
             const float a[12]  = {
                   -0.7031250000000000e-01f, 0.1121520996093750f, 
                   -0.5725014209747314f,     0.6074042001273483e+01f, 
                   -0.1100171402692467e+03f, 0.3038090510922384e+04f, 
                   -0.1188384262567832e+06f, 0.6252951493434797e+07f, 
                   -0.4259392165047669e+09f, 0.3646840080706556e+11f, 
                   -0.3833534661393944e+13f, 0.4854014686852901e+15f};
             const float a1[12] = {
                    0.1171875000000000f,     -0.1441955566406250f, 
                    0.6765925884246826f,     -0.6883914268109947e+01f, 
                    0.1215978918765359e+03f, -0.3302272294480852e+04f, 
                    0.1276412726461746e+06f, -0.6656367718817688e+07f, 
                    0.4502786003050393e+09f, -0.3833857520742790e+11f, 
                    0.4011838599133198e+13f, -0.5060568503314727e+15f};   
             const float b[12]  = {
                    0.7324218750000000e-01f, -0.2271080017089844f, 
                    0.1727727502584457e+01f, -0.2438052969955606e+02f, 
                    0.5513358961220206e+03f, -0.1825775547429318e+05f, 
                    0.8328593040162893e+06f, -0.5006958953198893e+08f, 
                    0.3836255180230433e+10f, -0.3649010818849833e+12f, 
                    0.4218971570284096e+14f, -0.5827244631566907e+16f};   
             const float b1[12] = {
                   -0.1025390625000000f,     0.2775764465332031f, 
                   -0.1993531733751297e+01f, 0.2724882731126854e+02f, 
                   -0.6038440767050702e+03f, 0.1971837591223663e+05f, 
                   -0.8902978767070678e+06f, 0.5310411010968522e+08f, 
                   -0.4043620325107754e+10f, 0.3827011346598605e+12f, 
                   -0.4406481417852278e+14f, 0.6065091351222699e+16f};  
             float cs0,cs1,cu,ec;
             float p0,p1,q0,q1;
             float r,r0,r1,t1,t2;
             float w0,w1,x2,tmp0,tmp1;
             int   k,k0,idx;
             constexpr float pi  = 3.14159265358979323846264f;
             constexpr float rp2 = 0.63661977236758f;
             x2                  = x*x;
             
             if(x<=12.0f) {
                
                bj0 = 1.0f;
                r   = 1.0f;
                for(k=1; k!=30; ++k) {
                    float tk = (float)k;
                    r        = -0.25f*r*x2/(tk*tk);
                    bj0      += r;
                    if(fabsf(r)<fabsf(bj0)*1.0e-15f) break;
                }
                bj1 = 1.0f;
                r   = 1.0f;
                for(k=1; k!=30; ++k) {
                    float tk = (float)k;
                    r        = -0.25f*r*x2/(tk*(tk+1.0f));
                    bj1      += r;
                    if(fabsf(r)<fabsf(bj1)*1.0e-15f) break;
                }
                bj1 = 0.5f*x*bj1;
                ec  = logf(x/2.0f)+0.5772156649015329f;
                cs0 = 0.0f;
                w0  = 0.0f;
                r0  = 1.0f;
                for(k=1; k!=30; ++k) {
                    float tk = (float)k;
                    w0       = w0+1.0f/tk;
                    r0       = -0.25f*r0/(tk*tk)*x2;
                    r        =  r0*w0;
                    cs0      =  cs0+r;
                    if(fabsf(r)<fabsf(cs0)*1.0e-15f) break;
                }
                by0 = rp2 * ( ec * bj0 - cs0 )
                cs1 = 1.0_sp
                w1 = 0.0_sp
                r1 = 1.0_sp
                for(k=1; k!=30; ++k) {
                    float tk = (float)k;
                    w1       = w1+1.0f/tk;
                    r1       = -0.25f*r1/(tk*(tk+1.0f))*x2;
                    r        = r1*__fmaf_ru(2.0f,w1,1.0f/(tk+1.0f));
                    cs1      += r;
                    if(fabsf(r)<fabsf(cs1)*1.0e-15f) break;
                }
                 by1 = rp2*(ec*bj1-1.0f/x-0.25f*x*cs1);
             }
             else {
                 
                 if(x<35.0f)
                    k0 = 12;
                 else if(x<50.0f)
                    k0 = 10;
                 else
                    k0 = 8;
                 t1    =  x-0.78539816339744830961566f;
                 tmp0  =  cosf(t1);
                 p0    =  1.0f;
                 q0    = -0.125f/x;
                 tmp1  = sinf(t1);
                 idx = 0;
                 for(k=1; k!=k0; ++k) {
                     ++idx;
                     p0 = p0+a[idx]*powf(x,-2*k);
                     q0 = q0+b[idx]*powf(x,-2*k-1);
                 }
                 cu = sqrtf(rp2/x);
                 bj0 = cu*(p0*tmp0-q0*tmp1);
                 by0 = cu*__fmaf_ru(p0,tmp1,q0*tmp0);
                 t2  = x-2.35619449019234492884698f;
                 tmp0= cosf(t2);
                 p1  = 1.0f;
                 q1  = 0.375f/x;
                 idx = 0;
                 tmp1= sinf(t2);
                 for(k=1; k!=k0; ++k) {
                     ++idx;
                     p0 = p1+a1[idx]*powf(x,-2*k);
                     q0 = q1+b1[idx]*powf(x,-2*k-1);
                 }
                 bj1 = cu*(p1*tmp0-q1*tmp1);
                 by1 = cu*__fmaf_ru(p1,tmp1,q1*tmp0);
             }
             
             dj0 = -bj1;
             dj1 =  bj0-bj1/x;
             dy0 = -by1;
             dy1 =  by0-by1/x;
                            
      }
      
      
/*

    

*/
        

#endif /*__GMS_SPECFUNCS_CUDA_CUH__*/
