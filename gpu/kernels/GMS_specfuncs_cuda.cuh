
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
          
             

#endif /*__GMS_SPECFUNCS_CUDA_CUH__*/
