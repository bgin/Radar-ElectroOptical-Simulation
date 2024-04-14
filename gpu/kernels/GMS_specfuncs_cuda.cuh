
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
              bi = xq*(pir*vk1+2.0f/sr3*vi1); //! pir * vk1 + 2.0 * invsr3 * vii
              ad = -xa/sr3*pir*vk2;
              bd = xa*(pir*vk2+2.0f/sr3*vi2);
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
                    r  = r*x/(3.0f*k)*x/(3.0f*k-1.0f)*x;
                    fx = fx+r;
                    if ( fabsf(r)<fabsf(fx)*eps) break;
                }
                gx = x;
                r  = x;
                #pragma unroll
                for(k=1;k!=40;++k) {
                    r = r*x/(3.0f*k)*x/(3.0f*k+1.0f)*x;
                    gx = gx+r;
                    if(fabsf(r)<fabsf(gx)*eps) break;
                }
                 ai = c1 * fx - c2 * gx;
                 bi = sr3 * ( c1 * fx + c2 * gx );
                 df = 0.5f * x * x;
                 r = df;
                 #pragma unroll
                for(k=1;k!=40;++k) {
                    r = r * x / ( 3.0f * k ) * x / ( 3.0f * k + 2.0f ) * x;
                    df = df + r; 
                    if (fabsf( r ) < fabsf( df ) * eps ) break; 
                }
                dg = 1.0f;
                r  = 1.0f;
                #pragma unroll
                for(k=1;k!=40;++k) {
                    r = r * x / ( 3.0f * k ) * x / ( 3.0f * k - 2.0f ) * x;
                    dg = dg + r;
                    if ( fabsf( r ) < fabsf( dg ) * eps ) break;
                }
                 ad = c1 * df - c2 * dg;
                 bd = sr3 * ( c1 * df + c2 * dg );
             }
             else {
                 xe = xa * xq / 1.5f;
                 xr1 = 1.0 / xe;
                 xar = 1.0 / xq;
                 xf = sqrtf( xar );
                 rp = 0.5641895835477563f;
                 r = 1.0f;
                 #pragma unroll
                 for(k=1;k!=40;++k) {
                     r = r * ( 6.0 * k - 1.0 ) * 
                     0.00462962962962962962963f * ( 6.0 * k - 3.0 ) /
                     k * ( 6.0 * k - 5.0 ) / ( 2.0 * k - 1.0 );
                     ck[k] = r;
                     dk[k] = - ( 6.0 * k + 1.0 ) / ( 6.0 * k - 1.0 ) * ck[k];
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
                    ai = 0.5 * rp * xf * xp1 * sai;
                    bi = rp * xf / xp1 * sbi;
                    ad = -0.5 * rp / xf * xp1 * sad;
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
                    bi = rp * xf * ( xcs * ssa + xss * ssb );
                    ad = -rp / xf * ( xcs * sda + xss * sdb );
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
                           float t1 = 4.0f * kidx + 1.0f;
                           rq = -0.78125e-2f * rq * ( vv - t0*t0 ) *
                           ( vv - t1*t1 ) / ( k * ( 2.0 * k + 1.0 ) * x2 );
                           qx += rq;
                       }
                       qx = 0.125f * ( vv - 1.0f ) * qx / x;
                       xk = x - ( 0.5f * l * 0.3333333333333333f + 0.25f ) * pi;
                       a0 = sqrtf( rp2 / x );
                       ck = cosf( xk );
                       sk = sinf( xk );
                       if(l == 1) {
                           vj1 = a0 * ( px * ck - qx * sk );
                           vy1 = a0 * ( px * sk + qx * ck );
                       }
                       else {
                           vj2 = a0 * ( px * ck - qx * sk );
                           vy2 = a0 * ( px * sk + qx * ck )
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
                        gr = gr*z+gk;
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

#endif /*__GMS_SPECFUNCS_CUDA_CUH__*/
