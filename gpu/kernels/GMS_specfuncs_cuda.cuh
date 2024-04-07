
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
                            float & ai,
                            float & bi,
                            float & ad,
                            float & bd) {
            
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
              bi = xq*(pir*vk1+2.0f/sr3*vi1); //! pir * vk1 + 2.0_sp * invsr3 * vii
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
                             float & ai,
                             float & bi,
                             float & ad,
                             float & bd) {
                             
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
                 xr1 = 1.0_sp / xe;
                 xar = 1.0_sp / xq;
                 xf = sqrtf( xar );
                 rp = 0.5641895835477563f;
                 r = 1.0f;
                 #pragma unroll
                 for(k=1;k!=40;++k) {
                     r = r * ( 6.0_sp * k - 1.0_sp ) * 
                     0.00462962962962962962963f * ( 6.0_sp * k - 3.0_sp ) /
                     k * ( 6.0_sp * k - 5.0_sp ) / ( 2.0_sp * k - 1.0_sp );
                     ck[k] = r;
                     dk[k] = - ( 6.0_sp * k + 1.0_sp ) / ( 6.0_sp * k - 1.0_sp ) * ck[k];
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
                    ai = 0.5_sp * rp * xf * xp1 * sai;
                    bi = rp * xf / xp1 * sbi;
                    ad = -0.5_sp * rp / xf * xp1 * sad;
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
                             float & vj1,
                             float & vj2,
                             float & vy1,
                             float & vy2,
                             float & vi1,
                             float & vi2,
                             float & vk1,
                             float & vk2) {
          
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
                      a0 = powf(0.5_sp * x, vl);
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
                           ( vv - t1*t1 ) / ( k * ( 2.0_sp * k + 1.0_sp ) * x2 );
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







#endif /*__GMS_SPECFUNCS_CUDA_CUH__*/
