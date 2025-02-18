
#ifndef __GMS_MATRIX_COMPUTATIONS_H__
#define __GMS_MATRIX_COMPUTATIONS_H__


namespace file_info {

   const unsigned int gGMS_MATRIX_COMPUTATIONS_MAJOR = 1U;
   const unsigned int gGMS_MATRIX_COMPUTATIONS_MINOR = 0U;
   const unsigned int gGMS_MATRIX_COMPUTATIONS_MICRO = 0U;
   const unsigned int gGMS_MATRIX_COMPUTATIONS_FULLVER =
     1000U*gGMS_MATRIX_COMPUTATIONS_MAJOR+100U*gGMS_MATRIX_COMPUTATIONS_MINOR+
     10U*gGMS_MATRIX_COMPUTATIONS_MICRO;
   const char * const pgGMS_MATRIX_COMPUTATIONS_CREATE_DATE = "17-04-2020 11:10 +00200 (FRI 17 APR 2020 GMT+2)";
   const char * const pgGMS_MATRIX_COMPUTATIONS_BUILD_DATE  = __DATE__ ":" __TIME__;
   const char * const pgGMS_MATRIX_COMPUTATIONS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const pgGMS_MATRIX_COMPUTATIONS_SYNOPSIS    = "Matrix real and complex helper functions";
}

#include <cstdint>
#include <complex>
#include "GMS_config.h"
#include "GMS_indices.h"

namespace gms {

       namespace math {

               /*
                    Matrix exponential computation --  unrolled version 1
                */
		__ATTR_HOT__
		__ATTR_ALIGN__(16)
		
		void exp4x4m_cmplxr4v1( const std::complex<float> * __restrict __ATTR_ALIGN__(64) L,
		                      const std::complex<float> * __restrict __ATTR_ALIGN__(64) Q,
				      const std::complex<float> * __restrict __ATTR_ALIGN__(64) INVQ,
				      const float z,
				      std::complex<float> * __restrict __ATTR_ALIGN__(64) result); 

		/*
                        Complex 4x4 matrix multiplication (single precision)
                 */
		 __ATTR_HOT__
		 __ATTR_ALIGN__(16)
		
		 void mul4x4m_cmplxr4(const std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		                      const std::complex<float> * __restrict __ATTR_ALIGN__(64) b,
				      std::complex<float> * __restrict __ATTR_ALIGN__(64) c); 
             /*
                  Helper function to multiply 3 complex matrices 4x4
              */
	      __ATTR_ALWAYS_INLINE__
	      static inline
	      void mul4x4m_cmplxr4_helper(const std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
	                                  const std::complex<float> * __restrict __ATTR_ALIGN__(64) b,
					  const std::complex<float> * __restrict __ATTR_ALIGN__(64) c,
					  std::complex<float> * __restrict __ATTR_ALIGN__(64) result) {
                  __ATTR_ALIGN__(64) std::complex<float> d[16];
		  mul4x4m_cmplxr4(a,b,d);
		  mul4x4m_cmplxr4(d,c,result);
	       }
             /*
                 The exponential of complex matrix 4x4 version 2
              */
	    __ATTR_HOT__
	    __ATTR_ALIGN__(16)
	  
	    void exp4x4m_cmplxr4v2(const std::complex<float> * __restrict __ATTR_ALIGN__(64) L,
	                           const std::complex<float> * __restrict __ATTR_ALIGN__(64) Q,
				   const std::complex<float> * __restrict __ATTR_ALIGN__(64) INVQ,
				   const float z,
				  float * __restrict __ATTR_ALIGN__(64) result); 
		   /*
                        4x4 real matrix extinction
                    */
		   __ATTR_ALWAYS_INLINE__
		   static inline
		   void extinct_m4x4r4(const std::complex<float> * __restrict __ATTR_ALIGN__(32) M,
		                       float * __restrict __ATTR_ALIGN__(64) K) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                             M = (const std::complex<float>*)__builtin_assume_aligned(M,32);
			     K = (float*)__builtin_assume_aligned(K,64);
#elif defined __ICC || defined __INTEL_COMPILER
                             __assume_aligned(M,32);
			     __assume_aligned(K,64);
#endif
                             K[0] = -2.0f*M[0].real();
			     K[1] = 0.0f;
			     K[2] = -2.0f*M[1].real();
			     K[3] = -2.0f*M[1].imag();
			     K[4] = 0.0f;
			     K[5] = -2.0f*M[3].real();
			     K[6] = -2.0f*M[2].real();
			     K[7] = -2.0f*M[2].imag();
			     K[8] = -M[2].real();
			     K[9] = -M[1].real();
			     K[10] = -(M[0].real()+M[3].real());
			     K[11] = -(M[0].imag()-M[1].imag());
			     K[12] = -2.0f*M[2].real();
			     K[13] = -2.0f*M[1].real();
			     K[14] = -(M[0].imag()-M[1].imag());
			     K[15] = -(M[0].real()+M[3].real());
		    }

		   /*
                         Eigenvalue of complex 4x4 matrix
                    */
		   
		   __ATTR_HOT__
		   __ATTR_ALIGN__(16)
		   void eigen4x4_cmplxr4(const std::complex<float> * __restrict __ATTR_ALIGN__(64) m,
		                         std::complex<float> * __restrict __ATTR_ALIGN__(64) l,
					 std::complex<float> * __restrict __ATTR_ALIGN__(64) q,
					 std::complex<float> * __restrict __ATTR_ALIGN__(64) invq,
					 const float freq,
					 const float wlength,
					 const float k0); 
		   /*
                         Complex matrix 4x4 inversion (using unpotimized Linpack routines
                         'cgeco' and 'cgedi'
                    */
                   __ATTR_ALWAYS_INLINE__
		   static inline
		   void invm4x4_cmplxr4(const std::complex<float> * __restrict __ATTR_ALIGN__(64) in,
		                        std::complex<float> * __restrict __ATTR_ALIGN__(64) out) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        in  =  (const std::complex<float>*)__builtin_assume_aligned(in,64);
			out =  (std::complex<float>*)__builtin_assume_aligned(out,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(in,64);
			__assume_aligned(out,64);
#endif
                        __ATTR_ALIGN__(32) std::complex<float> sv[4] = {};
			__ATTR_ALIGN__(16) std::complex<float> d[2]  = {};
			__ATTR_ALIGN__(16) int32_t iv[4] = {};
			float rc;
			int32_t job;
			// Small loop manually unrolled
			out[0] = in[0];
		        out[1] = in[1];
			out[2] = in[2];
			out[3] = in[3];
			out[4] = in[4];
			out[5] = in[5];
			out[6] = in[6];
			out[7] = in[7];
			out[8] = in[8];
			out[9] = in[9];
			out[10] = in[10];
			out[11] = in[11];
			out[12] = in[12];
			out[13] = in[13];
			out[14] = in[14];
			out[15] = in[15];
			cgeco(&out[0],4,4,&ipvt[0]);
			job = 1;
			cgedi(&out[0],4,4,&ipvt[0],&det[0],job);
		   }
		   

		  
//****************************************************************************80
//
//  Purpose:
//
//    CGECO factors a complex matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, CGEFA is slightly faster.
//
//    To solve A*X = B, follow CGECO by CGESL.
//
//    To compute inverse(A)*C, follow CGECO by CGESL.
//
//    To compute determinant(A), follow CGECO by CGEDI.
//
//    To compute inverse(A), follow CGECO by CGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.

                  
	           __ATTR_HOT__
		   __ATTR_ALIGN__(16)
	           float cgeco(std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
	                       const int32_t lda,
		               const int32_t n,
		               int32_t * __restrict __ATTR_ALIGN__(64) ipvt); 
//  Purpose:
//
//    CGEDI computes the determinant and inverse of a matrix.
//
//  Discussion:
//
//    The matrix must have been factored by CGECO or CGEFA.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if CGECO has set 0.0 < RCOND or CGEFA has set
//    INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science		 

                   __ATTR_HOT__
		   __ATTR_ALIGN__(16)
		
		   void cgedi(std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		              const int32_t lda,
			      const int32_t n,
			      int32_t * __restrict __ATTR_ALIGN__(64) ipvt,
			      std::complex<float> det[2],
			      const int32_t job);
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//

                  __ATTR_HOT__
		  __ATTR_ALIGN__(64)
		 
                  int32_t cgefa(std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		                const int32_t lda,
				const int32_t n,
				int32_t * __restrict __ATTR_ALIGN__(64) ipvt); 
		     
		  
		  
		   static inline
		   float r4_max(const float x,
		                const float y) {
                       if(y < x) {
                           return (x);
		       }
		       else {
                           return (y);
		       }
		    }
//****************************************************************************80
//
//  Purpose:
//
//    SCASUM takes the sum of the absolute values of a vector.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
		 
		   static inline
		   float scasum(const int32_t n,
		                const std::complex<float> * __restrict __ATTR_ALIGN__(64) x,
				const int32_t incx) {
#if defined __GNUC__ || !defined __INTEL_COMPILER
                        x = (const std::complex<float>*)__builtin_assume_aligned(x,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(x,64);
#endif
                        float value;
			value = 0.0f;
			if ( n <= 0 || incx <= 0 )
                           {
                              return value;
                        }
                        if ( incx == 1 )
			
                           {
                               for ( i = 0; i < n; i++ )
                                   {
                                     value = value + std::fabs ( x[i].real() )
                                                   + std::fabs ( x[i].imag() );
                                   }
                           }
                        else
                              {
                                   ix = 0;
                                   for ( i = 0; i < n; i++ )
                                       {
                                          value = value + std::fabs (  x[ix].real() )
                                                        + std::fabs (  x[ix].imag() );
                                          ix = ix + incx;
                                       }
                              }
                                return value;
		    }

//Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.		    

		 
		    static inline
		    float cabs1(const std::complex<float> z) {
                        float value;
			value = 0.0f;
			value = std::fabs(z.real())+std::fabs(z.imag());
			return (value);
		    }

// Modified:
//
//    11 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.

                 
		   static inline
		   void cswap(const int32_t n,
		              std::complex<float> * __restrict __ATTR_ALIGN__(64) cx,
			      const int32_t incx,
			      std::complex<float> * __restrict __ATTR_ALIGN__(64) cy,
			      const int32_t incy) {
			      if(n <= 0) {
                                  return;
                              }
#if defined __GNUC__ && !defined __INTEL_COMPILER
                              cx = (std::complex<float>*)__builtin_assume_aligned(cx,64);
			      cy = (std::complex<float>*)__builtin_assume_aligned(cy,64);
#elif defined __ICC || defined __INTEL_COMPILER
                              __assume_aligned(cx,64);
			      __assume_aligned(cy,64);
#endif
                             std::complex<float>ctemp;
                             int i;
                             int ix;
                             int iy;
                             if(incx == 1 && incy == 1) {
                                for(i = 0; i < n; i++) {
                                    ctemp = cx[i];
                                    cx[i] = cy[i];
                                    cy[i] = ctemp;
                                 }
                             }
                              else {
                                   if(0 <= incx) {
                                      ix = 0;
                                   }
                                   else {
                                      ix = (-n + 1 )*incx;
                                   }
                                   if(0 <= incy) {
                                      iy = 0;
                                   }
                                   else {
                                   iy = ( -n + 1 ) * incy;
                                  }
                              for(i = 0; i < n; i++) {
                                  ctemp = cx[ix];
                                  cx[ix] = cy[iy];
                                  cy[iy] = ctemp;
                                  ix = ix + incx;
                                  iy = iy + incy;
                              }
                         } 
		   }

//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.		    

                    
		    static inline
		    int32_t icamax(const int32_t n,
		                   std::complex<float> * __restrict __ATTR_ALIGN__(64) x,
				   const int32_t incx) {
			   if( n < 1 || incx  <=  0 ) {
                                 return (0);
                           }
                           if ( n == 1 ) {
                                 return (1);
                           }
#if defined __GNUC__ && !defined __INTEL_COMPILER
                          x = (std::complex<float>*)__builtin_assume_aligned(x,64);
#elif defined __ICC || defined __INTEL_COMPILER
                          __assume_aligned(x,64);
#endif
                          int i;
                          int ix;
                          float smax;
                          int value;

                          value = 0;

                          if ( incx != 1 ){
  
                              ix = 0;
                              smax = cabs1 ( x[0] );
                              ix = ix + incx;
                              for(i = 1; i < n; i++ ){
                                    const float t0 = cabs1(x[ix]);
                                    if ( smax < t0 ){
                                         value = i + 1;
                                         smax = t0;
                                    }
                                    ix = ix + incx;
                              }
                         }
                          else
                              {
                               smax = cabs1 ( x[0] );
                               for ( i = 1; i < n; i++ ) {
			            const float t0 = cabs1(x[i]);
                                    if ( smax < t0 ){
                                         value = i + 1;
                                         smax = t0;
                                     }
                               }
                           }
                           return value;
		      }


//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> Z1, Z2, the arguments.
//
//    Output, complex <float> CSIGN1,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
                 
		   static inline
		   std::complex<float> csign1(const std::complex<float> z1,
		                              const std::complex<float> z2) {
                        if(cabs1(z2) == 0.0f) {
                           return (std::complex<float>(0.0f,0.0f));
			}
			else {
                           return (cabs1(z1)*(z2/cabs1(z2)));
			}
		   }

//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
                 
		   static inline
		   void csscal(const int32_t n,
		               const float sa,
			       std::complex<float> * __restrict __ATTR_ALIGN__(64) cx,
			       const int32_t incx) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                               cx = (std::complex<float>*)__builtin_assume_aligned(cx,64);
#elif defined __ICC || defined __INTEL_COMPILER
                               __assume_aligned(cx,64);
#endif
                               if(n <= 0 || incx <= 0) {
                                  return;
			       }
			      if( incx == 1) {
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(64)
#pragma vector
#endif
                                  for(int32_t i = 0; i != n; ++i) {
                                      const std::complex<float> t0 = sa * cx[i];
				      cx[i] = t0;
				  }
			       }
			       else {
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(64)
#pragma vector
#endif
                                  for(int32_t i = 0; i != n; ++i) {
                                      const std::complex<float> t0 = sa * cx[i*incx];
				      cx[i*incx] = t0;
				  }
			       }
		   }

//
//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.

                 
		   static inline
		   std::complex<float> cdotc(const int32_t n,
		                             const std::complex<float> * __restrict __ATTR_ALIGN__(64) cx,
					     const int32_t incx,
					     const std::complex<float> * __restrict __ATTR_ALIGN__(64) cy,
					     const int32_t incy) {
			if(n <= 0) {
                           return (std::complex<float>(0.0f,0.0f));
			}
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        cx = (const std::complex<float>*)__builtin_assume_aligned(cx,64);
			cy = (const std::complex<float>*)__builtin_assume_aligned(cy,64);
#elif defined __ICC || defined __INTEL_COMPILER
			__assume_aligned(cx,64);
			__assume_aligned(cy,64);
#endif
			std::complex<float> value;
			int32_t ix,iy;
			value = {0.0f,0.0f};
			if(incx == 1 && incy == 1) {
#pragma code_align(64)
#pragma vector always
                           for(int32_t i = 0; i != n; ++i) {
                               value = value + std::conj(cx[i]) * cy[i];
			   }
			}
			else {
                           if(0 <= incx) {
                              ix = 0;
			   }
			   else {
                              ix = (-n + 1) * incx;
			   }
			   if(0 <= incy) {
                              iy = 0;
			   }
			   else {
                              iy = (-n + 1) * incy;
			   }
			   for(int32_t  i = 0; i != n; ++i) {
                               value = value + std::conj(cx[ix]) * cy[iy];
			       ix += incx;
			       iy += incy
			   }
			}
			return (value);
		   }

//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.

                 
		   static inline
		   void caxpy(const int32_t n,
		              const std::complex<float> ca,
			      const std::complex<float> * __restrict __ATTR_ALIGN__(64) cx,
			      const int32_t incx,
			      std::complex<float> * __restrict __ATTR_ALIGN__(64) cy,
			      const int32_t incy) {
			      if(n <= 0) {
                                 return;
			      }
			      if(cabs1(ca) == 0.0f) {
                                 return;
			      }
#if defined __GNUC__ && !defined __INTEL_COMPILER
		             cx = (const std::complex<float>*)__builtin_assume_aligned(cx,64);
			     cy = (std::complex<float>*)__builtin_assume_aligned(cy,64);
#elif defined __ICC || defined __INTEL_COMPILER
                             assume_aligned(cx,64);
			     assume_aligned(cy,64);
#endif
                             int32_t i,ix,iy;
			     if(incx != 1 || incy != 1) {
                                if(incx <= 0) {
                                   ix = 0;
				}
				else {
                                   ix = (-n+1)*incx;
				}
				if(incy <= 0) {
                                   iy = 0;
				}
				else {
                                   iy = (-n+1)*incy;
				}
#pragma code_align(64)
#pragma vector always
				for(i = 0; i != n; ++i) {
                                    std::complex<float> t0 = cy[iy] + ca * cx[ix];
				    cy[iy] = t0;
				    ix += incx;
				    iy += incy;
				}
			     }
			     else {
#pragma code_align(64)
#pragma vector always
                                for(i = 0; i != n; +++i) {
                                    std::complex<float> t0 = cy[i] + ca * cx[i];
				    cy[i] = t0;
				}
			     }
		   }

   } // math

} // gms








#endif /*__GMS_MATRIX_COMPUTATIONS_H__*/
