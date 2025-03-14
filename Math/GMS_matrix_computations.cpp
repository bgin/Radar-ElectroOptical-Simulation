

#include "GMS_matrix_computations.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"



               /*
                    Matrix exponential computation --  unrolled version 1
                */
		
		void gms::math::exp4x4m_cmplxr4v1( const std::complex<float> * __restrict __ATTR_ALIGN__(64) L,
		                      const std::complex<float> * __restrict __ATTR_ALIGN__(64) Q,
				      const std::complex<float> * __restrict __ATTR_ALIGN__(64) INVQ,
				      const float z,
				      std::complex<float> * __restrict __ATTR_ALIGN__(64) result) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        L      = (const std::complex<float>*)__builtin_assume_aligned(L,64);
			Q      = (const std::complex<float>*)__builtin_assume_aligned(Q,64);
			INVQ   = (const std::complex<float>*)__builtin_assume_aligned(INVQ,64);
			result = (std::complex<float>*)__builtin_assume_aligned(result,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(L,64);
			__assume_aligned(Q,64);
			__assume_aligned(INVQ,64);
			__assume_aligned(result,64);
#endif
                        __ATTR_ALIGNED__(64) std::complex<float> t0[16] = {};
			t0[0]  = std::exp(L[0]*z);
			t0[5]  = std::exp(L[1]*z);
			t0[10] = std::exp(L[2]*z);
			t0[15] = std::exp(L[3]*z);
	        }

		/*
                        Complex 4x4 matrix multiplication (single precision)
                 */
		
		 void gms::math::mul4x4m_cmplxr4(const std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		                      const std::complex<float> * __restrict __ATTR_ALIGN__(64) b,
				      std::complex<float> * __restrict __ATTR_ALIGN__(64) c) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                      a   = (const std::complex<float>*)__builtin_assume_aligned(a,64);
		      b   = (const std::complex<float>*)__builtin_assume_aligned(b,64);
		      c   = (std::complex<float>*)__builtin_assume_aligned(c,64);
#elif defined __ICC || defined __INTEL_COMPILER
                      __assume_aligned(a,64);
		      __assume_aligned(b,64);
		      __assume_aligned(c,64);
#endif
                      c = {};
                      for(int32_t i = 0; i != 4; i += 2) {
                          for(int32_t j = 0; j != 4; j += 2) {
			      std::complex<float> s0 = {0.0f,0.0f};
			      std::complex<float> s1 = {0.0f,0.0f};
			      std::complex<float> s2 = {0.0f,0.0f};
			      std::complex<float> s3 = {0.0f,0.0f};
			      for(int32_t k = 0; k != 4; ++k) {
                                  std::complex<float> r1 = a[Ix2D(j,4,k)];
				  std::complex<float> r2 = b[Ix2D(k,4,i)];
				  s0 = s0 + r1 * r2;
				  std::complex<float> r3 = a[Ix2D(j+1,4,k)];
				  s1 = s1 + r3 * r2;
				  std::complex<float> r4 = b[Ix2D(k+1,4,i)];
				  s2 = s2 + r1 * r4;
				  s3 = s3 + r3 * r4;
			      }
			      c[Ix2D(i,4,j)]     = s0;
			      c[Ix2D(i+1,4,j)]   = s1;
			      c[Ix2D(i,4,j+1)]   = s2;
			      c[Ix2D(i+1,4,j+1)] = s3;
			  }
		      }
                 }
            
             /*
                 The exponential of complex matrix 4x4 version 2
              */
	   
	    void gms::math::exp4x4m_cmplxr4v2(const std::complex<float> * __restrict __ATTR_ALIGN__(64) L,
	                           const std::complex<float> * __restrict __ATTR_ALIGN__(64) Q,
				   const std::complex<float> * __restrict __ATTR_ALIGN__(64) INVQ,
				   const float z,
				  float * __restrict __ATTR_ALIGN__(64) result) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        L      = (const std::complex<float>*)__builtin_assume_aligned(L,64);
			Q      = (const std::complex<float>*)__builtin_assume_aligned(Q,64);
			INVQ   = (const std::complex<float>*)__builtin_assume_aligned(INVQ,64);
			result = (float*)__builtin_assume_aligned(result,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(L,64);
			__assume_aligned(Q,64);
			__assume_aligned(INVQ,64);
			__assume_aligned(result,64);
#endif
#if defined __AVX512F__
			__ATTR_ALIGN__(64) std::complex<float> Diag[16] = {};
                        __ATTR_ALIGN__(64) float Qre[16]    = {};
			__ATTR_ALIGN__(64) float Qim[16]    = {};
			__ATTR_ALIGN__(64) float INVQre[16] = {};
			__ATTR_ALIGN__(64) float INVQim[16] = {};
			__ATTR_ALIGN__(64) float Lre[16]    = {};
			__ATTR_ALIGN__(64) float Lim[16]    = {};
			__ATTR_ALIGN__(64) float FRre[16]   = {};
			__ATTR_ALIGN__(64) float FRim[16]   = {};
			__ATTR_ALIGN__(64) float FIre[16]   = {};
			__ATTR_ALIGN__(64) float FIim[16]   = {};
			__ATTR_ALIGN__(64) float Diff[16]   = {};
			__ATTR_ALIGN__(16) float Sum[16]    = {};
#else
			__ATTR_ALIGN__(32) std::complex<float> Diag[16] = {};
                        __ATTR_ALIGN__(32) float Qre[16]    = {};
			__ATTR_ALIGN__(32) float Qim[16]    = {};
			__ATTR_ALIGN__(32) float INVQre[16] = {};
			__ATTR_ALIGN__(32) float INVQim[16] = {};
			__ATTR_ALIGN__(32) float Lre[16]    = {};
			__ATTR_ALIGN__(32) float Lim[16]    = {};
			__ATTR_ALIGN__(32) float FRre[16]   = {};
			__ATTR_ALIGN__(32) float FRim[16]   = {};
			__ATTR_ALIGN__(32) float FIre[16]   = {};
			__ATTR_ALIGN__(32) float FIim[16]   = {};
			__ATTR_ALIGN__(32) float Diff[16]   = {};
			__ATTR_ALIGN__(32) float Sum[16]    = {};
#endif
			for(int32_t i = 0; i != 16; ++i) {
                            Qre[i]    = Q[i].real();
			    Qim[i]    = Q[i].imag();
			    INVQre[i] = INVQ[i].real();
			    INVQim[i] = INVQ[i].imag();
			}
                        Diag[0] = std::exp(L[0]*z);
			Lre[0]  = Diag[0].real();
			Lim[0]  = Diag[0].imag();
			Diag[5] = std::exp(L[1]*z);
			Lre[5]  = Diag[5].real();
			Lim[5]  = Diag[5].imag();
			Diag[10] = std::exp(L[10]*z);
			Lre[10]  = Diag[10].real();
			Lim[10]  = Diag[10].imag();
			Diag[15] = std::exp(L[15]*z);
			Lre[15]  = Diag[15].real();
			Lim[15]  = Diag[15].imag();
#if defined __AVX512F__
                      
			__m512 zmm0,zmm1;
		      
			zmm0 = _mm512_setr_ps(Lre[15],Lre[10],Lre[5],Lre[0],
			                      Lre[15],Lre[10],Lre[5],Lre[0],
					      Lre[15],Lre[10],Lre[5],Lre[0],
					      Lre[15],Lre[10],Lre[5],Lre[0]);
			zmm1 = _mm512_setr_ps(Lim[15],Lim[10],Lim[5],Lim[0],
			                      Lim[15],Lim[10],Lim[5],Lim[0],
					      Lim[15],Lim[10],Lim[5],Lim[0],
					      Lim[15],Lim[10],Lim[5],Lim[0]);
			const __m512 qre = _mm512_load_ps(&INVQre[0]);
			const __m512 qim = _mm512_load_ps(&INVQim[0]);
			_mm512_store_ps(&FRre[0],_mm512_mul_ps(qre,zmm0));
			_mm512_store_ps(&FRim[0],_mm512_mul_ps(qim,zmm1));
			_mm512_store_ps(&FIre[0],_mm512_mul_ps(qre,zmm0));
			_mm512_store_ps(&FIim[0],_mm512_mul_ps(qim,zmm1));
			_mm512_store_ps(&Diff[0],_mm512_sub_ps(_mm512_load_ps(&FRre[0]),
			                                       _mm512_load_ps(&FIim[0])));
			_mm512_store_ps(&Sum[0], _mm512_add_ps(_mm512_load_ps(&FRim[0]),
			                                       _mm512_load_ps(&FIre[0])));
#else
                       	__m256 ymm0,ymm1;
			ymm0 = _mm256_setr_ps(Lre[15],Lre[10],Lre[5],Lre[0],
					      Lre[15],Lre[10],Lre[5],Lre[0]);
			ymm1 = _mm256_setr_ps(Lim[15],Lim[10],Lim[5],Lim[0],
					      Lim[15],Lim[10],Lim[5],Lim[0]);
		   
			const __m256 qrel = _mm256_load_ps(&QINVre[0]);
			const __m256 qreh = _mm256_load_ps(&QINVre[8]);
			const __m256 qiml = _mm256_load_ps(&QINVim[0]);
			const __m256 qimh = _mm256_load_ps(&QINVim[8]);
			_mm256_store_ps(&FRre[0],_mm256_mul_ps(qrel,ymm0));
			_mm256_store_ps(&FRre[8],_mm256_mul_ps(qreh,ymm0));
			_mm256_store_ps(&FRim[0],_mm256_mul_ps(qiml,ymm1));
			_mm256_store_ps(&FRim[8],_mm256_mul_ps(qimh,ymm1));
			_mm256_store_ps(&FIre[0],_mm256_mul_ps(qrel,ymm0));
			_mm256_store_ps(&FIre[8],_mm256_mul_ps(qreh,ymm0));
			_mm256_store_ps(&FIim[0],_mm256_mul_ps(qiml,ymm1));
			_mm256_store_ps(&FIim[8],_mm256_mul_ps(qimh,ymm1));
			_mm256_store_ps(&Diff[0],_mm256_sub_ps(_mm256_load_ps(&FRre[0]),
							       _mm256_load_ps(&FIim[0])));
			_mm256_store_ps(&Sum[0], _mm256_add_ps(_mm256_load_ps(&FRim[0]),
							       _mm256_load_ps(&FIre[0])));
			_mm256_store_ps(&Diff[8],_mm256_sub_ps(_mm256_load_ps(&FRre[8]),
							       _mm256_load_ps(&FIim[8])));
			_mm256_store_ps(&Sum[8], _mm256_add_ps(_mm256_load_ps(&FRim[8]),
							       _mm256_load_ps(&FIre[8])));
#endif
			// This should be vectorized.
			result = {};
		        for(int32_t i = 0; i != 4; i += 2) {
                           for(int32_t j = 0; j != 4; j += 2) {
                              float s0 = 0.0f;
			      float s1 = 0.0f;
			      float s2 = 0.0f;
			      float s3 = 0.0f;
			      for(int32_t k = 0; k != 4; ++k) {
                                  float r1 = Qre[Ix2D(j,4,k)];
				  float r2 = Diff[Ix2D(k,4,i)];
				  float r3 = Qim[Ix2D(j,4,k)];
				  float r4 = Sum[Ix2D(k,4,i)];
				  s0 = s0 + r1*r2-r3*r4;
				  float r5 = Qre[Ix2D(j+1,4,k)];
				  float r6 = Diff[Ix2D(k,4,i+1)];
				  float r7 = Qim[Ix2D(j+1,4,k)];
				  s1 = s1 + r5*r2-r7*r4;
				  float r8 = Sum[Ix2D(k,4,i+1)];
				  s2 = s2 + r1*r6-r3*r8;
				  s3 = s3 + r5*r6-r7*r8;
			      }
			      result[Ix2D(i,4,j)]     = s0;
			      result[Ix2D(i+1,4,j)]   = s1;
			      result[Ix2D(i,4,j+1)]   = s2;
			      result[Ix2D(i+1,4,j+1)] = s3; 
			   }
			}
		
		     }

		
		   /*
                         Eigenvalue of complex 4x4 matrix
                    */
		   
		   void gms::math::eigen4x4_cmplxr4(const std::complex<float> * __restrict __ATTR_ALIGN__(64) m,
		                         std::complex<float> * __restrict __ATTR_ALIGN__(64) l,
					 std::complex<float> * __restrict __ATTR_ALIGN__(64) q,
					 std::complex<float> * __restrict __ATTR_ALIGN__(64) invq,
					 const float freq,
					 const float wlength,
					 const float k0) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        m = (const std::complex<float>*)__builtin_assume_aligned(m,64);
			l = (std::complex<float>*)__builtin_assume_aligned(l,64);
			q = (std::complex<float>*)__builtin_assume_aligned(q,64);
			invq = (std::complex<float>*)__builtin_assume_aligned(invq,64);
#elif defined __ICC || defined __INTEL_COMPILER
                        __assume_aligned(m,64);
			__assume_aligned(l,64);
			__assume_aligned(q,64);
			__assume_aligned(invq,64);
#endif
                       std::complex<float> evv[2];
		       std::complex<float> evh[2];
		       std::complex<float> ev;
		       std::complex<float> eh;
                       std::complex<float> j;
		       std::complex<float> cw1;
		       std::complex<float> cw2;
		       std::complex<float> diff;
		       std::complex<float> sum;
		       std::complex<float> t0;
		       std::complex<float> t1;
		       std::complex<float> t2;
		       float tol,w1,w2,c1,c2;
		       bool b1,b2,b3,b4,b5,b6;
		       j = {0.0f,1.0f};
		       tol = 0.001f;
		       b1 = false;
		       b1 = std::cabs(m[2]) == 0.0f;
		       b2 = false;
		       b2 = std::cabs(m[0]) != 0.0f;
		       if(b1) {
                          c1 = 0.0f;
		       }else if(b2) {
                          c1 = std::cabs(m[2]/m[0]);
		       }else {
                          c1 = tol + 1.0f;
		       }
		       b3 = false;
		       b3 = std::cabs(m[1]) == 0.0f;
		       b4 = false;
		       b4 = std::cabs(m[3]) != 0.0f;
		       if(b3) {
                          c2 = 0.0f;
		       }else if(b4) {
                          c2 = std::cabs(m[1]/m[3]);
		       }else {
                          c2 = tol + 1.0f;
		       }
		       b5 = c1 < tol;
		       b6 = c2 < tol;
		       if(b5 && b6) {
                          ev = k0 - j*m[0];
			  eh = k0 - j*m[3];
			  evv[0] = {1.0f,0.0f};
			  evv[1] = {0.0f,0.0f};
			  evh[0] = {0.0f,0.0f};
			  evh[1] = {1.0f,0.0f};
			  cw1 = -(m[0]+m[3]).real();
			  cw2 =  (m[0]+m[3]).imag();
			  l[0] = -2.0f*m[0].real();
			  l[1] = {cw1,-cw2};
			  l[2] = {cw1,cw2};
			  l[3]  = -2.0f*m[3].real();
			  q[0]  = {1.0f,0.0f};
			  q[1]  = {0.0f,0.0f};
			  q[2]  = {0.0f,0.0f};
			  q[3]  = {0.0f,0.0f};
			  q[4]  = {0.0f,0.0f};
			  q[5]  = {0.0f,0.0f};
			  q[6]  = {1.0f,0.0f};
			  q[7]  = {0.0f,-1.0f};
			  q[8]  = {0.0f,0.0f};
			  q[9]  = {0.0f,0.0f};
			  q[10] = {1.0f,0.0f};
			  q[11] = {0.0f,1.0f};
			  q[12] = {0.0f,0.0f};
			  q[13] = {1.0f,0.0f};
			  q[14] = {0.0f,0.0f};
			  q[15] = {0.0f,0.0f};

			  invq[0]  = {1.0f,0.0f};
			  invq[1]  = {0.0f,0.0f};
			  invq[2]  = {0.0f,0.0f};
			  invq[3]  = {0.0f,0.0f};
			  invq[4]  = {0.0f,0.0f};
			  invq[5]  = {0.0f,0.0f};
			  invq[6]  = {0.0f,0.0f};
			  invq[7]  = {1.0f,0.0f};
			  invq[8]  = {0.0f,0.0f};
			  invq[9]  = {0.5f,0.0f};
			  invq[10] = {0.5f,0.0f};
			  invq[11] = {0.0f,0.0f};
			  invq[12] = {0.0f,0.0f};
			  invq[13] = {0.0f,0.5f};
			  invq[14] = {0.0f,-0.5f};
			  invq[15] = {0.0f,0.0f};
		       }
		       else {
                              diff = m[0]-m[3];
			      sum  = m[0]+m[3];
			      std::complex<float> arg0 = 4.0f*m[1]*m[2];
			      t0 = std::sqrt(diff*diff+arg0);
			      ev = k0-(j/2.0f)*(sum+t0);
			      eh = k0-(j/2.0f)*(sum-t0);
			      t1 = 2.0f*m[1]/(diff+t0);
			      t2 = 2.0f*m[2]/(-diff-t0);
			      ev = k0-j*(m[0]+m[3]+r)/2.0f;
			      eh = k0-j*(m[0]+m[3]-r)/2.0f;
			      evv[0] = {1.0f,0.0f};
			      evv[1] = ev;
			      evh[0] = eh;
			      evh[1] = {1.0f,0.0f};
			      l[0] = 2.0f*ev.real();
			      l[1] = j*(std::conj(eh)-ev);
			      l[2] = j*(std::conj(ev)-eh);
			      l[3] = 2.0*eh.imag();
			      w1 = std::abs(ev);
			      w2 = std::abs(eh);
			      cw1 = std::conj(ev);
			      cw2 = std::conj(eh);
			      q[0] = {1.0f,0.0f};
			      q[1] = w1*w2;
			      q[2] = 2.0f*t1.real();
			      q[3] = -2.0f*t1.imag();
			      q[4] = cw2;
			      q[5] = t1;
			      q[6] = 1.0f+t1*cw2;
			      q[7] = -j*(1.0f-t1*cw2);
			      q[8] = t2;
			      q[9] = cw1;
			      q[10] = 1.0f+t2*cw1;
			      q[11] = j*(1.0f-t2*cw1);
			      q[12] = cw2*cw2;
			      q[13] = {1.0f,0.0f};
			      q[14] = 2.0f*t2.real();
			      q[15] = 2.0f*t2.imag();
			      invm4x4_cmplxr4(&q[0],&invq[0]);
		       }
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

                 
	           float gms::math::cgeco(std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
	                       const int32_t lda,
		               const int32_t n,
		               int32_t * __restrict __ATTR_ALIGN__(64) ipvt) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                        a    = (std::complex<float>*)__builtin_assume_aligned(a,64);
			ipvt = (int32_t*)__builtin_assume_aligned(ipvt,64);
#elif defined __ICC || defineed __INTEL_COMPILER
                        __assume_aligned(a,64);
			__assume_aligned(ipvt,64);
#endif
                        std::complex<float> ek,t,wk,wkm;
			std::complex<float> * z;
			float anorm,s,sm,ynorm;
			int32_t i,j,k,l;
		        z = gms::common::gms_cmplxr4_emalloca(static_cast<size_t>(n),64);
			//
                        //  Compute the 1-norm of A.
                       //
		       anorm = 0.0f;
		       for(j = 0; j != n; ++j) {
		           const float t0 = scasum(n,a+0+j*lda,1);
                           anorm = r4_max(anorm,t0);
		       }
		       // Factor
		       cgefa(a,lda,n,ipvt);
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and hermitian(A)*Y = E.
//
//  Hermitian(A) is the conjugate transpose of A.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where hermitian(U)*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(U)*W = E.
//                    
                      ek = {1.0f,0.0f};
		      z = {};
		      for(k = 1; k <= n; ++k) {
                          if(cabs1(z[k-1]) != 0.0f) {
                             ek = csign1(ek,-z[k-1]);
			  }
			  const float t0 = cabs1(a[k-1+(k-1)*lda]);
			  const float t1 = cabs1(ek-z[k-1]);
			  if(t0 < t1) {
                             s = t0/t1;
			     csscal(n,s,z,1);
			     ek = std::complex<float>(s,0.0f)*ek;
			  }
			  wk = ek - z[k-1];
			  wkm = -ek - z[k-1];
			  s = cabs1(wk);
			  sm = cabs1(wkm);
			  const std::complex<float> t2 = std::conj(t0);
			  if(t0 != 0.0f) {
			     wk = wk / t2;
			     wkm = wkm / t2;
			  }
			  else {
                              wk = {1.0f,0.0f};
			      wkm = {1.0f,0.0f};
			  }
			  for(j = k+1; j <= n; ++j) {
			      const std::complex<float> t3 = a[k-1+(j-1)*lda];
			      const std::complex<float> t4 = std::conj(t3);
                              sm = sm + cabs1(z[j-1]+wkm*t4);
			      const std::complex<float> t5 = z[j-1]+wk*t4;
			      z[j-1] = t5;
			      s = s + cabs1(z[j-1]);
			  }
			  if(s < sm) {
                             t = wkm-wk;
			     wk = wkm;
			     for(j = k+1; j <= n; ++j) {
                                 std::complex<float> t6 = z[j-1]+t*std::conj(a[k-1+(j-1)*lda]);
				 z[j-1] = t6;
			     }
			  }
			  z[k-1] = wk;
		      }
		      s = 1.0f/scasum(n,z,1);
		      csscal(n,s,z,1);
		      for(k = n; 1 <= k; --k) {
                          if(k < n) {
                             const std::complex<float> t7 = z[k-1] + cdotc(n-k,a+k+(k-1)*lda,1,z+k,1);
			     z[k-1] = t7;
			  }
			  float t8 = cabs1(z[k-1]);
			  if(1.0f < t8) {
                             s = 1.0f/t8;
			     csscal(n,s,z,1);
			  }
			  l = ipvt[k-1];
			  t = z[l-k];
			  z[l-1] = z[k-1];
			  z[k-1] = t;
		      }

		      s = 1.0f / scasum(n,z,1);
		      csscal(n,s,z,1);
		      ynorm = 1.0f;
		      for(k = 1; k != n; ++k) {
                          l = ipvt[k-1];
			  t = z[l-1];
			  z[l-1] = z[k-1];
			  z[k-1] = t;
			  if(k < n) {
                             caxpy(n-k,t,a+k+(k-1)*lda,1,z+k,1);
			  }
			  float t9 = cabs1(z[k-1]);
			  if(1.0f < t9) {
                             s = .0f/t9;
			     csscal(n,s,z,1);
			     ynorm *= s;
			  }
		      }

		      s = 1.0f/scasum(n,z,1);
		      csscal(n,s,z,1);
		      ynorm *= s;

		      for(k = n; 1 <= k, k--) {
                          const float t10 = cabs1(a[k-1+(k-1)*lda;
			  const float t11 = z[k-1];
			  if(t10 < t11) {
                             s = t10/t11;
			     csscal(n,s,z,1);
			     ynorm *= s;
			  }
			  if(t10 != 0.0f) {
                             const std::complex<float> c0 = z[k-1]/a[k-1+(k-1)*lda];
			     z[k-1] = c0;
			  }
			  else {
                             z[k-1] = {1.0f,0.0f};
			  }
			  t = -z[k-1];
			  caxpy(k-1,t,a+0+(k-1)*lda,1,z,1);
		      }

		      s = 1.0f/scasum(n,z,1);
		      csscal(n,s,z,1);
		      ynorm *= s;
		      if(anorm != 0.0f) {
                         rcond = ynorm/anorm;
		      }
		      else {
                          rcond = 0.0f;
		      }
		      _mm_free(z);
		      return (rcond);
	         }

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

                  
		   void gms::math::cgedi(std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		              const int32_t lda,
			      const int32_t n,
			      int32_t * __restrict __ATTR_ALIGN__(64) ipvt,
			      std::complex<float> det[2],
			      const int32_t job) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
			    a = (std::complex<float>*)__builtin_assume_aligned(a,64);
			    ipvt = (int32_t*)__builtin_assume_aligned(ipvt,64);
#elif defined __ICC || defined __INTEL_COMPILER
			    __assume_aligned(a,64);
			    __assume_aligned(ipvt,64);
#endif
			    int i;
                            int j;
                            int k;
                            int l;
			    std::complex<float>t;
			    std::complex<float> *work = NULL;
//
//  Compute the determinant.
//
                            if(job/10 != 0){
                               det[0] = std::complex<float>(1.0f,0.0f );
                               det[1] = std::complex<float>(0.0f,0.f0 );
                               for(i = 1; i <= n; i++) {
				   if(ipvt[i-1] != i){
                                      det[0] = -det[0];
                                  }
                                  det[0] = a[i-1+(i-1)*lda] * det[0];
                                  if(cabs1(det[0]) == 0.0f ){
                                     break;
                                  }
                                  while(cabs1(det[0]) < 1.0f ){
				        det[0] = det[0] * std::complex<float>(10.0f,0.0f);
                                        det[1] = det[1] - std::complex<float>(1.0f,0.0f);
                                   }
                                  while(10.0f <= cabs1(det[0])){
                                       det[0] = det[0] /std::complex<float>(10.0f,0.0f);
                                       det[1] = det[1] + std::complex<float>(1.0f,0.0f);
                                   }
                               }
                           }
//
//  Compute inverse(U).
//
                           if((job % 10) != 0) {
                               work = gms::common::gms_cmplxr4_emalloca(static_cast<size_t>(n),64);
                               for(k = 1; k <= n; k++) {
			           std::complex<float> t0 = std::complex<float>(1.0,0.0)/a[k-1+(k-1)*lda];
                                   a[k-1+(k-1)*lda] = t0;
                                   t = -a[k-1+(k-1)*lda];
                                   cscal ( k-1, t, a+0+(k-1)*lda, 1 );
                                   for(j = k+1; j <= n; j++) {
                                       t = a[k-1+(j-1)*lda];
                                       a[k-1+(j-1)*lda] = std::complex<float>(0.0f,0.0f );
                                       caxpy(k,t,a+0+(k-1)*lda,1,a+0+(j-1)*lda,1);
                                   }
                               }
//
//  Form inverse(U) * inverse(L).
//
                              for(k = n-1; 1 <= k; k--) {
                                  for(i = k+1; i <= n; i++) {
                                      work[i-1] = a[i-1+(k-1)*lda];
                                      a[i-1+(k-1)*lda] =std::complex<float>(0.0f,0.0f);
                                    }
                                   for(j = k+1; j <= n; j++) {
                                       t = work[j-1];
                                       caxpy(n, t, a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1);
                                   }
                                   l = ipvt[k-1];
                                   if(l != k) {
                                      cswap(n,a+0+(k-1)*lda, 1,a+0+(l-1)*lda, 1);
                                   }
                                }
                                _mm_free(work);
                            }
		    	}

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

                
                  int32_t gms::math::cgefa(std::complex<float> * __restrict __ATTR_ALIGN__(64) a,
		                const int32_t lda,
				const int32_t n,
				int32_t * __restrict __ATTR_ALIGN__(64) ipvt) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
			    a = (std::complex<float>*)__builtin_assume_aligned(a,64);
			    ipvt = (int32_t*)__builtin_assume_aligned(ipvt,64);
#elif defined __ICC || defined __INTEL_COMPILER
			    __assume_aligned(a,64);
			    __assume_aligned(ipvt,64);
#endif
			     int info;
                             int j;
                             int k;
                             int l;
			    std::complex<float>t;
//
//  Gaussian elimination with partial pivoting.
//
                            info = 0;
                            for(k = 1; k <= n-1; k++) {
                  
//
//  Find L = pivot index.
//
                                 l = icamax ( n-k+1, a+(k-1)+(k-1)*lda, 1 ) + k - 1;
                                 ipvt[k-1] = l;
//
//  Zero pivot implies this column already triangularized.
//
                                 if(cabs1(a[l-1+(k-1)*lda]) == 0.0f){
                                      info = k;
                                      continue;
                                 }
//
//  Interchange if necessary.
//
                                 if (l != k) {
                                    t   = a[l-1+(k-1)*lda];
                                    a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
                                    a[k-1+(k-1)*lda] = t;
                                 }
//
//  Compute multipliers
//
                                 t = -std::complex<float>(1.0f, 0.0f)/a[k-1+(k-1)*lda];
                                 cscal(n-k,t,a+k+(k-1)*lda, 1 );
//
//  Row elimination with column indexing
//
                                 for(j = k+1; j <= n; j++) {
                                       t = a[l-1+(j-1)*lda];
                                       if(l != k) {
                                             a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
                                             a[k-1+(j-1)*lda] = t;
                                       }
                                       caxpy ( n-k, t, a+k+(k-1)*lda, 1, a+k+(j-1)*lda, 1 );
                                 }

                              }
                              ipvt[n-1] = n;
                              if cabs1(a[n-1+(n-1)*lda]) == 0.0f){
                                   info = n;
                              }
                              return info;
		     }

		     
		  


