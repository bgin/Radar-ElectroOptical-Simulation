
#include <omp.h>
#include "GMS_rotation_kernels_avx512_omp.h"
#include "GMS_rotation_kernels_avx512.hpp"


      void
      gms::math::
      q4x16_to_rmat9x16_zmm16r4_omp(const __m512 * __restrict q_x,
		                    const __m512 * __restrict q_y,
				    const __m512 * __restrict q_z,
				    const __m512 * __restrict q_w,
				    DCM9x16 * __restrict mat,
				    const int32_t n) {
            
            int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	    
#pragma omp parallel for schedule(static,4) private(i) default(none) \
                         shared(n,q_x,q_y,q_z,q_w,mat)
	    for(i = 0; i != n; ++i) {
                mat[i] = q4x16_to_rmat9x16_zmm16r4(q_x[i],
		                                   q_y[i],
					           q_z[i],
					           q_w[i]);
					
	    }

      }


      void
      gms::math::
      q4x8_to_rmat9x8_zmm8r8_omp(const __m512d * __restrict q_x,
		                 const __m512d * __restrict q_y,
				 const __m512d * __restrict q_z,
				 const __m512d * __restrict q_w,
				 DCM9x8 * __restrict mat,
				 const int32_t n) {

            int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	    
#pragma omp parallel for schedule(static,4) private(i) default(none) \
                         shared(n,q_x,q_y,q_z,q_w,mat)
	    for(i = 0; i != n; ++i) {
                mat[i] = q4x8_to_rmat9x8_zmm8r8(q_x[i],
		                                q_y[i],
				                q_z[i],
				                q_w[i]);
				      
	    }	    
      }


      void
      gms::math::
      rand_sphere_rm9x16_zmm16r4_omp(const __m512 * __restrict vr1,
		                     const __m512 * __restrict vr2,
				     const __m512 * __restrict vr3,
				     DCM9x16 * __restrict mat,
				     const int32_t n) {

            int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	    
#pragma omp parallel for schedule(static,4) private(i) default(none) \
                         shared(n,vr1,vr2,vr3,mat)
	    for(i = 0; i != n; ++i) {
                mat[i] = random_sphere_rm9x16_zmm16r4(vr1[i],
		                                      vr2[i],
						      vr3[i]);
	    }
      }


      void
      gms::math::
      rand_sphere_rm9x8_zmm8r8_omp(const __m512d * __restrict vr1,
		                   const __m512d * __restrict vr2,
				   const __m512d * __restrict vr3,
				   DCM9x8 * __restrict mat,
				   const int32_t n) {
            int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	    
#pragma omp parallel for schedule(static,4) private(i) default(none) \
                         shared(n,vr1,vr2,vr3,mat)
	    for(i = 0; i != n; ++i) {
                mat[i] = random_sphere_rm9x8_zmm8r8(vr1[i],
		                                    vr2[i],
						    vr3[i]);
	    }

     }


      void
      gms::math::
      urand_q4x16_a_zmm16r4_omp(const __m512 * __restrict vrx,
		                const __m512 * __restrict vry,
				const __m512 * __restrict vrz,
				const __m512 * __restrict vrw,
				float * __restrict __ATTR_ALIGN__(64) q_x,
				float * __restrict __ATTR_ALIGN__(64) q_y,
				float * __restrict __ATTR_ALIGN__(64) q_z,
				float * __restrict __ATTR_ALIGN__(64) q_w,
				const int32_t n) {
             const int32_t stride = 16;
             int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(stride,n,vrx,vry,vrz,vrw,q_x,q_y,q_z,q_w)
	    for(i = 0; i != n; ++i) {
                  urand_q4x16_a_zmm16r4(vrx[i],
		                        vry[i],
					vrz[i],
					vrw[i],
					q_x[j],
					q_y[j],
					q_z[j],
					q_w[j]);
	          j += stride;
	    }
      }


       void
       gms::math::
       urand_q4x16_u_zmm16r4_omp(const __m512 * __restrict vrx,
		                 const __m512 * __restrict vry,
				 const __m512 * __restrict vrz,
				 const __m512 * __restrict vrw,
				 float * __restrict q_x, //  length of float arrays must be multiplicity of 16
				 float * __restrict q_y, //  length of float arrays must be multiplicity of 16
				 float * __restrict q_z, //  length of float arrays must be multiplicity of 16
				 float * __restrict q_w, //  length of float arrays must be multiplicity of 16
				 const int32_t n) {
	     if(__builtin_expect((n%16)!=0,1)) {return;}
             const int32_t stride = 16;
             int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(stride,n,vrx,vry,vrz,vrw,q_x,q_y,q_z,q_w)
	    for(i = 0; i != n; ++i) {
                  urand_q4x16_u_zmm16r4(vrx[i],
		                        vry[i],
					vrz[i],
					vrw[i],
					q_x[j],
					q_y[j],
					q_z[j],
					q_w[j]);
	          j += stride;
	    }
      }


       void
       gms::math::
       urand_q4x8_a_zmm8r8_omp( const __m512d * __restrict vrx,
		                const __m512d * __restrict vry,
				const __m512d * __restrict vrz,
				const __m512d * __restrict vrw,
				double * __restrict  __ATTR_ALIGN__(64) q_x,
				double * __restrict  __ATTR_ALIGN__(64) q_y,
				double * __restrict  __ATTR_ALIGN__(64) q_z,
				double * __restrict  __ATTR_ALIGN__(64) q_w,
				const int32_t n) {

             if(__builtin_expect((n%8)!=0,1)) {return;}
	     const int32_t stride = 8;
	     int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(stride,n,vrx,vry,vrz,vrw,q_x,q_y,q_z,q_w)
	    for(i = 0; i != n; ++i) {
                  urand_q4x8_a_zmm8r8(  vrx[i],
		                        vry[i],
					vrz[i],
					vrw[i],
					q_x[j],
					q_y[j],
					q_z[j],
					q_w[j]);
	          j += stride;
	    }
       }


       void
       gms::math::
       urand_q4x8_u_zmm8r8_omp( const __m512d * __restrict vrx,
		                const __m512d * __restrict vry,
				const __m512d * __restrict vrz,
				const __m512d * __restrict vrw,
				double * __restrict  q_x,
				double * __restrict  q_y,
				double * __restrict  q_z,
				double * __restrict  q_w,
				const int32_t n) {

             if(__builtin_expect((n%8)!=0,1)) {return;}
	     const int32_t stride = 8;
	     int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(stride,n,vrx,vry,vrz,vrw,q_x,q_y,q_z,q_w)
	    for(i = 0; i != n; ++i) {
                  urand_q4x8_u_zmm8r8(  vrx[i],
		                        vry[i],
					vrz[i],
					vrw[i],
					q_x[j],
					q_y[j],
					q_z[j],
					q_w[j]);
	          j += stride;
	    }
       }


       void
       gms::math::
       q4x16_to_ea3x16_a_zmm16r4_omp(const __m512 * __restrict q_x,
		                     const __m512 * __restrict q_y,
				     const __m512 * __restrict q_w,
				     const __m512 * __restrict q_z,
				     float * __restrict  __ATTR_ALIGN__(64) alpha,
				     float * __restrict  __ATTR_ALIGN__(64) beta,
				     float * __restrict  __ATTR_ALIGN__(64) gamma,
				     const int32_t n) {

             if(__builtin_expect((n%16)!=0,1)) {return;}
	     const int32_t stride = 16;
	     int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 q4x16_to_ea3x16_a_zmm16r4( q_x[i],
		                            q_y[i],
					    q_w[i],
					    q_z[i],
					    alpha[j],
					    beta[j],
					    gamma[j]);
		  j += stride;
	    }
       }


       void
       gms::math::
       q4x16_to_ea3x16_u_zmm16r4_omp(const __m512 * __restrict q_x,
		                     const __m512 * __restrict q_y,
				     const __m512 * __restrict q_w,
				     const __m512 * __restrict q_z,
				     float * __restrict  alpha,
				     float * __restrict  beta,
				     float * __restrict  gamma,
				     const int32_t n) {

             if(__builtin_expect((n%16)!=0,1)) {return;}
	     const int32_t stride = 16;
	     int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 q4x16_to_ea3x16_u_zmm16r4( q_x[i],
		                            q_y[i],
					    q_w[i],
					    q_z[i],
					    alpha[j],
					    beta[j],
					    gamma[j]);
		  j += stride;
	    }
       }


       void
       gms::math::
       q4x8_to_ea3x8_a_zmm8r8_omp(const __m512d * __restrict q_x,
		                  const __m512d * __restrict q_y,
				  const __m512d * __restrict q_z,
				  const __m512d * __restrict q_w,
				  double * __restrict __ATTR_ALIGN__(64) alpha,
				  double * __restrict __ATTR_ALIGN__(64) beta,
				  double * __restrict __ATTR_ALIGN__(64) gamma,
				  const int32_t n) {
				  
             if(__builtin_expect((n%8)!=0,1)) {return;}
	     const int32_t stride = 8;
	     int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 q4x8_to_ea3x8_a_zmm8r8(    q_x[i],
		                            q_y[i],
					    q_w[i],
					    q_z[i],
					    alpha[j],
					    beta[j],
					    gamma[j]);
		  j += stride;
	    }

      }


       void
       gms::math::
       q4x8_to_ea3x8_u_zmm8r8_omp(const __m512d * __restrict q_x,
		                  const __m512d * __restrict q_y,
				  const __m512d * __restrict q_z,
				  const __m512d * __restrict q_w,
				  double * __restrict alpha,
				  double * __restrict beta,
				  double * __restrict gamma,
				  const int32_t n) {
				  
             if(__builtin_expect((n%8)!=0,1)) {return;}
	     const int32_t stride = 8;
	     int32_t i,j;
	     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 q4x8_to_ea3x8_u_zmm8r8(    q_x[i],
		                            q_y[i],
					    q_w[i],
					    q_z[i],
					    alpha[j],
					    beta[j],
					    gamma[j]);
		  j += stride;
	    }

      }


       void
       gms::math::
       q4x16_to_ax4x16_a_zmm16r4_omp(const __m512 * __restrict q_x,
		                     const __m512 * __restrict q_y,
				     const __m512 * __restrict q_z,
				     const __m512 * __restrict q_w,
				     float * __restrict  __ATTR_ALIGN__(64) ax1,
				     float * __restrict  __ATTR_ALIGN__(64) ax2,
				     float * __restrict  __ATTR_ALIGN__(64) ax3,
				     float * __restrict  __ATTR_ALIGN__(64) ax4,
				     const int32_t n) {

              if(__builtin_expect((n%16)!=0,1)) {return;}
	      const int32_t stride = 16;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,ax1,ax2,ax3,ax4)
	    for(i = 0; i != n; ++i) {
                 q4x16_to_ax4x16_a_zmm16r4(q_x[i],
		                           q_y[i],
					   q_z[i],
					   q_w[i],
					   ax1[j],
					   ax2[j],
					   ax3[j],
					   ax4[j]);
			 j += stride;
	    }
       }


       void
       gms::math::
       q4x16_to_ax4x16_u_zmm16r4_omp(const __m512 * __restrict q_x,
		                     const __m512 * __restrict q_y,
				     const __m512 * __restrict q_z,
				     const __m512 * __restrict q_w,
				     float * __restrict   ax1,
				     float * __restrict   ax2,
				     float * __restrict   ax3,
				     float * __restrict   ax4,
				     const int32_t n) {

              if(__builtin_expect((n%16)!=0,1)) {return;}
	      const int32_t stride = 16;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,ax1,ax2,ax3,ax4)
	    for(i = 0; i != n; ++i) {
                 q4x16_to_ax4x16_u_zmm16r4(q_x[i],
		                           q_y[i],
					   q_z[i],
					   q_w[i],
					   ax1[j],
					   ax2[j],
					   ax3[j],
					   ax4[j]);
			 j += stride;
	    }
       }


       void
       gms::math::
       q4x8_to_ax4x8_a_zmm8r8_omp(   const __m512d * __restrict q_x,
		                     const __m512d * __restrict q_y,
				     const __m512d * __restrict q_z,
				     const __m512d * __restrict q_w,
				     double * __restrict   __ATTR_ALIGN__(64) ax1,
				     double * __restrict   __ATTR_ALIGN__(64) ax2,
				     double * __restrict   __ATTR_ALIGN__(64) ax3,
				     double * __restrict   __ATTR_ALIGN__(64) ax4,
				     const int32_t n) {

              if(__builtin_expect((n%8)!=0,1)) {return;}
	      const int32_t stride = 8;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,ax1,ax2,ax3,ax4)
	    for(i = 0; i != n; ++i) {
                 q4x8_to_ax4x8_a_zmm8r8(q_x[i],
		                           q_y[i],
					   q_z[i],
					   q_w[i],
					   ax1[j],
					   ax2[j],
					   ax3[j],
					   ax4[j]);
			 j += stride;
	    }
       }


       void
       gms::math::
       q4x8_to_ax4x8_u_zmm8r8_omp(   const __m512d * __restrict q_x,
		                     const __m512d * __restrict q_y,
				     const __m512d * __restrict q_z,
				     const __m512d * __restrict q_w,
				     double * __restrict    ax1,
				     double * __restrict    ax2,
				     double * __restrict    ax3,
				     double * __restrict    ax4,
				     const int32_t n) {

              if(__builtin_expect((n%8)!=0,1)) {return;}
	      const int32_t stride = 8;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,ax1,ax2,ax3,ax4)
	    for(i = 0; i != n; ++i) {
                 q4x8_to_ax4x8_u_zmm8r8(   q_x[i],
		                           q_y[i],
					   q_z[i],
					   q_w[i],
					   ax1[j],
					   ax2[j],
					   ax3[j],
					   ax4[j]);
			 j += stride;
	    }
       }


       void
       gms::math::
       q4x16_to_rv4x16_a_zmm16r4_omp(const __m512 * __restrict q_x,
		                     const __m512 * __restrict q_y,
				     const __m512 * __restrict q_z,
				     const __m512 * __restrict q_w,
				     float * __restrict __ATTR_ALIGN__(64) r_x,
				     float * __restrict __ATTR_ALIGN__(64) r_y,
				     float * __restrict __ATTR_ALIGN__(64) r_z,
				     float * __restrict __ATTR_ALIGN__(64) r_w,
				     const int32_t n) {

              if(__builtin_expect((n%16)!=0,1)) {return;}
              const int32_t stride = 16;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,r_x,r_y,r_z,r_w)
	      for(i = 0; i != n; ++i) {
                   q4x16_to_rv4x16_a_zmm16r4(q_x[i],
		                             q_y[i],
					     q_z[i],
					     q_w[i],
					     r_x[j],
					     r_y[j],
					     r_z[j],
					     r_w[j]);
			  j += stride;
	     }
	}


       void
       gms::math::
       q4x16_to_rv4x16_u_zmm16r4_omp(const __m512 * __restrict q_x,
		                     const __m512 * __restrict q_y,
				     const __m512 * __restrict q_z,
				     const __m512 * __restrict q_w,
				     float * __restrict  r_x,
				     float * __restrict  r_y,
				     float * __restrict  r_z,
				     float * __restrict  r_w,
				     const int32_t n) {

              if(__builtin_expect((n%16)!=0,1)) {return;}
              const int32_t stride = 16;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,r_x,r_y,r_z,r_w)
	      for(i = 0; i != n; ++i) {
                   q4x16_to_rv4x16_u_zmm16r4(q_x[i],
		                             q_y[i],
					     q_z[i],
					     q_w[i],
					     r_x[j],
					     r_y[j],
					     r_z[j],
					     r_w[j]);
			  j += stride;
	     }
	}


       void
       gms::math::
       q4x8_to_rv4x8_a_zmm8r8_omp(   const __m512d * __restrict q_x,
		                     const __m512d * __restrict q_y,
				     const __m512d * __restrict q_z,
				     const __m512d * __restrict q_w,
				     double * __restrict __ATTR_ALIGN__(64) r_x,
				     double * __restrict __ATTR_ALIGN__(64) r_y,
				     double * __restrict __ATTR_ALIGN__(64) r_z,
				     double * __restrict __ATTR_ALIGN__(64) r_w,
				     const int32_t n) {

              if(__builtin_expect((n%8)!=0,1)) {return;}
              const int32_t stride = 8;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,r_x,r_y,r_z,r_w)
	      for(i = 0; i != n; ++i) {
                   q4x8_to_rv4x8_a_zmm8r8(q_x[i],
		                             q_y[i],
					     q_z[i],
					     q_w[i],
					     r_x[j],
					     r_y[j],
					     r_z[j],
					     r_w[j]);
			  j += stride;
	     }
	}


       void
       gms::math::
       q4x8_to_rv4x8_u_zmm8r8_omp(   const __m512d * __restrict q_x,
		                     const __m512d * __restrict q_y,
				     const __m512d * __restrict q_z,
				     const __m512d * __restrict q_w,
				     double * __restrict r_x,
				     double * __restrict r_y,
				     double * __restrict r_z,
				     double * __restrict r_w,
				     const int32_t n) {

              if(__builtin_expect((n%8)!=0,1)) {return;}
              const int32_t stride = 8;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,q_x,q_y,q_z,q_w,r_x,r_y,r_z,r_w)
	      for(i = 0; i != n; ++i) {
                   q4x8_to_rv4x8_u_zmm8r8(q_x[i],
		                             q_y[i],
					     q_z[i],
					     q_w[i],
					     r_x[j],
					     r_y[j],
					     r_z[j],
					     r_w[j]);
			  j += stride;
	     }
	}


        void
	rmat9x16_to_ea3x16_a_zmm16r4_omp(const DCM9x16 * __restrict mat,
		                         float * __restrict __ATTR_ALIGN__(64) alpha,
					 float * __restrict __ATTR_ALIGN__(64) beta,
					 float * __restrict __ATTR_ALIGN__(64) gamma,
					 const int32_t n) {

	      if(__builtin_expect((n%16)!=0,1)) {return;}
              const int32_t stride = 16;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,mat,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 rmat9x16_to_ea3x16_a_zmm16r4(mat[i],
		                              alpha[j],
					      beta[j],
					      gamma[j]);
		     j += stride;
	    }

       }


        void
	rmat9x16_to_ea3x16_u_zmm16r4_omp(const DCM9x16 * __restrict mat,
		                         float * __restrict  alpha,
					 float * __restrict  beta,
					 float * __restrict  gamma,
					 const int32_t n) {

	      if(__builtin_expect((n%16)!=0,1)) {return;}
              const int32_t stride = 16;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,mat,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 rmat9x16_to_ea3x16_u_zmm16r4(mat[i],
		                              alpha[j],
					      beta[j],
					      gamma[j]);
		     j += stride;
	    }

       }


       
        void
	rmat9x8_to_ea3x8_a_zmm8r8_omp(const DCM9x8 * __restrict mat,
		                      double * __restrict __ATTR_ALIGN__(64) alpha,
				      double * __restrict __ATTR_ALIGN__(64) beta,
				      double * __restrict __ATTR_ALIGN__(64) gamma,
					 const int32_t n) {

	      if(__builtin_expect((n%8)!=0,1)) {return;}
              const int32_t stride = 8;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,mat,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 rmat9x8_to_ea3x8_a_zmm8r8(mat[i],
		                              alpha[j],
					      beta[j],
					      gamma[j]);
		     j += stride;
	    }

       }


        void
	rmat9x8_to_ea3x8_u_zmm8r8_omp(const DCM9x8 * __restrict mat,
		                      double * __restrict  alpha,
				      double * __restrict  beta,
				      double * __restrict  gamma,
					 const int32_t n) {

	      if(__builtin_expect((n%8)!=0,1)) {return;}
              const int32_t stride = 8;
	      int32_t i,j;
	      j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif	     
#pragma omp parallel for schedule(static,4) private(i) firstprivate(j)  \
            default(none) shared(n,stride,mat,alpha,beta,gamma)
	    for(i = 0; i != n; ++i) {
                 rmat9x8_to_ea3x8_u_zmm8r8(mat[i],
		                              alpha[j],
					      beta[j],
					      gamma[j]);
		     j += stride;
	    }

       }


       


	


	


	


	

       


       



       

       


       

      
