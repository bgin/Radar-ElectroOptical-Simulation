
#ifndef __GMS_AVXVECF32_H__
#define __GMS_AVXVECF32_H__


namespace file_info{


       const unsigned int gGMS_AVXVECF32_MAJOR = 1;
       const unsigned int gGMS_AVXVECF32_MINOR = 0;
       const unsigned int gGMS_AVXVECF32_MICRO = 0;
       const unsigned int gGMS_AVXVECF32_FULLVER =
          1000U*gGMS_AVXVECF32_MAJOR+100U*gGMS_AVXVECF32_MINOR+10U*gGMS_AVXVECF32_MICRO;
       const char * const pgGMS_AVXVECF32_CREATION_DATE = "22-12-2019 12:09 +00200 (SUN 22 DEC 2019 GMT+2)";
       const char * const pgGMS_AVXVECF32_BUILD_DATE    = __DATE__ " " __TIME__;
       const char * const pgGMS_AVXVECF32_AUTHOR        = "Programmer: Bernard Gingold, contact; beniekg@gmail.com";
       const char * const pgGMS_AVXVECF32_SYNOPSIS      = "Wrapper class around __m256 data type.";
}

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {
       namespace math {

                 struct AVXVec8 {


		     __m256 m_v8 __ATTR_ALIGN__(32) ;

                     AVXVec8() {
                       m_v8 = _mm256_setzero_ps();
		     }
		     
		     
		     AVXVec8(const float* __restrict __ATTR_ALIGN__(32) v){
                        v = (const float*)__builtin_assume_aligned(32,v);
			m_v8 = _mm256_load_ps(&v[0]);
		     }
		     

		     AVXVec8(  const float s0,
		               const float s1,
			       const float s2,
			       const float s3,
			       const float s4,
			       const float s5,
			       const float s6,
			       const float s7) {
                         m_v8 = _mm256_setr_ps(s0,s1,s2,s3,s4,s5,s6,s7);
		     }
		     
		     

		     AVXVec8(const float s) {
                        m_v8 = _mm256_set1_ps(s);
		     }

		     AVXVec8(const __m256 v) {
                        m_v8 = v;
		     }

		     AVXVec8(const __m256i v) {
                        m_v8 = _mm256_castsi256_ps(v);
		     }

		     AVXVec8(const AVXVecF32 &x) {
                        m_v8 = x.m_v8;
		     }

		     AVXVec8(const __m128 lo,
		             const __m128 hi) {
                        m_v8 = _mm256_insertf128_ps(m_v8,lo,0);
			m_v8 = _mm256_insertf128_ps(m_v8,hi,1);
		     }

		    ~AVXVec8() = default;

		    __m128 lo_part() const {
                      return (_mm256_extractf128_ps(m_v8,0));
		    }

		    __m128 hi_part() const {
                      return (_mm256_extractf128_ps(m_v8,1));
		    }

		    //
		    // Load-store functions
		    //
		    // Address aligned on 32-byte boundary
                    AVXVec8   & load_a(const float * __restrict __ATTR_ALIGN__(32) addr) {
                          addr = (const float*)__builtin_assume_aligned(32,addr);
			  m_v8 = _mm256_load_ps(&addr[0]);
			  return (*this);
		    }

		    AVXVec8   & load_u(const float * __restrict  addr) {
                          m_v8 = _mm256_loadu_ps(&addr[0]);
			  return (*this);
		    }
		    
		    // Address argument should be aligned on 32-byte boundary
		    void store_a(float * __restrict __ATTR_ALIGN__(32) addr) const {
                          addr = (float*)__builtin_assume_aligned(32,addr);
			  _mm256_store_ps(&addr[0],m_v8);
		    }

		    void store_u(float * __restrict addr ) const {
                          _mm256_store_ps(&addr[0],m_v8);
		    }

		    void stream_store(float * __restrict addr ) const {
                         _mm256_stream_ps(&addr[0],m_v8);
			 _mm_sfence();
		    }

		    float extract_scalar(const uint32_t idx) const {
                         float t[8] __ATTR_ALIGN__(32) = {};
			 store_a(t);
			 return (t[idx&7]);
		    }
		    
                    // Inserts a single float in location 0...7
		    AVXVec8 & insert_scalar(const uint32_t pos
		                            const float val) {
                         float mem[8] __ATTR_ALIGN__(32) = {};
			 store_a(&mem[0]);
			 mem[pos & 7] = val;
			 load_a(&mem[0]);
			 return (*this);
		    }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		           AVXVec8 & operator=(const AVXVec8 &x)   {
                               if(this == &x) return (*this);
			       m_v8 = x.m_v8
			       return (*this);
		    }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		           AVXVec8 & operator=(const AVXVec8 x) {
                               if(this == &x) return (*this);
			       m_v8 = x.m_v8
			       return (*this);
		    }

		    operator __m256 () const { return (m_v8);}

		    float operator[](const uint32_t pos) const {
                          return (reinterpret_cast<const float*>(&m_v8)[pos];)
		    }

	   } __attribute__ ((aligned(64)));

	   	//
		// global(namespace) static functions
		//

		//
		// Extract __m128 part. Value of second parameter
		// must be 0 or 1 only.
		__ATTR_HOT__  __ATTR_ALIGN__(16)
                static inline __m128  extract_vec4(AVXVec8 v,
		                                  const uint32_t pos){
			 return (_mm256_extractf128_ps(v,pos));		  
		      }
	 	                             
                __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 select_vec8(const AVXVec8 a,
		                                  const AVXVec8 b,
						  const __m256 pred) {
                         return (_mm256_blendv_ps(a,b,pred));
		      }

		//
		//	Arithmetic and mathematical operations
		//
	        __ATTR_HOT__  __ATTR_ALIGN__(16)
                static inline AVXVec8 max(const AVXVec8 x,
		                          const AVXVec8 y) {
                         return (_mm256_max_ps(x,y));
		      }
		      
                __ATTR_HOT__  __ATTR_ALIGN__(16) 
		static inline AVXVec8 min(const AVXVec8 x,
		                          const AVXVec8 y) {
                         return (_mm256_min_ps(x,y));
		      }
		      
                __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 abs(const AVXVec8 x)  {
                        sttaic  const __m256 mask = _mm256_set1_ps(0x7FFFFFFF);
			 return (_mm256_and_ps(x,mask));
		      }

		__ATTR_HOT__  __ATTR_ALIGN__(16)      
		static inline AVXVec8 sqrt(const AVXVec8 x) {
                         return (_mm256_sqrt_ps(x));
	  	     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)     
		static inline AVXVec8 ceil(const AVXVec8 x)  {
                         return (_mm256_ceil_ps(x));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 floor(const AVXVec8 x)   {
                         return (_mm256_floor_ps(x));
	  	     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)     
		static inline AVXVec8 round(const AVXVec8 x,
		                            int round)  {
                         return (_mm256_round_ps(x,round));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 addsub(const AVXVec8 x,
		                             const AVXVec8 y) {
                         return (_mm256_addsub_ps(x,y));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)     
                static inline AVXVec8 dot(const AVXVec8 x,
		                          const AVXVec8 y,
					  const int32_t imm8) {
                         return (_mm256_dp_ps(x,y,imm8));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)     
		static inline AVXVec8 hadd(const AVXVec8 x,
		                           const AVXVec8 y) {
                         return (_mm256_hadd_ps(x,y));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 hsub(const AVXVec8 x,
		                           const AVXVec8 y) {
                         return (_mm256_hsub_ps(x,y));
		     }

	            
                __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 rsqrt(const AVXVec8 x) {
                         return (__m256_rsqrt_ps(x));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
	        static inline AVXVec8 rcp(const AVXVec8 x) {
                         return (_mm256_rcp_ps(x));
		     }

		__ATTR_HOT__  __ATTR_ALIGN__(16)     
		static inline AVXVec8 i32gather(float const* __restrict base_addr,
		                                __m256i vindex , const int scale)   {
                         return (_mm256_i32gather_ps(&base_addr[0],vindex,scale));
		    }				

		__ATTR_HOT__  __ATTR_ALIGN__(16)    
		static inline AVXVec8 mask_i32gather(AVXVec8 src,
		                                     float const* __restrict base_addr,
						     __m256i vindex , __m256 mask,
						     const int scale)     {
                          return (_mm256_mask_i32gather_ps(src,&base_addr[0],vindex,mask,scale));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 mask_load(float const* __restrict addr,
		                                __m256i mask) {
                          return (_mm256_maskload_ps(addr,mask));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
                static inline AVXVec8 loadu2_m128(float const* __restrict hiaddr,
		                                  float const* __restrict loaddr) {
                          return (_mm256_loadu2_m128(&hiaddr[0],&loaddr[0]));
		    }
		// Return true if all sign  bits are 1

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
	        static inline int    testc(const AVXVec8 x,
		                           const AVXVec8 y) {
                          return (_mm256_testc_ps(x,y));
		    }

		// Return true if all sign bit is zero
		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline bool    testz(const AVXVec8 x,
		                            const AVXVec8 y) {
                          return (_mm256_testz_ps(x,y));
		    }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 permute(const AVXVecx x,
		                              int imm8)   {
                          return (_mm256_permute_ps(x,imm8));
		    }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 unary_minus(const AVXVec8 x){
                          return (_mm256_sub_ps(_mm256_setzero_ps(),x));
	 	    }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 sin(const AVXVec8 x) {
                          return (_mm256_sin_ps(x));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 cos(const AVXVec8 x) {
                          return (_mm256_cos_ps(x));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 sinh(const AVXVec8 x) {
                          return (_mm256_sinh_ps(x));
		    }

		__ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 cosh(const AVXVec8 x) {
                          return (_mm256_cosh_ps(x));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 tan(const AVXVec8 x) {
                          return (_mm256_tan_ps(x));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)    
		static inline AVXVec8 tanh(const AVXVec8 x) {
                          return (_mm256_tanh_ps(x));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)   
		static inline AVXVec8 asin(const AVXVec8 x) {
                          return (_mm256_asin_ps(x));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)   
		static inline AVXVec8 asinh(const AVXVec8 x) {
                          return (_mm256_asinh_ps(x));
		    }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)    
		static inline AVXVec8 acos(const AVXVec8 x) {
                          return (_mm256_acos_ps(x));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 acosh(const AVXVec8 x) {
                          return (_mm256_acosh_ps(x));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 atan(const AVXVec8 x) {
                          return (_mm256_atan_ps(x));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 atanh(const AVXVec8 x) {
                          return (_mm256_atanh_ps(x));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 atan2(const AVXVec8 x,
		                            const AVXVec8 y) {
                          return (_mm256_atan2_ps(x,y));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
	        static inline AVXVec8 exp(const AVXVec8 x) {
                          return (_mm256_exp_ps(x));
		    }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 log10(const AVXVec8 x) {
                          return (_mm256_log10_ps(x));
		   }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 log(const AVXVec8 x)  {
                          return (_mm256_log_ps(x));
		   }

	           __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 pow(const AVXVec8 x,
		                          const AVXVec8 y) {
                          return (_mm256_pow_ps(x,y));
		   }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)  
	        static inline AVXVec8 fmadd(const AVXVec8 a,
		                            const AVXVec8 b,
					    const AVXVec8 c)  {
                          return (_mm256_fmadd_ps(a,b,c));
		   }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
	        static inline AVXVec8 fmaddsub(const AVXVec8 a,
		                               const AVXVec8 b,
					       const AVXVec8 c) {
                          return (_mm256_fmaddsub_ps(a,b,c));
		   }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 fmsub(const AVXVec8 a,
		                            const AVXVec8 b,
					    const AVXVec8 c)  {
                           return (_mm256_fmsub_ps(a,b,c));
		   }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 fmsubadd(const AVXVec8 a,
		                               const AVXVec8 b,
					       const AVXVec8 c) {
                           return (_mm256_fmsubadd_ps(a,b,c));
		   }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 fnmadd(const AVXVec8 a,
		                             const AVXVec8 b,
					     const AVXVec8 c) {
                           return (_mm256_fnmadd_ps(a,b,c));
		   }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 fnmsub(const AVXVec8 a,
		                             const AVXVec8 b,
					     const AVXVec8 c) {
                           return (_mm256_fnmsub_ps(a,b,c));
		   }

	          __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline int v8_lt_zero(const AVXVec8 x) {
                           __m256 vzmask = _mm256_setzero_ps();
			   static const __m256 vzero = _mm256_set1_ps(0.f);
			   vzmask = _mm256_cmp_ps(vzero,x,_CMP_LT_OQ);
			   return (_mm256_testc_ps(vzmask,vzero));
		   }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline bool  v8_bt_ones(const AVXVec8 x) {
                          __m256 vmaskgt = _mm256_setzero_ps();
			  __m256 vmasklt = _mm256_setzero_ps();
			  static const __m256 one    = _mm256_set1_ps(1.f);
			  static const __m256 negone = _mm256_set1_ps(-1.0f);
			  static const __m256 zero   = _mm256_setzero_ps();
			  vmaskgt = _mm256_cmp_ps(one,x,_CMP_GT_OQ);
			  vmasklt = _mm256_cmp_ps(negone,x,_CMP_LT_OQ);
			  if(_mm256_testc_ps(vmaskgt,zero) &&
			     _mm256_testc_ps(vmasklt,zero))
			     return (true);
			  else
			     return (false);
		   }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		  static inline int v8_ne_zero(const AVXVec8 x) {
                         __m256 vmaskz = _mm256_setzero_ps();
			 static const __m256 zero = _mm256_setzero_ps();
			 vmaskz = _mm256_cmp_ps(zero,x,_CMP_EQ_OQ);
			 return (_mm256_testc_ps(vmaskz,zero));
		  }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline int v8_ne_negone(const AVXVec8 x) {
                         _mm256 vmaskz = _mm256_setzero_ps();
			 static const __m256 negone = _mm256_set1_ps(-1.f);
			 vmaskz = _mm256_cmp_ps(negone,x,_CMP_LT_OQ);
			 return (_mm256_testc_ps(vmaskz,zero));
		 }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline bool v8_denormal_present(const AVXVec8 x) {
                         using namespace std;
			 __m256 xmaskz = _mm256_setzero_ps();
			 __m256 ymaskz = _mm256_setzero_ps();
			 static const __m256 zero = _mm256_setzero_ps();
			 static const __m256 tiny = _mm256_set1_ps(numeric_limits<float>::min());
			 constexpr uint32_t absmsk = 0x7FFFFFFFF;
			 ymask = _mm256_cmp_ps(zero,x,_CMP_NEQ_OQ);
			 __m256 vabs = _mm256_and_ps(x,_mm256_set1_ps(absmsk));
			 xmask = _mm256_cmp_ps(vabs,tiny,_CMP_LT_OQ);
			 if(_mm256_testc_ps(ymask,zero) &&
			    _mm256_testc_ps(xmask,zero))
			    return (true);
			 else
			    return (false);
		 }

		// C = A+B, vector + vector
		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator+(const AVXVec8 x, 
						const AVXVec8 y) {
                       return (_mm256_add_ps(x,y));
		 }

		// C = A+B, vector + scalar
		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator+(const AVXVec8 x, 
						  const float s) {
                       return (_mm256_add_ps(x,_mm256_set1_ps(s)));
		 }

		  __ATTR_HOT__  __ATTR_ALIGN__(16)
                static inline AVXVec8 operator+(const AVXVec8 x,
		                                const __m256 y) {
                       return (_mm256_add_ps(x,y));
		 }
		 
		// C = A+B, scalar + vector
		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator+(const float s, 
						  const AVXVec8 x)  {
                       return (_mm256_add_ps(_mm256_set1_ps(s),x));
		 }

		 __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator+(const __m256 x,
		                                const AVXVec8 y) {
                       return (_mm256_add_ps(x,y));
		 }

		// A = A+B, vector + vector (in-place)
		  __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator+=(AVXVec8 x,
						 const AVXVec8 y) {
                       x = x+y
		       return (x);
		 }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline AVXVec8 operator+=(AVXVec8 x,
		                                  const float s) {
                       x = x+AVXVec8{s};
		       return (x);
		 }

		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline AVXVec8 operator+=(AVXVec8 x,
		                                  const __m256 y) {
                       x = x+y;
		       return (x);
		 }

		// A = A+1
		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator++(AVXVec8 x) {
                       x = x+AVXVec8{-1.f};
		       return (x);
		 }

	    // C = A-B, vector - vector
	           __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator-(const AVXVec8 x,
						const AVXVec8 y) {
                        return (_mm256_sub_ps(x,y));
		 }

		// C = A-B, vector - scalar
		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator-(const AVXVec8 x, 
						  const float s) {
                        return (_mm256_sub_ps(x,_mm256_set1_ps(s))); 
		 }

		// C = A-B, scalar - vector
		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator-(const float s, 
						const AVXVec8 x)   {
                        return (_mm256_sub_ps(_mm256_set1_ps(s),x));
		 }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline AVXVec8 operator-(const AVXVec8 x,
		                                 const __m256 y) {
                        return (_mm256_sub_ps(x,y));
		 }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline AVXVec8 operator-(const __m256 x,
		                                 const AVXVec8 y) {
                        return (_mm256_sub_ps(x,y));
		 }

		// A = A-B, vector - vector (in-place)
		   __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator-=(AVXVec8 x, 
						 const AVXVec8 y) {
                      x = x-y;
		      return (x);
		 }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		static inline AVXVec8 operator-=(AVXVec8 x,
		                                 const float s) {
                      x = x-AVXVec8{s};
		      return (x);
		 }

		    __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline AVXVec8 operator-=(const float s,
		                                  AVXVec8 x) {
                      x = AVXVec8{s}-x;
		      return (x);
		 }

		     __ATTR_HOT__  __ATTR_ALIGN__(16)
		 static inline AVXVec8 operator-=(AVXVec8 x,
		                                  const __m256 y) {
                      x = x-y;
		      return (x);
		 }

		

			// A = A-1.0L
		   __ATTR_HOT__  __ATTR_ALIGN__(16)	
		static inline AVXVec8 operator--(AVXVec8 x)  {
                     x = x-AVXVec8{-1.f};
		     return (x);
		}

		// C = A*B, vector * vector
		static inline AVXVec8 operator*(const AVXVec8, 
						  const AVXVec8)       

	    // C = A*B, vector * scalar
		static inline AVXVec8 operator*(const AVXVec8,
						  const float)         

		// C = A*B, scalar * vector
		static inline AVXVec8 operator*(const float,
						const AVXVec8)        

		// A = A*B, vector * vector (in-place)
		static inline AVXVec8 operator*=(AVXVec8, 
						 const AVXVec8)
						 
		// C = A/B, vector / vector
		static inline AVXVec8 operator/(const AVXVec8, 
						  const AVXVec8)       

		// C = A/B, vector / scalar
		static inline AVXVec8 operator/(const AVXVec8
						  const float)          

		// C = A/B, scalar / vector
		static inline AVXVec8 operator/(const float, 
						  const AVXVec8)        

		// A = A/B, vector / vector (in-place)
		static inline AVXVec8 operator/=(AVXVec8, 
						 const AVXVec8)         

	        //

		// C = A==B, vector == vector, C is of type __m256d
		static inline __m256 operator==(const AVXVec8, 
						 const AVXVec8)
						 
		static inline __m256 operator==(const AVXVec8
						 const float )          

		static inline __m256 operator==(const float,
						 const AVXVec8)         

		// C = A != b, vector != vector, C is of type __m256d
		static inline __m256 operator!=(const AVXVec8, 
						 const AVXVec8)          

		// C = A>B, vector > vector, C is of type __m256d
		static inline __m256 operator>(const AVXVec8, 
						const AVXVec8)
						
		// C = A<B, vector < vector, C is of type __m256d
		static inline __m256 operator<(const AVXVec8, 
						const AVXVec8)           

		// C = A>=B, vector >= B, C is of type __m256d
		static inline __m256 operator>=(const AVXVec8, 
						 const AVXVec8)           

		// C = A<=B, vector <= vector, C is of type __m256d
		static inline __m256 operator<=(const AVXVec8, 
						 const AVXVec8)          

			// C = A&B, vector & vector
		static inline AVXVec8 operator&(const AVXVec8,
						  const AVXVec8)           

		
		// A = A&B, vector & vector (in-place)
		static inline AVXVec8 operator&=(AVXVec8,
						   const AVXVec8)         


		// C = A | B, vector | vector
		static inline AVXVec8 operator|(const AVXVec8,
						  const AVXVec8)          

		
		// A = A | B, vector | vector
		static inline AVXVec8 operator|=(AVXVec8,
						   const AVXVec8)        

		// C = A ^ B, vector ^ vector
		static inline AVXVec8 operator^(const AVXVec8,
						  const AVXVec8)        

		// A = A ^ B, vector ^ vector
		static inline AVXVec8 operator^=(AVXVec8,
						   const AVXVec8)        			 





    }
}





#endif  /*__GMS_AVXVECF32_H__*/
