
#include "Config.fpp"

module  mod_avx512c16f32


 !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         'mod_avx512c16f32
 !          
 !          Purpose:
 !                      This module contains a decomposed to real and imaginary
 !                      parts a complex vector of 16 elements (complex numbers)
 !                      This representation nicely fits into 2 AVX-512 ZMMx
 !                      registers.
 !          History:
 !                        Date: 03-02-2020
 !                        Time: 11:57 GMT+2
 !                        
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                      Bernard Gingold
 !          
 !                 
 !          References:
 !         
 !                          Own implementation (ported from similar C++ version)
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
 !==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.

     use mod_kinds, only : int1,int4,sp
     use mod_vectypes, only : ZMM16r4_t
     use, intrinsic :: ISO_C_BINDING
     implicit none
     !=====================================================59
     !  File and module information:
     !  version,creation and build date, author,description
     !=====================================================59

     ! Major version
     integer(kind=int4), parameter :: MOD_AVX512C16F32_MAJOR = 1
     ! Minor version
     integer(kind=int4), parameter :: MOD_AVX512C16F32_MINOR = 0
     ! Micro version
     integer(kind=int4), parameter :: MOD_AVX512C16F32_MICRO = 0
     ! Full version
     integer(kind=int4), parameter :: MOD_AVX512C16F32_FULLVER = &
          1000*MOD_AVX512C16F32_MAJOR+100*MOD_AVX512C16F32_MINOR+10*MOD_AVX512C16F32_MICRO
     ! Module creation date
     character(*),       parameter :: MOD_AVX512C16F32_CREATION_DATE = "03-02-2020 11:57 +00200 (MON 03 FEB 2020 GMT+2)"
     ! Module build date
     character(*),       parameter :: MOD_AVX512C16F32_BUILD_DATE    = __DATE__ " " __TIME__
     ! Module author info
     character(*),       parameter :: MOD_AVX512C16F32_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com"
     ! Short description
     character(*),       parameter :: MOD_AVX512C16F32_SYNOPSIS      = "Packed complex vector of 16 elements (single precision)"

     type, public :: AVX512c16f32_t
        public
        sequence
        real(kind=sp), dimension(0:15) :: re
        real(kind=sp), dimension(0:15) :: im
     end type AVX512c16f32_t

     interface operator (+)
         module procedure c16_add_c16
         module procedure c16_add_c1
         module procedure c16_add_v16
         module procedure c16_add_s1
         module procedure c1_add_c16
         module procedure v16_add_c16
         module procedure s1_add_c16
     end interface operator (+)

     interface operator (-)
         module procedure c16_sub_c16
         module procedure c16_sub_c1
         module procedure c16_sub_v16
         module procedure c16_sub_s1
         module procedure c1_sub_c16
         module procedure v16_sub_c16
         module procedure s1_sub_c16
     end interface operator (-)

     interface operator (*)
         module procedure c16_mul_c16
         module procedure c16_mul_c1
         module procedure c16_mul_v16
         module procedure c16_mul_s1
         module procedure c1_mul_c16
         module procedure v16_mul_c16
         module procedure s1_mul_c16
     end interface operator (*)

     interface operator (/)
         module procedure c16_div_c16
         module procedure c16_div_c1
         module procedure c16_div_v16
         module procedure c16_div_s1
         module procedure c1_div_c16
         module procedure v16_div_c16
         module procedure s1_div_c16
     end interface operator (/)

#if (USE_INTRINSIC_VECTOR_COMPARE) == 1

      interface operator (==)
         module procedure c16_eq_c16
         module procedure c16_eq_c1
       
         module procedure c1_eq_c16
        
      end interface operator (==)

      interface operator (/=)
         module procedure c16_neq_c16
         module procedure c16_neq_c1
       
         module procedure c2_neq_c16
      
      end interface operator (/=)

      interface operator (>)
         module procedure c16_gt_c16
         module procedure c16_gt_c1
        
         module procedure c1_gt_c16
     
      end interface operator (>)

      interface operator (<)
         module procedure c16_lt_c16
         module procedure c16_lt_c1
        
         module procedure c1_lt_c16
        
      end interface operator (<)

      interface operator (>=)
         module procedure c16_ge_c16
         module procedure c16_ge_c1
       
         module procedure c1_ge_c16
        
      end interface operator (>=)

      interface operator (<=)
         module procedure c16_le_c16
         module procedure c16_le_c1
        
         module procedure c1_le_c16
        
      end interface operator (<=)
#else
     interface operator (==)
         module procedure c16_eq_c16
         module procedure c16_eq_c1
       
         module procedure c1_eq_c16
        
      end interface operator (==)

      interface operator (/=)
         module procedure c16_neq_c16
         module procedure c16_neq_c1
       
         module procedure c1_neq_c16
      
      end interface operator (/=)

      interface operator (>)
         module procedure c16_gt_c16
         module procedure c16_gt_c2
        
         module procedure c2_gt_c16
     
      end interface operator (>)

      interface operator (<)
         module procedure c16_lt_c16
         module procedure c16_lt_c1
        
         module procedure c1_lt_c16
        
      end interface operator (<)

      interface operator (>=)
         module procedure c16_ge_c16
         module procedure c16_ge_c1
       
         module procedure c1_ge_c16
        
      end interface operator (>=)

      interface operator (<=)
         module procedure c16_le_c16
         module procedure c16_le_c1
        
         module procedure c1_le_c16
        
      end interface operator (<=)

#endif

    contains

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER     
      pure function default_init() result(iq) !GCC$ ATTRIBUTES cold :: default_init !GCC$ ATTRIBUTES vectorcall :: default_init !GCC$ ATTRIBUTES inline :: default_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: default_init
      pure function default_init() result(iq)
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: default_init
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: default_init
      
#endif
        use mod_vecconsts, only : v16_n0
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = v16_n0
        iq.im = v16_n0
      end function default_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function array_init(re,im) result(iq) !GCC$ ATTRIBUTES cold :: array_init !GCC$ ATTRIBUTES vectorcall :: array_init !GCC$ ATTRIBUTES inline :: array_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: array_init
      pure function array_init(re,im) result(iq)
         !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: array_init
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: array_init
#endif
        real(kind=sp), dimension(0:15), intent(in) :: re
        real(kind=sp), dimension(0:15), intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = re
        iq.im = im
      end function array_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function complex1_init(c) result(iq) !GCC$ ATTRIBUTES cold :: complex1_init !GCC$ ATTRIBUTES vectorcall :: complex1_init !GCC$ ATTRIBUTES inline :: complex1_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: complex1_init
      pure function complex1_init(c) result(iq)
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: complex1_init
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: complex_init
#endif
        complex(kind=sp),     intent(in) :: c
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = real(c,kind=sp)
        iq.im = aimag(c,kind=sp)
      end function complex1_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function complex2x16_init(c) result(iq) !GCC$ ATTRIBUTES cold :: complex2x16_init !GCC$ ATTRIBUTES vectorcall :: complex2x16_init !GCC$ ATTRIBUTES inline :: complex2x16_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: complex2x16_init
      pure function complex2x16_init(c) result(iq)
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: complex2x16_init
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: complex2x16_init
#endif
        complex(kind=sp), dimension(0:15),  intent(in) :: c
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = real(c,kind=sp)
        iq.im = aimag(c,kind=sp)
      end function complex2x16_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function zmm16r42x_init(v1,v2) result(iq) !GCC$ ATTRIBUTES cold :: zmm16r42x_init !GCC$ ATTRIBUTES vectorcall :: zmm16r42x_init !GCC$ ATTRIBUTES inline :: zmm16r42x_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: zmm16r42x_init
      pure function zmm16r42x_init(v1,v2) result(iq)
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: zmm16r42x_init
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: zmm16r42x_init
#endif
        type(ZMM16r4_t),    intent(in) :: v1
        type(ZMM16r4_t),    intent(in) :: v2
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = v1
        iq.im = v2
      end function zmm16r42x_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function zmm16r41x_init(v) result(iq) !GCC$ ATTRIBUTES cold :: zmm16r41x_init !GCC$ ATTRIBUTES vectorcall :: zmm16r41x_init !GCC$ ATTRIBUTES inline :: zmm16r41x_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: zmm16r41x_init
      pure function zmm16r41x_init(v) result(iq)
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: zmm16r41x_init
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: zmm16r41x_init
#endif
        type(ZMM16r4_t),    intent(in) :: v
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = v
        iq.im = v16_n0
      end function zmm16r41x_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function r41x_init(s) result(iq) !GCC$ ATTRIBUTES cold :: r41x_init !GCC$ ATTRIBUTES inline :: r41x_init
#elif defined __ICC || defined __INTEL_COMPILER
      pure function r41x_init(s) result(iq)
        !DIR$ ATTRIBUTES INLINE :: r41x_init
        !DIR$M ATTRIBUTES CODE_ALIGN : 16 :: r41x_init
#endif
        real(kind=sp),    intent(in) :: s
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = s
        iq.im = v16_n0
      end function r41x_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function copy_init(rhs) result(iq) !GCC$ ATTRIBUTES cold :: copy_init !GCC$ ATTRIBUTES vectorcall :: copy_init !GCC$ ATTRIBUTES inline :: copy_init
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: copy_init
      pure function copy_init(rhs) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: copy_init
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: copy_init
#endif
        type(AVX512c16f32_t),    intent(in) :: rhs
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
         iq = rhs
      end function copy_init

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_add_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_add_c16 !GCC$ ATTRIBUTES vectorcall :: c16_add_c16 !GCC$ ATTRIBUTES inline :: c16_add_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_add_c16
      pure function c16_add_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_add_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_add_c16
#endif
        type(AVX512c16f32_t),    intent(in) :: x
        type(AVX512c16f32_t),    intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re+y.re
        iq.im = x.im+y.im
      end function c16_add_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_add_c1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_add_c1 !GCC$ ATTRIBUTES vectorcall :: c16_add_c1 !GCC$ ATTRIBUTES inline :: c16_add_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_add_c1
      pure function c16_add_c1(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_add_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_add_c1
#endif
        type(AVX512c16f32_t),    intent(in) :: x
        complex(kind=sp),        intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re+real(y,kind=sp)
        iq.im = x.im+aimag(y,kind=sp)
      end function c16_add_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_add_v16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_add_v16 !GCC$ ATTRIBUTES vectorcall :: c16_add_v16 !GCC$ ATTRIBUTES inline :: c16_add_v16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c16_add_v16
      pure function c16_add_v16(x,y) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_add_v16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_add_v16
#endif
        use mod_vecconsts, only : v16_n0
        type(AVX512c16f32_t),    intent(in) :: x
        type(ZMM16r4_t),         intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re+y.v
        iq.im = v16_n0
      end function c16_add_v16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_add_s1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_add_s1 !GCC$ ATTRIBUTES vectorcall :: c16_add_s1 !GCC$ ATTRIBUTES inline :: c16_add_s1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_add_s1
      pure function c16_add_s1(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_add_s1
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_add_s1
#endif
        use mod_veccconts, only : v16_n0
        type(AVX512c16f32_t),    intent(in) :: x
        real(kind=sp),           intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re+y
        iq.im = v16_n0
      end function c16_add_s1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_add_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c1_add_c16 !GCC$ ATTRIBUTES vectorcall :: c1_add_c16 !GCC$ ATTRIBUTES inline :: c1_add_c16
#elif defined __ICC || defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES INLINE :: c1_add_c16
      pure function c1_add_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_add_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c1_add_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = real(x,kind=sp)+y.re
        iq.im = aimag(x,kind=sp)+y.im
      end function c1_add_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function v16_add_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: v16_add_c16 !GCC$ ATTRIBUTES vectorcall :: v16_add_c16 !GCC$ ATTRIBUTES inline :: v16_add_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: v16_add_c16
      pure function v16_add_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: v16_add_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: v16_add_c16
#endif
        use mod_vecconsts, only : v16_n0
        type(ZMM16r4_t),      intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.v+y.re
        iq.im = v16_n0
      end function v16_add_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function s1_add_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: s1_add_c16 !GCC$ ATTRIBUTES vectorcall :: s1_add_c16 !GCC$ ATTRIBUTES inline :: s1_add_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: s1_add_c16
      pure function s1_add_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: s1_add_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: s1_add_c16
#endif
         use mod_vecconsts, only : v16_n0
        real(kind=sp),         intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = s1
        iq.im = v16_n0
      end function s1_add_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_sub_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_sub_c16 !GCC$ ATTRIBUTES vectorcall :: c16_sub_c16 !GCC$ ATTRIBUTES inline :: c16_sub_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_sub_c16
      pure function c16_sub_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_sub_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_sub_c16
#endif
        type(AVX512c16f32_t),     intent(in) :: x
        type(AVX512c16f32_t),     intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re-y.re
        iq.im = x.im-y.im
      end function c16_sub_c16
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER        
      pure function c16_sub_c1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_sub_c1 !GCC$ ATTRIBUTES vectorcall :: c16_sub_c1 !GCC$ ATTRIBUTES inline :: c16_sub_c1
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c16_sub_c1
      pure function c16_sub_c1(x,y) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_sub_c1
         !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_sub_c1
#endif
        type(AVX512c16f32_t),    intent(in) :: x
        complex(kind=sp),        intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re-real(y,kind=sp)
        iq.im = x.im-aimag(y,kind=sp)
      end function c16_sub_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_sub_v16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_sub_v16 !GCC$ ATTRIBUTES vectorcall :: c16_sub_v16 !GCC$ ATTRIBUTES inline :: c16_sub_v16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_sub_v16
      pure function c16_sub_v16(x,y) result(iq)
            use mod_vecconsts, only : v16_n0
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_sub_v16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_sub_v16
#endif
        type(AVX512c16f32_t),     intent(in) :: x
        type(ZMM16r4_t),          intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re-y.v
        iq.im = v16_n0
      end function c16_sub_v16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_sub_s1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_sub_s1 !GCC$ ATTRIBUTES inline :: c16_sub_s1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_sub_s1
      pure function c16_sub_s1(x,y) result(iq)
        use mod_vecconsts, only : v16_n0
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_sub_s1
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_sub_s1
#endif
        type(AVX512c16f32_t),    intent(in) :: x
        real(kind=sp),           intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re-s1
        iq.im = v16_n0
      end function c16_sub_s1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_sub_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c1_sub_c16 !GCC$ ATTRIBUTES vectorcall :: c1_sub_c16 !GCC$ ATTRIBUTES inline :: c1_sub_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_sub_c16
      pure function c1_sub_c16(x,y) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_sub_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c1_sub_c16  
#endif
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = real(x,kind=sp)-y.re
        iq.im = aimag(x,kind=sp)-y.im
      end function c1_sub_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function v16_sub_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: v16_sub_c16 !GCC$ ATTRIBUTES vectorcall :: v16_sub_c16 !GCC$ ATTRIBUTES inline :: v16_sub_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: v16_sub_c16
        pure function v16_sub_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: v16_sub_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: v16_sub_c16
#endif
          use mod_vecconsts, only : v16_n0
          type(ZMM16r4_t),      intent(in) :: x
          type(AVX512c16f32_t), intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.v-y.re
        iq.im = v16_n0
      end function v16_sub_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function s1_sub_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: s1_sub_c16 !GCC$ ATTRIBUTES vectorcall :: s1_sub_c16 !GCC$ ATTRIBUTES inline :: s1_sub_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: s1_sub_c16
      pure function s1_sub_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: s1_sub_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: s1_sub_c16
        
#endif
         use mod_vecconsts, only : v16_n0
        real(kind=sp),           intent(in) :: x
        type(AVX512c16f32_t),    intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x-y.re
        iq.im = v16_n0
      end function s1_sub_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_mul_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_mul_c16 !GCC$ ATTRIBUTES vectorcall :: c16_mul_c16 !GCC$ ATTRIBUTES inline :: c16_mul_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_mul_c16
      pure function c16_mul_c16(x,y) result(uq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_mul_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_mul_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: zmm0,zmm1,zmm2,zmm3
        type(ZMM16r4_t), automatic :: zmm0,zmm1,zmm2,zmm3
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: zmm0 !GCC$ ATTRIBUTES aligned(64) :: zmm0
        type(ZMM16r4_t), automatic :: zmm1 !GCC$ ATTRIBUTES aligned(64) :: zmm1
        type(ZMM16r4_t), automatic :: zmm2 !GCC$ ATTRIBUTES aligned(64) :: zmm2
        type(ZMM16r4_t), automatic :: zmm3 !GCC$ ATTRIBUTES aligned(64) :: zmm3
#endif
       zmm0.v = x.re*y.re
       zmm1.v = x.im*y.im
       iq.re  = zmm0.v+zmm1.v
       zmm2.v = x.im*y.re
       zmm3.v = x.re*y.im
       iq.im  = zmm2.v-zmm3.v
     end function c16_mul_c16
    
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function c16_mul_c1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_mul_c1 !GCC$ ATTRIBUTES vectorcall :: c16_mul_c1 !GCC$ ATTRIBUTES inline :: c16_mul_c1
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: c16_mul_c1
     pure function c16_mul_c1(x,y) result(iq)
       !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_mul_c1
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_mul_c1
#endif
       type(AVX512c16f32_t),     intent(in) :: x
       complex(kind=sp),         intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: zmm0,zmm1,zmm2,zmm3
        type(ZMM16r4_t), automatic :: zmm0,zmm1,zmm2,zmm3
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: zmm0 !GCC$ ATTRIBUTES aligned(64) :: zmm0
        type(ZMM16r4_t), automatic :: zmm1 !GCC$ ATTRIBUTES aligned(64) :: zmm1
        type(ZMM16r4_t), automatic :: zmm2 !GCC$ ATTRIBUTES aligned(64) :: zmm2
        type(ZMM16r4_t), automatic :: zmm3 !GCC$ ATTRIBUTES aligned(64) :: zmm3
#endif
       zmm0.v = x.re*real(y,kind=dp)
       zmm1.v = x.im*aimag(y,kind=dp)
       iq.re  = zmm0.v+zmm1.v
       zmm2.v = x.im*real(y,kind=dp)
       zmm3.v = x.re*aimag(y,kind=dp)
       iq.im  = zmm2.v-zmm3.v
     end function c16_mul_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function c16_mul_v16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_mul_v16 !GCC$ ATTRIBUTES vectorcall :: c16_mul_v16 !GCC$ ATTRIBUTES inline :: c16_mul_v16
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: c16_mul_v16
     pure function c16_mul_v16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_mul_v16
       !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_mul_v16
#endif
       type(AVX512c16f32_t),     intent(in) :: x
       type(ZMM16r4_t),          intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re*y.v
        iq.im = x.im*y.v
     end function c16_mul_v16
       
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function c16_mul_s1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_mul_s1 !GCC$ ATTRIBUTES vectorcall :: c16_mul_s1 !GCC$ ATTRIBUTES inline :: c16_mul_s1
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: c16_mul_s1
     pure function c16_mul_s1(x,y) result(iq)
       !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_mul_s1
       !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_mul_s1
#endif
       type(AVX512c16f32_t), intent(in) :: x
       real(kind=sp),        intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re*y
        iq.im = x.im*y
      end function c16_mul_s1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_mul_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c1_mul_c16 !GCC$ ATTRIBUTES vectorcall :: c1_mul_c16 !GCC$ ATTRIBUTES inline :: c1_mul_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_mul_c16
      pure function c1_mul_c16(x,y) result(iq)
       !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_mul_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_mul_c16
#endif
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: zmm0,zmm1,zmm2,zmm3
        type(ZMM16r4_t), automatic :: zmm0,zmm1,zmm2,zmm3
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: zmm0 !GCC$ ATTRIBUTES aligned(64) :: zmm0
        type(ZMM16r4_t), automatic :: zmm1 !GCC$ ATTRIBUTES aligned(64) :: zmm1
        type(ZMM16r4_t), automatic :: zmm2 !GCC$ ATTRIBUTES aligned(64) :: zmm2
        type(ZMM16r4_t), automatic :: zmm3 !GCC$ ATTRIBUTES aligned(64) :: zmm3
#endif
       zmm0.v = real(x,kind=sp)*y.re
       zmm1.v = aimag(x,kind=sp)*y.im
       iq.re  = zmm0.v+zmm1.v
       zmm2.v = real(x,kind=sp)*y.im
       zmm3.v = aimag(x,kind=sp)*y.re
       iq.im  = zmm2.v-zmm3.v
     end function c1_mul_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function v16_mul_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: v16_mul_c16 !GCC$ ATTRIBUTES vectorcall :: v16_mul_c16 !GCC$ ATTRIBUTES inline :: v16_mul_c16
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: v16_mul_c16
     pure function v16_mul_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: v16_mul_c16
       !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: v16_mul_c16
#endif
       type(ZMM16r4_t),       intent(in) :: x
       type(AVX512c16f32_t),  intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.v*y.re
        iq.im = x.v*y.im
      end function v16_mul_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function s1_mul_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: s1_mul_c16 !GCC$ ATTRIBUTES vectorcall :: s1_mul_c16 !GCC$ ATTRIBUTES inline :: s1_mul_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: s1_mul_c16
        pure function s1_mul_c16(x,y) result(iq)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: s1_mul_c16
          !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: s1_mul_c16
#endif
          real(kind=sp),        intent(in) :: x
          type(AVX512c16f32_t), intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x*y.re
        iq.im = x*y.im
      end function s1_mul_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER     
      pure function c16_div_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_div_c16 !GCC$ ATTRIBUTES vectorcall :: c16_div_c16 !GCC$ ATTRIBUTES inline :: c16_div_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_div_c16
      pure function c16_div_c16(x,y) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_div_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_div_c16
#endif
        type(AVX512c16f32_t),    intent(in) :: x
        type(AVX512c16f32_t),    intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: zmm0,zmm1,zmm2,zmm3,den
        type(ZMM16r4_t), automatic :: zmm0,zmm1,zmm2,zmm3,den
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: zmm0 !GCC$ ATTRIBUTES aligned(64) :: zmm0
        type(ZMM16r4_t), automatic :: zmm1 !GCC$ ATTRIBUTES aligned(64) :: zmm1
        type(ZMM16r4_t), automatic :: zmm2 !GCC$ ATTRIBUTES aligned(64) :: zmm2
        type(ZMM16r4_t), automatic :: zmm3 !GCC$ ATTRIBUTES aligned(64) :: zmm3
        type(ZMM16r4_t), automatic :: den  !GCC$ ATTRIBUTES aligned(64) :: den
#endif
#if (USE_SAFE_COMPLEX_DIVISION) == 1
       iq = cdiv_smith(x,y)
#else
       zmm0.v = x.re*y.re
       zmm1.v = x.im*y.im
       zmm2.v = x.im*y.re
       zmm3.v = x.re*y.im
       den.v  = (y.re*y.re)+(y.im*y.im)
       iq.re  = (zmm0.v+zmm1.v)/den.v
       iq.im  = (zmm2.v-zmm3.v)/den.v
#endif        
     end function c16_div_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function c16_div_c1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_div_c1 !GCC$ ATTRIBUTES vectorcall :: c16_div_c1 !GCC$ ATTRIBUTES inline c16_div_c1
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: c16_div_c1
     pure function c16_div_c1(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_div_c1
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_div_c1
#endif
       type(AVX512c16f32_t),    intent(in) :: x
       complex(kind=sp),        intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: zmm0,zmm1,zmm2,zmm3,den
        type(ZMM16r4_t), automatic :: zmm0,zmm1,zmm2,zmm3,den
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: zmm0 !GCC$ ATTRIBUTES aligned(64) :: zmm0
        type(ZMM16r4_t), automatic :: zmm1 !GCC$ ATTRIBUTES aligned(64) :: zmm1
        type(ZMM16r4_t), automatic :: zmm2 !GCC$ ATTRIBUTES aligned(64) :: zmm2
        type(ZMM16r4_t), automatic :: zmm3 !GCC$ ATTRIBUTES aligned(64) :: zmm3
        type(ZMM16r4_t), automatic :: den  !GCC$ ATTRIBUTES aligned(64) :: den
#endif
       zmm0.v = x.re*real(y,kind=sp)
       zmm1.v = x.im*aimag(y,kind=sp)
       zmm2.v = x.im*real(y,kind=sp)
       zmm3.v = x.re*aimag(y,kind=sp)
       den.v  = (real(y,kind=sp)*real(y,kind=sp))+ &
            (aimag(y,kind=sp)*aimag(y,kind=sp))
       iq.re = (zmm0.v+zmm1.v)/den.v
       iq.im = (zmm2.v-zmm3.v)/den.v
     end function c16_div_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function c16_div_v16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_div_v16 !GCC$ ATTRIBUTES vectorcall :: c16_div_v16 !GCC$ ATTRIBUTES inline :: c16_div_v16
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: c16_div_v16
     pure function c16_div_v16(x,y) result(iq)
       !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_div_v16
       !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_div_v16
#endif
       type(AVX512c16f32_t),     intent(in) :: x
       type(ZMM16r4_t),          intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re/y.v
        iq.im = x.im/y.v
      end function c16_div_v16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_div_s1(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c16_div_s1 !GCC$ ATTRIBUTES vectorcall :: c16_div_s1 !GCC$ ATTRIBUTES inline :: c16_div_s1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_div_s1
      pure function c16_div_s1(x,y) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_div_s1
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: c16_div_s1
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        real(kind=sp),             intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x.re/y
        iq.im = x.im/y
      end function c16_div_s1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_div_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: c1_div_c16 !GCC$ ATTRIBUTES vectorcall :: c1_div_c16 !GCC$ ATTRIBUTES inline :: c1_div_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_div_c16
      pure function c1_div_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_div_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_div_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: zmm0,zmm1,zmm2,zmm3,den
        type(ZMM16r4_t), automatic :: zmm0,zmm1,zmm2,zmm3,den
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: zmm0 !GCC$ ATTRIBUTES aligned(64) :: zmm0
        type(ZMM16r4_t), automatic :: zmm1 !GCC$ ATTRIBUTES aligned(64) :: zmm1
        type(ZMM16r4_t), automatic :: zmm2 !GCC$ ATTRIBUTES aligned(64) :: zmm2
        type(ZMM16r4_t), automatic :: zmm3 !GCC$ ATTRIBUTES aligned(64) :: zmm3
        type(ZMM16r4_t), automatic :: den  !GCC$ ATTRIBUTES aligned(64) :: den
#endif
        real(kind=sp), automatic :: r,i
        r = real(x,kind=sp)
        i = aimag(x,kind=sp)
        zmm0.v = r*y.re
        zmm1.v = i*y.im
        zmm2.v = i*y.re
        zmm3.v = r*y.im
        den.v  = (y.re*y.re)+(y.im*y.im)
        iq.re  = (zmm0.v+zmm1.v)/den.v
        iq.im  = (zmm2.v-zmm3.v)/den.v
      end function c1_div_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function v16_div_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: v16_div_c16 !GCC$ ATTRIBUTES vectorcall :: v16_div_c16 !GCC$ ATTRIBUTES inline :: v16_div_c16
#elif defined __ICC || defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES INLINE :: v16_div_c16
      pure function v16_div_c16(x,y) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: v16_div_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: v16_div_c16
#endif
        type(ZMM16r4_t),       intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x/y.re
        iq.im = x/y.im
      end function v16_div_c16
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function s1_div_c16(x,y) result(iq) !GCC$ ATTRIBUTES hot :: s1_div_c16 !GCC$ ATTRIBUTES vectorcall :: s1_div_c16 !GCC$ ATTRIBUTES inline :: s1_div_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: s1_div_c16
        pure function s1_div_c16(x,y) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: s1_div_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: s1_div_c16
#endif
          real(kind=sp),         intent(in) :: x
          type(AVX512c16f32_t),  intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = x/y.re
        iq.im = x/y.im
      end function s1_div_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function conjugate(x) result(iq) !GCC$ ATTRIBUTES hot :: conjugate !GCC$ ATTRIBUTES vectorcall :: conjugate !DIR$ ATTRIBUTES inline :: conjugate
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: conjugate
      pure function conjugate(x) result(iq)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: conjugate
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: conjugate
#endif
        use mod_vecconsts, only : v16_n0
        type(AVX512c16f32_t),    intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = v16_n0-x.re
        iq.im = x.im
      end function conjugate
        
#if (USE_INTRINSIC_VECTOR_COMPARE) == 1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_eq_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_eq_c16 !GCC$ ATTRIBUTES vectorcall :: c16_eq_c16 !GCC$ ATTRIBUTES inline :: c16_eq_c16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function c16_eq_c16(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_eq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 ::c16_eq_c16
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        type(AVX512c16f32_t),      intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,0)
        lim.zmm = x.im
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,0)
      end function c16_eq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_eq_c1(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_eq_c1 !GCC$ ATTRIBUTES vectorcall :: c16_eq_c1 !GCC$ ATTRIBUTES inline :: c16_eq_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_eq_c1
      pure function c16_eq_c1(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_eq_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_eq_c1
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        type(AVX512c16f32_t),     intent(in) :: x
        complex(kind=sp),         intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = real(y,kind=sp)
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,0)
        lim.zmm = x.im
        rim.zmm = aimag(y,kind=sp)
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,0)
      end function c16_eq_c1
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_eq_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c1_eq_c16 !GCC$ ATTRIBUTES vectorcall :: c1_eq_c16 !GCC$ ATTRIBUTES inline :: c1_eq_c16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c1_eq_c16
      pure function c1_eq_c16(x,y) result(mmask16)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_eq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_eq_c16
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = real(x,kind=sp)
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,0)
        lim.zmm = aimag(x,kind=sp)
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,0)
      end function c1_eq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_neq_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_neq_c16 !GCC$ ATTRIBUTES vectorcall :: c16_neq_c16 !GCC$ ATTRIBUTES inline :: c16_neq_c16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function c16_neq_c16(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_neq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_neq_c16
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        type(AVX512c16f32_t),      intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,12)
        lim.zmm = x.im
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,12)
      end function c16_neq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_neq_c1(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_neq_c1 !GCC$ ATTRIBUTES vectorcall :: c16_neq_c1 !GCC$ ATTRIBUTES inline :: c16_neq_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_neq_c1
      pure function c16_neq_c1(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_neq_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_neq_c1
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        type(AVX512c16f32_t),     intent(in) :: x
        complex(kind=sp),         intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = real(y,kind=sp)
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,12)
        lim.zmm = x.im
        rim.zmm = aimag(y,kind=sp)
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,12)
      end function c16_neq_c1
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_neq_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c1_neq_c16 !GCC$ ATTRIBUTES vectorcall :: c1_neq_c16 !GCC$ ATTRIBUTES inline :: c1_neq_c16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c1_neq_c16
      pure function c1_neq_c16(x,y) result(mmask16)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_neq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_neq_c16
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = real(x,kind=sp)
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,12)
        lim.zmm = aimag(x,kind=sp)
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,12)
      end function c1_neq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_gt_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_gt_c16 !GCC$ ATTRIBUTES vectorcall :: c16_gt_c16 !GCC$ ATTRIBUTES inline :: c16_gt_c16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function c16_gt_c16(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_gt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_gt_c16
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        type(AVX512c16f32_t),      intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,30)
        lim.zmm = x.im
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,30)
      end function c16_gt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_gt_c1(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_gt_c1 !GCC$ ATTRIBUTES vectorcall :: c16_gt_c1 !GCC$ ATTRIBUTES inline :: c16_gt_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_gt_c1
      pure function c16_gt_c1(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_gt_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_gt_c1
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        type(AVX512c16f32_t),     intent(in) :: x
        complex(kind=sp),         intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = real(y,kind=sp)
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,30)
        lim.zmm = x.im
        rim.zmm = aimag(y,kind=sp)
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,30)
      end function c16_gt_c1
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_gt_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c1_gt_c16 !GCC$ ATTRIBUTES vectorcall :: c1_gt_c16 !GCC$ ATTRIBUTES inline :: c1_gt_c16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c1_gt_c16
      pure function c1_gt_c16(x,y) result(mmask16)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_gt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_gt_c16
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = real(x,kind=sp)
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,30)
        lim.zmm = aimag(x,kind=sp)
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,30)
      end function c1_gt_c16
      
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_lt_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_lt_c16 !GCC$ ATTRIBUTES vectorcall :: c16_lt_c16 !GCC$ ATTRIBUTES inline :: c16_lt_c16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function c16_lt_c16(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_lt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_lt_c16
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        type(AVX512c16f32_t),      intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,17)
        lim.zmm = x.im
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,17)
      end function c16_lt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_lt_c1(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_lt_c1 !GCC$ ATTRIBUTES vectorcall :: c16_lt_c1 !GCC$ ATTRIBUTES inline :: c16_lt_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_lt_c1
      pure function c16_lt_c1(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_lt_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_lt_c1
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        type(AVX512c16f32_t),     intent(in) :: x
        complex(kind=sp),         intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = real(y,kind=sp)
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,17)
        lim.zmm = x.im
        rim.zmm = aimag(y,kind=sp)
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,17)
      end function c16_lt_c1
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_lt_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c1_lt_c16 !GCC$ ATTRIBUTES vectorcall :: c1_lt_c16 !GCC$ ATTRIBUTES inline :: c1_lt_c16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c1_eq_c16
      pure function c1_lt_c16(x,y) result(mmask16)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_lt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_lt_c16
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = real(x,kind=sp)
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,17)
        lim.zmm = aimag(x,kind=sp)
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,17)
      end function c1_lt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_ge_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_ge_c16 !GCC$ ATTRIBUTES vectorcall :: c16_ge_c16 !GCC$ ATTRIBUTES inline :: c16_ge_c16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function c16_ge_c16(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_ge_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_ge_c16
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        type(AVX512c16f32_t),      intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,29)
        lim.zmm = x.im
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,29)
      end function c16_ge_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_ge_c1(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_ge_c1 !GCC$ ATTRIBUTES vectorcall :: c16_ge_c1 !GCC$ ATTRIBUTES inline :: c16_ge_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_ge_c1
      pure function c16_ge_c1(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_ge_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_ge_c1
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        type(AVX512c16f32_t),     intent(in) :: x
        complex(kind=sp),         intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = real(y,kind=sp)
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,29)
        lim.zmm = x.im
        rim.zmm = aimag(y,kind=sp)
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,29)
      end function c16_ge_c1
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_ge_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c1_ge_c16 !GCC$ ATTRIBUTES vectorcall :: c1_ge_c16 !GCC$ ATTRIBUTES inline :: c1_ge_c16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c1_ge_c16
      pure function c1_ge_c16(x,y) result(mmask16)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_ge_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_ge_c16
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = real(x,kind=sp)
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,29)
        lim.zmm = aimag(x,kind=sp)
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,29)
      end function c1_ge_c16
      
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_le_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_le_c16 !GCC$ ATTRIBUTES vectorcall :: c16_le_c16 !GCC$ ATTRIBUTES inline :: c16_le_c16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function c16_le_c16(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_le_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_le_c16
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        type(AVX512c16f32_t),      intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,18)
        lim.zmm = x.im
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,18)
      end function c16_le_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_le_c1(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c16_le_c1 !GCC$ ATTRIBUTES vectorcall :: c16_le_c1 !GCC$ ATTRIBUTES inline :: c16_le_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_le_c1
      pure function c16_le_c1(x,y) result(mmask16)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_le_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_le_c1
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        type(AVX512c16f32_t),     intent(in) :: x
        complex(kind=sp),         intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = x.re
        rre.zmm = real(y,kind=sp)
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,18)
        lim.zmm = x.im
        rim.zmm = aimag(y,kind=sp)
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,18)
      end function c16_le_c1
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_le_c16(x,y) result(mmask16) !GCC$ ATTRIBUTES hot :: c1_le_c16 !GCC$ ATTRIBUTES vectorcall :: c1_le_c16 !GCC$ ATTRIBUTES inline :: c1_le_c16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: c1_le_c16
      pure function c1_le_c16(x,y) result(mmask16)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_le_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_le_c16
#endif
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
        complex(kind=sp),     intent(in) :: x
        type(AVX512c16f32_t), intent(in) :: y
        integer(c_short), dimension(0:1) :: mmask16
#if defined __INTEL_COMPILER 
        
        !DIR$ ATTRIBUTES ALIGN : 64 :: lre,lim,rre,rim
        type(v16f32), automatic :: lre,lim,rre,rim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        
        type(v16f32), automatic :: lre !GCC$ ATTRIBUTES aligned(64) :: lre
        type(v16f32), automatic :: lim !GCC$ ATTRIBUTES aligned(64) :: lim
        type(v16f32), automatic :: rre !GCC$ ATTRIBUTES aligned(64) :: rre
        type(v16f32), automatic :: rim !GCC$ ATTRIBUTES aligned(64) :: rim
        
#endif
        mmask16 = 0
        lre.zmm = real(x,kind=sp)
        rre.zmm = y.re
        mmask16(0) = v16f32_cmp_ps_mask(lre,rre,18)
        lim.zmm = aimag(x,kind=sp)
        rim.zmm = y.im
        mmask16(1) = v16f32_cmp_ps_mask(lim,rim,18)
      end function c1_le_c16
      

#else

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_eq_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_eq_c16 !GCC$ ATTRIBUTES vectorcall :: c16_eq_c16 !GCC$ ATTRIBUTES inline :: c16_eq_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_eq_c16
      pure function c16_eq_c16(x,y) result(bres)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_eq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_eq_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre = .false.
        mre = (x.re == y.re)
        mim = .false.
        mim = (x.im == y.im)
        bres = .false.
        bres(0) = all(mre)
        bres(1) = all(mim)
      end function c16_eq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_eq_c1(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_eq_c1 !GCC$ ATTRIBUTES vectorcall :: c16_eq_c1 !GCC$ ATTRIBUTES inline :: c16_eq_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_eq_c1
      pure function c16_eq_c1(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_eq_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_eq_c1  
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        complex(kind=sp),          intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre  = .false.
        mim  = .false.
        bres = .false.
        mre = x.re == real(y,kind=sp)
        bres(0) = all(mre)
        mim = x.im == aimag(y,kind=sp)
        bres(1) = all(mim)
      end function c16_eq_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_eq_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c1_eq_c16 !GCC$ ATTRIBUTES vectorcall :: c1_eq_c16 !GCC$ ATTRIBUTES inline :: c1_eq_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_eq_c16
      pure function c1_eq_c16(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_eq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_eq_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
       mre  = .false.
       mim  = .false.
       bres = .false.
       mre = real(x,kind=sp)  == y.re
       bres(0) = all(mre)
       mim = aimag(x,kind=sp) == y.im
       bres(1) = all(mim)
     end function c1_eq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_neq_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_neq_c16 !GCC$ ATTRIBUTES vectorcall :: c16_neq_c16 !GCC$ ATTRIBUTES inline :: c16_neq_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_neq_c16
      pure function c16_eq_c16(x,y) result(bres)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_neq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_neq_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre = .false.
        mre = (x.re /= y.re)
        mim = .false.
        mim = (x.im /= y.im)
        bres = .false.
        bres(0) = all(mre)
        bres(1) = all(mim)
      end function c16_neq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_neq_c1(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_neq_c1 !GCC$ ATTRIBUTES vectorcall :: c16_neq_c1 !GCC$ ATTRIBUTES inline :: c16_neq_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_neq_c1
      pure function c16_eq_c1(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_neq_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_neq_c1  
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        complex(kind=sp),          intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre  = .false.
        mim  = .false.
        bres = .false.
        mre = (x.re /= real(y,kind=sp))
        bres(0) = all(mre)
        mim = (x.im /= aimag(y,kind=sp))
        bres(1) = all(mim)
      end function c16_neq_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_neq_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c1_neq_c16 !GCC$ ATTRIBUTES vectorcall :: c1_neq_c16 !GCC$ ATTRIBUTES inline :: c1_neq_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_neq_c16
      pure function c1_eq_c16(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_neq_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_neq_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
       mre  = .false.
       mim  = .false.
       bres = .false.
       mre = (real(x,kind=sp)  /= y.re)
       bres(0) = all(mre)
       mim = (aimag(x,kind=sp) /= y.im)
       bres(1) = all(mim)
     end function c1_neq_c16


#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_gt_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_gt_c16 !GCC$ ATTRIBUTES vectorcall :: c16_gt_c16 !GCC$ ATTRIBUTES inline :: c16_gt_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_gt_c16
      pure function c16_gt_c16(x,y) result(bres)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_gt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_gt_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre = .false.
        mre = (x.re > y.re)
        mim = .false.
        mim = (x.im > y.im)
        bres = .false.
        bres(0) = all(mre)
        bres(1) = all(mim)
      end function c16_gt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_gt_c1(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_gt_c1 !GCC$ ATTRIBUTES vectorcall :: c16_gt_c1 !GCC$ ATTRIBUTES inline :: c16_gt_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_gt_c1
      pure function c16_gt_c1(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_gt_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_gt_c1  
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        complex(kind=sp),          intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre  = .false.
        mim  = .false.
        bres = .false.
        mre = (x.re > real(y,kind=sp))
        bres(0) = all(mre)
        mim = (x.im > aimag(y,kind=sp))
        bres(1) = all(mim)
      end function c16_gt_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_gt_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c1_gt_c16 !GCC$ ATTRIBUTES vectorcall :: c1_gt_c16 !GCC$ ATTRIBUTES inline :: c1_gt_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_gt_c16
      pure function c1_gt_c16(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_gt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_gt_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
       mre  = .false.
       mim  = .false.
       bres = .false.
       mre = (real(x,kind=sp) > y.re)
       bres(0) = all(mre)
       mim = (aimag(x,kind=sp) > y.im)
       bres(1) = all(mim)
     end function c1_gt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_lt_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_lt_c16 !GCC$ ATTRIBUTES vectorcall :: c16_lt_c16 !GCC$ ATTRIBUTES inline :: c16_lt_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_lt_c16
      pure function c16_lt_c16(x,y) result(bres)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_lt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_lt_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre = .false.
        mre = (x.re < y.re)
        mim = .false.
        mim = (x.im < y.im)
        bres = .false.
        bres(0) = all(mre)
        bres(1) = all(mim)
      end function c16_lt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_lt_c1(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_lt_c1 !GCC$ ATTRIBUTES vectorcall :: c16_lt_c1 !GCC$ ATTRIBUTES inline :: c16_lt_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_lt_c1
      pure function c16_lt_c1(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_lt_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_lt_c1  
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        complex(kind=sp),          intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre  = .false.
        mim  = .false.
        bres = .false.
        mre = (x.re < real(y,kind=sp))
        bres(0) = all(mre)
        mim = (x.im < aimag(y,kind=sp))
        bres(1) = all(mim)
      end function c16_lt_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_lt_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c1_lt_c16 !GCC$ ATTRIBUTES vectorcall :: c1_lt_c16 !GCC$ ATTRIBUTES inline :: c1_lt_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_lt_c16
      pure function c1_lt_c16(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_lt_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_lt_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
       mre  = .false.
       mim  = .false.
       bres = .false.
       mre = (real(x,kind=sp) < y.re)
       bres(0) = all(mre)
       mim = (aimag(x,kind=sp) < y.im)
       bres(1) = all(mim)
     end function c1_lt_c16


#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_ge_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_ge_c16 !GCC$ ATTRIBUTES vectorcall :: c16_ge_c16 !GCC$ ATTRIBUTES inline :: c16_ge_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_ge_c16
      pure function c16_ge_c16(x,y) result(bres)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_ge_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_ge_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre = .false.
        mre = (x.re >= y.re)
        mim = .false.
        mim = (x.im >= y.im)
        bres = .false.
        bres(0) = all(mre)
        bres(1) = all(mim)
      end function c16_ge_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_ge_c1(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_ge_c1 !GCC$ ATTRIBUTES vectorcall :: c16_ge_c1 !GCC$ ATTRIBUTES inline :: c16_ge_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_ge_c1
      pure function c16_ge_c1(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_ge_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_ge_c1  
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        complex(kind=sp),          intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre  = .false.
        mim  = .false.
        bres = .false.
        mre = (x.re >= real(y,kind=sp))
        bres(0) = all(mre)
        mim = (x.im >= aimag(y,kind=sp))
        bres(1) = all(mim)
      end function c16_ge_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_ge_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c1_ge_c16 !GCC$ ATTRIBUTES vectorcall :: c1_ge_c16 !GCC$ ATTRIBUTES inline :: c1_ge_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_ge_c16
      pure function c1_eq_c16(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_ge_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_ge_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
       mre  = .false.
       mim  = .false.
       bres = .false.
       mre = (real(x,kind=sp)  >= y.re)
       bres(0) = all(mre)
       mim = (aimag(x,kind=sp) >= y.im)
       bres(1) = all(mim)
     end function c1_eq_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_le_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_le_c16 !GCC$ ATTRIBUTES vectorcall :: c16_le_c16 !GCC$ ATTRIBUTES inline :: c16_le_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_le_c16
      pure function c16_le_c16(x,y) result(bres)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_le_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_le_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre = .false.
        mre = (x.re <= y.re)
        mim = .false.
        mim = (x.im <= y.im)
        bres = .false.
        bres(0) = all(mre)
        bres(1) = all(mim)
      end function c16_le_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c16_le_c1(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c16_le_c1 !GCC$ ATTRIBUTES vectorcall :: c16_le_c1 !GCC$ ATTRIBUTES inline :: c16_le_c1
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c16_le_c1
      pure function c16_le_c1(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c16_le_c1
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c16_le_c1  
#endif
        type(AVX512c16f32_t),      intent(in) :: x
        complex(kind=sp),          intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
        logical(kind=int1), dimension(0:1) :: bres
        mre  = .false.
        mim  = .false.
        bres = .false.
        mre = (x.re <= real(y,kind=sp))
        bres(0) = all(mre)
        mim = (x.im <= aimag(y,kind=sp))
        bres(1) = all(mim)
      end function c16_le_c1

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function c1_le_c16(x,y) result(bres) !GCC$ ATTRIBUTES hot :: c1_le_c16 !GCC$ ATTRIBUTES vectorcall :: c1_le_c16 !GCC$ ATTRIBUTES inline :: c1_le_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: c1_le_c16
      pure function c1_le_c16(x,y) result(bres)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: c1_le_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: c1_le_c16
#endif
        complex(kind=sp),      intent(in) :: x
        type(AVX512c16f32_t),  intent(in) :: y
#if defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: mre,mim
        logical(kind=int4), dimension(0:15) :: mre,mim
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        logical(kind=int4), dimension(0:15) :: mre !GCC$ ATTRIBUTES aligned(64) :: mre
        logical(kind=int4), dimension(0:15) :: mim !GCC$ ATTRIBUTES aligned(64) :: mim
#endif
       mre  = .false.
       mim  = .false.
       bres = .false.
       mre = (real(x,kind=sp)  <= y.re)
       bres(0) = all(mre)
       mim = (aimag(x,kind=sp) <= y.im)
       bres(1) = all(mim)
     end function c1_le_c16




        
      

#endif

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function polar(rho,theta) result(iq) !GCC$ ATTRIBUTES hot :: polar !GCC$ ATTRIBUTES vectorcall :: polar !GCC$ ATTRIBUTES inline :: polar
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: polar
       pure function polar(rho,theta)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: polar
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: polar
#endif
         type(ZMM16r4_t),     intent(in) :: rho
         type(ZMM16r4_t),     intent(in) :: theta
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t) :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t) :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = rho*cos(theta.v)
        iq.im = rho*sin(theta.v)
      end function polar

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function carg_c16(x) result(arg) !GCC$ ATTRIBUTES hot :: carg_c16 !GCC$ ATTRIBUTES vectorcall :: carg_c16 !GCC$ inline :: carg_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: carg_c16
      pure function carg_c16(x)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: carg_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: carg_c16
#endif
        type(AVX512c16f32_t),  intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: arg
        type(ZMM16r4_t) :: arg
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(ZMM16r4_t)  :: arg !GCC$ ATTRIBUTES aligned(64) :: arg
#endif
        arg.v = atan2(x.re,x.im)
      end function carg_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function carg_2xv16(re,im) result(arg) !GCC$ ATTRIBUTES hot :: carg_2xv16 !GCC$ ATTRIBUTES vectorcall :: carg_2xv16 !GCC$ ATTRIBUTES inline :: carg_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: carg_2xv16
        pure function carg_2xv16(re,im) result(arg)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: carg_2xv16
          !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: carg_2xv16
#endif
          type(ZMM16r4_t),     intent(in) :: re
          type(ZMM16r4_t),     intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: arg
        type(ZMM16r4_t) :: arg
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(ZMM16r4_t)  :: arg !GCC$ ATTRIBUTES aligned(64) :: arg
#endif
        arg.v = atan2(re,im)
      end function carg_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function csin_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: csin_c16 !GCC$ ATTRIBUTES vectorcall :: csin_c16 !GCC$ ATTRIBUTES inline :: csin_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: csin_c16
        pure function csin_c16(x) result(iq)
           !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: csin_c16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: csin_c16
#endif
          type(AVX512c16f32_t),    intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = sin(x.re)*cosh(x.im)
        iq.im = cos(x.re)*sinh(x.im)
      end function csin_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function csin_2xv16(re,im)  result(iq) !GCC$ ATTRIBUTES hot :: csin_2xv16 !GCC$ ATTRIBUTES vectorcall :: csin_2xv16 !GCC$ ATTRIBUTES inline :: csin_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: csin_2xv16
      pure function csin_2xv16(re,im) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: csin_2xv16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: csin_2xv16
#endif
        type(ZMM16r4_t),    intent(in) :: re
        type(ZMM16r4_t),    intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = sin(re.v)*cosh(im.v)
        iq.im = cos(re.v)*sinh(im.v)
      end function csin_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function csinh_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: csin_c16 !GCC$ ATTRIBUTES vectorcall :: csin_c16 !GCC$ ATTRIBUTES inline :: csin_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: csinh_c16
      pure function csinh_c16(x) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: csinh_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: csinh_c16
#endif
        type(AVX512c16f32_t),    intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = sinh(x.re)*cos(x.im)
        iq.im = cosh(x.re)*sin(x.im)
      end function csinh_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function csinh_2xv16(re,im) result(iq) !GCC$ ATTRIBUTES hot :: csinh_2xv16 !GCC$ ATTRIBUTES vectorcall :: csinh_2xv16 !GCC$ ATTRIBUTES inline :: csinh_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: csinh_2xv16
      pure function csinh_2xv16(re,im) result(iq)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: csinh_2xv16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: csinh_2xv16
#endif
        type(ZMM16r4_t),     intent(in) :: re
        type(ZMM16r4_t),     intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = sinh(re.v)*cos(im.v)
        iq.im = cosh(re.v)*sin(im.v)
      end function csinh_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ccos_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: ccos_c16 !GCC$ ATTRIBUTES vectorcall :: ccos_c16 !GCC$ ATTRIBUTES inline :: ccos_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: ccos_c16
      pure function ccos_c16(x) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ccos_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ccos_c16
#endif
        type(AVX512c16f32_t),   intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = cos(x.re)*cosh(x.im)
        iq.im = sin(x.re)*sinh(x.im)
      end function ccos_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ccos_2xv16(re,im) result(iq) !GCC$ ATTRIBUTES hot :: ccos_2xv16 !GCC$ ATTRIBUTES vectorcall :: ccos_2xv16 !GCC$ ATTRIBUTES inline :: ccos_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES INLINE :: ccos_2xv16
      pure function ccos_2xv16(re,im) result(iq)
        !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ccos_2xv16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ccos_2xv16
#endif
        type(ZMM16r4_t),    intent(in) :: re
        type(ZMM16r4_t),    intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = cos(re)*cosh(im)
        iq.im = sin(re)*sinh(im)
      end function ccos_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ccosh_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: ccosh_c16 !GCC$ ATTRIBUTES vectorcall :: ccosh_c16 !GCC$ ATTRIBUTES inline :: ccosh_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: ccosh_c16
      pure function ccosh_c16(x) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ccosh_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ccosh_c16
#endif
        type(AVX512c16f32_t),     intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = cosh(x.re)*cos(x.im)
        iq.im = sinh(x.re)*sin(x.im)
      end function ccosh_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ccosh_2xv16(re,im) result(iq) !GCC$ ATTRIBUTES hot :: ccosh_2xv16 !GCC$ ATTRIBUTES vectorcall :: ccosh_2xv16 !GCC$ ATTRIBUTES inline :: ccosh_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: ccosh_2xv16
      pure function ccos_2xv16(re,im) result(iq)
           !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ccosh_2xv16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ccosh_2xv16
#endif
        type(ZMM16r4_t),    intent(in) :: re
        type(ZMM16r4_t),    intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = cosh(re.v)*cos(im.v)
        iq.im = sinh(re.v)*sin(im.v)
      end function ccos_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cexp_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: cexp_c16 !GCC$ ATTRIBUTES vectorcall :: cexp_c16 !GCC$ ATTRIBUTES inline :: cexp_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: cexp_c16
      pure function cexp_c16(x) result(iq)
            !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cexp_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cexp_c16
#endif
          type(AVX512c16f32_t),   intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = exp(x.re)*cos(x.im)
        iq.im = exp(x.re)*cos(x.im)
      end function cexp_c16
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ctan_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: ctan_c16 !GCC$ ATTRIBUTES vectorcall :: ctan_c16 !GCC$ ATTRIBUTES inline :: ctan_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: ctan_c16
      pure function ctan_c16(x) result(iq)
              !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ctan_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ctan_c16
#endif
        type(AVX512c16f32_t),    intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq = csin_c16(x)/ccos_c16(x)
      end function ctan_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ctan_2xv16(re,im) result(iq) !GCC$ ATTRIBUTES hot :: ctan_2xv16 !GCC$ ATTRIBUTES vectorcall :: ctan_2xv16 !GCC$ ATTRIBUTES inline :: ctan_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE ::  ctan_2xv16
      pure function ctan_2xv16(re,im) result(iq)
                !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ctan_2xv16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ctan_2xv16
#endif
        type(ZMM16r4_t),     intent(in) :: re
        type(ZMM16r4_t),     intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq = csin_2xv16(re,im)/ccos_2xv16(re,im)
      end function ctan_2xv16
        
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function ctanh_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: ctanh_c16 !GCC$ ATTRIBUTES vectorcall :: ctanh_c16 !GCC$ ATTRIBUTES inline :: ctanh_c16
#elif defined __ICC || defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES INLINE :: ctanh_c16
        pure function ctanh_c16(x) result(iq)
             !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: ctanh_c16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ctanh_c16
#endif
          type(AVX512c16f32_t),     intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq = csinh_c16(x)/ccosh_c16(x)
      end function ctanh_c16
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cexp_2xv16(re,im) result(iq) !GCC$ ATTRIBUTES hot :: cexp_2xv16 !GCC$ ATTRIBUTES vectorcall :: cexp_2xv16 !GCC$ ATTRIBUTES inline :: cexp_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
      pure function cexp_2xv16(re,im) result(iq)
           !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cexp_2xv16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cexp_2xv16
#endif
        type(ZMM16r4_t),     intent(in) :: re
        type(ZMM16r4_t),     intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
#endif
        iq.re = exp(re.v)*cos(im.v)
        iq.im = exp(re.v)*sin(im.v)
      end function cexp_2xv16
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cabs_c16(x) result(val) !GCC$ ATTRIBUTES hot :: cabs_c16 !GCC$ ATTRIBUTES vectorcall :: cabs_c16 !GCC$ ATTRIBUTES inline :: cabs_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: cabs_c16
        pure function cabs_c16(x) result(val)
           !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cabs_c16
          !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: cabs_c16
#endif
          type(AVX512c16f32_t),  intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: val
        type(ZMM16r4_t)  :: val
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(ZMM16r4_t)  :: val !GCC$ ATTRIBUTES aligned(64) :: val
#endif
          val.v = sqrt(x.re*x.re+x.im*x.im)
      end function cabs_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cabs_2xv16(re,im) result(val) !GCC$ ATTRIBUTES hot :: cabs_2xv16 !GCC$ ATTRIBUTES vectorcall :: cabs_2xv16 !GCC$ ATTRIBUTES inline :: cabs_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: cabs_2xv16
        pure function cabs_2xv16(re,im) result(val)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cabs_2xv16
          !DIR$ ATTRIBUTES CODE_ALIGN : 16 :: cabs_2xv16
#endif
          type(ZMM16r4_t),    intent(in) :: re
          type(ZMM16r4_t),    intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: val
        type(ZMM16r4_t)  :: val
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(ZMM16r4_t)  :: val !GCC$ ATTRIBUTES aligned(64) :: val
#endif
        val.v = sqrt(re*re+im*im)
      end function cabs_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cpow_c16(c16,n) result(iq) !GCC$ ATTRIBUTES hot :: cpow_c16 !GCC$ ATTRIBUTES vectorcall :: cpow_c16 !GCC$ ATTRIBUTES inline :: cpow_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: cpow_c16
      pure function cpow_c16(x,n) result(iq)
         !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cpow_c16
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cpow_c16
#endif
        type(AVX512c16f32_t),    intent(in) :: x
        real(kind=sp),           intent(in) :: n
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: r,theta,pow,trig
        type(ZMM16r4_t), automatic :: r,theta,pow,trig
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: r !GCC$ ATTRIBUTES aligned(64) :: r
        type(ZMM16r4_t), automatic :: theta !GCC$ ATTRIBUTES aligned(64) :: theta
        type(ZMM16r4_t), automatic :: pow !GCC$ ATTRIBUTES aligned(64) :: pow
        type(ZMM16r4_t), automatic :: trig !GCC$ ATTRIBUTES aligned(64) :: trig
#endif
        r.v = sqrt(x.re*x.re+x.im*x.im)
        pow.v = r.v**n
        theta.v = atan(x.im/x.re)
        trig.v = theta.v*n
        iq.re = pow.v*cos(trig.v)
        iq.im = pow.v*sin(trig.v)
      end function cpow_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cpow_2xv16(re,im,n) result(iq)  !GCC$ ATTRIBUTES hot :: cpow_2xv16 !GCC$ ATTRIBUTES vectorcall :: cpow_2xv16 !GCC$ ATTRIBUTES inline :: cpow_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES INLINE :: cpow_2xv16
        pure function cpow_2xv16(re,im,n)
             !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cpow_2xv16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cpow_2xv16
#endif
          type(ZMM16r4_t),   intent(in) :: re
          type(ZMM16r4_t),   intent(in) :: im
          real(kind=sp),     intent(in) :: n
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: r,theta,pow,trig
        type(ZMM16r4_t), automatic :: r,theta,pow,trig
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: r !GCC$ ATTRIBUTES aligned(64) :: r
        type(ZMM16r4_t), automatic :: theta !GCC$ ATTRIBUTES aligned(64) :: theta
        type(ZMM16r4_t), automatic :: pow !GCC$ ATTRIBUTES aligned(64) :: pow
        type(ZMM16r4_t), automatic :: trig !GCC$ ATTRIBUTES aligned(64) :: trig
#endif
        r.v = sqrt(re*re+im*im)
        pow.v = r.v**n
        theta.v = atan(im/re)
        trig.v = theta.v*n
        iq.re = pow.v*cos(trig.v)
        iq.im = pow.v*sin(trig.v)
      end function cpow_2xv16
        
#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function clog_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: clog_c16 !GCC$ ATTRIBUTES vectorcall :: clog_c16 !GCC$ ATTRIBUTES inline :: clog_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: clog_c16
        pure function clog_c16(x) result(iq)
    !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: clog_c16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: clog_c16
#endif
          type(AVX512c16f32_t),  intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: t0
        type(ZMM16r4_t), automatic :: t0
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: t0 !GCC$ ATTRIBUTES aligned(64) :: t0
       
#endif
       t0 = cabs_c16(x)
       iq.re = log(t0.v)
       iq.im = carg_c16(x)
     end function clog_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
     pure function clog_2xv16(re,im) result(iq)  !GCC$ ATTRIBUTES hot :: clog_2xv16 !GCC$ ATTRIBUTES vectorcall :: clog_2xv16 !GCC$ ATTRIBUTES inline :: clog_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: clog_2xv16
       pure function clog_2xv16(re,im) result(iq)
          !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: clog_2xv16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: clog_2xv16
#endif
         type(ZMM16r4_t),  intent(in)  :: re
         type(ZMM16r4_t),  intent(in)  :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: t0
        type(ZMM16r4_t), automatic :: t0
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: t0 !GCC$ ATTRIBUTES aligned(64) :: t0
       
#endif
        t0 = cabs_2xv16(re,im)
        iq.re = log(t0.v)
        iq.im = carg_2xv16(re,im)
      end function clog_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function csqrt_c16(x) result(iq) !GCC$ ATTRIBUTES hot :: csqrt_c16 !GCC$ ATTRIBUTES vectorcall :: csqrt_c16 !GCC$ ATTRIBUTES inline :: csqrt_c16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: csqrt_c16
        pure function csqrt_c16(x) result(iq)
            !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: csqrt_c16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: csqrt_c16
#endif
          use mod_vecconsts, only : v16_1over2
          type(AVX512c16f32_t),     intent(in) :: x
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: t0,t1,t2
        type(ZMM16r4_t), automatic :: t0,t1,t2
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: t0 !GCC$ ATTRIBUTES aligned(64) :: t0
        type(ZMM16r4_t), automatic :: t1 !GCC$ ATTRIBUTES aligned(64) :: t1
        type(ZMM16r4_t), automatic :: t2 !GCC$ ATTRIBUTES aligned(64) :: t2
       
#endif
        t0 = cabs_c16(x)
        t1.v = v16_1over2.v*(t0.v+c16.re)
        iq.re = sqrt(t1.v)
        t2.v = v16_1over2*(t0.v-c16.re)
        iq.im = sqrt(t2.v)
      end function csqrt_c16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function csqrt_2xv16(re,im) result(iq) !GCC$ ATTRIBUTES hot :: csqrt_2xv16 !GCC$ ATTRIBUTES vectorcall :: csqrt_2xv16 !GCC$ ATTRIBUTES inline :: csqrt_2xv16
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: csqrt_2xv16
        pure function csqrt_2xv16(re,im) result(iq)
            !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: csqrt_2xv16
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: csqrt_2xv16
#endif
          use mod_vecconsts, only : v16_1over2
          type(ZMM16r4_t),  intent(in) :: re
          type(ZMM16r4_t),  intent(in) :: im
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: t0,t1,t2
        type(ZMM16r4_t), automatic :: t0,t1,t2
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: t0 !GCC$ ATTRIBUTES aligned(64) :: t0
        type(ZMM16r4_t), automatic :: t1 !GCC$ ATTRIBUTES aligned(64) :: t1
        type(ZMM16r4_t), automatic :: t2 !GCC$ ATTRIBUTES aligned(64) :: t2
       
#endif
        t0 = cabs_2xv16(re,im)
        t1.v = v16_1over2*(t0.v+re.v)
        iq.re = sqrt(t1.v)
        t2.v = v8_1over2*(t0.v-re.v)
      end function csqrt_2xv16

#if defined __GFORTRAN__ && !defined __INTEL_COMPILER
      pure function cdiv_smith(x,y) result(iq) !GCC$ ATTRIBUTES hot :: cdiv_smith !GCC$ ATTRIBUTES vectorcall :: cdiv_smith !GCC$ ATTRIBUTES inline :: cdiv_smith
#elif defined __ICC || defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: cdiv_smith
        pure function cdiv_smith(x,y) result(iq)
             !DIR$ ATTRIBUTES VECTOR:PROCESSOR(skylake_avx512) :: cdiv_smith
          !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cdiv_smith
#endif
#if  (USE_INTRINSIC_VECTOR_COMPARE) == 1
        use, intrinsic :: ISO_C_BINDING
        use mod_avx512_bindings, only : v16f32, v16f32_cmp_ps_mask
#endif
        type(AVX512c16f32_t),   intent(in) :: x
        type(AVX512c16f32_t),   intent(in) :: y
#if defined __INTEL_COMPILER 
        !DIR$ ATTRIBUTES ALIGN : 64 :: iq
        type(AVX512c16f32_t)  :: iq
        !DIR$ ATTRIBUTES ALIGN : 64 :: ratio,denom
        type(ZMM16r4_t), automatic :: ratio,denom
#elif defined __GFORTRAN__ && !defined __INTEL_COMPILER
        type(AVX512c16f32_t)  :: iq !GCC$ ATTRIBUTES aligned(64) :: iq
        type(ZMM16r4_t), automatic :: ratio !GCC$ ATTRIBUTES aligned(64) :: ratio
        type(ZMM16r4_t), automatic :: denom !GCC$ ATTRIBUTES aligned(64) :: denom
       
       
#endif
        logical(kind=int4), dimension(0:15) :: bres
#if   (USE_INTRINSIC_VECTOR_COMPARE) == 1
        !DIR$ ATTRIBUTES ALIGN : 64 :: tre,tim
        type(v16f32) :: tre,tim
        integer(c_short), automatic :: mask_gte
        integer(c_short), parameter :: all_ones = Z'FFFF'
#endif
#if   (USE_INTRINSIC_VECTOR_COMPARE) == 1
        tre.zmm = abs(y.re)
        tim.zmm = abs(y.im)
        mask_gte = v16f32_cmp_ps_mask(tre,tim,29);
#else
        bres = abs(y.re) >= abs(y.im)
#endif
        re_part = v16_n0
        im_part = v16_n0
#if   (USE_INTRINSIC_VECTOR_COMPARE) == 1
        if(mask_gte == all_ones) then
#elif
        if(all(bres)) then
#endif
           ratio.v   = y.im/y.re
           denom.v   = y.re+(ratio.v*y.im)
           iq.re     = (x.re+x.im*ratio.v)/denom.v
           iq.im     = (x.im-x.re*ratio.v)/denom.v
        else
           ratio.v   = y.re/y.im
           denom.v   = y.im+ratio.v*y.re
           iq.re     = (x.re*ratio.v+x.im)/denom.v
           iq.im     = (x.im*ratio.v-x.re)/denom.v
        end if       
      end function cdiv_smith
        
      
end module mod_avx512c16f32
