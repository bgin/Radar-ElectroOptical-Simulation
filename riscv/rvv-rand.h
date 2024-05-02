
#ifndef __RVV_RAND_H__
#define __RVV_RAND_H__

#include <riscv_vector.h>

/*Various statistical CDF and PDF functional approximations and random
  samples belonging to these functions.
  This is used as a random input data generator for unit and performance tests
  of pixman-rvv implementation.
*/

/*
    !*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static void  
rvv_i4_uniform(int * __restrict__    seed,
               int * __restrict__    randi,
               int                   n) {
    
        int * __restrict__    ps = seed;
        int * __restrict__    pr = randi;
        size_t                vn = (size_t)n;
        size_t vl;
        for(vn > 0;vn -= vl;ps += vl;pr += vl) {
            vl = __riscv_vsetvl_e32m1(vn);
            register vint32m1_t v0 = __riscv_vle32_v_i32m1(ps,vl);
            register vint32m1_t vk = __riscv_vdiv_vx_i32m1(v0,127773,vl);
            register vint32m1_t v1 = __riscv_vmul_vx_i32m1(k,2836,vl);
            register vint32m1_t v2 = __riscv_vsub_vv_i32m1(v0,__riscv_vmul_vx_i32m1(vk,127773,vl),vl);
            v0                     = __riscv_vsub_vv_i32m1(__riscv_vmul_vx_i32m1(v2,16807,vl),v1,vl);
            vbool32_t           vb = __riscv_vmslt_vx_i32m1_b32(v0,0,vl);
            if(__riscv_vfirst_m_b32()) v0 = __riscv_vadd_vx_i32m1(v0,2147483647,vl);
            __riscv_vse32_v_i32m1(pr,v0,vl);
        }
}


__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static void  
rvv_r4_uniform(const int * __restrict__ randi,
               float     * __restrict__ randr,
               int                      n) {
        
        const int * __restrict__ pri = randi;
        float     * __restrict__ prr = randr;
        size_t                   vn  = (size_t)n;
        size_t                   vl;
        for(vn > 0;vn -= vl;pri += vl;prr += vl) {
            vl = __riscv_vsetvl_e32m1(vn);
            register vint32m1_t   iv0 = __riscv_vle32_v_i32m1(pri,vl);
            register vfloat32m1_t fv0 = __riscv_vreinterpret_v_i32m1_f32m1(iv0);
            register vfloat32m1_t fv1 = __riscv_vfmul_vf_v32m1(fv0,4.656612875e-10f,vl);
            __riscv_vse32_v_f32m1(prr,fv1,vl);
        }
}


/*
    
     !*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!

*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static void
rvv_i4_uniform_ab(const int * __restrict__ a,
                  const int * __restrict__ b,
                  int * __restrict__ seed,
                  int * __restrict__ randi,
                  const int          n) {
        
        const vfloat32m1_t vone     = __riscv_vfmv_v_f_f32m1(1.0f,4);
        const int * __restrict__ pa = a;
        const int * __restrict__ pb = b;
        int       * __restrict__ ps = seed;
        int       * __restrict__ pr = randi;
        size_t                   vn = (size_t)n;
        size_t                   vl;
        for(vn > 0;vn -= vl;pa += vl;pb += vl;ps += vl;pr += vl) {
            vl = __riscv_vsetvl_e32m1(vn);
            register vint32m1_t v0   = __riscv_vle32_v_i32m1(ps,vl);
            register vint32m1_t va   = __riscv_vle32_v_i32m1(pa,vl);
            register vint32m1_t vk   = __riscv_vdiv_vx_i32m1(v0,127773,vl);
            register vint32m1_t vb   = __riscv_vle32_v_i32m1(pb,vl);
            register vint32m1_t v1   = __riscv_vmul_vx_i32m1(k,2836,vl);
            register vfloat32m1_t vfa= __riscv_vreinterpret_v_i32m1_f32m1(__riscv_vmin_vv_i32m1(va,vb,vl));
            register vint32m1_t v2   = __riscv_vsub_vv_i32m1(v0,__riscv_vmul_vx_i32m1(vk,127773,vl),vl);
            register vfloat32m1_t vfb= __riscv_reinterpret_v_i32m1_f32m1(__riscv_vmax_vv_i32m1(va,vb,vl));
            v0                       = __riscv_vsub_vv_i32m1(__riscv_vmul_vx_i32m1(v2,16807,vl),v1,vl);
            vbool32_t           vb   = __riscv_vmslt_vx_i32m1_b32(v0,0,vl);
            register vfloat32m1_t vf1= __riscv_vfsub_vf_f32m1(vfa,0.5f,vl);
            if(__riscv_vfirst_m_b32()) v0 = __riscv_vadd_vx_i32m1(v0,2147483647,vl);
            register vfloat32m1_t vf2= __riscv_vfadd_vf_f32m1(vfb,0.5f,vl);
            register vfloat32m1_t vr = __riscv_vfmul_vf_f32m1(__riscv_vreinterpret_v_i32m1_f32m1(v0),4.656612875e-10f,vl);
            vr                       = __riscv_vfmadd_vv_f32m1(__riscv_vfsub_vv_f32m1(vone,vr,vl),
                                                               vf1,__riscv_vfmul_vv_f32m1(vr,vf2),vl);
            register vint32m1_t   val= __riscv_vreinterpret_v_f32m1_i32m1(vr);
            val                      = __riscv_vmax_vv_i32m1(val,__riscv_vmin_vv_i32m1(va,vb,vl),vl);
            val                      = __riscv_vmin_vv_i32m1(val,__riscv_vmax_vv_i32m1(va,vb,vl),vl);
            __riscv_vse32_v_i32m1(pr,val,vl);
        }        
}


/*
   Various xorshift rvv-vectorized simple implementations.
   Relying on scalar C-version as implemented by: https://en.wikipedia.org/wiki/Xorshift
*/

/* The state must be initialized to non-zero */
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static void
rvv_xorshift32(unsigned int * __restrict__ state,
               int                         n) {
    
        unsigned int * __restrict__ ps = state;
        size_t                      vn = (size_t)n;
        size_t                      vl;
        for(vn > 0;vn -= vl;ps += vl) {
            vl = __riscv_vsetvl_e32m1(vn);
            register vuint32m1_t vi0 = __riscv_vle32_v_u32m1(ps,vl);
            register vuint32m1_t vx  = vi0;
            vx                       = __riscv_vxor_vv_u32m1(vx,__riscv_vsll_vx_u32m1(vx,13,vl),vl);
            vx                       = __riscv_vxor_vv_u32m1(vx,__riscv_vsrl_vx_u32m1(vx,17,vl),vl);
            vx                       = __riscv_vxor_vv_u32m1(vx,__riscv_vsll_vx_u32m1(vx,5,vl),vl);
            __riscv_vse32_v_u32m1(ps,vx,vl);
        }
}


/*Non-linear xorshift variants*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static void
rvv_xorwow(unsigned int * __restrict__ x0,
           unsigned int * __restrict__ x1,
           unsigned int * __restrict__ x2,
           unsigned int * __restrict__ x3,
           unsigned int * __restrict__ x4,
           unsigned int * __restrict__ ct,
           unsigned int * __restrict__ randi,
           int                         n) {

        unsigned int * __restrict__ px0 = x0;
        unsigned int * __restrict__ px1 = x1;
        unsigned int * __restrict__ px2 = x2;
        unsigned int * __restrict__ px3 = x3;
        unsigned int * __restrict__ px4 = x4;
        unsigned int * __restrict__ pct = ct;
        unsigned int * __restrict__ pr  = randi;
        size_t                      vn  = (size_t)n;
        size_t                      vl;
        for(vn>0;vn-=vl;px0+=vl;px1+=vl;px2+=vl;px3+=vl;px4+=vl;ct+=vl;pr+=vl) {
            vl = __riscv_vsetvl_e32m1(vn);
            register vuint32m1_t vt3= __riscv_vle32_v_u32m1(px3,vl);
            register vuint32m1_t vt2= __riscv_vle32_v_u32m1(px2,vl);
            register vuint32m1_t vt1= __riscv_vle32_v_u32m1(px1,vl);
            register vuint32m1_t vt = __riscv_vle32_v_u32m1(px4,vl);
            register vuint32m1_t vs = __riscv_vle32_v_u32m1(px0,vl);
            __riscv_vse32_v_u32m1(px4,vt3,vl);
            __riscv_vse32_v_u32m1(px3,vt2,vl);
            __riscv_vse32_v_u32m1(px2,vt1,vl);
            __riscv_vse32_v_u32m1(px1,vs,vl);
            vt                      = __riscv_vxor_vv_u32m1(vt,__riscv_vsrl_vx_u32m1(vt,2,vl),vl);
            vt                      = __riscv_vxor_vv_u32m1(vt,__riscv_vsll_vx_u32m1(vt,1,vl),vl);
            vt                      = __riscv_vxor_vv_u32m1(vt,__riscv_vxor_vv_u32m1(vs,__riscv_vsll_vx_u32m1(vs,4,vl),vl),vl);
            register vuint32m1_t vc = __riscv_vle32_v_u32m1(pct,vl);
            __riscv_vse32_v_u32m1(px0,vt,vl);
            register vuint32m1_t vc2=__riscv_vadd_vx_u32m1(vc,362437,vl);
            __riscv_vse32_v_u32m1(pct,vc2,vl);
            __riscv_vse32_v_u32m1(pr,__riscv_vadd_vv_u32m1(vt,vc2,vl),vl);
        }
}

/*
    LFSR113, LFSR258: adapted from this paper: http://lomont.org/papers/2008/Lomont_PRNG_2008.pdf
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static vuint32m1_t
rvv_lfsr113(vuint32m1_t * __restrict__ z1,
            vuint32m1_t * __restrict__ z2,
            vuint32m1_t * __restrict__ z3,
            vuint32m1_t * __restrict__ z4,
            size_t                     vl) {
    
        vuint32m1_t b;
        b     = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(*(z1),6,vl),*(z1),vl),13,vl);
        *(z1) = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(*(z1),4294967294,vl),18,vl),b,vl);
        b     = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(*(z2),2,vl),*(z2),vl),27,vl);
        *(z2) = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(*(z2),4294967288,vl),2,vl),b,vl);
        b     = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(*(z3),13,vl),*(z3),vl),21,vl);
        *(z3) = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(*(z3),4294967280,vl),7,vl),b,vl);
        b     = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(*(z4),3,vl),*(z4),vl),12,vl);
        *(z4) = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(*(z4),4294967168,vl),13,vl),b,vl);
        return (__riscv_vxor_vv_u32m1(__riscv_vxor_vv_u32m1(*(z1),*(z2),vl),
                                      __riscv_vxor_vv_u32m1(*(z3),*(z4),vl),vl));
}


__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static void
rvv_lfsr113_v2(uint32_t * __restrict__ z1,
               uint32_t * __restrict__ z2,
               uint32_t * __restrict__ z3,
               uint32_t * __restrict__ z4,
               uint32_t * __restrict__ randi,
               int                     vl) {

        uint32_t * __restrict__ pz1 = z1;
        uint32_t * __restrict__ pz2 = z2;
        uint32_t * __restrict__ pz3 = z3;
        uint32_t * __restrict__ pz4 = z4;
        uint32_t * __restrict__ pr  = randi;
        size_t                  vn  = (size_t)n;
        size_t                  vl;
        for(vn > 0;vn-=vl;pz1+=vl;pz2+=vl;pz3+=vl;pz4+=vl;pr+=vl) {
            vl = __riscv_vsetvl_e32m1(vn);
            vuint32m1_t vz1 = __riscv_vle32_v_u32m1(pz1,vl);
            vuint32m1_t b   = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(vz1,6,vl),vz1,vl),13,vl);
            vz1             = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(vz1,4294967294,vl),18,vl),b,vl);
            vuint32m1_t vz2 = __riscv_vle32_v_u32m1(pz2,vl);
            b               = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(vz2,2,vl),vz2,vl),27,vl);
            vz2             = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(vz2,4294967288,vl),2,vl),b,vl);
            vuint32m1_t vz3 = __riscv_vle32_v_u32m1(pz3,vl);
            b               = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(vz3,13,vl),vz3,vl),21,vl);
            vz3             = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(vz3,4294967280,vl),7,vl),b,vl);
            vuint32m1_t vz4 = __riscv_vle32_v_u32m1(pz4,vl);
            b               = __riscv_vsrl_vx_u32m1(__riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(vz4,3,vl),vz4,vl),12,vl);
            vz4             = __riscv_vxor_vv_u32m1(__riscv_vsll_vx_u32m1(__riscv_vand_vx_u32m1(vz4,4294967168,vl),13,vl),b,vl);
            register vuint32m1_t vres = __riscv_vxor_vv_u32m1(__riscv_vxor_vv_u32m1(vz1,vz2,vl),
                                        __riscv_vxor_vv_u32m1(vz3,vz4,vl),vl);
            __riscv_vse32_v_u32m1(pr,vres,vl);
        }
}


/* The original code from http://www.burtleburtle.net/bob/rand/smallprng.html */

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
__attribute__((riscv_vector_cc))
inline 
static vuint32m1_t
rvv_rot(const vuint32m1_t x,
    const uint32_t k,
    size_t         vl) {
  
    const uint32_t val = (uint32_t)32-k;
    register vuint32m1_t vt0 = __riscv_vsrl_vx_u32m1(x,val,vl);
    register vuint32m1_t vt1 = __riscv_vsll_vx_u32m1(x,k,vl);
    register vuint32m1_t vt2 = __riscv_vor_vv_u32m1(vt0,vt1,vl);
    return (vt2);
}

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline 
static vuint32m1_t
rvv_ranval(vuint32m1_t * __restrict__ a,
           vuint32m1_t * __restrict__ b,
           vuint32m1_t * __restrict__ c,
           vuint32m1_t * __restrict__ d,
           size_t                     vl) {
        
    register vuint32m1_t e = __riscv_vsub_vv_u32m1(*(a),rvv_rot(*(b),27,vl),vl);
    *(a)                   = __riscv_vxor_vv_u32m1(*(b),rvv_rot(*(c),17,vl),vl);
    *(b)                   = __riscv_vadd_vv_u32m1(*(c),*(d),vl);
    *(c)                   = __riscv_vadd_vv_u32m1(*(d),e,vl);
    *(d)                   = __riscv_vadd_vv_u32m1(e,*(a),vl);
    return (*(d));
}

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline static
void rvv_raninit(  vuint32m1_t * __restrict__ a,
                   vuint32m1_t * __restrict__ b,
                   vuint32m1_t * __restrict__ c,
                   vuint32m1_t * __restrict__ d,
                   const vuint32m1_t          seed, 
                   size_t                     vl) {
        
        *(a) = __riscv_vmv_v_x_u32m1(0xf1ea5eed,vl);
        *(b) = *(c) = *(d) = seed;
        int32_t i;
        for(i = 0; i<20; ++i) (void)rvv_ranval(a,b,c,d,vl);
}


/* This is xoshiro128++ 1.0, one of our 32-bit all-purpose, rock-solid
   generators. It has excellent speed, a state size (128 bits) that is
   large enough for mild parallelism, and it passes all tests we are aware
   of.

   For generating just single-precision (i.e., 32-bit) floating-point
   numbers, xoshiro128+ is even faster.

   The state must be seeded so that it is not everywhere zero. */

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline static
vuint32m1_t rvv_xorshiro128_pp(vuint32m1_t * __restrict__ s0,
                     vuint32m1_t * __restrict__ s1,
                     vuint32m1_t * __restrict__ s2,
                     vuint32m1_t * __restrict__ s3,
                     size_t                     vl) {
        
        register vuint32m1_t x0 = *(s0);
        register vuint32m1_t x1 = *(s1);
        register vuint32m1_t x2 = *(s2);
        register vuint32m1_t x3 = *(s3);
        register vuint32m1_t result = __riscv_vadd_vv_u32m1(rvv_rot(__riscv_vadd_vv_u32m1(x0,x3,vl),7,vl),x0,vl);
        register vuint32m1_t  t = __riscv_vsll_vx_u32m1(x1,9,vl);
        x2                      = __riscv_vxor_vv_u32m1(x2,x0,vl);
        x1                      = __riscv_vxor_vv_u32m1(x1,x2,vl);
        x3                      = __riscv_vxor_vv_u32m1(x3,x1,vl);
        x0                      = __riscv_vxor_vv_u32m1(x0,x3,vl);
        x2                      = __riscv_vxor_vv_u32m1(x2,t,vl);
        x3                      = rvv_rot(x3,11,vl);
        return (result);
}


/* This is xoshiro128** 1.1, one of our 32-bit all-purpose, rock-solid
   generators. It has excellent speed, a state size (128 bits) that is
   large enough for mild parallelism, and it passes all tests we are aware
   of.

   Note that version 1.0 had mistakenly s[0] instead of s[1] as state
   word passed to the scrambler.

   For generating just single-precision (i.e., 32-bit) floating-point
   numbers, xoshiro128+ is even faster.

   The state must be seeded so that it is not everywhere zero. */


__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
inline static
vuint32m1_t rvv_xorshiro128_ss(vuint32m1_t * __restrict__ s0,
                               vuint32m1_t * __restrict__ s1,
                               vuint32m1_t * __restrict__ s2,
                               vuint32m1_t * __restrict__ s3,
                               size_t                     vl) {

        register vuint32m1_t x0     = *(s0);
        register vuint32m1_t x1     = *(s1);
        register vuint32m1_t x2     = *(s2);
        register vuint32m1_t x3     = *(s3);
        register vuint32m1_t result = __riscv_vmul_vx_u32m1(rvv_rot(__riscv_vmul_vx_u32m1(x1,5,vl),7,vl),9,vl);
        register vuint32m1_t t      = __riscv_vsll_vx_u32m1(x1,9,vl);
        x2                          = __riscv_vxor_vv_u32m1(x2,x0,vl);
        x3                          = __riscv_vxor_vv_u32m1(x3,x1,vl);
        x1                          = __riscv_vxor_vv_u32m1(x1,x2,vl);
        x0                          = __riscv_vxor_vv_u32m1(x0,x3,vl);
        x2                          = __riscv_vxor_vv_u32m1(x2,t,vl);
        x3                          = rvv_rot(x3,11,vl);
        return (result);
}


// https://en.wikipedia.org/wiki/Permuted_congruential_generator



__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
__attribute__((riscv_vector_cc))
inline 
static vuint32m1_t
rvv_pcg32_fast(const vuint64m1_t seed,
               size_t vl) {
    
    static vuint64m1_t   mcg_state  = __riscv_vadd_vx_u64m1(__riscv_vadd_vv_u64m1(seed,seed,vl),1,vl);
    register vuint64m1_t multiplier = __riscv_vmv_v_x_u64m1(6364136223846793005u,vl);
    register vuint64m1_t x          = mcg_state;
    register vuint32m1_t count      = __riscv_vreinterpret_v_u64m1_u32m1(__riscv_vsrl_vx_u64m1(x,61,vl));
    mcg_state                       = __riscv_vmul_vv_u64m1(x,multiplier,vl);
    x                               = __riscv_vxor_vv_u64m1(x,__riscv_vsrl_vx_u64m1(x,22,vl),vl);
    return (__riscv_vreinterpret_v_u64m1_u32m1(__riscv_vsrl_vx_u64m1(x,22+count,vl)));
}



    























#endif /*__RVV_RAND_H__*/