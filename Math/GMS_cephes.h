
#ifndef __GMS_CEPHES_H__
#define __GMS_CEPHES_H__

/*
    Slightly adapted version of CEPHES library.
    The rationale: removed the overhead of GLIBC calls for math.h declared
    mathematical functions (elementary and transcendental).
*/

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "GSM_config.h"

#if !defined (DENORMAL)
    #define DENORMAL 1
#endif

#define EXPMSK 0x807f
#define MEXP 255
#define NBITS 24

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_ceilf( const float x ){
     float y = 0.0f;
     y = floorf( (float )x );
     if( y < x )
	 y += 1.0;
     return(y);
}

/* Bit clearing masks: */

static unsigned short bmask[] = {
       0xffff,
       0xfffe,
       0xfffc,
       0xfff8,
       0xfff0,
       0xffe0,
       0xffc0,
       0xff80,
       0xff00,
       0xfe00,
       0xfc00,
       0xf800,
       0xf000,
       0xe000,
       0xc000,
       0x8000,
       0x0000,
    };

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_floorf( float x ) {
unsigned short *p;
union
  {
    float y;
    unsigned short i[2];
  } u;
int e;
u.y = x;
/* find the exponent (power of 2) */
p = &u.i[1];
e = (( *p >> 7) & 0xff) - 0x7f;
p -= 1;
if( e < 0 ){
    if( u.y < 0 )
	 return( -1.0 );
    else
	 return( 0.0 );
}
e = (NBITS -1) - e;
/* clean out 16 bits at a time */
while( e >= 16 ) {
	*p++ = 0;
	 e -= 16;
}

/* clear the remaining bits */
if( e > 0 )
	*p &= bmask[e];

if( (x < 0) && (u.y != x) )
	u.y -= 1.0;

   return(u.y);
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_frexpf( float x, int *pw2 ) {
union
  {
    float y;
    unsigned short i[2];
  } u;
int i, k;
short *q;
u.y = x;
q = &u.i[1];
/* find the exponent (power of 2) */

i  = ( *q >> 7) & 0xff;
if( i == 0 ) {
    if( u.y == 0.0 ) {
	*pw2 = 0;
	 return(0.0);
    }
/* Number is denormal or zero */
#if DENORMAL
/* Handle denormal number. */
	do{
		
	   u.y *= 2.0;
	   i -= 1;
	   k  = ( *q >> 7) & 0xff;
	}
	while( k == 0 );
	i = i + k;
#else
	*pw2 = 0;
	return( 0.0 );
#endif /* DENORMAL */
	}
i -= 0x7e;
*pw2 = i;
*q &= 0x807f;	/* strip all exponent bits */
*q |= 0x3f00;	/* mantissa between 0.5 and 1 */
return( u.y );
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_ldexpf( float x, int pw2 ) {
union
  {
    float y;
    unsigned short i[2];
  } u;
short *q;
int e;
u.y = x;
q = &u.i[1];
while( (e = ( *q >> 7) & 0xff) == 0 ) {
	if( u.y == (float )0.0 ) {
		return( 0.0 );
	}
/* Input is denormal. */
	if( pw2 > 0 ) {
	    u.y *= 2.0;
	    pw2 -= 1;
	}
	if( pw2 < 0 ) {
	    if( pw2 < -24 )
		return( 0.0 );
	    u.y *= 0.5;
	    pw2 += 1;
	}
	if( pw2 == 0 )
		return(u.y);
	}
e += pw2;

/* Handle overflow */
if( e > MEXP )
	{
	return( 3.4028235e+38 );
	}

*q &= 0x807f;

/* Handle denormalized results */
if( e < 1 )
	{
#if DENORMAL
	if( e < -24 )
		return( 0.0 );
	*q |= 0x80; /* Set LSB of exponent. */
	/* For denormals, significant bits may be lost even
	   when dividing by 2.  Construct 2^-(1-e) so the result
	   is obtained with only one multiplication.  */
	u.y *= ldexpf(1.0f, e - 1);
	return(u.y);
#else
	return( 0.0 );
#endif
	}
*q |= (e & 0xff) << 7;
return(u.y);
}

/* Return 1 if the sign bit of x is 1, else 0.  */

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
int ceph_signbitf(const float x) {
union
	{
	float f;
	short s[4];
	int i;
	} u;

u.f = x;

if( sizeof(int) == 4 ){
    return( u.i < 0 );
}
else{
     return( u.s[1] < 0 );
   }
	
}

/* Return 1 if x is a number that is Not a Number, else return 0.  */

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
int ceph_isnanf(const float x) {
union
	{
	float f;
	unsigned short s[2];
	unsigned int i;
	} u;

u.f = x;
if( sizeof(int) == 4 ) {
	
     if( ((u.i & 0x7f800000) == 0x7f800000)
	    && ((u.i & 0x007fffff) != 0) )
		return 1;

	return(0);
}
else
    { /* size int not 4 */

	if( (u.s[1] & 0x7f80) == 0x7f80) {
	     if( ((u.s[1] & 0x007f) | u.s[0]) != 0 )
			return(1);
	}

	return(0);
     } /* size int not 4 */

}

/* Return 1 if x is not infinite and is not a NaN.  */
__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
int ceph_isfinitef(const float x) {
union{
	
	float f;
	unsigned short s[2];
	unsigned int i;
} u;
u.f = x;
if( sizeof(int) == 4 ) {
	if( (u.i & 0x7f800000) != 0x7f800000)
	   return 1;

	return(0);
}
else
    {

	if( (u.s[1] & 0x7f80) != 0x7f80)
	     return 1;

	return(0);
     }

}

/* Single precision inverse hyperbolic cosine
 * test interval: [1.0, 1.5]
 * trials: 10000
 * peak relative error: 1.7e-7
 * rms relative error: 5.0e-8
 *
 * Copyright (C) 1989 by Stephen L. Moshier.  All rights reserved.
 */

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float ceph_acoshf(const double xx) {
float x, z;
x = xx;
if( x < 1.0 ) {
    return(0.0);
}
if( x > 1500.0 )
   return( ceph_logf(x) + 1.44269504088896341  );

z = x - 1.0;
if( z < 0.5 ) {
    z =
	(((( 1.7596881071E-3 * z
	  - 7.5272886713E-3) * z
	  + 2.6454905019E-2) * z
	  - 1.1784741703E-1) * z
	  + 1.4142135263E0) * ceph_sqrtf( z );
}
else
    {
	z = ceph_sqrtf( z*(x+1.0) );
	z = ceph_logf(x + z);
}
return( z );
}

/* Single precision circular arcsine
 * test interval: [-0.5, +0.5]
 * trials: 10000
 * peak relative error: 6.7e-8
 * rms relative error: 2.5e-8
 */
__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float ceph_asinf(const double xx) {
float a, x, z;
int sign, flag;
x = xx;
if( x > 0 ) {
  sign = 1;
  a = x;
}
else
    {
  sign = -1;
  a = -x;
}

if( a > 1.0 ) {
    return( 0.0 );
}

if( a < 1.0e-4 ) {
    z = a;
    goto done;
}

if( a > 0.5 ) {
    z = 0.5 * (1.0 - a);
    x = sqrtf( z );
    flag = 1;
}
else
    {
    x = a;
    z = x * x;
    flag = 0;
}
z =
(((( 4.2163199048E-2 * z
  + 2.4181311049E-2) * z
  + 4.5470025998E-2) * z
  + 7.4953002686E-2) * z
  + 1.6666752422E-1) * z * x
  + x;

if( flag != 0 ) {
    z = z + z;
    z = 1.5707963267948966192 - z;
}
done:
if( sign < 0 )
	z = -z;
return( z );
}


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float ceph_acosf(const float x ) {
if( x < -1.0 )
	goto domerr;

if( x < -0.5) 
	return( 3.14159265358979323846
                - 2.0 * ceph_asinf(ceph_sqrtf(0.5*(1.0+x)) ) );
if( x > 1.0 ) {
domerr:
	return( 0.0 );
}
if( x > 0.5 )
	return( 2.0 * ceph_asinf(  ceph_sqrtf(0.5*(1.0-x) ) ) );
	
return(1.5707963267948966192 - ceph_asinf(x) );
}

/* Single precision circular arcsine
 * test interval: [-tan(pi/8), +tan(pi/8)]
 * trials: 10000
 * peak relative error: 7.7e-8
 * rms relative error: 2.9e-8
 */
 
__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float ceph_atanf(const float xx ) {
float x, y, z;
int sign;
x = xx;

/* make argument positive and save the sign */
if( xx < 0.0 ) {	
    sign = -1;
    x = -xx;
}
else
    {
     sign = 1;
     x = xx;
}
/* range reduction */
if( x > 2.414213562373095 ) {  /* tan 3pi/8 */
    y = 1.5707963267948966192;
    x = -( 1.0/x );
}
else if( x > 0.4142135623730950 ) { /* tan pi/8 */
	y = 0.7853981633974483096;
	x = (x-1.0)/(x+1.0);
}
else
    y = 0.0;

z = x * x;
y +=
((( 8.05374449538e-2 * z
  - 1.38776856032E-1) * z
  + 1.99777106478E-1) * z
  - 3.33329491539E-1) * z * x
  + x;
if( sign < 0 )
	y = -y;
return( y );
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float atan2f(const float y,
             const float x ) {
float z, w;
int code;
code = 0;
if( x < 0.0 )
	code = 2;
if( y < 0.0 )
	code |= 1;
if( x == 0.0 ) {
      if( code & 1 ){
	  return( -1.5707963267948966192 );
      }
      if( y == 0.0 )
	     return( 0.0 );
      return(1.5707963267948966192 );
}

if( y == 0.0 ) {
	 if( code & 2 )
		return(3.14159265358979323846);
	 return( 0.0 );
}

switch( code ) {
	
	default:
#if ANSIC
	case 0:
	case 1: w = 0.0; break;
	case 2: w = 3.14159265358979323846; break;
	case 3: w = -3.14159265358979323846; break;
}
z = ceph_atanf( y/x );
return( w + z );
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float atanhf(const double xx) {
float x, z;
x = xx;
if( x < 0 )
	z = -x;
else
	z = x;
if( z >= 1.0 ) {
	if( x == 1.0 )
		return(3.4028234663852885981170418348451692544e38);
	if( x == -1.0 )
		return( -3.4028234663852885981170418348451692544e38);
   	return(3.4028234663852885981170418348451692544e38);
}

if( z < 1.0e-4 )
	return(x);

if( z < 0.5 ) {
	z = x * x;
	z =
	(((( 1.81740078349E-1 * z
	  + 8.24370301058E-2) * z
	  + 1.46691431730E-1) * z
	  + 1.99782164500E-1) * z
	  + 3.33337300303E-1) * z * x
	  + x;
}
else
    {
	z = 0.5 * ceph_logf( (1.0+x)/(1.0-x) );
}
return( z );
}



__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float cbrtf(const double xx) {
int e, rem, sign;
float x, z;
x = xx;
if( x == 0 )
	return( 0.0 );
if( x > 0 )
	sign = 1;
else{
    sign = -1;
    x = -x;
}

z = x;
/* extract power of 2, leaving
 * mantissa between 0.5 and 1
 */
x = ceph_frexpf( x, &e );

/* Approximate cube root of number between .5 and 1,
 * peak relative error = 9.2e-6
 */
x = (((-0.13466110473359520655053  * x
      + 0.54664601366395524503440 ) * x
      - 0.95438224771509446525043 ) * x
      + 1.1399983354717293273738  ) * x
      + 0.40238979564544752126924;

/* exponent divided by 3 */
if( e >= 0 ) {
    rem = e;
    e /= 3;
    rem -= 3*e;
    if( rem == 1 )
	x *= 1.25992104989487316477;
    else if( rem == 2 )
		x *= 1.58740105196819947475;
}
/* argument less than 1 */
else {
	
   e = -e;
   rem = e;
   e /= 3;
   rem -= 3*e;
   if( rem == 1 )
       x /= 1.25992104989487316477;
   else if( rem == 2 )
       x /= 1.58740105196819947475;
       e = -e;
}
/* multiply by power of 2 */
x = ceph_ldexpf( x, e );
/* Newton iteration */
x -= ( x - (z/(x*x)) ) * 0.333333333333;
if( sign < 0 )
	x = -x;
return(x);
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline 
float ceph_coshf(const double xx) {
float x, y;
x = xx;
if( x < 0 )
	x = -x;
if( x > 88.72283905206835)
    return( 3.4028234663852885981170418348451692544e38);
y = ceph_expf(x);
y = y + 1.0/y;
return( 0.5*y );
}




#endif /*__GMS_CEPHES_H__*/
