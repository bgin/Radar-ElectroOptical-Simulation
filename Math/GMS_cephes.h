
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

static int sgngamf;

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


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_expf(const double xx) {
float x, z;
int n;
x = xx;
if( x > 88.72283905206835) {
    return(3.4028234663852885981170418348451692544e38);
}
if( x < -103.278929903431851103) {
    return(0.0);
}

/* Express e**x = e**g 2**n
 *   = e**g e**( n loge(2) )
 *   = e**( g + n loge(2) )
 */
z = ceph_floorf(1.44269504088896341 * x + 0.5 ); /* floor() truncates toward -infinity. */
x -= z * 0.693359375;
x -= z * -2.12194440e-4;
n = z;

z = x * x;
/* Theoretical peak relative error in [-0.5, +0.5] is 4.2e-9. */
z =
((((( 1.9875691500E-4  * x
   + 1.3981999507E-3) * x
   + 8.3334519073E-3) * x
   + 4.1665795894E-2) * x
   + 1.6666665459E-1) * x
   + 5.0000001201E-1) * z
   + x
   + 1.0;

/* multiply by power of 2 */
x = ceph_ldexpf( z, n );
return( x );
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_expx2f(const float x) {
  float u, u1, m;
  if (x < 0.0)
    x = -x;

  /* Represent x as an exact multiple of 1/32 plus a residual.  */
  m = .03125f * ceph_floorf(32.0f * x + 0.5f);
  x -= m;
  /* x**2 = m**2 + 2mf + f**2 */
  u = m * m;
  u1 = 2 * m * x  +  x * x;

  if ((u+u1) > 88.72283905206835 )
    return (3.4028234663852885981170418348451692544e38);

  /* u is exact, u1 is small.  */
  u = ceph_expf(u) * ceph_expf(u1);
  return(u);
}


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float logf( const double xx ) {
register float y;
float x, z, fe;
int e;
x = xx;
fe = 0.0;
/* Test for domain */
if( x <= 0.0 ) {
    if( x == 0.0 )
	   ;
    else
       return(-103.278929903431851103);
}
x = ceph_frexpf( x, &e );
if( x < 0.707106781186547524) {
    e -= 1;
    x = x + x - 1.0; /*  2x - 1  */
}	
else{
    x = x - 1.0;
}
z = x * x;
/* 3.4e-9 */
/*
p = logfcof;
y = *p++ * x;
for( i=0; i<8; i++ )
	{
	y += *p++;
	y *= x;
	}
y *= z;
*/

y =
(((((((( 7.0376836292E-2 * x
- 1.1514610310E-1) * x
+ 1.1676998740E-1) * x
- 1.2420140846E-1) * x
+ 1.4249322787E-1) * x
- 1.6668057665E-1) * x
+ 2.0000714765E-1) * x
- 2.4999993993E-1) * x
+ 3.3333331174E-1) * x * z;

if( e ){
   fe = e;
   y += -2.12194440e-4 * fe;
}

y +=  -0.5 * z;  /* y - 0.5 x^2 */
z = x + y;   /* ... + x  */
if( e )
	z += 0.693359375 * fe;
return( z );
}

#define ceph_fabsf(x) ( (x) < 0 ? -(x) : (x) )

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float betaf(const  float aa, const float bb ) {
float a, b, y;
int sign;
sign = 1;
a = aa;
b = bb;
if( a <= 0.0 ) {
     if( a == ceph_floorf(a) )
		goto over;
}
if( b <= 0.0 ) {
     if( b == ceph_floorf(b) )
		goto over;
}

y = a + b;
if( ceph_fabsf(y) > 34.84425627277176174) {
	y = ceph_lgamf(y);
	sign *= sgngamf; /* keep track of the sign */
	y = ceph_lgamf(b) - y;
	sign *= sgngamf;
	y = ceph_lgamf(a) + y;
	sign *= ceph_sgngamf;
	if( y > 88.72283905206835) {
		
over:
	    return( sign * 3.4028234663852885981170418348451692544e38 );
	 }
      return( sign * ceph_expf(y) );
}

y = ceph_gammaf(y);
if( y == 0.0 )
	goto over;

if( a > b ) {
    y =  ceph_gammaf(a)/y;
    y *= ceph_gammaf(b);
}
else{
    y = ceph_gammaf(b)/y;
    y *= ceph_gammaf(a);
}
return(y);
}


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_chbevlf( float x, float * __restrict array, int n ) {
float b0, b1, b2, *p;
int i;
p = array;
b0 = *p++;
b1 = 0.0;
i = n - 1;
do{
   	b2 = b1;
	b1 = b0;
	b0 = x * b1  -  b2  + *p++;
}
while( --i );
return( 0.5*(b0-b2) );
}

/*
   DESCRIPTION:
 *
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 *
 * This is accomplished using the inverse gamma integral
 * function and the relation
 *
 *    x/2 = igami( df/2, y );
*/



__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_chdtrcf(const float dff, const float xx){
float df, x;
df = dff;
x = xx;
if( (x < 0.0) || (df < 1.0) ) {
    return(0.0);
}
return( ceph_igamcf( 0.5*df, 0.5*x ) );
}


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float chdtrf(const float dff, const float xx) {
float df, x;
df = dff;
x = xx;
if( (x < 0.0) || (df < 1.0) ) {
	return(0.0);
}
return( ceph_igamf( 0.5*df, 0.5*x ) );
}


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float chdtrif( const float dff, const float yy ) {
float y, df, x;
y = yy;
df = dff;
if( (y < 0.0) || (y > 1.0) || (df < 1.0) ) {
	return(0.0);
}
x = ceph_igamif( 0.5 * df, y );
return( 2.0 * x );
}

/*
DESCRIPTION:
 *
 * Approximates the integral
 *
 *                             x
 *                             -
 *                      2     | |        2
 *  dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
 *                          | |
 *                           -
 *                           0
 *
 * Three different rational approximations are employed, for
 * the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,10        50000       4.4e-7      6.3e-8
*/



//extern float PIF, MACHEPF;



__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float dawsnf( const float xxx ) {
/* Dawson's integral, interval 0 to 3.25 */
const static float AN[10] = {
 1.13681498971755972054E-11,
 8.49262267667473811108E-10,
 1.94434204175553054283E-8,
 9.53151741254484363489E-7,
 3.07828309874913200438E-6,
 3.52513368520288738649E-4,
-8.50149846724410912031E-4,
 4.22618223005546594270E-2,
-9.17480371773452345351E-2,
 9.99999999999999994612E-1,
};
const static float AD[11] = {
 2.40372073066762605484E-11,
 1.48864681368493396752E-9,
 5.21265281010541664570E-8,
 1.27258478273186970203E-6,
 2.32490249820789513991E-5,
 3.25524741826057911661E-4,
 3.48805814657162590916E-3,
 2.79448531198828973716E-2,
 1.58874241960120565368E-1,
 5.74918629489320327824E-1,
 1.00000000000000000539E0,
};

/* interval 3.25 to 6.25 */
const static float BN[11] = {
 5.08955156417900903354E-1,
-2.44754418142697847934E-1,
 9.41512335303534411857E-2,
-2.18711255142039025206E-2,
 3.66207612329569181322E-3,
-4.23209114460388756528E-4,
 3.59641304793896631888E-5,
-2.14640351719968974225E-6,
 9.10010780076391431042E-8,
-2.40274520828250956942E-9,
 3.59233385440928410398E-11,
};
const static float BD[10] = {
/*  1.00000000000000000000E0,*/
-6.31839869873368190192E-1,
 2.36706788228248691528E-1,
-5.31806367003223277662E-2,
 8.48041718586295374409E-3,
-9.47996768486665330168E-4,
 7.81025592944552338085E-5,
-4.55875153252442634831E-6,
 1.89100358111421846170E-7,
-4.91324691331920606875E-9,
 7.18466403235734541950E-11,
};

/* 6.25 to infinity */
const static float CN[5] = {
-5.90592860534773254987E-1,
 6.29235242724368800674E-1,
-1.72858975380388136411E-1,
 1.64837047825189632310E-2,
-4.86827613020462700845E-4,
};
const static float CD[5] = {
/* 1.00000000000000000000E0,*/
-2.69820057197544900361E0,
 1.73270799045947845857E0,
-3.93708582281939493482E-1,
 3.44278924041233391079E-2,
-9.73655226040941223894E-4,
};

float xx, x, y;
int sign;
xx = xxx;
sign = 1;
if( xx < 0.0 ) {
	
      sign = -1;
      xx = -xx;
}

if( xx < 3.25 ) {
      x = xx*xx;
      y = xx * polevlf( x, AN, 9 )/polevlf( x, AD, 10 );
      return( sign * y );
}
x = 1.0/(xx*xx);
if( xx < 6.25 ){
       y = 1.0/xx + x * polevlf( x, BN, 10) / (p1evlf( x, BD, 10) * xx);
       return( sign * 0.5 * y );
}

if( xx > 1.0e9 )
	return( (sign * 0.5)/xx );

/* 6.25 to infinity */
y = 1.0/xx + x * polevlf( x, CN, 4) / (p1evlf( x, CD, 5) * xx);
return( sign * 0.5 * y );
}

/*DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *                phi
 *                 -
 *                | |
 *                |                   2
 * E(phi\m)  =    |    sqrt( 1 - m sin t ) dt
 *                |
 *              | |    
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ellief( const float phia, const float ma ) {
float phi, m, a, b, c, e, temp;
float lphi, t;
int d, mod;

phi = phia;
m = ma;
if( m == 0.0 )
	return( phi );
if( m == 1.0 )
	return( ceph_sinf(phi) );
lphi = phi;
if( lphi < 0.0 )
	lphi = -lphi;
a = 1.0;
b = 1.0 - m;
b = ceph_sqrtf(b);
c = ceph_sqrtf(m);
d = 1;
e = 0.0;
t = ceph_tanf( lphi );
mod = (lphi +1.5707963267948966192 )/3.141592653589793238;

while( fabsf(c/a) > 5.9604644775390625E-8)
	{
	temp = b/a;
	lphi = lphi + ceph_atanf(t*temp) + mod * 3.141592653589793238;
	mod = (lphi + 1.5707963267948966192)/3.141592653589793238;
	t = t * ( 1.0 + temp )/( 1.0 - temp * t * t );
	c = 0.5 * ( a - b );
	temp = ceph_sqrtf( a * b );
	a = 0.5 * ( a + b );
	b = temp;
	d += d;
	e += c * ceph_sinf(lphi);
	}

b = 1.0 - m;
temp = ellpef(b)/ellpkf(b);
temp *= (ceph_atanf(t) + mod * 3.141592653589793238)/(d * a);
temp += e;
if( phi < 0.0 )
	temp = -temp;
return( temp );
}

/*
    SYNOPSIS:
 *
 * float phi, m, y, ellikf();
 *
 * y = ellikf( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi\m)  =    |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ellikf(const  float phia, const float ma ){
float phi, m, a, b, c, temp;
float t;
int d, mod, sign;
phi = phia;
m = ma;
if( m == 0.0 )
	return( phi );
if( phi < 0.0 ){
	
	phi = -phi;
	sign = -1;
}
else{
	sign = 0;
}
a = 1.0;
b = 1.0 - m;
if( b == 0.0 )
	return(ceph_logf(ceph_tanf(0.5*(1.5707963267948966192 + phi))));
b = ceph_sqrtf(b);
c = ceph_sqrtf(m);
d = 1;
t = ceph_tanf( phi );
mod = (phi + 1.5707963267948966192)/3.141592653589793238;
while(fabsf(c/a) > 5.9604644775390625E-8) {
	temp = b/a;
	phi = phi + ceph_atanf(t*temp) + mod * PIF;
	mod = (phi + 1.5707963267948966192)/3.141592653589793238;
	t = t * ( 1.0 + temp )/( 1.0 - temp * t * t );
	c = ( a - b )/2.0;
	temp = ceph_sqrtf( a * b );
	a = ( a + b )/2.0;
	b = temp;
	d += d;
	}

temp = (ceph_atanf(t) + mod * 3.141592653589793238)/(d * a);
if( sign < 0 )
	temp = -temp;
return( temp );
}

/*
 SYNOPSIS:
 *
 * float m1, y, ellpef();
 *
 * y = ellpef( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *            pi/2
 *             -
 *            | |                 2
 * E(m)  =    |    sqrt( 1 - m sin t ) dt
 *          | |    
 *           -
 *            0
 *
 * Where m = 1 - m1, using the approximation
 *
 *      P(x)  -  x log x Q(x).
 *
 * Though there are no singularities, the argument m1 is used
 * rather than m for compatibility with ellpk().
 *
 * E(1) = 1; E(0) = pi/2
*/


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ellpef( const float xx){
const static float P[] = {
  1.53552577301013293365E-4,
  2.50888492163602060990E-3,
  8.68786816565889628429E-3,
  1.07350949056076193403E-2,
  7.77395492516787092951E-3,
  7.58395289413514708519E-3,
  1.15688436810574127319E-2,
  2.18317996015557253103E-2,
  5.68051945617860553470E-2,
  4.43147180560990850618E-1,
  1.00000000000000000299E0
};
const static float Q[] = {
  3.27954898576485872656E-5,
  1.00962792679356715133E-3,
  6.50609489976927491433E-3,
  1.68862163993311317300E-2,
  2.61769742454493659583E-2,
  3.34833904888224918614E-2,
  4.27180926518931511717E-2,
  5.85936634471101055642E-2,
  9.37499997197644278445E-2,
  2.49999999999888314361E-1
};

float x;
x = xx;
if( (x <= 0.0) || (x > 1.0) ){
	  if( x == 0.0 )
		return( 1.0 );
          return( 0.0 );
}
return( polevlf(x,P,10) - ceph_logf(x) * (x * polevlf(x,Q,9)) );
}

/*
 	Jacobian Elliptic Functions
 *
 *
 *
 * SYNOPSIS:
 *
 * float u, m, sn, cn, dn, phi;
 * int ellpj();
 *
 * ellpj( u, m, _&sn, _&cn, _&dn, _&phi );
 *
 *
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
int ellpjf( float uu, float mm,
            float * __restrict sn,
	    float * __restrict cn,
	    float * __restrict dn,
	    float * __restrict ph ) {
 __ATTR_ALIGN__(32) float a[10];
 __ATTR_ALIGN__(32) float c[10];
float u, m, ai, b, phi, t, twon;
int i;
u = uu;
m = mm;
/* Check for special cases */
if( m < 0.0 || m > 1.0 ){
     return(-1);
}
if( m < 1.0e-5 ) {
	t = ceph_sinf(u);
	b = ceph_cosf(u);
	ai = 0.25 * m * (u - t*b);
	*sn = t - ai*b;
	*cn = b + ai*t;
	*ph = u - ai;
	*dn = 1.0 - 0.5*m*t*t;
	return(0);
}
if( m >= 0.99999 ) {
	ai = 0.25 * (1.0-m);
	b = ceph_coshf(u);
	t = ceph_tanhf(u);
	phi = 1.0/b;
	twon = b * ceph_sinhf(u);
	*sn = t + ai * (twon - u)/(b*b);
	*ph = 2.0*ceph_atanf(ceph_expf(u)) - 1.5707963267948966192 + ai*(twon - u)/b;
	ai *= t * phi;
	*cn = phi - ai * (twon - u);
	*dn = phi + ai * (twon + u);
	return(0);
}

/*	A. G. M. scale		*/
a[0] = 1.0;
b = ceph_sqrtf(1.0 - m);
c[0] = ceph_sqrtf(m);
twon = 1.0;
i = 0;
while(ceph_fabsf((c[i]/a[i])) > 5.9604644775390625E-8) {
	
	if( i > 8 ) 
	     break;
	ai = a[i];
	++i;
	c[i] = 0.5 * ( ai - b );
	t = ceph_sqrtf( ai * b );
	a[i] = 0.5 * ( ai + b );
	b = t;
	twon += twon;
}


/* backward recurrence */
phi = twon * a[i] * u;
do{
	
	t = c[i] * ceph_sinf(phi) / a[i];
	b = phi;
	phi = 0.5 * (ceph_asinf(t) + phi);
}
while( --i );
t = ceph_sinf(phi);
*sn = t;
*cn = ceph_cosf(phi);
/* Thanks to Hartmut Henkel for reporting a bug here:  */
*dn = ceph_sqrtf(1.0 - m * t * t);
*ph = phi;
return(0);
}

/*
    SYNOPSIS:
 *
 * float m1, y, ellpkf();
 *
 * y = ellpkf( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ellpkf(const float xx){
 __ATTR_ALIGN__(32) const float P[] =
{
 1.37982864606273237150E-4,
 2.28025724005875567385E-3,
 7.97404013220415179367E-3,
 9.85821379021226008714E-3,
 6.87489687449949877925E-3,
 6.18901033637687613229E-3,
 8.79078273952743772254E-3,
 1.49380448916805252718E-2,
 3.08851465246711995998E-2,
 9.65735902811690126535E-2,
 1.38629436111989062502E0
};

 __ATTR_ALIGN__(32) const float Q[] =
{
 2.94078955048598507511E-5,
 9.14184723865917226571E-4,
 5.94058303753167793257E-3,
 1.54850516649762399335E-2,
 2.39089602715924892727E-2,
 3.01204715227604046988E-2,
 3.73774314173823228969E-2,
 4.88280347570998239232E-2,
 7.03124996963957469739E-2,
 1.24999999999870820058E-1,
 4.99999999999999999821E-1
};
float x;
x = xx;
if( (x < 0.0) || (x > 1.0) ) {
	
	return( 0.0 );
}

if( x > 5.9604644775390625E-8){
    return( polevlf(x,P,10) - ceph_logf(x) * polevlf(x,Q,10) );
}
else
    {
	if( x == 0.0 ) {
	   return(3.4028234663852885981170418348451692544e38);
       }
	else{
	     return(1.3862943611198906188E0 - 0.5 * ceph_logf(x));
	}
    }
}

/*
   SYNOPSIS:
 *
 * float x, y, exp2f();
 *
 * y = exp2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 2 raised to the x power.
 *
 * Range reduction is accomplished by separating the argument
 * into an integer k and fraction f such that
 *     x    k  f
 *    2  = 2  2.
 *
 * A polynomial approximates 2**x in the basic range [-0.5, 0.5].
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_exp2f( const float xx ){
__ATTR_ALIGN__(32) const float P[] = {
 1.535336188319500E-004,
 1.339887440266574E-003,
 9.618437357674640E-003,
 5.550332471162809E-002,
 2.402264791363012E-001,
 6.931472028550421E-001
};
float x, px;
int i0;
x = xx;
if( x > 127.0){
	return(3.4028234663852885981170418348451692544e38);
}
if( x < -127.0) {
	return(0.0);
}
/* The following is necessary because range reduction blows up: */
if( x == 0 )
	return(1.0);

/* separate into integer and fractional parts */
px = ceph_floorf(x);
i0 = px;
x = x - px;
if( x > 0.5 ) {
     i0 += 1;
     x -= 1.0;
}

/* rational approximation
 * exp2(x) = 1.0 +  xP(x)
 */
px = 1.0 + x * polevlf( x, P, 5 );
/* scale by power of 2 */
px = ceph_ldexpf( px, i0 );
return(px);
}

/*
   SYNOPSIS:
 *
 * float x, y, exp10f();
 *
 * y = exp10f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 10 raised to the x power.
 *
 * Range reduction is accomplished by expressing the argument
 * as 10**x = 2**n 10**f, with |f| < 0.5 log10(2).
 * A polynomial approximates 10**f.
 *
*/
__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_exp10f(const float xx) {
 __ATTR_ALIGN__(32) const float P[] = {
 2.063216740311022E-001,
 5.420251702225484E-001,
 1.171292686296281E+000,
 2.034649854009453E+000,
 2.650948748208892E+000,
 2.302585167056758E+000
};
float x, px, qx;
short n;
x = xx;
if( x > 38.230809449325611792) {
	return(3.4028234663852885981170418348451692544e38);
}
if( x < -38.230809449325611792){ /* Would like to use MINLOG but can't */
	return(0.0);
}
/* The following is necessary because range reduction blows up: */
if( x == 0 )
	return(1.0);

/* Express 10**x = 10**g 2**n
 *   = 10**g 10**( n log10(2) )
 *   = 10**( g + n log10(2) )
 */
px = x * 3.32192809488736234787e0;
qx = ceph_floorf( px + 0.5 );
n = qx;
x -= qx * 3.00781250000000000000E-1;
x -= qx * 2.48745663981195213739E-4;

/* rational approximation for exponential
 * of the fractional part:
 * 10**x - 1  =  2x P(x**2)/( Q(x**2) - P(x**2) )
 */
px = 1.0 + x * polevlf( x, P, 5 );
/* multiply by power of 2 */
x = ceph_ldexpf( px, n );
return(x);
}

/*
 SYNOPSIS:
 *
 * int n;
 * float x, y, expnf();
 *
 * y = expnf( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float expnf( const int n, const float xx ) {
float x, ans, r, t, yk, xk;
float pk, pkm1, pkm2, qk, qkm1, qkm2;
float psi, z;
int i, k;
constexpr float big = 16777216.0;
x = xx;
if( n < 0 )
	goto domerr;
if( x < 0 ) {
	
domerr:
	return(3.4028234663852885981170418348451692544e38);
}
if( x > 88.72283905206835)
	return( 0.0 );

if( x == 0.0 ) {
      if( n < 2 ) {
	  return(3.4028234663852885981170418348451692544e38);
	}
	else
		return( 1.0/(n-1.0) );
	}

if( n == 0 )
	return( ceph_expf(-x)/x );

/*							expn.c	*/
/*		Expansion for large n		*/

if( n > 5000 ) {
	
	xk = x + n;
	yk = 1.0 / (xk * xk);
	t = n;
	ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
	ans = yk * (ans + t * (t  -  2.0 * x));
	ans = yk * (ans + t);
	ans = (ans + 1.0) * ceph_expf( -x ) / xk;
	goto done;
}

if( x > 1.0 )
	goto cfrac;
 
/*							expn.c	*/

/*		Power series expansion		*/

psi = -0.57721566490153286060 - ceph_logf(x);
for( i=1; i<n; i++ )
	psi = psi + 1.0/i;

z = -x;
xk = 0.0;
yk = 1.0;
pk = 1.0 - n;
if( n == 1 )
	ans = 0.0;
else
	ans = 1.0/pk;
do{
	xk += 1.0;
	yk *= z/xk;
	pk += 1.0;
	if( pk != 0.0 ) {
	      ans += yk/pk;
	}
	if( ans != 0.0 )
		t = ceph_fabsf(yk/ans);
	else
		t = 1.0;
}
while(t > 5.9604644775390625E-8);
k = xk;
t = n;
r = n - 1;
ans = (ceph_powf(z, r) * psi / ceph_gammaf(t)) - ans;
goto done;

/*							expn.c	*/
/*		continued fraction		*/
cfrac:
k = 1;
pkm2 = 1.0;
qkm2 = x;
pkm1 = 1.0;
qkm1 = x + n;
ans = pkm1/qkm1;

do{
	
    k += 1;
    if( k & 1 ) {
	yk = 1.0;
	xk = n + (k-1)/2;
    }
    else
	{
	yk = x;
	xk = k/2;
    }
     pk = pkm1 * yk  +  pkm2 * xk;
     qk = qkm1 * yk  +  qkm2 * xk;
     if( qk != 0 ) {
	 r = pk/qk;
	 t = ceph_fabsf( (ans - r)/r );
	 ans = r;
     }
     else
	t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
if( ceph_fabsf(pk) > big ) {
      pkm2 *= 5.9604644775390625E-8;
      pkm1 *= 5.9604644775390625E-8;
      qkm2 *= 5.9604644775390625E-8;
      qkm1 *= 5.9604644775390625E-8;
    }
}
while( t > MACHEPF );
ans *= expf( -x );
done:
return( ans );
}

/*
   *
 * SYNOPSIS:
 *
 * double x, y, expx2f();
 *
 * y = expx2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes y = exp(x*x) while suppressing error amplification
 * that would ordinarily arise from the inexactness of the argument x*x.
 * 
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float expx2f (const float x) {
  float u, u1, m;
  if (x < 0)
    x = -x;
  /* Represent x as an exact multiple of 1/32 plus a residual.  */
  m = .03125f * ceph_floorf(32.0f * x + 0.5f);
  x -= m;
  /* x**2 = m**2 + 2mf + f**2 */
  u = m * m;
  u1 = 2 * m * x  +  x * x;
  if ((u+u1) > 88.72283905206835)
    return (3.4028234663852885981170418348451692544e38);

  /* u is exact, u1 is small.  */
  u = ceph_expf(u) * ceph_expf(u1);
  return(u);
}



/** SYNOPSIS:
 *
 * int N;
 * float x, y, coef[N+1], polevlf[];
 *
 * y = polevlf( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float polevlf( float xx, float * __restrict coef, int N ) {
float ans, x;
float *p;
int i;
x = xx;
p = coef;
ans = *p++;
/*
for( i=0; i<N; i++ )
	ans = ans * x  +  *p++;
*/

i = N;
do{
	ans = ans * x  +  *p++;
}
while( --i );
return( ans );
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float p1evlf( float xx,
              float * __restrict coef,
	      int N ) {
float ans, x;
float *p;
int i;
x = xx;
p = coef;
ans = x + *p++;
i = N-1;
do{
	ans = ans * x  + *p++;
}
while( --i );

return( ans );
}



#endif /*__GMS_CEPHES_H__*/
