
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

/*
  SYNOPSIS:
 *
 * float x, S, C;
 * void fresnlf();
 *
 * fresnlf( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
void fresnlf(const  float xxa,
             float * __restrict ssa,
	     float *__restrict cca ){
/* S(x) for small x */
__ATTR_ALIGN__(32) const float sn[7] = {
 1.647629463788700E-009,
-1.522754752581096E-007,
 8.424748808502400E-006,
-3.120693124703272E-004,
 7.244727626597022E-003,
-9.228055941124598E-002,
 5.235987735681432E-001
};

/* C(x) for small x */
__ATTR_ALIGN__(32) const float cn[7] = {
 1.416802502367354E-008,
-1.157231412229871E-006,
 5.387223446683264E-005,
-1.604381798862293E-003,
 2.818489036795073E-002,
-2.467398198317899E-001,
 9.999999760004487E-001
};


/* Auxiliary function f(x) */
__ATTR_ALIGN__(32) const float fn[8] = {
-1.903009855649792E+012,
 1.355942388050252E+011,
-4.158143148511033E+009,
 7.343848463587323E+007,
-8.732356681548485E+005,
 8.560515466275470E+003,
-1.032877601091159E+002,
 2.999401847870011E+000
};

/* Auxiliary function g(x) */
__ATTR_ALIGN__(32) const float gn[8] = {
-1.860843997624650E+011,
 1.278350673393208E+010,
-3.779387713202229E+008,
 6.492611570598858E+006,
-7.787789623358162E+004,
 8.602931494734327E+002,
-1.493439396592284E+001,
 9.999841934744914E-001
};
float f, g, cc, ss, c, s, t, u, x, x2;
/*debug double t1;*/
x = xxa;
x = ceph_fabsf(x);
x2 = x * x;
if( x2 < 2.5625 ) {
	t = x2 * x2;
	ss = x * x2 * polevlf( t, sn, 6);
	cc = x * polevlf( t, cn, 6);
	goto done;
}

if( x > 36974.0 ) {
	
	cc = 0.5;
	ss = 0.5;
	goto done;
}


/*		Asymptotic power series auxiliary functions
 *		for large argument
 */
	x2 = x * x;
	t = 3.141592653589793238 * x2;
	u = 1.0/(t * t);
	t = 1.0/t;
	f = 1.0 - u * polevlf( u, fn, 7);
	g = t * polevlf( u, gn, 7);

	t = 1.5707963267948966192 * x2;
	c = ceph_cosf(t);
	s = ceph_sinf(t);
	t = 3.141592653589793238 * x;
	cc = 0.5  +  (f * s  -  g * c)/t;
	ss = 0.5  -  (f * c  +  g * s)/t;

done:
if( xxa < 0.0 ) {
	cc = -cc;
	ss = -ss;
}

*cca = cc;
*ssa = ss;
return 0;
}

/*
    SYNOPSIS:
 *
 * float x, y, gammaf();
 * extern int sgngamf;
 *
 * y = gammaf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngamf.
 * This same variable is also filled in by the logarithmic
 * gamma function lgam().
 *
 * Arguments between 0 and 10 are reduced by recurrence and the
 * function is approximated by a polynomial function covering
 * the interval (2,3).  Large arguments are handled by Stirling's
 * formula. Negative arguments are made positive using
 * a reflection formula.  
 *
*/

/* Gamma function computed by Stirling's formula,
 * sqrt(2 pi) x^(x-.5) exp(-x) (1 + 1/x P(1/x))
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float stirf(const float xx) {
const float STIR[] = {
-2.705194986674176E-003,
 3.473255786154910E-003,
 8.333331788340907E-002
};
float x, y, w, v;
x = xx;
w = 1.0/x;
w = 1.0 + w * polevlf( w, STIR, 2 );
y = ceph_expf( -x );
if( x > 26.77) {
	 /* Avoid overflow in pow() */
	v = ceph_powf( x, 0.5 * x - 0.25 );
	y *= v;
	y *= v;
}
else{
    y = ceph_powf( x, x - 0.5 ) * y;
}
y = 2.50662827463100050242 * y * w;
return( y );
}


__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float ceph_gammaf(const float xx ) {
/* gamma(x+2), 0 < x < 1 */
__ATTR_ALIGN__(32) const float P[] = {
 1.536830450601906E-003,
 5.397581592950993E-003,
 4.130370201859976E-003,
 7.232307985516519E-002,
 8.203960091619193E-002,
 4.117857447645796E-001,
 4.227867745131584E-001,
 9.999999822945073E-001
};
float p, q, x, z, nz;
int i, direction, negative;
x = xx;
sgngamf = 1;
negative = 0;
nz = 0.0;
if( x < 0.0 ) {
	negative = 1;
	q = -x;
	p = ceph_floorf(q);
	if( p == q )
		goto goverf;
	i = p;
	if((i & 1) == 0 )
		sgngamf = -1;
	nz = q - p;
	if( nz > 0.5 ) {
	        p += 1.0;
		nz = q - p;
	}
	nz = q * ceph_sinf( 3.141592653589793238 * nz );
	if( nz == 0.0 ) {
		
goverf:
	     return(sgngamf * 3.4028234663852885981170418348451692544e38);
        }
	if( nz < 0 )
		nz = -nz;
	x = q;
	}
if( x >= 10.0 ) {
    z = stirf(x);
}
if( x < 2.0 )
	direction = 1;
else
	direction = 0;
z = 1.0;
while( x >= 3.0 ) {
        x -= 1.0;
	z *= x;
}
/*
while( x < 0.0 )
	{
	if( x > -1.E-4 )
		goto small;
	z *=x;
	x += 1.0;
	}
*/
while( x < 2.0 ) {
	
	if( x < 1.e-4 )
		goto small;
	z *=x;
	x += 1.0;
}

if( direction )
	z = 1.0/z;

if( x == 2.0 )
	return(z);

x -= 2.0;
p = z * polevlf( x, P, 7 );

gdone:

if( negative ) {
    p = sgngamf * 3.141592653589793238/(nz * p );
}
return(p);

small:
if( x == 0.0 ){
	return(3.4028234663852885981170418348451692544e38);
}
else{
	p = z / ((1.0 + 0.5772156649015329 * x) * x);
	goto gdone;
    }
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float lgamf(const float xx ){
/* log gamma(x+2), -.5 < x < .5 */
__ATTR_ALIGN__(32) const float B[] = {
 6.055172732649237E-004,
-1.311620815545743E-003,
 2.863437556468661E-003,
-7.366775108654962E-003,
 2.058355474821512E-002,
-6.735323259371034E-002,
 3.224669577325661E-001,
 4.227843421859038E-001
};

/* log gamma(x+1), -.25 < x < .25 */
__ATTR_ALIGN__(32) const float C[] = {
 1.369488127325832E-001,
-1.590086327657347E-001,
 1.692415923504637E-001,
-2.067882815621965E-001,
 2.705806208275915E-001,
-4.006931650563372E-001,
 8.224670749082976E-001,
-5.772156501719101E-001
};
float p, q, w, z, x;
float nx, tx;
int i, direction;
sgngamf = 1;
x = xx;
if( x < 0.0 ) {
	q = -x;
	w = ceph_lgamf(q); /* note this modifies sgngam! */
	p = ceph_floorf(q);
	if( p == q )
		goto loverf;
	i = p;
	if( (i & 1) == 0 )
		sgngamf = -1;
	else
		sgngamf = 1;
	z = q - p;
	if( z > 0.5 ) {
		p += 1.0;
		z = p - q;
	}
	z = q * ceph_sinf(3.141592653589793238 * z);
	if(z == 0.0)
		goto loverf;
	z = -ceph_logf( 0.318309886183790671538*z ) - w;
	return( z );
	}

if( x < 6.5 ) {
	direction = 0;
	z = 1.0;
	tx = x;
	nx = 0.0;
	if(x >= 1.5) {
	       while( tx > 2.5 ){
			nx -= 1.0;
			tx = x + nx;
			z *=tx;
		}
		x += nx - 2.0;
iv1r5:
		p = x * polevlf( x, B, 7 );
		goto cont;
	}
	if( x >= 1.25 ){
		z *= x;
		x -= 1.0; /* x + 1 - 2 */
		direction = 1;
		goto iv1r5;
	}
	if( x >= 0.75 ) {
	        x -= 1.0;
		p = x * polevlf( x, C, 7 );
		q = 0.0;
		goto contz;
	}
	while( tx < 1.5 ){
	       if( tx == 0.0 )
			goto loverf;
		z *=tx;
		nx += 1.0;
		tx = x + nx;
	}
	direction = 1;
	x += nx - 2.0;
	p = x * polevlf( x, B, 7 );

cont:
	if( z < 0.0 ){
		sgngamf = -1;
		z = -z;
	}
	else
	    {
		sgngamf = 1;
	}
	q = ceph_logf(z);
	if( direction )
		q = -q;
contz:
	return( p + q );
	}
if( x > 2.035093e36){
	
loverf:
	return(sgngamf * 3.4028234663852885981170418348451692544e38);
}

/* Note, though an asymptotic formula could be used for x >= 3,
 * there is cancellation error in the following if x < 6.5.  */
q =  0.91893853320467274178 - x;
q += ( x - 0.5 ) * ceph_logf(x);
if( x <= 1.0e4 ){
	z = 1.0/x;
	p = z * z;
	q += ((    6.789774945028216E-004 * p
		 - 2.769887652139868E-003 ) * p
		+  8.333316229807355E-002 ) * z;
}
return( q );
}

/*
	Gauss hypergeometric function   F
 *	                               2 1
 *
 *
 * SYNOPSIS:
 *
 * float a, b, c, x, y, hyp2f1f();
 *
 * y = hyp2f1f( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 *  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
 *                           2 1
 *
 *           inf.
 *            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
 *   =  1 +   >   -----------------------------  x   .
 *            -         c(c+1)...(c+k) (k+1)!
 *          k = 0
 *
 *  Cases addressed are
 *	Tests and escapes for negative integer a, b, or c
 *	Linear transformation if c - a or c - b negative integer
 *	Special case c = a or c = b
 *	Linear transformation for  x near +1
 *	Transformation for x < -0.5
 *	Psi function expansion if x > 0.5 and c - a - b integer
 *      Conditionally, a recurrence on c to make c-a-b > 0
 *
*/

#if !defined(ceph_roundf)
    #define ceph_roundf(x) (ceph_floorf((x)+(float )0.5))
#endif


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hyp2f1f(const float aa,
              const float bb,
	      const float cc,
	      const float xx ){
float a, b, c, x;
float d, d1, d2, e;
float p, q, r, s, y, ax;
float ia, ib, ic, id, err;
int flag, i, aid;
a = aa;
b = bb;
c = cc;
x = xx;
err = 0.0;
ax = ceph_fabsf(x);
s = 1.0 - x;
flag = 0;
ia = ceph_roundf(a); /* nearest integer to a */
ib = ceph_roundf(b);
if( a <= 0 )
	{
	if( ceph_fabsf(a-ia) < 1.0e-5)		/* a is a negative integer */
		flag |= 1;
	}

if( b <= 0 )
	{
	if( ceph_fabsf(b-ib) < 1.0e-5)		/* b is a negative integer */
		flag |= 2;
	}

if( ax < 1.0 )
	{
	if( ceph_fabsf(b-c) < 1.0e-5)		/* b = c */
		{
		y = ceph_powf( s, -a );	/* s to the -a power */
		goto hypdon;
		}
	if( ceph_fabsf(a-c) < 1.0e-5)		/* a = c */
		{
		y = ceph_powf( s, -b );	/* s to the -b power */
		goto hypdon;
		}
	}

if( c <= 0.0 )
	{
	ic = ceph_roundf(c); 	/* nearest integer to c */
	if( ceph_fabsf(c-ic) < 1.0e-5)		/* c is a negative integer */
		{
		/* check if termination before explosion */
		if( (flag & 1) && (ia > ic) )
			goto hypok;
		if( (flag & 2) && (ib > ic) )
			goto hypok;
		goto hypdiv;
		}
	}

if( flag )			/* function is a polynomial */
	goto hypok;

if( ax > 1.0 )			/* series diverges	*/
	goto hypdiv;

p = c - a;
ia = ceph_roundf(p);
if( (ia <= 0.0) && (ceph_fabsf(p-ia) < 1.0e-5) )	/* negative int c - a */
	flag |= 4;

r = c - b;
ib = ceph_roundf(r); /* nearest integer to r */
if( (ib <= 0.0) && (ceph_fabsf(r-ib) < 1.0e-5) )	/* negative int c - b */
	flag |= 8;

d = c - a - b;
id = ceph_roundf(d); /* nearest integer to d */
q = ceph_fabsf(d-id);

if( ceph_fabsf(ax-1.0) < 1.0e-5)			/* |x| == 1.0	*/
	{
	if( x > 0.0 )
		{
		if( flag & 12 ) /* negative int c-a or c-b */
			{
			if( d >= 0.0 )
				goto hypf;
			else
				goto hypdiv;
			}
		if( d <= 0.0 )
			goto hypdiv;
		y = ceph_gammaf(c)*ceph_gammaf(d)/(ceph_gammaf(p)*ceph_gammaf(r));
		goto hypdon;
		}

	if( d <= -1.0 )
		goto hypdiv;
	}

/* Conditionally make d > 0 by recurrence on c
 * AMS55 #15.2.27
 */
if( d < 0.0 )
	{
/* Try the power series first */
	y = hyt2f1f( a, b, c, x, &err );
	if( err < 1.0e-5)
		goto hypdon;
/* Apply the recurrence if power series fails */
	err = 0.0;
	aid = 2 - id;
	e = c + aid;
	d2 = hyp2f1f(a,b,e,x);
	d1 = hyp2f1f(a,b,e+1.0,x);
	q = a + b + 1.0;
	for( i=0; i<aid; i++ )
		{
		r = e - 1.0;
		y = (e*(r-(2.0*e-q)*x)*d2 + (e-a)*(e-b)*x*d1)/(e*r*s);
		e = r;
		d1 = d2;
		d2 = y;
		}
	goto hypdon;
	}


if( flag & 12 )
	goto hypf; /* negative integer c-a or c-b */

hypok:
y = hyt2f1f( a, b, c, x, &err );

hypdon:
if( err > 1.0e-5)
	{
//	mtherr( "hyp2f1", PLOSS );
/*	printf( "Estimated err = %.2e\n", err );*/
	}
return(y);

/* The transformation for c-a or c-b negative integer
 * AMS55 #15.3.3
 */
hypf:
y = ceph_powf( s, d ) * hys2f1f( c-a, c-b, c, x, &err );
goto hypdon;

/* The alarm exit */
hypdiv:
return(3.4028234663852885981170418348451692544e38);
}

/* Apply transformations for |x| near 1
 * then call the power series
 */


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hyt2f1f( const float aa,
               const float bb,
	       const float cc,
	       const float xx,
	       float *loss ) {
float a, b, c, x;
float p, q, r, s, t, y, d, err, err1;
float ax, id, d1, d2, e, y1;
int i, aid;

a = aa;
b = bb;
c = cc;
x = xx;
err = 0.0;
s = 1.0 - x;
if( x < -0.5 )
	{
	if( b > a )
		y = ceph_powf( s, -a ) * hys2f1f( a, c-b, c, -x/s, &err );

	else
		y = ceph_powf( s, -b ) * hys2f1f( c-a, b, c, -x/s, &err );

	goto done;
	}



d = c - a - b;
id = ceph_roundf(d);	/* nearest integer to d */

if( x > 0.8 )
{

if( ceph_fabsf(d-id) > 1.0e-5) /* test for integer c-a-b */
	{
/* Try the power series first */
	y = hys2f1f( a, b, c, x, &err );
	if( err < 1.0e-5)
		goto done;
/* If power series fails, then apply AMS55 #15.3.6 */
	q = hys2f1f( a, b, 1.0-d, s, &err );	
	q *= ceph_gammaf(d) /(ceph_gammaf(c-a) * ceph_gammaf(c-b));
	r = ceph_powf(s,d) * hys2f1f( c-a, c-b, d+1.0, s, &err1 );
	r *= ceph_gammaf(-d)/(ceph_gammaf(a) * ceph_gammaf(b));
	y = q + r;

	q = ceph_fabsf(q); /* estimate cancellation error */
	r = ceph_fabsf(r);
	if( q > r )
		r = q;
	err += err1 + (5.9604644775390625E-8*r)/y;

	y *= ceph_gammaf(c);
	goto done;
	}	
else
	{
/* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12 */
	if( id >= 0.0 )
		{
		e = d;
		d1 = d;
		d2 = 0.0;
		aid = id;
		}
	else
		{
		e = -d;
		d1 = 0.0;
		d2 = d;
		aid = -id;
		}

	ax = ceph_logf(s);

	/* sum for t = 0 */
	y = psif(1.0) + psif(1.0+e) - psif(a+d1) - psif(b+d1) - ax;
	y /= ceph_gammaf(e+1.0);

	p = (a+d1) * (b+d1) * s / ceph_gammaf(e+2.0);	/* Poch for t=1 */
	t = 1.0;
	do
		{
		r = psif(1.0+t) + psif(1.0+t+e) - psif(a+t+d1)
			- psif(b+t+d1) - ax;
		q = p * r;
		y += q;
		p *= s * (a+t+d1) / (t+1.0);
		p *= (b+t+d1) / (t+1.0+e);
		t += 1.0;
		}
	while( fabsf(q/y) > EPS );


	if( id == 0.0 )
		{
		y *= ceph_gammaf(c)/(ceph_gammaf(a)*ceph_gammaf(b));
		goto psidon;
		}

	y1 = 1.0;

	if( aid == 1 )
		goto nosum;

	t = 0.0;
	p = 1.0;
	for( i=1; i<aid; i++ )
		{
		r = 1.0-e+t;
		p *= s * (a+t+d2) * (b+t+d2) / r;
		t += 1.0;
		p /= t;
		y1 += p;
		}


nosum:
	p = ceph_gammaf(c);
	y1 *= ceph_gammaf(e) * p / (ceph_gammaf(a+d1) * ceph_gammaf(b+d1));
	y *= p / (ceph_gammaf(a+d2) * ceph_gammaf(b+d2));
	if( (aid & 1) != 0 )
		y = -y;

	q = powf( s, id );	/* s to the id power */
	if( id > 0.0 )
		y *= q;
	else
		y1 *= q;

	y += y1;
psidon:
	goto done;
	}
}


/* Use defining power series if no special cases */
y = hys2f1f( a, b, c, x, &err );

done:
*loss = err;
return(y);
}

/* Defining power series expansion of Gauss hypergeometric function */

#if 1

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hys2f1f( const float aa,
               const float bb,
	       const float cc,
	       const float xx,
	       float *loss ) {
int i;
float a, b, c, x;
float f, g, h, k, m, s, u, umax;


a = aa;
b = bb;
c = cc;
x = xx;
i = 0;
umax = 0.0;
f = a;
g = b;
h = c;
k = 0.0;
s = 1.0;
u = 1.0;

do
	{
	if( ceph_fabsf(h) < 1.0e-5)
		return(3.4028234663852885981170418348451692544e38);
	m = k + 1.0;
	u = u * ((f+k) * (g+k) * x / ((h+k) * m));
	s += u;
	k = ceph_fabsf(u);  /* remember largest term summed */
	if( k > umax )
		umax = k;
	k = m;
	if( ++i > 10000 ) /* should never happen */
		{
		*loss = 1.0;
		return(s);
		}
	}
while( ceph_fabsf(u/s) > 5.9604644775390625E-8);

/* return estimated relative error */
*loss = (5.9604644775390625E-8*umax)/ceph_fabsf(s) + (5.9604644775390625E-8*i);

return(s);
}


#else /* 0 */




__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hys2f1f( const float aa,
                      const float bb,
		      const float cc,
		      const float xx,
		      float *loss ) {
int i;
double a, b, c, x;
double f, g, h, k, m, s, u, umax;

a = aa;
b = bb;
c = cc;
x = xx;
i = 0;
umax = 0.0;
f = a;
g = b;
h = c;
k = 0.0;
s = 1.0;
u = 1.0;

do
	{
	if( ceph_fabsf(h) < 1.0e-5)
		{
		*loss = 1.0;
		return(3.4028234663852885981170418348451692544e38);
		}
	m = k + 1.0;
	u = u * ((f+k) * (g+k) * x / ((h+k) * m));
	s += u;
	k = ceph_fabsf(u);  /* remember largest term summed */
	if( k > umax )
		umax = k;
	k = m;
	if( ++i > 10000 ) /* should never happen */
		{
		*loss = 1.0;
		return(s);
		}
	}
while( ceph_fabsf(u/s) > 5.9604644775390625E-8);

/* return estimated relative error */
*loss = (5.9604644775390625E-8*umax)/ceph_fabsf(s) + (5.9604644775390625E-8*i);

return(s);
}
#endif

/*
   	Confluent hypergeometric function
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, b, x, y, hypergf();
 *
 * y = hypergf( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the confluent hypergeometric function
 *
 *                          1           2
 *                       a x    a(a+1) x
 *   F ( a,b;x )  =  1 + ---- + --------- + ...
 *  1 1                  b 1!   b(b+1) 2!
 *
 * Many higher transcendental functions are special cases of
 * this power series.
 *
 * As is evident from the formula, b must not be a negative
 * integer or zero unless a is an integer with 0 >= a > b.
 *
 * The routine attempts both a direct summation of the series
 * and an asymptotic expansion.  In each case error due to
 * roundoff, cancellation, and nonconvergence is estimated.
 * The result with smaller estimated error is returned.
*/

 __ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hypergf( float aa, float bb, float xx ){
float a, b, x, asum, psum, acanc, pcanc, temp;
a = aa;
b = bb;
x = xx;
/* See if a Kummer transformation will help */
temp = b - a;
if( ceph_fabsf(temp) < 0.001 * ceph_fabsf(a) )
	return( ceph_expf(x) * hypergf( temp, b, -x )  );

psum = hy1f1pf( a, b, x, &pcanc );
if( pcanc < 1.0e-6 )
	goto done;


/* try asymptotic series */

asum = hy1f1af( a, b, x, &acanc );


/* Pick the result with less estimated error */

if( acanc < pcanc )
	{
	pcanc = acanc;
	psum = asum;
	}

done:
if( pcanc > 1.0e-3 )
    ;

return( psum );
}

/* Power series summation for confluent hypergeometric function		*/



__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hy1f1pf( const float aa,
               const float bb,
	       const float xx,
	       float *err ) {
float a, b, x, n, a0, sum, t, u, temp;
float an, bn, maxt, pcanc;
a = aa;
b = bb;
x = xx;
/* set up for power series summation */
an = a;
bn = b;
a0 = 1.0;
sum = 1.0;
n = 1.0;
t = 1.0;
maxt = 0.0;

while( t > 5.9604644775390625E-8)
	{
	if( bn == 0 )			/* check bn first since if both	*/
		{
	
		return(3.4028234663852885981170418348451692544e38 );	/* an and bn are zero it is	*/
		}
	if( an == 0 )			/* a singularity		*/
		return( sum );
	if( n > 200 )
		goto pdone;
	u = x * ( an / (bn * n) );

	/* check for blowup */
	temp = ceph_fabsf(u);
	if( (temp > 1.0 ) && (maxt > (3.4028234663852885981170418348451692544e38/temp)) )
		{
		pcanc = 1.0;	/* estimate 100% error */
		goto blowup;
		}

	a0 *= u;
	sum += a0;
	t = ceph_fabsf(a0);
	if( t > maxt )
		maxt = t;
/*
	if( (maxt/fabsf(sum)) > 1.0e17 )
		{
		pcanc = 1.0;
		goto blowup;
		}
*/
	an += 1.0;
	bn += 1.0;
	n += 1.0;
	}

pdone:

/* estimate error due to roundoff and cancellation */
if( sum != 0.0 )
	maxt /= ceph_fabsf(sum);
maxt *= 5.9604644775390625E-8; 	/* this way avoids multiply overflow */
pcanc = ceph_fabsf( 5.9604644775390625E-8* n  +  maxt );

blowup:

*err = pcanc;

return( sum );
}

/*							hy1f1a()	*/
/* asymptotic formula for hypergeometric function:
 *
 *        (    -a                         
 *  --    ( |z|                           
 * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
 *        (  --                           
 *        ( |  (b-a)                      
 *
 *
 *                                x    a-b                     )
 *                               e  |x|                        )
 *                             + -------- 2f0( b-a, 1-a, 1/x ) )
 *                                --                           )
 *                               |  (a)                        )
 */


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hy1f1af( const
                      const float aa,
		      const float bb,
		      const float xx,
		      float *err ) {
float a, b, x, h1, h2, t, u, temp, acanc, asum, err1, err2;

a = aa;
b = bb;
x = xx;
if( x == 0 )
	{
	acanc = 1.0;
	asum = 3.4028234663852885981170418348451692544e38;
	goto adone;
	}
temp = ceph_logf( ceph_fabsf(x) );
t = x + temp * (a-b);
u = -temp * a;

if( b > 0 )
	{
	temp = ceph_lgamf(b);
	t += temp;
	u += temp;
	}

h1 = hyp2f0f( a, a-b+1, -1.0/x, 1, &err1 );

temp = ceph_expf(u) / ceph_gammaf(b-a);
h1 *= temp;
err1 *= temp;

h2 = hyp2f0f( b-a, 1.0-a, 1.0/x, 2, &err2 );

if( a < 0 )
	temp = ceph_expf(t) / ceph_gammaf(a);
else
	temp = ceph_expf( t - ceph_lgamf(a) );

h2 *= temp;
err2 *= temp;

if( x < 0.0 )
	asum = h1;
else
	asum = h2;

acanc = ceph_fabsf(err1) + ceph_fabsf(err2);


if( b < 0 )
	{
	temp = ceph_gammaf(b);
	asum *= temp;
	acanc *= ceph_fabsf(temp);
	}


if( asum != 0.0 )
	acanc /= ceph_fabsf(asum);

acanc *= 30.0;	/* fudge factor, since error of asymptotic formula
		 * often seems this much larger than advertised */

adone:


*err = acanc;
return( asum );
}

__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float hyp2f0f(const float aa,
              const float bb,
	      const float xx,
	      const int type,
	      const float *err){
float a, b, x, a0, alast, t, tlast, maxt;
float n, an, bn, u, sum, temp;

a = aa;
b = bb;
x = xx;
an = a;
bn = b;
a0 = 1.0;
alast = 1.0;
sum = 0.0;
n = 1.0;
t = 1.0;
tlast = 1.0e9;
maxt = 0.0;

do
	{
	if( an == 0 )
		goto pdone;
	if( bn == 0 )
		goto pdone;

	u = an * (bn * x / n);

	/* check for blowup */
	temp = ceph_fabsf(u);
	if( (temp > 1.0 ) && (maxt > (3.4028234663852885981170418348451692544e38/temp)) )
		goto error;

	a0 *= u;
	t = ceph_fabsf(a0);

	/* terminating condition for asymptotic series */
	if( t > tlast )
		goto ndone;

	tlast = t;
	sum += alast;	/* the sum is one term behind */
	alast = a0;

	if( n > 200 )
		goto ndone;

	an += 1.0;
	bn += 1.0;
	n += 1.0;
	if( t > maxt )
		maxt = t;
	}
while( t > 5.9604644775390625E-8);


pdone:	/* series converged! */

/* estimate error due to roundoff and cancellation */
*err = ceph_fabsf(  5.9604644775390625E-8 * (n + maxt)  );

alast = a0;
goto done;

ndone:	/* series did not converge */

/* The following "Converging factors" are supposed to improve accuracy,
 * but do not actually seem to accomplish very much. */

n -= 1.0;
x = 1.0/x;

switch( type )	/* "type" given as subroutine argument */
{
case 1:
	alast *= ( 0.5 + (0.125 + 0.25*b - 0.5*a + 0.25*x - 0.25*n)/x );
	break;

case 2:
	alast *= 2.0/3.0 - b + 2.0*a + x - n;
	break;

default:
	;
}

/* estimate error due to roundoff, cancellation, and nonconvergence */
*err = 5.9604644775390625E-8* (n + maxt)  +  ceph_fabsf( a0 );


done:
sum += alast;
return( sum );

/* series blew up: */
error:
*err = 3.4028234663852885981170418348451692544e38;
return( sum );
}

/*
    	Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i0();
 *
 * y = i0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float i0f( float x ) {
__ATTR_ALIGN__(64) const float A[] =
{
-1.30002500998624804212E-8f,
 6.04699502254191894932E-8f,
-2.67079385394061173391E-7f,
 1.11738753912010371815E-6f,
-4.41673835845875056359E-6f,
 1.64484480707288970893E-5f,
-5.75419501008210370398E-5f,
 1.88502885095841655729E-4f,
-5.76375574538582365885E-4f,
 1.63947561694133579842E-3f,
-4.32430999505057594430E-3f,
 1.05464603945949983183E-2f,
-2.37374148058994688156E-2f,
 4.93052842396707084878E-2f,
-9.49010970480476444210E-2f,
 1.71620901522208775349E-1f,
-3.04682672343198398683E-1f,
 6.76795274409476084995E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

__ATTR_ALIGN__(32) const float B[] =
{
 3.39623202570838634515E-9f,
 2.26666899049817806459E-8f,
 2.04891858946906374183E-7f,
 2.89137052083475648297E-6f,
 6.88975834691682398426E-5f,
 3.36911647825569408990E-3f,
 8.04490411014108831608E-1f
};
float y;

if( x < 0 )
	x = -x;
if( x <= 8.0f )
	{
	y = 0.5f*x - 2.0f;
	return( ceph_expf(x) * chbevlf( y, A, 18 ) );
	}

return(  ceph_expf(x) * chbevlf( 32.0f/x - 2.0f, B, 7 ) / ceph_sqrtf(x) );
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float i0ef( float x ) {
__ATTR_ALIGN__(64) const float A[] =
{
-1.30002500998624804212E-8f,
 6.04699502254191894932E-8f,
-2.67079385394061173391E-7f,
 1.11738753912010371815E-6f,
-4.41673835845875056359E-6f,
 1.64484480707288970893E-5f,
-5.75419501008210370398E-5f,
 1.88502885095841655729E-4f,
-5.76375574538582365885E-4f,
 1.63947561694133579842E-3f,
-4.32430999505057594430E-3f,
 1.05464603945949983183E-2f,
-2.37374148058994688156E-2f,
 4.93052842396707084878E-2f,
-9.49010970480476444210E-2f,
 1.71620901522208775349E-1f,
-3.04682672343198398683E-1f,
 6.76795274409476084995E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

__ATTR_ALIGN__(32) const float B[] =
{
 3.39623202570838634515E-9f,
 2.26666899049817806459E-8f,
 2.04891858946906374183E-7f,
 2.89137052083475648297E-6f,
 6.88975834691682398426E-5f,
 3.36911647825569408990E-3f,
 8.04490411014108831608E-1f
};
float y;
if( x < 0 )
	x = -x;
if( x <= 8.0f )
	{
	y = 0.5f*x - 2.0f;
	return( chbevlf( y, A, 18 ) );
	}

return(  chbevlf( 32.0f/x - 2.0f, B, 7 ) / ceph_sqrtf(x) );
}

/*
   	Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i1f();
 *
 * y = i1f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float i1f(const float xx) {
/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */

__ATTR_ALIGN__(64) const float A[] =
{
 9.38153738649577178388E-9f,
-4.44505912879632808065E-8f,
 2.00329475355213526229E-7f,
-8.56872026469545474066E-7f,
 3.47025130813767847674E-6f,
-1.32731636560394358279E-5f,
 4.78156510755005422638E-5f,
-1.61760815825896745588E-4f,
 5.12285956168575772895E-4f,
-1.51357245063125314899E-3f,
 4.15642294431288815669E-3f,
-1.05640848946261981558E-2f,
 2.47264490306265168283E-2f,
-5.29459812080949914269E-2f,
 1.02643658689847095384E-1f,
-1.76416518357834055153E-1f,
 2.52587186443633654823E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */

 __ATTR_ALIGN__(32) const float B[] =
{
-3.83538038596423702205E-9f,
-2.63146884688951950684E-8f,
-2.51223623787020892529E-7f,
-3.88256480887769039346E-6f,
-1.10588938762623716291E-4f,
-9.76109749136146840777E-3f,
 7.78576235018280120474E-1f,
 0.0F
};
float x, y, z;
x = xx;
z = ceph_fabsf(x);
if( z <= 8.0f )
	{
	y = 0.5f*z - 2.0f;
	z = chbevlf( y, A, 17 ) * z * ceph_expf(z);
	}
else
	{
	z = ceph_expf(z) * chbevlf( 32.0f/z - 2.0f, B, 7 ) / ceph_sqrtf(z);
	}
if( x < 0.0f )
	z = -z;
return( z );
}

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float i1ef( const float xx ) {
float x, y, z;
__ATTR_ALIGN__(64) const float A[] =
{
 9.38153738649577178388E-9f,
-4.44505912879632808065E-8f,
 2.00329475355213526229E-7f,
-8.56872026469545474066E-7f,
 3.47025130813767847674E-6f,
-1.32731636560394358279E-5f,
 4.78156510755005422638E-5f,
-1.61760815825896745588E-4f,
 5.12285956168575772895E-4f,
-1.51357245063125314899E-3f,
 4.15642294431288815669E-3f,
-1.05640848946261981558E-2f,
 2.47264490306265168283E-2f,
-5.29459812080949914269E-2f,
 1.02643658689847095384E-1f,
-1.76416518357834055153E-1f,
 2.52587186443633654823E-1f
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */

 __ATTR_ALIGN__(32) const float B[] =
{
-3.83538038596423702205E-9f,
-2.63146884688951950684E-8f,
-2.51223623787020892529E-7f,
-3.88256480887769039346E-6f,
-1.10588938762623716291E-4f,
-9.76109749136146840777E-3f,
 7.78576235018280120474E-1f,
 0.0F
};
x = xx;
z = ceph_fabsf(x);
if( z <= 8.0f )
	{
	y = 0.5f*z - 2.0f;
	z = chbevlf( y, A, 17 ) * z;
	}
else
	{
	z = chbevlf( 32.0f/z - 2.0f, B, 7 ) / sqrtf(z);
	}
if( x < 0.0f )
	z = -z;
return( z );
}

/*
   	Incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamf();
 *
 * y = igamf( a, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
*/

__ATTR_PURE__
__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
__ATTR_ALIGN__(32)
static inline
float igamcf( const float aa, const float xx ) {
float a, x, ans, c, yc, ax, y, z;
float pk, pkm1, pkm2, qk, qkm1, qkm2;
float r, t;
const float big =  16777216.0f;
a = aa;
x = xx;
if( (x <= 0) || ( a <= 0) )
	return( 1.0 );

if( (x < 1.0) || (x < a) )
	return( 1.0 - ceph_igamf(a,x) );

ax = a * ceph_logf(x) - x - ceph_lgamf(a);
if( ax < -88.72283905206835 )
	{

	return( 0.0 );
	}
ax = ceph_expf(ax);

/* continued fraction */
y = 1.0 - a;
z = x + y + 1.0;
c = 0.0;
pkm2 = 1.0;
qkm2 = x;
pkm1 = x + 1.0;
qkm1 = z * x;
ans = pkm1/qkm1;

do
	{
	c += 1.0;
	y += 1.0;
	z += 2.0;
	yc = y * c;
	pk = pkm1 * z  -  pkm2 * yc;
	qk = qkm1 * z  -  qkm2 * yc;
	if( qk != 0 )
		{
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
	if( eph_fabsf(pk) > big )
		{
		pkm2 *= 5.9604644775390625E-8;
		pkm1 *= 5.9604644775390625E-8;
		qkm2 *= 5.9604644775390625E-8;
		qkm1 *= 5.9604644775390625E-8;
		}
	}
while( t > 5.9604644775390625E-8);

return( ans * ax );
}

/* left tail of incomplete gamma function:
 *
 *          inf.      k
 *   a  -x   -       x
 *  x  e     >   ----------
 *           -     -
 *          k=0   | (a+k+1)
 *
 */


float igamf( const float aa, const float xx ){
float a, x, ans, ax, c, r;
a = aa;
x = xx;
if( (x <= 0) || ( a <= 0) )
	return( 0.0 );

if( (x > 1.0) && (x > a ) )
	return( 1.0 - ceph_igamcf(a,x) );

/* Compute  x**a * exp(-x) / gamma(a)  */
ax = a * ceph_logf(x) - x - ceph_lgamf(a);
if( ax < -88.72283905206835;)
	{

	return( 0.0 );
	}
ax = ceph_expf(ax);

/* power series */
r = a;
c = 1.0;
ans = 1.0;

do
	{
	r += 1.0;
	c *= x/r;
	ans += c;
	}
while( c/ans > 5.9604644775390625E-8 );

return( ans * ax/a );
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
float * __restrict p;
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
