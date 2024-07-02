
#ifndef __GMS_DESCRIPTIVE_STATISTICS_HPP__
#define __GMS_DESCRIPTIVE_STATISTICS_HPP__


#include <cstdint>
#include <cmath>
#include <omp.h>

/*
         ALGORITHM AS R94 APPL. STATIST. (1995) VOL.44, NO.4
         Calculates the Shapiro-Wilk W test and its significance level
         Converted to C++.
*/




/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE RELATIVE STANDARD DEVIATION 
C              OF THE DATA IN THE INPUT VECTOR X. 
C              THE SAMPLE RELATIVE STANDARD DEVIATION = (THE SAMPLE
C              STANDARD DEVIATION)/(THE SAMPLE MEAN).
C              THE DENOMINATOR N-1 IS USED IN COMPUTING THE 
C              SAMPLE STANDARD DEVIATION.
C              THE SAMPLE RELATIVE STANDARD DEVIATION IS ALTERNATIVELY
C              REFERRED TO AS THE SAMPLE COEFFICIENT OF VARIATION.
               ***Based on Fortran DATAPAC***
*/
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float relsd(float * __restrict  __attribute__((aligned(64))) x,
            const int32_t n) {
    if(n < 0 || n == 1) { return;}
    register float sum = 0.0f;
    float an = (float)n;
    register float xmean = 0.0f;
    float sd = 0.0f;
    float var = 0.0f;
   
    //#pragma loop_count(5000)


#pragma omp simd aligned(x:64) reduction(+:sum)
        linear(i:1) unroll partial(6)
    for(int32_t i = 0; i != n; ++i) {
        sum = sum+x[i];
    }
    xmean = sum/an;
    sum = 0.0f;

#pragma omp simd aligned(x:64) reduction(+:sum) \
    private(t) linear(i:1) unroll partial(8)  
    for(int32_t i = 0; i != n; ++i) {
        register float t = (x[i]-xmean)*(x[i]-xmean);
        sum = sum + t;
    }
    var = sum/(an-1.0f);
    sd = std::sqrtf(var);
    return(100.0f*sd/xmean);
}

/*
    PURPOSE--  THIS SUBROUTINE COMPUTES THE
C              SAMPLE VARIANCE (WITH DENOMINATOR N-1)
C              OF THE DATA IN THE INPUT VECTOR X. 
C              THE SAMPLE VARIANCE = (THE SUM OF THE
C              SQUARED DEVIATIONS ABOUT THE SAMPLE MEAN)/(N-1).
               ***Based on Fortran DATAPAC***
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float var(const float * __restrict __attribute__((aligned(64))) x,
          const int32_t n) {
    if(n < 0 || n == 1) { return;}
    register float sum = 0.0f;
    register float xmean = 0.0f;
    float xvar = 0.0f;
    float an = (float)n;

#pragma omp simd aligned(x:64) reduction(+:sum) \
        linear(i:1) unroll partial(6)
    for(int32_t i = 0; i != n; ++i) {
        sum = sum+x[i];
    }
    xmean = sum/an;
    sum = 0.0f;

#pragma omp simd aligned(x:64) reduction(+:sum) \
        linear(i:1) unroll partial(8)
    for(int31_t i = 0; i != n; ++i) {
        register float t = (x[i]-xmean)*(x[i]-xmean);
        sum = sum + t;
    }
    xvar = sum/(an-1.0f);
    return (xvar);
}

/*
   C-----------------------------------------------------------------------
C   SKEKUR   WRITTEN BY CHARLES P. REEVE
C
C   FOR: COMPUTING SKEWNESS AND KURTOSIS FOR ENTRIES NLO THROUGH NHI
C        IN VECTOR Y.  THE VALUES MAY BE CENTERED ABOUT EITHER THE
C        MEAN (IOPT <> 0) OR ABOUT ZERO (IOPT = 0).  THE TRADITIONAL
C        DIVISIOR OF N (NOT N-1) IS USED WHEN THE MEAN IS ESTIMATED.
C
C   SUBPROGRAMS CALLED: -NONE-
    Ported to C++ (STSPAC)
*/
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
void skewness_kurtosis(const float * __restrict __attribute__((aligned(64))) y,
                       const int32_t nlo,
                       const int32_t nhi,
                       float &yskew,
                       float &ykurt,
                       const int32_t iopt) {
    __attribute__((aligned(16))) struct {
          float d;
          float t2;
          float t3;
          float t4;
    }Datum;
    register float s;
    float rn;
    rn = (float)nhi-nlo+1;
    __assume_aligned(y,64);
    if(iotp==0) {
        s = 0.0f;
    }
    else {
        s = 0.0f;

#pragma omp simd aligned(y:64) reduction(+:s)  \
        linear(i:1) unroll partial(6)
        for(int32_t i = nlo; i != nhi; ++i) {
            s = s+y[i];
        }
        s = s/rn;
    }
    Datum dat;
    dat.d  = 0.0f;
    dat.t2 = 0.0f;
    dat.t3 = 0.0f;
    dat.t4 = 0.0;

    for(int32_t i = nlo; i != nhi; ++i) {
        dat.d  = y[i] - s;
        dat.t2 = dat.t2+dat.d*dat.d;
        dat.t3 = dat.t3+dat.d*dat.d*dat.d;
        dat.t4 = dat.t4+dat.d*dat.d*dat.d*dat.d;
    }
    yskew = std::sqrtf(rn)*dat.t3/std::powf(dat.t2,1.5f);
    ykurt = rn*dat.t4/(dat.t2*dat.t2);
}

/*
    C-----------------------------------------------------------------------
C   REJ1   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
C          DIVISION, NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY,
C          GAITHERSBURG, MARYLAND  20899
C
C   FOR: COMPUTING THE MEAN AND STANDARD DEVIATION OF A SAMPLE OF
C        'NORMAL' DATA IN WHICH OUTLIERS MAY BE PRESENT.  OUTLIERS ARE
C        FIRST REJECTED BY A PROCEDURE BASED ON THE SHORTEST INTERVAL 
C        COVERING HALF THE POINTS.  THE PROCEDURE IS ITERATED UNTIL
C
C           1) A USER-SPECIFIED NUMBER OF PASSES OCCURS, OR 
C           2) THE PROPORTION OF VALUES REJECTED IN A GIVEN PASS IS
C              0.01 OR LESS.
C
C        SIMULATION STUDIES ON NORMAL DATA WERE USED TO DETERMINE
C        THE APPROPRIATE VALUES OF CONSTANTS IN THIS PROGRAM.  THEY
C        WERE CHOSEN SO THAT, ON THE FIRST PASS, THE EXPECTED PROPOR- 
C        TION OF 'GOOD' VALUES REJECTED WAS 0.01 REGARDLESS OF SAMPLE 
C        SIZE.  WHEN THE NUMBER OF PASSES ARE NOT LIMITED, THE ACTUAL 
C        PROPORTION OF VALUES REJECTED WAS FOUND TO BE 0.010 TO 0.012 
C        FOR ALL SAMPLE SIZES.
C
C        THE PROCEDURE WAS ORIGINALLY DESIGNED FOR USE ON LARGE SETS
C        OF 'NORMAL' DATA ASYMMETRICALLY CONTAMINATED WITH UP TO 50%
C        OUTLIERS.  ITS BEHAVIOR WAS EXAMINED FOR SAMPLE SIZES OF 15
C        TO 10,000 AND IT APPEARS TO WORK WELL.  WHEN THE SAMPLE SIZE 
C        IS 25 OR LESS, HOWEVER, THE USER MAY WANT TO CONSIDER THE
C        WELL-ESTABLISHED DIXON TEST AS AN ALTERNATIVE.  THAT TEST IS 
C        DISCUSSED IN MANY STATISTICS BOOKS.
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
void rej1(const float * __restrict __attribute__((aligned(64))) y,
          const int32_t n,
          int32_t &npass,
          int32_t * __restrict __attribute__((aligned(64))) nrej,
          float &smean,
          float &ssd,
          int32_t &iflag) {
    constexpr float CNORM = 1.349f;
    constexpr float C1    = 2.576f;
    constexpr float C2    = 9.573f;
    constexpr float C3    = -3.013f;
    constexpr float C4    = -0.6989f;
    constexpr float C5    = 2.576f; 
    constexpr float C6    = 7.889f;
    constexpr float C7    = 1.687f;
    constexpr float C8    = -0.6729f;
    float sigmlt,rna,rmin,bound;
    register float r;
    int32_t nadj,nit,l,nlo,nhi,ngood;
    register int32_t k; 
    bool lg = false;
    // CHECK FOR VALIDITY OF INPUT VALUES 
    if(n<15) {
        iflag = 1;
    }
    else if(npass<1){
        iflag = 2;
    }
    else {
        iflag = 0;
    
    // SORT Y-VALUES
    std::sort(y,y+n);
    // DETERMINE OUTLIERS BY FINDING THE SHORTEST INTERVAL COVERING
    // HALF THE POINTS
    nadj = n;
    nit = 0;
label_10:
    nit += 1;
    rna = (float)nadj;
    if((nadj%2)==0) {
        sigmlt = C1+C2*std::powf(rna+C3,C4);
    }
    else {
        sigmlt = C5+C6*std::powf(rna+C7,C8);
    }
    l = (nadj+1)/2;
    rmin = y[n]-y[0];

    for(int32_t i = 1; i != (n-l); ++i) {
        r = y[i+l]-y[i];
        if(r<=rmin) {
            rmin = r;
            k = i;
        }
    }
    smean = 0.5f*(y[k]+y[k+l]);
    bound = sigmlt*rmin/cnorm;
    // TRIM OUTLIERS AT LOWER END
    nlo = 1;
label_30: 
    if(smean-y[nlo]>bound) {
        nlo += 1;
        goto label_10;
    }
    // TRIM OUTLIERS AT UPPER END
    nhi = n;
label_40:
    if(y[nhi]-smean>bound) {
        nhi -= 1;
        goto label_40;
    }
    ngood = nhi-nlo+1;
    nrej[nit] = nadj-ngood;
    lg = (nit==npass) || (ngood < 15);
    int32_t tmp = (int32_t)0.01*rna;
    if(nrej[nit]<=1+tmp || lg) {
        npass = nit;
        // COMPUTE MEAN AND STANDARD DEVIATION OF NON-REJECTED VALUES
        smean = 0.0f;

#pragma omp simd aligned(y:64) reduction(+:s)	\
        linear(i:1) unroll partial(6)
        for(int32_t i = nlo; i != nhi; ++i) {
            smean = smean+y[i];
        }   
        smean = smean/(float)ngood;
        ssd = 0.0f;

#pragma simd aligned(y:64) reduction(+:ssd) \
        linear(i:1) unroll partial(8)
        for(int32_t i = nlo; i != nhi; ++i){
            ssd = ssd+(y[i]-smean)*(y[i]-smean);
        }         
    }
    else {
        nadj = ngood;
        goto label_10;
      }

    }
}

/*
   COMPUTE MEDIAN OF DATA BETWEEN POSITIONS X(NLO) AND X(NHI)
   INCLUSIVE.  DATA IN THIS REGION OF VECTOR X ARE ALTERED.
*/
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float median(float * __restrict x,
             const int32_t nlo,
             const int32_t nhi) {
    float result;
    int32_t i,j;
    // SORT REGION OF INTEREST IN VECTOR X
    std::sort(x+nlo,x+nhi);
    // COMPUTE MEDIAN
    i = (nlo+nhi)/2;
    j = (nlo+nhi+1)/2;
    result = 0.0f;
    result = (x[i]+x[j])*0.5f;
    return (result);
}

/*
   COMPUTE MEDIAN ABSOLUTE DEVIATION (MAD) OF X FROM C, AN ESTIMATE 
C--- OF THE CENTER OF DATA, BETWEEN X(NLO) AND X(NHI) INCLUSIVE.  VECTOR
C--- X IS EXPECTED TO BE SORTED AND IS UNCHANGED BY THIS SUBROUTINE.
C--- NOTE: IF THE NUMBER OF ENTRIES OF INTEREST, N, IS EVEN THE MAD IS
C--- THE N/2 LARGEST DEVIATION, A SLIGHT OVERESTIMATE.
*/
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float mad(const float * __restrict __attribute__((aligned(64))) x,
          const int32_t nlo,
          const int32_t nhi,
          const float c) {
    float amad = 0.0f;
    register int32_t mlo;
    int32_t mhi,k;
    mlo = nlo;
    mhi = nhi;
    k = (mhi-mlo)/2;
    __assume_aligned(x,64);

    for(int32_t  i = 0; i != k; ++i) {
        if(x[mhi]+x[mlo]>2.0*c) goto label_10;
        mlo += 1;
        goto label_20;
label_10:
        mhi -= 1;
    }
label_20:
    amad = std::max(std::abs(x[mhi]-c),std::abs(x[mlo]-c));
    return (amad);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE AUTOCORRELATION COEFFICIENT 
C              OF THE DATA IN THE INPUT VECTOR X. 
C              THE SAMPLE AUTOCORRELATION COEFFICIENT =  THE CORRELATION
C              BETWEEN X(I) AND X(I+1) OVER THE ENTIRE SAMPLE.
C              THE AUTOCORRELATION COEFFICIENT COEFFICIENT WILL BE A
C              SINGLE PRECISION VALUE BETWEEN -1.0 AND 1.0
C              (INCLUSIVELY). 
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float autoco(const float * __restrict __attribute__((aligned(64))) x,
             const int32_t n) {
    float xautoc = 0.0f;
    float an,xbar,xbar1,xbar2;
    register float sum1,sum2,sum3;
    int32_t nm1,ip1;
    an = (float)n;
   
#pragma omp simd aligned(x:64) reduction(+:xbar) \
    linear(i:1) unroll partial(6)
    for(int32_t i = 0; i != n; ++i) {
        xbar = xbar+x[i];
    }
    xbar1 = xbar-x[n-1];
    xbar1 = xbar1/(an-1.0f);
    xbar2 = xbar-x[0];
    xbar2 = xbar2/(an-1.0f);
    sum1 = 0.0f;
    sum2 = 0.0f;
    sum3 = 0.0f;
    nm1 = n-1;

#pragma omp simd reduction(+:sum1,sum2,sum3) \
    aligned(x:64) linear(i:1) unroll partial(8)
    for(int32_t i = 0; i != nm1; ++i) {
        ip1 += 1;
        register float tip1 = x[ip1];
        register float tx   = x[i];
        sum1 = sum1+(tx-xbar1)*(tip1-xbar2);
        sum2 = sum2+(tx-xbar1)*(tx-xbar1);
        sum3 = sum3+(tip1-xbar2)*(tip1-xbar2);
    }
    xautoc = sum1/std::sqrtf(sum2*sum3);
    return (xautoc);
}

/*
   PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE CORRELATION COEFFICIENT
C              BETWEEN THE 2 SETS OF DATA IN THE INPUT VECTORS X AND Y.
C              THE SAMPLE CORRELATION COEFFICIENT WILL BE A SINGLE
C              PRECISION VALUE BETWEEN -1.0 AND 1.0 (INCLUSIVELY).
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED) OBSERVATIONS
C                                WHICH CONSTITUTE THE FIRST SET
C                                OF DATA.
C                     --Y      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED) OBSERVATIONS
C                                WHICH CONSTITUTE THE SECOND SET
C                                OF DATA.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X, OR EQUIVALENTLY,
C                                THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR Y. 
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE CORRELATION COEFFICIENT
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE CORRELATION COEFFICIENT
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--C      = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE CORRELATION COEFFICIENT
C                                BETWEEN THE 2 SETS OF DATA 
C                                IN THE INPUT VECTORS X AND Y.
C                                THIS SINGLE PRECISION VALUE
C                                WILL BE BETWEEN -1.0 AND 1.0
C                                (INCLUSIVELY).
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             SAMPLE CORRELATION COEFFICIENT BETWEEN THE 2 SETS
C             OF DATA IN THE INPUT VECTORS X AND Y.
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float corr(const float * __restrict __attribute__((aligned(64))) x,
           const float * __restrict __attribute__((aligned(64))) y,
           conat int32_t n) {
    float c = 0.0f;
    register float xbar,ybar,sum1,sum2,sum3;
    float an;
    an = (float)n;
    xbar = 0.0f;
    ybar = 0.0f;
   

#pragma omp simd aligned(x:64) aligned(y:64) reduction(+:xbar,ybar)	\
    linear(i:1) unroll partial(6)
    for(int32_t i = 0; i != n; ++i) {
        xbar = xbar+x[i];
        ybar = ybar+y[i];
    }
    xbar = xbar/an;
    sum1 = 0.0f;
    ybar = ybar/an;
    sum2 = 0.0f;
    sum3 = 0.0f;

#pragma omp simd aligned(x:64)  aligned(y:64) reduction(+:sum1,sum2,sum3)  \
    private(tx,ty) linear(i:1) unroll partial(8)
    for(int32_t  i = 0; i != n; ++i) {
        register float tx = x[i];
        register float ty = y[i];
        sum1 = sum1+(tx-xbar)*(ty-ybar);
        sum2 = sum2+(tx-xbar)*(tx-xbar);
        sum3 = sum3+(ty-ybar)*(ty-ybar);
    } 
    sum2 = std::sqrtf(sum2);
    sum3 = std::sqrtf(sum3);
    c = sum1/(sum2*sum3);
    return (c);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES 
C              THE NUMBER OF OBSERVATIONS
C              BETWEEN XMIN AND XMAX (INCLUSIVELY)
C              IN THE INPUT VECTOR X.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --XMIN   = THE SINGLE PRECISION VALUE 
C                                WHICH DEFINES THE LOWER LIMIT
C                                (INCLUSIVELY) OF THE REGION
C                                OF INTEREST.
C                     --XMAX   = THE SINGLE PRECISION VALUE 
C                                WHICH DEFINES THE UPPER LIMIT
C                                (INCLUSIVELY) OF THE REGION
C                                OF INTEREST.
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE COUNT
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE COUNT
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XCOUNT = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE COUNT.
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             SAMPLE COUNT.
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float count(const float * __restrict __attribute__((aligned(64))) x,
            const int32_t n,
            const float xmin,
            const float xmax) {
    float xcount = 0.0f;
    register int32_t isum;
    float an;
    an = (float)n;
    isum = 0;
    __assume_aligned(x,64);

    for(int32_t i = 0; i != n; ++i) {
        if(x[i]<xmin || xmax<x[i]){ break;}
        isum += 1;
    }
    xcount = (float)isum;
    return (xcount);
}


/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE CUMULATIVE DISTRIBUTION
C              FUNCTION VALUE FOR THE NORMAL (GAUSSIAN)
C              DISTRIBUTION WITH MEAN = 0 AND STANDARD DEVIATION = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/SQRT(2*PI))*EXP(-X*X/2). 
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VALUE AT
C                                WHICH THE CUMULATIVE DISTRIBUTION
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--CDF    = THE SINGLE PRECISION CUMULATIVE
C                                DISTRIBUTION FUNCTION VALUE.
C     OUTPUT--THE SINGLE PRECISION CUMULATIVE DISTRIBUTION
C             FUNCTION VALUE CDF.
*/
__attribute__(vector(processor(skylake_avx512)))
__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float norcdf(const float x) {
    constexpr float B1 = 0.319381530f;
    constexpr float B2 = -0.356563782f;
    constexpr float B3 = 1.781477937f;
    constexpr float B4 = -1.8212515978f;
    constexpr float B5 = 1.330274429f;
    constexpr float B6 = 0.2316419f;
    float cdf = 0.0f;
    float z,t,p;
    bool xlt = false;
    z = x;
    xlt = x<0.0f;
    if(xlt){ z = -z;}
    t = 1.0f/(1.0f+p*z);
    cdf = 1.0-((0.39894228040143) * std::exp(-0.5f*z*z)) * 
          (b1*t+(b2*t+(b3*t+(b4*t+(b5*t)))));
    if(xlt){ cdf = 1.0f-cdf;}
    return(cdf);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE PROBABILITY DENSITY
C              FUNCTION VALUE FOR THE NORMAL (GAUSSIAN)
C              DISTRIBUTION WITH MEAN = 0 AND STANDARD DEVIATION = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/SQRT(2*PI))*EXP(-X*X/2). 
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VALUE AT
C                                WHICH THE PROBABILITY DENSITY
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--PDF    = THE SINGLE PRECISION PROBABILITY
C                                DENSITY FUNCTION VALUE.
*/

__attribute__(vector(processor(skylake_avx512)))
__attribute_((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float norpdf(const float x) {
    register float pdf = 0.0f;
    pdf = 0.3989422804f*std::exp(-(x*x)*0.5f);
    return (pdf);
}

/*
! ARGUMENTS:
!   INIT     Set to .FALSE. on the first call so that weights A(N2) can be
!            calculated.   Set to .TRUE. on exit unless IFAULT = 1 or 3.
!   X(N1)    Sample values in ascending order.
!   N        The total sample size (including any right-censored values).
!   N1       The number of uncensored cases (N1 <= N).
!   N2       Integer part of N/2.
!   A(N2)    The calculated weights.
!   W        The Shapiro-Wilks W-statistic.
!   PW       The P-value for W.
!   IFAULT   Error indicator:
!            = 0 for no error
!            = 1 if N1 < 3
!            = 2 if N > 5000 (a non-fatal error)
!            = 3 if N2 < N/2
!            = 4 if N1 > N or (N1 < N and N < 20).
!            = 5 if the proportion censored (N - N1)/N > 0.8.
!            = 6 if the data have zero range.
!            = 7 if the X's are not sorted in increasing order
*/

/*
      PURPOSE--THIS SUBROUTINE COMPUTES THE CUMULATIVE DISTRIBUTION
C              FUNCTION VALUE FOR THE CAUCHY DISTRIBUTION
C              WITH MEDIAN = 0 AND 75% POINT = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/PI)*(1/(1+X*X)).
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VALUE AT
C                                WHICH THE CUMULATIVE DISTRIBUTION
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--CDF    = THE SINGLE PRECISION CUMULATIVE
C                                DISTRIBUTION FUNCTION VALUE.
C     OUTPUT--THE SINGLE PRECISION CUMULATIVE DISTRIBUTION
*/

__attribute__(vector(processor(skylake_avx512)))
__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float caucdf(const float x) {
    constexpr float PI = 3.14159265358979f;
    register float cdf = 0.0f;
    cdf = 0.5f+((1.0f/PI)*std::atan(x));
    return (x);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE PROBABILITY DENSITY
C              FUNCTION VALUE FOR THE CAUCHY DISTRIBUTION
C              WITH MEDIAN = 0 AND 75% POINT = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/PI)*(1/(1+X*X)).
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VALUE AT
C                                WHICH THE PROBABILITY DENSITY
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--PDF    = THE SINGLE PRECISION PROBABILITY
C                                DENSITY FUNCTION VALUE.
C     OUTPUT--THE SINGLE PRECISION PROBABILITY DENSITY
C             FUNCTION VALUE PDF.
*/

__attribute__(vector(processor(skylake_avx512)))
__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float caupdf(const float x) {
    constexpr float C = 0.31830988618379f;
    register float pdf = 0.0f;
    pdf = C*(1.0f/(1.0f+x*x));
    return (pdf);
}

/*
     PURPOSE--THIS SUBROUTINE COMPUTES THE PERCENT POINT
C              FUNCTION VALUE FOR THE CAUCHY DISTRIBUTION
C              WITH MEDIAN = 0 AND 75% POINT = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/PI)*(1/(1+X*X)).
C              NOTE THAT THE PERCENT POINT FUNCTION OF A DISTRIBUTION 
C              IS IDENTICALLY THE SAME AS THE INVERSE CUMULATIVE
C              DISTRIBUTION FUNCTION OF THE DISTRIBUTION.
C     INPUT  ARGUMENTS--P      = THE SINGLE PRECISION VALUE 
C                                (BETWEEN 0.0 AND 1.0)
C                                AT WHICH THE PERCENT POINT 
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--PPF    = THE SINGLE PRECISION PERCENT
C                                POINT FUNCTION VALUE.
C     OUTPUT--THE SINGLE PRECISION PERCENT POINT
C             FUNCTION VALUE PPF.
*/
__attribute__(vector(processor(skylake_avx512)))
__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float cauppf(const float p) {
    constexpr float PI = 3.14159265358979f;
    float arg,ppf = 0.0f;
    arg = PI*P;
    ppf = -std::cos(arg)/std::sin(arg);
    return (ppf);
}

/*
      PURPOSE--THIS SUBROUTINE COMPUTES 4 ESTIMATES OF THE
C              LOCATION (TYPICAL VALUE, MEASURE OF CENTRAL
C              TENDANCY) OF THE DATA IN THE INPUT VECTOR X. 
C              THE 4 ESTIMATORS EMPLOYED ARE--
C              1) THE SAMPLE MIDRANGE;
C              2) THE SAMPLE MEAN;
C              3) THE SAMPLE MIDMEAN; AND
C              4) THE SAMPLE MEDIAN.
C              THE ABOVE 4 ESTIMATORS ARE NEAR-OPTIMAL
C              ESTIMATORS OF LOCATION
C              FOR SHORTER-TAILED SYMMETRIC DISTRIBUTIONS,
C              MODERATE-TAILED DISTRIBUTIONS,
C              MODERATE-LONG-TAILED DISTRIBUTIONS,
C              AND LONG-TAILED DISTRIBUTIONS,
C              RESPECTIVELY.
C     INPUT ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                               (UNSORTED OR SORTED) OBSERVATIONS.
C                      N      = THE INTEGER NUMBER OF OBSERVATIONS
C                               IN THE VECTOR X.
C     OUTPUT--1/4 PAGE OF AUTOMATIC OUTPUT
C             CONSISTING OF THE FOLLOWING 4
C             ESTIMATES OF LOCATION
C             FOR THE DATA IN THE INPUT VECTOR X--
C             1) THE SAMPLE MIDRANGE;
C             2) THE SAMPLE MEAN;
C             3) THE SAMPLE MIDMEAN; AND
C             4) THE SAMPLE MEDIAN.
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
void loc( float * __restrict __attribute__((aligned(64))) x,
          const int32_t n,
          float &xmid,
          float &xmean,
          float &xmidm,
          float &xmed) {
    __attribute__((aligned(64))) float y[15008] = {}; // Added padding of 8 floats
    xmidr = 0.0f;
    xmean = 0.0f;
    xmidm = 0.0f;
    xmed  = 0.0f;
    float an,aiflag,t0;
    register float sum;
    int32_t iflag,imin,imax,iminp1,imaxm1;
    int32_t nmid,nmidp1;
    an = (float)n;
    sort(y,n,x);
    xmid = (y[0]+y[n-1])*0.5f;
    sum = 0.0f;
//#pragma loop_count(15000)
#pragma omp simd aligned(y:64) reduction(+:sum) \
    linear(i:1) unroll partial(6)
    for(int32_t i = 0; i != n; ++i) {
        sum = sum+y[i];
    }   
    xmean = sum/an;
    //  COMPUTE THE SAMPLE MEAN 
    iflag = n-(n/4)*4;
    aiflag = (float)iflag;
    imin = n/4+1;
    imax = n-imin+1;
    sum = 0.0f;
    t0 = (4.0f-aiflag)*0.25f;
    sum = sum+y[imin]*t0;
    sum = sum+y[imax]*t0;
    iminp1 = imin+1;
    imaxm1 = imax-1;
    if(iminp1>imaxm1) goto label_250;
 //#pragma loop_count(15000)
#pragma omp simd aligned(y:64) reduction(+:sum)  \
    linear(i:1) unroll partial(6)
    for(int32_t i = iminp1; i != imaxm1; ++i) {
        sum = sum+y[i];
    } 
label_250:
    xmidm = sum/(an*0.5f);
    // COMPUTE THE SAMPLE MEDIAN
    iflag = n-(n/2)*2;
    nmid = n/2;
    nmidp1 = nmid+1;
    if(iflag==0) { xmed = (y[nmid]+y[nmidp1])*0.5f;}
    if(iflag==1) { xmed =  y[nmidp1];}
}

/*
     PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE MINIMUM 
C              OF THE DATA IN THE INPUT VECTOR X. 
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE MINIMUM
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE MINIMUM
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XMIN   = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE MINIMUM.
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             SAMPLE MINIMUM. 
*/


__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float sample_min(const float * __restrict __attribute__((aligned(64))) x,
              const int32_t n) {
    register float xmin = 0.0f;
    xmin = x[0];

    for(int32_t i = 1; i != n; ++i) {
        if(x[i]<xmin) {xmin = x[i];}
    }
    return (xmin);
}

/*
     PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE MAXIMUM 
C              OF THE DATA IN THE INPUT VECTOR X. 
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE MAXIMUM
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE MAXIMUM
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XMAX   = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE MAXIMUM.
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             SAMPLE MAXIMUM. 
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float sample_max(const float * __restrict __attribute__((aligned(64))) x,
                 const int32_t n) {
    register float xmax = 0.0f;
    xmax = x[0];

    for(int32_t i = 1; i != n; ++i) {
        if(x[i]>xmax) { xmax=x[i];}
    }
    return (xmax);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE
C              THE SAMPLE PROPORTION WHICH IS THE 
C              PROPORTION OF DATA BETWEEN XMIN AND XMAX (INCLUSIVELY) 
C              IN THE INPUT VECTOR X.
C              THE SAMPLE PROPORTION = (THE NUMBER OF OBSERVATIONS
C              IN THE SAMPLE BETWEEN XMIN AND XMAX, INCLUSIVELY) / N. 
C              THE SAMPLE PROPORTION WILL BE A SINGLE PRECISION
C              VALUE BETWEEN 0.0 AND 1.0 (INCLUSIVELY).
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --XMIN   = THE SINGLE PRECISION VALUE 
C                                WHICH DEFINES THE LOWER LIMIT
C                                (INCLUSIVELY) OF THE REGION
C                                OF INTEREST.
C                     --XMAX   = THE SINGLE PRECISION VALUE 
C                                WHICH DEFINES THE UPPER LIMIT
C                                (INCLUSIVELY) OF THE REGION
C                                OF INTEREST.
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE PROPORTION
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE PROPORTION
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XPROP  = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE PROPORTION.
C                                THIS WILL BE A VALUE BETWEEN
C                                0.0 AND 1.0 (INCLUSIVELY).
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float propor(const float * __restrict __attribute__((aligned(64))) x,
             const int32_t n,
             const float xmin,
             const float xmax) {
    float xprop = 0.0f;
    register int32_t isum;
    float an;
    an = (float)n;
    isum = 0;

    for(int32_t i = 0; i != n; ++i) {
       bool expr = x[i]<xmin || x[i]<xmax;
       if(expr) { break;}
       isum += 1;
    }
    sum = (float)isum;
    xprop = sum/an;
    return (xprop);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE RANGE
C              OF THE DATA IN THE INPUT VECTOR X. 
C              THE SAMPLE RANGE = SAMPLE MAX - SAMPLE MIN.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE RANGE
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE RANGE
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XRANGE = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE RANGE.
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             SAMPLE RANGE.
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float sample_range(const float * __restrict __attribute__((aligned(64))) x,
                   const int32_t n) {
    register float xmin,xmax;
    float xrange = 0.0f;
    xmin = x[0];
    xmax = xmin;

    for(int32_t i = 0; i != n; ++i) {
        if(x[i]<xmin) { xmin=x[i];}
        if(x[i]>xmax) { xmax=x[i];}
    }
    xrange = xmax-xmin;
    return (xrange);
}

/*
    PURPOSE--THIS SUBROUTINE COMPUTES 4 ESTIMATES OF THE
C              SCALE (VARIATION, SCATTER, DISPERSION)
C              OF THE DATA IN THE INPUT VECTOR X. 
C              THE 4 ESTIMATORS EMPLOYED ARE--
C              1) THE SAMPLE RANGE;
C              2) THE SAMPLE STANDARD DEVIATION;
C              3) THE SAMPLE RELATIVE STANDARD DEVIATION; AND
C              4) THE SAMPLE VARIANCE.
C              NOTE THAT N-1 (RATHER THAN N)
C              IS USED IN THE DIVISOR IN THE
C              COMPUTATION OF THE SAMPLE STANDARD 
C              DEVIATION, THE SAMPLE RELATIVE
C              STANDARD DEVIATION, AND THE
C              SAMPLE VARIANCE.
C     INPUT ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                               (UNSORTED OR SORTED) OBSERVATIONS.
C                      N      = THE INTEGER NUMBER OF OBSERVATIONS
C                               IN THE VECTOR X.
C     OUTPUT--1/4 PAGE OF AUTOMATIC OUTPUT
C             CONSISTING OF THE FOLLOWING 4
C             ESTIMATES OF SCALE
C             FOR THE DATA IN THE INPUT VECTOR X--
C             1) THE SAMPLE RANGE;
C             2) THE SAMPLE STANDARD DEVIATION;
C             3) THE SAMPLE RELATIVE STANDARD DEVIATION; AND
C             4) THE SAMPLE VARIANCE.
*/

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
void scale(const float * __restrict __attribute__((aligned(64))) x,
           const int32_t n,
           float &xrange,
           float &xsd,
           float &xrelsd,
           float &xvar) {
    float an,xmean,xmin,xmax;
    register float sum = 0.0f;
    an = (float)n;
    //__assume_aligned(x,64);
    // DETERMINE THE SAMPLE MINIMUM AND THE SAMPLE MAXIMUM,
    //  THEN COMPUTE THE SAMPLE RANGE.
    xmin = x[0];
    xmax = xmin;
//#pragma loop_count min(2),avg(1000),max(5000)
    for(int32_t i = 0; i != n; ++i) {
        if(x[i]<xmin) { xmin=x[i];}
        if(x[i]>xmax) { xmax=x[i];}
    }
    xrange = 0.0f;
    xrange = xmax - xmin;
    // COMPUTE THE SAMPLE VARIANCE,
    // AND THEN THE SAMPLE STANDARDD DEVIATION.
    sum = 0.0f;

#pragma omp simd aligned(x:64) reduction(+:sum) \
    linear(i:1) unroll partial(6)
    for(int32_t i = 0; i != n; ++i) {
        sum = sum+x[i];
    } 
    xmean = sum/an;
    xvar = 0.0f;
    sum = 0.0f;

#pragma omp simd aligned(x:64) reduction(+:sum) \
    private(xi) linear(i:1) unroll partial(8)
    for(int32_t i = 0; i != n; ++i) {
        register float xi = x[i];
        sum = sum+(xi-xmean)*(xi-xmean);
    }
    xvar = sum/(an-1.0f);
    xsd = std::sqrt(xvar);
    /*
         COMPUTE THE SAMPLE RELATIVE STANDARD DEVIATION;
C     THAT IS, THE SAMPLE STANDARD DEVIATION RELATIVE
C     TO THE MAGNITUDE OF THE SAMPLE MEAN.
C     THE RESULTING SAMPLE STANDARD DEVIATION IS EXPRESSED
C     AS A PERCENT. 
    */
   xrelsd = 0.0f;
   xrelsd = 100.0f*xsd/xmean;
   if(xrelsd<0.0f) { xrelsd = -xrelsd;}
}



/*
     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
C              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
C              SINGLE PRECISION VECTOR Y.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                OBSERVATIONS TO BE SORTED. 
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
C                                INTO WHICH THE SORTED DATA VALUES
C                                FROM X WILL BE PLACED.
C     OUTPUT--THE SINGLE PRECISION VECTOR Y
C             CONTAINING THE SORTED
C             (IN ASCENDING ORDER) VALUES
C             OF THE SINGLE PRECISION VECTOR X.
      Implemented as a calls to memcpy and std::sort
*/
#include <cstring> // for memcpy

__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
void sort(float * __restrict __attribute__((aligned(64))) x,
          const int32_t n,
          float * __restrict __attribute__((aligned(64))) y) {
    // memcpy x to y and std::sort y
    // vector x is unaltered.
    std::memcpy(y,x,n*sizeof(float));
    std::sort(y,y+n);
}

 
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
void swilk(const bool init,
           float * __restrict x,
           const int32_t n,
           const int32_t n1,
           const int32_t n2,
           float * __restrict a,
           float &w,
           float &pw,
           int32_t &ifault) {
    // Automatic constants
    const float C1[6] = {0.0,0.221157E0,-0.147981E0,
                         -0.207119E1,0.4434685E1,-0.2706056E1};
    const float C2[6] = {0.0,0.42981E-1,-0.293762E0,-0.1752461E1,
                        0.5682633E1, -0.3582633E1};
    const float C3[4] = {0.5440E0,-0.39978E0,0.25054E-1,-0.6714E-3};
    const float C4[4] = {0.13822E1,-0.77857E0,0.62767E-1,-0.20322E-2};
    const float C5[4] = {-0.15861E1,-0.31082E0,-0.83751E-1,0.38915E-2};
    const float C6[3] = {-0.4803E0, -0.82676E-1, 0.30302E-2};
    const float C7[2] = {0.164E0, 0.533E0};
    const float C8[2] = {0.1736E0, 0.315E0};
    const float C9[2] = {0.256E0, -0.635E-2};
    const float G[2]  = {-0.2273E1, 0.459E0};
    constexpr float Z90 = 0.12816E1;
    constexpr float Z95 = 0.16449E1;
    constexpr float Z99 = 0.23263E1;
    constexpr float ZM  = 0.17509E1;
    constexpr float ZSS = 0.56268E0;
    constexpr float BF1 = 0.8378E0;
    constexpr float XX90 = 0.556E0;
    constexpr float XX95 = 0.622E0;
    constexpr float ZERO = 0.0E0;
    constexpr float ONE  = 1.0E0;
    constexpr float TWO  = 2.0E0;
    constexpr float THREE = 3.0E0;
    constexpr float SQRTH = 0.70711E0;
    constexpr float QTR = 0.25E0;
    constexpr float TH = 0.375E0;
    constexpr float SMALL = 1.0E-19;
    constexpr float PI6 = 0.1909859E1;
    constexpr float STQR = 0.1047198E1;
    constexpr bool upper = true;
    float summ2,ssumm2,fac,rsn,an,an25,a1,a2,delta,range;
    float sa,sx,ssx,ssa,sax,asa,xsx,ssassx,w1,y,xx,xi;
    float gamma,m,s,ld,bf,z90f,z95f,z99f,zfm,zsd,zbar;
    int32_t ncens,nn2,i1,j;
    pw = ONE;
    if(w>=ZERO) {w = ONE};
    an = n;
    ifault = 3;
    nn2 = n/2;
    if(n2<nn2) { return;}
    ifault = 1;
    if(n<3) { return;}
     __assume_aligned(a,64);
     __assume_aligned(x,64);
    // If INIT is false, calculates coefficients for the test
    if(!init) {
       if(n==3) {
            a[0] = SQRTH;
        }
        else {
            an25 = an+QTR;
            summ2 = ZERO;
	    //#pragma loop_count min(3),avg(1000),max(5000)
            for(int32_t i = 1; i != n2; ++i) {
                a[i] = ppnd((i-th)/an25);
                summ2 += a[i] * a[i];
            }
            summ2 *= TWO;
            ssumm2 = std::sqrt(summ2);
            rsn = ONE/std::sqrt(an);
            a1 = poly(C1,6,rsn)-a[0]/ssumm2;
            //         Normalize coefficients
            if(n>5) {
                i1 = 3;
                a2 = -a[1]/ssumm2+poly(C2,6,rsn);
                const t0 = a[0]*a[0];
                const t1 = a[1]*a[1];
                fac = std::sqrt((summ2-TWO*t0-TWO*t1)/
                         (ONE-TWO*t0-TWO*t1) );
                a[0] = a1;
                a[1] = a2;
            }
            else {
                i1 = 2;
                fac = std::sqrt((summ2-TWO*a[0]*a[0])/
                                (ONE-TWO*a[1]*a[1]));
                a[0] = a1;
            }
#pragma omp simd aligned(a:64) private(t) \
            linear(i:1) unroll partial(2)
            for(int32_t i = i1; i != nn2; ++i) {
                const t = -a[i];
                a[i] = t/fac;
            }
        }
        init = true;
    }
    if(n1<3) {return;}
    ncens = n-n1;
    ifault = 4;
    const bool b0 = ncens>0 && n<20;
    if(ncens<0 || b0) { return;}
    ifault = ;
    delta = (float)ncens/an;
    if(delta>0.8f) {return;}
    //  If W input as negative, calculate significance level of -W
    if(w<ZERO) {
        w1 = ONE+w;
        ifault = 0;
        goto label_70;
    }
    //  Check for zero range
    ifault = 6;
    range = x[n1]-x[0];
    if(range<SMALL) {return;}
    // Check for correct sort order on range - scaled X
    ifault = 7;
    xx = x[0]/range;
    sx = xx;
    sa = -a[0];
    j = n-1;
//#pragma vector aligned
//#pragma vector always
//#pragma loop_count min(2),avg(1000),max(5000)    
    for(int32_t i = 1; i != n2; ++i) {
        xi = x[i]/range;
        sx += xi;
        if(i!=j) { 
            sa = sa+copysignf((float)1,(float)i-j)
                 *a[std::min(i,j)]
        };
        xx = xi;
        j = j-1;
    }
    ifault = 0;
    if(n>5000) { ifault=2;}
    //
    // Calculate W statistic as squared correlation
    //        between data and coefficients
    //
    sa = sa/n1;
    sx = sx/n1;
    ssa = ZERO;
    ssx = ZERO;
    sax = ZERO;
    j = n;
//#pragma loop_count min(2),avg(1000),max(5000)  
//#pragma vector aligned
//#pragma vector always
//#pragma loop_count min(2),avg(1000),max(5000)
#pragma omp simd aligned(a:64) private(asa,xsx)
    linear(i:1) reduction(+:ssa,ssx,sax) unroll partial(8)
    for(int32_t i = 1; i != n1; ++i) {
        if(i!=j) {
            asa = std::copysignf((float)1,(float)i-j)*
                  a[std::min(i,j)]-sa;
        }
        else {
            asa = -sa;
        }
        xsx = x[i]/range;
        ssa = ssa+asa*asa;
        ssx = ssx+xsx*xsx;
        sax = sax+asa*xsx;
        j = j-1;
    }
    //
    //  W1 equals (1-W) claculated to avoid excessive rounding error
    //       for W very near 1 (a potential problem in very large samples)
    // 
    ssassx = std::sqrtf(ssa*ssx);
    w1 = (ssassx-sax)*(ssassx+sax)/(ssa*ssx);
label_70:
    w = ONE-w1;
    //
    // Calculate significance level for W (exact for N=3)
    //
    if(n==3) {
        pw = PI6*(std::asin(std::sqrt(w))-STQR);
        return;
    }
    y = std::log(w1);
    xx = std::log(an);
    m = ZERO;
    S = ONE;
    if(n<=11) {
        gamma = poly(G,2,an);
        if(y>=gamma) {
            pw=small;
            return;
        }
        y = -std::log(gamma-y);
        m = poly(C3,4,an);
        s = std::exp(poly(C4,4,an));
    }
    else {
        m = poly(C5,4,xx);
        s = std::exp(poly(C6,3,xx));
    }
    if(ncens>0) {
        // Censoring by proportion NCENS/N.  Calculate mean and sd
        //       of normal equivalent deviate of W.
        ld = -log(delta);
        bf = ONE+xx*BF1;
        z90f = Z90+bf*std::pow(poly(C7,2,std::pow(XX90,xx)),ld);
        z95f = Z95+bf*std::pow(poly(C8,2,std::pow(XX95,xx)),ld);
        z99f = z99+bf*std::pow(poly(C9,2,xx),ld);
        //  Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
        // pseudo-mean and pseudo-sd of z as the slope and intercept
        zfm = (z90f+z95f+z99f)/THREE;
        zsd = (Z90*(z90f-zfm)+Z95*(z95f-zfm)+Z99*(z99f-zfm))/zss;
        zbar = zfm-zsd*zm;
        m = m+zbar*s;
        s = s*zsd;
    }
    pw = alnorm((y-m)/s,upper);
} 


__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
double ppnd(double p){
    constexpr double split = 4.19999999999999984456877655248E-1;
    constexpr double A0    = 2.50662823883999985596915394126E0;
    constexpr double A1    = -1.86150006252900013237194798421E1;
    constexpr double A2    =  4.13911977353400004631112096831E1;
    constexpr double A3    = -2.54410604963699995550996391103E1;
    constexpr double B1    = -8.47351093089999984897531248862E0;
    constexpr double B2    = 2.308336743742999885853350861E1;
    constexpr double B3    = -2.1062241018259999947304095258E1;
    constexpr double B4    = 3.13082909833000000432434717368E0;
    constexpr double C0    = -2.78718931138000014513522728521E0;
    constexpr double C1    = -2.2979647913400000902583997231E0;
    constexpr double C2    = 4.85014127135000006063592081773E0;
    constexpr double C3    = 2.32121276858000014087224371906E0;
    constexpr double D1    = 3.54388924762000012691487427219E0;
    constexpr double D2    = 1.63706781897000008818565675028E0;
    constexpr double ZERO  = 0.0;
    constexpr double ONE   = 1.0;
    constexpr double HALF  = 0.5;
    double q,r,result;
    q = p-HALF;
    if(std::fabs(q)>split) { goto label10;}
    // 0.08 < P < 0.92
    r = q*q;
    result = q*(((A3*r+A2)*r+A1)*r+A0)/((((B4*r+B3)*r+B2)*r+B1)*r+ONE);
    return (result);
    // P < 0.08 OR P > 0.92, SET R = MIN(P,1-P)
label10:
    r = p;
    if(q>ZERO) {r = ONE-p;}
    if(r<=ZERO) {goto label20;}
    r = std::sqrt(-std::log(r));
    result = (((C3*r+C2)*r+C1)*r+C0)/((D2*r+D1)*r+ONE);
    if(q<ZERO) {result = -result};
    return (result);
label20:
    result = 0.0;
    return (result);
}

__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
float poly(const float * __restrict c,
           const int32_t nord,
           const float x) {
    float p,result;
    int32_t n2,j;
    result = c[0];
    if(nord==1) { return;}
    p = x*c[nord];
    if(nord==2) {goto label20;}
    n2 = nord-2;
    j = n2+1;
    for(int32_t i = 0; i != n2; ++i) {
        p = (p+c[j])*x;
        j = j-1;
    }
label20:
    result += p;
    return (result;)
}

__attribute__((pure))
__attribute__((hot))
__attribute__((aligned(32)))
__attribute__((always_inline))
static inline
double alnorm(const double x,
              const bool upper) {
    /*
     Evaluates the tail area of the standardised normal curve
c       from x to infinity if upper is .true. or
c       from minus infinity to x if upper is .false.
    */
   constexpr double zero = 0.0e+0;
   constexpr double one  = 1.0e+0;
   constexpr double half = 0.5e+0;
   constexpr double ltone = 7.0e+0;
   constexpr double utzero = 18.66e+0;
   constexpr double con = 1.28e+00;
   constexpr double p = 3.98942280443999985894976134659E-1;
   constexpr double q = 3.99903485040000006289773182289E-1;
   constexpr double r =  3.98942280385000014319984984468E-1;
   constexpr double a1 = 5.75885480458000031944720831234E0;
   constexpr double a2 = 2.62433121678999992099079463514E0;
   constexpr double a3 = 5.92885724437999961367040668847E0;
   constexpr double b1 = -2.98213557806999993715635355329E1;
   constexpr double b2 = 4.86959930692000000362895661965E1;
   constexpr double c1 = -3.80520000000000021871322027772E-8;
   constexpr double c2 = 3.98064793999999973986986256946E-4;
   constexpr double c3 = -1.51679116634999999746469256934E-1;
   constexpr double c4 = 4.83859128080000022720241759089E0;
   constexpr double c5 = 7.4238092402699995542292299433E-1;
   constexpr double c6 = 3.99019417011000010475640920049E0;
   constexpr double d1 = 1.00000615301999995487847172626E0;
   constexpr double d2 = 1.98615381364000009867254448181E0;
   constexpr double d3 = 5.29330324926000006513504558825E0;
   constexpr double d4 = -1.51508972450999994663334291545E1;
   constexpr double d5 = 3.07899330340000005890033207834E1;
   double con,z,y,result;
   bool up;
   z = x;
   if(z>=zero) {goto label_10;}
   up = !up;
   z = -z;
label_10:
   if(z<=ltone || (up && z <= utzero)) {goto label_20;}
   result = 0.0;
   goto label_40;
label_20:
   y = half*z*z;
   if(z>con) {goto label_30;}
   result = half-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
   goto label_40;
label_30:
   result = std::exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
label_40:
   if(!up) {result = result-one;}
   return (result);
}








#endif /* __GMS_DESCRIPTIVE_STATISTICS_HPP__*/
