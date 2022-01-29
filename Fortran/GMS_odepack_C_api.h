
#ifndef __GMS_ODEPACK_C_API_H__
#define __GMS_ODEPACK_C_API_H__ 290120221022



namespace file_info {


  const unsigned int GMS_ODEPACK_C_API_MAJOR = 1;
  const unsigned int GMS_ODEPACK_C_API_MINOR = 0;
  const unsigned int GMS_ODEPACK_C_API_MICRO = 0;
  const unsigned int GMS_ODEPACK_C_API_FULLVER =
    1000U*GMS_ODEPACK_C_API_MAJOR+100U*GMS_ODEPACK_C_API_MINOR+10U*GMS_ODEPACK_C_API_MICRO;
  const char * const GMS_ODEPACK_C_API_CREATE_DATE = "29-01-2022 10:22 +00200 (SAT 29 JAN 2022 GMT+2)";
  const char * const GMS_ODEPACK_C_API_BUILD_DATE  = __DATE__":"__TIME__;
  const char * const GMS_ODEPACK_C_API_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const GMS_ODEPACK_C_API_SYNOPSIS    = "C-wrappers to ODEPACK Library.";
}








extern "C" {

/*
typedef void(*F)(int *, double *, double *, double *);
typedef void(*JAC)(int *, double *, double *, int *, int *, double *, int *);
typedef void(*PJAC)(int *, double *, double *, int *, double *, double *, double *, double *, int *, F,JAC);
typedef void(*SLVS)(double *, int *, double *, double *);
typedef void(*G)(int *, double *, double *, int *, double *);
typedef void(*PSOL)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int * ); 
*/


/*
      SUBROUTINE DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/



void DLSODE(void(*F)(int *, double *, double *, double *),
            int *    __restrict,
	    double * __restrict,
	    double * __restrict,
	    int *    __restrict,
	    double * __restrict,
	    double * __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    double * __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    void(*JAC)(int *, double *, double *, int *, int *, double *, int *),
	    int *    __restrict);


/*
      SUBROUTINE DLSODES (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void DLSODES(void(*F)(int *, double *, double *, double *),
             int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     void(*JAC)(int *, double *, double *, int *, int *, double *, int *),
	     int *      __restrict);


/*

      SUBROUTINE DLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void DLSODA( void(*F)(int *, double *, double *, double *),
             int *        __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     void(*JAC)(int *, double *, double *, int *, int *, double *, int *),
	     int *      __restrict);


/*

      SUBROUTINE DLSODAR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT,
     2            G, NG, JROOT)
      EXTERNAL F, JAC, G
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT,
     1   NG, JROOT
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1   JROOT(NG)
*/


void DLSODAR( void(*F)(int *, double *, double *, double *),
             int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     void(*JAC)(int *, double *, double *, int *, int *, double *, int *),
	     int *      __restrict,
	     void(*G)(int *, double *, double *, int *, double *),
	     int *      __restrict,
	     int *      __restrict);


/*

      SUBROUTINE DLSODPK (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL, MF)
      EXTERNAL F, JAC, PSOL
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void DLSODPK(void(*F)(int *, double *, double *, double *),
             int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     void(*JAC)(int *, double *, double *, int *, int *, double *, int *),
	     void(*PSOL)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int * ),
	     int *      __restrict);


/*

     SUBROUTINE DLSODKR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL,
     2            MF, G, NG, JROOT)
      EXTERNAL F, JAC, PSOL, G
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF,
     1        NG, JROOT
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          JROOT(*)
*/


void DLSODKR(void(*F)(int *, double *, double *, double *),
             int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     double *   __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     int *      __restrict,
	     void(*JAC)(int *, double *, double *, int *, int *, double *, int *),
	     void(*PSOL)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int * ),
	     int *      __restrict,
	     void(*G)(int *, double *, double *, int *, double *),
	     int *      __restrict,
	     int *      __restrict);


/*
     
       SUBROUTINE DLSODI (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
*/


void DLSODI(void(*RES)(int *, double *, double *, double *, double *, int *),
            void(*ADDA)(int *, double *, double *, int *, int *, double *, int *),
	    void(*JAC)(int *, double *, double *, double *, int *, int *, double *, int *),
	    int *     __restrict,
	    double *  __restrict,
	    double *  __restrict,
	    double *  __restrict,
	    double *  __restrict,
	    int *     __restrict,
	    double *  __restrict,
	    double *  __restrict,
	    int *     __restrict,
	    int *     __restrict,
	    int *     __restrict,
	    double *  __restrict,
	    int *     __restrict,
	    int *     __restrict,
	    int *     __restrict,
	    int *     __restrict);



/*
 
      SUBROUTINE DLSOIBT (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
*/


void DLSOIBT(void(*RES)(int *, double *, double *, double *, double *, int *),
             void(*ADDA)(int *, double *, double *, int *, int *, double *, double *, double *),
	     void(*JAC)(int *, double *, double *, double *, int *, int *, double *, double *, double *),
	     int *     __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     int *     __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     double *  __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict);



/*

      
      SUBROUTINE DLSODIS (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
*/


void DLSODIS(void(*RES)(int *, double *, double *, double *, double *, int *),
             void(*ADDA)(int *, double *, double *, int *, int *, int *, double *),
	     void(*JAC)(int *, double *, double *, double *, int *, int *, int *, double *),
	     int *     __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     int *     __restrict,
	     double *  __restrict,
	     double *  __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     double *  __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict);



/*
       
           ODEPACK Single-Precision Interfaces

*/



/*

      SUBROUTINE SLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void SLSODE(void(*F)(int *, float *, float *, float *),
            int *   __restrict,
	    float * __restrict,
	    float * __restrict,
	    int *   __restrict,
	    float * __restrict,
	    float * __restrict,
	    int *   __restrict,
	    int *   __restrict,
	    int *   __restrict,
	    float * __restrict,
	    int *   __restrict,
	    int *   __restrict,
	    int *   __restrict,
	    void(*JAC)(int *, float *, float *, int *, int *, float *, int *),
	    int *   __restrict);


/*

      SUBROUTINE SLSODES (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void SLSODES(void(*F)(int *, float *, float *, float *),
             int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *      __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     void(*JAC)(int *, float *, float *, int *, int *, float *, int *),
	     int *     __restrict);



/*

         SUBROUTINE SLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
         EXTERNAL F, JAC
         INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT
         REAL Y, T, TOUT, RTOL, ATOL, RWORK
         DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void SLSODA( void(*F)(int *, float *, float *, float *),
             int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     void(*JAC)(int *, float *, float *, int *, int *, float *, int *),
	     int *     __restrict);


/*
  
      SUBROUTINE SLSODAR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT,
     2            G, NG, JROOT)
      EXTERNAL F, JAC, G
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT,
     1   NG, JROOT
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1   JROOT(NG)
*/


void SLSODAR( void(*F)(int *, float *, float *, float *),
             int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     void(*JAC)(int *, float *, float *, int *, int *, float *, int *),
	     int *     __restrict,
	     void(*G)(int *, float *, float *, int *, float *),
	     int *     __restrict,
	     int *     __restrict);



/*

      SUBROUTINE SLSODPK (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL, MF)
      EXTERNAL F, JAC, PSOL
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/



void SLSODPK(void(*F)(int *, float *, float *, float *),
             int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     void(*JAC)(int *, float *, float *, int *, int *, float *, int *),
	     void(*PSOL)(int *, float *, float *, float *, float *, float *, float *, int *, float *, int *, int * ),
	     int *     __restrict);



/*

     SUBROUTINE SLSODKR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL,
     2            MF, G, NG, JROOT)
      EXTERNAL F, JAC, PSOL, G
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF,
     1        NG, JROOT
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          JROOT(*)
*/


void SLSODKR(void(*F)(int *, float *, float *, float *),
             int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     float *   __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     int *     __restrict,
	     void(*JAC)(int *, float *, float *, int *, int *, float *, int *),
	     void(*PSOL)(int *, float *, float *, float *, float *, float *, float *, int *, float *, int *, int * ),
	     int *     __restrict,
	     void(*G)(int *, float *, float *, int *, float *),
	     int *     __restrict,
	     int *     __restrict);
	     

/*

      SUBROUTINE SLSODI (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
*/


void SLSODI(void(*RES)(int *, float *, float *, float *, float *, int *),
            void(*ADDA)(int *, float *, float *, int *, int *, float *, int *),
	    void(*JAC)(int *, float *, float *, float *, int *, int *, float *, int *),
	    int *    __restrict,
	    float *  __restrict,
	    float *  __restrict,
	    float *  __restrict,
	    float *  __restrict,
	    int *    __restrict,
	    float *  __restrict,
	    float *  __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    float *  __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    int *    __restrict,
	    int *    __restrict);



/*

       SUBROUTINE SLSOIBT (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
*/


void SLSOIBT(void(*RES)(int *, float *, float *, float *, float *, int *),
             void(*ADDA)(int *, float *, float *, int *, int *, float *, float *, float *),
	     void(*JAC)(int *, float *, float *, float *, int *, int *, float *, float *, float *),
	     int *    __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     int *    __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     float *  __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     int *    __restrict);



/*

       SUBROUTINE SLSODIS (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
*/



void SLSODIS(void(*RES)(int *, float *, float *, float *, float *, int *),
             void(*ADDA)(int *, float *, float *, int *, int *, int *, float *),
	     void(*JAC)(int *, float *, float *, float *, int *, int *, int *, float *),
	     int *    __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     int *    __restrict,
	     float *  __restrict,
	     float *  __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     float *  __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     int *    __restrict,
	     int *    __restrict);


  










































} // End of extern "C"
#endif


#endif /*__GMS_ODEPACK_F77_H__*/
