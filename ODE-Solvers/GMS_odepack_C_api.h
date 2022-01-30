
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



void DLSODE_(void(*F)(int *    __restrict, 
                      double * __restrict, 
                      double * __restrict, 
                      double * __restrict),
            int *              __restrict,
	    double *           __restrict,
	    double *           __restrict,
	    int *              __restrict,
	    double *           __restrict,
	    double *           __restrict,
	    int *              __restrict,
	    int *              __restrict,
	    int *              __restrict,
	    double *           __restrict,
	    int *              __restrict,
	    int *              __restrict,
	    int *              __restrict,
	    void(*JAC)(int *   __restrict , 
                       double * __restrict, 
                       double * __restrict, 
                       int *    __restrict, 
                       int *    __restrict, 
                       double   __restrict*, 
                       int *    __restrict),
	    int *    __restrict);


/*
      SUBROUTINE DLSODES (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void DLSODES_(void(*F)(int *      __restrict, 
                       double *   __restrict, 
                       double *   __restrict, 
                       double *   __restrict),
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
	     void(*JAC)(int *     __restrict, 
                        double *  __restrict, 
                        double *  __restrict, 
                        int *     __restrict, 
                        int *     __restrict, 
                        double *  __restrict, 
                        int *     __restrict),
	     int *      __restrict);


/*

      SUBROUTINE DLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
*/


void DLSODA_( void(*F)(int * __restrict, double * __restrict, double * __restrict, double * __restrict),
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
	     void(*JAC)(int * __restrict, double * __restrict, double * __restrict, 
                        int * __restrict, int * __restrict, double * __restrict, int * __restrict),
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


void DLSODAR_( void(*F)(int * __restrict, double * __restrict, 
                        double * __restrict, double * __restrict),
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
	     void(*JAC)(int *, double * __restrict, double * __restrict, int * __restrict, 
                        int * __restrict, double * __restrict, int * __restrict),
	     int *      __restrict,
	     void(*G)(int * __restrict, double * __restrict, double * __restrict, 
                      int * __restrict, double * __restrict),
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


void DLSODPK_(void(*F)(int * __restrict, double * __restrict, double * __restrict, double * __restrict),
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
	     void(*JAC)(int * __restrict, double * __restrict, double * __restrict, 
                        int * __restrict, int * __restrict, double * __restrict, int * __restrict),
	     void(*PSOL)(int * __restrict, double * __restrict, double * __restrict, double * __restrict, 
                         double * __restrict, double * __restrict, double * __restrict, int * __restrict, 
                         double * __restrict, int * __restrict, int * __restrict),
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


void DLSODKR_(void(*F)(int * __restrict, double * __restrict, double * __restrict, double * __restrict),
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
	     void(*JAC)(int * __restrict, double * __restrict, double * __restrict, 
                        int * __restrict, int * __restrict, double * __restrict, int * __restrict),
	     void(*PSOL)(int * __restrict, double * __restrict, double * __restrict, double * __restrict, 
                         double * __restrict, double * __restrict, double * __restrict, int * __restrict, 
                         double * __restrict, int * __restrict, int * __restrict),
	     int *      __restrict,
	     void(*G)(int * __restrict, double * __restrict, double * __restrict, int * __restrict, double * __restrict),
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


void DLSODI_(void(*RES)(int * __restrict, double * __restrict, double * __restrict, 
                        double * __restrict, double * __restrict, int * __restrict),
            void(*ADDA)(int * __restrict, double * __restrict, double * __restrict, 
                        int * __restrict, int * __restrict, double * __restrict, int * __restrict),
	    void(*JAC)(int * __restrict, double * __restrict, double * __restrict, 
                       double * __restrict, int * __restrict, int * __restrict, double * __restrict, int * __restrict),
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


void DLSOIBT_(void(*RES)(int *, double *, double *, double *, double *, int *),
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


void DLSODIS_(void(*RES)(int *, double *, double *, double *, double *, int *),
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


void SLSODE_(void(*F)(int *, float *, float *, float *),
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


void SLSODES_(void(*F)(int *, float *, float *, float *),
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


void SLSODA_( void(*F)(int *, float *, float *, float *),
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


void SLSODAR_( void(*F)(int *, float *, float *, float *),
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



void SLSODPK_(void(*F)(int *, float *, float *, float *),
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


void SLSODKR_(void(*F)(int *, float *, float *, float *),
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


void SLSODI_(void(*RES)(int *, float *, float *, float *, float *, int *),
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


void SLSOIBT_(void(*RES)(int *, float *, float *, float *, float *, int *),
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



void SLSODIS_(void(*RES)(int *, float *, float *, float *, float *, int *),
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



/*
is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
      CALL DSRCOM(RSAV,ISAV,JOB)  Saves and restores the contents of the
C                               internal COMMON blocks used by DLSODE
C                               (see Part 3 below).  RSAV must be a
C                               real array of length 218 or more, and
C                               ISAV must be an integer array of length
C                               37 or more.  JOB = 1 means save COMMON
C                               into RSAV/ISAV.  JOB = 2 means restore
C                               COMMON from same.  DSRCOM is useful if
C                               one is interrupting a run and restarting
C                               later, or alternating between two or
C                               more problems solved with DLSODE.
*/
void DSRCOM_(double * __restrict, int * __restrict, int * __restrict);

/*
   computes an interpolated value of the y vector at t = TOUT.

   Detailed instructions for using DINTDY
C     --------------------------------------
C     The form of the CALL is:
C
C           CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C     The input parameters are:
C
C     T          Value of independent variable where answers are
C                desired (normally the same as the T last returned by
C                DLSODE).  For valid results, T must lie between
C                TCUR - HU and TCUR.  (See "Optional Outputs" above
C                for TCUR and HU.)
C     K          Integer order of the derivative desired.  K must
C                satisfy 0 <= K <= NQCUR, where NQCUR is the current
C                order (see "Optional Outputs").  The capability
C                corresponding to K = 0, i.e., computing y(t), is
C                already provided by DLSODE directly.  Since
C                NQCUR >= 1, the first derivative dy/dt is always
C                available with DINTDY.
C     RWORK(21)  The base address of the history array YH.
C     NYH        Column length of YH, equal to the initial value of NEQ.
C
C     The output parameters are:
C
C     DKY        Real array of length NEQ containing the computed value
C                of the Kth derivative of y(t).
C     IFLAG      Integer flag, returned as 0 if K and T were legal,
C                -1 if K was illegal, and -2 if T was illegal.
C                On an error return, a message is also written.
*/
void DINTDY_(double * __restrict, int * __restrict, double * __restrict,
             int * __restrict, double * __restrict, int * __restrict);

  










































} // End of extern "C"
#endif


#endif /*__GMS_ODEPACK_F77_H__*/
