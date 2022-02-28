

#ifndef __GMS_CULSODA_ODE_CUH__
#define __GMS_CULSODA_ODE_CUH__

#include <stdint.h>

/*

     List of test ODE for CUDA lsoda integrator.
*/


struct ODE_1 {

    __device__ void operator()(int32_t *  __restrict__ neq,
                               double  *  __restrict__ t,
                               double  *  __restrict__ y,
                               double  *  __restrict__ ydot) {
        
            ydot[0] = 1.0E4*y[1]*y[2]-0.04E0*y[0];
            ydot[2] = 3.0E7*y[1]*y[1];
            ydot[1] = 1.0*(ydot[0]+ydot[2]);
    }                 
};

/*

! DEMONSTRATION PROGRAM FOR THE DVODE_F90 PACKAGE.

! The following is a simple example problem, with the coding
! needed for its solution by DVODE_F90. The problem is from
! chemical kinetics, and consists of the following three rate
! equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
*/
struct ODE_2 {
 
   __device__ void operator()( int32_t *  __restrict__ neq,
                               double  *  __restrict__ t,
                               double  *  __restrict__ y,
                               double  *  __restrict__ ydot) {
    
          ydot[0] = -.04E0*y[0] + 1.0E4*y[1]*y[2];
          ydot[2] = 3.E7*y[1]*y[1];
          ydot[1] = -ydot[0]-ydot[2];
   } 

};


struct JAC_2 {

    __device__ void operator()(int32_t * __restrict__ neq,
                               double  * __restrict__ t,
                               double  * __restrict__ y,
                               int32_t * __restrict__ ml,
                               int32_t * __restrict__ mu,
                               double  * __restrict__ pd,
                               int32_t * __restrict__ nrpd) {
         pd[0] = -0.04;
         pd[1] = 1.0E4*y[2];
         pd[2] = 1.0E4*y[1];
         pd[3] = 0.04;
         pd[5] = -pd[2];
         pd[7] = 6.0E7*y[1];
         pd[4] = -pd[1]-pd[7];
    }
};


/*
       SUBROUTINE F2(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(2), YDOT(2)

        YDOT(1) = Y(2)
        YDOT(2) = 100.0D0*(1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
        RETURN
      END SUBROUTINE F2

      SUBROUTINE JAC2(NEQ,T,Y,ML,MU,PD,NROWPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(2), PD(NROWPD,2)

        PD(1,1) = 0.0D0
        PD(1,2) = 1.0D0
        PD(2,1) = -200.0D0*Y(1)*Y(2) - 1.0D0
        PD(2,2) = 100.0D0*(1.0D0-Y(1)*Y(1))
        RETURN
      END SUBROUTINE JAC2
*/


struct ODE_3 {

     __device__ void operator()(int32_t *  __restrict__ neq,
                                double  *  __restrict__ t,
                                double  *  __restrict__ y,
                                double  *  __restrict__ ydot) {
         
           ydot[0] = ydot[1];
           ydot[1] = 100.0E0*(1.0E0-y[0]*y[0])*y[1]-y[0];

     }
};


struct JAC_3 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
          
            pd[0] = 0.0E0;
            pd[1] = 1.0E0;
            pd[2] = -200.0E0*y[0]*y[1]-1.0E0;
            pd[3] = 100.0E0*(1.0E0-y[0]*y[1]);
      }
};


/*
      SUBROUTINE F1(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(1), YDOT(1)

        YDOT(1) = ((2.0D0*LOG(Y(1))+8.0D0)/T-5.0D0)*Y(1)
        RETURN
      END SUBROUTINE F1
*/


struct ODE_3 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
           
           ydot[0] = ((2.0E0*__logf(y[0])+8.0E0)/t-5.0E0)*y[0];
      }
};


/*
      SUBROUTINE F1(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(2), YDOT(2)

        YDOT(1) = Y(2)
        YDOT(2) = 3.0D0*(1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
        RETURN
      END SUBROUTINE F1

      SUBROUTINE JAC1(NEQ,T,Y,ML,MU,PD,NROWPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(2), PD(NROWPD,2)

        PD(1,1) = 0.0D0
        PD(1,2) = 1.0D0
        PD(2,1) = -6.0D0*Y(1)*Y(2) - 1.0D0
        PD(2,2) = 3.0D0*(1.0D0-Y(1)*Y(1))
        RETURN
      END SUBROUTINE JAC1
*/


struct ODE_4 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
          
           ydot[0] = y[1];
           ydot[1] = 3.0E0*(1.0E0-y[0]*y[0])*y[1]-y[0];
     }
};


struct JAC_4 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
          
           pd[0] = 0.0E0;
           pd[1] = 1.0E0;
           pd[2] = -6.0E0*y[0]*y[1]-1.0E0;
           pd[3] = 3.0E0*(1.0E0-y[0]*y[0]);
     }
};





#endif /*__GMS_CULSODA_ODE_CUH__*/
