#ifndef __GMS_FFTPACK51_DECLARATIONS_H__
#define __GMS_FFTPACK51_DECLARATIONS_H__



/*
 extern "C"  Wrapper declarations of corresponding F77 FFTPack functions
@author: Bernard Gingold
@version:  1.0  26/10/2015

*/

#if defined __cplusplus

extern "C"
{
	
	/*
	@brief SUBROUTINE RFFT1I (N, WSAVE, LENSAV, IER)
     INTEGER    N, LENSAV, IER
     REAL       WSAVE(LENSAV)
	*/
	void   RFFT1I(int *, float *, int *, int *);

	/*
	@brief  SUBROUTINE RFFT1B (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

           INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
           REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void    RFFT1B(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief SUBROUTINE RFFT1F (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

           INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
           REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void    RFFT1F(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	brief  SUBROUTINE RFFT2I (L, M, WSAVE, LENSAV, IER)
    INTEGER    L, M, LENSAV, IER
    REAL       WSAVE(LENSAV)
	*/
	void    RFFT2I(int *, int *, float *, int *);

	/*
	@brief SUBROUTINE RFFT2B (LDIM, L, M, R, WSAVE, LENSAV, WORK, LENWRK, IER)
    INTEGER    LDIM, L, M, LENSAV, LENWRK, IER
    REAL       R(LDIM,M), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void    RFFT2B(int *, int *, int *, float *, float *, int *, float *, int *, int *);

	/*
	@brief SUBROUTINE RFFT2F (LDIM, L, M, R, WSAVE, LENSAV, WORK, LENWRK, IER)
     INTEGER    LDIM, L, M,  LENSAV, LENWRK, IER
     REAL       R(LDIM,M), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void    RFFT2F(int *, int *, int *, float *, float *, int *, float *, int *, int *);

	/*
	@brief SUBROUTINE RFFTMI (N, WSAVE, LENSAV, IER)
    INTEGER    N, LENSAV, IER
    REAL       WSAVE(LENSAV)
	*/
	void    RFFTMI(int *, float *, int *, int *);

	/*
	@brief  SUBROUTINE RFFTMB (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV,
  1                  WORK, LENWRK, IER)

   INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
   REAL       R(LENR), WSAVE(LENSAV)     ,WORK(LENWRK)
	*/
	void    RFFTMB(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief  SUBROUTINE RFFTMF (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV,
         1                  WORK, LENWRK, IER)

         INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
         REAL       R(LENR), WSAVE(LENSAV)     ,WORK(LENWRK)
	*/
	void    RFFTMF(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief SUBROUTINE COST1I (N, WSAVE, LENSAV, IER)
     INTEGER    N, LENSAV, IER
     REAL       WSAVE(LENSAV)
	*/
	void    COSTI(int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COST1B (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

             INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void     COST1B(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief   SUBROUTINE COST1F (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)
             INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK) 
	*/
	void     COST1F(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COSTMI (N, WSAVE, LENSAV, IER)
              INTEGER    N, LENSAV, IER
              REAL       WSAVE(LENSAV)
	*/
	void      COSTMI(int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE COSTMB (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
               1                   WORK, LENWRK, IER)

             INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void       COSTMB(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE COSTMF (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
               1                   WORK, LENWRK, IER)

               INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
               REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)

	*/
	void       COSTMF(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE SINT1I (N, WSAVE, LENSAV, IER)
     INTEGER    N, LENSAV, IER
     REAL       WSAVE(LENSAV)
	*/
	void        SINT1I(int *, float *, int *, int *);

	/*
	@brief      SUBROUTINE SINT1B (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

                INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
                REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void        SINT1B(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE SINT1F (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

               INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
               REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void        SINT1F(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE SINTMI (N, WSAVE, LENSAV, IER)
               INTEGER    N, LENSAV, IER
               REAL       WSAVE(LENSAV)
	*/
	void       SINTM1(int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE SINTMB (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
    1                   WORK, LENWRK, IER)

             INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      SINTMB(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE SINTMF (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
1                   WORK, LENWRK, IER)

             INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      SINTMF(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COSQ1I (N, WSAVE, LENSAV, IER)
    INTEGER    N, LENSAV, IER
    REAL       WSAVE(LENSAV)
	*/
	void      COSQ1I(int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COSQ1B (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

              INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
              REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      COSQ1B(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COSQ1F (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

             INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      COSQ1F(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COSQMI (N, WSAVE, LENSAV, IER)
              INTEGER    N, LENSAV, IER
              REAL       WSAVE(LENSAV)
	*/
	void      COSQMI(int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE COSQMB (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
                   WORK, LENWRK, IER)

             INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      COSQMB(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE COSQMF (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
                   WORK, LENWRK, IER)

               INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
               REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      COSQMF(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief     SUBROUTINE SINQ1I (N, WSAVE, LENSAV, IER)
               INTEGER    N, LENSAV, IER
               REAL       WSAVE(LENSAV)
	*/
	void      SINQ1I(int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE SINQ1B (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

    INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
    REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void       SINQ1B(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief    SUBROUTINE SINQ1F (N, INC, R, LENR, WSAVE, LENSAV, WORK, LENWRK, IER)

             INTEGER    N, INC, LENR, LENSAV, LENWRK, IER
             REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void      SINQ1F(int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	brief    SUBROUTINE SINQMI (N, WSAVE, LENSAV, IER)
             INTEGER    N, LENSAV, IER
             REAL       WSAVE(LENSAV)
	*/
	void     SINQMI(int *, float *, int *, int *);

	/*
	@brief   SUBROUTINE SINQMB (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
1                   WORK, LENWRK, IER)

            INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
            REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void     SINQMB(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	/*
	@brief  SUBROUTINE SINQMF (LOT, JUMP, N, INC, R, LENR, WSAVE, LENSAV, 
1                   WORK, LENWRK, IER)

            INTEGER    LOT, JUMP, N, INC, LENR, LENSAV, LENWRK, IER
            REAL       R(LENR), WSAVE(LENSAV), WORK(LENWRK)
	*/
	void     SINQMF(int *, int *, int *, int *, float *, int *, float *, int *, float *, int *, int *);

	
}
      
#endif

#endif  /*__GMS_FFTPACK51_DECLARATIONS_H__*/
