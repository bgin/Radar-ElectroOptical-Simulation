!** DCFOD
PURE SUBROUTINE DCFOD(Meth,Elco,Tesco)
       use mod_kinds, only : i4,dp
  !> Subsidiary to DDEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (CFOD-S, DCFOD-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   DCFOD defines coefficients needed in the integrator package DDEBDF
  !
  !***
  ! **See also:**  DDEBDF
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !
  INTEGER(i4), INTENT(IN) :: Meth
  REAL(dp), INTENT(OUT) :: Elco(13,12), Tesco(3,12)
  !
  INTEGER(i4) :: i, ib, nq, nqm1, nqp1
  REAL(dp) :: agamq, fnq, fnqm1, pc(12), pint, ragq, rq1fac, rqfac, tsign, xpin
  !     ------------------------------------------------------------------
  !      DCFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
  !      NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
  !      GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
  !      THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH =
  !      2.  (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
  !      DCFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
  !      AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
  !
  !      THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
  !      THE COEFFICIENTS EL(I), 1 <= I <= NQ+1, FOR THE METHOD OF
  !      ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A
  !      GENERATING POLYNOMIAL, I.E.,
  !          L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
  !      FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
  !          DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) =
  !      0.  FOR THE BDF METHODS, L(X) IS GIVEN BY
  !          L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
  !      WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
  !
  !      THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
  !      LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
  !      AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
  !      SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
  !      NQ + 1 IF K = 3.
  !     ------------------------------------------------------------------

  !* FIRST EXECUTABLE STATEMENT  DCFOD
  IF( Meth==2 ) THEN
    !
    pc(1) = 1._dp
    rq1fac = 1._dp
    DO nq = 1, 5
      !           ------------------------------------------------------------
      !            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
      !                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ).
      !            INITIALLY, P(X) = 1.
      !           ------------------------------------------------------------
      fnq = nq
      nqp1 = nq + 1
      !           FORM COEFFICIENTS OF P(X)*(X+NQ).
      !           ------------------------------------
      pc(nqp1) = 0._dp
      DO ib = 1, nq
        i = nq + 2 - ib
        pc(i) = pc(i-1) + fnq*pc(i)
      END DO
      pc(1) = fnq*pc(1)
      !           STORE COEFFICIENTS IN ELCO AND TESCO.
      !           --------------------------------
      DO i = 1, nqp1
        Elco(i,nq) = pc(i)/pc(2)
      END DO
      Elco(2,nq) = 1._dp
      Tesco(1,nq) = rq1fac
      Tesco(2,nq) = nqp1/Elco(1,nq)
      Tesco(3,nq) = (nq+2)/Elco(1,nq)
      rq1fac = rq1fac/fnq
    END DO
  ELSE
    !
    Elco(1,1) = 1._dp
    Elco(2,1) = 1._dp
    Tesco(1,1) = 0._dp
    Tesco(2,1) = 2._dp
    Tesco(1,2) = 1._dp
    Tesco(3,12) = 0._dp
    pc(1) = 1._dp
    rqfac = 1._dp
    DO nq = 2, 12
      !           ------------------------------------------------------------
      !            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
      !                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ-1).
      !            INITIALLY, P(X) = 1.
      !           ------------------------------------------------------------
      rq1fac = rqfac
      rqfac = rqfac/nq
      nqm1 = nq - 1
      fnqm1 = nqm1
      nqp1 = nq + 1
      !           FORM COEFFICIENTS OF P(X)*(X+NQ-1).
      !           ----------------------------------
      pc(nq) = 0._dp
      DO ib = 1, nqm1
        i = nqp1 - ib
        pc(i) = pc(i-1) + fnqm1*pc(i)
      END DO
      pc(1) = fnqm1*pc(1)
      !           COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X).
      !           -----------------------
      pint = pc(1)
      xpin = pc(1)/2._dp
      tsign = 1._dp
      DO i = 2, nq
        tsign = -tsign
        pint = pint + tsign*pc(i)/i
        xpin = xpin + tsign*pc(i)/(i+1)
      END DO
      !           STORE COEFFICIENTS IN ELCO AND TESCO.
      !           --------------------------------
      Elco(1,nq) = pint*rq1fac
      Elco(2,nq) = 1._dp
      DO i = 2, nq
        Elco(i+1,nq) = rq1fac*pc(i)/i
      END DO
      agamq = rqfac*xpin
      ragq = 1._dp/agamq
      Tesco(2,nq) = ragq
      IF( nq<12 ) Tesco(1,nqp1) = ragq*rqfac/nqp1
      Tesco(3,nqm1) = ragq
    END DO
  END IF
  ! ----------------------- END OF SUBROUTINE DCFOD -----------------------
END SUBROUTINE DCFOD
