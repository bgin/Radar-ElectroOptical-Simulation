

module mol_solvers


!============================================================
! This module contains various subroutines adapted from
! the codes of William E. Schiesser.
! Reference:
!               https://www.lehigh.edu/~wes1/books/
! This is rather a low priority module for solving
! PDE at least in case of elliptic types.
!===========================================================

  
  public
  implicit none

  !================================================================!
  !           Finite-Difference Scheme collection                  !
  !================================================================!

contains
  
      SUBROUTINE DSS002(XL,XU,N,U,UX)
          !dir$ optimize:3
          !dir$ attributes forceinline :: DSS002
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS002
#if 0
C...
C...  SUBROUTINE DSS002 COMPUTES THE FIRST DERIVATIVE, U , OF A
C...                                                    X
C...  VARIABLE U OVER THE SPATIAL DOMAIN XL LE X LE XU
C...
C...  ARGUMENT LIST
C...
C...     XL      LOWER BOUNDARY VALUE OF X (INPUT)
C...
C...     XU      UPPER BOUNDARY VALUE OF X (INPUT)
C...
C...     N       NUMBER OF GRID POINTS IN THE X DOMAIN INCLUDING THE
C...             BOUNDARY POINTS (INPUT)
C...
C...     U       ONE-DIMENSIONAL ARRAY CONTAINING THE VALUES OF U AT
C...             THE N GRID POINT POINTS FOR WHICH THE DERIVATIVE IS
C...             TO BE COMPUTED (INPUT)
C...
C...     UX      ONE-DIMENSIONAL ARRAY CONTAINING THE NUMERICAL
C...             VALUES OF THE DERIVATIVES OF U AT THE N GRID POINTS
C...             (OUTPUT)
C...
C...  SUBROUTINE DSS002 COMPUTES THE FIRST DERIVATIVE, U , OF A
C...                                                    X
C...  VARIABLE U OVER THE SPATIAL DOMAIN XL LE X LE XU FROM THE
C...  CLASSICAL THREE-POINT, SECOND-ORDER FINITE DIFFERENCE APPROXI-
C...  TIONS
C...
C...                                       2
C...  U1  = (1/2DX)(-3U1 + 4U2 - U3) + O(DX ) (LEFT BOUNDARY,     (1)
C...    X                                         X = XL)
C...
C...                                   2
C...  UI  = (1/2DX)(UI+1 - UI-1) + O(DX ) (INTERIOR POINT,        (2)
C...    X                                   X NE XL, XU)
C...
C...                                          2
C...  UN  = (1/2DX)(3UN - 4UN-1 + UN-2) + O(DX ) (RIGHT BOUNDARY, (3)
C...    X                                            X = XU)
C...
C...  EQUATIONS (1) TO (3) APPLY OVER A GRID IN X WITH CORRESPONDING
C...  VALUES OF THE FUNCTION U(X) REPRESENTED AS
C...
C...   U1      U2       U3         UI        UN-2      UN-1    UN
C...
C...  X=XL  X=XL+DX  X=XL+2DX ... X=XI ... X=XU-2DX  X=XU-DX  X=XU
C...
C...  THE ORIGIN OF EQUATIONS (1) TO (3) IS OUTLINED BELOW.
C...
C...  CONSIDER THE FOLLOWING POLYNOMIAL IN X OF ARBITRARY ORDER
C...
C...                                     2             3
C...  U(X) = A0 + A1(X - X0) + A2(X - X0)  + A3(X - X0)  + ....   (4)
C...
C...  WE SEEK THE VALUES OF THE COEFFICIENTS A0, A1, A2, ... FOR A
C...  PARTICULAR FUNCTION U(X).  IF X = X0 IS SUBSTITUTED IN EQUATION
C...  (4), WE HAVE IMMEDIATELY A0 = U(X0).  NEXT, IF EQUATION (4) IS
C...  DIFFERENTIATED WITH RESPECT TO X,
C...
C...                                                   2
C...  DU(X)/DX = U (X) = A1 + 2A2(X - X0) + 3A3(X - X0)  + ...    (5)
C...              X
C...
C...  AGAIN, WITH X = X0, A1 = DU(X0)/DX = U (X0).  DIFFERENTIATION
C...                                        X
C...  OF EQUATION (5) IN TURN GIVES
C...
C...  D2U(X)/DX2 = U  (X) = 2A2 + 6A3(X - X0) + ...
C...                2X
C...
C...  AND FOR X = X0, A2 = U  (X0)/2F (2F = 1*2, I.E., 2 FACTORIAL).
C...                        2X
C...
C...  WE CAN CONTINUE THIS PROCESS OF DIFFERENTIATION FOLLOWED BY THE
C...  SUBSTITUTION X = X0 TO OBTAIN THE SUCCESSIVE COEFFICIENTS IN
C...  EQUATION (4), A3, A4, ...  FINALLY, SUBSTITUTION OF THESE CO-
C...  EFFICIENTS IN EQUATION (4) GIVES
C...
C...                                                 2
C...  U(X) = U(X0) + U (X0)(X - X0) + U  (X0)(X - X0)  +
C...                  X       1F       2X       2F
C...                                                              (6)
C...                                3                  4
C...                 U  (X0)(X - X0)  + U  (X0)(X - X0)  + ...
C...                  3X       3F        4X       4F
C...
C...  THE CORRESPONDENCE BETWEEN EQUATION (6) AND THE WELL-KNOWN
C...  TAYLOR SERIES SHOULD BE CLEAR.  THUS THE EXPANSION OF A
C...  FUNCTION, U(X), AROUND A NEIGHBORING POINT X0 IN TERMS OF U(X0)
C...  AND THE DERIVATIVES OF U(X) AT X = X0 IS EQUIVALENT TO APPROXI-
C...  MATING U(X) NEAR X0 BY A POLYNOMIAL.
C...
C...  EQUATION (6) IS THE STARTING POINT FOR THE DERIVATION OF THE
C...  CLASSICAL FINITE DIFFERENCE APPROXIMATIONS OF DERIVATIVES SUCH
C...  AS THE THREE-POINT FORMULAS OF EQUATIONS (1), (2) AND (3).  WE
C...  WILL NOW CONSIDER THE DERIVATION OF THESE THREE-POINT FORMULAS
C...  IN A STANDARD FORMAT WHICH CAN THEN BE EXTENDED TO HIGHER
C...  MULTI-POINT FORMULAS IN OTHER SUBROUTINES, E.G., FIVE-POINT
C...  FORMULAS IN SUBROUTINE DSS004.
C...
C...  THREE-POINT FORMULAS
C...
C...     (1)  LEFT END, POINT I = 1
C...
C...  IF EQUATION (6) IS WRITTEN AROUND THE POINT X = XL FOR X = XL +
C...  DX AND X = XL + 2DX, FOR WHICH THE CORRESPONDING VALUES OF U(X)
C...  ARE U1, U2 AND U3 (U1 AND U2 ARE SEPARATED WITH RESPECT TO X BY
C...  DISTANCE DX AS ARE U2 AND U3, I.E., WE ASSUME A UNIFORM GRID
C...  SPACING, DX, FOR INDEPENDENT VARIABLE X)
C...
C...                                2            3
C...  U2 = U1 + U1 ( DX) + U1  ( DX)  + U1  ( DX)  + ...          (7)
C...              X  1F      2X  2F       3X  3F
C...
C...                                2            3
C...  U3 = U1 + U1 (2DX) + U1  (2DX)  + U1  (2DX)  + ...          (8)
C...              X  1F      2X  2F       3X  3F
C...
C...  WE CAN NOW TAKE A LINEAR COMBINATION OF EQUATIONS (7) AND (8)
C...  BY FIRST MULTIPLYING EQUATION (7) BY A CONSTANT, A, AND EQUA-
C...  TION (8) BY CONSTANT B
C...
C...                                  2           3
C...  A(U2 = U1 + U1 ( DX) + U1  ( DX) + U1  ( DX) + ...)         (9)
C...                X  1F      2X  2F      3X  3F
C...
C...                                  2           3
C...  B(U3 = U1 + U1 (2DX) + U1  (2DX) + U1  (2DX) + ...)        (10)
C...                X  1F      2X  2F      3X  3F
C...
C...  CONSTANTS A AND B ARE THEN SELECTED SO THAT THE COEFFICIENTS OF
C...  THE U1  TERMS SUM TO ONE (SINCE WE ARE INTERESTED IN OBTAINING
C...        X
C...  A FINITE DIFFERENCE APPROXIMATION FOR THIS FIRST DERIVATIVE).
C...  ALSO, WE SELECT A AND B SO THAT THE COEFFICIENTS OF THE U1
C...                                                            2X
C...  TERMS SUM TO ZERO IN ORDER TO DROP OUT THE CONTRIBUTION OF THIS
C...  SECOND DERIVATIVE (THE BASIC IDEA IS TO DROP OUT AS MANY OF THE
C...  DERIVATIVES AS POSSIBLE IN THE TAYLOR SERIES BEYOND THE DERI-
C...  VATIVE OF INTEREST, IN THIS CASE U1 , IN ORDER TO PRODUCE A
C...                                     X
C...  FINITE DIFFERENCE APPROXIMATION FOR THE DERIVATIVE OF MAXIMUM
C...  ACCURACY).  IN THIS CASE WE HAVE ONLY TWO CONSTANTS, A AND B,
C...  TO SELECT SO WE CAN DROP OUT ONLY THE SECOND DERIVATIVE, U1  ,
C...                                                             2X
C...  IN THE TAYLOR SERIES (IN ADDITION TO RETAINING THE FIRST DERI-
C...  VATIVE).  THIS PROCEDURE LEADS TO TWO LINEAR ALGEBRAIC EQUA-
C...  TIONS IN THE TWO CONSTANTS
C...
C...  A + 2B = 1
C...
C...  A + 4B = 0
C...
C...  SOLUTION OF THESE EQUATIONS FOR A AND B GIVES
C...
C...  A = 2, B = -1/2
C...
C...  SOLUTION OF EQUATIONS (9) AND (10) FOR U1  WITH THESE VALUES OF
C...  A AND B GIVES EQUATION (1)               X
C...
C...                                      2
C...  U1 = (1/2DX)(-3U1 + 4U2 - U3) + O(DX )                      (1)
C...    X
C...               2
C...  THE TERM O(DX ) INDICATES A PRINCIPAL ERROR TERM DUE TO TRUNCA-
C...                                                2
C...  TION OF THE TAYLOR SERIES WHICH IS OF ORDER DX .  THIS TERM IN
C...                    2
C...  FACT EQUALS U1  DX /3F, WHICH IS EASILY OBTAINED IN DERIVING
C...                3X
C...  EQUATION (1).
C...
C...  THIS SAME BASIC PROCEDURE CAN NOW BE APPLIED TO THE DERIVATION
C...  OF EQUATIONS (2) AND (3).
C...
C...     (2)  INTERIOR POINT I
C...
C...                                    2           3
C...  A(UI-1 = UI + UI (-DX) + UI  (-DX) + UI  (-DX) + ...)
C...                  X  1F      2X  2F      3X  3F
C...
C...                                    2           3
C...  B(UI+1 = UI + UI ( DX) + UI  ( DX) + UI  ( DX) + ...)
C...                  X  1F      2X  2F      3X  3F
C...
C...  -A + B = 1
C...
C...   A + B = 0
C...
C...  A = -1/2, B = 1/2
C...                                   2
C...  UI  = (1/2DX)(UI+1 - UI-1) + O(DX )                         (2)
C...    X
C...
C...     (3)  RIGHT END, POINT I = N
C...
C...                                      2            3
C...  A(UN-2 = UN + UN (-2DX) + UN  (-2DX) + UN  (-2DX) + ...)
C...                  X   1F      2X   2F      3X   3F
C...
C...                                      2            3
C...  B(UN-1 = UN + UN ( -DX) + UN  ( -DX) + UN  ( -DX) + ...)
C...                  X   1F      2X   2F      3X   3F
C...
C...  -2A - B = 1
C...
C...   4A + B = 0
C...
C...   A = -2, B = 1/2
C...                                          2
C...  UN  = (1/2DX)(3UN - 4UN-1 + UN-2) + O(DX )                  (3)
C...    X
C...
C...  THE WEIGHTING COEFFICIENTS FOR EQUATIONS (1), (2) AND (3) CAN
C...  BE SUMMARIZED AS
C...
C...          -3   4  -1
C...
C...     1/2  -1   0   1
C...
C...           1  -4   3
C...
C...  WHICH ARE THE COEFFICIENTS REPORTED BY BICKLEY FOR N = 2, M =
C...  1, P = 0, 1, 2 (BICKLEY, W. G., FORMULAE FOR NUMERICAL DIFFER-
C...  ENTIATION, MATH. GAZ., VOL. 25, 1941).
C...
C...  EQUATIONS (1), (2) AND (3) CAN NOW BE PROGRAMMED TO GENERATE
C...  THE DERIVATIVE U (X) OF FUNCTION U(X) (ARGUMENTS U AND UX OF
C...                  X
C...  SUBROUTINE DSS002, RESPECTIVELY).
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION     DX,  R2FDX,      U,     UX,     XL,     XU
       !               DFLOAT
!C...
      DIMENSION U(N),UX(N)
!C...
!C...  COMPUTE THE SPATIAL INCREMENT
      DX=(XU-XL)/DFLOAT(N-1)
      R2FDX=1.D+00/(2.D+00*DX)
      NM1=N-1
!C...
!C...  EQUATION (1) (NOTE - THE RHS OF THE FINITE DIFFERENCE APPROXI-
!C...  TIONS, EQUATIONS (1), (2) AND (3) HAVE BEEN FORMATTED SO THAT
!C...  THE NUMERICAL WEIGHTING COEFFICIENTS CAN BE MORE EASILY ASSOCI-
!C...  ATED WITH THE BICKLEY MATRIX LISTED ABOVE)
      UX(1)=R2FDX*(-3.D+00*U(1)+4.D+00*U(2)-1.D+00*U(3)) 
      
!C...
!C...  EQUATION (2)
      !dir$ assume_aligned UX:64,U:64
      !dir$ vector aligned
      !dir$ vector vectorlength(8)
      !dir$ ivdep
      !dir$ vector always
      DO  I=2,NM1
          UX(I)=R2FDX*(-1.D+00*U(I-1)+0.D+00*U(I)+1.D+00*U(I+1))  
      END DO
!C...
!C...  EQUATION (3)
      UX(N)=R2FDX*(1.D+00*U(N-2)-4.D+00*U(N-1)+3.D+00*U(  N))
                
      
    END SUBROUTINE


    
   SUBROUTINE DSS004(XL,XU,N,U,UX)
          !dir$ optimize:3
          !dir$ attributes forceinline :: DSS004
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS004
#if 0
C...
C...  SUBROUTINE DSS004 COMPUTES THE FIRST DERIVATIVE, U , OF A
C...                                                    X
C...  VARIABLE U OVER THE SPATIAL DOMAIN XL LE X LE XU FROM CLASSICAL
C...  FIVE-POINT, FOURTH-ORDER FINITE DIFFERENCE APPROXIMATIONS
C...
C...  ARGUMENT LIST
C...
C...     XL      LOWER BOUNDARY VALUE OF X (INPUT)
C...
C...     XU      UPPER BOUNDARY VALUE OF X (INPUT)
C...
C...     N       NUMBER OF GRID POINTS IN THE X DOMAIN INCLUDING THE
C...             BOUNDARY POINTS (INPUT)
C...
C...     U       ONE-DIMENSIONAL ARRAY CONTAINING THE VALUES OF U AT
C...             THE N GRID POINT POINTS FOR WHICH THE DERIVATIVE IS
C...             TO BE COMPUTED (INPUT)
C...
C...     UX      ONE-DIMENSIONAL ARRAY CONTAINING THE NUMERICAL
C...             VALUES OF THE DERIVATIVES OF U AT THE N GRID POINTS
C...             (OUTPUT)
C...
C...  THE MATHEMATICAL DETAILS OF THE FOLLOWING TAYLOR SERIES (OR
C...  POLYNOMIALS) ARE GIVEN IN SUBROUTINE DSS002.
C...
C...  FIVE-POINT FORMULAS
C...
C...     (1)  LEFT END, POINT I = 1
C...
C...                                   2            3            4
C...  A(U2 = U1 + U1  ( DX) + U1  ( DX)  + U1  ( DX)  + U1  ( DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U1  ( DX)  + U1  ( DX)  + U1  ( DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...                                   2            3            4
C...  B(U3 = U1 + U1  (2DX) + U1  (2DX)  + U1  (2DX)  + U1  (2DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U1  (2DX)  + U1  (2DX)  + U1  (2DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...                                   2            3            4
C...  C(U4 = U1 + U1  (3DX) + U1  (3DX)  + U1  (3DX)  + U1  (3DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U1  (3DX)  + U1  (3DX)  + U1  (3DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...                                   2            3            4
C...  D(U5 = U1 + U1  (4DX) + U1  (4DX)  + U1  (4DX)  + U1  (4DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U1  (4DX)  + U1  (4DX)  + U1  (4DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...  CONSTANTS A, B, C AND D ARE SELECTED SO THAT THE COEFFICIENTS
C...  OF THE U1  TERMS SUM TO ONE AND THE COEFFICIENTS OF THE U1  ,
C...           X                                                2X
C...  U1   AND U1   TERMS SUM TO ZERO
C...    3X       4X
C...
C...  A +   2B +   3C +   4D = 1
C...
C...  A +   4B +   9C +  16D = 0
C...
C...  A +   8B +  27C +  64D = 0
C...
C...  A +  16B +  81C + 256D = 0
C...
C...  SIMULTANEOUS SOLUTION FOR A, B, C AND D FOLLOWED BY THE SOLU-
C...  TION OF THE PRECEDING TAYLOR SERIES, TRUNCATED AFTER THE U
C...                                                            4X
C...  TERMS, FOR U1  GIVES THE FOLLOWING FIVE-POINT APPROXIMATION
C...               X
C...                                                         4
C...  U1  = (1/12DX)(-25U1 + 48U2 - 36U3 + 16U4 - 3U5) + O(DX )   (1)
C...    X
C...
C...     (2)  INTERIOR POINT, I = 2
C...
C...                                   2            3            4
C...  A(U1 = U2 + U2  (-DX) + U2  (-DX)  + U2  (-DX)  + U2  (-DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U2  (-DX)  + U2  (-DX)  + U2  (-DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...                                   2            3            4
C...  B(U3 = U2 + U2  ( DX) + U2  ( DX)  + U2  ( DX)  + U2  ( DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U2  ( DX)  + U2  ( DX)  + U2  ( DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...                                   2            3            4
C...  C(U4 = U2 + U2  (2DX) + U2  (2DX)  + U2  (2DX)  + U2  (2DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U2  (2DX)  + U2  (2DX)  + U2  (2DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...                                   2            3            4
C...  D(U5 = U2 + U2  (3DX) + U2  (3DX)  + U2  (3DX)  + U2  (3DX)
C...                X   1F      2X  2F       3X  3F       4X  4F
C...
C...                      5            6            7
C...           + U2  (3DX)  + U2  (3DX)  + U2  (3DX)  + ...)
C...               5X  5F       6X  6F       7X  7F
C...
C...  -A +   B +  2C +  3D = 1
C...
C...   A +   B +  4C +  9D = 0
C...
C...  -A +   B +  8C + 27D = 0
C...
C...   A +   B + 16C + 81D = 0
C...
C...  SIMULTANEOUS SOLUTION FOR A, B, C AND D FOLLOWED BY THE SOLU-
C...  TION OF THE PRECEDING TAYLOR SERIES, TRUNCATED AFTER THE U
C...                                                            4X
C...  TERMS, FOR U1  GIVES THE FOLLOWING FIVE-POINT APPROXIMATION
C...               X
C...                                                        4
C...  U2  = (1/12DX)(-3U1 - 10U2 + 18U3 -  6U4 +  U5) + O(DX )    (2)
C...    X
C...
C...     (3)  INTERIOR POINT I, I NE 2, N-1
C...
C...                                        2             3
C...  A(UI-2 = UI + UI  (-2DX)  + UI  (-2DX)  + UI  (-2DX)
C...                  X    1F       2X   2F       3X   3F
C...
C...                          4             5             6
C...              + UI  (-2DX)  + UI  (-2DX)  + UI  (-2DX)  + ...)
C...                  4X   4F       5X   5F       6X   6F
C...
C...                                        2             3
C...  B(UI-1 = UI + UI  ( -DX)  + UI  ( -DX)  + UI  ( -DX)
C...                  X    1F       2X   2F       3X   3F
C...
C...                          4             5             6
C...              + UI  ( -DX)  + UI  ( -DX)  + UI  ( -DX)  + ...)
C...                  4X   4F       5X   5F       6X   6F
C...
C...                                        2             3
C...  C(UI+1 = UI + UI  (  DX)  + UI  (  DX)  + UI  (  DX)
C...                  X    1F       2X   2F       3X   3F
C...
C...                          4             5             6
C...              + UI  (  DX)  + UI  (  DX)  + UI  (  DX)  + ...)
C...                  4X   4F       5X   5F       6X   6F
C...
C...                                        2             3
C...  D(UI+2 = UI + UI  ( 2DX)  + UI  ( 2DX)  + UI  ( 2DX)
C...                  X    1F       2X   2F       3X   3F
C...
C...                          4             5             6
C...              + UI  ( 2DX)  + UI  ( 2DX)  + UI  ( 2DX)  + ...)
C...                  4X   4F       5X   5F       6X   6F
C...
C...   -2A -   B +   C +  2D = 1
C...
C...    4A +   B +   C +  4D = 0
C...
C...   -8A -   B +   C +  8D = 0
C...
C...   16A +   B +   C + 16D = 0
C...
C...  SIMULTANEOUS SOLUTION FOR A, B, C AND D FOLLOWED BY THE SOLU-
C...  TION OF THE PRECEDING TAYLOR SERIES, TRUNCATED AFTER THE U
C...                                                            4X
C...  TERMS, FOR U1  GIVES THE FOLLOWING FIVE-POINT APPROXIMATION
C...               X
C...                                                          4
C...  UI  = (1/12DX)(UI-2 - 8UI-1 + 0UI + 8UI+1 - UI+2) + O(DX )  (3)
C...    X
C...
C...     (4)  INTERIOR POINT, I = N-1
C...
C...                                              2               3
C...  A(UN-4 = UN-1 + UN-1  (-3DX)  + UN-1  (-3DX)  + UN-1  (-3DX)
C...                      X    1F         2X   2F         3X   3F
C...
C...                       4               5               6
C...         + UN-1  (-3DX)  + UN-1  (-3DX)  + UN-1  (-3DX)  + ...
C...               4X   4F         5X   5F         6X   6F
C...
C...                                              2               3
C...  B(UN-3 = UN-1 + UN-1  (-2DX)  + UN-1  (-2DX)  + UN-1  (-2DX)
C...                      X    1F         2X   2F         3X   3F
C...
C...                       4               5               6
C...         + UN-1  (-2DX)  + UN-1  (-2DX)  + UN-1  (-2DX)  + ...
C...               4X   4F         5X   5F         6X   6F
C...
C...                                              2               3
C...  C(UN-2 = UN-1 + UN-1  ( -DX)  + UN-1  (- -X)  + UN-1  ( -DX)
C...                      X    1F         2X   2F         3X   3F
C...
C...                       4               5               6
C...         + UN-1  ( -DX)  + UN-1  ( -DX)  + UN-1  ( -DX)  + ...
C...               4X   4F         5X   5F         6X   6F
C...
C...                                              2               3
C...  D(UN   = UN-1 + UN-1  (  DX)  + UN-1  (  DX)  + UN-1  (  DX)
C...                      X    1F         2X   2F         3X   3F
C...
C...                       4               5               6
C...         + UN-1  (  DX)  + UN-1  (  DX)  + UN-1  (  DX)  + ...
C...               4X   4F         5X   5F         6X   6F
C...
C...  -3A -  2B -   C +   D = 1
C...
C...   9A +  4B +   C +   D = 0
C...
C... -27A -  8B -   C +   D = 0
C...
C...  81A + 16B +   C +   D = 0
C...
C...  SIMULTANEOUS SOLUTION FOR A, B, C AND D FOLLOWED BY THE SOLU-
C...  TION OF THE PRECEDING TAYLOR SERIES, TRUNCATED AFTER THE U
C...                                                            4X
C...  TERMS, FOR U1  GIVES THE FOLLOWING FIVE-POINT APPROXIMATION
C...               X
C...                                                                4
C...  UN-1  = (1/12DX)(-UN-4 + 6UN-3 - 18UN-2 + 10UN-1 + 3UN) + O(DX )
C...      X
C...                                                              (4)
C...
C...    (5)  RIGHT END, POINT I = N
C...
C...                                       2             3
C...  A(UN-4 = UN + UN (-4DX)  + UN  (-4DX)  + UN  (-4DX)
C...                  X   1F       2X   2F       3X   3F
C...
C...                         4             5             6
C...             + UN  (-4DX)  + UN  (-4DX)  + UN  (-4DX)  + ...)
C...                 4X   4F       5X   5F       6X   6F
C...
C...                                       2             3
C...  B(UN-3 = UN + UN (-3DX)  + UN  (-3DX)  + UN  (-3DX)
C...                  X   1F       2X   2F       3X   3F
C...
C...                         4             5             6
C...             + UN  (-3DX)  + UN  (-3DX)  + UN  (-3DX)  + ...)
C...                 4X   4F       5X   5F       6X   6F
C...
C...                                       2             3
C...  C(UN-2 = UN + UN (-2DX)  + UN  (-2DX)  + UN  (-2DX)
C...                  X   1F       2X   2F       3X   3F
C...
C...                         4             5             6
C...             + UN  (-2DX)  + UN  (-2DX)  + UN  (-2DX)  + ...)
C...                 4X   4F       5X   5F       6X   6F
C...
C...                                       2             3
C...  D(UN-1 = UN + UN ( -DX)  + UN  ( -DX)  + UN  ( -DX)
C...                  X   1F       2X   2F       3X   3F
C...
C...                         4             5             6
C...             + UN  ( -DX)  + UN  ( -DX)  + UN  ( -DX)  + ...)
C...                 4X   4F       5X   5F       6X   6F
C...
C...   -4A -  3B -  2C -   D = 1
C...
C...   16A +  9B +  4C +   D = 0
C...
C...  -64A - 27B -  8C -   D = 0
C...
C...  256A + 81B + 16C +   D = 0
C...
C...  SIMULTANEOUS SOLUTION FOR A, B, C AND D FOLLOWED BY THE SOLU-
C...  TION OF THE PRECEDING TAYLOR SERIES, TRUNCATED AFTER THE U
C...                                                            4X
C...  TERMS, FOR U1  GIVES THE FOLLOWING FIVE-POINT APPROXIMATION
C...               X
C...                                                                4
C...  UN  = (1/12DX)(3UN-4 - 16UN-3 + 36UN-2 - 48UN-1 + 25UN) + O(DX )
C...    X
C...                                                              (5)
C...
C...  THE WEIGHTING COEFFICIENTS FOR EQUATIONS (1) TO (5) CAN BE
C...  SUMMARIZED AS
C...
C...             -25   48  -36   16   -3
C...
C...              -3  -10   18   -6    1
C...
C...       1/12    1   -8    0    8   -1
C...
C...              -1    6  -18   10    3
C...
C...               3  -16   36  -48   25
C...
C...  WHICH ARE THE COEFFICIENTS REPORTED BY BICKLEY FOR N = 4, M =
C...  1, P = 0, 1, 2, 3, 4 (BICKLEY, W. G., FORMULAE FOR NUMERICAL
C...  DIFFERENTIATION, MATH. GAZ., VOL. 25, 1941.  NOTE - THE BICKLEY
C...  COEFFICIENTS HAVE BEEN DIVIDED BY A COMMON FACTOR OF TWO).
C...
C...  EQUATIONS (1) TO (5) CAN NOW BE PROGRAMMED TO GENERATE THE
C...  DERIVATIVE U (X) OF FUNCTION U(X) (ARGUMENTS U AND UX OF SUB-
C...              X
C...  ROUTINE DSS004 RESPECTIVELY).
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION     DX,R4FDX,U,UX,XL,XU
     !1                 DFLOAT

      DIMENSION U(N),UX(N)

!C...  COMPUTE THE SPATIAL INCREMENT
      DX=(XU-XL)/DFLOAT(N-1)
      R4FDX=1.D+00/(12.D+00*DX)
      NM2=N-2
!C...
!C...  EQUATION (1) (NOTE - THE RHS OF EQUATIONS (1), (2), (3), (4)
!C...  AND (5) HAVE BEEN FORMATTED SO THAT THE NUMERICAL WEIGHTING
!C...  COEFFICIENTS CAN BE MORE EASILY ASSOCIATED WITH THE BICKLEY
!C...  MATRIX ABOVE)
      UX(1)=R4FDX*(-25.D+00*U(1)+48.D+00*U(2)-36.D+00*U(3)+  &   
                    16.D+00*U(4)-3.D+00*U(5))
  
       
!C...
!C...  EQUATION (2)
      UX(  2)=R4FDX*            &
        ( -3.D+00      *U(  1)  &
          -10.D+00     *U(  2)  &
          +18.D+00     *U(  3)  &
          -6.D+00      *U(  4)  &
          +1.D+00      *U(  5))
!C...
!C...  EQUATION (3)
      !dir$ assume_aligned U:64,UX:64
      !dir$ vector aligned
      !dir$ ivdep
      !dir$ vector always
      DO  I=3,NM2
      UX(  I)=R4FDX*          &
      ( +1.D+00     *U(I-2)   &
       -8.D+00      *U(I-1)   &
       +0.D+00      *U(  I)   &
       +8.D+00      *U(I+1)   &
       -1.D+00      *U(I+2))
      END DO
!C...
!C...  EQUATION (4)
      UX(N-1)=R4FDX*          &
      ( -1.D+00      *U(N-4)  &
       +6.D+00      *U(N-3)   &
      -18.D+00      *U(N-2)   &
      +10.D+00      *U(N-1)   &
       +3.D+00      *U(N  ))
C...
C...  EQUATION (5)
      UX(  N)=R4FDX*
      ( +3.D+00      *U(N-4)  &
      -16.D+00      *U(N-3)   &
      +36.D+00      *U(N-2)   &
      -48.D+00      *U(N-1)   &
      +25.D+00      *U(N  ))
     
    END SUBROUTINE


    
    SUBROUTINE DSS012(XL,XU,N,U,UX,V)
          !dir$ optimize:3
          !dir$ attributes forceinline :: DSS012
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS012
#if 0
C...
C...  SUBROUTINE DSS012 IS AN APPLICATION OF FIRST-ORDER DIRECTIONAL
C...  DIFFERENCING IN THE NUMERICAL METHOD OF LINES.  IT IS INTENDED
C...  SPECIFICALLY FOR THE ANALYSIS OF CONVECTIVE SYSTEMS MODELLED BY
C...  FIRST-ORDER HYPERBOLIC PARTIAL DIFFERENTIAL EQUATIONS WITH THE
C...  SIMPLEST FORM
C...
C...                            U  + V*U  = 0                        (1)
C...                             T      X
C...
C...  THE FIRST FIVE PARAMETERS, XL, XU, N, U AND UX, ARE THE SAME
C...  AS FOR SUBROUTINES DSS002 TO DSS010 AS DEFINED IN THOSE ROUTINES.
C...  THE SIXTH PARAMETER, V, MUST BE PROVIDED TO DSS012 SO THAT THE
C...  DIRECTION OF FLOW IN EQUATION (1) CAN BE USED TO SELECT THE
C...  APPROPRIATE FINITE DIFFERENCE APPROXIMATION FOR THE FIRST-ORDER
C...  SPATIAL DERIVATIVE IN EQUATION (1), U .  THE CONVENTION FOR THE
C...  SIGN OF V IS                         X
C...
C...     FLOW LEFT TO RIGHT                 V GT 0
C...     (I.E., IN THE DIRECTION            (I.E., THE SIXTH ARGUMENT IS
C...     OF INCREASING X)                   POSITIVE IN CALLING DSS012)
C...
C...     FLOW RIGHT TO LEFT                 V LT 0
C...     (I.E., IN THE DIRECTION            (I.E., THE SIXTH ARGUMENT IS
C...     OF DECREASING X)                   NEGATIVE IN CALLING DSS012)
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION     DX,      U,     UX,      V,     XL,     XU
     !1                 DFLOAT
!C...
      DIMENSION U(N),UX(N)
!C!...
!C...  COMPUTE THE SPATIAL INCREMENT, THEN SELECT THE FINITE DIFFERENCE
!C...  APPROXIMATION DEPENDING ON THE SIGN OF V IN EQUATION (1).  THE
!C...  ORIGIN OF THE FINITE DIFFERENCE APPROXIMATIONS USED BELOW IS GIVEN
!C...  AT THE END OF SUBROUTINE DSS012.
      DX=(XU-XL)/DFLOAT(N-1)
      IF(V.LT.0.D+00)GO TO 10
!C...
!C...     (1)  FINITE DIFFERENCE APPROXIMATION FOR POSITIVE V
              UX(1)=(U(2)-U(1))/DX
              !dir$ assume_aligned UX:64,U:64
              !dir$ vector aligned
              !dir$ ivdep
              !dir$ vector always
              DO 1 I=2,N
                   UX(I)=(U(I)-U(I-1))/DX
1             CONTINUE
              RETURN
!C...
!C...     (2)  FINITE DIFFERENCE APPROXIMATION FOR NEGATIVE V
10            NM1=N-1
              !dir$ assume_aligned UX:64,U:64
              !dir$ vector aligned
              !dir$ ivdep
              !dir$ vector always
              DO 2 I=1,NM1
                   UX(I)=(U(I+1)-U(I))/DX
2             CONTINUE
              UX(N)=(U(N)-U(N-1))/DX
              RETURN
#if 0
C...
C...  THE BACKWARD DIFFERENCES IN SECTION (1) ABOVE ARE BASED ON THE
C...  TAYLOR SERIES
C...
C...                                  2           3
C...  UI-1 = UI + UI (-DX) + UI  (-DX) + UI  (-DX) + ...
C...                X  1F      2X  2F      3X  3F
C...
C...                                          2
C...  IF THIS SERIES IS TRUNCATED AFTER THE DX  TERM AND THE RESULTING
C...  EQUATION SOLVED FOR U ,  WE OBTAIN IMMEDIATELY
C...                       X
C...
C...  UI  = (UI - UI-1)/DX + O(DX)
C...    X
C...
C...  WHICH IS THE FIRST-ORDER BACKWARD DIFFERENCE USED IN DO LOOP 1.
C...  THE DERIVATIVE U1  IS COMPUTED BY USING THE POINT TO THE RIGHT OF
C...                   X
C...  U1, I.E., U2, SINCE THIS IS THE ONLY POINT AVAILABLE IF FICTITIOUS
C...  POINTS TO THE LEFT OF U1 ARE TO BE AVOIDED.
C...
C...  THE FORWARD DIFFERENCES IN SECTION (2) ABOVE ARE BASED ON THE
C...  TAYLOR SERIES
C...
C...                                  2           3
C...  UI+1 = UI + UI ( DX) + UI  ( DX) + UI  ( DX) + ...
C...                X  1F      2X  2F      3X  3F
C...
C...                                          2
C...  IF THIS SERIES IS TRUNCATED AFTER THE DX  TERM AND THE RESULTING
C...  EQUATION SOLVED FOR U ,  WE OBTAIN IMMEDIATELY
C...                       X
C...
C...  UI  = (UI+1 - UI)/DX + O(DX)
C...    X
C...
C...  WHICH IS THE FIRST-ORDER FORWARD DIFFERENCE USED IN DO LOOP 2.
C...  THE DERIVATIVE UN  IS COMPUTED BY USING THE POINT TO THE LEFT OF
C...                   X
C...  UN (UN-1), SINCE THIS IS THE ONLY POINT AVAILABLE IF FICTITIOUS
C...  POINTS TO THE RIGHT OF UN ARE TO BE AVOIDED.
#endif

   END SUBROUTINE


   
   SUBROUTINE DSS018(XL,XU,N,U,UX,V)
          !dir$ optimize:3
          !dir$ attributes forceinline :: DSS018
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS018
#if 0
C...
C...  SUBROUTINE DSS018 IS AN APPLICATION OF THIRD-ORDER DIRECTIONAL
C...  DIFFERENCING IN THE NUMERICAL METHOD OF LINES.  IT IS INTENDED
C...  SPECIFICALLY FOR THE ANALYSIS OF CONVECTIVE SYSTEMS MODELLED BY
C...  FIRST-ORDER HYPERBOLIC PARTIAL DIFFERENTIAL EQUATIONS AS DIS-
C...  CUSSED IN SUBROUTINE DSS012.  THE COEFFICIENTS OF THE FINITE
C...  DIFFERENCE APPROXIMATIONS USED HEREIN ARE TAKEN FROM BICKLEY, W.
C...  G., FORMULAE FOR NUMERICAL DIFFERENTIATION, THE MATHEMATICAL
C...  GAZETTE, PP. 19-27, 1941, N = 3, M = 1, P = 0, 1, 2, 3.  THE
C...  IMPLEMENTATION IS THE **FOUR-POINT BIASED UPWIND FORMULA** OF
C...  M. B. CARVER AND H. W. HINDS, THE METHOD OF LINES AND THE
C...  ADVECTION EQUATION, SIMULATION, VOL. 31, NO. 2, PP. 59-69,
C...  AUGUST, 1978
C...
#endif
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N),UX(N)
!C...
!C...  COMPUTE THE COMMON FACTOR FOR EACH FINITE DIFFERENCE APPROXIMATION
!C...  CONTAINING THE SPATIAL INCREMENT, THEN SELECT THE FINITE DIFFER-
!C...  ENCE APPROXIMATION DEPENDING ON THE SIGN OF V (SIXTH ARGUMENT).
      DX=(XU-XL)/DFLOAT(N-1)
      R3FDX=1.D+00/(6.D+00*DX)
      IF(V.LT.0.D+00)GO TO 10
!C...
!C...     (1)  FINITE DIFFERENCE APPROXIMATION FOR POSITIVE V
      UX(  1)=R3FDX*           &
      ( -11.D+00     *U(  1)   &
       +18.D+00     *U(  2)   &
        -9.D+00     *U(  3)   &
        +2.D+00     *U(  4))
      UX(  2)=R3FDX*
     (  -2.D+00     *U(  1)   &
        -3.D+00     *U(  2)   &
        +6.D+00     *U(  3)   &
        -1.D+00     *U(  4))
      NM1=N-1
      !dir$ asssume_aligned UX:64,U:64
      !dir$ vector aligned
      !dir$ ivdep
      !dir$ vector vectorlength(8)
      !dir$ vector always
      DO 1 I=3,NM1
      UX(  I)=R3FDX*           &
      (  +1.D+00     *U(I-2)   &
         -6.D+00     *U(I-1)   &
         +3.D+00     *U(I  )   &
         +2.D+00     *U(I+1))
1     CONTINUE
      UX(  N)=R3FDX*           &
      (  -2.D+00      *U(N-3)  &
         +9.D+00      *U(N-2)  &
         -18.D+00     *U(N-1)  &
         +11.D+00     *U(N  ))
      RETURN
!C...
!C...     (2)  FINITE DIFFERENCE APPROXIMATION FOR NEGATIVE V
10    UX(  1)=R3FDX*           &
     ( -11.D+00     *U(  1)    &
       +18.D+00     *U(  2)    &
        -9.D+00     *U(  3)    &
        +2.D+00     *U(  4))
      NM2=N-2
      !dir$ asssume_aligned UX:64,U:64
      !dir$ vector aligned
      !dir$ ivdep
      !dir$ vector vectorlength(8)
      !dir$ vector always
      DO 2 I=2,NM2
      UX(  I)=R3FDX*         &
     (  -2.D+00     *U(I-1)  &
        -3.D+00     *U(I  )  &
        +6.D+00     *U(I+1)  &
        -1.D+00     *U(I+2))
2     CONTINUE
      UX(N-1)=R3FDX*         &
     (  +1.D+00     *U(N-3)  &
        -6.D+00     *U(N-2)  &
        +3.D+00     *U(N-1)  &
        +2.D+00     *U(N  ))
      UX(  N)=R3FDX*         &
     (  -2.D+00     *U(N-3)  &
        +9.D+00     *U(N-2)  &
       -18.D+00     *U(N-1)  &
       +11.D+00     *U(N  ))
    
    END SUBROUTINE


   
     SUBROUTINE DSS020(XL,XU,N,U,UX,V)
          !dir$ optimize:3
          !dir$ attributes forceinline :: DSS020
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS020
#if 0
C...
C...  SUBROUTINE DSS020 IS AN APPLICATION OF FOURTH-ORDER DIRECTIONAL
C...  DIFFERENCING IN THE NUMERICAL METHOD OF LINES.  IT IS INTENDED
C...  SPECIFICALLY FOR THE ANALYSIS OF CONVECTIVE SYSTEMS MODELLED BY
C...  FIRST-ORDER HYPERBOLIC PARTIAL DIFFERENTIAL EQUATIONS AS DIS-
C...  CUSSED IN SUBROUTINE DSS012.  THE COEFFICIENTS OF THE FINITE
C...  DIFFERENCE APPROXIMATIONS USED HEREIN ARE TAKEN FROM BICKLEY, W.
C...  G., FORMULAE FOR NUMERICAL DIFFERENTIATION, THE MATHEMATICAL
C...  GAZETTE, PP. 19-27, 1941, N = 4, M = 1, P = 0, 1, 2, 3, 4.  THE
C...  IMPLEMENTATION IS THE **FIVE-POINT BIASED UPWIND FORMULA** OF
C...  M. B. CARVER AND H. W. HINDS, THE METHOD OF LINES AND THE
C...  ADVECTION EQUATION, SIMULATION, VOL. 31, NO. 2, PP. 59-69,
C...  AUGUST, 1978
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION      DX,  R4FDX,      U,     UX,      V,     XL
     !1                      XU, DFLOAT
!C...
      DIMENSION U(N),UX(N)
!C...
!C...  COMPUTE THE COMMON FACTOR FOR EACH FINITE DIFFERENCE APPROXIMATION
!C...  CONTAINING THE SPATIAL INCREMENT, THEN SELECT THE FINITE DIFFER-
!C...  ENCE APPROXIMATION DEPENDING ON THE SIGN OF V (SIXTH ARGUMENT).
      DX=(XU-XL)/DFLOAT(N-1)
      R4FDX=1.D+00/(12.D+00*DX)
      IF(V.LT.0.D+00)GO TO 10
!C...
!C...     (1)  FINITE DIFFERENCE APPROXIMATION FOR POSITIVE V
      UX(  1)=R4FDX*         &
      ( -25.D+00    *U(  1)  &
       +48.D+00     *U(  2)  &
       -36.D+00     *U(  3)  &
       +16.D+00     *U(  4)  &
        -3.D+00     *U(  5))
      UX(  2)=R4FDX*         &
      (  -3.D+00    *U(  1)  &
       -10.D+00     *U(  2)  &
       +18.D+00     *U(  3)  &
        -6.D+00     *U(  4)  &
        +1.D+00     *U(  5))
      UX(  3)=R4FDX*         &
      (  +1.D+00    *U(  1)  &
        -8.D+00     *U(  2)  &
        +0.D+00     *U(  3)  &
        +8.D+00     *U(  4)  &
        -1.D+00     *U(  5))
      NM1=N-1
      !dir$ asssume_aligned UX:64,U:64
      !dir$ vector aligned
      !dir$ ivdep
      !dir$ vector vectorlength(8)
      !dir$ vector always
      DO 1 I=4,NM1
      UX(  I)=R4FDX*          &
      (  -1.D+00     *U(I-3)  &
         +6.D+00     *U(I-2)  &
        -18.D+00     *U(I-1)  &
        +10.D+00     *U(I  )  &
         +3.D+00     *U(I+1))
1     CONTINUE
      UX(  N)=R4FDX*          &
     (  +3.D+00     *U(N-4)   &
       -16.D+00     *U(N-3)   &
       +36.D+00     *U(N-2)   &
       -48.D+00     *U(N-1)   &
       +25.D+00     *U(N  ))
      RETURN
C...
C...     (2)  FINITE DIFFERENCE APPROXIMATION FOR NEGATIVE V
10    UX(  1)=R4FDX*          &
     ( -25.D+00     *U(  1)   &
       +48.D+00     *U(  2)   &
       -36.D+00     *U(  3)   &
       +16.D+00     *U(  4)   &
        -3.D+00     *U(  5))
      NM3=N-3
      !dir$ asssume_aligned UX:64,U:64
      !dir$ vector aligned
      !dir$ ivdep
      !dir$ vector vectorlength(8)
      !dir$ vector always
      DO 2 I=2,NM3
      UX(  I)=R4FDX*         &
     (  -3.D+00     *U(I-1)  &
       -10.D+00     *U(I  )  &
       +18.D+00     *U(I+1)  &
        -6.D+00     *U(I+2)  &
        +1.D+00     *U(I+3))
2     CONTINUE
      UX(N-2)=R4FDX*         &
     (  +1.D+00     *U(N-4)  &
        -8.D+00     *U(N-3)  &
        +0.D+00     *U(N-2)  &
        +8.D+00     *U(N-1)  &
        -1.D+00     *U(N  ))
      UX(N-1)=R4FDX*         &
     (  -1.D+00     *U(N-4)  &
        +6.D+00     *U(N-3)  &
       -18.D+00     *U(N-2)  &
       +10.D+00     *U(N-1)  &
        +3.D+00     *U(N  ))
      UX(  N)=R4FDX*         &
     (  +3.D+00     *U(N-4)  &
       -16.D+00     *U(N-3)  &
       +36.D+00     *U(N-2)  &
       -48.D+00     *U(N-1)  &
       +25.D+00     *U(N  ))
     
    END SUBROUTINE


    
    SUBROUTINE DSS034(XL,XU,N1,N2,ND,U2D,UX2D,V)
          !dir$ optimize:3
          !dir$ attributes forceinline :: DSS034
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS034
#if 0
C...
C...  SUBROUTINE DSS034 COMPUTES A PARTIAL DERIVATIVE OVER A TWO-
C...  DIMENSIONAL DOMAIN USING EITHER FIVE-POINT CENTERED OR FIVE-
C...  POINT BIASED UPWIND APPROXIMATIONS.  IT IS INTENDED PRIMARILY
C...  FOR THE NUMERICAL METHOD OF LINES (NMOL) NUMERICAL INTEGRATION
C...  OF PARTIAL DIFFERENTIAL EQUATIONS (PDES) IN TWO DIMENSIONS.
C...  DSS034 IS RECOMMENDED PARTICULARLY FOR CONVECTIVE-DIFFUSION
C...  PROBLEMS WHICH REQUIRE COMBINATIONS OF BIASED UPWIND AND CENTERED
C...  APPROXIMATIONS FOR FIRST AND SECOND-ORDER SPATIAL DERIVATIVES IN
C...  PDES.
C...
C...  ARGUMENT LIST
C...
C...     XL        LOWER VALUE OF THE INDEPENDENT VARIABLE FOR WHICH
C...               THE PARTIAL DERIVATIVE IS TO BE COMPUTED (INPUT)
C...
C...     XU        UPPER VALUE OF THE INDEPENDENT VARIABLE FOR WHICH
C...               THE PARTIAL DERIVATIVE IS TO BE COMPUTED (INPUT)
C...
C...     N1        NUMBER OF GRID POINTS FOR THE FIRST INDEPENDENT
C...               VARIABLE (INPUT)
C...
C...     N2        NUMBER OF GRID POINTS FOR THE SECOND INDEPENDENT
C...               VARIABLE (INPUT)
C...
C...     ND        NUMBER OF THE INDEPENDENT VARIABLE FOR WHICH THE
C...               PARTIAL DERIVATIVE IS TO BE COMPUTED (INPUT)
C...
C...     U2D       TWO-DIMENSIONAL ARRAY CONTAINING THE DEPENDENT VARI-
C...               ABLE WHICH IS TO BE DIFFERENTIATED WITH RESPECT TO
C...               INDEPENDENT VARIABLE ND (INPUT)
C...
C...     UX2D      TWO-DIMENSIONAL ARRAY CONTAINING THE PARTIAL DERI-
C...               VATIVE OF THE DEPENDENT VARIABLE WITH RESPECT TO
C...               INDEPENDENT VARIABLE ND (OUTPUT)
C...
C...     V         VARIABLE TO SELECT EITHER THE FIVE-POINT CENTERED
C...               OR FIVE-POINT BIASED UPWIND APPROXIMATION FOR THE
C...               PARTIAL DERIVATIVE.  V EQ 0 CALLS THE FIVE-POINT
C...               CENTERED APPROXIMATION.  V NE 0 CALLS THE FIVE-POINT
C...               BIASED UPWIND APPROXIMATION (INPUT)
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION   UX1D,   UX2D,    U1D,    U2D,      V,     XL
     !1                     XU
!C...
!C...  THE FOLLOWING TWO-DIMENSIONAL ARRAYS CONTAIN THE DEPENDENT
!C...  VARIABLE (U2D) AND ITS PARTIAL DERIVATIVE (UX2D)
      DIMENSION   U2D(N1,N2), UX2D(N1,N2)
#if 0
C...
C...  THE FOLLOWING ONE-DIMENSIONAL ARRAYS CONTAIN THE DEPENDENT
C...  VARIABLE (U1D) AND ITS PARTIAL DERIVATIVE (UX1D).  IN EACH
C...  CASE, ONE OF THE INDEPENDENT VARIABLES IS CONSTANT AND THE
C...  OTHER INDEPENDENT VARIABLE VARIES OVER ITS TOTAL INTERVAL.
C...  THESE ARRAYS ARE USED FOR TEMPORARY STORAGE IN CALLING THE
C...  ONE-DIMENSIONAL ROUTINES DSS004 AND DSS020.
C...
C...  NOTE THAT THE ARRAYS HAVE ABSOLUTE DIMENSIONS AND MAY THERE-
C...  FORE HAVE TO BE INCREASED IN SIZE.  HOWEVER, WITH A SIZE
C...  OF 51, THE TWO-DIMENSIONAL PROBLEM COULD HAVE A GRID OF
C...  51 X 51 POINTS, THEREBY GENERATING AN APPROXIMATING ODE
C...  SYSTEM WITH A MULTIPLE OF 51 X 51 EQUATIONS, DEPENDING ON
C...  THE NUMBER OF SIMULTANEOUS PDES.  THIS IS A VERY LARGE ODE
C...  PROBLEM, AND THEREFORE THE FOLLOWING ABSOLUTE DIMENSIONING
C...  IS CONSIDERED ADEQUATE FOR MOST PROBLEMS.
#endif
      ! Changed to 1024x1024
      !DIMENSION      U1D(51),    UX1D(51)
      DOUBLE PRECISION, DIMENSION(1024) :: U1D
      DOUBLE PRECISION, DIMENSION(1024) :: UX1D
      !dir$ attributes align : 64 :: U1D
      !dir$ attributes align : 64 :: UX1D
!C...
!C...  GO TO STATEMENT 2 IF THE PARTIAL DERIVATIVE IS TO BE COMPUTED
!C...  WITH RESPECT TO THE SECOND INDEPENDENT VARIABLE
      IF(ND.EQ.2)GO TO 2
!C...
!C...  ******************************************************************
!C...
!C...  THE PARTIAL DERIVATIVE IS TO BE COMPUTED WITH RESPECT TO THE
!C...  FIRST INDEPENDENT VARIABLE DEFINED OVER AN INTERVAL CONSISTING
!C...  OF N1 GRID POINTS.  COMPUTE THE PARTIAL DERIVATIVE AT THE N1 X
!C...  N2 GRID POINTS VIA NESTED DO LOOPS 10, 11 AND 12
      DO 10 J=1,N2
!C...
!C...  TRANSFER THE DEPENDENT VARIABLE IN THE TWO-DIMENSIONAL ARRAY U2D
!C...  TO THE ONE-DIMENSIONAL ARRAY U1D SO THAT SUBROUTINES DSS004 AND
!C...  DSS020 CAN BE USED TO CALCULATE THE PARTIAL DERIVATIVE
      !dir$ assume_aligned U1D:64,U2D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 11 I=1,N1
            U1D(I)=U2D(I,J)
11    CONTINUE
!C...
!C...  IF V EQ 0, A FIVE-POINT CENTERED APPROXIMATION IS USED FOR THE
!C...  PARTIAL DERIVATIVE
      IF(V.EQ.0.D+00)CALL DSS004(XL,XU,N1,U1D,UX1D)
!C...
!C...  IF V NE 0, A FIVE-POINT BIASED UPWIND APPROXIMATION IS USED FOR
!C...  THE PARTIAL DERIVATIVE
      IF(V.NE.0.D+00)CALL DSS020(XL,XU,N1,U1D,UX1D,V)
!C...
!C...  RETURN THE PARTIAL DERIVATIVE IN THE ONE-DIMENSIONAL ARRAY UX1D
!C...  TO THE TWO-DIMENSIONAL ARRAY UX2D
      !dir$ assume_aligned UX2D:64,UX1D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 12 I=1,N1
            UX2D(I,J)=UX1D(I)
12    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE AT A PARTICULAR VALUE OF THE SECOND INDE-
!C...  PENDENT VARIABLE HAS BEEN CALCULATED.  REPEAT THE CALCULATION FOR
!C...  THE NEXT VALUE OF THE SECOND INDEPENDENT VARIABLE
10    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE HAS BEEN CALCULATED OVER THE ENTIRE N1 X
!C...  N2 GRID.  THEREFORE RETURN TO THE CALLING PROGRAM WITH THE PARTIAL
!C...  DERIVATIVE IN THE TWO-DIMENSIONAL ARRAY UX2D
      RETURN

!C...  THE PARTIAL DERIVATIVE IS TO BE COMPUTED WITH RESPECT TO THE
!C...  SECOND INDEPENDENT VARIABLE DEFINED OVER AN INTERVAL CONSISTING
!C...  OF N2 GRID POINTS.  COMPUTE THE PARTIAL DERIVATIVE AT THE N1 X
!C...  N2 GRID POINTS VIA NESTED DO LOOPS 20 AND 21
2     DO 20 I=1,N1
!C...
!C...  TRANSFER THE DEPENDENT VARIABLE IN THE TWO-DIMENSIONAL ARRAY U2D
!C...  TO THE ONE-DIMENSIONAL ARRAY U1D SO THAT SUBROUTINES DSS004 AND
!C...  DSS020 CAN BE USED TO CALCULATE THE PARTIAL DERIVATIVE
      !dir$ assume_aligned U1D:64,U2D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 21 J=1,N2
            U1D(J)=U2D(I,J)
21    CONTINUE
!C...
!C...  IF V EQ 0, A FIVE-POINT CENTERED APPROXIMATION IS USED FOR THE
!C!...  PARTIAL DERIVATIVE
      IF(V.EQ.0.D+00)CALL DSS004(XL,XU,N2,U1D,UX1D)
!C...
!C...  IF V NE 0, A FIVE-POINT BIASED UPWIND APPROXIMATION IS USED FOR
!C...  THE PARTIAL DERIVATIVE
      IF(V.NE.0.D+00)CALL DSS020(XL,XU,N2,U1D,UX1D,V)
!C...
!C...  RETURN THE PARTIAL DERIVATIVE IN THE ONE-DIMENSIONAL ARRAY UX1D
!C...  TO THE TWO-DIMENSIONAL ARRAY UX2D
      !dir$ assume_aligned UX2D:64,UX1D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 22 J=1,N2
             UX2D(I,J)=UX1D(J)
22    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE AT A PARTICULAR VALUE OF THE FIRST INDE-
!C...  PENDENT VARIABLE HAS BEEN CALCULATED.  REPEAT THE CALCULATION FOR
!C...  THE NEXT VALUE OF THE FIRST INDEPENDENT VARIABLE
20    CONTINUE

     
   END SUBROUTINE


   
   SUBROUTINE DSS036(XL,XU,N1,N2,N3,ND,U3D,UX3D,V)
            !dir$ optimize:3
          !dir$ attributes forceinline :: DSS036
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS036
#if 0
C...
C...  SUBROUTINE DSS036 COMPUTES A PARTIAL DERIVATIVE OVER A THREE-
C...  DIMENSIONAL DOMAIN USING EITHER FIVE-POINT CENTERED OR FIVE-
C...  POINT BIASED UPWIND APPROXIMATIONS.  IT IS INTENDED PRIMARILY
C...  FOR THE NUMERICAL METHOD OF LINES (NMOL) NUMERICAL INTEGRATION
C...  OF PARTIAL DIFFERENTIAL EQUATIONS (PDES) IN THREE DIMENSIONS.
C...
C...  SUBROUTINE DSS036 IS CALLED IN ESSENTIALLY THE SAME WAY AS
C...  SUBROUTINE DSS034.  THE ONLY DIFFERENCE IS AN ADDITIONAL ARGU-
C...  MENT, N3, TO DEFINE THE NUMBER OF GRID POINTS IN THE THIRD
C...  DIMENSION.  THE COMMENTS IN DSS034 SHOULD THEREFORE BE USEFUL
C...  IN UNDERSTANDING THE OPERATION OF DSS036.  IN PARTICULAR,
C...  DSS036 CALLS SUBROUTINES DSS004 AND DSS020 TO IMPLEMENT THE
C...  FIVE-POINT CENTERED APPROXIMATION AND FIVE-POINT BIASED UPWIND
C...  APPROXIMATION OF THE PARTIAL DERIAVTIVE, RESPECTIVELY.
C...
C...  ARGUMENT LIST
C...
C...     XL        LOWER VALUE OF THE INDEPENDENT VARIABLE FOR WHICH
C...               THE PARTIAL DERIVATIVE IS TO BE COMPUTED (INPUT)
C...
C...     XU        UPPER VALUE OF THE INDEPENDENT VARIABLE FOR WHICH
C...               THE PARTIAL DERIVATIVE IS TO BE COMPUTED (INPUT)
C...
C...     N1        NUMBER OF GRID POINTS FOR THE FIRST INDEPENDENT
C...               VARIABLE (INPUT)
C...
C...     N2        NUMBER OF GRID POINTS FOR THE SECOND INDEPENDENT
C...               VARIABLE (INPUT)
C...
C...     N3        NUMBER OF GRID POINTS FOR THE THIRD INDEPENDENT
C...               VARIABLE (INPUT)
C...
C...     ND        NUMBER OF THE INDEPENDENT VARIABLE FOR WHICH THE
C...               PARTIAL DERIVATIVE IS TO BE COMPUTED (INPUT)
C...
C...     U3D       THREE-DIMENSIONAL ARRAY CONTAINING THE DEPENDENT
C...               VARIABLE WHICH IS TO BE DIFFERENTIATED WITH RESPECT
C...               TO INDEPENDENT VARIABLE ND (INPUT)
C...
C...     UX3D      THREE-DIMENSIONAL ARRAY CONTAINING THE PARTIAL DERI-
C...               VATIVE OF THE DEPENDENT VARIABLE WITH RESPECT TO
C...               INDEPENDENT VARIABLE ND (OUTPUT)
C...
C...     V         VARIABLE TO SELECT EITHER THE FIVE-POINT CENTERED
C...               OR FIVE-POINT BIASED UPWIND APPROXIMATION FOR THE
C...               PARTIAL DERIVATIVE.  V EQ 0 CALLS THE FIVE-POINT
C...               CENTERED APPROXIMATION.  V NE 0 CALLS THE FIVE-POINT
C...               BIASED UPWIND APPROXIMATION (INPUT)
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION   UX1D,UX3D,U1D,U3D,V,XL,XU
                          
!C...
!C...  THE FOLLOWING THREE-DIMENSIONAL ARRAYS CONTAIN THE DEPENDENT
!C...  VARIABLE (U3D) AND ITS PARTIAL DERIVATIVE (UX3D)
      DOUBLE PRECISION, DIMENSION(N1,N2,N3) ::  U3D,UX3D
#if 0
C...
C...  THE FOLLOWING ONE-DIMENSIONAL ARRAYS CONTAIN THE DEPENDENT
C...  VARIABLE (U1D) AND ITS PARTIAL DERIVATIVE (UX1D).  IN EACH
C...  CASE, ONE OF THE INDEPENDENT VARIABLES IS CONSTANT AND THE
C...  OTHER TWO INDEPENDENT VARIABLES VARY OVER THEIR TOTAL INTERVALS.
C...  THESE ARRAYS ARE USED FOR TEMPORARY STORAGE IN CALLING THE
C...  ONE-DIMENSIONAL ROUTINES DSS004 AND DSS020.
C...
C...  NOTE THAT THE ARRAYS HAVE ABSOLUTE DIMENSIONS AND MAY THERE-
C...  FORE HAVE TO BE INCREASED IN SIZE.  HOWEVER, WITH A SIZE
C...  OF 51, THE THREE-DIMENSIONAL PROBLEM COULD HAVE A GRID OF
C...  51 X 51 X 51 POINTS, THEREBY GENERATING AN APPROXIMATING ODE
C...  SYSTEM WITH A MULTIPLE OF 51 X 51 X 51 EQUATIONS, DEPENDING ON
C...  THE NUMBER OF SIMULTANEOUS PDES.  THIS IS A VERY LARGE ODE
C...  PROBLEM, AND THEREFORE THE FOLLOWING ABSOLUTE DIMENSIONING
C...  IS CONSIDERED ADEQUATE FOR MOST PROBLEMS.
#endif
      DOUBLE PRECISION,DIMENSION(1024) :: U1D,UX1D
       dir$ attributes align : 64 :: U1D,UX1D
!C...
!C...  GO TO STATEMENT 2 IF THE PARTIAL DERIVATIVE IS TO BE COMPUTED
!C...  WITH RESPECT TO THE SECOND INDEPENDENT VARIABLE
      IF(ND.EQ.2)GO TO 2
!C...
!C...  GO TO STATEMENT 3 IF THE PARTIAL DERIVATIVE IS TO BE COMPUTED
!C...  WITH RESPECT TO THE THIRD INDEPENDENT VARIABLE
      IF(ND.EQ.3)GO TO 3
!C...
!C...  ******************************************************************
!C...
!C...  THE PARTIAL DERIVATIVE IS TO BE COMPUTED WITH RESPECT TO THE
!C...  FIRST INDEPENDENT VARIABLE DEFINED OVER AN INTERVAL CONSISTING
!C...  OF N1 GRID POINTS.  COMPUTE THE PARTIAL DERIVATIVE AT THE N1 X
!C...  N2 X N3 GRID POINTS VIA NESTED DO LOOPS 10, 11, 12 AND 13
      DO 10 J=1,N2
      DO 11 K=1,N3
!C...
!C...  TRANSFER THE DEPENDENT VARIABLE IN THE THREE-DIMENSIONAL ARRAY U3D
!C...  TO THE ONE-DIMENSIONAL ARRAY U1D SO THAT SUBROUTINES DSS004 AND
!C...  DSS020 CAN BE USED TO CALCULATE THE PARTIAL DERIVATIVE
       !dir$ assume_aligned U1D:64,U3D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 12 I=1,N1
             U1D(I)=U3D(I,J,K)
12    CONTINUE
!C...
!C...  IF V EQ 0, A FIVE-POINT CENTERED APPROXIMATION IS USED FOR THE
!C...  PARTIAL DERIVATIVE
      IF(V.EQ.0.D+00)CALL DSS004(XL,XU,N1,U1D,UX1D)
!C...
!C...  IF V NE 0, A FIVE-POINT BIASED UPWIND APPROXIMATION IS USED FOR
!C...  THE PARTIAL DERIVATIVE
      IF(V.NE.0.D+00)CALL DSS020(XL,XU,N1,U1D,UX1D,V)
!C...
!C...  RETURN THE PARTIAL DERIVATIVE IN THE ONE-DIMENSIONAL ARRAY UX1D
!C...  TO THE THREE-DIMENSIONAL ARRAY UX3D
      !dir$ assume_aligned UX3D:64,UX1D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 13 I=1,N1
            UX3D(I,J,K)=UX1D(I)
13    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE AT PARTICULAR VALUES OF THE SECOND AND
!C...  THIRD INDEPENDENT VARIABLE HAS BEEN CALCULATED.  REPEAT THE
!C...  CALCULATION FOR THE OTHER VALUES OF THE SECOND AND THIRD
!C...  INDEPENDENT VARIABLES
11    CONTINUE
10    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE HAS BEEN CALCULATED OVER THE ENTIRE N1 X
!C...  N2 X N3 GRID.  THEREFORE RETURN TO THE CALLING PROGRAM WITH THE
!C...  PARTIAL DERIVATIVE IN THE THREE-DIMENSIONAL ARRAY UX3D
      RETURN
!C...
!C...  ******************************************************************
!C...
!C...  THE PARTIAL DERIVATIVE IS TO BE COMPUTED WITH RESPECT TO THE
!C...  SECOND INDEPENDENT VARIABLE DEFINED OVER AN INTERVAL CONSISTING
!C...  OF N2 GRID POINTS.  COMPUTE THE PARTIAL DERIVATIVE AT THE N1 X
!C...  N2 X N3 GRID POINTS VIA NESTED DO LOOPS 20, 21, 22 AND 23
2     DO 20 I=1,N1
      DO 21 K=1,N3
!C...
!C...  TRANSFER THE DEPENDENT VARIABLE IN THE THREE-DIMENSIONAL ARRAY U3D
!C...  TO THE ONE-DIMENSIONAL ARRAY U1D SO THAT SUBROUTINES DSS004 AND
!C...  DSS020 CAN BE USED TO CALCULATE THE PARTIAL DERIVATIVE
       !dir$ assume_aligned U1D:64,U3D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 22 J=1,N2
            U1D(J)=U3D(I,J,K)
22    CONTINUE
!C...
!C...  IF V EQ 0, A FIVE-POINT CENTERED APPROXIMATION IS USED FOR THE
!C...  PARTIAL DERIVATIVE
      IF(V.EQ.0.D+00)CALL DSS004(XL,XU,N2,U1D,UX1D)
!C...
!C...  IF V NE 0, A FIVE-POINT BIASED UPWIND APPROXIMATION IS USED FOR
!C...  THE PARTIAL DERIVATIVE
      IF(V.NE.0.D+00)CALL DSS020(XL,XU,N2,U1D,UX1D,V)
!C...
!C...  RETURN THE PARTIAL DERIVATIVE IN THE ONE-DIMENSIONAL ARRAY UX1D
!C...  TO THE THREE-DIMENSIONAL ARRAY UX3D
       !dir$ assume_aligned UX3D:64,UX1D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 23 J=1,N2
            UX3D(I,J,K)=UX1D(J)
23    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE AT PARTICULAR VALUES OF THE FIRST AND
!C...  THIRD INDEPENDENT VARIABLE HAS BEEN CALCULATED.  REPEAT THE
!C...  CALCULATION FOR THE OTHER VALUES OF THE FIRST AND THIRD
!C...  INDEPENDENT VARIABLES
21    CONTINUE
20    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE HAS BEEN CALCULATED OVER THE ENTIRE N1 X
!C...  N2 X N3 GRID.  THEREFORE RETURN TO THE CALLING PROGRAM WITH THE
!C...  PARTIAL DERIVATIVE IN THE THREE-DIMENSIONAL ARRAY UX3D
      RETURN
!C...
!C...  ******************************************************************
!C...
!C...  THE PARTIAL DERIVATIVE IS TO BE COMPUTED WITH RESPECT TO THE
!C...  THIRD INDEPENDENT VARIABLE DEFINED OVER AN INTERVAL CONSISTING
!C...  OF N3 GRID POINTS.  COMPUTE THE PARTIAL DERIVATIVE AT THE N1 X
!C...  N2 X N3 GRID POINTS VIA NESTED DO LOOPS 30, 31, 32 AND 33
3     DO 30 I=1,N1
      DO 31 J=1,N2
!C...
!C...  TRANSFER THE DEPENDENT VARIABLE IN THE THREE-DIMENSIONAL ARRAY U3D
!C...  TO THE ONE-DIMENSIONAL ARRAY U1D SO THAT SUBROUTINES DSS004 AND
!C...  DSS020 CAN BE USED TO CALCULATE THE PARTIAL DERIVATIVE
        !dir$ assume_aligned U1D:64,U3D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 32 K=1,N3
             U1D(K)=U3D(I,J,K)
32    CONTINUE
!C...
!C...  IF V EQ 0, A FIVE-POINT CENTERED APPROXIMATION IS USED FOR THE
!C...  PARTIAL DERIVATIVE
      IF(V.EQ.0.D+00)CALL DSS004(XL,XU,N3,U1D,UX1D)
!C...
!C...  IF V NE 0, A FIVE-POINT BIASED UPWIND APPROXIMATION IS USED FOR
!C...  THE PARTIAL DERIVATIVE
      IF(V.NE.0.D+00)CALL DSS020(XL,XU,N3,U1D,UX1D,V)
!C...
!C...  RETURN THE PARTIAL DERIVATIVE IN THE ONE-DIMENSIONAL ARRAY UX1D
!C...  TO THE THREE-DIMENSIONAL ARRAY UX3D
      !dir$ assume_aligned UX3D:64,UX1D:64
      !dir$ vector aligned
      !dir$ unroll(16)
      !dir$ vector always
      DO 33 K=1,N3
             UX3D(I,J,K)=UX1D(K)
33    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE AT PARTICULAR VALUES OF THE FIRST AND
!C...  SECOND INDEPENDENT VARIABLE HAS BEEN CALCULATED.  REPEAT THE
!C...  CALCULATION FOR THE OTHER VALUES OF THE FIRST AND SECOND
!C...  INDEPENDENT VARIABLES
31    CONTINUE
30    CONTINUE
!C...
!C...  THE PARTIAL DERIVATIVE HAS BEEN CALCULATED OVER THE ENTIRE N1 X
!C...  N2 X N3 GRID.  THEREFORE RETURN TO THE CALLING PROGRAM WITH THE
!C...  PARTIAL DERIVATIVE IN THE THREE-DIMENSIONAL ARRAY UX3D
    
    END SUBROUTINE


   
    SUBROUTINE DSS042(XL,XU,N,U,UX,UXX,NL,NU)
            !dir$ optimize:3
          !dir$ attributes forceinline :: DSS042
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS042
#if 0
C...
C...  SUBROUTINE DSS042 COMPUTES A SECOND-ORDER APPROXIMATION OF A
C...  SECOND-ORDER DERIVATIVE, WITH OR WITHOUT THE NORMAL DERIVATIVE
C...  AT THE BOUNDARY.
C...
C...     DOUBLE PRECISION
C...
C...  ARGUMENT LIST
C...
C...     XL      LEFT VALUE OF THE SPATIAL INDEPENDENT VARIABLE (INPUT)
C...
C...     XU      RIGHT VALUE OF THE SPATIAL INDEPENDENT VARIABLE (INPUT)
C...
C...     N       NUMBER OF SPATIAL GRID POINTS, INCLUDING THE END
C...             POINTS (INPUT)
C...
C...     U       ONE-DIMENSIONAL ARRAY OF THE DEPENDENT VARIABLE TO BE
C...             DIFFERENTIATED (INPUT)
C...
C...     UX      ONE-DIMENSIONAL ARRAY OF THE FIRST DERIVATIVE OF U.
C...             THE END VALUES OF UX, UX(1) AND UX(N), ARE USED IN
C...             NEUMANN BOUNDARY CONDITIONS AT X = XL AND X = XU,
C...             DEPENDING ON THE ARGUMENTS NL AND NU (SEE THE DE-
C...             SCRIPTION OF NL AND NU BELOW)
C...
C...     UXX     ONE-DIMENSIONAL ARRAY OF THE SECOND DERIVATIVE OF U
C...             (OUTPUT)
C...
C...     NL      INTEGER INDEX FOR THE TYPE OF BOUNDARY CONDITION AT
C...             X = XL (INPUT).  THE ALLOWABLE VALUES ARE
C...
C...                1 - DIRICHLET BOUNDARY CONDITION AT X = XL
C...                    (UX(1) IS NOT USED)
C...
C...                2 - NEUMANN BOUNDARY CONDITION AT X = XL
C...                    (UX(1) IS USED)
C...
C...     NU      INTEGER INDEX FOR THE TYPE OF BOUNDARY CONDITION AT
C...             X = XU (INPUT).  THE ALLOWABLE VALUES ARE
C...
C...                1 - DIRICHLET BOUNDARY CONDITION AT X = XU
C...                    (UX(N) IS NOT USED)
C...
C...                2 - NEUMANN BOUNDARY CONDITION AT X = XU
C...                    (UX(N) IS USED)
C...
C...  TYPE REAL VARIABLES AS DOUBLE PRECISION
#endif
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION,DIMENSION(N) :: U, UX, UXX
      !dir$ attributes align : 64 :: U,UX,UXX
#if 0
C...
C...  THE FOLLOWING DERIVATION IS FOR A SET OF SECOND-ORDER, FOUR-POINT
C...  APPROXIMATIONS FOR A SECOND DERIVATIVE THAT CAN BE USED AT THE
C...  BOUNDARIES OF A SPATIAL DOMAIN.  THESE APPROXIMATIONS HAVE THE
C...  FEATURES
C...
C...     (1)  ONLY INTERIOR AND BOUNDARY POINTS ARE USED (I.E., NO
C...          FICTITIOUS POINTS ARE USED)
C...
C...     (2)  THE NORMAL DERIVATIVE AT THE BOUNDARY IS INCLUDED AS PART
C...          OF THE APPROXIMATION FOR THE SECOND DERIVATIVE
C...
C...     (3)  APPROXIMATIONS FOR THE BOUNDARY CONDITIONS ARE NOT USED.
C...
C...  THE DERIVATION IS BY PROFESSOR GILBERT A. STENGLE, DEPARTMENT OF
C...  MATHEMATICS, LEHIGH UNIVERSITY, BETHLEHEM, PA 18015, AND WAS DONE
C...  ON DECEMBER 7, 1985.
C...
C...  FOR AN APPROXIMATION AT THE LEFT BOUNDARY, INVOLVING THE POINTS
C...  I, I+1, I+2 AND I+3, CONSIDER THE FOLLOWING TAYLOR SERIES EXPAN-
C...  SIONS
C...
C...                  UX(I)( DX)   UXX(I)( DX)**2   UXXX(I)( DX)**3
C...  U(I+1) = U(I) + ---------- + -------------- + --------------- +...
C...                       1             2                 6
C...
C...
C...                  UX(I)(2DX)   UXX(I)(2DX)**2   UXXX(I)(2DX)**3
C...  U(I+2) = U(I) + ---------- + -------------- + --------------- +...
C...                       1             2                 6
C...
C...  IF WE NOW FORM THE FOLLOWING LINEAR COMBINATION, INVOLVING CON-
C...  STANTS A, B, C AND D TO BE DETERMINED, AND USE THE PRECEDING TWO
C...  TAYLOR SERIES,
C...
C...     A*U(I) + B*UX(I) + C*U(I+1) + D*U(I+2)
C...
C...  WE HAVE
C...
C...     A*U(I) + B*UX(I) + C*U(I+1) + D*U(I+2) =
C...
C...     (A + B + C + D)*U(I) +
C...
C...     (B + DX*C + 2*DX*D)*UX(I) +
C...
C...     (C*(DX**2)/2 + D*((2*DX)**2)/2)*UXX(I) +
C...
C...     (C*(DX**3)/6 + D*((2*DX)**3)/6)*UXXX(I) + O(DX**4)
C...
C...  THE THIRD DERIVATIVE, UXXX(I), CAN BE DROPPED BY TAKING
C...
C...     C = -8*D
C...
C...  THE SECOND DERIVATIVE, UXX(I), CAN BE RETAINED BY TAKING
C...
C...     (DX**2)(C/2 + 2*D) = 1
C...
C...  WHICH, WHEN COMBINED WITH THE PRECEDING RESULT GIVES
C...
C...     D = -1/(2*(DX**2))
C...
C...     C = 4/(DX**2)
C...
C...  THE FIRST DERIVATIVE, UX(I), CAN BE DROPPED BY TAKING
C...
C...     B + DX*C + 2*DX*D = 0
C...
C...  OR
C...
C...     B = -DX*C - 2*DX*D = -4/DX - 2*DX*(-1/(2*(DX**2))) = -3/DX
C...
C...  FINALLY, U(I), CAN BE DROPPED BY TAKING
C...
C...     A = - C - D = 8*D - D = -7*D = -7/(2*(DX**2))
C...
C...  IF WE NOW SOLVE FOR THE DERIVATIVE OF INTEREST, UXX(I),
C...
C...     UXX(I) = -7/(2(DX**2))*U(I) - 3/DX*UX(I)
C...
C...              + 8/(DX**2)*U(I+1) - 1/(2*(DX**2))U(I+2) + O(DX**2)
C...
C...       = (1/(2*(DX**2)))*(-U(I+2) + 8*U(I+1) - 7*U(I) - 6*DX*UX(I))
C...
C...         + O(DX**2)
C...
C...  WHICH IS THE FOUR-POINT, SECOND-ORDER APPROXIMATION FOR THE SECOND
C...  DERIVATIVE, UXX(I), INCLUDING THE FIRST DERIVATIVE, UX(I).
C...
C...  FOUR CHECKS OF THIS APPROXIMATION CAN EASILY BE MADE FOR U(I) =
C...  1, U(I) = X, U(I) = X**2 AND U(I) = X**3
C...
C...     UXX(I) = (1/(2*(DX**2)))*(-1 + 8*1 - 7*1 - 6*DX*0) = 0
C...
C...     UXX(I) = (1/(2*(DX**2)))*(-(X + 2*DX) + 8*(X + DX)
C...
C...              -7*X - 6*DX*1) = 0
C...
C...     UXX(I) = (1/(2*(DX**2)))*(-(X + 2*DX)**2 + 8*(X + DX)**2
C...
C...            - 7*(X**2) - 6*DX*(2*X))
C...
C...             = (-  X**2 -  4*X*DX - 4*DX**2
C...
C...               + 8*X**2 + 16*X*DX + 8*DX**2
C...
C...               - 7*X**2 - 12*X*DX)/(2*(DX**2)) = 2
C...
C...     UXX(I) = (1/(2*(DX**2)))*(-(X + 2*DX)**3 + 8*(X + DX)**3
C...
C...            - 7*(X**3) - 6*DX*(3*X**2))
C...
C...            = (1/(2*(DX**2)))*(- X**3 - 6*DX*X**2 - 12*X*DX**2
C...
C...            - 8*DX**3 + 8*X**3 + 24*DX*X**2 + 24*X*DX**2 + 8*DX**3
C...
C...            - 7*X**3 - 18*DX*X**2)
C...
C...            = (1/(2*(DX**2)))*(12*X*DX**2) = 6*X
C...
C...  THE PRECEDING APPROXIMATION FOR UXX(I) CAN BE APPLIED AT THE
C...  LEFT BOUNDARY VALUE OF X BY TAKING I = 1.  AN APPROXIMATION AT
C...  THE RIGHT BOUNDARY IS OBTAINED BY TAKING DX = -DX AND REVERSING
C...  THE SUBSCRIPTS IN THE PRECEDING APPROXIMATION, WITH I = N
C...
C...     UXX(I)
C...
C...       = (1/(2*(DX**2)))*(-U(I-2) + 8*U(I-1) - 7*U(I) + 6*DX*UX(I))
C...
C...         + O(DX**2)
C...
C...  TO OBTAIN APPROXIMATIONS OF THE SECOND DERVIAVTIVE WHICH DO NOT
C...  INVOLVE THE FIRST DERIVATIVE, WE TAKE AS THE LINEAR COMBINATION
C...
C...     A*U(I) + B*U(I+1) + C*U(I+2) + D*U(I+3) =
C...
C...  WE HAVE
C...
C...     A*U(I) + B*U(I+1) + C*U(I+2) + D*U(I+3) =
C...
C...     (A + B + C + D)*U(I)+
C...
C...     (DX*B + 2*DX*C + 4*DX*D)*UX(I)+
C...
C...     (B*(DX**2)/2 + C*((2*DX)**2)/2 + D*((3*DX)**2)/2)*UXX(I) +
C...
C...     (B*(DX**3)/6 + C*((2*DX)**3)/6 + D*((3*DX)**3)/6)*UXX(I) +
C...
C...     O(DX**4)
C...
C...  THE THIRD DERIVATIVE, UXXX(I), CAN BE DROPPED BY TAKING
C...
C...     B + 8*C + 27*D = 0
C...
C...  THE SECOND DERIVATIVE, UXX(I), CAN BE RETAINED BY TAKING
C...
C...     (DX**2)*(B/2 + 2*C + (9/2)*D) = 1
C...
C...  THE FIRST DERIVATIVE CAN BE DROPPED BY TAKING
C...
C...     B + 2*C + 3*D = 0
C...
C...  SOLUTION OF THE PRECEDING EQUATIONS FOR C AND D BY ELIMINATION OF
C...  B GIVES
C...
C...     6*C + 24*D = 0
C...
C...     4*C + 18*D = -2/(DX**2)
C...
C...  THEN, ELIMINATING C, GIVES
C...
C...     (18 - 16)*D = -2/(DX**2)
C...
C...  OR
C...
C...     D = -1/(DX**2)
C...
C...     C = (24/6)/(DX**2) = 4/(DX**2)
C...
C...     B = -8/(DX**2) + 3/(DX**2) = -5/(DX**2)
C...
C...  U(I) CAN BE DROPPED BY TAKING
C...
C...     A + B + C + D = 0
C...
C...  OR
C...
C...     A = (5 - 4 + 1)/(DX**2) = 2/(DX**2)
C...
C...  IF WE NOW SOLVE FOR THE DERIVATIVE OF INTEREST, UXX(I),
C...
C...     UXX(I) = (1/DX**2)*(2*U(I) - 5*U(I+1) + 4*U(I+2) - 1*U(I+3))
C...
C...            + O(DX**2)
C...
C...  WHICH IS THE FOUR-POINT, SECOND-ORDER APPROXIMATION FOR THE SECOND
C...  DERIVATIVE, UXX(I), WITHOUT THE FIRST DERIVATIVE, UX(I).
C...
C...  FOUR CHECKS OF THIS APPROXIMATION CAN EASILY BE MADE FOR U(I) =
C...  1, U(I) = X, U(I) = X**2 AND U(I) = X**3
C...
C...     UXX(I) = (1/DX**2)*(2 - 5 + 4 - 1) = 0
C...
C...     UXX(I) = (1/DX**2)*(2*X - 5*(X + DX) + 4*(X + 2*DX)
C...
C...              - 1*(X + 3*DX)) = 0
C...
C...     UXX(I) = (1/DX**2)*(2*X**2 - 5*(X + DX)**2 + 4*(X + 2*DX)**2
C...
C...              - 1*(X + 3*DX)**2) = 2
C...
C...     UXX(I) = (1/DX**2)*(2*X**3 - 5*(X + DX)**3 + 4*(X + 2*DX)**3
C...
C...            - 1*(X + 3*DX)**3)
C...
C...             = (1/DX**2)*(2*X**3 - 5*X**3 - 15*X*DX**2
C...
C...             - 15*DX*X**2 - 5*DX**3 + 4*X**3 + 24*DX*X**2
C...
C...             + 48*X*DX**2 + 32*DX**3 - X**3 - 9*DX*X**2
C...
C...             - 27*X*DX**2 - 27DX**3)
C...
C...             = (1/DX**2)*(6*X*DX**2) = 6*X
C...
C...  THE PRECEDING APPROXIMATION FOR UXX(I) CAN BE APPLIED AT THE
C...  LEFT BOUNDARY VALUE OF X BY TAKING I = 1.  AN APPROXIMATION AT
C...  THE RIGHT BOUNDARY IS OBTAINED BY TAKING DX = -DX AND REVERSING
C...  THE SUBSCRIPTS IN THE PRECEDING APPROXIMATION, WITH I = N
C...
C...     UXX(I) = (1/DX**2)*(2*U(I) - 5*U(I-1) + 4*U(I-2) - 1*U(I-3))
C...
C...            + O(DX**2)
C...
C...  GRID SPACING
#endif
      DX=(XU-XL)/DFLOAT(N-1)
!C...
!C...  CALCULATE UXX AT THE LEFT BOUNDARY, WITHOUT UX
      IF(NL.EQ.1)THEN
      UXX(1)=((     2.D0)*U(  1) &
             +(    -5.D0)*U(  2) &
             +(     4.D0)*U(  3) &
             +(    -1.D0)*U(  4))/(DX**2)
!C...
!C...  CALCULATE UXX AT THE LEFT BOUNDARY, INCLUDING UX
      ELSE IF(NL.EQ.2)THEN
      UXX(1)=((    -7.D0)*U(  1)  &
             +(     8.D0)*U(  2)  &
             +(    -1.D0)*U(  3))/(2.D0*DX**2) &
             +(    -6.D0)*UX( 1) /(2.D0*DX)
      END IF
!C...
!C...  CALCULATE UXX AT THE RIGHT BOUNDARY, WITHOUT UX
      IF(NU.EQ.1)THEN
      UXX(N)=((     2.D0)*U(N  )   &
              +(    -5.D0)*U(N-1)  &
              +(     4.D0)*U(N-2)  &
              +(    -1.D0)*U(N-3))/(DX**2)
!C...
!C...  CALCULATE UXX AT THE RIGHT BOUNDARY, INCLUDING UX
      ELSE IF(NU.EQ.2)THEN
      UXX(N)=((    -7.D0)*U(N  )  &
            +(     8.D0)*U(N-1)   &
            +(    -1.D0)*U(N-2))/(2.D0*DX**2) &
            +(     6.D0)*UX(N ) /(2.D0*DX)
      END IF
!C...
!C...  CALCULATE UXX AT THE INTERIOR GRID POINTS
      !dir$ assume_aligned UXX:64,U:64
      !dir$ ivdep
      !dir$ vector aligned
      !dir$ vector always
      DO 1 I=2,N-1
            UXX(I)=(U(I+1)-2.D0*U(I)+U(I-1))/DX**2
1     CONTINUE
    
   END SUBROUTINE

   
   SUBROUTINE DSS044(XL,XU,N,U,UX,UXX,NL,NU)
             !dir$ optimize:3
          !dir$ attributes forceinline :: DSS044
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: DSS044
#if 0
C...
C...  SUBROUTINE DSS044 COMPUTES A FOURTH-ORDER APPROXIMATION OF A
C...  SECOND-ORDER DERIVATIVE, WITH OR WITHOUT THE NORMAL DERIVATIVE
C...  AT THE BOUNDARY.
C...
C...     SINGLE PRECISION
C...
C...  ARGUMENT LIST
C...
C...     XL      LEFT VALUE OF THE SPATIAL INDEPENDENT VARIABLE (INPUT)
C...
C...     XU      RIGHT VALUE OF THE SPATIAL INDEPENDENT VARIABLE (INPUT)
C...
C...     N       NUMBER OF SPATIAL GRID POINTS, INCLUDING THE END
C...             POINTS (INPUT)
C...
C...     U       ONE-DIMENSIONAL ARRAY OF THE DEPENDENT VARIABLE TO BE
C...             DIFFERENTIATED (INPUT)
C...
C...     UX      ONE-DIMENSIONAL ARRAY OF THE FIRST DERIVATIVE OF U.
C...             THE END VALUES OF UX, UX(1) AND UX(N), ARE USED IN
C...             NEUMANN BOUNDARY CONDITIONS AT X = XL AND X = XU,
C...             DEPENDING ON THE ARGUMENTS NL AND NU (SEE THE DE-
C...             SCRIPTION OF NL AND NU BELOW)
C...
C...     UXX     ONE-DIMENSIONAL ARRAY OF THE SECOND DERIVATIVE OF U
C...             (OUTPUT)
C...
C...     NL      INTEGER INDEX FOR THE TYPE OF BOUNDARY CONDITION AT
C...             X = XL (INPUT).  THE ALLOWABLE VALUES ARE
C...
C...                1 - DIRICHLET BOUNDARY CONDITION AT X = XL
C...                    (UX(1) IS NOT USED)
C...
C...                2 - NEUMANN BOUNDARY CONDITION AT X = XL
C...                    (UX(1) IS USED)
C...
C...     NU      INTEGER INDEX FOR THE TYPE OF BOUNDARY CONDITION AT
C...             X = XU (INPUT).  THE ALLOWABLE VALUES ARE
C...
C...                1 - DIRICHLET BOUNDARY CONDITION AT X = XU
C...                    (UX(N) IS NOT USED)
C...
C...                2 - NEUMANN BOUNDARY CONDITION AT X = XU
C...                    (UX(N) IS USED)
C...
C...  TYPE REAL VARIABLES AS DOUBLE PRECISION
#endif
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION,DIMENSION(N) :: U,UX,UXX
      !dir$ attributes align : 64 :: U,UX,UXX
#if 0
C...
C...  THE FOLLOWING DERIVATION WAS COMPLETED BY W. E. SCHIESSER, DEPTS
C...  OF CHE AND CSEE, LEHIGH UNIVERSITY, BETHLEHEM, PA 18015, USA, ON
C...  DECEMBER 15, 1986.  ADDITIONAL DETAILS ARE GIVEN IN SUBROUTINE
C...  DSS042.
C...
C...  ******************************************************************
C...
C...  (1)  UXX AT THE INTERIOR POINTS 3, 4,..., N-2
C...
C...  TO DEVELOP A SET OF FOURTH-ORDER CORRECT DIFFERENTIATION FORMULAS
C...  FOR THE SECOND DERIVATIVE UXX, WE CONSIDER FIRST THE INTERIOR
C...  GRID POINTS AT WHICH A SYMMETRIC FORMULA CAN BE USED.
C...
C...  IF WE CONSIDER A FORMULA OF THE FORM
C...
C...     A*U(I-2) + B*U(I-1) + E*U(I) + C*U(I+1) + D*U(I+2)
C...
C...  TAYLOR SERIES EXPANSIONS OF U(I-2), U(I-1), U(I+1) AND U(I+2)
C...  CAN BE SUBSTITUTED INTO THIS FORMULA.  WE THEN CONSIDER THE
C...  LINEAR ALBEGRAIC EQUATIONS RELATING A, B, C AND D WHICH WILL
C...  RETAIN CERTAIN TERMS, I.E., UXX, AND DROP OTHERS, E.G., UXXX,
C...  UXXXX AND UXXXXX.
C...
C...  THUS, FOR GRID POINTS 3, 4,..., N-2
C...
C...     TO RETAIN UXX
C...
C...        4*A + B + C +  4*D = 2                              (1)
C...
C...     TO DROP UXXX
C...
C...       -8*A - B + C +  8*D = 0                              (2)
C...
C...     TO DROP UXXXX
C...
C...       16*A + B + C + 16*D = 0                              (3)
C...
C...     TO DROP UXXXXX
C...
C...      -32*A - B + C + 32*D = 0                              (4)
C...
C...  EQUATIONS (1) TO (4) CAN BE SOLVED FOR A, B, C AND D.  IF EQUA-
C...  TION (1) IS ADDED TO EQUATION (2)
C...
C...        -4*A + 2*C + 12*D = 2                               (5)
C...
C...  IF EQUATION (1) IS SUBTRACTED FROM EQUATION (3)
C...
C...        12*A + 12*D = -2                                    (6)
C...
C...  IF EQUATION (1) IS ADDED TO EQUATION (4)
C...
C...        -28*A + 2*C + 36*D = 2                              (7)
C...
C...  EQUATIONS (5) TO (7) CAN BE SOLVED FOR A, C AND D.  IF EQUATION
C...  (5) IS SUBTRACTED FROM EQUATION (7), AND THE RESULT COMBINED
C...  WITH EQUATION (6)
C...
C...         12*A + 12*D = -2                                   (6)
C...
C...        -24*A + 24*D = 0                                    (8)
C...
C...  EQUATIONS (6) AND (8) CAN BE SOLVED FOR A AND D.  FROM (8), A
C...  = D.  FROM EQUATION (6), A = -1/12 AND D = -1/12.  THEN, FROM
C...  EQUATION (5), C = 4/3, AND FROM EQUATION (1), B = 4/3.
C...
C...  THE FINAL DIFFERENTIATION FORMULA IS THEN OBTAINED AS
C...
C...     (-1/12)*U(I-2) +   (4/3)*U(I-1) +
C...
C...       (4/3)*U(I+1) + (-1/12)*U(I+2)
C...
C...     (-1/12 + 4/3 - 1/12 + 4/3)*U(I) + UXX(I)*(DX**2) + O(DX**6)
C...
C...  OR
C...
C...     UXX(I) = (1/(12*DX**2))*(-1*U(I-2) + 16*U(I-1)
C...
C...                  - 30*U(I) + 16*U(I+1) -  1*U(I+2)         (9)
C...
C...                  + O(DX**4)
C...
C...  NOTE THAT THE UX TERM DROPS OUT, I.E., THE BASIC EQUATION IS
C...
C...        -2*A - B + C + 2*D =
C...
C...        -2*(-1/12) - (4/3) + (4/3) + 2*(-1/12) = 0
C...
C...  EQUATION (9) WAS OBTAINED BY DROPPING ALL TERMS IN THE UNDERLYING
C...  TAYLOR SERIES UP TO AND INCLUDING THE FIFTH DERIVATIVE, UXXXXX.
C...  THUS, EQUATION (9) IS EXACT FOR POLYNOMIALS UP TO AND INCLUDING
C...  FIFTH ORDER.  THIS CAN BE CHECKED BY SUBSTITUTING THE FUNCTIONS
C...  1, X, X**2, X**3, X**4 AND X**5 IN EQUATION (9) AND COMPUTING THE
C...  CORRESPONDING DERIVATIVES FOR COMPARISON WITH THE KNOWN SECOND
C...  DERIVATIVES.  THIS IS DONE FOR 1 MERELY BY SUMMING THE WEIGHTING
C...  COEFFICIENTS IN EQUATION (9), WHICH SHOULD SUM TO ZERO, I.E.,
C...  -1 + 16 - 30 + 16 -1 = 0.
C...
C...  FOR THE REMAINING FUNCTIONS, THE ALGEBRA IS RATHER INVOLVED, BUT
C...  THESE FUNCTIONS CAN BE CHECKED NUMERICALLY, I.E., NUMERICAL VALUES
C...  OF X**2, X**3, X**4 AND X**5 CAN BE SUBSTITUTED IN EQUATION (9)
C...  AND THE COMPUTED DERIVATIVES CAN BE COMPARED WITH THE KNOW NUMERI-
C...  CAL SECOND DERIVATIVES.  THIS IS NOT A PROOF OF CORRECTNESS OF
C...  EQUATION (9), BUT WOULD LIKELY DETECT ANY ERRORS IN EQUATION (9).
C...
C...  CHECK THE VALUES OF A, B, C AND D IN EQUATIONS (1) TO (4)
C...  CALL CHECK(1)
C...
C...  ******************************************************************
C...
C...  (2)  UXX AT THE INTERIOR POINTS I = 2 AND N-1
C...
C...  FOR GRID POINT 2, WE CONSIDER A FORMULA OF THE FORM
C...
C...     A*U(I-1) + F*U(I) + B*U(I+1) + C*U(I+2) + D*U(I+3) + E*U(I+4)
C...
C...  TAYLOR SERIES EXPANSIONS OF U(I-1), U(I+1), U(I+2), U(I+3) AND
C...  U(I+4) WHEN SUBSTITUTED INTO THIS FORMULA GIVE LINEAR ALGEBRAIC
C...  EQUATIONS RELATING A, B, C, D AND E.
C...
C...     TO DROP UX
C...
C...        -A + B + 2*C + 3*D + 4*E = 0                       (10)
C...
C...     TO RETAIN UXX
C...
C...        A + B + 4*C + 9*D + 16*E = 2                       (11)
C...
C...     TO DROP UXXX
C...
C...       -A + B + 8*C + 27*D + 64*E = 0                      (12)
C...
C...     TO DROP UXXXX
C...
C...        A + B + 16*C + 81*D + 256*E = 0                    (13)
C...
C...     TO DROP UXXXXX
C...
C...       -A + B + 32*C + 243*D + 1024*E = 0                  (14)
C...
C...  EQUATIONS (11), (12), (13) AND (14) CAN BE SOLVED FOR A, B, C,
C...  D AND E.  IF EQUATION (10) IS ADDED TO EQUATION (11)
C...
C...     2*B + 6*C + 12*D +20*E = 2                            (15)
C...
C...  IF EQUATION (10) IS SUBTRACTED FROM EQUATION (12)
C...
C...     6*C + 24*D + 60*E = 0                                 (16)
C...
C...  IF EQUATION (10) IS ADDED TO EQUATION (13)
C...
C...     2*B + 18*C + 84*D + 260*E = 0                         (17)
C...
C...  IF EQUATION (10) IS SUBTRACTED FROM EQUATION (14)
C...
C...     30*C + 240*D + 1020*E = 0                             (18)
C...
C...  EQUATIONS (15), (16), (17) AND (18) CAN BE SOLVED FOR B, C, D
C...  AND E.
C...
C...     6*C + 24*D + 60*E = 0                                 (16)
C...
C...  IF EQUATION (15) IS SUBTRACTED FROM EQUATION (17)
C...
C...     12*C + 72*D + 240*E = -2                              (19)
C...
C...     30*C + 240*D + 1020*E = 0                             (18)
C...
C...  EQUATIONS (16), (18) AND (19) CAN BE SOLVED FOR C, D AND E.  IF
C...  TWO TIMES EQUATION (16) IS SUBTRACTED FROM EQUATION (19),
C...
C...     24*D + 120*E = -2                                     (20)
C...
C...  IF FIVE TIMES EQUATION (16) IS SUBTRACTED FROM EQUATION (18),
C...
C...     120*D + 720*E = 0                                     (21)
C...
C...  EQUATIONS (20) AND (21) CAN BE SOLVED FOR D AND E.  FROM (21),
C...  E = (-1/6)*D.  SUBSTITUTION IN EQUATION (20) GIVES D = -1/2.
C...  THUS, E = 1/12.  FROM EQUATION (16), C = 7/6.  FROM EQUATION
C...  (15), B = -1/3.  FROM EQUATION (10), A = 5/6.
C...
C...  THE FINAL DIFFERENTIATION FORMULA IS THEN OBTAINED AS
C...
C...  (5/6)*U(I-1) + (-1/3)*U(I+1) + (7/6)*U(I+2) + (-1/2)*U(I+3)
C...
C...  + (1/12)*U(I+4) = (5/6 - 1/3 + 7/6 - 1/2 + 1/12)*U(I)
C...
C...  + UXX*(DX**2) + O(DX**6)
C...
C...  OR
C...
C...  UXX(I) = (1/12*DX**2)*(10*U(I-1) - 15*U(I) - 4*U(I+1)
C...                                                           (22)
C...         + 14*U(I+2) - 6*U(I+3) + 1*U(I+4)) + O(DX**4)
C...
C...  EQUATION (22) WILL BE APPLIED AT I = 2 AND N-1.  THUS
C...
C...  UXX(2) = (1/12*DX**2)*(10*U(1) - 15*U(2) - 4*U(3)
C...                                                           (23)
C...         + 14*U(4) - 6*U(5) + 1*U(6)) + O(DX**4)
C...
C...  UXX(N-1) = (1/12*DX**2)*(10*U(N) - 15*U(N-1) - 4*U(N-2)
C...                                                           (24)
C...           + 14*U(N-3) - 6*U(N-4) + 1*U(N-5)) + O(DX**4)
C...
C...  CHECK THE VALUES OF A, B, C, D AND E IN EQUATIONS (10) TO (14)
C...  CALL CHECK(2)
C...
C...  ******************************************************************
C...
C...  (3)  UXX AT THE BOUNDARY POINTS 1 AND N
C...
C...  FINALLY, FOR GRID POINT 1, AN APPROXIMATION WITH A NEUMANN BOUND-
C...  ARY CONDITON OF THE FORM
C...
C...     A*U(I+1) + B*U(I+2) + C*U(I+3) + D*U(I+4) + E*UX(I) + F*U(I)
C...
C...  WILL BE USED.  THE CORRESPONDING ALGEBRAIC EQUATIONS ARE
C...
C...     TO DROP UX
C...
C...        A + 2*B + 3*C + 4*D + E = 0                        (25)
C...
C...     TO RETAIN UXX
C...
C...        A + 4*B + 9*C + 16*D = 2                           (26)
C...
C...     TO DROP UXXX
C...
C...        A + 8*B + 27*C + 64*D = 0                          (27)
C...
C...     TO DROP UXXXX
C...
C...        A + 16*B + 81*C + 256*D = 0                        (28)
C...
C...     TO DROP UXXXXX
C...
C...        A + 32*B + 243*C + 1024*D = 0                      (29)
C...
C...  EQUATIONS (25) TO (29) CAN BE SOLVED FOR A, B, C, D AND E.  IF
C...
C...  EQUATION (26) IS SUBTRACTED FROM EQUATIONS (27), (28) AND (29),
C...
C...     4*B + 18*C + 48*D = -2                                (30)
C...
C...     12*B + 72*C + 240*D = -2                              (31)
C...
C...     28*B + 234*C + 1008*D = -2                            (32)
C...
C...  EQUATIONS (30), (31) AND (32) CAN BE SOLVED FOR B, C AND D
C...
C...     18*C + 96*D = 4                                       (33)
C...
C...     108*C + 672*D = 12                                    (34)
C...
C...  EQUATIONS (3) AND (34) CAN BE SOLVED FOR C AND D, C = 8/9, D =
C...  -1/8.
C...
C...  FROM EQUATION (30), B = -3.  FROM EQUATION (26), A = 8.  FROM
C...  EQUATION (25), E = -25/6.
C...
C...  THE FINAL DIFFERENTIATION FORMULA IS THEN OBTAINED AS
C...
C...  8*U(I+1) - 3*U(I+2) + (8/9)*U(I+3) - (1/8)*U(I+4)
C...
C...  - (25/6)*UX(I)*DX
C...
C...  = (8 - 3 + (8/9) - (1/8))*U(I) + UXX*(DX**2) + O(DX**6)
C...
C...  OR
C...
C...  UXX(I) = (1/12*DX**2)*((-415/6)*U(I) + 96*U(I+1) - 36*U(I+2)
C...                                                                (35)
C...  + (32/3)*U(I+3) - (3/2)*U(I+4) - 50*UX(I)*DX) + O(DX**4)
C...
C...  EQUATION (35) WILL BE APPLIED AT I = 1 AND I = N
C...
C...  UXX(1) = (1/12*DX**2)*((-415/6)*U(1) + 96*U(2) - 36*U(3)
C...                                                                (36)
C...  + (32/3)*U(4) - (3/2)*U(5) - 50*UX(1)*DX) + O(DX**4)
C...
C...  UXX(N) = (1/12*DX**2)*((-415/6)*U(N) + 96*U(N-1) - 36*U(N-2)
C...                                                                (37)
C...  + (32/3)*U(N-3) - (3/2)*U(N-4) + 50*UX(N)*DX) + O(DX**4)
C...
C...  CHECK THE VALUES OF A, B, C AND D IN EQUATIONS (25) TO (29)
C...  CALL CHECK(3)
C...
C...  ALTERNATIVELY, FOR GRID POINT 1, AN APPROXIMATION WITH A DIRICHLET
C...  BOUNDARY CONDITION OF THE FORM
C...
C...  A*U(I+1) + B*U(I+2) + C*U(I+3) + D*U(I+4) + E*U(I+5) + F*U(I)
C...
C...  CAN BE USED.  THE CORRESPONDING ALGEBRAIC EQUATIONS ARE
C...
C...     TO DROP UX
C...
C...        A + 2*B + 3*C + 4*D + 5*E = 0                      (38)
C...
C...     TO RETAIN UXX
C...
C...        A + 4*B + 9*C + 16*D + 25*E = 2                    (39)
C...
C...     TO DROP UXXX
C...
C...        A + 8*B + 27*C + 64*D + 125*E = 0                  (40)
C...
C...     TO DROP UXXXX
C...
C...        A + 16*B + 81*C + 256*D + 625*E = 0                (41)
C...
C...     TO DROP UXXXXX
C...
C...        A + 32*B + 243*C + 1024*D + 3125*E = 0             (42)
C...
C...  EQUATIONS (38), (39), (40), (41) AMD (42) CAN BE SOLVED FOR A,
C...  B, C, D AND E.
C...
C...        2*B + 6*C + 12*D + 20*E = 2                        (43)
C...
C...        6*B + 24*C + 60*D + 120*E = 0                      (44)
C...
C...        14*B + 78*C + 252*D + 620*E = 0                    (45)
C...
C...        30*B + 240*C + 1020*D + 3120*E = 0                 (46)
C...
C...  EQUATIONS (43), (44), (45) AND (46) CAN BE SOLVED FOR B, C, D
C...  AND E
C...
C...        6*C + 24*D + 60*E = -6                             (47)
C...
C...        36*C + 168*D + 480*E = -14                         (48)
C...
C...        150*C + 840*D + 2820*E = -30                       (49)
C...
C...  EQUATIONS (47), (48) AND (49) CAN BE SOLVED FOR C, D AND E
C...
C...        24*D + 120*E = 22                                  (50)
C...
C...        240*D + 1320*E = 120                               (51)
C...
C...  FROM EQUATIONS (50) AND (51), D = 61/12, E = -5/6.  FROM EQUATION
C...  (47), C = -13.  FROM EQUATION (43), B = 107/6.  FROM EQUATION
C...  (38), A = -77/6.
C...
C...  THE FINAL DIFFERENTIATION FORMULA IS THEN OBTAINED AS
C...
C...  (-77/6)*U(I+1) + (107/6)*U(I+2) - 13*U(I+3) + (61/12)*U(I+4)
C...
C...  - (5/6)*U(I+5) = (-77/6 + 107/6 - 13 + 61/12 - 5/6)*U(I) +
C...
C...  UXX(I)*(DX**2) + O(DX**6)
C...
C...  OR
C...
C...  UXX(I) = (1/12*DX**2)*(45*U(I) - 154*U(I+1) + 214*U(I+2)
C...                                                                (52)
C...         - 156*U(I+3) + 61*U(I+4) - 10*U(I+5)) + O(DX**4)
C...
C...  EQUATION (52) WILL BE APPLIED AT I = 1 AND I = N
C...
C...  UXX(1) = (1/12*DX**2)*(45*U(1) - 154*U(2) + 214*U(3)
C...                                                                (53)
C...         - 156*U(4) + 61*U(5) - 10*U(6)) + O(DX**4)
C...
C...  UXX(N) = (1/12*DX**2)*(45*U(N) - 154*U(N-1) + 214*U(N-2)
C...                                                                (54)
C...         -156*U(N-3) + 61*U(N-4) - 10*U(N-5)) + O(DX**4)
C...
C...  CHECK THE VALUES OF A, B, C, D, AND E IN EQUATIONS (38) TO (42)
C...  CALL CHECK(4)
C...
C...  ******************************************************************
C...
C...  GRID SPACING
#endif
      DX=(XU-XL)/DFLOAT(N-1)
!C...
!C...  1/(12*DX**2) FOR SUBSEQUENT USE
      R12DXS=1./(12.0D0*DX**2)
!C...
!C...  UXX AT THE LEFT BOUNDARY
!C...
!C...     WITHOUT UX (EQUATION (53))
         IF(NL.EQ.1)THEN
         UXX(1)=R12DXS*                &
                      (   45.0D0*U(1)  &
                        -154.0D0*U(2)  &
                        +214.0D0*U(3)  &
                        -156.0D0*U(4)  &
                         +61.0D0*U(5)  &
                         -10.0D0*U(6))
!C...
!C...     WITH UX (EQUATION (36))
         ELSE IF(NL.EQ.2)THEN
         UXX(1)=R12DXS*                     &
                      (-415.0D0/6.0D0*U(1)  &
                             +96.0D0*U(2)   &
                              -36.0D0*U(3)  &
                        +32.0D0/3.0D0*U(4)  &
                         -3.0D0/2.0D0*U(5)  &
                              -50.0D0*UX(1)*DX)
         END IF
!C...
!C...  UXX AT THE RIGHT BOUNDARY
!C...
!C...     WITHOUT UX (EQUATION (54))
         IF(NU.EQ.1)THEN
         UXX(N)=R12DXS*                   &
                      (   45.0D0*U(N  )   &
                        -154.0D0*U(N-1)   &
                        +214.0D0*U(N-2)   &
                        -156.0D0*U(N-3)   &
                         +61.0D0*U(N-4)   &
                         -10.0D0*U(N-5))
!C...
!C...     WITH UX (EQUATION (37))
         ELSE IF(NU.EQ.2)THEN
         UXX(N)=R12DXS*                       &
                      (-415.0D0/6.0D0*U(N  )  &
                              +96.0D0*U(N-1)  &
                              -36.0D0*U(N-2)  &
                        +32.0D0/3.0D0*U(N-3)  &
                         -3.0D0/2.0D0*U(N-4)  &
                              +50.0D0*UX(N  )*DX)
         END IF
!C...
!C...  UXX AT THE INTERIOR GRID POINTS
!C...
!C...     I = 2 (EQUATION (23))
         UXX(2)=R12DXS*                 &
                      (   10.0D0*U(1)   &
                         -15.0D0*U(2)   &
                          -4.0D0*U(3)   &
                         +14.0D0*U(4)   &
                          -6.0D0*U(5)   &
                          +1.0D0*U(6))
!C...
!C...     I = N-1 (EQUATION (24))
         UXX(N-1)=R12DXS*                 &
                        ( 10.0D0*U(N  )   &
                         -15.0D0*U(N-1)   &
                          -4.0D0*U(N-2)   &
                         +14.0D0*U(N-3)   &
                          -6.0D0*U(N-4)   &
                          +1.0D0*U(N-5))
!C...
!C...     I = 3, 4,..., N-2 (EQUATION (9))
           !dir$ assume_aligned UXX:64,U:64
           !dir$ vector aligned
           !dir$ ivdep
           !dir$ vector always
         DO 1 I=3,N-2
              UXX(I)=R12DXS*             &
                      (   -1.0D0*U(I-2)  &
                         +16.0D0*U(I-1)  &
                         -30.0D0*U(I  )  &
                         +16.0D0*U(I+1)  &
                          -1.0D0*U(I+2))
1        CONTINUE
      
      END SUBROUTINE


 

   
!=======================================================================!

   SUBROUTINE DDRIV2 (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK, &
        LENW,IWORK,LENIW,G)
#if 0
C***BEGIN PROLOGUE  DDRIV2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  871105   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of DDRIV2 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  DDRIV2 uses double precision arithmetic.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C  I.  ABSTRACT  .......................................................
C
C    The function of DDRIV2 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.
C    DDRIV2 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C       (REMEMBER--To run DDRIV2 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
C
C    The user should use parameter names in the call sequence of DDRIV2
C    for those quantities whose value may be altered by DDRIV2.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of differential equations.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F and
C             G.  Thus parameters required by F and G can be stored in
C             this array in components N+1 and above.  (Note: Changes
C             by the user to the first N components of this array will
C             take effect only after a restart, i.e., after setting
C             MSTATE to +1(-1).)
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C             This computes YDOT = F(Y,T), the right hand side of the
C             differential equations.  Here Y is a vector of length at
C             least N.  The actual length of Y is determined by the
C             user's declaration in the program which calls DDRIV2.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  DDRIV2.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV2, he should set N to zero.
C             DDRIV2 will signal this by returning a value of MSTATE
C             equal to +6(-6).  Altering the value of N in F has no
C             effect on the value of N in the call sequence of DDRIV2.
C
C    TOUT   = (Input) The point at which the solution is desired.
C
C    MSTATE = An integer describing the status of integration.  The user
C             must initialize MSTATE to +1 or -1.  If MSTATE is
C             positive, the routine will integrate past TOUT and
C             interpolate the solution.  This is the most efficient
C             mode.  If MSTATE is negative, the routine will adjust its
C             internal step to reach TOUT exactly (useful if a
C             singularity exists beyond TOUT.)  The meaning of the
C             magnitude of MSTATE:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of MSTATE should be tested by the
C                  user.  Unless DDRIV2 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, DDRIV2 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV2
C                  again.  Other than an error in problem setup, the
C                  most likely cause for this condition is trying to
C                  integrate a stiff set of equations with the non-stiff
C                  integrator option. (See description of MINT below.)
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV2
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling DDRIV2 again.
C               6  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               7  (Output)(Unsuccessful) N has been set to zero in
C                  FUNCTION G.  See description of G below.
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  DDRIV2 attempts
C             to find the value of T at which one of the equations
C             changes sign.  DDRIV2 can find at most one root per
C             equation per internal integration step, and will then
C             return the solution either at TOUT or at a root, whichever
C             occurs first in the direction of integration.  The index
C             of the equation whose root is being reported is stored in
C             the sixth element of IWORK.
C             NOTE: NROOT is never altered by this program.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  EPS = 0 is allowed.  On output, the adjusted
C             relative accuracy if the input value was too small.  The
C             value of EPS should be set as large as is reasonable,
C             because the amount of work done by DDRIV2 increases as
C             EPS decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest physically
C             meaningful value for the solution.  This is used inter-
C             nally to compute an array YWT(I) = MAX(ABS(Y(I)), EWT).
C             One step error estimates divided by YWT(I) are kept less
C             than EPS.  Setting EWT to zero provides pure relative
C             error control.  However, setting EWT smaller than
C             necessary can adversely affect the running time.
C
C    MINT   = (Input) The integration method flag.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.
C             MINT may not be changed without restarting, i.e., setting
C             the magnitude of MSTATE to 1.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The length of WORK should be at least
C               16*N + 2*NROOT + 204         if MINT is 1, or
C               N*N + 10*N + 2*NROOT + 204   if MINT is 2, or
C               N*N + 17*N + 2*NROOT + 204   if MINT is 3,
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to DDRIV2.
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               21      if MINT is 1, or
C               N+21    if MINT is 2 or 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to DDRIV2.
C
C    G      = A double precision FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to
C             obtain the value of each of the NROOT equations for which
C             a root is desired.  G is of the form:
C                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT)
C                   DOUBLE PRECISION Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls DDRIV2.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.  Normally a return from G passes
C             control back to  DDRIV2.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV2, he should set N to zero.
C             DDRIV2 will signal this by returning a value of MSTATE
C             equal to +7(-7).  In this case, the index of the equation
C             being evaluated is stored in the sixth element of IWORK.
C             Altering the value of N in G has no effect on the value of
C             N in the call sequence of DDRIV2.
C
C***LONG DESCRIPTION
C
C  III.  OTHER COMMUNICATION TO THE USER  ..............................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERROR.  That program will
C       terminate the user's run if it detects a probable problem setup
C       error, e.g., insufficient storage allocated by the user for the
C       WORK array.  Messages are written on the standard error message
C       file.  At installations which have this error handling package
C       the user should determine the standard error handling file from
C       the local documentation.  Otherwise the short but serviceable
C       routine, XERROR, available with this package, can be used.  That
C       program writes on logical unit 6 to transmit messages.  A
C       complete description of XERROR is given in the Sandia
C       Laboratories report SAND78-1189 by R. E. Jones.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken since last initialization.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. On any return from DDRIV2 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    B. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV2.
C
C    C. When the routine G is not required, difficulties associated with
C       an unsatisfied external can be avoided by using the name of the
C       routine which calculates the right hand side of the differential
C       equations in place of G in the call sequence of DDRIV2.
C
C  V.  USAGE  ..........................................................
C
C               PROGRAM SAMPLE
C               EXTERNAL F
C               PARAMETER(MINT = 1, NROOT = 0, N = ...,
C              8          LENW = 16*N + 2*NROOT + 204, LENIW = 21)
C                                           N is the number of equations
C               DOUBLE PRECISION EPS, EWT, T, TOUT, WORK(LENW), Y(N)
C               INTEGER IWORK(LENIW)
C               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C               T = 0.                           Initial point
C               DO 10 I = 1,N
C          10     Y(I) = ...                     Set initial conditions
C               TOUT = T
C               EWT = ...
C               MSTATE = 1
C               EPS = ...
C          20   CALL DDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
C              8             MINT, WORK, LENW, IWORK, LENIW, F)
C                                          Last argument is not the same
C                                          as F if rootfinding is used.
C               IF (MSTATE .GT. 2) STOP
C               WRITE(6, 100) TOUT, (Y(I), I=1,N)
C               TOUT = TOUT + 1.
C               IF (TOUT .LE. 10.) GO TO 20
C          100  FORMAT(...)
C               END (Sample)
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDRIV3,XERROR
C***END PROLOGUE  DDRIV2
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
#endif
      EXTERNAL F, G
      DOUBLE PRECISION EPS, EWT, EWTCOM(1), G, HMAX, T, TOUT, &
          WORK(*), Y(*)
      INTEGER IWORK(*)
      CHARACTER MSG*81
      PARAMETER(IMPL = 0, MXSTEP = 1000)
!C***FIRST EXECUTABLE STATEMENT  DDRIV2
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(MSG, '(''DDRIV21FE Illegal input. Improper value for '', ''the integration method flag,'', I8)') MINT
        CALL XERROR(MSG(1:81), 81, 21, 2)
        RETURN
      END IF
      IF (MSTATE .GE. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      EWTCOM(1) = EWT
      IF (EWT .NE. 0.D0) THEN
        IERROR = 3
      ELSE
        IERROR = 2
      END IF
      IF (MINT .EQ. 1) THEN
        MITER = 0
        MXORD = 12
      ELSE IF (MINT .EQ. 2) THEN
        MITER = 2
        MXORD = 5
      ELSE IF (MINT .EQ. 3) THEN
        MITER = 2
        MXORD = 12
      END IF
      HMAX = 2.D0*ABS(TOUT - T)
      CALL DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,   &
                 IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,  &
                  LENW, IWORK, LENIW, F, F, NDE, MXSTEP, G, F)
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
    END SUBROUTINE
    
      SUBROUTINE DDRIV3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR, &
        MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,  &
        FA,NDE,MXSTEP,G,USERS)
!C***BEGIN PROLOGUE  DDRIV3
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C   From the book "Numerical Methods and Software"
!C      by D. Kahaner, C. Moler, S. Nash
!C         Prentice Hall 1988
!C***END PROLOGUE  DDRIV3
      EXTERNAL F, JACOBN, FA, G, USERS
      DOUBLE PRECISION AE, BIG, EPS, EWT(*), G, GLAST, H, HMAX, HSIGN,  &
          NROUND, RE, D1MACH, SIZE, DNRM2, SUM, T, TLAST, TOUT, TROOT, &
          UROUND, WORK(*), Y(*)
      INTEGER IWORK(*)
      LOGICAL CONVRG
      CHARACTER MSG*205
      PARAMETER(NROUND = 20.D0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,                        &
               IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,              &
               IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,           &
               ITOUT = 167, ITQ = 168, ITREND = 204, IYH = 205,          &
               INDMXR = 1, INQUSD = 2, INSTEP = 3, INFE = 4, INJE = 5,   &
               INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,           &
               IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,          &
               INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,         &
               IMTR = 18, IMXRDS = 19, IMXORD = 20, INDPRT = 21,         &
               INDPVT = 22)
!C***FIRST EXECUTABLE STATEMENT  DDRIV3
      NPAR = N
      UROUND = D1MACH (4)
      IF (NROOT .NE. 0) THEN
        AE = D1MACH(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.D0) THEN
        WRITE(MSG, '(''DDRIV36FE Illegal input.  EPS,'', D16.8, '', is negative.'')') EPS
        CALL XERROR(MSG(1:60), 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(MSG, '(''DDRIV37FE Illegal input.  Number of equations,'', I8, '', is not positive.'')') N
        CALL XERROR(MSG(1:72), 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(MSG, '(''DDRIV314FE Illegal input.  Maximum order,'', I8, '', is not positive.'')') MXORD
      
        CALL XERROR(MSG(1:67), 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.           &
       (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))                 &
       .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.                         &
       (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.             &
       ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.           &
       (IMPL .EQ. 2 .AND. MINT .EQ. 1) .OR.                               &
       (NSTATE .LT. 1 .OR. NSTATE .GT. 10)) THEN
        WRITE(MSG, '(''DDRIV39FE Illegal input.  Improper value for '',''NSTATE(MSTATE), MINT, MITER or IMPL.'')')
       
        CALL XERROR(MSG(1:81), 81, 9, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR. &
       MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
        WRITE(MSG, '(''DDRIV310FE Illegal input.  Insufficient '',   &
       ''storage allocated for the IWORK array.  Based on the '')')
        WRITE(MSG(94:), '(''value of the input parameters involved, '', &
       ''the required storage is'', I8)') LIWCHK
        CALL XERROR(MSG(1:164), 164, 10, 2)
        RETURN
      END IF
!C                                                Allocate the WORK array
!C                                         IYH is the index of YH in WORK
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N
!C                                             IDFDY is the index of DFDY
!C
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3)  THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2)  THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5)  THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF
!C                                               IYWT is the index of YWT
      ISAVE1 = IYWT + N
!C                                           ISAVE1 is the index of SAVE1
      ISAVE2 = ISAVE1 + N
!C                                           ISAVE2 is the index of SAVE2
      IGNOW = ISAVE2 + N
!C                                             IGNOW is the index of GNOW
      ITROOT = IGNOW + NROOT
!C                                           ITROOT is the index of TROOT
      IFAC = ITROOT + NROOT
!C                                               IFAC is the index of FAC
      IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. MINT .EQ. 3) THEN
        IA = IFAC + N
      ELSE
        IA = IFAC
      END IF
!C                                                   IA is the index of A
      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN
        LENCHK = IA - 1
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*N
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*N
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN
        LENCHK = IA - 1 + N
      END IF
      IF (LENW .LT. LENCHK) THEN
        WRITE(MSG, '(''DDRIV38FE Illegal input.  Insufficient '',    &
       ''storage allocated for the WORK array.  Based on the '')')
        WRITE(MSG(92:), '(''value of the input parameters involved, '', &
       ''the required storage is'', I8)') LENCHK
        CALL XERROR(MSG(1:162), 162, 8, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        MATDIM = 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        MATDIM = N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MATDIM = 2*ML + MU + 1
      END IF
      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN
        NDECOM = N
      ELSE IF (IMPL .EQ. 2) THEN
        NDECOM = NDE
      END IF
      IF (NSTATE .EQ. 1) THEN
!C                                                  Initialize parameters
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
          IWORK(IMXORD) = MIN(MXORD, 12)
        ELSE IF (MINT .EQ. 2) THEN
          IWORK(IMXORD) = MIN(MXORD, 5)
        END IF
        IWORK(IMXRDS) = MXORD
        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN
          IWORK(IMNT) = MINT
          IWORK(IMTR) = MITER
          IWORK(IMNTLD) = MINT
          IWORK(IMTRLD) = MITER
        ELSE IF (MINT .EQ. 3) THEN
          IWORK(IMNT) = 1
          IWORK(IMTR) = 0
          IWORK(IMNTLD) = IWORK(IMNT)
          IWORK(IMTRLD) = IWORK(IMTR)
          IWORK(IMTRSV) = MITER
        END IF
        WORK(IHMAX) = HMAX
        H = (TOUT - T)*(1.D0 - 4.D0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.D0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.D0
        WORK(IHUSED) =0.D0
        WORK(IAVGRD) = 0.D0
        IWORK(INDMXR) = 0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        IWORK(INROOT) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
!C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
        IF (T .EQ. TOUT) RETURN
        GO TO 180
      END IF
!C                                             On a continuation, check
!C                                             that output points have
!C                                             been or will be overtaken.
      IF (IWORK(ICNVRG) .EQ. 1) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      END IF
      T = WORK(IT)
      H = WORK(IH)
      HSIGN = WORK(IHSIGN)
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180
!C
!C                                   IWORK(IJROOT) flags unreported
!C                                   roots, and is set to the value of
!C                                   NTASK when a root was last selected.
!C                                   It is set to zero when all roots
!C                                   have been reported.  IWORK(INROOT)
!C                                   contains the index and WORK(ITOUT)
!C                                   contains the value of the root last
!C                                   selected to be reported.
!C                                   IWORK(INRTLD) contains the value of
!C                                   NROOT and IWORK(INDTRT) contains
!C                                   the value of ITROOT when the array
!C                                   of roots was last calculated.
      IF (NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
!C                                      TOUT has just been reported.
!C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
              T = TROOT
              NSTATE = 5
              GO TO 580
            END IF
!C                                         A root has just been reported.
!C                                         Select the next root.
          ELSE
            TROOT = T
            IROOT = 0
            DO 50 I = 1,IWORK(INRTLD)
              JTROOT = IWORK(INDTRT) + I - 1
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN
!C
!C                                              Check for multiple roots.
!C
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND. &
               I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                  GO TO 60
                END IF
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = TROOT
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
                  NSTATE = 5
                  T = TROOT
                  GO TO 580
                END IF
              END IF
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN
!C
!C                                     If there are no more roots, or the
!C                                     user has altered TOUT to be less
!C                                     than a root, set IJROOT to zero.
!C
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN
                IWORK(IJROOT) = 0
              ELSE
                CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
                NSTATE = 5
                T = TROOT
                GO TO 580
              END IF
            END IF
          END IF
        END IF
      END IF
!C
      IF (NTASK .EQ. 1) THEN
        NSTATE = 2
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
!C                                                      Check if TOUT has
!C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''DDRIV32WRN With NTASK='', I1, '' on input, '',  &
         ''T,'', D16.8, '', was beyond TOUT,'', D16.8, ''.  Solution'', &
         '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG(1:124), 124, 2, 0)
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          NSTATE = 2
          GO TO 580
        END IF
!C                                   Determine if TOUT has been overtaken
!C
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          NSTATE = 2
          GO TO 560
        END IF
!C                                             If there are no more roots
!C                                             to report, report T.
        IF (NSTATE .EQ. 5) THEN
          NSTATE = 2
          GO TO 560
        END IF
        NSTATE = 2
!C                                                       See if TOUT will
!C                                                       be overtaken.
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''DDRIV32WRN With NTASK='', I1, '' on input, '',   &
         ''T,'', D16.8, '', was beyond TOUT,'', D16.8, ''.  Solution'',  &
         '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG(1:124), 124, 2, 0)
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
!C                         Implement changes in MINT, MITER, and/or HMAX.
!C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND. &
       MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. WORK(IHMAX)) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. WORK(IH)) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = H
        END IF
        WORK(IHMAX) = HMAX
      END IF

 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
        JYH = IYH + I - 1
 190    Y(I) = WORK(JYH)
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          JGNOW = IGNOW + I - 1
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
 200     CONTINUE
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
 230      WORK(JYWT) = 1.D0
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
          JYWT = I + IYWT - 1
 250      WORK(JYWT) = EWT(I)
        GO TO 410
      END IF
!C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. 0.D0) GO TO 290
          JYWT = I + IYWT - 1
 280      WORK(JYWT) = ABS(Y(I))
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (NPAR, T, Y, WORK(ISAVE2))
          IF (NPAR .EQ. 0) THEN
            NSTATE = 6
            RETURN
          END IF
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1),      &
                      WORK(ISAVE2), T, H, WORK(IEL), IMPL, NPAR,   &
                      NDECOM, IFLAG)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (NPAR, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),  &
                         WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 9
              RETURN
            END IF
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF (WORK(JA) .EQ. 0.D0) GO TO 690
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. 0.D0) THEN
            WORK(JYWT) = ABS(Y(J))
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              JSAVE2 = J + ISAVE2 - 1
              WORK(JYWT) = ABS(H*WORK(JSAVE2))
            ELSE
              JHYP = J + IYH + N - 1
              WORK(JYWT) = ABS(WORK(JHYP))
            END IF
          END IF
          IF (WORK(JYWT) .EQ. 0.D0) WORK(JYWT) = UROUND
 360      CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
          JYWT = I + IYWT - 1
 380      WORK(JYWT) = MAX(EWT(1), ABS(Y(I)))
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
          JYWT = I + IYWT - 1
 400      WORK(JYWT) = MAX(EWT(I), ABS(Y(I)))
      END IF
!C
 410  DO 420 I = 1,N
        JYWT = I + IYWT - 1
        JSAVE2 = I + ISAVE2 - 1
 420    WORK(JSAVE2) = Y(I)/WORK(JYWT)
      SUM = DNRM2(N, WORK(ISAVE2), 1)/SQRT(DBLE(N))
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.D0 + 10.D0*UROUND)
        WRITE(MSG, '(''DDRIV34REC At T,'', D16.8, '', the requested '',   &
       ''accuracy, EPS, was not obtainable with the machine '',          &
       ''precision.  EPS has been increased to'')') T
        WRITE(MSG(137:), '(D16.8)') EPS
        CALL XERROR(MSG(1:152), 152, 4, 1)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        WRITE(MSG, '(''DDRIV35WRN At T,'', D16.8, '', the step size,'',  &
       D16.8, '', is smaller than the roundoff level of T.  '')') T, H
        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',  &
       ''change in the right hand side of the differential '',        &
       ''equations.'')')
        CALL XERROR(MSG(1:205), 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
          WRITE(MSG, '(''DDRIV33WRN At T,'', D16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', D16.8)')
     8    T, MXSTEP, TOUT
          CALL XERROR(MSG(1:103), 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL DDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND, USERS,
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV)
C
      CALL DDSTP (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, NPAR,
     8           NDECOM, WORK(IYWT), UROUND, USERS,  WORK(IAVGH),
     8           WORK(IAVGRD), WORK(IH), WORK(IHUSED), IWORK(IJTASK),
     8           IWORK(IMNTLD), IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH),
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IFAC),
     8            WORK(IHOLD), IWORK(INDPVT), JSTATE, IWORK(INQ),
     8            IWORK(INWAIT), WORK(IRC), WORK(IRMAX), WORK(ISAVE1),
     8            WORK(ISAVE2), WORK(ITQ), WORK(ITREND), MINT,
     8            IWORK(IMTRSV), IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      GO TO (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
          IF (GLAST*WORK(JGNOW) .GT. 0.D0) THEN
            WORK(JTROOT) = T + H
          ELSE
            IF (WORK(JGNOW) .EQ. 0.D0) THEN
              WORK(JTROOT) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.D0) THEN
                WORK(JTROOT) = T + H
              ELSE
                IF (ABS(WORK(IHUSED)) .GE. UROUND*ABS(T)) THEN
                  TLAST = T - WORK(IHUSED)
                  IROOT = I
                  TROOT = T
                  CALL DDZRO (AE, G, H, NPAR, IWORK(INQ), IROOT, RE, T,
     8                        WORK(IYH), UROUND,  TROOT, TLAST,
     8                        WORK(JGNOW), GLAST,  Y)
                  DO 480 J = 1,N
  480               Y(J) = WORK(IYH + J -1)
                  IF (NPAR .EQ. 0) THEN
                    IWORK(INROOT) = I
                    NSTATE = 7
                    RETURN
                  END IF
                  WORK(JTROOT) = TROOT
                ELSE
                  WORK(JTROOT) = T
                  IROOT = I
                END IF
              END IF
            END IF
          END IF
 500      CONTINUE
        IF (IROOT .EQ. 0) THEN
          IWORK(IJROOT) = 0
!C                                                  Select the first root
        ELSE
          IWORK(IJROOT) = NTASK
          IWORK(INRTLD) = NROOT
          IWORK(INDTRT) = ITROOT
          TROOT = T + H
          DO 510 I = 1,NROOT
            JTROOT = ITROOT + I - 1
            IF (WORK(JTROOT)*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = WORK(JTROOT)
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = TROOT
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
            NSTATE = 5
            T = TROOT
            GO TO 580
          END IF
        END IF
      END IF
!C                               Test for NTASK condition to be satisfied
      NSTATE = 2
      IF (NTASK .EQ. 1) THEN
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260
        CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        GO TO 580
!C                               TOUT is assumed to have been attained
!C                               exactly if T is within twenty roundoff
!C                               units of TOUT, relative to max(TOUT, T).
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
          GO TO 260
        END IF
      END IF
!C                                      All returns are made through this
!C                                      section.  IMXERR is determined.
 560  DO 570 I = 1,N
        JYH = I + IYH - 1
 570    Y(I) = WORK(JYH)
 580  IF (CONVRG) THEN
        IWORK(ICNVRG) = 1
      ELSE
        IWORK(ICNVRG) = 0
      END IF
      IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.D0
      IMXERR = 1
      IWORK(INDMXR) = IMXERR
      DO  590 I = 1,N
!C                                            SIZE = ABS(ERROR(I)/YWT(I))
        JYWT = I + IYWT - 1
        JERROR = I + ISAVE1 - 1
        SIZE = ABS(WORK(JERROR)/WORK(JYWT))
        IF (BIG .LT. SIZE) THEN
          BIG = SIZE
          IMXERR = I
          IWORK(INDMXR) = IMXERR
        END IF
 590    CONTINUE
      RETURN
!C
 660  NSTATE = JSTATE
      RETURN
!C                                        Fatal errors are processed here
!C
 670  WRITE(MSG, '(''DDRIV311FE At T,'', D16.8, '', the attempted '',  &
       ''step size has gone to zero.  Often this occurs if the '',    &
       ''problem setup is incorrect.'')') T
      CALL XERROR(MSG(1:129), 129, 11, 2)
      RETURN
!C
 680  WRITE(MSG, '(''DDRIV312FE At T,'', D16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL XERROR(MSG(1:165), 165, 12, 2)
      RETURN
!C
 690  WRITE(MSG, '(''DDRIV313FE At T,'', D16.8, '', while solving'',  &
       '' A*YDOT = F, A is singular.'')') T
      CALL XERROR(MSG(1:74), 74, 13, 2)
      RETURN
   END SUBROUTINE
   
      SUBROUTINE DDNTP (H,K,N,NQ,T,TOUT,YH,Y)
!C***BEGIN PROLOGUE  DDNTP
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***END PROLOGUE  DDNTP
      DOUBLE PRECISION FACTOR, H, R, T, TOUT, Y(*), YH(N,*)
!C***FIRST EXECUTABLE STATEMENT  DDNTP
      IF (K .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        KUSED = MIN(K, NQ)
        FACTOR = 1.D0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*DBLE(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        DO 80 JJ = KUSED+1,NQ
          J = K + 1 + NQ - JJ
          FACTOR = 1.D0
          DO 60 KK = 1,KUSED
 60         FACTOR = FACTOR*DBLE(J-KK)
          DO 70 I = 1,N
 70         Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80       CONTINUE
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
   END SUBROUTINE
   
      SUBROUTINE DDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)
!C***BEGIN PROLOGUE  DDZRO
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***ROUTINES CALLED  DDNTP
!C***END PROLOGUE  DDZRO
      DOUBLE PRECISION A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC, &
          H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)
!C***FIRST EXECUTABLE STATEMENT  DDZRO
      ER = 4.D0*UROUND
      RW = MAX(RE, ER)
      IC = 0
      ACBS = ABS(B - C)
      A = C
      FA = FC
      KOUNT = 0
!C                                                    Perform interchange
 10   IF (ABS(FC) .LT. ABS(FB)) THEN
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
      END IF
      CMB = 0.5D0*(C - B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AE
!C                                                Test stopping criterion
      IF (ACMB .LE. TOL) RETURN
      IF (KOUNT .GT. 50) RETURN
!C                                    Calculate new iterate implicitly as
!C                                    B + P/Q, where we arrange P .GE. 0.
!C                         The implicit form is used to prevent overflow.
      P = (B - A)*FB
      Q = FA - FB
      IF (P .LT. 0.D0) THEN
        P = -P
        Q = -Q
      END IF
!C                          Update A and check for satisfactory reduction
!C                          in the size of our bounding interval.
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.D0*ACMB .GE. ACBS) THEN
!C                                                                 Bisect
          B = 0.5D0*(C + B)
          GO TO 20
        END IF
        IC = 0
      END IF
      ACBS = ACMB
!C                                            Test for too small a change
      IF (P .LE. ABS(Q)*TOL) THEN
!C                                                 Increment by tolerance
        B = B + SIGN(TOL, CMB)
!C                                               Root ought to be between
!C                                               B and (C + B)/2.
      ELSE IF (P .LT. CMB*Q) THEN
!C                                                            Interpolate
        B = B + P/Q
      ELSE
!C                                                                 Bisect
        B = 0.5D0*(C + B)
      END IF
!C                                             Have completed computation
!C                                             for new iterate B.
 20   CALL DDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (N .EQ. 0) RETURN
      IF (FB .EQ. 0.D0) RETURN
      KOUNT = KOUNT + 1
!C
!C             Decide whether next step is interpolation or extrapolation
!C
      IF (SIGN(1.0D0, FB) .EQ. SIGN(1.0D0, FC)) THEN
        C = A
        FC = FA
      END IF
      GO TO 10
    END SUBROUTINE DDZRO
    
    SUBROUTINE DDSTP (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,      &
        MITER,ML,MU,N,NDE,YWT,UROUND,USERS,AVGH,AVGORD,H,HUSED,JTASK,   &
        MNTOLD,MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,FAC, &
        HOLD,IPVT,JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,  &
        MTRSV,MXRDSV)
!C***BEGIN PROLOGUE  DDSTP
!!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***ROUTINES CALLED  DDNTL,DDPST,DDCOR,DDPSC,DDSCL,DNRM2
!C***END PROLOGUE  DDSTP
      EXTERNAL F, JACOBN, FA, USERS
      DOUBLE PRECISION A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,  &
          BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,    &
          ERDN, ERUP, ETEST, FAC(*), H, HMAX, HN, HOLD, HS, HUSED,     &
          NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM,  &
          SAVE1(*), SAVE2(*), DNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD, &
          UROUND, Y(*), YH(N,*), YWT(*), Y0NRM
      INTEGER IPVT(*)
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3D0, BIAS2 = 1.2D0, BIAS3 = 1.4D0, MXFAIL = 3, &
               MXITER = 3, MXTRY = 50, RCTEST = .3D0, RMFAIL = 2.D0,     &
               RMNORM = 10.D0, TRSHLD = 1.D0)
      DATA IER /.FALSE./
!C***FIRST EXECUTABLE STATEMENT  DDSTP
      NSV = N
      BND = 0.D0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,            &
                   MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,          &
                   UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,     &
                   YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,      &
                   RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
      END IF
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL DDPSC (1, N, NQ,  YH)
      EVALJC = ((ABS(RC - 1.D0) .GT. RCTEST) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
!C
 110  ITER = 0
      DO 115 I = 1,N
 115    Y(I) = YH(I,1)
      CALL F (N, T, Y, SAVE2)
      IF (N .EQ. 0) THEN
        JSTATE = 6
        GO TO 430
      END IF
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,              &
                   MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND,       &
                   NFE, NJE,  A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG,          &
                   BND, JSTATE)
        IF (N .EQ. 0) GO TO 430
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.D0
      END IF
      DO 125 I = 1,N
 125    SAVE1(I) = 0.D0
!C                      Up to MXITER corrector iterations are taken.
!C                      Convergence is tested by requiring the r.m.s.
!C                      norm of changes to be less than EPS.  The sum of
!C                      the corrections is accumulated in the vector
!C                      SAVE1(I).  It is approximately equal to the L-th
!C                      derivative of Y multiplied by
!C                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus
!C                      proportional to the actual errors to the lowest
!C                      power of H present (H**L).  The YH array is not
!C                      altered in the correction loop.  The norm of the
!C                      iterate difference is stored in D.  If
!C                      ITER .GT. 0, an estimate of the convergence rate
!C                      constant is stored in TREND, and this is used in
!C                      the convergence test.
!C
 130  CALL DDCOR (DFDY, EL, FA, H, IMPL, IPVT, MATDIM, MITER, ML,        &
                 MU, N, NDE, NQ, T, USERS, Y, YH, YWT,  EVALFA, SAVE1,  &
                 SAVE2,  A, D, JSTATE)
        IF (N .EQ. 0) GO TO 430
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = DNRM2(N, SAVE1, 1)
          DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
          Y0NRM = DNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          NUMER = DNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.D0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
          IF (DENOM .NE. 0.D0)  &
         BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9D0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.D0*TREND, 1.D0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          GO TO 430
        END IF
        NFE = NFE + 1
        GO TO 130
      END IF
!C                     The corrector iteration failed to converge in
!C                     MXITER tries.  If partials are involved but are
!C                     not up to date, they are reevaluated for the next
!C                     try.  Otherwise the YH array is retracted to its
!C                     values before prediction, and H is reduced, if
!C                     possible.  If not, a no-convergence exit is taken.
      IF (CONVRG) THEN
        EVALJC = .TRUE.
        EVALFA = .FALSE.
        GO TO 110
      END IF
 160  T = TOLD
      CALL DDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3D0
      ELSE
        RH = .9D0*(EPS/CTEST)**(.2D0)
      END IF
      IF (RH*H .EQ. 0.D0) GO TO 400
      CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
!C                          The corrector has converged.  CONVRG is set
!C                          to .TRUE. if partial derivatives were used,
!C                          to indicate that they may need updating on
!C                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER .NE. 0)
      DO 180 I = 1,NDE
 180    SAVE2(I) = SAVE1(I)/YWT(I)
      ETEST = DNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(DBLE(NDE)))
!C
!C                           The error test failed.  NFAIL keeps track of
!C                           multiple failures.  Restore T and the YH
!C                           array to their previous values, and prepare
!C                           to try the step again.  Compute the optimum
!C                           step size for this or one lower order.
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL DDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.D0/(BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
          IF (NQ .GT. 1) THEN
            DO 190 I = 1,NDE
 190          SAVE2(I) = YH(I,NQ+1)/YWT(I)
            ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
            RH1 = 1.D0/MAX(1.D0, BIAS1*(ERDN/EPS)**(1.D0/DBLE(NQ)))
            IF (RH2 .LT. RH1) THEN
              NQ = NQ - 1
              RC = RC*EL(1,NQ)/EL(1,NQ+1)
              RH = RH1
            ELSE
              RH = RH2
            END IF
          ELSE
            RH = RH2
          END IF
          NWAIT = NQ + 2
          IF (RH*H .EQ. 0.D0) GO TO 400
          CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          GO TO 100
        END IF
!C                Control reaches this section if the error test has
!C                failed MXFAIL or more times.  It is assumed that the
!C                derivatives that have accumulated in the YH array have
!C                errors of the wrong order.  Hence the first derivative
!C                is recomputed, the order is set to 1, and the step is
!C                retried.
        NFAIL = 0
        JTASK = 2
        DO 215 I = 1,N
 215      Y(I) = YH(I,1)
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,          &
                   MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,        &
                   UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,   &
                   YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,    &
                   RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        RMAX = RMNORM
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
!C                          After a successful step, update the YH array.
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = (DBLE(NSTEP-1)*AVGH + H)/DBLE(NSTEP)
      AVGORD = (DBLE(NSTEP-1)*AVGORD + DBLE(NQ))/DBLE(NSTEP)
      DO 230 J = 1,NQ+1
        DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
      DO 235 I = 1,N
 235    Y(I) = YH(I,1)
!C                                          If ISWFLG is 3, consider
!C                                          changing integration methods.
!C
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.D0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/DBLE(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,  &
           (ETEST/(EPS*EL(NQ+1,1)))**(1.D0/DBLE(NQ+1)))
            IF (HS .GT. 1.2D0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.D0
              RMAX = RMNORM
              TREND = 1.D0
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/DBLE(NQ+1)))
            HN = ABS(H)/MAX(UROUND, &
           (ETEST*EL(NQ+1,1)/EPS)**(1.D0/DBLE(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.D0
              CONVRG = .FALSE.
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          END IF
        END IF
      END IF
      IF (SWITCH) THEN
        MINT = 2
        MNTOLD = MINT
        MITER = MTRSV
        MTROLD = MITER
        MAXORD = MIN(MXRDSV, 5)
        NQ = MIN(NQ, MAXORD)
        RC = 0.D0
        RMAX = RMNORM
        TREND = 1.D0
        CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
!C                           Consider changing H if NWAIT = 1.  Otherwise
!C                           decrease NWAIT by 1.  If NWAIT is then 1 and
!C                           NQ.LT.MAXORD, then SAVE1 is saved for use in
!C                           a possible order increase on the next step.
!C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
        IF (RH.GT.TRSHLD) CALL DDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
      ELSE IF (NWAIT .GT. 1) THEN
        NWAIT = NWAIT - 1
        IF (NWAIT .EQ. 1 .AND. NQ .LT. MAXORD) THEN
          DO 250 I = 1,NDE
 250        YH(I,MAXORD+1) = SAVE1(I)
        END IF
!C             If a change in H is considered, an increase or decrease in
!C             order by one is considered also.  A change in H is made
!C             only if it is by a factor of at least TRSHLD.  Factors
!C             RH1, RH2, and RH3 are computed, by which H could be
!C             multiplied at order NQ - 1, order NQ, or order NQ + 1,
!C             respectively.  The largest of these is determined and the
!C             new order chosen accordingly.  If the order is to be
!C             increased, we compute one additional scaled derivative.
!C             If there is a change of order, reset NQ and the
!C             coefficients.  In any case H is reset according to RH and
!C             the YH array is rescaled.
      ELSE
        IF (NQ .EQ. 1) THEN
          RH1 = 0.D0
        ELSE
          DO 270 I = 1,NDE
 270        SAVE2(I) = YH(I,NQ+1)/YWT(I)
          ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
          RH1 = 1.D0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.D0/DBLE(NQ)))
        END IF
        RH2 = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.D0
        ELSE
          DO 290 I = 1,NDE
 290        SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
          ERUP = DNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(DBLE(NDE)))
          RH3 = 1.D0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.D0/DBLE(NQ+2)))
        END IF
        IF (RH1 .GT. RH2 .AND. RH1 .GE. RH3) THEN
          RH = RH1
          IF (RH .LE. TRSHLD) GO TO 380
          NQ = NQ - 1
          RC = RC*EL(1,NQ)/EL(1,NQ+1)
        ELSE IF (RH2 .GE. RH1 .AND. RH2 .GE. RH3) THEN
          RH = RH2
          IF (RH .LE. TRSHLD) GO TO 380
        ELSE
          RH = RH3
          IF (RH .LE. TRSHLD) GO TO 380
          DO 360 I = 1,N
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/DBLE(NQ+1)
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.D0) RH = MIN(RH, 1.D0/(2.D0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
        RMAX = RMNORM
 380    NWAIT = NQ + 2
      END IF
!C               All returns are made through this section.  H is saved
!C               in HOLD to allow the caller to change H on the next step
      JSTATE = 1
      HOLD = H
      RETURN
!C
 400  JSTATE = 2
      HOLD = H
      DO 405 I = 1,N
 405    Y(I) = YH(I,1)
      RETURN
!C
 410  JSTATE = 3
      HOLD = H
      RETURN
!C
 420  JSTATE = 4
      HOLD = H
      RETURN
!C
 430  T = TOLD
      CALL DDPSC (-1, NSV, NQ,  YH)
      DO 435 I = 1,NSV
 435    Y(I) = YH(I,1)
 440  HOLD = H
      RETURN
   END SUBROUTINE
   
      SUBROUTINE DDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,   &
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,UROUND,USERS,Y,YWT,H,MNTOLD,   &
     8   MTROLD,NFE,RC,YH,A,CONVRG,EL,FAC,IER,IPVT,NQ,NWAIT,RH,RMAX,   &
     8   SAVE2,TQ,TREND,ISWFLG,JSTATE)
!C***BEGIN PROLOGUE  DDNTL
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***ROUTINES CALLED  DDCST,DDSCL,DGEFA,DGESL,DGBFA,DGBSL,DNRM2
!C***END PROLOGUE  DDNTL
      DOUBLE PRECISION A(MATDIM,*), EL(13,12), EPS, FAC(*), H, HMAX,      &
          HOLD, OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX,   &
          SMIN, DNRM2, SUM, SUM0, T, TQ(3,12), TREND, UROUND, Y(*),      &
          YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.D0)
!C***FIRST EXECUTABLE STATEMENT  DDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) THEN
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RMAX = RMINIT
        END IF
        RC = 0.D0
        CONVRG = .FALSE.
        TREND = 1.D0
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          RETURN
        END IF
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N, &
                       NDE, IFLAG)
            IF (N .EQ. 0) THEN
              JSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
            DO 150 I = 1,NDE
              IF (A(I,1) .EQ. 0.D0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = 0.D0
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/YWT(I)
        SUM = DNRM2(NDE, SAVE1, 1)
        SUM0 = 1.D0/MAX(1.D0, ABS(T))
        SMAX = MAX(SUM0, SUM)
        SMIN = MIN(SUM0, SUM)
        SUM = SMAX*SQRT(1.D0 + (SMIN/SMAX)**2)/SQRT(DBLE(NDE))
        H = SIGN(MIN(2.D0*EPS/SUM, ABS(H)), H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
        IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. ISWFLG .EQ. 3) THEN
          DO 20 I = 1,N
 20         FAC(I) = SQRT(UROUND)
        END IF
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.D0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          RH = H/HOLD
          CALL DDSCL (HMAX, N, NQ, RMAX,  HOLD, RC, RH, YH)
        END IF
      END IF
   END SUBROUTINE
   
      SUBROUTINE DDPST (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,   &
        NQ,SAVE2,T,USERS,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,FAC,IER,IPVT,    &
        SAVE1,ISWFLG,BND,JSTATE)
!C***BEGIN PROLOGUE  DDPST
!C***REFER TO  DDRIV3
!C  Subroutine DDPST is called to reevaluate the partials.
!C  If MITER is 1, 2, 4, or 5, the matrix
!C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
!C  decomposition, with the results also stored in DFDY.
!C***ROUTINES CALLED  DGEFA,DGBFA,DNRM2
!C***DATE WRITTEN   790601   (YYMMDD)
!C***REVISION DATE  870401   (YYMMDD)
!C***CATEGORY NO.  I1A2,I1A1B
!C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
!C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
!C***END PROLOGUE  DDPST
      DOUBLE PRECISION A(MATDIM,*), BL, BND, BP, BR, BU, DFDY(MATDIM,*),  & 
          DFDYMX, DIFF, DY, EL(13,12), FAC(*), FACMAX, FACMIN, FACTOR,   &
          H, SAVE1(*), SAVE2(*), SCALE, DNRM2, T, UROUND, Y(*),          &
          YH(N,*), YJ, YS, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER
      PARAMETER(FACMAX = .5D0)
!C***FIRST EXECUTABLE STATEMENT  DDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BU = UROUND**(.25D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          DO 170 J = 1,N
            YS = MAX(ABS(YWT(J)), ABS(Y(J)))
 120        DY = FAC(J)*YS
            IF (DY .EQ. 0.D0) THEN
              IF (FAC(J) .LT. FACMAX) THEN
                FAC(J) = MIN(100.D0*FAC(J), FACMAX)
                GO TO 120
              ELSE
                DY = YS
              END IF
            END IF
            IF (NQ .EQ. 1) THEN
              DY = SIGN(DY, SAVE2(J))
            ELSE
              DY = SIGN(DY, YH(J,3))
            END IF
            DY = (Y(J) + DY) - Y(J)
            YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            Y(J) = YJ
            FACTOR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
!C                                                                 Step 1
            DIFF = ABS(SAVE2(1) - SAVE1(1))
            IMAX = 1
            DO 150 I = 2,N
              IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                IMAX = I
                DIFF = ABS(SAVE2(I) - SAVE1(I))
              END IF
 150          CONTINUE
!C                                                                 Step 2
            IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT. 0.D0) THEN
              SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
!C                                                                 Step 3
              IF (DIFF .GT. BU*SCALE) THEN
                FAC(J) = MAX(FACMIN, FAC(J)*.1D0)
              ELSE IF (BR*SCALE .LE. DIFF .AND. DIFF .LE. BL*SCALE) THEN
                FAC(J) = MIN(FAC(J)*10.D0, FACMAX)
!C                                                                 Step 4
              ELSE IF (DIFF .LT. BR*SCALE) THEN
                FAC(J) = MIN(BP*FAC(J), FACMAX)
              END IF
            END IF
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
          NFE = NFE + N
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        END IF
        CALL DGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BU = UROUND**(.25D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 290 K = J,N,MW
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
 280          DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) THEN
                IF (FAC(K) .LT. FACMAX) THEN
                  FAC(K) = MIN(100.D0*FAC(K), FACMAX)
                  GO TO 280
                ELSE
                  DY = YS
                END IF
              END IF
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              DFDY(MW,K) = Y(K)
 290          Y(K) = Y(K) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            DO 330 K = J,N,MW
              Y(K) = DFDY(MW,K)
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
              DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) DY = YS
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              FACTOR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO 300 I = I1,I2
                I3 = K + I - MW
 300            DFDY(I,K) = FACTOR*(SAVE1(I3) - SAVE2(I3))
!C                                                                 Step 1
              IMAX = MAX(1, K - MU)
              DIFF = ABS(SAVE2(IMAX) - SAVE1(IMAX))
              I1 = IMAX
              I2 = MIN(K + ML, N)
              DO 310 I = I1+1,I2
                IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                  IMAX = I
                  DIFF = ABS(SAVE2(I) - SAVE1(I))
                END IF
 310            CONTINUE
!C                                                                 Step 2
              IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT.0.D0) THEN
                SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
!C                                                                 Step 3
                IF (DIFF .GT. BU*SCALE) THEN
                  FAC(K) = MAX(FACMIN, FAC(K)*.1D0)
                ELSE IF (BR*SCALE .LE.DIFF .AND. DIFF .LE.BL*SCALE) THEN
                  FAC(K) = MIN(FAC(K)*10.D0, FACMAX)
!C                                                                 Step 4
                ELSE IF (DIFF .LT. BR*SCALE) THEN
                  FAC(K) = MIN(BP*FAC(K), FACMAX)
                END IF
              END IF
 330          CONTINUE
 340        CONTINUE
          NFE = NFE + J2
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.D0
          DO 345 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 345 I = I1,I2
 345          DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))
          BND = 0.D0
          IF (DFDYMX .NE. 0.D0) THEN
            DO 350 J = 1,N
              I1 = MAX(ML+1, MW+1-J)
              I2 = MIN(MW+N-J, MW+ML)
              DO 350 I = I1,I2
 350            BND = BND + (DFDY(I,J)/DFDYMX)**2
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 380 I = I1,I2
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        END IF
        CALL DGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
      END IF
    END SUBROUTINE
    
      SUBROUTINE DDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,  &
       NDE,NQ,T,USERS,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D,JSTATE)
!C***BEGIN PROLOGUE  DDCOR
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***ROUTINES CALLED  DGESL,DGBSL,DNRM2
!C***END PROLOGUE  DDCOR
      DOUBLE PRECISION A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,  &
           SAVE1(*), SAVE2(*), DNRM2, T, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL EVALFA
!C***FIRST EXECUTABLE STATEMENT  DDCOR
      IF (MITER .EQ. 0) THEN
        DO 100 I = 1,N
 100      SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
        D = DNRM2(N, SAVE1, 1)/SQRT(DBLE(N))
        DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 150 I = 1,N
 150        SAVE2(I) = H*SAVE2(I)
          DO 160 J = 1,N
            DO 160 I = 1,N
 160          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL DGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        DO 200 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 200      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 250 I = 1,N
 250        SAVE2(I) = H*SAVE2(I)
          MW = ML + 1 + MU
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
              I3 = I + J - MW
 260          SAVE2(I3) = SAVE2(I3) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL DGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        DO 300 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 300      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
        DO 320 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 320      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      END IF
    END SUBROUTINE
    
      SUBROUTINE DDCST (MAXORD,MINT,ISWFLG,EL,TQ)
!C***BEGIN PROLOGUE  DDCST
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***ROUTINES CALLED  (NONE)
!C***END PROLOGUE  DDCST
      DOUBLE PRECISION EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
!C***FIRST EXECUTABLE STATEMENT  DDCST
      FACTRL(1) = 1.D0
      DO 10 I = 2,MAXORD
 10     FACTRL(I) = DBLE(I)*FACTRL(I-1)
!C                                             COMPUTE ADAMS COEFFICIENTS
      IF (MINT .EQ. 1) THEN
        GAMMA(1) = 1.D0
        DO 40 I = 1,MAXORD+1
          SUM = 0.D0
          DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/DBLE(I-J+2)
 40       GAMMA(I+1) = SUM
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        EL(2,2) = 1.D0
        EL(3,2) = 1.D0
        DO 60 J = 3,MAXORD
          EL(2,J) = FACTRL(J-1)
          DO 50 I = 3,J
 50         EL(I,J) = DBLE(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60       EL(J+1,J) = 1.D0
        DO 80 J = 2,MAXORD
          EL(1,J) = EL(1,J-1) + GAMMA(J)
          EL(2,J) = 1.D0
          DO 80 I = 3,J+1
 80         EL(I,J) = EL(I,J)/(DBLE(I-1)*FACTRL(J-1))
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.D0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.D0/GAMMA(J+1)
 100      TQ(3,J) = -1.D0/GAMMA(J+2)
!C                                              COMPUTE GEAR COEFFICIENTS
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        DO 130 J = 2,MAXORD
          EL(1,J) = FACTRL(J)
          DO 120 I = 2,J
 120        EL(I,J) = DBLE(J)*EL(I,J-1) + EL(I-1,J-1)
 130      EL(J+1,J) = 1.D0
        SUM = 1.D0
        DO 150 J = 2,MAXORD
          SUM = SUM + 1.D0/DBLE(J)
          DO 150 I = 1,J+1
 150        EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.D0/FACTRL(J-1)
          TQ(2,J) = DBLE(J+1)/EL(1,J)
 170      TQ(3,J) = DBLE(J+2)/EL(1,J)
      END IF
!C                          Compute constants used in the stiffness test.
!C                          These are the ratio of TQ(2,NQ) for the Gear
!C                          methods to those for the Adams methods.
      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.D0
          DO 190 I = 1,MXRD
            SUM = 0.D0
            DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/DBLE(I-J+2)
 190        GAMMA(I+1) = SUM
        END IF
        SUM = 1.D0
        DO 200 I = 2,MXRD
          SUM = SUM + 1.D0/DBLE(I)
 200      EL(1+I,1) = -DBLE(I+1)*SUM*GAMMA(I+1)
      END IF
    END SUBROUTINE
    
      SUBROUTINE DDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH)
!C***BEGIN PROLOGUE  DDSCL
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***END PROLOGUE  DDSCL
      DOUBLE PRECISION H, HMAX, RC, RH, RMAX, R1, YH(N,*)
!C***FIRST EXECUTABLE STATEMENT  DDSCL
      IF (H .LT. 1.D0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.D0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
    END SUBROUTINE
    
    SUBROUTINE DDPSC (KSGN,N,NQ,YH)
      use omp_lib
!C***BEGIN PROLOGUE  DDPSC
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!!C***END PROLOGUE  DDPSC
      DOUBLE PRECISION YH(N,*)
!C***FIRST EXECUTABLE STATEMENT  DDPSC
      IF (KSGN .GT. 0) THEN
         !dir$ assume_aligned YH:64
        DO 10 J1 = 1,NQ
          DO 10 J2 = J1,NQ
             J = NQ - J2 + J1
             !dir$ vector aligned
             !dir$ ivdep
             !$omp simd reduction(+:YH)
            DO 10 I = 1,N
 10           YH(I,J) = YH(I,J) + YH(I,J+1)
      ELSE
               !dir$ assume_aligned YH:64
        DO 30 J1 = 1,NQ
          DO 30 J2 = J1,NQ
             J = NQ - J2 + J1
              !dir$ vector aligned
             !dir$ ivdep
             !$omp simd reduction(+:YH)
            DO 30 I = 1,N
 30           YH(I,J) = YH(I,J) - YH(I,J+1)
      END IF
    END SUBROUTINE
    
      subroutine f (n, t, y, yp)
      double precision alfa,t,y(*),yp(*)
      common /const/ alfa, impl, miter
      save
      data istart /0/
      if (istart.eq.0) then
        istart = 1
        n = 0
        return
      end if
      yp(1) = 1.d0 + alfa*(y(2) - y(1)) - y(1)*y(3)
      yp(2) = alfa*(y(1) - y(2)) - y(2)*y(3)
      if (impl.eq.0 .or. impl.eq.1) then
        yp(3) = 1.d0 - y(3)*(y(1) + y(2))
      else if (impl.eq.2) then
        yp(3) = y(1) + y(2) - y(3)
      end if
    end subroutine f
    
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
!C***BEGIN PROLOGUE  DGBFA
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C   From the book "Numerical Methods and Software"
!C      by D. Kahaner, C. Moler, S. Nash
!C         Prentice Hall 1988
!C***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      DOUBLE PRECISION ABD(LDA,1)
!C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
!C
!C***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
!C
!C     ZERO INITIAL FILL-IN COLUMNS
!C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
!C
!C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
!C
!C        ZERO NEXT FILL-IN COLUMN
!C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
!C
!C        FIND L = PIVOT INDEX
!C
         LM = MIN0(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
!C
!C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
!C
!C           INTERCHANGE IF NECESSARY
!C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
!C
!C           COMPUTE MULTIPLIERS
!C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
!C
!C           ROW ELIMINATION WITH COLUMN INDEXING
!C
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
    END SUBROUTINE
    
      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
!C***BEGIN PROLOGUE  DGBSL
!C   THIS PROLOGUE HAS BEEN OMITTED FOR REASONS OF SPACE
!C   FOR A COMPLETE COPY OF THIS SUBROUTINE CONTACT THE AUTHORS
!C    From the book "Numerical Methods and Software"
!C       by D. Kahaner, C. Moler, S. Nash
!C          Prentice Hall 1988
!C***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      DOUBLE PRECISION ABD(LDA,1),B(1)
!C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
!C***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!C
!C        JOB = 0 , SOLVE  A * X = B
!C        FIRST SOLVE L*Y = B
!C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
!C
!C        NOW SOLVE  U*X = Y
!C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!C
!C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!C        FIRST SOLVE  TRANS(U)*Y = B
!C
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
!C
!C        NOW SOLVE TRANS(L)*X = Y
!C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
    END SUBROUTINE
    

    SUBROUTINE DGEFS(A,LDA,N,V,ITASK,IND,WORK,IWORK,RCOND)
#if 0
C***BEGIN PROLOGUE  DGEFS
C***DATE WRITTEN   800326   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  DGEFS solves a GENERAL double precision
C            NXN system of linear equations.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C    Subroutine DGEFS solves a general NxN system of double
C    precision linear equations using LINPACK subroutines DGECO
C    and DGESL.  That is, if A is an NxN double precision matrix
C    and if X and B are double precision N-vectors, then DGEFS
C    solves the equation
C
C                          A*X=B.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices U and L using partial pivoting.  These
C    factors and the pivoting information are used to find the
C    solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option to only solve (ITASK.GT.1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, N and IWORK must not have been altered by the user follow-
C    ing factorization (ITASK=1).  IND will not be changed by DGEFS
C    in this case.
C
C  Argument Description ***
C
C    A      DOUBLE PRECISION(LDA,N)
C             on entry, the doubly subscripted array with dimension
C               (LDA,N) which contains the coefficient matrix.
C             on return, an upper triangular matrix U and the
C               multipliers necessary to construct a matrix L
C               so that A=L*U.
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  The first N elements of
C             the array A are the elements of the first column of
C             the matrix A.  N must be greater than or equal to 1.
C             (terminal error message IND=-2)
C    V      DOUBLE PRECISION(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             If ITASK=1, the matrix A is factored and then the
C               linear equation is solved.
C             If ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A and IWORK.
C             If ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  see error message corresponding to IND below.
C    WORK   DOUBLE PRECISION(N)
C             a singly subscripted array of dimension at least N.
C    IWORK  INTEGER(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than 1.
C    IND=-3  terminal   ITASK is less than 1.
C    IND=-4  terminal   The matrix A is computationally singular.
C                         A solution has not been computed.
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the matrix
C                         A may be poorly scaled.
C
C               Note-  The above terminal(*fatal*) error messages are
C                      designed to be handled by XERRWV in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERROR.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C***REFERENCES  SUBROUTINE DGEFS WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C                 THE LINPACK SUBROUTINES USED BY DGEFS ARE DESCRIBED IN
C                 DETAIL IN THE *LINPACK USERS GUIDE* PUBLISHED BY
C                 THE SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS
C                 (SIAM) DATED 1979.
C***ROUTINES CALLED  D1MACH,DGECO,DGESL,XERROR,XERRWV
C***END PROLOGUE  DGEFS
C
#endif
      INTEGER LDA,N,ITASK,IND,IWORK(N)
      DOUBLE PRECISION A(LDA,N),V(N),WORK(N),D1MACH
      DOUBLE PRECISION RCOND
!C***FIRST EXECUTABLE STATEMENT  DGEFS
      IF (LDA.LT.N)  GO TO 101
      IF (N.LE.0)  GO TO 102
      IF (ITASK.LT.1) GO TO 103
      IF (ITASK.GT.1) GO TO 20
!C
!C     FACTOR MATRIX A INTO LU
      CALL DGECO(A,LDA,N,IWORK,RCOND,WORK)
!C
!C     CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
      IF (RCOND.EQ.0.0D0)  GO TO 104
!C
!C     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
      IND=-IDINT(DLOG10(D1MACH(4)/RCOND))
!C
!C     CHECK FOR IND GREATER THAN ZERO
      IF (IND.GT.0)  GO TO 20
      IND=-10
      CALL XERROR( 'DGEFS ERROR (IND=-10) -- SOLUTION MAY HAVE NO SIGNIF
     1ICANCE',58,-10,0)
!C
!C     SOLVE AFTER FACTORING
   20 CALL DGESL(A,LDA,N,IWORK,V,0)
      RETURN
!C
!C     IF LDA.LT.N, IND=-1, TERMINAL XERRWV MESSAGE
  101 IND=-1
      CALL XERRWV( 'DGEFS ERROR (IND=-1) -- LDA=I1 IS LESS THAN N=I2',
     148,-1,1,2,LDA,N,0,0,0)
      RETURN
!C
!C     IF N.LT.1, IND=-2, TERMINAL XERRWV MESSAGE
  102 IND=-2
      CALL XERRWV( 'DGEFS ERROR (IND=-2) -- N=I1 IS LESS THAN 1',
     143,-2,1,1,N,0,0,0,0)
      RETURN
!C
!C     IF ITASK.LT.1, IND=-3, TERMINAL XERRWV MESSAGE
  103 IND=-3
      CALL XERRWV( 'DGEFS ERROR (IND=-3) -- ITASK=I1 IS LESS THAN 1',
     147,-3,1,1,ITASK,0,0,0,0)
      RETURN
!C
!C     IF SINGULAR MATRIX, IND=-4, TERMINAL XERRWV MESSAGE
  104 IND=-4
      CALL XERRWV( 'DGEFS ERROR (IND=-4) -- SINGULAR MATRIX A - NO SOLUT
     1ION',55,-4,1,0,0,0,0,0,0)
      RETURN
!C
      END SUBROUTINE

      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
!C***BEGIN PROLOGUE  DGESL
!C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
!C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
!C     From the book "Numerical Methods and Software"
!C          by  D. Kahaner, C. Moler, S. Nash
!C               Prentice Hall 1988
!C***ROUTINES CALLED  DAXPY,DDOT
!C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
!C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
!C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
!C
!C        JOB = 0 , SOLVE  A * X = B
!C        FIRST SOLVE  L*Y = B
!C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
!C
!C        NOW SOLVE  U*X = Y
!C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
!C
!C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!C        FIRST SOLVE  TRANS(U)*Y = B
!C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
!C
!C        NOW SOLVE TRANS(L)*X = Y
!C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END

      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
!C***BEGIN PROLOGUE  DGECO
!C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
!C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
!C     From the book "Numerical Methods and Software"
!C          by  D. Kahaner, C. Moler, S. Nash
!C               Prentice Hall 1988
!C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL
!C***END PROLOGUE  DGECO
      INTEGER LDA,N,IPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
!C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
!C
!C     COMPUTE 1-NORM OF A
!C
!C***FIRST EXECUTABLE STATEMENT  DGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
!C
!C     FACTOR
!C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
!C
!C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!C     OVERFLOW.
!C
!C     SOLVE TRANS(U)*W = E
!C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!C
!C     SOLVE TRANS(L)*Y = W
!C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
!C
      YNORM = 1.0D0
!C
!C     SOLVE L*V = Y
!C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!C
!C     SOLVE  U*Z = V
!C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
!C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
!C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END SUBROUTINE

      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
C***BEGIN PROLOGUE  DGEFA
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END


    SUBROUTINE EIGEN(N,A,WR,WI,Z,RW,IW)
          !dir$ optimize:3
#if 0
C...
C...  SUBROUTINE EIGEN COMPUTES THE TEMPORAL EIGENVALUES OF AN NTH ORDER
C...  ODE SYSTEM WHOSE COEFFICIENT MATRIX IS STORED IN ARRAY A, AND
C...  OPTIONALLY, THE ASSOCIATED EIGENVECTORS.  EIGEN IN TURN CALLS
C...  EISPACK ROUTINES BALANC, ELMHES, ELTRAN, HQR2, AND OPTIONALLY,
C...  BALBAK IF THE EIGENVECTORS ARE TO BE COMPUTED (EISPACK IS A
C...  STANDARD LIBRARY OF ROUTINES FOR EIGENVALUE ANALYSIS OF LINEAR
C...  ALGEBRAIC SYSTEMS).
C...
C...  SUBROUTINES EIGEN, IN ALL CASES, COMPUTES AND PRINTS THE EIGEN-
C...  VALUES.  OPTIONALLY, IT COMPUTES, BUT DOES NOT PRINT, THE EIGEN-
C...  VECTORS (WHICH CAN BE PRINTED FROM THE CALLING PROGRAM)
C...
C...  ARGUMENT LIST
C...
C...     N       ORDER OF THE SQUARE MATRIX FOR WHICH THE EIGENVALUES,
C...             AND OPTIONALLY, THE EIGENVECTORS ARE TO BE COMPUTED
C...             (INPUT)
C...
C...             N LT 0, NO EIGENVECTORS
C...
C...             N GT 0, EIGENVECTORS ARE COMPUTED
C...
C...     A       MATRIX FOR WHICH EIGENVALUES AND EIGENVECTORS ARE TO
C...             BE COMPUTED (INPUT)
C...
C...     WR      ARRAY CONTAINING THE REAL PARTS OF THE N EIGENVALUES
C...             (OUTPUT)
C...
C...     WI      ARRAY CONTAINING THE IMAGINARY PARTS OF THE N EIGEN-
C...             VALUES (INPUT)
C...
C...     Z       ARRAY CONTAINING THE N EIGENVECTORS (OUTPUT)
C...
C...     RW,IW   REAL AND INTEGER WORK ARRAYS (NOT USED FOR INPUT OR
C...             OUTPUT)
C...
C...  THE FOLLOWING COMMENTS WERE TAKEN FROM THE EISPACK DOCUMENTATION:
C...
C...  THE EIGENVALUES OF A REAL GENERAL MATRIX ARE EITHER REAL OR
C...  COMPLEX CONJUGATE PAIRS.  WR(I) AND WI(I) CONTAIN THE REAL
C...  AND IMAGINARY PARTS OF THE I-TH EIGENVALUE.  THE EIGENVALUES ARE
C...  UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS OF EIGENVALUES
C...  APPEAR CONSECUTIVELY WITH THE EIGENVALUE HAVING THE POSITIVE
C...  IMAGINARY PART FIRST.
C...
C...  IF THE J-TH EIGENVALUE IS REAL, THEN THE J-TH COLUMN OF Z
C...  CONTAINS ITS EIGENVECTOR.
C...
C...  IF THE J-TH EIGENVALUE IS COMPLEX WITH POSITIVE IMAGINARY PART,
C...  THEN THE J-TH AND J+1-TH COLUMNS OF Z CONTAIN THE REAL AND IMAG-
C...  INARY PARTS OF ITS EIGENVECTOR.
C...
C...  THE CONJUGATE OF THIS VECTOR IS THE EIGENVECTOR FOR THE CONJUGATE
C...  EIGENVALUE.
C...
C...  THE EIGENVECTORS ARE NOT NORMALIZED.
C...
C...  VARIABLE DIMENSION THE ARRAYS USED BY THE EISPACK ROUTINES
#endif
      REAL A(N,N),   WR(N),   WI(N), Z(N,N), RW(N)
      INTEGER  NM      LOW,     UPP,  ERROR, IW(N)
!C...
!C...  COMMON/IO/ CONTAINS THE INPUT/OUTPUT UNIT (DEVICE) NUMBERS
      COMMON/IO/       NI,        NO
!C...
!C...  SET VARIABLE IVEC DEPENDING ON WHETHER EIGENVECTORS ARE TO BE
!C...  COMPUTED
!C...
!C...     NO EIGENVECTORS
         IF(N.LT.0)THEN
            N=-N
            IVEC=0
!C...
!C...     EIGENVECTORS WILL BE COMPUTED
         ELSE
            IVEC=1
         END IF
!C...
!C...  ROW DIMENSION OF THE ARRAYS USED BY THE EISPACK ROUTINES
      NM=N
!C...
!C...  COMPUTE THE EIGENVALUES
      CALL BALANC(NM,N,A,LOW,UPP,RW)
      CALL ELMHES(NM,N,LOW,UPP,A,IW)
      CALL ELTRAN(NM,N,LOW,UPP,A,IW,Z)
      CALL HQR2  (NM,N,LOW,UPP,A,WR,WI,Z,ERROR)
!C...
!C...  PRINT THE EIGENVALUES
      WRITE(NO,900)ERROR,LOW,UPP
      WRITE(NO,901)
      DO 1 I=1,N
      IF(ABS(WI(I)).GT.1.0E-10)GO TO 2
!C...
!C...  THE EIGENVALUE IS REAL
      WRITE(NO,902)I,WR(I),WI(I)
      GO TO 1
!C...
!C...  THE EIGENVALUE IS COMPLEX, SO THE DAMPING COEFFICIENT AND NATURAL
!C...  FREQUENCY ARE COMPUTED FOR PRINTING
2     RATIO=WI(I)/WR(I)
      ZD=SQRT(1.0E+00/(1.0E+00+RATIO**2))
      WN=ABS(WR(I))/ZD
      WRITE(NO,903)I,WR(I),WI(I),ZD,WN
1     CONTINUE
!C...
!C...  NO EIGENVECTORS
      IF(IVEC.EQ.0)THEN
         RETURN
!C...
!C...  EIGENVECTORS ARE COMPUTED
      ELSE
         CALL BALBAK(NM,N,LOW,UPP,RW,N,Z)
      END IF
900   FORMAT(1H1,' ERROR = ',I3,' LOW = ',I3,' UPP = ',I3,//)
901   FORMAT(4X,'I',11X,'REAL',11X,'IMAG',14X,'Z',13X,'WN',/)
902   FORMAT(I5,2F15.3)
903   FORMAT(I5,4F15.3)
      RETURN
      END
!C
!C     ------------------------------------------------------------------
!C
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
          !dir$ optimize:3
          !dir$ attributes forceinline :: BALANC
          !dir$ attributes optimization_parameter: "target_arch=skylake-avx512" :: BALANC
     
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
!C     REAL ABS
      LOGICAL NOCONV
#if 0
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT-
C
C        A CONTAINS THE BALANCED MATRIX,
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
C
C                **********
#endif 
      RADIX = 2.

      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
!C     ********** IN-LINE PROCEDURE FOR ROW AND
!C                COLUMN EXCHANGE **********
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
!C
       !dir$ assume_aligned A:64
       !dir$ vector aligned
       !dir$ vector always
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
!C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
!C
   50 GO TO (80,130), IEXC
!C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!C                AND PUSH THEM DOWN **********
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
!C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
!C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0) GO TO 120
  110    CONTINUE
!C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
!C
      GO TO 140
!C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!C                AND PUSH THEM LEFT **********
  130 K = K + 1
!C
  140 DO 170 J = K, L
!C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0) GO TO 170
  150    CONTINUE
!C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
!C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
      DO 180 I = K, L
  180 SCALE(I) = 1.0
!C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
  190 NOCONV = .FALSE.
!C
      DO 270 I = K, L
         C = 0.0
         R = 0.0
!C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
  200    CONTINUE

         G = R / RADIX
         F = 1.0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
!C     ********** NOW BALANCE **********
  240    IF ((C + R) / F .GE. 0.95 * S) GO TO 270
         G = 1.0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.

         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G

         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F

  270 CONTINUE

      IF (NOCONV) GO TO 190

  280 LOW = K
      IGH = L
      RETURN
!     ********** LAST CARD OF BALANC **********
      END

      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)

      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
!C     REAL ABS
      INTEGER INT(IGH)
#if 0
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N,
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT-
C
C        A CONTAINS THE HESSENBERG MATRIX.  THE MULTIPLIERS
C          WHICH WERE USED IN THE REDUCTION ARE STORED IN THE
C          REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX,
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
#endif
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200

      DO 180 M = KP1, LA
         MM1 = M - 1
         X = 0.0
         I = M

         DO 100 J = M, IGH
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE

         INT(M) = I
         IF (I .EQ. M) GO TO 130
!C    ********** INTERCHANGE ROWS AND COLUMNS OF A **********
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE

         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
!C    ********** END INTERCHANGE **********
  130    IF (X .EQ. 0.0) GO TO 180
         MP1 = M + 1
!C
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
!C
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
!C
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
!C
  160    CONTINUE
!C
  180 CONTINUE
!C
  200 RETURN
!C    ********** LAST CARD OF ELMHES **********
      END

      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
!C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
!C     REAL SQRT,ABS,SIGN
!C     INTEGER MIN0
      LOGICAL NOTLAS
      COMPLEX Z3
!C     COMPLEX CMPLX
      REAL T3(2)
      EQUIVALENCE (Z3,T3(1))
#if 0
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N,
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        H HAS BEEN DESTROYED,
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N,
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
C     PROCEDURE CDIV BY COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
#endif
      MACHEP = 2.**(-47)

      IERR = 0
!C     ********** STORE ROOTS ISOLATED BY BALANC **********
      DO 50 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0
   50 CONTINUE
!C
      EN = IGH
      T = 0.0
!C     ********** SEARCH FOR NEXT EIGENVALUES **********
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
!C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!C                FOR L=EN STEP -1 UNTIL LOW DO -- **********
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         IF (ABS(H(L,L-1)) .LE. MACHEP * (ABS(H(L-1,L-1))
     X      + ABS(H(L,L)))) GO TO 100
   80 CONTINUE
!C     ********** FORM SHIFT **********
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
!C     ********** FORM EXCEPTIONAL SHIFT **********
      T = T + X
!C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
!C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75 * S
      Y = X
      W = -0.4375 * S * S
  130 ITS = ITS + 1
!C     ********** LOOK FOR TWO CONSECUTIVE SMALL
!C                SUB-DIAGONAL ELEMENTS.
!C                FOR M=EN-2 STEP -1 UNTIL L DO -- **********
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         IF (ABS(H(M,M-1)) * (ABS(Q) + ABS(R)) .LE. MACHEP * ABS(P)
     X    * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
!C
  150 MP2 = M + 2
!C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0
  160 CONTINUE
!C     ********** DOUBLE QR STEP INVOLVING ROWS L TO EN AND
!C                COLUMNS M TO EN **********
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
!C     ********** ROW MODIFICATION **********
         DO 210 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
!C
         J = MIN0(EN,K+3)
!C     ********** COLUMN MODIFICATION **********
         DO 230 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
!C     ********** ACCUMULATE TRANSFORMATIONS **********
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            IF (.NOT. NOTLAS) GO TO 240
            P = P + ZZ * Z(I,K+2)
            Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K) = Z(I,K) - P
  250    CONTINUE
!C
  260 CONTINUE
!C
      GO TO 70
!C     ********** ONE ROOT FOUND **********
  270 H(EN,EN) = X  +  T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0
      EN = NA
      GO TO 60
!C     ********** TWO ROOTS FOUND **********
  280 P = (Y - X) / 2.0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0) GO TO 320
!C     ********** REAL PAIR **********
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0
      WI(EN) = 0.0
      X = H(EN,NA)
      R = SQRT(X*X+ZZ*ZZ)
      P = X / R
      Q = ZZ / R
!C     ********** ROW MODIFICATION **********
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
!C     ********** COLUMN MODIFICATION **********
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
!C     ********** ACCUMULATE TRANSFORMATIONS **********
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
!C
      GO TO 330
!C     ********** COMPLEX PAIR **********
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
!C     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
!C                VECTORS OF UPPER TRIANGULAR FORM **********
  340 NORM = 0.0
      K = 1

      DO 360 I = 1, N

         DO 350 J = K, N
  350    NORM = NORM + ABS(H(I,J))

         K = I
  360 CONTINUE

      IF (NORM .EQ. 0.0) GO TO 1001
!C     ********** FOR EN=N STEP -1 UNTIL 1 DO -- **********
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
!C     ********** REAL VECTOR **********
  600    M = EN
         H(EN,EN) = 1.0
         IF (NA .EQ. 0) GO TO 800
!C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = H(I,EN)
            IF (M .GT. NA) GO TO 620
!C
            DO 610 J = M, NA
  610       R = R + H(I,J) * H(J,EN)
!C
  620       IF (WI(I) .GE. 0.0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0) GO TO 640
            T = W
            IF (W .EQ. 0.0) T = MACHEP * NORM
            H(I,EN) = -R / T
            GO TO 700
!C     ********** SOLVE REAL EQUATIONS **********
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 700
  650       H(I+1,EN) = (-S - Y * T) / ZZ
  700    CONTINUE
!C     ********** END REAL VECTOR **********
         GO TO 800
!C     ********** COMPLEX VECTOR **********
  710    M = NA
!C     ********** LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
!C                EIGENVECTOR MATRIX IS TRIANGULAR **********
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    Z3 = CMPLX(0.0,-H(NA,EN)) / CMPLX(H(NA,NA)-P,Q)
         H(NA,NA) = T3(1)
         H(NA,EN) = T3(2)
  730    H(EN,NA) = 0.0
         H(EN,EN) = 1.0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
!C
         DO 790 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0
            SA = H(I,EN)
!C
            DO 760 J = M, NA
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
!C
            IF (WI(I) .GE. 0.0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 790
  770       M = I
            IF (WI(I) .NE. 0.0) GO TO 780
            Z3 = CMPLX(-RA,-SA) / CMPLX(W,Q)
            H(I,NA) = T3(1)
            H(I,EN) = T3(2)
            GO TO 790
!C     ********** SOLVE COMPLEX EQUATIONS **********
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0 * Q
            IF (VR .EQ. 0.0 .AND. VI .EQ. 0.0) VR = MACHEP * NORM
     X       * (ABS(W) + ABS(Q) + ABS(X) + ABS(Y) + ABS(ZZ))
            Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA) / CMPLX(VR,VI)
            H(I,NA) = T3(1)
            H(I,EN) = T3(2)
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN)) / CMPLX(ZZ,Q)
            H(I+1,NA) = T3(1)
            H(I+1,EN) = T3(2)
  790    CONTINUE
!C     ********** END COMPLEX VECTOR **********
  800 CONTINUE
!C     ********** END BACK SUBSTITUTION.
!C                VECTORS OF ISOLATED ROOTS **********
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
!C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
!C
  840 CONTINUE
!C     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!C                VECTORS OF ORIGINAL FULL MATRIX.
!                FOR J=N STEP -1 UNTIL LOW DO -- **********
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)

         DO 880 I = LOW, IGH
            ZZ = 0.0

            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)

            Z(I,J) = ZZ
  880 CONTINUE

      GO TO 1001
!C     ********** SET ERROR -- NO CONVERGENCE TO AN
!C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = EN
 1001 RETURN
!C     ********** LAST CARD OF HQR2 **********
      END
!C
!C     ------------------------------------------------------------------
!C
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
!C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
#if 0
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE ACCUMULATES THE STABILIZED ELEMENTARY
C     SIMILARITY TRANSFORMATIONS USED IN THE REDUCTION OF A
C     REAL GENERAL MATRIX TO UPPER HESSENBERG FORM BY  ELMHES.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N,
C
C        A CONTAINS THE MULTIPLIERS WHICH WERE USED IN THE
C          REDUCTION BY  ELMHES  IN ITS LOWER TRIANGLE
C          BELOW THE SUBDIAGONAL,
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION BY  ELMHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  ELMHES.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z TO IDENTITY MATRIX **********
#endif
      DO 80 I = 1, N

         DO 60 J = 1, N
   60    Z(I,J) = 0.0

         Z(I,I) = 1.0
   80 CONTINUE

      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
!C     ********** FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1

         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)

         I = INT(MP)
         IF (I .EQ. MP) GO TO 140

         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0
  130    CONTINUE

         Z(I,MP) = 1.0
  140 CONTINUE

  200 RETURN
!C     ********** LAST CARD OF ELTRAN **********
      END

      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)

      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),Z(NM,M)
      REAL S
#if 0
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  BALANC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC,
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  BALANC,
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
#endif
      IF (IGH .EQ. LOW) GO TO 120

      DO 110 I = LOW, IGH
         S = SCALE(I)
!C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!C                IF THE FOREGOING STATEMENT IS REPLACED BY
!C                S=1.0/SCALE(I). **********
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
!C
  110 CONTINUE
!C     ********- FOR I=LOW-1 STEP -1 UNTIL 1,
!C               IGH+1 STEP 1 UNTIL N DO -- **********
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
!C
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
!C
  140 CONTINUE
!C
     

   END SUBROUTINE


  
   SUBROUTINE RKF45(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
#if 0
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
C     WRITTEN BY H.A.WATTS AND L.F.SHAMPINE
C                   SANDIA LABORATORIES
C                  ALBUQUERQUE,NEW MEXICO
C
C    RKF45 IS PRIMARILY DESIGNED TO SOLVE NON-STIFF AND MILDLY STIFF
C    DIFFERENTIAL EQUATIONS WHEN DERIVATIVE EVALUATIONS ARE INEXPENSIVE.
C    RKF45 SHOULD GENERALLY NOT BE USED WHEN THE USER IS DEMANDING
C    HIGH ACCURACY.
C
C ABSTRACT
C
C    SUBROUTINE  RKF45  INTEGRATES A SYSTEM OF NEQN FIRST ORDER
C    ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C             DY(I)/DT = F(T,Y(1),Y(2),...,Y(NEQN))
C              WHERE THE Y(I) ARE GIVEN AT T .
C    TYPICALLY THE SUBROUTINE IS USED TO INTEGRATE FROM T TO TOUT BUT IT
C    CAN BE USED AS A ONE-STEP INTEGRATOR TO ADVANCE THE SOLUTION A
C    SINGLE STEP IN THE DIRECTION OF TOUT.  ON RETURN THE PARAMETERS IN
C    THE CALL LIST ARE SET FOR CONTINUING THE INTEGRATION. THE USER HAS
C    ONLY TO CALL RKF45 AGAIN (AND PERHAPS DEFINE A NEW VALUE FOR TOUT).
C    ACTUALLY, RKF45 IS AN INTERFACING ROUTINE WHICH CALLS SUBROUTINE
C    RKFS FOR THE SOLUTION.  RKFS IN TURN CALLS SUBROUTINE  FEHL WHICH
C    COMPUTES AN APPROXIMATE SOLUTION OVER ONE STEP.
C
C    RKF45  USES THE RUNGE-KUTTA-FEHLBERG (4,5)  METHOD DESCRIBED
C    IN THE REFERENCE
C    E.FEHLBERG , LOW-ORDER CLASSICAL RUNGE-KUTTA FORMULAS WITH STEPSIZE
C                 CONTROL , NASA TR R-315
C
C    THE PERFORMANCE OF RKF45 IS ILLUSTRATED IN THE REFERENCE
C    L.F.SHAMPINE,H.A.WATTS,S.DAVENPORT, SOLVING NON-STIFF ORDINARY
C                 DIFFERENTIAL EQUATIONS-THE STATE OF THE ART ,
C                 SANDIA LABORATORIES REPORT SAND75-0182 ,
C                 TO APPEAR IN SIAM REVIEW.
C
C
C    THE PARAMETERS REPRESENT-
C      F -- SUBROUTINE F(T,Y,YP) TO EVALUATE DERIVATIVES YP(I)=DY(I)/DT
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED
C      Y(*) -- SOLUTION VECTOR AT T
C      T -- INDEPENDENT VARIABLE
C      TOUT -- OUTPUT POINT AT WHICH SOLUTION IS DESIRED
C      RELERR,ABSERR -- RELATIVE AND ABSOLUTE ERROR TOLERANCES FOR LOCAL
C            ERROR TEST. AT EACH STEP THE CODE REQUIRES THAT
C                 ABS(LOCAL ERROR) .LE. RELERR*ABS(Y) + ABSERR
C            FOR EACH COMPONENT OF THE LOCAL ERROR AND SOLUTION VECTORS
C      IFLAG -- INDICATOR FOR STATUS OF INTEGRATION
C      WORK(*) -- ARRAY TO HOLD INFORMATION INTERNAL TO RKF45 WHICH IS
C            NECESSARY FOR SUBSEQUENT CALLS. MUST BE DIMENSIONED
C            AT LEAST  3+6*NEQN
C      IWORK(*) -- INTEGER ARRAY USED TO HOLD INFORMATION INTERNAL TO
C            RKF45 WHICH IS NECESSARY FOR SUBSEQUENT CALLS. MUST BE
C            DIMENSIONED AT LEAST  5
C
C
C  FIRST CALL TO RKF45
C
C    THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR THE ARRAYS
C    IN THE CALL LIST  -      Y(NEQN) , WORK(3+6*NEQN) , IWORK(5)  ,
C    DECLARE F IN AN EXTERNAL STATEMENT, SUPPLY SUBROUTINE F(T,Y,YP) AND
C    INITIALIZE THE FOLLOWING PARAMETERS-
C
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED.  (NEQN .GE. 1)
C      Y(*) -- VECTOR OF INITIAL CONDITIONS
C      T -- STARTING POINT OF INTEGRATION , MUST BE A VARIABLE
C      TOUT -- OUTPUT POINT AT WHICH SOLUTION IS DESIRED.
C            T=TOUT IS ALLOWED ON THE FIRST CALL ONLY, IN WHICH CASE
C            RKF45 RETURNS WITH IFLAG=2 IF CONTINUATION IS POSSIBLE.
C      RELERR,ABSERR -- RELATIVE AND ABSOLUTE LOCAL ERROR TOLERANCES
C            WHICH MUST BE NON-NEGATIVE. RELERR MUST BE A VARIABLE WHILE
C            ABSERR MAY BE A CONSTANT. THE CODE SHOULD NORMALLY NOT BE
C            USED WITH RELATIVE ERROR CONTROL SMALLER THAN ABOUT 1.E-8 .
C            TO AVOID LIMITING PRECISION DIFFICULTIES THE CODE REQUIRES
C            RELERR TO BE LARGER THAN AN INTERNALLY COMPUTED RELATIVE
C            ERROR PARAMETER WHICH IS MACHINE DEPENDENT. IN PARTICULAR,
C            PURE ABSOLUTE ERROR IS NOT PERMITTED. IF A SMALLER THAN
C            ALLOWABLE VALUE OF RELERR IS ATTEMPTED, RKF45 INCREASES
C            RELERR APPROPRIATELY AND RETURNS CONTROL TO THE USER BEFORE
C            CONTINUING THE INTEGRATION.
C      IFLAG -- +1,-1  INDICATOR TO INITIALIZE THE CODE FOR EACH NEW
C            PROBLEM. NORMAL INPUT IS +1. THE USER SHOULD SET IFLAG=-1
C            ONLY WHEN ONE-STEP INTEGRATOR CONTROL IS ESSENTIAL. IN THIS
C            CASE, RKF45 ATTEMPTS TO ADVANCE THE SOLUTION A SINGLE STEP
C            IN THE DIRECTION OF TOUT EACH TIME IT IS CALLED. SINCE THIS
C            MODE OF OPERATION RESULTS IN EXTRA COMPUTING OVERHEAD, IT
C            SHOULD BE AVOIDED UNLESS NEEDED.
C
C
C  OUTPUT FROM RKF45
C
C      Y(*) -- SOLUTION AT T
C      T -- LAST POINT REACHED IN INTEGRATION.
C      IFLAG = 2 -- INTEGRATION REACHED TOUT. INDICATES SUCCESSFUL RETUR
C                   AND IS THE NORMAL MODE FOR CONTINUING INTEGRATION.
C            =-2 -- A SINGLE SUCCESSFUL STEP IN THE DIRECTION OF TOUT
C                   HAS BEEN TAKEN. NORMAL MODE FOR CONTINUING
C                   INTEGRATION ONE STEP AT A TIME.
C            = 3 -- INTEGRATION WAS NOT COMPLETED BECAUSE RELATIVE ERROR
C                   TOLERANCE WAS TOO SMALL. RELERR HAS BEEN INCREASED
C                   APPROPRIATELY FOR CONTINUING.
C            = 4 -- INTEGRATION WAS NOT COMPLETED BECAUSE MORE THAN
C                   3000 DERIVATIVE EVALUATIONS WERE NEEDED. THIS
C                   IS APPROXIMATELY 500 STEPS.
C            = 5 -- INTEGRATION WAS NOT COMPLETED BECAUSE SOLUTION
C                   VANISHED MAKING A PURE RELATIVE ERROR TEST
C                   IMPOSSIBLE. MUST USE NON-ZERO ABSERR TO CONTINUE.
C                   USING THE ONE-STEP INTEGRATION MODE FOR ONE STEP
C                   IS A GOOD WAY TO PROCEED.
C            = 6 -- INTEGRATION WAS NOT COMPLETED BECAUSE REQUESTED
C                   ACCURACY COULD NOT BE ACHIEVED USING SMALLEST
C                   ALLOWABLE STEPSIZE. USER MUST INCREASE THE ERROR
C                   TOLERANCE BEFORE CONTINUED INTEGRATION CAN BE
C                   ATTEMPTED.
C            = 7 -- IT IS LIKELY THAT RKF45 IS INEFFICIENT FOR SOLVING
C                   THIS PROBLEM. TOO MUCH OUTPUT IS RESTRICTING THE
C                   NATURAL STEPSIZE CHOICE. USE THE ONE-STEP INTEGRATOR
C                   MODE.
C            = 8 -- INVALID INPUT PARAMETERS
C                   THIS INDICATOR OCCURS IF ANY OF THE FOLLOWING IS
C                   SATISFIED -   NEQN .LE. 0
C                                 T=TOUT  AND  IFLAG .NE. +1 OR -1
C                                 RELERR OR ABSERR .LT. 0.
C                                 IFLAG .EQ. 0  OR  .LT. -2  OR  .GT. 8
C      WORK(*),IWORK(*) -- INFORMATION WHICH IS USUALLY OF NO INTEREST
C                   TO THE USER BUT NECESSARY FOR SUBSEQUENT CALLS.
C                   WORK(1),...,WORK(NEQN) CONTAIN THE FIRST DERIVATIVES
C                   OF THE SOLUTION VECTOR Y AT T. WORK(NEQN+1) CONTAINS
C                   THE STEPSIZE H TO BE ATTEMPTED ON THE NEXT STEP.
C                   IWORK(1) CONTAINS THE DERIVATIVE EVALUATION COUNTER.
C
C
C  SUBSEQUENT CALLS TO RKF45
C
C    SUBROUTINE RKF45 RETURNS WITH ALL INFORMATION NEEDED TO CONTINUE
C    THE INTEGRATION. IF THE INTEGRATION REACHED TOUT, THE USER NEED ONL
C    DEFINE A NEW TOUT AND CALL RKF45 AGAIN. IN THE ONE-STEP INTEGRATOR
C    MODE (IFLAG=-2) THE USER MUST KEEP IN MIND THAT EACH STEP TAKEN IS
C    IN THE DIRECTION OF THE CURRENT TOUT. UPON REACHING TOUT (INDICATED
C    BY CHANGING IFLAG TO 2),THE USER MUST THEN DEFINE A NEW TOUT AND
C    RESET IFLAG TO -2 TO CONTINUE IN THE ONE-STEP INTEGRATOR MODE.
C
C    IF THE INTEGRATION WAS NOT COMPLETED BUT THE USER STILL WANTS TO
C    CONTINUE (IFLAG=3,4 CASES), HE JUST CALLS RKF45 AGAIN. WITH IFLAG=3
C    THE RELERR PARAMETER HAS BEEN ADJUSTED APPROPRIATELY FOR CONTINUING
C    THE INTEGRATION. IN THE CASE OF IFLAG=4 THE FUNCTION COUNTER WILL
C    BE RESET TO 0 AND ANOTHER 3000 FUNCTION EVALUATIONS ARE ALLOWED.
C
C    HOWEVER,IN THE CASE IFLAG=5, THE USER MUST FIRST ALTER THE ERROR
C    CRITERION TO USE A POSITIVE VALUE OF ABSERR BEFORE INTEGRATION CAN
C    PROCEED. IF HE DOES NOT,EXECUTION IS TERMINATED.
C
C    ALSO,IN THE CASE IFLAG=6, IT IS NECESSARY FOR THE USER TO RESET
C    IFLAG TO 2 (OR -2 WHEN THE ONE-STEP INTEGRATION MODE IS BEING USED)
C    AS WELL AS INCREASING EITHER ABSERR,RELERR OR BOTH BEFORE THE
C    INTEGRATION CAN BE CONTINUED. IF THIS IS NOT DONE, EXECUTION WILL
C    BE TERMINATED. THE OCCURRENCE OF IFLAG=6 INDICATES A TROUBLE SPOT
C    (SOLUTION IS CHANGING RAPIDLY,SINGULARITY MAY BE PRESENT) AND IT
C    OFTEN IS INADVISABLE TO CONTINUE.
C
C    IF IFLAG=7 IS ENCOUNTERED, THE USER SHOULD USE THE ONE-STEP
C    INTEGRATION MODE WITH THE STEPSIZE DETERMINED BY THE CODE OR
C    CONSIDER SWITCHING TO THE ADAMS CODES DE/STEP,INTRP. IF THE USER
C    INSISTS UPON CONTINUING THE INTEGRATION WITH RKF45, HE MUST RESET
C    IFLAG TO 2 BEFORE CALLING RKF45 AGAIN. OTHERWISE,EXECUTION WILL BE
C    TERMINATED.
C
C    IF IFLAG=8 IS OBTAINED, INTEGRATION CAN NOT BE CONTINUED UNLESS
C    THE INVALID INPUT PARAMETERS ARE CORRECTED.
C
C    IT SHOULD BE NOTED THAT THE ARRAYS WORK,IWORK CONTAIN INFORMATION
C    REQUIRED FOR SUBSEQUENT INTEGRATION. ACCORDINGLY, WORK AND IWORK
C    SHOULD NOT BE ALTERED.
C
C
#endif
      INTEGER NEQN,IFLAG,IWORK(5)
      DOUBLE PRECISION Y(NEQN),T,TOUT,RELERR,ABSERR,WORK(1)
!C     IF COMPILER CHECKS SUBSCRIPTS, CHANGE WORK(1) TO WORK(3+6*NEQN)
!C
      EXTERNAL F
!C
      INTEGER K1,K2,K3,K4,K5,K6,K1M
!C
!C
!C     COMPUTE INDICES FOR THE SPLITTING OF THE WORK ARRAY
!C
      K1M=NEQN+1
      K1=K1M+1
      K2=K1+NEQN
      K3=K2+NEQN
      K4=K3+NEQN
      K5=K4+NEQN
      K6=K5+NEQN
!C
!C     THIS INTERFACING ROUTINE MERELY RELIEVES THE USER OF A LONG
!C     CALLING LIST VIA THE SPLITTING APART OF TWO WORKING STORAGE
!C     ARRAYS. IF THIS IS NOT COMPATIBLE WITH THE USERS COMPILER,
!C     HE MUST USE RKFS DIRECTLY.
!C
      CALL RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),WORK(K1M),   &
               WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),WORK(K6),   &
               WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),IWORK(4),IWORK(5))

      RETURN
      END SUBROUTINE

      SUBROUTINE RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,YP,H,F1,F2,F3, &
                     F4,F5,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,KFLAG)
#if 0
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
C
C     RKFS INTEGRATES A SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL
C     EQUATIONS AS DESCRIBED IN THE COMMENTS FOR RKF45 .
C     THE ARRAYS YP,F1,F2,F3,F4,AND F5 (OF DIMENSION AT LEAST NEQN) AND
C     THE VARIABLES H,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,AND KFLAG ARE USED
C     INTERNALLY BY THE CODE AND APPEAR IN THE CALL LIST TO ELIMINATE
C     LOCAL RETENTION OF VARIABLES BETWEEN CALLS. ACCORDINGLY, THEY
C     SHOULD NOT BE ALTERED. ITEMS OF POSSIBLE INTEREST ARE
C         YP - DERIVATIVE OF SOLUTION VECTOR AT T
C         H  - AN APPROPRIATE STEPSIZE TO BE USED FOR THE NEXT STEP
C         NFE- COUNTER ON THE NUMBER OF DERIVATIVE FUNCTION EVALUATIONS
C
C
#endif
      LOGICAL HFAILD,OUTPUT

      INTEGER  NEQN,IFLAG,NFE,KOP,INIT,JFLAG,KFLAG
      DOUBLE PRECISION  Y(NEQN),T,TOUT,RELERR,ABSERR,H,YP(NEQN), &
       F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),SAVRE,      &
       SAVAE

      EXTERNAL F

      DOUBLE PRECISION  A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,REMIN,RER,S, &
       SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK

      INTEGER  K,MAXNFE,MFLAG

      DOUBLE PRECISION  DABS,DMAX1,DMIN1,DSIGN
!C
!C     REMIN IS THE MINIMUM ACCEPTABLE VALUE OF RELERR.  ATTEMPTS
!C     TO OBTAIN HIGHER ACCURACY WITH THIS SUBROUTINE ARE USUALLY
!C     VERY EXPENSIVE AND OFTEN UNSUCCESSFUL.
!C
      DATA REMIN/1.D-12/
!C
!C
!C     THE EXPENSE IS CONTROLLED BY RESTRICTING THE NUMBER
!C     OF FUNCTION EVALUATIONS TO BE APPROXIMATELY MAXNFE.
!C     AS SET, THIS CORRESPONDS TO ABOUT 500 STEPS.
!C
      DATA MAXNFE/3000/
!C
!C
!C     CHECK INPUT PARAMETERS
!C
!C
      IF (NEQN .LT. 1) GO TO 10
      IF ((RELERR .LT. 0.0D0)  .OR.  (ABSERR .LT. 0.0D0)) GO TO 10
      MFLAG=IABS(IFLAG)
      IF ((MFLAG .EQ. 0) .OR. (MFLAG .GT. 8)) GO TO 10
      IF (MFLAG .NE. 1) GO TO 20
!C
!C     FIRST CALL, COMPUTE MACHINE EPSILON
!C
      EPS = 1.0D0
    5 EPS = EPS/2.0D0
      EPSP1 = EPS + 1.0D0
      IF (EPSP1 .GT. 1.0D0) GO TO 5
      U26 = 26.0D0*EPS
      GO TO 50
!C
!C     INVALID INPUT
   10 IFLAG=8
      RETURN
!C
!C     CHECK CONTINUATION POSSIBILITIES
!C
   20 IF ((T .EQ. TOUT) .AND. (KFLAG .NE. 3)) GO TO 10
      IF (MFLAG .NE. 2) GO TO 25
!C
!C     IFLAG = +2 OR -2
      IF ((KFLAG .EQ. 3) .OR. (INIT .EQ. 0)) GO TO 45
      IF (KFLAG .EQ. 4) GO TO 40
      IF ((KFLAG .EQ. 5)  .AND.  (ABSERR .EQ. 0.0D0)) GO TO 30
      IF ((KFLAG .EQ. 6)  .AND.  (RELERR .LE. SAVRE)  .AND. &
         (ABSERR .LE. SAVAE)) GO TO 30
      GO TO 50
!C
!C     IFLAG = 3,4,5,6,7 OR 8
   25 IF (IFLAG .EQ. 3) GO TO 45
      IF (IFLAG .EQ. 4) GO TO 40
      IF ((IFLAG .EQ. 5) .AND. (ABSERR .GT. 0.0D0)) GO TO 45
!C
!C     INTEGRATION CANNOT BE CONTINUED SINCE USER DID NOT RESPOND TO
!C     THE INSTRUCTIONS PERTAINING TO IFLAG=5,6,7 OR 8
   30 STOP
!C
!C     RESET FUNCTION EVALUATION COUNTER
   40 NFE=0
      IF (MFLAG .EQ. 2) GO TO 50
!C
!C     RESET FLAG VALUE FROM PREVIOUS CALL
   45 IFLAG=JFLAG
      IF (KFLAG .EQ. 3) MFLAG=IABS(IFLAG)
!C!
!C     SAVE INPUT IFLAG AND SET CONTINUATION FLAG VALUE FOR SUBSEQUENT
!C     INPUT CHECKING
   50 JFLAG=IFLAG
      KFLAG=0
!C
!C     SAVE RELERR AND ABSERR FOR CHECKING INPUT ON SUBSEQUENT CALLS
      SAVRE=RELERR
      SAVAE=ABSERR
!C
!C     RESTRICT RELATIVE ERROR TOLERANCE TO BE AT LEAST AS LARGE AS
!C     2*EPS+REMIN TO AVOID LIMITING PRECISION DIFFICULTIES ARISING
!C     FROM IMPOSSIBLE ACCURACY REQUESTS
!C
      RER=2.0D0*EPS+REMIN
      IF (RELERR .GE. RER) GO TO 55
!C
!C     RELATIVE ERROR TOLERANCE TOO SMALL
      RELERR=RER
      IFLAG=3
      KFLAG=3
      RETURN
!C
   55 DT=TOUT-T
!C
      IF (MFLAG .EQ. 1) GO TO 60
      IF (INIT .EQ. 0) GO TO 65
      GO TO 80
!C
!C     INITIALIZATION --
!C                       SET INITIALIZATION COMPLETION INDICATOR,INIT
!C                       SET INDICATOR FOR TOO MANY OUTPUT POINTS,KOP
!C                       EVALUATE INITIAL DERIVATIVES
!C                       SET COUNTER FOR FUNCTION EVALUATIONS,NFE
!C                       ESTIMATE STARTING STEPSIZE
!C
   60 INIT=0
      KOP=0
!C
      A=T
      CALL F(A,Y,YP)
      NFE=1
      IF (T .NE. TOUT) GO TO 65
      IFLAG=2
      RETURN
!C
C
   65 INIT=1
      H=DABS(DT)
      TOLN=0.
      DO 70 K=1,NEQN
        TOL=RELERR*DABS(Y(K))+ABSERR
        IF (TOL .LE. 0.) GO TO 70
        TOLN=TOL
        YPK=DABS(YP(K))
        IF (YPK*H**5 .GT. TOL) H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF (TOLN .LE. 0.0D0) H=0.0D0
      H=DMAX1(H,U26*DMAX1(DABS(T),DABS(DT)))
      JFLAG=ISIGN(2,IFLAG)
!C
!C
!C     SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
!C
   80 H=DSIGN(H,DT)
!C
!C     TEST TO SEE IF RKF45 IS BEING SEVERELY IMPACTED BY TOO MANY
!C     OUTPUT POINTS
!C
      IF (DABS(H) .GE. 2.0D0*DABS(DT)) KOP=KOP+1
      IF (KOP .NE. 100) GO TO 85
!C
!C     UNNECESSARY FREQUENCY OF OUTPUT
      KOP=0
      IFLAG=7
      RETURN
!C
   85 IF (DABS(DT) .GT. U26*DABS(T)) GO TO 95
!C
!C     IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
!C
      DO 90 K=1,NEQN
   90   Y(K)=Y(K)+DT*YP(K)
      A=TOUT
      CALL F(A,Y,YP)
      NFE=NFE+1
      GO TO 300
!C
!C
!C     INITIALIZE OUTPUT POINT INDICATOR
!C
   95 OUTPUT= .FALSE.
!C!
!C     TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
!C     SCALE THE ERROR TOLERANCES
!C
      SCALE=2.0D0/RELERR
      AE=SCALE*ABSERR
!C
!C
!C     STEP BY STEP INTEGRATION
!C
  100 HFAILD= .FALSE.
!C
!C     SET SMALLEST ALLOWABLE STEPSIZE
!C
      HMIN=U26*DABS(T)
!C
!C     ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
!C     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE AND
!C     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
C
!      DT=TOUT-T
      IF (DABS(DT) .GE. 2.0D0*DABS(H)) GO TO 200
      IF (DABS(DT) .GT. DABS(H)) GO TO 150
!C
!C     THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
!C     OUTPUT POINT
!C
      OUTPUT= .TRUE.
      H=DT
      GO TO 200
!C
  150 H=0.5D0*DT
#if 0
C
C
C
C     CORE INTEGRATOR FOR TAKING A SINGLE STEP
C
C     THE TOLERANCES HAVE BEEN SCALED TO AVOID PREMATURE UNDERFLOW IN
C     COMPUTING THE ERROR TOLERANCE FUNCTION ET.
C     TO AVOID PROBLEMS WITH ZERO CROSSINGS,RELATIVE ERROR IS MEASURED
C     USING THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE
C     BEGINNING AND END OF A STEP.
C     THE ERROR ESTIMATE FORMULA HAS BEEN GROUPED TO CONTROL LOSS OF
C     SIGNIFICANCE.
C     TO DISTINGUISH THE VARIOUS ARGUMENTS, H IS NOT PERMITTED
C     TO BECOME SMALLER THAN 26 UNITS OF ROUNDOFF IN T.
C     PRACTICAL LIMITS ON THE CHANGE IN THE STEPSIZE ARE ENFORCED TO
C     SMOOTH THE STEPSIZE SELECTION PROCESS AND TO AVOID EXCESSIVE
C     CHATTERING ON PROBLEMS HAVING DISCONTINUITIES.
C     TO PREVENT UNNECESSARY FAILURES, THE CODE USES 9/10 THE STEPSIZE
C     IT ESTIMATES WILL SUCCEED.
C     AFTER A STEP FAILURE, THE STEPSIZE IS NOT ALLOWED TO INCREASE FOR
C     THE NEXT ATTEMPTED STEP. THIS MAKES THE CODE MORE EFFICIENT ON
C     PROBLEMS HAVING DISCONTINUITIES AND MORE EFFECTIVE IN GENERAL
C     SINCE LOCAL EXTRAPOLATION IS BEING USED AND EXTRA CAUTION SEEMS
C     WARRANTED.
C
C
C     TEST NUMBER OF DERIVATIVE FUNCTION EVALUATIONS.
C     IF OKAY,TRY TO ADVANCE THE INTEGRATION FROM T TO T+H
C
#endif
  200 IF (NFE .LE. MAXNFE) GO TO 220
!C
!C     TOO MUCH WORK
      IFLAG=4
      KFLAG=4
      RETURN
!C
!C     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
!C
  220 CALL FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE=NFE+5
!C
!C     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR ESTIMATES
!C     AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE ERROR IS
!C     MEASURED WITH RESPECT TO THE AVERAGE OF THE MAGNITUDES OF THE
!C     SOLUTION AT THE BEGINNING AND END OF THE STEP.
!C
      EEOET=0.0D0
      DO 250 K=1,NEQN
        ET=DABS(Y(K))+DABS(F1(K))+AE
        IF (ET .GT. 0.0D0) GO TO 240
!C
!C       INAPPROPRIATE ERROR TOLERANCE
        IFLAG=5
        RETURN
!C
  240   EE=DABS((-2090.0D0*YP(K)+(21970.0D0*F3(K)-15048.0D0*F4(K)))+
     1                        (22528.0D0*F2(K)-27360.0D0*F5(K)))
  250   EEOET=DMAX1(EEOET,EE/ET)
!C
      ESTTOL=DABS(H)*EEOET*SCALE/752400.0D0
!C
      IF (ESTTOL .LE. 1.0D0) GO TO 260
!C
!C
!C     UNSUCCESSFUL STEP
!C                       REDUCE THE STEPSIZE , TRY AGAIN
!C                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
!C
      HFAILD= .TRUE.
      OUTPUT= .FALSE.
      S=0.1D0
      IF (ESTTOL .LT. 59049.0D0) S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF (DABS(H) .GT. HMIN) GO TO 200
!C
!C     REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG=6
      KFLAG=6
      RETURN
!C
!C
!C     SUCCESSFUL STEP
!C                        STORE SOLUTION AT T+H
!C                        AND EVALUATE DERIVATIVES THERE
!C
  260 T=T+H
      DO 270 K=1,NEQN
  270   Y(K)=F1(K)
      A=T
      CALL F(A,Y,YP)
      NFE=NFE+1
!C
!C
!C                       CHOOSE NEXT STEPSIZE
!C                       THE INCREASE IS LIMITED TO A FACTOR OF 5
!C                       IF STEP FAILURE HAS JUST OCCURRED, NEXT
!C                          STEPSIZE IS NOT ALLOWED TO INCREASE
!C
      S=5.0D0
      IF (ESTTOL .GT. 1.889568D-4) S=0.9D0/ESTTOL**0.2D0
      IF (HFAILD) S=DMIN1(S,1.0D0)
      H=DSIGN(DMAX1(S*DABS(H),HMIN),H)

      IF (OUTPUT) GO TO 300
      IF (IFLAG .GT. 0) GO TO 100
     ONE-STEP MODE
      IFLAG=-2
      RETURN
   
  300 T=TOUT
      IFLAG=2
      RETURN

      END SUBROUTINE

      SUBROUTINE FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
#if 0
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
C    FEHL INTEGRATES A SYSTEM OF NEQN FIRST ORDER
C    ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C             DY(I)/DT=F(T,Y(1),---,Y(NEQN))
C    WHERE THE INITIAL VALUES Y(I) AND THE INITIAL DERIVATIVES
C    YP(I) ARE SPECIFIED AT THE STARTING POINT T. FEHL ADVANCES
C    THE SOLUTION OVER THE FIXED STEP H AND RETURNS
C    THE FIFTH ORDER (SIXTH ORDER ACCURATE LOCALLY) SOLUTION
C    APPROXIMATION AT T+H IN ARRAY S(I).
C    F1,---,F5 ARE ARRAYS OF DIMENSION NEQN WHICH ARE NEEDED
C    FOR INTERNAL STORAGE.
C    THE FORMULAS HAVE BEEN GROUPED TO CONTROL LOSS OF SIGNIFICANCE.
C    FEHL SHOULD BE CALLED WITH AN H NOT SMALLER THAN 13 UNITS OF
C    ROUNDOFF IN T SO THAT THE VARIOUS INDEPENDENT ARGUMENTS CAN BE
C    DISTINGUISHED.
C
C
#endif
      INTEGER  NEQN
      DOUBLE PRECISION  Y(NEQN),T,H,YP(NEQN),F1(NEQN),F2(NEQN), &
       F3(NEQN),F4(NEQN),F5(NEQN),S(NEQN)

      DOUBLE PRECISION  CH
      INTEGER  K

      CH=H/4.0D0
      DO 221 K=1,NEQN
  221   F5(K)=Y(K)+CH*YP(K)
      CALL F(T+CH,F5,F1)

      CH=3.0D0*H/32.0D0
      DO 222 K=1,NEQN
  222   F5(K)=Y(K)+CH*(YP(K)+3.0D0*F1(K))
      CALL F(T+3.0D0*H/8.0D0,F5,F2)

      CH=H/2197.0D0
      DO 223 K=1,NEQN
  223   F5(K)=Y(K)+CH*(1932.0D0*YP(K)+(7296.0D0*F2(K)-7200.0D0*F1(K)))
      CALL F(T+12.0D0*H/13.0D0,F5,F3)

      CH=H/4104.0D0
      DO 224 K=1,NEQN
  224   F5(K)=Y(K)+CH*((8341.0D0*YP(K)-845.0D0*F3(K))+ &
                                 (29440.0D0*F2(K)-32832.0D0*F1(K)))
      CALL F(T+H,F5,F4)

      CH=H/20520.0D0
      DO 225 K=1,NEQN
  225   F1(K)=Y(K)+CH*((-6080.0D0*YP(K)+(9295.0D0*F3(K)-
     1         5643.0D0*F4(K)))+(41040.0D0*F1(K)-28352.0D0*F2(K)))
      CALL F(T+H/2.0D0,F1,F5)
!C
!C     COMPUTE APPROXIMATE SOLUTION AT T+H
!C
      CH=H/7618050.0D0
      DO 230 K=1,NEQN
  230   S(K)=Y(K)+CH*((902880.0D0*YP(K)+(3855735.0D0*F3(K)- &
             1371249.0D0*F4(K)))+(3953664.0D0*F2(K)+   &
             277020.0D0*F5(K)))

      RETURN
    END SUBROUTINE


     
    SUBROUTINE JMAP(N,A,Y,YOLD,F,FOLD)
#if 0
C...
C...  ******************************************************************
C...
C...  SUBROUTINE JMAP MAPS THE JACOBIAN MATRIX OF AN NTH-ORDER SYSTEM OF
C...  ALGEBRAIC EQUATIONS OR FIRST-ORDER ORDINARY DIFFERENTIAL EQUATIONS
C...
C...  DOUBLE PRECISION VERSION
C...
C...  ******************************************************************
C...
C...  COPYRIGHT - LEHIGH UNIVERSITY, 1980
C...
C...  AUTHORS
C...
C...     G. R. DISSINGER
C...           AND
C...     W. E. SCHIESSER
C...     MOUNTAINTOP, BLDG A(111), RM D307
C...     LEHIGH UNIVERSITY
C...     BETHLEHEM, PA  18015
C...
C...     215/758-4264
C...
C...  ARGUMENT LIST
C...
C...     N       NUMBER OF ALGEBRAIC OR FIRST-ORDER ORDINARY DIFFEREN-
C...             TIAL EQUATIONS (ODES) FOR WHICH THE JACOBIAN MATRIX IS
C...             TO BE MAPPED, I.E., THE ORDER OF THE ALGEBRAIC OR ODE
C...             SYSTEM (INPUT)
C...
C...     A       TWO-DIMENSIONAL ARRAY CONTAINING THE JACOBIAN MATRIX
C...             OF THE NTH-ORDER ALGEBRAIC OR ODE SYSTEM (OUTPUT)
C...
C...     Y       ONE-DIMENSIONAL ARRAY OF N DEPENDENT VARIABLES SET TO
C...             INITIAL VALUES BY A CALL TO SUBROUTINE INITAL IN JMAP
C...             (INPUT TO SUBROUTINE JMAP VIA SUBROUTINE INITAL)
C...
C...     YOLD    ONE-DIMENSIONAL WORK ARRAY TO STORE THE N DEPENDENT
C...             VARIABLES IN Y TEMPORARILY SO THAT WHEN THE INDIVIDUAL
C...             Y*S ARE PERTURBED IN COMPUTING THE JACOBIAN MATRIX ELE-
C...             MENTS, THE Y*S CAN THEN BE RESTORED TO THEIR ORIGINAL
C...             VALUES
C...
C...     F       ONE-DIMENSIONAL ARRAY OF N DERIVATIVES OF THE DEPENDENT
C...             VARIABLES SET BY A CALL TO SUBROUTINE DERV IN JMAP
C...             (INPUT TO SUBROUTINE JMAP VIA SUBROUTINE DERV)
C...
C...     FOLD    ONE-DIMENSIONAL WORK ARRAY TO STORE THE N DERIVATIVES
C...             IN F SO THAT A FINITE DIFFERENCE APPROXIMATION OF THE
C...             INDIVIDUAL ELEMENTS OF THE JACOBIAN MATRIX CAN BE COM-
C...             PUTED
C...
C...  FROM THIS POINT ON, THE DISCUSSION IS IN TERMS OF FIRST-ORDER
C...  ODES, BUT IT CAN AS WELL BE PRESENTED IN TERMS OF LINEAR OR NON-
C...  ALGEBRAIC EQUATIONS, OR TRANSCENDENTAL EQUATIONS, OR A COMBINA-
C...  TION OF ALL THREE TYPES OF EQUATIONS.  THE ONLY REQUIREMENT IS
C...  THAT INITIAL CONDITIONS ABOUT WHICH THE JACOBIAN MATRIX IS TO BE
C...  COMPUTED MUST BE DEFINED IN SUBROUTINE INITAL, AND THE EQUATIONS
C...  THEMSELVES MUST BE PROGRAMMED IN SUBROUTINE DERV.
C...
C...  SUBROUTINES INITAL AND DERV DEFINING THE ODE SYSTEM ARE SUPPLIED
C...  BY THE USER, AS WELL AS A CALLING PROGRAM FOR JMAP.  THE MATHE-
C...  MATICAL CONCEPTS ARE EXPLAINED FURTHER IN THE FOLLOWING COMMENTS.
C...  THE CALCULATION OF THE INDIVIDUAL ELEMENTS OF THE JACOBIAN MATRIX
C...  BY A FINITE DIFFERENCE APPROXIMATION IS PROGRAMMED AS STATEMENT 14
C...  WITHIN DO LOOP 7.
C...
C...  THE SYSTEM OF ODES FOR WHICH THE JACOBIAN MATRIX IS MAPPED IS
C...
C...     DY1/DT = F1(Y1,Y2,...,YN)
C...
C...     DY2/DT = F2(Y1,Y2,...,YN)                                   (1)
C...       .              .
C...       .              .
C...       .              .
C...     DYN/DT = FN(Y1,Y2,...,YN)
C...
C...  WHICH CAN BE SUMMARIZED IN VECTOR FORM AS
C...
C...      -      - -
C...     DY/DT = F(Y)                                                (2)
C...
C...  WHERE
C...
C...     -                T
C...     Y = (Y1,Y2,..,YN)
C...
C...     -                T
C...     F = (F1,F2,..,FN)
C...
C...                              -
C...  SINCE THE DERIVATIVE VECTOR F IS IN GENERAL A NONLINEAR FUNCTION
C...                                   -
C...  OF THE DEPENDENT VARIABLE VECTOR Y, A TAYLOR SERIES EXPANSION
C...  TRUNCATED AFTER LINEAR TERMS GIVES A LINEARZED APPROXIMATION OF
C...  THE ORIGINAL SYSTEM
C...
C...      -      -
C...     DY/DT = JY                                                  (3)
C...
C...        -
C...  WHERE J IS THE JACOBIAN MATRIX OF THE ORIGINAL SYSTEM, I.E.,
C...
C...         ...         ...
C...         .F         F  .
C...         . 11        1N.
C...         .             .
C...     -   .    F        .
C...     J = .     IJ      .                                         (4)
C...         .             .
C...         .F         F  .
C...         . N1        NN.
C...         ...         ...
C...
C...  F   IS THE PARTIAL DERIVATIVE OF F  WITH RESPECT TO Y .  THUS THE
C...   IJ                               I                  J
C...  JACOBIAN MATRIX IS SQUARE (N X N).
C...                                                  -
C...  SUBROUTINE JMAP PRINTS A TWO-DIMENSIONAL MAP OF J WITH THE NUMBERS
C...  0 TO 9 INDICATING THE RELATIVE ORDER-OF-MAGNITUDE OF THE INDIVID-
C...  UAL ELEMENTS OF THE MATRIX.  THE VALUES OF THE ROW SUBSCRIPT I ARE
C...  PRINTED DOWN THE LEFT SIDE OF THE MAP AND THE VALUES OF THE COLUMN
C...  SUBSCRIPT J ARE PRINTED ACROSS THE TOP.  THE MAP IS PRINTED IN
C...  SECTIONS 99 COLUMNS WIDE.  THUS IF THE DIFFERENTIAL EQUATION
C...  SYSTEM IS GREATER THAN 99TH-ORDER, SUCCESSIVE SECTIONS OF THE MAP
C...  WILL BE PRINTED VERTICALLY.  THESE CAN THEN BE JOINED TOGETHER TO
C...  MAKE UP THE COMPLETE MAP.
C...
C...  THE N X N PARTIAL DERIVATIVES IN THE JACOBIAN MATRIX ARE COMPUTED
C...  APPROXIMATELY BY A SIMPLE DIFFERENCING PROCEDURE.  THE INITIAL
C...            -
C...  VALUES OF Y REQUIRED TO START THE CALCULATION ARE OBTAINED BY A
C...                                        -
C...  CALL TO SUBROUTINE DERV.  VALUES OF F ARE COMPUTED BY A SERIES
C...  OF CALLS TO SUBROUTINE DERV.  ALTHOUGH THESE SUBROUTINE NAMES
C...  PERTAIN SPECIFICALLY TO DSS/2, JMAP CAN EASILY BE ADAPTED FOR USE
C...  WITH ANY INITIAL-VALUE ODE INTEGRATION SYSTEM.
C...
C...  TYPE SELECTED REAL VARIABLES AS DOUBLE PRECISION (NOTE THIS DOES
C...  NOT INCLUDE ARRAY A FOR TWO REASONS (1) THE ODE JACOBIAN MATRIX
C...  STORED IN ARRAY A DOES NOT HAVE TO BE OF HIGH ACCURACY TO BE
C...  MAPPED BY SUBROUTINE JMAP, AND IN FACT, THE ELEMENTS OF THE ODE
C...  JACOBIAN MATRIX ARE NOT OF HIGH ACCURACY SINCE THEY ARE COMPUTED
C...  BY FIRST ORDER DIFFERENCES, AND (2) SINCE ARRAY IS OF SIZE N X N
C...  WHERE N IS THE NUMBER OF ODES, THIS SIZE INCREASES VERY RAPIDLY
C...  WITH N, AND THE STORAGE PROBLEM FOR ARRAY A IS MADE MORE SEVERE
C...  BY DECLARING ARRAY A AS DOUBLE PRECISION
#endif
      DOUBLE PRECISION      Y,   YOLD,      F,   FOLD,  FRACT,  YOLDJ
!C...
!C...  COMMON/IO/ CONTAINS THE INPUT/OUTPUT UNIT (DEVICE) NUMBERS
      COMMON/IO/       NI,        NO
!C...
!C...  VARIABLE DIMENSIONS FOR THE ARRAYS PASSED TO SUBROUTINE JMAP
      DIMENSION A(N,N),Y(N),YOLD(N),F(N),FOLD(N)
!C...
!C...  ABSOLUTE DIMENSIONS FOR THE ARRAYS USED IN SUBROUTINE JMAP
      DIMENSION SYM(11),SYMBOL(100),N1(100),N10(100),N100(100)
!C...
!C...  SET THE SCALE FACTOR, FRACTIONAL CHANGE USED IN THE CALCULATION
!C...  AND PROCESSING OF THE JACOBIAN MATRIX
      DATA SCALE,FRACT/1.0,0.001D+00/
!C...
!C...  DEFINE THE SYMBOLS USED IN PRINTING THE JACOBIAN MAP
      DATA SYM/1H ,1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
!C...
!C...  SET THE ARRAYS USED TO PRINT THE UNITS (N1), TENS (N10) AND
!C...  HUNDREDS (N100) PLACES IN THE COLUMN SUBSCRIPTS OF THE JACOBIAN
!C...  MAP
      DATA N1(1),N1(2),N1(3),N1(4),N1(5),N1(6),N1(7),N1(8),N1(9),N1(10)/ &
          0,1,2,3,4,5,6,7,8,9/
      DATA N10/10*0,10*1,10*2,10*3,10*4,10*5,10*6,10*7,10*8,10*9/
      DO 1 I=11,100
      N1(I)=N1(I-10)
1     CONTINUE
!C...                                     -     -
!C...  SET THE BASE VALUES OF THE VECTORS Y AND F
      CALL DERV

      DO 2 I=1,N
      YOLD(I)=Y(I)
2     CONTINUE

      CALL DERV

      DO 3 I=1,N
      Y(I)=YOLD(I)
      FOLD(I)=F(I)
3     CONTINUE

      NSECT=N/100+1

      DO 4 NSEC=1,NSECT

      JL=(NSEC-1)*100
      JU=JL+99

      IF(N.LT.JU)JU=N

      IF(NSEC.NE.1)GO TO 11

      JL=1
      JU=99
      IF(N.LT.JU)JU=N
      JUP1=JU+1

      WRITE(NO,908)

      IF(JUP1.LT.11)GO TO 5

      WRITE(NO,901)(N10(JCOL),JCOL=11,JUP1)
5     WRITE(NO,902)( N1(JCOL),JCOL= 2,JUP1)
      GO TO 8

11    JD=JU-JL+1
      DO 12 JCOL=1,JD

      N100(JCOL)=NSEC-1
12    CONTINUE

      IF(NSEC.EQ.1)WRITE(NO,908)

      WRITE(NO,904)(N100(JCOL),JCOL=1,JD)
      WRITE(NO,904)( N10(JCOL),JCOL=1,JD)
      WRITE(NO,904)(  N1(JCOL),JCOL=1,JD)

8     DO 6 JCOL=JL,JU
                                                                IJ
      YOLDJ= Y(JCOL)
                                                  J
      Y(JCOL)=Y(JCOL)*(1.D+00+FRACT)
      IF(DABS(Y(JCOL)).LT.FRACT)Y(JCOL)=Y(JCOL)+FRACT
                                               I
      CALL DERV

      DO 7 IROW=1,N

      IF(Y(JCOL).NE.YOLDJ)GO TO 14
      YOLDJ=YOLDJ+FRACT

14    A(IROW,JCOL)=SNGL((F(IROW)-FOLD(IROW))/(Y(JCOL)-YOLDJ))
7     CONTINUE

      DO 13 I=1,N
      Y(I)=YOLD(I)
13    CONTINUE
6     CONTINUE

      DO 9 IROW=1,N
      DO 10 JCOL=JL,JU
      JS=JCOL-JL+1
      IF(ABS(A(IROW,JCOL)).LE.0.00001*SCALE)SYMBOL(JS)=SYM(1)
      IF(ABS(A(IROW,JCOL)).GT.0.00001*SCALE)SYMBOL(JS)=SYM(2)
      IF(ABS(A(IROW,JCOL)).GT.0.00010*SCALE)SYMBOL(JS)=SYM(3)
      IF(ABS(A(IROW,JCOL)).GT.0.00100*SCALE)SYMBOL(JS)=SYM(4)
      IF(ABS(A(IROW,JCOL)).GT.0.01000*SCALE)SYMBOL(JS)=SYM(5)
      IF(ABS(A(IROW,JCOL)).GT.0.10000*SCALE)SYMBOL(JS)=SYM(6)
      IF(ABS(A(IROW,JCOL)).GT.1.00000*SCALE)SYMBOL(JS)=SYM(7)
      IF(ABS(A(IROW,JCOL)).GT.10.0000*SCALE)SYMBOL(JS)=SYM(8)
      IF(ABS(A(IROW,JCOL)).GT.100.000*SCALE)SYMBOL(JS)=SYM(9)
      IF(ABS(A(IROW,JCOL)).GT.1000.00*SCALE)SYMBOL(JS)=SYM(10)
      IF(ABS(A(IROW,JCOL)).GT.10000.0*SCALE)SYMBOL(JS)=SYM(11)
10    CONTINUE

      JD=JU-JL+1
      IF(NSEC.EQ.1)WRITE(NO,906)IROW,(SYMBOL(JS),JS=1,JD)
      IF(NSEC.GT.1)WRITE(NO,903)IROW,(SYMBOL(JS),JS=1,JD)
9     CONTINUE

4     CONTINUE
      RETURN
901   FORMAT(1H ,/,20X,90I1)
902   FORMAT(11X,99I1)
903   FORMAT(I8,2X,100A1)
904   FORMAT(10X,100I1)
906   FORMAT(I8,3X,100A1)
907   FORMAT(1H ,//, &
      58H DIVISION BY ZERO WAS ATTEMPTED IN DO LOOP 7 OF SUBROUTINE, &
      18H JMAP WITH IROW = ,I4,12H AND JCOL = ,I4,8H SO THAT,/,      &
      58H PROGRAM EXECUTION IS TERMINATED.  CHECK THE PROGRAMMING I, &
      58HN SUBROUTINES INITAL AND DERV FOR A DEPENDENT VARIABLE    , /, &
      58H WHICH IS CONSTANT, I.E., Y(JCOL) = YOLDJ.                )
908   FORMAT(1H ,/, &
      58H DEPENDENT VARIABLE COLUMN INDEX J (FOR YJ) IS PRINTED HOR,     &
      58HIZONTALLY                                                 ,//,  &
      58H DERIVATIVE ROW INDEX I (FOR DYI/DT = FI(Y1,Y2,...,YJ,...,,     &
      58HYN) IS PRINTED VERTICALLY                                 ,//,  &
      58H JACOBIAN MATRIX ELEMENT IN THE MAP WITH INDICES I,J IS FO,     &
      58HR PFI/PYJ WHERE P DENOTES A PARTIAL DERIVATIVE            ,//)
      END SUBROUTINE

   



  





end module mol_solvers
