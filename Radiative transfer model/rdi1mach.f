c ---------------------------------------------------------------------
c  Fortran-90 versions of machine-constant routines R1MACH, D1MACH, I1MACH
c
c  {R,D,I}1MACH revisited: no more uncommenting DATA statements
c  
c  Presented at the IFIP WG 2.5 International Workshop on 
c  "Current Directions in Numerical Software and High Performance 
c  Computing", 19 - 20 October 1995, Kyoto, Japan. 
c  
c  The widely-used original routines were modified to use Fortran-90 
c  intrinsic functions.  This was not completely possible with I1MACH, 
c  which returns some parameters (logical unit numbers of standard
c  input, standard output, and standard error) that may require
c  user customization. 
c  
c  David Gay (dmg@bell-labs.com)
c  Eric Grosse (ehg@bell-labs.com)
c  Bell Laboratories
c  700 Mountain Avenue
c  Murray Hill, New Jersey 07974-0636
c  USA 
c  
c  References:
c  
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996.
c
c http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson
c    /d1mach.html  (THIS WEB SITE WORKED AS OF APR 2000)
c -------------------------------------------------------------------------


      REAL(4) FUNCTION R1MACH (I)
c
c   R1MACH can be used to obtain machine-dependent parameters for
c   single precision numbers.  The results for various values of I are:
c
c   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   R1MACH(3) = B**(-T), the smallest relative spacing.
c   R1MACH(4) = B**(1-T), the largest relative spacing.
c   R1MACH(5) = LOG10(B)
c
c   Assume single precision numbers are represented in the T-digit,
c   base-B form
c
c              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c     Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c     August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      REAL(4) :: B, X = 1.0_4

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          R1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          R1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          R1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          R1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          R1MACH = LOG10(B)
        CASE DEFAULT
          R1MACH = 0.0_4
      END SELECT

     
      END FUNCTION R1MACH


     REAL(8) FUNCTION D1MACH (I)
c
c   D1MACH can be used to obtain machine-dependent parameters for
c   double precision numbers.  The results for various values of I are:
c
c   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   D1MACH(3) = B**(-T), the smallest relative spacing.
c   D1MACH(4) = B**(1-T), the largest relative spacing.
c   D1MACH(5) = LOG10(B)
c
c   Assume double precision numbers are represented in the T-digit,
c   base-B form
c
c        sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      REAL(8) :: B, X = 1.0_8

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          D1MACH = 0.0_8
      END SELECT

      
      END FUNCTION D1MACH


      INTEGER(4) FUNCTION I1MACH (I)
c
c   I1MACH can be used to obtain machine-dependent parameters for the
c   local machine environment.  The results for various values of I are:
c
c   I/O unit numbers (**MAY REQUIRE USER CUSTOMIZATION**):
c     I1MACH( 1) = the standard input unit.
c     I1MACH( 2) = the standard output unit.
c     I1MACH( 3) = the standard punch unit (obsolete, will cause error)
c     I1MACH( 4) = the standard error message unit.
c                  (the error message unit is usually 0 in UNIX systems)
c
c   Words:
c     I1MACH( 5) = the number of bits per integer storage unit.
c     I1MACH( 6) = the number of characters per integer storage unit.
c                  (obsolete, will cause an error)
c
c   Integers:
c     assume integers are represented in the S-digit, base-A form
c
c          sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
c
c     where 0 <= X(I) < A for I=0,...,S-1.
c
c     I1MACH( 7) = A, the base.
c     I1MACH( 8) = S, the number of base-A digits.
c     I1MACH( 9) = A**S - 1, the largest magnitude.
c
c   Floating-Point Numbers:
c     Assume floating-point numbers are represented in the T-digit,
c     base-B form
c                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c     where 0 <= X(I) .LT. B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c     I1MACH(10) = B, the base.
c
c   Single-Precision:
c     I1MACH(11) = T, the number of base-B digits.
c     I1MACH(12) = EMIN, the smallest exponent E.
c     I1MACH(13) = EMAX, the largest exponent E.
c
c   Double-Precision:
c     I1MACH(14) = T, the number of base-B digits.
c     I1MACH(15) = EMIN, the smallest exponent E.
c     I1MACH(16) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   750101  DATE WRITTEN
c   960411  Modified for Fortran 90 (BE after suggestions by Eric Grosse)    
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(4) :: I
      REAL(4)  :: X_single  = 1.0_4
      REAL(8)  :: X_double = 1.0_8

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          I1MACH = -1
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          I1MACH = -1
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X_single)
        CASE (11)
          I1MACH = DIGITS(X_single)
        CASE (12)
          I1MACH = MINEXPONENT(X_single)
        CASE (13)
          I1MACH = MAXEXPONENT(X_single)
        CASE (14)
          I1MACH = DIGITS(X_double)
        CASE (15)
          I1MACH = MINEXPONENT(X_double)
        CASE (16)
          I1MACH = MAXEXPONENT(X_double) 
        CASE DEFAULT
          I1MACH = -1
      END SELECT

     
      END FUNCTION I1MACH

