
!================================================================================!
! Various helper subroutines for Kalman Filtering and Least-Squares Methods.
!
! Adapted from the book: 
! Advanced Kalman Filtering, Least-Squares and Modeling: A Practical Handbook
! Callable from C-side code.
! Author(s): Bruce P. Gibbs
!
! Publisher: Wiley, Year: 2011
!
! ISBN: 0470529709,9780470529706,9780470890035,9780470890042
!** Author: B. Gibbs, 12/2009
!
!
! The author grants the user a non-exclusive, worldwide, royalty-free copyright license to

! 1. reproduce and modify the software for your own purposes, and
! 2. to distribute the software provided that you give credit to the author,
!    do not make claims against the author or Wiley-Interscience,
!    and mark any modifications as your own.  If the software is incorporated in
!    a commercial product, you  agree to defend and indemnify this author and
!    Wiley-Interscience against any losses, damages and costs arising from claims,
!    lawsuits and other legal actions brought by a third party.

! The software is provided on an as is basis, without warranties or conditions
! of any kind, either express or implied (including any warranties or conditions
! of title, non-infringement, merchantability or fitness for a particular purpose).
! The user is solely responsible for determining the appropriateness of using and
! distributing the program and assumes all risks associated with its exercise of
! rights under this agreement.
!********************************************************************************
! @Modified by Bernard Gingold on 08-05-2022 03:59 GMT+2, beniekg@gmail.com
!********************************************************************************


 subroutine bicgstab (ier,  xe,  a,b,n,nd)
      !dir$ attributes code_align : 32 :: bicgstab
      !dir$ optimize : 3
      !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: bicgstab
      !dir$ attributes noinline :: bicgstab
      !dir$ attributes concurrency_safe :: bicgstab
      use omp_lib
!** preconditioned bi-conjugate gradient stabilized method for solving A * xe =b.
!** where A is n x n.  Note version assumes that matrix A is full.  For most problems of this
!** type matrix A is very sparse, so sparse matrix storage should be used for A and the
!** internal multiplications modified accordingly.
!** Reference: C.T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, Philadelphia (1995)


      implicit none

      integer(4),intent(in) :: n         !row and column dimension of A
      integer(4),intent(in) :: nd        !row dimension of A used in the calling routine (nd >= n)
      real(8),intent(in) :: a(nd,n)      !square n x n matrix
      real(8),intent(in) :: b(n)         !n x 1 right-hand-side vector
      real(8),intent(inout) :: xe(n)     !output n x 1 vector
      integer(4),intent(out) :: ier      !error code: 0 =no error, 1=maximum iterations (200) reached

      integer(4) i,j
      integer(4) :: imax=200
      real(8) d(n),hs(n,n),r(n),r0(n),rho0,rhol,rho,alpha,beta,w
      real(8) p(n),v(n),s(n),t(n)
      !dir$ attribute align : 64 :: d
      !dir$ attribute align : 64 :: hs
      !dir$ attribute align : 64 :: r
      !dir$ attribute align : 64 :: r0
      !dir$ attribute align : 64 :: p
      !dir$ attribute align : 64 :: v
      !dir$ attribute align : 64 :: s
      !dir$ attribute align : 64 :: t
      real(8) :: eps = 1.d-6             !tolerance test for convergence
      !_____________________________________

      ier =1
      !** compute pre-conditioning as diagonal
      !dir$ assume_aligned d:64
      !dir$ assume_aligned a:64
      !dir$ assume_aligned hs:64
      !dir$ ivdep
      !dir$ code_align(32)
       !$omp simd simdlen(8) linear(i:1)
      do i=1,n
        d(i)=sqrt(sum(a(:,i)**2))
        hs(:,i) =a(:,i)/d(i)
      enddo
      xe(:) =xe(:)/d(:)

      r(:) =b(:) -matmul(hs,xe)
      r0(:) =r(:)
      rho0 =sqrt(sum(b(:)**2))
      rhol =1.0_8
      alpha =1.0_8
      w =1.0_8
      v(:) =0.0_8
      p(:) =0.0_8
      rho =dot_product(r0,r)
      !dir$ assume_aligned p:64
      !dir$ assume_aligned v:64
      !dir$ assume_aligned s:64
      !dir$ assume_aligned r:64
      !dir$ assume_aligned t:64
      !dir$ assume_aligned xe:64
      !dir$ ivdep
      do i=1,imax
        beta =(rho/rhol)*(alpha/w)
        p(:) =r(:) +beta*(p(:)-w*v(:))
        v(:) =matmul(hs,p)
        alpha =rho/dot_product(r0,v)
        s(:) =r(:) -alpha*v(:)
        t(:) =matmul(hs,s)
        w =dot_product(t,s)/sqrt(sum(t(:)**2))
        rho =w*dot_product(r0,t)
        xe(:) =xe(:) +alpha*p(:) +w*s(:)
        r(:) =s(:) -w*t(:)
        if (sqrt(sum(r(:)**2)) < eps*rho0) exit
        if (i >= imax) then
          ier =1
          return
        endif
      enddo

      ier =0
      xe(:) =xe(:)*d(:)

      end subroutine bicgstab


      subroutine cfactor(ier, a, n,eps)
        !dir$ attributes code_align : 32 :: cfactor
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cfactor
        !dir$ attributes noinline :: cfactor
        !dir$ attributes concurrency_safe :: cfactor
!** subroutine cfactor computes the Cholesky factor of a symmetric
!** positive definite matrix A, i.e., A =U^T * U where U is upper triangular.
!** Matrix A is stored as upper triangular by columns and the output U is stored in A.
!** If matrix A is singular at row i, the same row of U will be set to zero and
!** ier will be set to the first row that is found to be singular.  cfactor also tests
!** for a loss of precision or singularity and returns an error code indicating
!** at which row the loss occurred.

!** Author: B. Gibbs, 12/2009

!**************************************
! The author grants the user a non-exclusive, worldwide, royalty-free copyright license to

! 1. reproduce and modify the software for your own purposes, and
! 2. to distribute the software provided that you give credit to the author,
!    do not make claims against the author or Wiley-Interscience,
!    and mark any modifications as your own.  If the software is incorporated in
!    a commercial product, you  agree to defend and indemnify this author and
!    Wiley-Interscience against any losses, damages and costs arising from claims,
!    lawsuits and other legal actions brought by a third party.

! The software is provided on an as is basis, without warranties or conditions
! of any kind, either express or implied (including any warranties or conditions
! of title, non-infringement, merchantability or fitness for a particular purpose).
! The user is solely responsible for determining the appropriateness of using and
! distributing the program and assumes all risks associated with its exercise of
! rights under this agreement.
!**************************************


      implicit none

      integer(4),intent(in) :: n               !Dimension of A
      real(8),intent(in) :: eps                !Tolerance for loss of precision, e.g. eps =1d-4
      real(8),intent(inout) :: a((n*(n+1))/2)  !Symmetric matrix to be factored, or output factor
!      real(8),intent(out) :: errm               !minimum remaining precision
      integer(4),intent(out) :: ier            !Error flag: ier = -1 = n < 0, ier = -i means matrix
                                                ! is singular at row i, ier = +i, loss of precision
                                                ! exceeds eps at row i

      integer(4) i,k,kpiv,ki,iflg,km1
      real(8) tol,dsum,dpiv,work,wmax
      real(8) errm   !minimum remaining precision
      !____________________________________________________

      iflg =0
!      if (ier == -2) iflg=1

      !**  test for invalid n
      if (n < 1) then
        ier =-1
        return
      endif

      ier =0
      kpiv =0          !initialize diagonal index
      wmax =0.0_8

      do k=1,n         !row index
        kpiv =kpiv+k   !index of diagonal element of row
        ki =kpiv
        km1 =k-1
        !**  calculate tolerance for testing loss of significance
        tol =abs(eps*a(kpiv))

        !**   start factorization-loop over k-th row
        do i=k,n    !i is column index, ki is index of (k,i) element
          if (k > 1) then
            dsum =dot_product(a(kpiv-km1:kpiv-1),a(ki-km1:ki-1))
          else
            dsum =0.0_8
          endif

          !**  compute difference for a(ki)
          dsum =a(ki)-dsum

          if (i == k) then
            !**  diagonal element: test for negative pivot element and for loss of significance
            wmax =max(wmax,abs(a(kpiv)/dsum))   !(a(kpiv)/dsum)**1.3 matches actual errors - WHY ??
            if (dsum <= tol) then
              if (dsum <= 0.d0) then
                write (6,'(/"MATRIX IS SINGULAR AT ROW ",I3)') i
                if (ier >= 0) ier =-i    !set for first row that is singular
                dpiv =0.0_8              !set dpiv to zero elements in row k
                a(kpiv) =0.0_8
                !dsum =1.d40  !when called by sinv, set diagonal to big number to get "pseudo-inverse" ?
                !return
              else
                work =log10(abs(a(kpiv)/dsum))
                write (6,100) k,work
  100           format(/'AT ROW',i5,',',f7.1,                           &
     &            ' DIGITS WERE LOST IN MATRIX FACTORIZATION')
                if (ier == 0) ier =k-1
              endif
            endif

            if (dsum > 0.0_8) then
              !** compute pivot element
              a(kpiv) =sqrt(dsum)
              dpiv =1.0_8/a(kpiv)
            endif

          else      !**  calculate terms in row
            a(ki) =dsum*dpiv
          endif

          ki =ki+i
        enddo       !____ end i loop
      enddo      !____ end k loop

!      errm =1.1d-16*sqrt(wmax)      !little difference between using max vs RSS
      errm =1.1d-16*wmax     !1.1d-16 is mantissa LSB (53) of IEEE S_floating on PC

     
      end subroutine cfactor

      subroutine cgnr (ier,xe,h,y,n,m,mmax)
          !dir$ attributes code_align : 32 :: cgnr
          !dir$ optimize : 3
          !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cgnr
          !dir$ attributes noinline :: cgnr
          !dir$ attributes concurrency_safe :: cgnr
          use omp_lib
!** preconditioned conjugate gradient method for solving the least squares normal equations
!** to minimize the residual.  The measurement equation is y =H*xe +r.
!** Least squares normal equations are xe = (H^T*H)^(-1) * H^T*y.
!** References:
!**    1) C.T. Kelley, Iterative Methods for Linear and Nonlinear Equations,
!**       SIAM, Philadelphia (1995),
!**    2) Å Björck, Numerical Methods for Least Squares Problems, SIAM, Philadelphia (1996)

!** Author: B. Gibbs, 12/2009

!**************************************
! The author grants the user a non-exclusive, worldwide, royalty-free copyright license to

! 1. reproduce and modify the software for your own purposes, and
! 2. to distribute the software provided that you give credit to the author,
!    do not make claims against the author or Wiley-Interscience,
!    and mark any modifications as your own.  If the software is incorporated in
!    a commercial product, you  agree to defend and indemnify this author and
!    Wiley-Interscience against any losses, damages and costs arising from claims,
!    lawsuits and other legal actions brought by a third party.

! The software is provided on an as is basis, without warranties or conditions
! of any kind, either express or implied (including any warranties or conditions
! of title, non-infringement, merchantability or fitness for a particular purpose).
! The user is solely responsible for determining the appropriateness of using and
! distributing the program and assumes all risks associated with its exercise of
! rights under this agreement.
!**************************************
      implicit none

      integer(4),intent(in) :: n         !column dimension of matrix H
      integer(4),intent(in) :: m         !number of actual measurements (rows) in H and y
      integer(4),intent(in) :: mmax      !row dimension of H
      real(8),intent(in) :: h(mmax,n)    !measurement partial matrix
      real(8),intent(in) :: y(m)         !measurement vector
      real(8),intent(inout) :: xe(n)     !state vector
      integer(4),intent(out) :: ier      !error return: 0 =OK, 1 =not converged in 100 iterations

      integer(4) i,j
      integer(4) :: imax =100            !max allowed iterations
      real(8) hs(m,n),r(m),q(m),rho0,tau,taul,alpha
      !dir$ attributes align : 64 :: hs
      !dir$ attributes align : 64 :: r
      !dir$ attributes align : 64 :: q
      real(8) d(n),p(n),s(n)
      !dir$ attributes align : 64 :: d
      !dir$ attributes align : 64 :: p
      !dir$ attributes align : 64 :: s
      real(8) :: eps =1.d-12
      !_____________________________________

      ier =0
      !** compute pre-conditioning as diagonal
      !dir$ assume_aligned d:64
      !dir$ assume_aligned h:64
      !dir$ assume_aligned hs:64
      !dir$ code_align(32)
      !$omp simd simdlen(8) linear(i:1)
      do i=1,n
        d(i) =sqrt(sum(h(:m,i)**2))
!        d(i) =1.d0      !### test (little difference)
        hs(:,i) =h(:m,i)/d(i)
      enddo
      xe(:) =xe(:)*d(:)

      r(:) =y(:)-matmul(hs,xe)
      s(:) =matmul(r(:),hs)
      p(:) =s(:)
      rho0 =sqrt(sum(s(:)**2))
      taul =1.0_8
      !dir$ loop_count(100)
      do i=1,imax
        tau =sum(s(:)**2)
        if (i > 1) then
          p(:) =s(:) +(tau/taul)*p(:)
        endif
        taul =tau
        q(:) =matmul(hs,p)
        alpha=tau/sum(q(:)**2)
        xe(:)=xe(:)+alpha*p(:)
        s(:)=s(:)-alpha*matmul(q,hs)
!       write (6,'("cgnr: ",2i3,2g13.5)') n,i,sqrt(sum(s(:)**2)),rho0
        if (sqrt(sum(s(:)**2)) < eps*rho0) exit
        if (i >= imax) then
          ier =1
          xe(:) =xe(:)/d(:)
          return
        endif
      enddo

      xe(:) =xe(:)/d(:)

     
      end subroutine cgnr


