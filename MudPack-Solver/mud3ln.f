c
c     file mud3ln.f
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                  copyright (c) 2008 by UCAR                   *
c     *                                                               *
c     *       University Corporation for Atmospheric Research         *
c     *                                                               *
c     *                      all rights reserved                      *
c     *                                                               *
c     *                     MUDPACK  version 5.0.1                    *
c     *                                                               *
c     *                 A Fortran Package of Multigrid                *
c     *                                                               *
c     *                Subroutines and Example Programs               *
c     *                                                               *
c     *      for Solving Elliptic Partial Differential Equations      *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *                         John Adams                            *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the National Center for Atmospheric Research          *
c     *                                                               *
c     *                Boulder, Colorado  (80307)  U.S.A.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the National Science Foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c ... purpose
c
c     mud3ln.f contains subroutines for line relaxation in the x and y
c     and z direction.  This file must be loaded with any of the real
c     3-d mudpack solvers except mud3sp.
c
      subroutine slxmd3(nx,ny,nz,phi,cof,tx,sum,nxa,nyc,nze)
       !dir$ attributes code_align : 32 :: slxmd3
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: slxmd3
       use omp_lib
c
c     x line relaxation thru red and then black points in the
c     (y,z) plane for periodic or nonperiodic x b.c.
c
      implicit none
      integer nx,ny,nz,i,ib,j,k
      integer nxa,nyc,nze,nper
      real phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8),tx(nx,ny,nz,*)
      real sum(ny,nz)
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nxa.ne.0) then
c
c     x direction not periodic
c     first solve for x lines thru red points in (y,z) plane
c
      !dir$ assume_aligned phi:64
      !dir$ assume_aligned cof:64
       
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
           do j=1,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do i=2,nx
          
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
c
c     backward sweep
c
	  do j=1,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
             i = nx-ib+1
             
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
c
c     end of k odd loop
c
	end do
c
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned cof:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
           do j=2,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
c
c     end of k even loop
c
	end do

	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned cof:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
           do j=2,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
            !dir$ fma
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
	end do

c
c
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned cof:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
           do j=1,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do j=1,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     x direction periodic
c
         do k=1,nz
           
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8) 
	  do j=1,ny
	    sum(j,k) = 0.0
	  end do
	end do
c
c      sweep x lines thru red (y,z) forward and back
c
      !dir$ assume_aligned phi:64
      !dir$ assume_aligned cof:64  
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k) SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
           do j=1,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do i=2,nx-2
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
           
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp simd reduction(+:sum)  
	    do j=1,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
c
c     backward sweep
c
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))
     +                      /tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
       !dir$ assume_aligned phi:64
       !dir$ assume_aligned cof:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
           do j=2,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx-2
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
             
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp simd reduction(+:sum)  
	    do j=2,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                     tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned cof:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
           do j=2,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx-2
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
            
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp simd reduction(+:sum)
	    do j=2,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                     tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-odd
c
       !dir$ assume_aligned cof:64
       !dir$ assume_aligned sum:64
       !dir$ assume_aligned phi:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
           do j=1,ny,2
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx-2
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
            
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !$omp simd reduction(+:sum)
	    do j=1,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                     tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      end

      subroutine slymd3(nx,ny,nz,phi,cof,ty,sum,nxa,nyc,nze)
        !dir$ attributes code_align : 32 :: slymd3
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: slymd3
c
c     y line relaxation thru red and then black points in the
c     (x,z) plane for periodic or nonperiodic y b.c.
c
      implicit none
      integer nx,ny,nz,i,j,jb,k
      integer nxa,nyc,nze,nper
      real phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8),ty(ny,nx,nz,*)
      real sum(nx,nz)
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nyc.ne.0) then
c
c     y direction not periodic
c     first solve for y lines thru red points in (x,z) plane
c
      !dir$ assume_aligned phi:64
      !dir$ assume_aligned cof:64
       !dir$ assume_aligned ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k) SHARED(phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=1,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do j=2,ny
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !$omp simd reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
             j = ny-jb+1
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
       end do
       !dir$ assume_aligned phi:64
       !dir$ assume_aligned cof:64
       !dir$ assume_aligned ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=2,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp simd reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
c
        !dir$ assume_aligned cof:64
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=2,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny
            !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp simd reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=1,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny
               !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp simd reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
         end do
            !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp simd reduction(/:phi)
	  do i=1,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
             j = ny-jb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
             !dir$ fma
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     y direction periodic
c
         do k=1,nz
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
	  do i=1,nx
	    sum(i,k) = 0.0
	  end do
	end do
c
c      sweep y lines thru red (x,z) forward and back
c
        !dir$ assume_aligned cof:64,sum:64,phi:64,ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k)  SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=1,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do j=2,ny-2
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp simd reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp simd reduction(-:sum)
	    do i=1,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
             j = ny-jb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
         !dir$ assume_aligned cof:64,sum:64,phi:64,ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k) SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=2,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny-2
                !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp simd reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
                !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !dir$ simd reduction(+:sum)
	    do i=2,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp simd reduction(-:phi)
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
             j = ny-jb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
        !dir$ assume_aligned cof:64,sum:64,phi:64,ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k), SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=2,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny-2
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp simd reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
               !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !$omp simd reduction(+:sum)
	    do i=2,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
             j = ny-jb+1
               !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
        !dir$ assume_aligned sum:64,phi:64,cof:64,ty:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,jb,k), SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
           do i=1,nx,2
             do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny-2
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !$omp reduction(+:sum)
	    do i=1,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
             !$omp reduction(-:phi)
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
          do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
             j = ny-jb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      return
      end

      subroutine slzmd3(nx,ny,nz,phi,cof,tz,sum,nxa,nyc,nze)
       !dir$ attributes code_align : 32 :: slzmd3
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: slzmd3
c
c     z line relaxation thru red and then black points in the
c     (x,y) plane for periodic or nonperiodic z b.c.
c
      implicit none
      integer nx,ny,nz,i,j,k,kb
      integer nxa,nyc,nze,nper
      real phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8),tz(nz,nx,ny,*)
      real sum(nx,ny)
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nze.ne.0) then
c
c     z direction not periodic
c     first solve for z lines thru red points in (x,y) plane
c
         !dir$ assume_aligned phi:64,cof:64,tz:64
!$OMP PARALLEL DO SCHEDULE(STTAIC,8) PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=1,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     forward sweep
c
	  do k=2,nz
            !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
            !$omp reduction(+:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c     backward sweep
c
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(+:phi)
	  do i=1,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
         
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=2,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(/:phi)
	  do i=2,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
         !dir$ assume_aligned phi:64,cof:64,tz:64
c
c     solve for z lines thru black points in (x,y)plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=2,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
         end do
          !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(/:phi)
	  do i=2,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        !dir$ assume_aligned phi:64,cof:64,tz:64
c    
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,kb) SHARED(phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=1,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz
             !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi) 
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(/:phi)
	  do i=1,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     z direction periodic
c
         do j=1,ny
             !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
          
	  do i=1,nx
	    sum(i,j) = 0.0
	  end do
	end do
c
c      sweep z lines thru red (x,y) forward and back
c
        !dir$ assume_aligned sum:64,phi:64,cof:64,tz:64
!$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=1,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     forward sweep
c
	  do k=2,nz-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(+:sum)
	    do i=1,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
         end do
          !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(-:phi)
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     sweep black z lines thru (x,y)
c
        !dir$ assume_aligned sum:64,phi:64,cof:64,tz:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,kb) SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=2,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz-2
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(+:sum)
	    do i=2,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(-:phi)
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
         end do
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
         do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
             k = nz-kb+1
            !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        !dir$ assume_aligned sum:64,phi:64,cof:64,tz:64
!$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=2,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(+:sum)
	    do i=2,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
         end do
          !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(-:phi)
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
         end do
          !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
         do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
        

        !dir$ assume_aligned sum:64,phi:64,cof:64,tz:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,kb) SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=1,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz-2
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
              !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(+:sum)
	    do i=1,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
         end do
          !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !$omp reduction(-:phi)
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
         end do
          !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
             k = nz-kb+1
           !dir$ vector aligned
           !dir$ vector always
           !dir$ unroll(8)
           !dir$ fma
           !$omp reduction(-:phi)
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      return
      end
