c
c     file cud24sp.f
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
c     cud24sp attempts to improve the second order approximation generated
c     by cud2sp to a fourth order approximation using difference corrections
c     within multigrid cycling
c
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cud2sp.f, cudcom.f
c
c
c
      subroutine cud24sp(work,phi,ierror)
      implicit none
      complex phi(*),work(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      common/cud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      integer ierror,nx,ny
      ierror = 0
      nx = nfx
      ny = nfy
      if (min0(nx,ny).lt.6) then
	ierror = 30
	return
      end if
      call cu24sp(nx,ny,phi,work)
      end

      subroutine cu24sp(nx,ny,phi,wk)
      implicit none
      integer nx,ny
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      complex phi(nx,ny),wk(*)
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      common/cud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      integer k,kb,ncx,ncy,ipf,irf,ipc,irc,icjc,i,j,ic,jc,jnx,ij,icx,icy
      irf = krbgn(ngrid)
      icx = kcxbgn(ngrid)
      icy = kcybgn(ngrid)
c
c     pass truncation error down (weighted) to all grid levels
c
      call ctr2sp(nx,ny,phi,wk(irf),wk(icx),wk(icy))
c
c     set phif in wk(ipf) to zero
c
      ipf = kpbgn(ngrid)
      do j=0,ny+1
	do i=0,nx+1
	  ij = j*(nx+2)+i
	  wk(ipf+ij) = (0.0,0.0)
	end do
      end do
      do kb=2,ngrid
	k = ngrid-kb+1
	nx = nxk(k+1)
	ny = nyk(k+1)
	ncx = nxk(k)
	ncy = nyk(k)
	ipf = kpbgn(k+1)
	ipc = kpbgn(k)
	irf = krbgn(k+1)
	irc = krbgn(k)
c
c     set phic in wk(ipc) to zero
c
	do jc=0,ncy+1
	  do ic=0,ncx+1
	    icjc = jc*(ncx+2)+ic
	    wk(ipc+icjc) = (0.0,0.0)
	  end do
	end do
c
c      full weighting of truncation error from k+1 to k
c
	call cres2(nx,ny,wk(irf),ncx,ncy,wk(irc),nxa,nxb,nyc,nyd)
      end do
c
c     execute one fmw(2,1) cycle for correction
c
      kcycle = 2
      iprer = 2
      ipost = 1
      intpol = 3
c
c     lift correction approximation from coarsest grid
c
      do k=1,ngrid-1
	kcur = k
	call kcycd2sp(wk)
	nx = nxk(k+1)
	ny = nyk(k+1)
	ipf = kpbgn(k+1)
	ipc = kpbgn(k)
	ncx = nxk(k)
	ncy = nyk(k)
	call cprolon2(ncx,ncy,wk(ipc),nx,ny,wk(ipf),nxa,nxb,nyc,nyd,
     +                intpol)
      end do
      kcur = ngrid
c
c     execute one w(2,1) cycle from the finest grid level
c
      call kcycd2sp(wk)
c
c     add fourth order correction term to solution
c
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      do j=1,ny
	jnx = j*(nx+2)
	do i=1,nx
	  ij = jnx+i+1
	  phi(i,j) = phi(i,j)+wk(ij)
	end do
      end do
      return
      end

      subroutine ctr2sp(nx,ny,u,rhs,cofx,cofy)
c
c     estimate truncation error using second order approximation in u
c
      implicit none
      integer nx,ny
      complex rhs(nx,ny),u(nx,ny),cofx(nx,3),cofy(ny,3)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      common/cud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      real dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      integer i,j,istart,ifinal,jstart,jfinal
      complex cxx,cx,cyy,cy,tx,ty,ux3,ux4,uy3,uy4
      do j=1,ny
	do i=1,nx
	  rhs(i,j) = (0.0,0.0)
	end do
      end do
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlxx = dlx*dlx
      dlyy = dly*dly
      tdlx3 = 2.*dlxx*dlx
      tdly3 = 2.*dlyy*dly
      dlx4 = dlxx*dlxx
      dly4 = dlyy*dlyy
c
c     set subscript limit to avoid specified (dirchlet) boundaries.
c
      istart=1
      ifinal=nx
      jstart=1
      jfinal=ny
      if(nxa.eq.1) istart=2
      if(nxb.eq.1) ifinal=nx-1
      if(nyc.eq.1) jstart=2
      if(nyd.eq.1) jfinal=ny-1
c
c     compute truncation on deep interior
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,ux3,ux4,uy3,uy4,cxx,cx,cyy,cy,tx,ty) &
!$OMP SHARED (u,cofx,cofy,nx,ny) &
!$OMP SHARED (dlx,dly,dlxx,dlyy,dlx4,dly4,tdlx3,tdly3)
      do j=3,ny-2
	do i=3,nx-2
	  ux3 = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
	  ux4 = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))
     +          /dlx4
	  uy3 = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
	  uy4 = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))
     +          /dly4
c
c     reset pde coefficients from discretized coefs
c
	  cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	  cx = dlx*(cofx(i,2)-cofx(i,1))
	  cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	  cy = dly*(cofy(j,2)-cofy(j,1))
	  tx = (cxx*ux4*0.5+cx*ux3)/6.0
	  ty = (cyy*uy4*0.5+cy*uy3)/6.0
	  rhs(i,j) = dlxx*tx+dlyy*ty
	end do
      end do
c
c     estimate truncation error at and near nonspecified boundaries
c
      do i=istart,ifinal
	do j=jstart,2
	  cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	  cx = dlx*(cofx(i,2)-cofx(i,1))
	  cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	  cy = dly*(cofy(j,2)-cofy(j,1))
	  call cpde2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
	  tx = (cxx*ux4*0.5+cx*ux3)/6.0
	  ty = (cyy*uy4*0.5+cy*uy3)/6.0
c
c     adjust truncation error terms at mixed boundaries
c
	  if(nxa.ne.0 .and. (i.eq.1 .or. i.eq.nx)) then
	    if (i.eq.1) then
	      cxx = 0.5*dlxx*cofx(i,2)
	      tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	    else
	      cxx = 0.5*dlxx*cofx(i,1)
	      tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	    end if
	  end if
	  if (nyc.ne.0 .and. j.eq.1) then
	    cyy = 0.5*dlyy*cofy(j,2)
	    ty=cyy*(uy4*0.25+uy3/dly)/3.0
	  end if
	  rhs(i,j) = dlxx*tx+dlyy*ty
	end do
	do j=ny-1,jfinal
	  cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	  cx = dlx*(cofx(i,2)-cofx(i,1))
	  cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	  cy = dly*(cofy(j,2)-cofy(j,1))
	  call cpde2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
	  tx = (cxx*ux4*0.5+cx*ux3)/6.0
	  ty = (cyy*uy4*0.5+cy*uy3)/6.0
	  if(nxa.ne.0 .and. (i.eq.1 .or. i.eq.nx)) then
	    if (i.eq.1) then
	      cxx = 0.5*dlxx*cofx(i,2)
	      tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	    else
	      cxx = 0.5*dlxx*cofx(i,1)
	      tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	    end if
	  end if
	  if (nyc.ne.0 .and. j.eq.ny) then
	    cyy = 0.5*dlyy*cofy(j,1)
	    ty=cyy*(uy4*0.25-uy3/dly)/3.0
	  end if
	  rhs(i,j) = dlxx*tx+dlyy*ty
	end do
      end do
      do j=3,ny-2
	do i=istart,2
	  cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	  cx = dlx*(cofx(i,2)-cofx(i,1))
	  cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	  cy = dly*(cofy(j,2)-cofy(j,1))
	  call cpde2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
	  tx = (cxx*ux4*0.5+cx*ux3)/6.0
	  ty = (cyy*uy4*0.5+cy*uy3)/6.0
	  if(nxa.ne.0 .and. i.eq.1) then
	    cxx = 0.5*dlxx*cofx(i,2)
	    tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	  end if
	  rhs(i,j) = dlxx*tx+dlyy*ty
	end do
	do i=nx-1,ifinal
	  cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	  cx = dlx*(cofx(i,2)-cofx(i,1))
	  cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	  cy = dly*(cofy(j,2)-cofy(j,1))
	  call cpde2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
	  tx = (cxx*ux4*0.5+cx*ux3)/6.0
	  ty = (cyy*uy4*0.5+cy*uy3)/6.0
	  if(nxa.ne.0 .and. i.eq.nx) then
	  cxx = 0.5*dlxx*cofx(i,1)
	  tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	  end if
	  rhs(i,j) = dlxx*tx+dlyy*ty
	end do
      end do
      return
      end
