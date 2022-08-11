c
c     file cud2sp.f
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
c     cud2sp attempts to produce a second order finite difference
c     approximation to the complex two dimensional separable elliptic
c     partial differential equation of the form:
c
c       cxx(x)*pxx + cx(x)*px + cex(x)*p(x,y) +
c
c       cyy(y)*pyy + cy(y)*py + cey(y)*p(x,y) = r(x,y)
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cudcom.f
c
c
      subroutine cud2sp(iparm,fparm,work,cfx,cfy,bndyc,rhs,phi,
     +                  mgopt,ierror)
      implicit none
      integer iparm,mgopt,ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,itx,ity,icx,icy
      dimension iparm(17),fparm(6),mgopt(4)
      complex work(*),phi(*),rhs(*)
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      common/cud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      external cfx,cfy,bndyc
      data int / 0 /
      save int
      ierror = 1
      intl = iparm(1)    ! set and check intl on all calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
	int = 1
	if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
c
c     set  arguments internally
c     these will not be rechecked if intl=1!
c
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      ixp = iparm(6)
      jyq = iparm(7)
      iex = iparm(8)
      jey = iparm(9)
      ngrid = max0(iex,jey)
      nfx = iparm(10)
      nfy = iparm(11)
      iguess = iparm(12)
      maxcy = iparm(13)
      method = iparm(14)
      nwork = iparm(15)
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
c       set defaults
	kcycle = 2
	iprer = 2
	ipost = 1
	intpol = 3
      else
	iprer = mgopt(2)
	ipost = mgopt(3)
	intpol = mgopt(4)
      end if
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      tolmax = fparm(5)
      if (intl .eq. 0) then  ! intialization call
c
c     check input arguments
c
	ierror = 2   ! check boundary condition flags
	if (max0(nxa,nxb,nyc,nyd).gt.2) return
	if (min0(nxa,nxb,nyc,nyd).lt.0) return
	if (nxa.eq.0.and.nxb.ne.0) return
	if (nxa.ne.0.and.nxb.eq.0) return
	if (nyc.eq.0.and.nyd.ne.0) return
	if (nyc.ne.0.and.nyd.eq.0) return
	ierror = 3   ! check grid sizes
	if (ixp.lt.2) return
	if (jyq.lt.2) return
	ierror = 4
	ngrid = max0(iex,jey)
	if (iex.lt.1) return
	if (jey.lt.1) return
	if (ngrid.gt.50) return
	ierror = 5
	if (nfx.ne.ixp*2**(iex-1)+1) return
	if (nfy.ne.jyq*2**(jey-1)+1) return
	ierror = 6
	if (iguess*(iguess-1).ne.0) return
	ierror = 7
	if (maxcy.lt.1) return
	ierror = 8
	if (method.lt.0 .or. method.gt.3) return
	ierror = 9
c       compute and test minimum work space
	isx = 0
	if (method.eq.1 .or. method.eq.3) then
	  if (nxa.ne.0) isx = 3
	  if (nxa.eq.0) isx = 5
	end if
	jsy = 0
	if (method.eq.2 .or. method.eq.3) then
	  if (nyc.ne.0) jsy = 3
	  if (nyc.eq.0) jsy = 5
	end if
	kps = 1
	do k=1,ngrid
c       set subgrid sizes
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kps = kps+(nx+2)*(ny+2)+(1+isx+jsy)*nx*ny+3*(nx+ny)
	end do
	iparm(16) = kps+(nfx+2)*(nfy+2)   ! exact minimum work space
	lwork = iparm(16)
	if (lwork .gt. nwork) return
	ierror = 10   ! check solution region
	if (xb.le.xa .or. yd.le.yc) return
	ierror = 11
	if (tolmax .lt. 0.0) return
	ierror = 12   ! multigrid parameters
	if (kcycle.lt.0) return
	if (min0(iprer,ipost).lt.1) return
	if ((intpol-1)*(intpol-3).ne.0) return
	if (max0(kcycle,iprer,ipost).gt.2) then
	  ierror = -5   ! inefficient multigrid cycling
	end if
	if (ierror .gt. 0) ierror = 0   ! no fatal errors
c
c     set work space pointers and discretize pde at each grid level
c
	iw = 1
	do kb=1,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kpbgn(k) = iw
	  krbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
	  kcxbgn(k) = krbgn(k)+nx*ny
	  kcybgn(k) = kcxbgn(k)+3*nx
	  ktxbgn(k) = kcybgn(k)+3*ny
	  ktybgn(k) = ktxbgn(k)+isx*nx*ny
	  iw = ktybgn(k)+jsy*nx*ny
	  icx = kcxbgn(k)
	  icy = kcybgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call discd2sp(nx,ny,work(icx),work(icy),work(itx),work(ity),
     +                  bndyc,cfx,cfy,work,ierror)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call cud2sp1(nx,ny,rhs,phi,cfx,cfy,bndyc,work)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end

      subroutine cud2sp1(nx,ny,rhsf,phif,cfx,cfy,bndyc,wk)
      implicit none
      integer nx,ny
      complex phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ipc,ir,irc
      integer ncx,ncy,jj,ij,i,j,iter
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      common/cud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      external cfx,cfy,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ir = krbgn(ngrid)
c
c     set phif,rhsf in wk and adjust right hand side
c
      call cswk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
c
c     no initial guess at finest grid level!
c
	do kb=2,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ir = krbgn(k+1)
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ipc = kpbgn(k)
	  irc = krbgn(k)
c
c     transfer down to all grid levels
c
	  call ctrsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,
     +                 wk(ipc),wk(irc))
	end do
c
c     adjust right hand side at all grid levels in case
c     rhs or specified b.c. in phi or gbdy changed
c
	do k=1,ngrid
	  nx = nxk(k)
	  ny = nyk(k)
	  ip = kpbgn(k)
	  ir = krbgn(k)
	  call adjcd2sp(nx,ny,wk(ip),wk(ir),bndyc,cfx,cfy)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcycd2sp(wk)
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ipc = kpbgn(k)
	  ncx = nxk(k)
	  ncy = nyk(k)
c
c     lift or prolong approximation from k to k+1
c
	  call cprolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,
     +                  nyc,nyd,intpol)
	end do
      else
c
c     adjust rhs at finest grid level only
c
	nx = nxk(ngrid)
	ny = nyk(ngrid)
	ip = kpbgn(ngrid)
	ir = krbgn(ngrid)
	call adjcd2sp(nx,ny,wk(ip),wk(ir),bndyc,cfx,cfy)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcycd2sp(wk)
	if (tolmax.gt.0.0) then
c
c      error control
c
	  relmax = 0.0
	  phmax = 0.0
	  do j=1,nfy
	    jj = j*(nfx+2)
	    do i=1,nfx
	      ij = jj+i+1
	      phmax = amax1(phmax,cabs(wk(ij)))
	      relmax = amax1(relmax,cabs(wk(ij)-phif(i,j)))
	      phif(i,j) = wk(ij)
	    end do
	  end do
c
c     set maximum relative difference and check for convergence
c
	  if (phmax.gt.0.0) relmax = relmax/phmax
	  if (relmax.le.tolmax) return
	end if
      end do
c
c     set final interate after maxcy cycles in phif
c
      do j=1,nfy
	jj = j*(nfx+2)
	do i=1,nfx
	  ij = jj+i+1
	  phif(i,j) = wk(ij)
	end do
      end do
      return
      end

      subroutine kcycd2sp(wk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      complex wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ipc,ir,irc,icx,icy,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      common/cud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ir = krbgn(klevel)
      icx = kcxbgn(klevel)
      icy = kcybgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
c
c     prerelax at current finest grid level
c
      do l=1,iprer
	call relcd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = krbgn(klevel-1)
      call rescd2sp(nx,ny,wk(ip),wk(ir),ncx,ncy,wk(ipc),wk(irc),
     +              wk(icx),wk(icy),wk(kps))
c
c    set counter for grid levels to zero
c
      do l = 1,kcur
	kount(l) = 0
      end do
c
c    set new grid level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c   kcycle control point
c
   10 continue
c
c      post relax when kcur revisited
c
      if (klevel .eq. kcur) go to 5
c
c   count hit at current level
c
      kount(klevel) = kount(klevel)+1
c
c   relax at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ir = krbgn(klevel)
      icx = kcxbgn(klevel)
      icy = kcybgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
	call relcd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle complete at klevel
c
	ipc = ip
	ncx = nx
	ncy = ny
	ip = kpbgn(klevel+1)
	ir = krbgn(klevel+1)
	icx = kcxbgn(klevel+1)
	icy = kcybgn(klevel+1)
	nx = nxk(klevel+1)
	ny = nyk(klevel+1)
c
c    inject correction to finer grid
c
	call ccor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +             intpol,wk(kps))
c
c    reset counter to zero
c
	kount(klevel) = 0
c
c     ascend to next higher level and set to postrelax there
c
	klevel = klevel+1
	nrel = ipost
	go to 10
      else
	if (klevel .gt. 1) then
c
c    kcycle not complete so descend unless at coarsest grid
c
	  ipc = kpbgn(klevel-1)
	  ncx = nxk(klevel-1)
	  ncy = nyk(klevel-1)
	  irc = krbgn(klevel-1)
	  call rescd2sp(nx,ny,wk(ip),wk(ir),ncx,ncy,wk(ipc),wk(irc),
     +                   wk(icx),wk(icy),wk(kps))
c
c     prerelax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 10
	else
c
c    postrelax at coarsest level
c
	  do l=1,ipost
	    call relcd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                    wk(ity),wk(kps))
	  end do
	  ipc = ip
	  ncx = nx
	  ncy = ny
	  ip = kpbgn(2)
	  ir = krbgn(2)
	  icx = kcxbgn(2)
	  icy = kcybgn(2)
	  nx = nxk(2)
	  ny = nyk(2)
c
c    inject correction to level 2
c
	  call ccor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +               intpol,wk(kps))
c
c     set to postrelax at level 2
c
	  nrel = ipost
	  klevel = 2
	  go to 10
	end if
      end if
    5 continue
c
c     post relax at current finest grid level
c
      nx = nxk(kcur)
      ny = nyk(kcur)
      ip = kpbgn(kcur)
      ir = krbgn(kcur)
      icx = kcxbgn(kcur)
      icy = kcybgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
	call relcd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                 wk(ity),wk(kps))
      end do
      return
      end

      subroutine discd2sp(nx,ny,cofx,cofy,tx,ty,bndyc,cfx,cfy,wk,ier)
c
c     discretize elliptic pde for cud2, set nonfatal errors
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy,im1,jm1,ier
      complex cofx(nx,3),cofy(ny,3),tx(nx,ny,*),ty(ny,nx,*)
      complex wk(*)
      real x,y,dlx,dlx2,dlxx,dly,dly2,dlyy,alfmax,cemax
      complex cxx,cyy,cx,cy,cex,cey,c1,c2,c3,alfa,gbdy
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      external bndyc,cfx,cfy
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      alfmax = 0.0
      cemax = 0.0
c
c     set x,y subscript limits for calls to cofx,cofy
c     (avoid specified boundaries)
c
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
c
c     compute discretization coefficients on interior and
c     nonspecified boundaries
c
      do j=jst,jfn
	y = yc+(j-1)*dly
	call cfy(y,cyy,cy,cey)
	cemax = amax1(cabs(cey),cemax)
c
c     flag hyperbolic pde if necessary
c
	if (cabs(cy)*dly.gt.2*cabs(cyy)) then
	  if (klevel.eq.ngrid) then
	    ier = -4
	  end if
	  cyy = cmplx(0.5*cabs(cy)*dly,0.0)
	end if
	c1 = cyy/dlyy-cy/dly2
	c2 = cyy/dlyy+cy/dly2
	c3 = cey-(c1+c2)
	cofy(j,1) = c1
	cofy(j,2) = c2
	cofy(j,3) = c3
      end do
      do i=ist,ifn
	x = xa+(i-1)*dlx
	call cfx(x,cxx,cx,cex)
	cemax = amax1(cabs(cex),cemax)
c
c      flag hyperbolic pde if necessary
c
	if (cabs(cx)*dlx.gt.2*cabs(cxx)) then
	  if (klevel.eq.ngrid) then
	    ier = -4
	  end if
	  cxx = cmplx(0.5*cabs(cx)*dlx,0.0)
	end if
	c1 = cxx/dlxx-cx/dlx2
	c2 = cxx/dlxx+cx/dlx2
	c3 = cex-(c1+c2)
	cofx(i,1) = c1
	cofx(i,2) = c2
	cofx(i,3) = c3
      end do
c
c     adjust discretization for mixed derivative b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	i = 1
	c1 = cofx(i,1)
	cofx(i,1) = (0.0,0.0)
	cofx(i,2) = cofx(i,2)+c1
	y = yc+dly
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,y,alfa,gbdy)
	alfmax = amax1(alfmax,cabs(alfa))
	cofx(i,3) = cofx(i,3)+dlx2*alfa*c1
      end if
      if (nxb.eq.2) then
	kbdy = 2
	i = nx
	y = yc+dly
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,y,alfa,gbdy)
	c2 = cofx(i,2)
	cofx(i,1) = cofx(i,1)+c2
	cofx(i,2) = (0.0,0.0)
	cofx(i,3) = cofx(i,3)-dlx2*alfa*c2
	alfmax = amax1(cabs(alfa),alfmax)
      end if
      if (nyc.eq.2) then
	kbdy = 3
	j = 1
	x = xa+dlx
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,alfa,gbdy)
	c1 = cofy(j,1)
	cofy(j,1) = (0.0,0.0)
	cofy(j,2) = cofy(j,2) + c1
	cofy(j,3) = cofy(j,3) + dly2*alfa*c1
	alfmax = amax1(cabs(alfa),alfmax)
      end if
      if (nyd.eq.2) then
	kbdy = 4
	j = ny
	x = xa+dlx
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,alfa,gbdy)
	c2 = cofy(j,2)
	cofy(j,2) = (0.0,0.0)
	cofy(j,1) = cofy(j,1) + c2
	cofy(j,3) = cofy(j,3) - dly2*c2*alfa
	alfmax = amax1(cabs(alfa),alfmax)
      end if
c
c     if detected then flag singular pde
c
	if (cemax.eq.0.0.and.alfmax.eq.0.0) then
	  if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	    if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
	      ier = -3
	   end if
	  end if
	end if
c
c     set and factor tridiagonal matrices for line relaxation(s) if flagged
c
      if (method.eq.1.or.method.eq.3) then
	if (nxa.ne.0) then
c
c    nonperiodic x line relaxation
c
	  do i=1,nx
	    im1 = max0(i-1,1)
	    do j=1,ny
	      tx(im1,j,1) = cofx(i,1)
	      tx(i,j,2) = cofx(i,3)+cofy(j,3)
	      tx(i,j,3) = cofx(i,2)
	    end do
	  end do
	  if (nxa.eq.1) then
	    do j=1,ny
	      tx(1,j,2) = (1.0,0.0)
	      tx(1,j,3) = (0.0,0.0)
	    end do
	  end if
	  if (nxb.eq.1) then
	    do j=1,ny
	      tx(nx-1,j,1) = (0.0,0.0)
	      tx(nx,j,2) = (1.0,0.0)
	    end do
	  end if
	  call cfactri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
	else
c
c     periodic x line relaxation
c
	  if (nx .gt. 3) then
c
c     set and factor iff nx > 3
c
	    do i=1,nx-1
	      do j=1,ny
		tx(i,j,1) = cofx(i,1)
		tx(i,j,2) = cofx(i,3)+cofy(j,3)
		tx(i,j,3) = cofx(i,2)
	      end do
	    end do
	    call cfactrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
     +                   tx(1,1,5),wk(kps))
	  end if
	end if
      end if

      if (method.eq.2.or.method.eq.3) then
	if (nyc.ne.0) then
c
c     nonperiodic y line relaxation
c
	  do j=1,ny
	    jm1 = max0(j-1,1)
	    do i=1,nx
	      ty(jm1,i,1) = cofy(j,1)
	      ty(j,i,2) = cofy(j,3)+cofx(i,3)
	      ty(j,i,3) = cofy(j,2)
	    end do
	  end do
	  if (nyc.eq.1) then
	    do i=1,nx
	      ty(1,i,2) = (1.0,0.0)
	      ty(1,i,3) = (0.0,0.0)
	    end do
	  end if
	  if (nyd.eq.1) then
	    do i=1,nx
	      ty(ny-1,i,1) = (0.0,0.0)
	      ty(ny,i,2) = (1.0,0.0)
	    end do
	  end if
	  call cfactri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
	else
c
c      periodic y line relaxation
c
	  if (ny .gt. 3) then
c
c     set and factor iff ny > 3
c
	    do j=1,ny-1
	      do i=1,nx
		ty(j,i,1) = cofy(j,1)
		ty(j,i,2) = cofy(j,3)+cofx(i,3)
		ty(j,i,3) = cofy(j,2)
	      end do
	    end do
	    call cfactrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                   ty(1,1,5),wk(kps))
	  end if
	end if
      end if
      return
      end

      subroutine adjcd2sp(nx,ny,phi,rhs,bndyc,cfx,cfy)
c
c     adjust righthand side for various boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy
      complex phi(0:nx+1,0:ny+1),rhs(nx,ny)
      real x,y,dlx,dlx2,dlxx,dly,dly2,dlyy
      complex cxx,cyy,cx,cy,cex,cey,c1,c2,alfa,gbdy
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fcud2/xa,xb,yc,yd,tolmax,relmax
      external bndyc,cfx,cfy
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlx2 = dlx+dlx
      dly2 = dly+dly
      dlxx = dlx*dlx
      dlyy = dly*dly
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
c
c     adjust right hand side at derivative boundaries
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
	call cfx(x,cxx,cx,cex)
	if (cabs(cx)*dlx.gt.2*cabs(cxx)) then
	  cxx = cmplx(0.5*cabs(cx)*dlx,0.0)
	end if
	c1 = cxx/dlxx-cx/dlx2
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)+dlx2*c1*gbdy
	end do
      end if
      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
	call cfx(x,cxx,cx,cex)
	if (cabs(cx)*dlx.gt.2*cabs(cxx)) then
	  cxx = cmplx(0.5*cabs(cx)*dlx,0.0)
	end if
	c2 = cxx/dlxx+cx/dlx2
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)-dlx2*c2*gbdy
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
	call cfy(y,cyy,cy,cey)
	if (cabs(cy)*dly.gt.2*cabs(cyy)) then
	  cyy = cmplx(0.5*cabs(cy)*dly,0.0)
	end if
	c1 = cyy/dlyy-cy/dly2
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)+dly2*c1*gbdy
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
	call cfy(y,cyy,cy,cey)
	if (cabs(cy)*dly.gt.2*cabs(cyy)) then
	  cyy = cmplx(0.5*cabs(cy)*dly,0.0)
	end if
	c2 = cyy/dlyy+cy/dly2
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)-dly2*c2*gbdy
	end do
      end if
c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  rhs(i,j) = phi(i,j)
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  rhs(i,j) = phi(i,j)
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  rhs(i,j) = phi(i,j)
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  rhs(i,j) = phi(i,j)
	end do
      end if
      return
      end

      subroutine rescd2sp(nx,ny,phi,rhs,ncx,ncy,phic,rhsc,cofx,cofy,
     +                    resf)
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc,ist,ifn,jst,jfn
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      complex rhs(nx,ny),rhsc(ncx,ncy),resf(nx,ny)
      complex phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      complex cofx(nx,3),cofy(ny,3)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = (0.0,0.0)
	end do
      end do
c
c     intialize residual to zero and set limits
c
      do j=1,ny
	do i=1,nx
	  resf(i,j) = (0.0,0.0)
	end do
      end do
      ist = 1
      if (nxa.eq.1) ist = 2
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 2
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(jst,jfn,ist,ifn,rhs,resf,cofx,cofy,phi)
!$OMP PRIVATE(i,j)
      do j=jst,jfn
	do i=ist,ifn
	  resf(i,j) =  rhs(i,j)-(
     +                 cofx(i,1)*phi(i-1,j)+
     +                 cofx(i,2)*phi(i+1,j)+
     +                 cofy(j,1)*phi(i,j-1)+
     +                 cofy(j,2)*phi(i,j+1)+
     +                 (cofx(i,3)+cofy(j,3))*phi(i,j))
	end do
      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call cres2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

      subroutine relcd2sp(nx,ny,phi,rhs,cofx,cofy,tx,ty,sum)
c
c     relaxation for cud2sp
c
      implicit none
      integer nx,ny
      complex phi(*),rhs(*),cofx(*),cofy(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
	call rcd2spp(nx,ny,phi,rhs,cofx,cofy)
      else if (method.eq.1) then           ! line x relaxation
	call slxcd2sp(nx,ny,phi,rhs,cofx,cofy,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slycd2sp(nx,ny,phi,rhs,cofx,cofy,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxcd2sp(nx,ny,phi,rhs,cofx,cofy,tx,sum)
	call slycd2sp(nx,ny,phi,rhs,cofx,cofy,ty,sum)
      end if
      return
      end

      subroutine rcd2spp(nx,ny,phi,rhs,cofx,cofy)
c
c     gauss-seidel red/black point relaxation
c
      implicit none
      integer nx,ny,i,j,ist,ifn,jst,jfn
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),rhs(nx,ny),cofx(nx,3),cofy(ny,3)
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = 1
      if (nxa.eq.1) ist = 3
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 3
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
c
c    periodic adjustment bypass block
c
      if (nxa*nyc.ne.0) then
c
c     relax on red grid points
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=ist,ifn,2
	  do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) - (
     +                  cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                  cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=2,ifn,2
	  do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
c
c     relax on black grid points
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=ist,ifn,2
	  do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=2,ifn,2
	  do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
	return
      end if
c
c    set periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on red grid points
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=ist,ifn,2
	do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=2,ifn,2
	do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
c
c    ensure periodic virtual boundary red points are set
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on black grid points
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=ist,ifn,2
	do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=2,ifn,2
	do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
c
c     final set of periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      return
      end

      subroutine slxcd2sp(nx,ny,phi,rhs,cofx,cofy,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j,ist,ifn,jst,jfn,jm1,jp1
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),cofx(nx,3),cofy(ny,3),tx(nx,ny,*)
      complex sum(ny),rhs(nx,ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call rcd2spp(nx,ny,phi,rhs,cofx,cofy)
	return
      end if
      jst = 3
      jfn = ny-1
      ist = 2
      if (nxa.ne.1) ist = 1
      ifn = nx-1
      if (nxb.ne.1) ifn = nx
      if (nxa.ne.0) then
c
c     non-periodic line relaxation in x direction
c     lines thru odd y points
c
	if (nyc.ne.1) then
	  jst = 1
	  j = 1
	  jm1 = ny-1
	  jp1 = 2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,jm1)+cofy(j,2)*phi(i,jp1))
	  end do
	end if
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(ist,ifn,phi,rhs,cofy,ny) PRIVATE(i,j)
	do j=3,ny-1,2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1))
	  end do
	end do
	if (mod(ny,2).ne.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,jm1)+cofy(j,2)*phi(i,jp1))
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(jst,jfn,phi,tx,nx) PRIVATE(i,ib,j)
      do j=jst,jfn,2
	do i=2,nx
	  phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	end do
c
c     backward sweep
c
	phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	do ib=2,nx
	  i = nx-ib+1
	  phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	end do
      end do
c
c     lines thru even y points
c
	jfn = ny-1
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(ist,ifn,cofy,rhs,phi,ny) PRIVATE(i,j)
	do j=2,ny-1,2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1))
	  end do
	end do
	if (mod(ny,2).eq.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(jfn,phi,tx,nx) PRIVATE(i,ib,j)
	do j=2,jfn,2
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
c
c     backward sweep
c
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
c
c     adjust for periodic y if ny even (only possible at coarsest level)
c     red/black line will not allow phi(i,1) = phi(i,ny)
c
	if (nyc.eq.0.and.mod(ny,2).eq.0) then
	  do i=1,nx
	    phi(i,ny) = phi(i,1)
	  end do
	end if
      else
c
c     line periodic relaxation in x direction
c
	do j=1,ny
	  sum(j) = (0.0,0.0)
	end do
c
c     set rhs and solve on lines thru odd y points
c
	jst = 3
	jfn = ny-1
	if (nyc.ne.1) then
	  jst = 1
	  j = 1
	  jm1 = ny-1
	  jp1 = 2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofy,rhs,phi,nx,ny) PRIVATE(i,j)
	do j=3,ny-1,2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,j-1)-cofy(j,2)*phi(i,j+1)
	  end do
	end do
c
c     set last y point line if odd and non-specified
c
	if (mod(ny,2).ne.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,jst,jfn,phi,tx,nx) PRIVATE(i,ib,j)
	do j=jst,jfn,2
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))
     +                   /tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic point
c
	do j=jst,jfn,2
	  phi(nx,j) = phi(1,j)
	end do
c
c     set rhs and solve on lines thru even y points
c
	jfn = ny-1
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cofy,phi,rhs,nx,ny) PRIVATE(i,j)
	do j=2,ny-1,2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,j-1)-cofy(j,2)*phi(i,j+1)
	  end do
	end do
c
c     set last y point if even and non-specified
c
	if (mod(ny,2).eq.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,jfn,phi,tx,nx) PRIVATE(i,ib,j)
	do j=2,jfn,2
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))
     +                  /tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic point
c
	do j=2,jfn,2
	  phi(nx,j) = phi(1,j)
	end do
c
c     adjust for periodic y if ny even (only possible at coarsest level)
c     red/black line will not allow phi(i,1) = phi(i,ny)
c
	if (nyc.eq.0.and.mod(ny,2).eq.0) then
	  do i=1,nx
	    phi(i,ny) = phi(i,1)
	  end do
	end if
      end if
c
c      final set of periodic and virtual x and y boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      return
      end

      subroutine slycd2sp(nx,ny,phi,rhs,cofx,cofy,ty,sum)
      implicit none
      integer nx,ny,i,j,jb,ist,ifn,jst,jfn,im1,ip1
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),rhs(nx,ny)
      complex cofx(nx,3),cofy(ny,3),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call rcd2spp(nx,ny,phi,rhs,cofx,cofy)
	return
      end if
      ist = 3
      ifn = nx-1
      jst = 2
      if (nyc.ne.1) jst = 1
      jfn = ny-1
      if (nyd.ne.1) jfn = ny
      if (nyc.ne.0) then
c
c     non-periodic case
c     set adjusted rhs and solve on odd x points
c
	if (nxa.ne.1) then
	  i = 1
	  ist = 1
	  im1 = nx-1
	  ip1 = 2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     x interior
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(jst,jfn,phi,rhs,cofx,nx) PRIVATE(i,j)
	do i=3,nx-1,2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
	if (mod(nx,2).ne.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(ist,ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=ist,ifn,2
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c     backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     repeat for even x points
c
	ifn = nx-1
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(jst,jfn,phi,rhs,cofx,nx) PRIVATE(i,j)
	do i=2,nx-1,2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
	if (mod(nx,2).eq.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=2,ifn,2
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c     backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     adjust x boundaries if nxa=0 and nx is even (only possible at coar
c     red/black line ordering above will not allow phi(1,j) = phi(nx,j)
c
	if (nxa.eq.0.and.mod(nx,2).eq.0) then
	  do j=1,ny
	    phi(nx,j) = phi(1,j)
	  end do
	end if
      else
c
c     line periodic relaxation in y direction
c
	do i=1,nx
	  sum(i) = (0.0,0.0)
	end do
c
c     set rhs and solve on odd x points
c
	ist = 3
	ifn = nx-1
	if (nxa.ne.1) then
	  ist = 1
	  i = 1
	  im1 = nx-1
	  ip1 = 2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(phi,rhs,cofx,nx,ny) PRIVATE(i,j)
	do i=3,nx-1,2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
c
c     set last x point if odd and non-specified
c
	if (mod(nx,2).ne.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(ist,ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=ist,ifn,2
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                 phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
c
c     set periodic point
c
	do i=ist,ifn,2
	  phi(i,ny) = phi(i,1)
	end do
c
c     set rhs and solve on even x points
c
	ifn = nx-1
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(phi,rhs,cofx,nx,ny) PRIVATE(i,j)
	do i=2,nx-1,2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
c
c     set last x line if even and non-specified
c
	if (mod(nx,2).eq.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=2,ifn,2
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
c
c     set periodic point
c
	do i=2,ifn,2
	  phi(i,ny) = phi(i,1)
	end do
c
c     adjust for periodic x if nx even (only possible at coarsest level)
c     red/black line will not allow phi(1,j) = phi(nx,j)
c
	if (nxa.eq.0.and.mod(nx,2).eq.0) then
	  do j=1,ny
	    phi(nx,j) = phi(1,j)
	  end do
	end if
      end if
c
c      set periodic and virtual x and y boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      return
      end
