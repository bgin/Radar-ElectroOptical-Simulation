c
c     file mud2sa.f
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
c     mud2sa attempts to produce a second order finite difference
c     approximation to the elliptic partial differential equation of the form:
c
c     (sigx(x,y)*px)x + (sigy(x,y)*py)y - xlmbda*p(x,y) = r(x,y)
c
c ... see documentation and test files provided in this distribution
c
c ... required mudpack files
c
c     mudcom.f
c
c
      subroutine mud2sa(iparm,fparm,work,sigx,sigy,xlmbda,bndyc,rhs,phi,
     +                  mgopt,ierror)
      implicit none
      integer iparm,mgopt,ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,ic,itx,ity
      dimension iparm(17),fparm(6),mgopt(4)
      real work(*),phi(*),rhs(*)
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      common/mud2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      real sigx,sigy,xlmbda
      external sigx,sigy,xlmbda,bndyc
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
	  kps = kps+(nx+2)*(ny+2)+nx*ny*(6+isx+jsy)
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
	  kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
	  ktxbgn(k) = kcbgn(k)+6*nx*ny
	  ktybgn(k) = ktxbgn(k)+isx*nx*ny
	  iw = ktybgn(k)+jsy*nx*ny
c         ip = kpbgn(k)
	  ic = kcbgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call dismd2sa(nx,ny,work(ic),work(itx),work(ity),
     +                bndyc,sigx,sigy,xlmbda,work,ierror)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call mud2sa1(nx,ny,rhs,phi,sigx,sigy,xlmbda,bndyc,work)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end

      subroutine mud2sa1(nx,ny,rhsf,phif,sigx,sigy,xlmbda,bndyc,wk)
      implicit none
      integer nx,ny
      real phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      common/mud2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      real sigx,sigy,xlmbda
      external sigx,sigy,xlmbda,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+5*nx*ny
c
c     set phif,rhsf in wk and adjust right hand side
c
      call swk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
c
c     no initial guess at finest grid level!
c
	do kb=2,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ir = kcbgn(k+1)+5*nx*ny
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ipc = kpbgn(k)
	  icc = kcbgn(k)
	  irc = icc+5*ncx*ncy
c
c     transfer down to all grid levels
c
	  call trsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,
     +                wk(ipc),wk(irc))
	end do
c
c     adjust right hand side at all grid levels in case
c     rhs or specified b.c. in phi or gbdy changed
c
	do k=1,ngrid
	  nx = nxk(k)
	  ny = nyk(k)
	  ip = kpbgn(k)
	  ic = kcbgn(k)
	  call adjmd2sa(nx,ny,wk(ip),wk(ic),bndyc,sigx,sigy,xlmbda)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcymd2sa(wk)
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ipc = kpbgn(k)
	  ncx = nxk(k)
	  ncy = nyk(k)
c
c     lift or prolong approximation from k to k+1
c
	  call prolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,
     +                 nyc,nyd,intpol)
	end do
      else
c
c     adjust rhs at finest grid level only
c
	nx = nxk(ngrid)
	ny = nyk(ngrid)
	ip = kpbgn(ngrid)
	ic = kcbgn(ngrid)
	call adjmd2sa(nx,ny,wk(ip),wk(ic),bndyc,sigx,sigy,xlmbda)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcymd2sa(wk)
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
	      phmax = amax1(phmax,abs(wk(ij)))
	      relmax = amax1(relmax,abs(wk(ij)-phif(i,j)))
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

      subroutine kcymd2sa(wk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      real wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      common/mud2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
c
c     prerelax at current finest grid level
c
      do l=1,iprer
	call relmd2sa(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+5*ncx*ncy
      call resmd2sa(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
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
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
	call relmd2sa(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle complete at klevel
c
	ipc = ip
	ip = kpbgn(klevel+1)
	ncx = nxk(klevel)
	ncy = nyk(klevel)
	nx = nxk(klevel+1)
	ny = nyk(klevel+1)
c
c    inject correction to finer grid
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
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
	  irc = kcbgn(klevel-1)+5*ncx*ncy
	  call resmd2sa(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     +                wk(kps))
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
	    call relmd2sa(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
	  end do
	  ipc = ip
	  ip = kpbgn(2)
	  ncx = nxk(1)
	  ncy = nyk(1)
	  nx = nxk(2)
	  ny = nyk(2)
c
c    inject correction to level 2
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
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
      ic = kcbgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
	call relmd2sa(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end

      subroutine dismd2sa(nx,ny,cof,tx,ty,bndyc,sigx,sigy,xlmbda,wk,ier)
c
c     discretize selfadjoint elliptic pde for mud2sa, set nonfatal errors
c     use conservative differencing
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy,l,im1,jm1,ier
      real cof(nx,ny,6),tx(nx,ny,*),ty(ny,nx,*)
      real wk(*),dlx,dlx2,dlxx,dly,dly2,dlyy,odlxx,odlyy,alfmax
      real x,y,c1,c2,c3,c4,c5,alfa,gbdy
      real sigmin,xlmin,xlmax,sgxm,sgxp,sgym,sgyp,xlm
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      real sigx,sigy,xlmbda
      external bndyc,sigx,sigy,xlmbda
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      odlxx = 1.0/dlxx
      odlyy = 1.0/dlyy
      alfmax = 0.0
c
c     set x,y subscript limits for calls to sigx,sigy,xlmbda,bndyc
c
      jst = 1
      jfn = ny
      ist = 1
      ifn = nx
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
c
c     compute discretization coefficients on interior and
c     nonspecified boundaries
c
      sigmin = 1.0
      xlmin = 1.0
      xlmax = 0.0
      do j=jst,jfn
	y = yc+(j-1)*dly
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  sgxm = sigx(x-0.5*dlx,y)
	  sgxp = sigx(x+0.5*dlx,y)
	  sgym = sigy(x,y-0.5*dly)
	  sgyp = sigy(x,y+0.5*dly)
	  xlm = xlmbda(x,y)
c
c     monitor possible nonfatal coefficient errors
c
	  sigmin = amin1(sigmin,sigx(x,y),sigy(x,y))
	  xlmin = amin1(xlmin,xlm)
	  xlmax = amax1(abs(xlm),xlmax)
	  c1 = sgxm*odlxx
	  c2 = sgxp*odlxx
	  c3 = sgym*odlyy
	  c4 = sgyp*odlyy
	  c5 = -xlm-(c1+c2+c3+c4)
	  cof(i,j,1) = c1
	  cof(i,j,2) = c2
	  cof(i,j,3) = c3
	  cof(i,j,4) = c4
	  cof(i,j,5) = c5
	end do
      end do
c
c     flag non-ellipticity
c
      if (sigmin.le.0.0) then
	ier = -2
      end if
c
c     set warning flag if xlmbda is less then zero at some point
c
      if (xlmin.lt.0.0) then
	if (ier.eq.0) then
	  ier = -4
	end if
      end if
c
c     adjust discretization for mixed derivative b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  c1 = cof(i,j,1)
	  cof(i,j,1) = 0.0
	  cof(i,j,2) = cof(i,j,2)+c1
	  cof(i,j,5) = cof(i,j,5)+dlx2*alfa*c1
	  alfmax = amax1(abs(alfa),alfmax)
	end do
      end if
      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  c2 = cof(i,j,2)
	  cof(i,j,1) = cof(i,j,1)+c2
	  cof(i,j,2) = 0.0
	  cof(i,j,5) = cof(i,j,5)-dlx2*alfa*c2
	  alfmax = amax1(abs(alfa),alfmax)
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  c3 = cof(i,j,3)
	  cof(i,j,3) = 0.0
	  cof(i,j,4) = cof(i,j,4)+c3
	  cof(i,j,5) = cof(i,j,5)+dly2*alfa*c3
	  alfmax = amax1(abs(alfa),alfmax)
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  c4 = cof(i,j,4)
	  cof(i,j,3) = cof(i,j,3)+c4
	  cof(i,j,4) = 0.0
	  cof(i,j,5) = cof(i,j,5)-dly2*c4*alfa
	  alfmax = amax1(abs(alfa),alfmax)
	end do
      end if
	if (xlmax.eq.0.0.and.alfmax.eq.0.0) then
	  if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	    if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
c
c     flag singular pde fatal error
c
	      ier = -3
	      return
	   end if
	  end if
	end if
c
c     set coefficient for specified boundaries
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
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
	      tx(im1,j,1) = cof(i,j,1)
	      tx(i,j,2) = cof(i,j,5)
	      tx(i,j,3) = cof(i,j,2)
	    end do
	  end do
	  call factri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
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
		tx(i,j,1) = cof(i,j,1)
		tx(i,j,2) = cof(i,j,5)
		tx(i,j,3) = cof(i,j,2)
	      end do
	    end do
	    call factrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
     +                  tx(1,1,5),wk(kps))
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
	      ty(jm1,i,1) = cof(i,j,3)
	      ty(j,i,2) = cof(i,j,5)
	      ty(j,i,3) = cof(i,j,4)
	    end do
	  end do
	  call factri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
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
		ty(j,i,1) = cof(i,j,3)
		ty(j,i,2) = cof(i,j,5)
		ty(j,i,3) = cof(i,j,4)
	      end do
	    end do
	    call factrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                  ty(1,1,5),wk(kps))
	  end if
	end if
      end if
      return
      end


      subroutine adjmd2sa(nx,ny,phi,cof,bndyc,sigx,sigy,xlmbda)
c
c     adjust righthand side in cof(i,j,6) for boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,sgxm,sgxp,sgym,sgyp
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6)
      real dlx,dlx2,dlxx,dly,dly2,dlyy,odlxx,odlyy
      real x,y,c1,c2,c3,c4,alfa,gbdy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      real sigx,sigy,xlmbda
      external bndyc,sigx,sigy,xlmbda
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlx2 = dlx+dlx
      dly2 = dly+dly
      dlxx = dlx*dlx
      dlyy = dly*dly
      odlxx = 1.0/dlxx
      odlyy = 1.0/dlyy
      jst = 1
      jfn = ny
      ist = 1
      ifn = nx
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
c
c     adjust at derivative boundaries
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  sgxm = sigx(x-0.5*dlx,y)
	  c1 = sgxm*odlxx
	  cof(i,j,6) = cof(i,j,6)+dlx2*c1*gbdy
	end do
      end if
      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  sgxp = sigx(x+0.5*dlx,y)
	  c2 = sgxp*odlxx
	  cof(i,j,6) = cof(i,j,6)-dlx2*c2*gbdy
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  sgym = sigy(x,y-0.5*dly)
	  c3 = sgym*odlyy
	  cof(i,j,6) = cof(i,j,6)+dly2*c3*gbdy
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  sgyp = sigy(x,y+0.5*dly)
	  c4 = sgyp*odlyy
	  cof(i,j,6) = cof(i,j,6)-dly2*c4*gbdy
	end do
      end if
c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      return
      end

      subroutine resmd2sa(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real rhsc(ncx,ncy),resf(nx,ny)
      real phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real cof(nx,ny,6)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = 0.0
	end do
      end do
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
	do i=1,nx
	  resf(i,j) =  cof(i,j,6)-(
     +             cof(i,j,1)*phi(i-1,j)+
     +             cof(i,j,2)*phi(i+1,j)+
     +             cof(i,j,3)*phi(i,j-1)+
     +             cof(i,j,4)*phi(i,j+1)+
     +             cof(i,j,5)*phi(i,j))
	end do
      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

      subroutine relmd2sa(nx,ny,phi,cof,tx,ty,sum)
c
c     relaxation for mud2sa
c
      implicit none
      integer nx,ny
      real phi(*),cof(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
	call relmd2sap(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
	call slxmd2sa(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slymd2sa(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxmd2sa(nx,ny,phi,cof,tx,sum)
	call slymd2sa(nx,ny,phi,cof,ty,sum)
      end if
      return
      end

      subroutine relmd2sap(nx,ny,phi,cof)
c
c     gauss-seidel red/black point relaxation
c
      implicit none
      integer nx,ny,i,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6)
c
c    periodic adjustment bypass block
c
      if (nxa*nyc.ne.0) then
c
c     relax on red grid points
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=1,nx,2
	  do j=1,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=2,nx,2
	  do j=2,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
c
c     relax on black grid points
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=1,nx,2
	  do j=2,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=2,nx,2
	  do j=1,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
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
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=1,nx,2
	do j=1,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=2,nx,2
	do j=2,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
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
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=1,nx,2
	do j=2,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=2,nx,2
	do j=1,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
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

      subroutine slxmd2sa(nx,ny,phi,cof,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),tx(nx,ny,*),sum(ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call relmd2sap(nx,ny,phi,cof)
	return
      end if
c
c     set periodic y virtual boundary if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if

      if (nxa.ne.0) then
c
c     x direction not periodic
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
c
c     sweep odd j lines
c
	do j=1,ny,2
	  do i=1,nx
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
c
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
c     sweep even j lines forward and back
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=2,ny,2
	  do i=1,nx
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
      else
c
c     x direction periodic
c
	do j=1,ny
	  sum(j) = 0.0
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c      sweep odd lines forward and back
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=1,ny,2
	  do i=1,nx-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
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
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic and virtual points for j odd
c
	do j=1,ny,2
	  phi(nx,j) = phi(1,j)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c     sweep even j lines
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=2,ny,2
	  do i=1,nx-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
c
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
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic and virtual points for j even
c
	do j=2,ny,2
	  phi(nx,j) = phi(1,j)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
c
c     set periodic y virtual boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if

      return
      end

      subroutine slymd2sa(nx,ny,phi,cof,ty,sum)
      implicit none
      integer nx,ny,i,j,jb
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call relmd2sap(nx,ny,phi,cof)
	return
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if

      if (nyc.ne.0) then
c
c     y direction not periodic
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
c
c     sweep odd x lines
c
	do i=1,nx,2
	  do j=1,ny
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
	  end do
c
c     forward sweep thru odd x lines
c
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     forward sweep even x lines
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
	  end do
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
      else
c
c     y direction periodic
c
	do i=1,nx
	  sum(i) = 0.0
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
c     forward sweep odd x lines
c
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,ty,sum,nx,ny) PRIVATE(i,j,jb)
	do i=1,nx,2
	  do j=1,ny-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
	  end do
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
c       set odd periodic and virtual y boundaries
c
	do i=1,nx,2
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
c     forward sweep even x lines
c
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)

	  end do
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
c       set even periodic and virtual y boundaries
c
	do i=2,nx,2
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if

      return
      end
