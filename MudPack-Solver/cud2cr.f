c
c     file cud2cr.f
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
c     cud2cr attempts to produce a second order finite difference
c     approximation to the 2-d complex nonseparable elliptic
c     partial differential equation with cross derivative
c
c       cxx(x,y)*pxx + cxy(x,y)*pxy + cyy(x,y)*pyy +
c
c       cx(x,y)*px + cy(x,y)*py + ce(x,y)*p(x,y) = r(x,y)
c
c     The coefficients cxx,cxy,cyy,cx,cy,ce and the right-hand side
c     r and the unknown solution function p are all complex valued
c     functions of the real independent variables x and y.
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cudcom.f
c
      subroutine cud2cr(iparm,fparm,work,coef,bndyc,rhs,phi,mgopt,
     +                  ierror)
      implicit none
      integer iparm,mgopt,ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,ic,itx,ity
      dimension iparm(17),fparm(6),mgopt(4)
      complex work(*),phi(*),rhs(*)
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      common/cud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      external coef,bndyc
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
	  kps = kps+(nx+2)*(ny+2)+nx*ny*(10+isx+jsy)
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
	  ktxbgn(k) = kcbgn(k)+10*nx*ny
	  ktybgn(k) = ktxbgn(k)+isx*nx*ny
	  iw = ktybgn(k)+jsy*nx*ny
	  ic = kcbgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call discd2cr(nx,ny,work(ic),work(itx),work(ity),
     +                  bndyc,coef,work,ierror)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call cud2cr1(nx,ny,rhs,phi,coef,bndyc,work)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end

      subroutine cud2cr1(nx,ny,rhsf,phif,coef,bndyc,wk)
      implicit none
      integer nx,ny
      complex phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      common/cud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      external coef,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+9*nx*ny
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
	  ir = kcbgn(k+1)+9*nx*ny
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ipc = kpbgn(k)
	  icc = kcbgn(k)
	  irc = icc+9*ncx*ncy
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
	  ic = kcbgn(k)
	  call adjcd2cr(nx,ny,wk(ip),wk(ic),bndyc,coef)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcycd2cr(wk)
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
	ic = kcbgn(ngrid)
	call adjcd2cr(nx,ny,wk(ip),wk(ic),bndyc,coef)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcycd2cr(wk)
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

      subroutine kcycd2cr(wk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      complex wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      common/cud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
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
	call relcd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+9*ncx*ncy
      call rescd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
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
	call relcd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
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
	  irc = kcbgn(klevel-1)+9*ncx*ncy
	  call rescd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
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
	    call relcd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
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
	call ccor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +             intpol,wk(kps))
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
	call relcd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end

      subroutine discd2cr(nx,ny,cf,tx,ty,bndyc,coef,wk,ier)
c
c     discretize elliptic pde for cud2cr, set nonfatal errors
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j,kbdy,l,im1,jm1,ier,jc
      complex wk(*),cf(nx,ny,10),tx(nx,ny,*),ty(ny,nx,*)
      real x,y,dlx,dlx2,dlxx,dly,dly2,dlyy,cmin,gammax,cemax
      real dxoy,dyox,dlxy,dlxy2,dlxy4
      real cxxr,cxxi,cyyr,cyyi
      complex cxx,cxy,cyy,cx,cy,ce,c1,c2,c3,c4,c5
      complex c6,c7,c8,c9
      complex alfaa,alfab,alfac,alfad,betaa,betab,betac,betad,det
      complex gamaa,gamab,gamac,gamad
      complex alfim1,alfi,alfip1,betim1,beti,betip1,gamim1,gami,gamip1
      complex alfjm1,alfj,alfjp1,betjm1,betj,betjp1,gamjm1,gamj,gamjp1
      complex gbdim1,gbdi,gbdip1,gbdj,gbdjm1,gbdjp1
      complex gbdya,gbdyb,gbdyc,gbdyd
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      complex zer
      external bndyc,coef
      zer = (0.0,0.0)
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlxy = dlx*dly
      dlxy2 = dlxy+dlxy
      dlxy4 = dlxy2+dlxy2
      dyox = dly/dlx
      dxoy = dlx/dly
      cmin = 1.0
      gammax = 0.0
      cemax = 0.0
c
c     compute discretization coefficients on interior and
c     nonspecified boundaries
c
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call coef(x,y,cxx,cxy,cyy,cx,cy,ce)
	  cxxr = cxx
	  cxxi = aimag(cxx)
	  cyyr = cyy
	  cyyi = aimag(cyy)
	  cmin = amin1(cmin,amax1(cxxr*cyyr,cxxi*cyyi))
	  cemax = amax1(cabs(ce),cemax)
c
c     flag hyperbolic pde
c
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    if (klevel.eq.ngrid) ier = -4
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    if (klevel.eq.ngrid) ier = -4
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c1=cxx/dlxx+cx/dlx2
	  c2=cxy/dlxy4
	  c3=cyy/dlyy+cy/dly2
	  c4=-c2
	  c5=cxx/dlxx-cx/dlx2
	  c6=c2
	  c7=cyy/dlyy-cy/dly2
	  c8=-c2
	  c9=ce-(c1+c3+c5+c7)
	  cf(i,j,1)=c1
	  cf(i,j,2)=c2
	  cf(i,j,3)=c3
	  cf(i,j,4)=c4
	  cf(i,j,5)=c5
	  cf(i,j,6)=c6
	  cf(i,j,7)=c7
	  cf(i,j,8)=c8
	  cf(i,j,9)=c9
	end do
      end do
c
c     adjust at cornors virtual-virtual points not used
c
      call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
      c4=-cxy/dlxy2
      cf(1,1,1)=cf(1,1,1)-c4
      cf(1,1,2)=zer
      cf(1,1,3)=cf(1,1,3)-c4
      cf(1,1,4)=c4
      cf(1,1,5)=cf(1,1,5)-c4
      cf(1,1,6)=zer
      cf(1,1,7)=cf(1,1,7)-c4
      cf(1,1,8)=c4
      cf(1,1,9)=cf(1,1,9)+2.0*c4
      call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
      c2=cxy/dlxy2
      cf(nx,1,1)=cf(nx,1,1)-c2
      cf(nx,1,2)=c2
      cf(nx,1,3)=cf(nx,1,3)-c2
      cf(nx,1,4)=zer
      cf(nx,1,5)=cf(nx,1,5)-c2
      cf(nx,1,6)=c2
      cf(nx,1,7)=cf(nx,1,7)-c2
      cf(nx,1,8)=zer
      cf(nx,1,9)=cf(nx,1,9)+2.0*c2
      call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
      c2=cxy/dlxy2
      cf(1,ny,1)=cf(1,ny,1)-c2
      cf(1,ny,2)=c2
      cf(1,ny,3)=cf(1,ny,3)-c2
      cf(1,ny,4)=zer
      cf(1,ny,5)=cf(1,ny,5)-c2
      cf(1,ny,6)=c2
      cf(1,ny,7)=cf(1,ny,7)-c2
      cf(1,ny,8)=zer
      cf(1,ny,9)=cf(1,ny,9)+2.0*c2
      call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
      i = nx
      c4=-cxy/dlxy2
      cf(nx,ny,1)=cf(nx,ny,1)-c4
      cf(nx,ny,2)=zer
      cf(nx,ny,3)=cf(nx,ny,3)-c4
      cf(nx,ny,4)=c4
      cf(nx,ny,5)=cf(nx,ny,5)-c4
      cf(nx,ny,6)=zer
      cf(nx,ny,7)=cf(nx,ny,7)-c4
      cf(nx,ny,8)=c4
      cf(nx,ny,9)=cf(nx,ny,9)+2.0*c4
c
c     adjust discretization for mixed derivative b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
	call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
	gammax = amax1(cabs(gamjm1),cabs(gamj),gammax)
	do j=2,ny-1
	  jc = j
	  y = yc+j*dly
	  call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
	  gammax = amax1(cabs(gamjp1),gammax)
c
c     check for illegal tangential derivative b.c.
c
	  if (cabs(alfjm1*alfj*alfjp1).eq.0.0) then
	    ier = 13
	    return
	  end if
	  c4=cf(1,jc,4)
	  c5=cf(1,jc,5)
	  c6=cf(1,jc,6)
	  cf(1,jc,1)=cf(1,jc,1)+c5
	  cf(1,jc,2)=cf(1,jc,2)+c4
	  cf(1,jc,3)=cf(1,jc,3)+c6*(-betjm1/alfjm1*dxoy)+
     +    c5*(betj/alfj*dxoy)+c4*(3.*betjp1/alfjp1*dxoy+dlx2*gamjp1/
     +    alfjp1)
	  cf(1,jc,4)=zer
	  cf(1,jc,5)=zer
	  cf(1,jc,6)=zer
	  cf(1,jc,7)=cf(1,jc,7)+c6*(-3.*betjm1/alfjm1*dxoy+dlx2*gamjm1/
     +    alfjm1)+c5*(-betj/alfj*dxoy)+c4*(betjp1/alfjp1*dxoy)
	  cf(1,jc,8)=cf(1,jc,8)+c6
	  cf(1,jc,9)=cf(1,jc,9)+c6*(4.*betjm1/alfjm1*dxoy)+
     +    c5*(dlx2*gamj/alfj)+c4*(-4.*betjp1/alfjp1*dxoy)
	  alfjm1 = alfj
	  betjm1 = betj
	  gamjm1 = gamj
	  gbdjm1 = gbdj
	  alfj = alfjp1
	  betj = betjp1
	  gamj = gamjp1
	  gbdj = gbdjp1
	end do
      end if

      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
	call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
	gammax = amax1(cabs(gamjm1),cabs(gamj),gammax)
	do j=2,ny-1
	  jc = j
	  y = yc+j*dly
	  call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
	  gammax = amax1(cabs(gamjp1),gammax)
c
c     check for illegal tangential derivative b.c.
c
	  if (cabs(alfjm1*alfj*alfjp1).eq.0.0) then
	    ier = 13
	    return
	  end if
	  c1=cf(nx,jc,1)
	  c2=cf(nx,jc,2)
	  c8=cf(nx,jc,8)
	  cf(nx,jc,1)=zer
	  cf(nx,jc,2)=zer
	  cf(nx,jc,3)=cf(nx,jc,3)-c8*(-betjm1/alfjm1*dxoy)-c1*(betj/
     +    alfj*dxoy)-c2*(3.0*betjp1/alfjp1*dxoy+dlx2*gamjp1/alfjp1)
	  cf(nx,jc,4)=cf(nx,jc,4)+c2
	  cf(nx,jc,5)=cf(nx,jc,5)+c1
	  cf(nx,jc,6)=cf(nx,jc,6)+c8
	  cf(nx,jc,7)=cf(nx,jc,7)+c8*(3.0*betjm1/alfjm1*dxoy-dlx2*gamjm1
     +    /alfjm1)+c1*(betj/alfj*dxoy)+c2*(-betjp1/alfjp1*dxoy)
	  cf(nx,jc,8)=zer
	  cf(nx,jc,9)=cf(nx,jc,9)-c8*(4.0*betjm1/alfjm1*dxoy)-
     +    c1*(dlx2*gamj/alfj)-c2*(-4.0*betjp1/alfjp1*dxoy)
	  alfjm1 = alfj
	  betjm1 = betj
	  gamjm1 = gamj
	  gbdjm1 = gbdj
	  alfj = alfjp1
	  betj = betjp1
	  gamj = gamjp1
	  gbdj = gbdjp1
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
	jc = 1
	call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
	call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
	gammax = amax1(cabs(gamim1),cabs(gami),gammax)
	do i=2,nx-1
	  x=xa+i*dlx
	  call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
	  gammax = amax1(cabs(gamip1),gammax)
c
c     check for illegal tangential derivative b.c.
c
	  if (cabs(betim1*beti*betip1).eq.0.0) then
	    ier = 13
	    return
	  end if
	  c6=cf(i,jc,6)
	  c7=cf(i,jc,7)
	  c8=cf(i,jc,8)
	  cf(i,jc,1)=cf(i,jc,1)+c6*(-alfim1/betim1*dyox)+c7*(alfi/beti*
     +    dyox)+c8*(3.0*alfip1/betip1*dyox+dly2*gamip1/betip1)
	  cf(i,jc,2)=cf(i,jc,2)+c8
	  cf(i,jc,3)=cf(i,jc,3)+c7
	  cf(i,jc,4)=cf(i,jc,4)+c6
	  cf(i,jc,5)=cf(i,jc,5)+c6*(-3.0*alfim1/betim1*dyox+dly2*gamim1/
     +    betim1)+c7*(-alfi/beti*dyox)+c8*(alfip1/betip1*dyox)
	  cf(i,jc,6)=zer
	  cf(i,jc,7)=zer
	  cf(i,jc,8)=zer
	  cf(i,jc,9)=cf(i,jc,9)+c6*(4.0*alfim1/betim1*dyox)+c7*dly2*gami
     +    /beti+c8*(-4.0*alfip1/betip1*dyox)
c     advance scalars for next pass
	   alfim1=alfi
	   betim1=beti
	   gamim1=gami
	   gbdim1 = gbdi
	   alfi=alfip1
	   beti=betip1
	   gami=gamip1
	   gbdi = gbdip1
	end do
      end if

      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
	jc = ny
	kbdy=4
	call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
	call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
	gammax = amax1(cabs(gamim1),cabs(gami),gammax)
	do i=2,nx-1
	  x=xa+i*dlx
	  call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
	  gammax = amax1(cabs(gamip1),gammax)
c
c     check for illegal tangential derivative b.c.
c
	  if (cabs(betim1*beti*betip1).eq.0.0) then
	    ier = 13
	    return
	  end if
	  c2=cf(i,jc,2)
	  c3=cf(i,jc,3)
	  c4=cf(i,jc,4)
	  cf(i,jc,1)=cf(i,jc,1)+c4*(alfim1/betim1*dyox)+c3*(-alfi/beti*
     +    dyox)+c2*(-3.*alfip1/betip1*dyox-dly2*gamip1/betip1)
c     set virtual point coefficients to zero
	  cf(i,jc,2)=zer
	  cf(i,jc,3)=zer
	  cf(i,jc,4)=zer
	  cf(i,jc,5)=cf(i,jc,5)+c4*(3.0*alfim1/betim1*dyox-dly2*gamim1/
     +    betim1)+c3*(alfi/beti*dyox)+c2*(-alfip1/betip1*dyox)
	  cf(i,jc,6)=cf(i,jc,6)+c4
	  cf(i,jc,7)=cf(i,jc,7)+c3
	  cf(i,jc,8)=cf(i,jc,8)+c2
	  cf(i,jc,9)=cf(i,jc,9)+c4*(-4.*alfim1/betim1*dyox)+
     +    c3*(-dly2*gami/beti)+c2*(4.*alfip1/betip1*dyox)
	  alfim1=alfi
	  betim1=beti
	  gamim1=gami
	  gbdim1 = gbdi
	  alfi=alfip1
	  beti=betip1
	  gami=gamip1
	  gbdi = gbdip1
	end do
      end if
c
c     flag singular pde
c
	if (cemax.eq.0.0.and.gammax.eq.0.0) then
	  if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	    if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
	      ier = -3
	   end if
	  end if
	end if
c
c     flag non-ellipticity
c
      if (cmin.le.0.0) then
	ier = -2
      end if
c
c     make cornor adjustments if necessary
c
      if (nyc.eq.2) then
	j = 1
	jc = 1
	if (nxa.eq.0) then
c
c     periodic-mixed at (xa,yc)
c
	  kbdy = 3
	  call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
	  call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
	  c7=cf(1,jc,7)
	  c8=cf(1,jc,8)
	  cf(1,jc,1)=cf(1,jc,1)+c7*(alfi/beti*dyox)+
     +    c8*(3.0*alfip1/betip1*dyox+dly2*gamip1/betip1)
	  cf(1,jc,2)=cf(1,jc,2)+c8
	  cf(1,jc,3)=cf(1,jc,3)+c7
	  cf(1,jc,5)=cf(1,jc,5)+c7*(-alfi/beti*dyox)+
     +    c8*(alfip1/betip1*dyox)
	  cf(1,jc,6)=zer
	  cf(1,jc,7)=zer
	  cf(1,jc,8)=zer
	  cf(1,jc,9)=cf(1,jc,9)+c8*(-4.0*alfip1/betip1*dyox)+c7*(
     +    dly2*gami/beti)
c
c     adjust periodic-mixed at (xb,yc)
c
	  call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
	  call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
	  c6=cf(nx,jc,6)
	  c7=cf(nx,jc,7)
	  cf(nx,jc,1)=cf(nx,jc,1)+c6*(-alfim1/betim1*dyox)+c7*(alfi/beti
     +    *dyox)
	  cf(nx,jc,3)=cf(nx,jc,3)+c7
	  cf(nx,jc,4)=cf(nx,jc,4)+c6
	  cf(nx,jc,5)=cf(nx,jc,5)+c6*(-3.0*alfim1/betim1*dyox+dly2*
     +    gamim1/betim1)+c7*(-alfi/beti*dyox)
	  cf(nx,jc,6)=zer
	  cf(nx,jc,7)=zer
	  cf(nx,jc,8)=zer
	  cf(nx,jc,9)=cf(nx,jc,9)+c6*(4.0*alfim1/betim1*dyox)+c7*
     +    (dly2*gami/beti)
	else if (nxa.eq.2) then
c
c     mixed-mixed at (xa,yc)
c
	  kbdy = 3
c     phase 1
	  call bndyc(kbdy,xa+dlx,alfac,betac,gamac,gbdyc )
	  kbdy = 1
	  call bndyc(kbdy,yc+dly,alfaa,betaa,gamaa,gbdya)
	  c4=cf(1,jc,4)
	  c8=cf(1,jc,8)
	  cf(1,jc,1)=cf(1,jc,1)+c8*(3.0*alfac/betac*dyox+dly2*gamac/
     +    betac)
	  cf(1,jc,2)=cf(1,jc,2)+c8 +c4
	  cf(1,jc,3)=cf(1,jc,3)+c4*(3.0*betaa/alfaa*dxoy+dlx2*gamaa
     +    /alfaa)
	  cf(1,jc,4)=zer
	  cf(1,jc,5)=cf(1,jc,5)+c8*(alfac/betac*dyox)
	  cf(1,jc,6)=zer
	  cf(1,jc,7)=cf(1,jc,7)+c4*(betaa/alfaa*dxoy)
	  cf(1,jc,8)=zer
	  cf(1,jc,9)=cf(1,jc,9)+c8*(-4.0*alfac/betac*dyox)+
     +    c4*(-4.0*betaa/alfaa*dxoy)
c     phase 2
	  c5=cf(1,jc,5)
	  c7=cf(1,jc,7)
	  kbdy = 3
	  call bndyc(kbdy,xa,alfac,betac,gamac,gbdyc)
	  kbdy = 1
	  call bndyc(kbdy,yc,alfaa,betaa,gamaa,gbdya)
	  det=alfaa*betac-betaa*alfac
	  if (cabs(det).eq.0.0) then
	    ier = 14
	    return
	  end if
	  cf(1,jc,1)=cf(1,jc,1)+c5
	  cf(1,jc,3)=cf(1,jc,3)+c7
	  cf(1,jc,5)=zer
	  cf(1,jc,7)=zer
	  cf(1,jc,9)=cf(1,jc,9)+c5*(dlx2*(gamaa*betac-betaa*gamac)/det)+
     +    c7*(dly2*(alfaa*gamac-gamaa*alfac)/det)
	end if

	if (nxb.eq.2) then
c
c     mixed-mixed at (xb,yc)
c     phase 1
c
	  kbdy=2
	  call bndyc(kbdy,yc+dly,alfab,betab,gamab,gbdyb)
	  kbdy=3
	  call bndyc(kbdy,xb-dlx,alfac,betac,gamac,gbdyc)
	  c2=cf(nx,jc,2)
	  c6=cf(nx,jc,6)
	  cf(nx,jc,1)=cf(nx,jc,1)+c6*(-alfac/betac*dyox)
	  cf(nx,jc,2)=zer
	  cf(nx,jc,3)=cf(nx,jc,3)+c2*(-3.0*betab/alfab*dxoy-dlx2*
     +    gamab/alfab)
	  cf(nx,jc,4)=cf(nx,jc,4)+c6+c2
	  cf(nx,jc,5)=cf(nx,jc,5)+c6*(-3.0*alfac/betac*dyox+dly2*
     +    gamac/betac)
	  cf(nx,jc,6)=zer
	  cf(nx,jc,7)=cf(nx,jc,7)-c2*(betab/alfab*dxoy)
	  cf(nx,jc,8)=zer
	  cf(nx,jc,9)=cf(nx,jc,9)+c6*(4.0*alfac/betac*dyox)+
     +    c2*( 4.0*betab/alfab*dxoy)
c     phase 2
	  kbdy=3
	  call bndyc(kbdy,xb,alfac,betac,gamac,gbdyc)
	  kbdy=2
	  call bndyc(kbdy,yc,alfab,betab,gamab,gbdyb)
	  det=betac*alfab-alfac*betab
	  if (cabs(det).eq.0.0) then
	    ier = 14
	    return
	  end if
	  c1=cf(nx,jc,1)
	  c7=cf(nx,jc,7)
	  cf(nx,jc,1)=zer
	  cf(nx,jc,3)=cf(nx,jc,3)+c7
	  cf(nx,jc,5)=cf(nx,jc,5)+c1
	  cf(nx,jc,7)=zer
	  cf(nx,jc,9)=cf(nx,jc,9)+c1*(dlx2*(betab*gamac-gamab*betac)/
     +    det)+c7*(dly2*(alfab*gamac-gamab*alfac)/det)
	end if
      end if

      if (nyd.eq.2) then
	j = ny
	jc = ny
	if (nxa.eq.0) then
c     periodic-mixed at (xa,yd) and (xb,yd)
	  kbdy=4
	  call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
	  call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
	  c2=cf(1,jc,2)
	  c3=cf(1,jc,3)
	  cf(1,jc,1)=cf(1,jc,1)+c3*(-alfi/beti*dyox)+
     +    c2*(-3.0*alfip1/betip1*dyox-dly2*gamip1/betip1)
	  cf(1,jc,2) = zer
	  cf(1,jc,3) = zer
	  cf(1,jc,4) = zer
	  cf(1,jc,5)=cf(1,jc,5)+c3*(alfi/beti*dyox)+
     +    c2*(-alfip1/betip1*dyox)
	  cf(1,jc,7)=cf(1,jc,7)+c3
	  cf(1,jc,8)=cf(1,jc,8)+c2
	  cf(1,jc,9) = cf(1,jc,9)+c3*(-dly2*gami/beti)+
     +    c2*(4.0*alfip1/betip1*dyox)
	  call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
	  call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
	  c3=cf(nx,jc,3)
	  c4=cf(nx,jc,4)
	  cf(nx,jc,1)=cf(nx,jc,1)+c4*(alfim1/betim1*dyox)+
     +    c3*(-alfi/beti*dyox)
	  cf(nx,jc,2)=zer
	  cf(nx,jc,3)=zer
	  cf(nx,jc,4)=zer
	  cf(nx,jc,5)=cf(nx,jc,5)+c4*(3.0*alfim1/betim1*dyox-dly2*
     +    gamim1/betim1)+c3*(alfi/beti*dyox)
	  cf(nx,jc,6)=cf(nx,jc,6)+c4
	  cf(nx,jc,7)=cf(nx,jc,7)+c3
	  cf(nx,jc,9)=cf(nx,jc,9)+c4*(-4.0*alfim1/betim1*dyox)+
     +    c3*(-dly2*gami/beti)
	else if (nxa.eq.2) then
c     mixed-mixed at (xa,yd)
c     phase 1
	  kbdy=4
	  call bndyc(kbdy,xa+dlx,alfad,betad,gamad,gbdyd)
	  kbdy=1
	  call bndyc(kbdy,yd-dly,alfaa,betaa,gamaa,gbdya)
	  c2=cf(1,jc,2)
	  c6=cf(1,jc,6)
	  cf(1,jc,1)=cf(1,jc,1)+c2*(-3.0*alfad/betad*dyox-dly2*
     +    gamad/betad)
	  cf(1,jc,2)=zer
	  cf(1,jc,3)=cf(1,jc,3)+c6*(-betaa/alfaa*dxoy)
	  cf(1,jc,4)=zer
	  cf(1,jc,5)=cf(1,jc,5)+c2*(-alfad/betad*dyox)
	  cf(1,jc,6)=zer
	  cf(1,jc,7)=cf(1,jc,7)+c6*(-3.0*betaa/alfaa*dxoy+dlx2*
     +    gamaa/alfaa)
	  cf(1,jc,8)=cf(1,jc,8)+c6+c2
	  cf(1,jc,9)=cf(1,jc,9)+c2*(4.0*alfad/betad*dyox)+
     +    c6*(4.0*betaa/alfaa*dxoy)
c     phase 2
	  kbdy=1
	  call bndyc(kbdy,yd,alfaa,betaa,gamaa,gbdya)
	  kbdy=4
	  call bndyc(kbdy,xa,alfad,betad,gamad,gbdyd)
	  det=alfad*betaa-betad*alfaa
	  if (cabs(det).eq.0.0) then
	    ier = 14
	    return
	  end if
	  c3=cf(1,jc,3)
	  c5=cf(1,jc,5)
	  cf(1,jc,1)=cf(1,jc,1)+c5
	  cf(1,jc,3)=zer
	  cf(1,jc,5)=zer
	  cf(1,jc,7)=cf(1,jc,7)+c3
	  cf(1,jc,9)=cf(1,jc,9)+c5*(dlx2*(betaa*gamad-gamaa*betad)/det)+
     +    c3*(dly2*(alfaa*gamad-gamaa*alfad)/det)
	end if

	if (nxb.eq.2) then
c     mixed-mixed at (xb,yd)
	  kbdy=4
	  call bndyc(kbdy,xb-dlx,alfad,betad,gamad,gbdyd)
	  kbdy=2
	  call bndyc(kbdy,yd-dly,alfab,betab,gamab,gbdyb)
	  c4=cf(nx,jc,4)
	  c8=cf(nx,jc,8)
	  cf(nx,jc,1)=cf(nx,jc,1)+c4*(alfad/betad*dyox)
	  cf(nx,jc,2) = zer
	  cf(nx,jc,3) = cf(nx,jc,3)+c8*(betab/alfab*dxoy)
	  cf(nx,jc,4) = zer
	  cf(nx,jc,5) = cf(nx,jc,5)+c4*(3.0*alfad/betad*dyox-dly2*
     +    gamad/betad)
	  cf(nx,jc,6) = cf(nx,jc,6)+c4+c8
	  cf(nx,jc,7) = cf(nx,jc,7)+c8*(3.0*betab/alfab*dxoy-dlx2*
     +    gamab/alfab)
	  cf(nx,jc,8) = zer
	  cf(nx,jc,9) = cf(nx,jc,9)+c4*(-4.0*alfad/betad*dyox)+
     +    c8*(-4.0*betab/alfab*dxoy)
c     phase 2
	  kbdy=4
	  call bndyc(kbdy,xb,alfad,betad,gamad,gbdyd)
	  kbdy=2
	  call bndyc(kbdy,yd,alfab,betab,gamab,gbdyb)
	  det=alfad*betab-betad*alfab
	  if (cabs(det).eq.0.0) then
	    ier = 14
	    return
	  end if
	  c1=cf(nx,jc,1)
	  c3=cf(nx,jc,3)
	  cf(nx,jc,1)=zer
	  cf(nx,jc,3)=zer
	  cf(nx,jc,5)=cf(nx,jc,5)+c1
	  cf(nx,jc,7)=cf(nx,jc,7)+c3
	  cf(nx,jc,9)=cf(nx,jc,9)+c1*(dlx2*(betad*gamab-gamad*betab)/
     +    det)+c3*(dly2*(alfab*gamad-gamab*alfad)/det)
	end if
      end if

      if (nxa.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xa,yc)
c
	kbdy = 1
	call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
	call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	c5 = cf(1,1,5)
	c4 = cf(1,1,4)
	cf(1,1,3) = cf(1,1,3)+c5*(betj/alfj*dxoy)+
     +  c4*(3.*betjp1/alfjp1*dxoy+dlx2*gamjp1/alfjp1)
	cf(1,1,2) = cf(1,1,2)+c4
	cf(1,1,1) = cf(1,1,1)+c5
	cf(1,1,7) = cf(1,1,7)+c5*(-betj/alfj*dxoy)+
     +  c4*(betjp1/alfjp1*dxoy)
	cf(1,1,6) = zer
	cf(1,1,5) = zer
	cf(1,1,4) = zer
	cf(1,1,9) = cf(1,1,9)+c4*(-4.*betjp1/alfjp1*dxoy)+c5*(dlx2*
     +  gamj/alfj)
c
c     adjust mixed-periodic at (xa,yd)
c
	call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
	c6=cf(1,ny,6)
	c5=cf(1,ny,5)
	cf(1,ny,3) = cf(1,ny,3)+c6*(-betjm1/alfjm1*dxoy)+c5*(betj/alfj*
     +  dxoy)
	cf(1,ny,1) = cf(1,ny,1)+c5
	cf(1,ny,8) = cf(1,ny,8)+c6
	cf(1,ny,7) = cf(1,ny,7)+c6*(-3.*betjm1/alfjm1*dxoy+dlx2*
     +  gamjm1/alfjm1)+c5*(-betj/alfj*dxoy)
	cf(1,ny,6) = zer
	cf(1,ny,5) = zer
	cf(1,ny,4) = zer
	cf(1,ny,9) = cf(1,ny,9)+c6*(4.*betjm1/alfjm1*dxoy)+c5*
     +  (dlx2*gamj/alfj)
      end if

      if (nxb.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xb,yc) and (xb,yd)
c
	kbdy = 2
	call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
	call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	c2 = cf(nx,1,2)
	c1 = cf(nx,1,1)
	cf(nx,1,3) = cf(nx,1,3)+c1*(-betj/alfj*dxoy)+
     +  c2*(-3.*betjp1/alfjp1*dxoy-dlx2*gamjp1/alfjp1)
	cf(nx,1,2) = zer
	cf(nx,1,1) = zer
	cf(nx,1,8) = zer
	cf(nx,1,7) = cf(nx,1,7)+c1*(betj/alfj*dxoy)+
     +  c2*(-betjp1/alfjp1*dxoy)
	cf(nx,1,5) = cf(nx,1,5)+c1
	cf(nx,1,4) = cf(nx,1,4)+c2
	cf(nx,1,9) = cf(nx,1,9)+c1*(-dlx2*gamj/alfj)+
     +  c2*(4.*betjp1/alfjp1*dxoy)
	call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
	c1 = cf(nx,ny,1)
	c8 = cf(nx,ny,8)
	cf(nx,ny,3) = cf(nx,ny,3)+c8*(betjm1/alfjm1*dxoy)+
     +  c1*(-betj/alfj*dxoy)
	cf(nx,ny,2) = zer
	cf(nx,ny,1) = zer
	cf(nx,ny,8) = zer
	cf(nx,ny,7) = cf(nx,ny,7)+c8*(3.*betjm1/alfjm1*dxoy-dlx2*
     +  gamjm1/alfjm1)+c1*(betj/alfj*dxoy)
	cf(nx,ny,6) = cf(nx,ny,6)+c8
	cf(nx,ny,5) = cf(nx,ny,5)+c1
	cf(nx,ny,9) = cf(nx,ny,9)+c8*(-4.*betjm1/alfjm1*dxoy)+
     +  c1*(-dlx2*gamj/alfj)
      end if
c
c     set coefficient for specified boundaries
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  do l=1,9
	    cf(i,j,l) = zer
	  end do
	  cf(i,j,9) = (1.0,0.0)
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  do l=1,9
	    cf(i,j,l) = zer
	  end do
	  cf(i,j,9) = (1.0,0.0)
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  do l=1,9
	    cf(i,j,l) = zer
	  end do
	  cf(i,j,9) = (1.0,0.0)
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  do l=1,9
	    cf(i,j,l) = zer
	  end do
	  cf(i,j,9) = (1.0,0.0)
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
	      tx(im1,j,1) = cf(i,j,5)
	      tx(i,j,2) = cf(i,j,9)
	      tx(i,j,3) = cf(i,j,1)
	    end do
	  end do
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
		tx(i,j,1) = cf(i,j,5)
		tx(i,j,2) = cf(i,j,9)
		tx(i,j,3) = cf(i,j,1)
	      end do
	    end do
	    call cfactrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
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
	      ty(jm1,i,1) = cf(i,j,7)
	      ty(j,i,2) = cf(i,j,9)
	      ty(j,i,3) = cf(i,j,3)
	    end do
	  end do
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
		ty(j,i,1) = cf(i,j,7)
		ty(j,i,2) = cf(i,j,9)
		ty(j,i,3) = cf(i,j,3)
	      end do
	    end do
	    call cfactrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                  ty(1,1,5),wk(kps))
	  end if
	end if
      end if
      return
      end

      subroutine adjcd2cr(nx,ny,phi,cf,bndyc,coef)
c
c     adjust righthand side in cf(i,j,10) for boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j,kbdy
      complex cf(nx,ny,10),phi(0:nx+1,0:ny+1)
      real x,y,dlx,dlx2,dlxx,dly,dly2,dlyy,dlxy,dlxy2,dlxy4,dxoy,dyox
      complex cxx,cxy,cyy,cx,cy,ce,c1,c2,c3,c4,c5
      complex c6,c7,c8
      complex alfaa,alfab,alfac,alfad,betaa,betab,betac,betad,det
      complex gamaa,gamab,gamac,gamad
      complex alfim1,alfi,alfip1,betim1,beti,betip1,gamim1,gami,gamip1
      complex alfjm1,alfj,alfjp1,betjm1,betj,betjp1,gamjm1,gamj,gamjp1
      complex gbdim1,gbdi,gbdip1,gbdj,gbdjm1,gbdjp1
      complex gbdya,gbdyb,gbdyc,gbdyd
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      external bndyc,coef
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlxy = dlx*dly
      dlxy2 = dlxy+dlxy
      dlxy4 = dlxy2+dlxy2
      dyox = dly/dlx
      dxoy = dlx/dly
c
c     adjust at derivative boundaries
c
      if (nyc.eq.2) then
	kbdy=3
	call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
	call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
	do i=2,nx-1
	  x=xa+i*dlx
	  call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
	  call coef(x-dlx,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c6=cxy/dlxy4
	  c7=cyy/dlyy-cy/dly2
	  c8=-c6
	  cf(i,1,10)=cf(i,1,10)+dly2*(c6*gbdim1/betim1+c7*gbdi/beti+
     +    c8*gbdip1/betip1)
	  betim1=beti
	  gbdim1=gbdi
	  beti=betip1
	  gbdi=gbdip1
	end do
	if (nxa.eq.0) then
	  kbdy=3
	  call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
	  call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
	  call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c7=cyy/dlyy-cy/dly2+cxy/dlxy2
	  c8=-cxy/dlxy2
	  cf(1,1,10)=cf(1,1,10)+dly2*(c8*gbdip1/betip1+c7*gbdi/beti)
	  call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
	  call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
	  call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c6=cxy/dlxy2
	  c7=cyy/dlyy-cy/dly2-c6
	  cf(nx,1,10)=cf(nx,1,10)+dly2*(c6*gbdim1/betim1+c7*gbdi/beti)
	else if (nxa.eq.2) then
	  kbdy=3
	  call bndyc(kbdy,xa+dlx,alfac,betac,gamac,gbdyc)
	  kbdy=1
	  call bndyc(kbdy,yc+dly,alfaa,betaa,gamaa,gbdya)
	  call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c4=-cxy/dlxy2
	  c8=c4
	  c5=cxx/dlxx-cx/dlx2-c4
	  c7=cyy/dlyy-cy/dly2-c4
	  c5=c5+c8*(alfac/betac*dyox)
	  c7=c7+c4*(betaa/alfaa*dxoy)
	  cf(1,1,10)=cf(1,1,10)+2.0*(c8*dly*gbdyc/betac+c4*dlx*gbdya/
     +    alfaa)
	  kbdy=3
	  call bndyc(kbdy,xa,alfac,betac,gamac,gbdyc)
	  kbdy=1
	  call bndyc(kbdy,yc,alfaa,betaa,gamaa,gbdya)
	  det=alfaa*betac-betaa*alfac
	  cf(1,1,10)=cf(1,1,10)-c5*(dlx2*(betaa*gbdyc-gbdya*betac)/det)-
     +    c7*(dly2*(alfac*gbdya-gbdyc*alfaa)/det)
	end if
	if (nxb.eq.2) then
c     correct for mixed-mixed at (xb,yc)
	  kbdy=3
	  call bndyc(kbdy,xb-dlx,alfac,betac,gamac,gbdyc)
	  kbdy=2
	  call bndyc(kbdy,yc+dly,alfab,betab,gamab,gbdyb)
	  call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c2=cxy/dlxy2
	  c6=c2
	  c1=cxx/dlxx+cx/dlx2-c2
	  c7=cyy/dlyy-cy/dly2-c2
	  c1=c1+c6*(-alfac/betac*dyox)
	  c7=c7-c2*(betab/alfab*dxoy)
	  cf(nx,1,10)=cf(nx,1,10)+c6*(dly2*gbdyc/betac)-
     +    c2*(dlx2*gbdyb/alfab)
c     phase 2
	  kbdy=3
	  call bndyc(kbdy,xb,alfac,betac,gamac,gbdyc)
	  kbdy=2
	  call bndyc(kbdy,yc,alfab,betab,gamab,gbdyb)
	  det=betac*alfab-alfac*betab
	  cf(nx,1,10)=cf(nx,1,10)+c1*(dlx2*(betab*gbdyc-gbdyb*betac)/
     +    det)-c7*(dly2*(alfac*gbdyb-gbdyc*alfab)/det)
	end if
      end if

      if (nxb.eq.2) then
c     mixed along x=xb
	kbdy=2
	call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
	do j=2,ny-1
	  y=yc+j*dly
	  call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
	  call coef(xb,y-dly,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c1=cxx/dlxx+cx/dlx2
	  c2=cxy/dlxy4
	  c8=-c2
	  cf(nx,j,10)=cf(nx,j,10)-dlx2*(c8*gbdjm1/alfjm1+c1*gbdj/alfj+
     +    c2*gbdjp1/alfjp1)
	  alfjm1=alfj
	  gbdjm1=gbdj
	  alfj=alfjp1
	  gbdj=gbdjp1
	end do
      end if

      if (nyd.eq.2) then
	kbdy=4
	call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
	call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
	do i=2,nx-1
	  x=xa+i*dlx
	  call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
	  call coef(x-dlx,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c2=cxy/dlxy4
	  c3=cyy/dlyy+cy/dly2
	  c4=-c2
	  cf(i,ny,10)=cf(i,ny,10)-dly2*(c4*gbdim1/betim1+c3*gbdi/beti+
     +    c2*gbdip1/betip1)
	  betim1=beti
	  gbdim1=gbdi
	  beti=betip1
	  gbdi=gbdip1
	end do
	if (nxa.eq.0) then
c     correct for periodic-mixed at (xa,yd), (xb,yd)
	  kbdy=4
	  call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
	  call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
	  call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c2=cxy/dlxy2
	  c3=cyy/dlyy+cy/dly2-c2
	  cf(1,ny,10)=cf(1,ny,10)-dly2*(c3*gbdi/beti+c2*gbdip1/betip1)
	  call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
	  call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
	  call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c4=-cxy/dlxy2
	  c3=cyy/dlyy+cy/dly2-c4
	  cf(nx,ny,10)=cf(nx,ny,10)-dly2*(c3*gbdi/beti+c4*gbdim1/betim1)
	else if (nxa.eq.2) then
c     correct for mixed-mixed at (xa,yd)
c     phase 1
	  call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c2=cxy/dlxy2
	  c6=c2
	  c3=cyy/dlyy+cy/dly2-c2
	  c5=cxx/dlxx-cx/dlx2-c2
	  kbdy=4
	  call bndyc(kbdy,xa+dlx,alfad,betad,gamad,gbdyd)
	  kbdy=1
	  call bndyc(kbdy,yd-dly,alfaa,betaa,gamaa,gbdya)
	  c3=c3+c6*(-betaa/alfaa*dxoy)
	  c5=c5+c2*(-alfad/betad*dyox)
	  cf(1,ny,10)=cf(1,ny,10)+c2*(-dly2*gbdyd/betad)+
     +    c6*(dlx2*gbdya/alfaa)
c     phase 2
	  kbdy=4
	  call bndyc(kbdy,xa,alfad,betad,gamad,gbdyd)
	  kbdy=1
	  call bndyc(kbdy,yd,alfaa,betaa,gamaa,gbdya)
	  det=alfad*betaa-betad*alfaa
	  cf(1,ny,10)=cf(1,ny,10)-c5*(dlx2*(betad*gbdya-gbdyd*betaa)/
     +    det)-c3*(dly2*(alfad*gbdya-gbdyd*alfaa)/det)
	end if
	if (nxb.eq.2) then
c     correct ofr mixed-mixed at (xb,yd)
	  call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c4=-cxy/dlxy2
	  c8=c4
	  c1=cxx/dlxx+cx/dlx2-c4
	  c3=cyy/dlyy+cy/dly2-c4
	  kbdy=4
	  call bndyc(kbdy,xb-dlx,alfad,betad,gamad,gbdyd)
	  kbdy=2
	  call bndyc(kbdy, yd-dly,alfab,betab,gamab,gbdyb)
	  c1=c1+c4*(alfad/betad*dyox)
	  c3=c3+c8*(betab/alfab*dxoy)
	  cf(nx,ny,10)=cf(nx,ny,10)+c4*(-dly2*gbdyd/betad)+
     +    c8*(-dlx2*gbdyb/alfab)
c     phase 2
	  kbdy=4
	  call bndyc(kbdy,xb,alfad,betad,gamad,gbdyd)
	  kbdy=2
	  call bndyc(kbdy,yd,alfab,betab,gamab,gbdyb)
	  det=alfad*betab-betad*alfab
	  cf(nx,ny,10)=cf(nx,ny,10)-c1*(dlx2*(betab*gbdyd-gbdyb*betad)
     +    /det)-c3*(dly2*(alfad*gbdyb-gbdyd*alfab)/det)
	end if
      end if

      if (nxa.eq.2) then
c     mixed along x=xa
	kbdy=1
	call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
	do j=2,ny-1
	  y=yc+j*dly
	  call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
	  call coef(xa,y-dly,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	  c4=-cxy/dlxy4
	  c5=cxx/dlxx-cx/dlx2
	  c6=-c4
	  cf(1,j,10)=cf(1,j,10)+dlx2*(c6*gbdjm1/alfjm1+c5*gbdj/alfj+
     +    c4*gbdjp1/alfjp1)
	  alfjm1=alfj
	  gbdjm1=gbdj
	  alfj=alfjp1
	  gbdj=gbdjp1
	end do
      end if

      if (nxa.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xa,yc)
c
	kbdy = 1
	call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
	call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	c4 = -cxy/dlxy2
	c5 = cxx/dlxx-cx/dlx2-c4
	cf(1,1,10) = cf(1,1,10)+dlx2*(c4*gbdjp1/alfjp1+c5*gbdj/alfj)
c
c     mixed-periodic at (xa,yd)
c
	call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
	call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	c6=cxy/dlxy2
	c5 = cxx/dlxx-cx/dlx2-c6
	cf(1,ny,10)=cf(1,ny,10)+dlx2*(c6*gbdjm1/alfjm1+c5*gbdj/alfj)
      end if

      if (nxb.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xb,yc) and (xb,yd)
c
	kbdy = 2
	call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
	call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	c2 = cxy/dlxy2
	c1 = cxx/dlxx+cx/dlx2-c2
	cf(nx,1,10) = cf(nx,1,10)-dlx2*(c1*gbdj/alfj+c2*gbdjp1/alfjp1)
	call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
	call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
	  if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	    cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	  end if
	  if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	    cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	  end if
	c8 = -cxy/dlxy2
	c1 = cxx/dlxx+cx/dlx2-c8
	cf(nx,ny,10) = cf(nx,ny,10)-dlx2*(c1*gbdj/alfj+c8*gbdjm1/alfjm1)
      end if
c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      return
      end

      subroutine rescd2cr(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
        use omp_lib
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      complex rhsc(ncx,ncy),resf(nx,ny)
      complex phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      complex cof(nx,ny,10)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = (0.0,0.0)
	end do
      end do
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
	do i=1,nx
	  resf(i,j) = cof(i,j,10)-(
     +                cof(i,j,1)*phi(i+1,j)+
     +                cof(i,j,2)*phi(i+1,j+1)+
     +                cof(i,j,3)*phi(i,j+1)+
     +                cof(i,j,4)*phi(i-1,j+1)+
     +                cof(i,j,5)*phi(i-1,j)+
     +                cof(i,j,6)*phi(i-1,j-1)+
     +                cof(i,j,7)*phi(i,j-1)+
     +                cof(i,j,8)*phi(i+1,j-1)+
     +                cof(i,j,9)*phi(i,j))
	end do
      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call cres2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

      subroutine relcd2cr(nx,ny,phi,cof,tx,ty,sum)
c
c     relaxation for cud2
c
      implicit none
      integer nx,ny
      complex phi(*),cof(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
	call relcd2crp(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
	call slxcd2cr(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slycd2cr(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxcd2cr(nx,ny,phi,cof,tx,sum)
	call slycd2cr(nx,ny,phi,cof,ty,sum)
      end if
      return
      end

      subroutine relcd2crp(nx,ny,phi,cof)
c
c     gauss-seidel four color point relaxation
c
      use omp_lib
      implicit none
      integer nx,ny,i,j,lcolor,i1,i2,i3,i4,it
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      i1 = 1
      i2 = 4
      i3 = 3
      i4 = 2
c
c     sweep four colored grid points
c
      do lcolor=1,4
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(i1,cof,phi,nx,ny) PRIVATE(i,j)
	do j=1,ny,4
	  do i=i1,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(i2,cof,phi,nx,ny) PRIVATE(i,j)
	do j=2,ny,4
	  do i=i2,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(i3,cof,phi,nx,ny) PRIVATE(i,j)
	do j=3,ny,4
	  do i=i3,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(i4,cof,phi,nx,ny) PRIVATE(i,j)
	do j=4,ny,4
	  do i=i4,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
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
c    permute (i1,i2,i3,i4) for next color
c
	it = i4
	i4 = i3
	i3 = i2
	i2 = i1
	i1 = it
      end do
      return
      end

      subroutine slxcd2cr(nx,ny,phi,cof,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      use omp_lib
      implicit none
      integer nx,ny,i,ib,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),cof(nx,ny,10),tx(nx,ny,*),sum(ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call relcd2crp(nx,ny,phi,cof)
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
c     x direction not periodic, sweep odd j lines
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
	do j=1,ny,2
	  do i=1,nx
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
	do j=2,ny,2
	  do i=1,nx
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
	  sum(j) = (0.0,0.0)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c      sweep odd lines forward and back
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=1,ny,2
	  do i=1,nx-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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

      subroutine slycd2cr(nx,ny,phi,cof,ty,sum)
      implicit none
      integer nx,ny,i,j,jb
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),cof(nx,ny,10),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call relcd2crp(nx,ny,phi,cof)
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
	do i=1,nx,2
	  do j=1,ny
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
	  sum(i) = (0.0,0.0)
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
c     forward sweep odd x lines
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=1,nx,2
	  do j=1,ny-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
!$OMP PARALLEL DO SCHEDULE(STATIC,8) SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
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
