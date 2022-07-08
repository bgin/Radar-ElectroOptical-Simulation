c
c     file muh2.f
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
c     muh2 attempts to produce a second order finite difference
c     approximation to the two dimensional nonseparable elliptic
c     partial differential equation of the form:
c
c       cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +
c
c       ce(x,y)*pe(x,y) = r(x,y)
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mudcom.f
c
      subroutine muh2(iparm,fparm,wk,iwk,coef,bndyc,rhs,phi,mgopt,
     +                ierror)
          !dir$ optimize:3 
          !dir$ attributes code_align : 32 :: muh2
      implicit none
      integer iparm(17),mgopt(4),iwk(*),ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm(6),xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,k,kb,nx,ny,ic,itx,ity,iw
      real wk(*),phi(*),rhs(*)
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      !dir$ attributes align : 64 :: /fmud2/
      common/mud2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      !dir$ attributes align : 64 :: /mud2c/
      integer ibeta,ialfa,izmat,idmat
      common/muh2c/ibeta,ialfa,izmat,idmat
      !dir$ attributes align : 64 :: /muh2c/
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
c       check periodic coarse grid(s)
	if (ixp.lt.3 .and. nxa.eq.0) return
	if (jyq.lt.3 .and. nyc.eq.0) return
	ierror = 4
	ngrid = max0(iex,jey)
	if (iex.lt.1) return
	if (jey.lt.1) return
	ierror = 13
	if (iex.eq.1 .and. jey.eq.1) then
	  if (maxcy.gt.1) return
	  if (iguess.eq.1) return
	end if
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
c
c     set subgrid sizes
c
	do k=1,ngrid
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kps = kps+(nx+2)*(ny+2)+nx*ny*(6+isx+jsy)
	end do
c
c     set pointers for direct at coarse grid
c
	nx = ixp+1
	ny = jyq+1
	ibeta = kps+1
	if (nyc .eq. 0) then
	  ialfa = ibeta + nx*nx*(ny-1)
	  izmat = ialfa+nx*nx*(ny-1)
	  idmat = izmat+nx*nx*(ny-2)
	  kps = idmat+nx*nx*(ny-2)
	else
	  ialfa = ibeta + nx*nx*ny
	  kps = ialfa+nx*nx*ny
	end if
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
c     set pointers and discretize pde at each grid level
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
	  ic = kcbgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call dismh2(nx,ny,wk(ic),wk(itx),wk(ity),
     +                bndyc,coef,wk,iwk,ierror)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call muh21(nx,ny,rhs,phi,coef,bndyc,wk,iwk)
c
c     return if singular pde detected
c
      if (ierror .eq. 14) return
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end

      subroutine muh21(nx,ny,rhsf,phif,coef,bndyc,wk,iwk)
          !dir$ optimize:3 
          !dir$ attributes code_align : 32 :: muh21
      implicit none
      integer nx,ny,iwk(*)
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
     +     kcycle,iprer,ipost,intpol,kps
     !dir$ attributes align : 64 :: /imud2/
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
     !dir$ attributes align : 64 :: /fmud2/
      common/mud2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +     nxk(50),nyk(50),isx,jsy
     !dir$ attributes align : 64 :: /mud2c/
      integer ibeta,ialfa,izmat,idmat
      common/muh2c/ibeta,ialfa,izmat,idmat
      external coef,bndyc
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
	  call adjmh2(nx,ny,wk(ip),wk(ic),bndyc,coef)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcymh2(wk,iwk)
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
	call adjmh2(nx,ny,wk(ip),wk(ic),bndyc,coef)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcymh2(wk,iwk)
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

      subroutine kcymh2(wk,iwk)
          !dir$ optimize:3 
          !dir$ attributes code_align : 32 :: kcymh2
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      real wk(*)
      integer iwk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +     kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      !dir$ attributes align : 64 :: /fmud2/
      common/mud2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +     nxk(50),nyk(50),isx,jsy
      !dir$ attributes align : 64 :: /mud2c/
      integer ibeta,ialfa,izmat,idmat
      common/muh2c/ibeta,ialfa,izmat,idmat
      !dir$ attributes align : 64 :: /muh2c/
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      if (kcur .eq. 1) then
c
c     solve at coarse level with direct method and return
c
	if (nyc .ne. 0) then
	  call dir2(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),iwk)
	  return
	else
	  call dir2p(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),wk(izmat),
     +               wk(idmat),iwk)
	  return
	end if
      end if
c
c     prerelax at current finest grid level > 1
c
      do l=1,iprer
	call relmh2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+5*ncx*ncy
      call resmh2(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
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
c   relax or solve directly at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      if (klevel.gt.1) then
	do l=1,nrel
	  call relmh2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
	end do
      else
c
c     use direct method at coarsest level
c
	if (nyc .ne. 0) then
	  call dir2(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),iwk)
	else
	  call dir2p(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),wk(izmat),
     +               wk(idmat),iwk)
	end if
c
c     insure direct method is not called again at coarse level
c
	kount(1) = kcycle+1
      end if
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
	  call resmh2(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     +                wk(kps))
c
c     prerelax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 10
	else
c
c    direct  at coarsest level takes place of postrelax
c
	  ip = kpbgn(1)
	  ic = kcbgn(1)
	  nx = nxk(1)
	  ny = nyk(1)
	  if (nyc .ne. 0) then
	    call dir2(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),iwk)
	  else
	    call dir2p(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),wk(izmat)
     +                 ,wk(idmat),iwk)
	  end if
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
	call relmh2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end

      subroutine dismh2(nx,ny,cof,tx,ty,bndyc,coef,wk,iwk,ier)
          !dir$ optimize:3 
          !dir$ attributes code_align : 32 :: dismh2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: dismh2
          !dir$ attributes optimization_parameter: 'g2s=on' :: dismh2
c
c     discretize elliptic pde for muh2, set nonfatal errors
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer iwk(*),nx,ny,ist,ifn,jst,jfn,i,j,kbdy,l,im1,jm1,ier
      real cof(nx,ny,6),tx(nx,ny,*),ty(ny,nx,*)
      real wk(*),dlx,dlx2,dlxx,dly,dly2,dlyy,cmin,alfmax,cemax
      real x,y,cxx,cyy,cx,cy,ce,c1,c2,c3,c4,c5,alfa,gbdy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +     kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      !dir$ attributes align : 64 :: /fmud2/
      integer ibeta,ialfa,izmat,idmat
      common/muh2c/ibeta,ialfa,izmat,idmat
      !dir$ attributes align : 64 :: /muh2c/
      external bndyc,coef
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      cmin = 1.0
      alfmax = 0.0
      cemax = 0.0
c
c     set x,y subscript limits for calls to coef,bndyc
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
      !dir$ assume_aligned cof:64
      do j=jst,jfn
         y = yc+(j-1)*dly
         !dir$ vector aligned
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call coef(x,y,cxx,cyy,cx,cy,ce)
	  cmin = amin1(cmin,cxx*cyy)
	  cemax = amax1(abs(ce),cemax)
c
c     flag hyperbolic pde
c
	  if (klevel.eq.ngrid) then
	    if (abs(cx)*dlx.gt.2*abs(cxx) .or.
     +          abs(cy)*dly.gt.2*abs(cyy)) then
	      ier = -4
	    end if
	  end if
	  cxx = amax1(cxx,abs(cx)*dlx*0.5)
	  cyy = amax1(cyy,abs(cy)*dly*0.5)
	  c1 = cxx/dlxx-cx/dlx2
	  c2 = cxx/dlxx+cx/dlx2
	  c3 = cyy/dlyy-cy/dly2
	  c4 = cyy/dlyy+cy/dly2
	  c5 = ce-(c1+c2+c3+c4)
	  cof(i,j,1) = c1
	  cof(i,j,2) = c2
	  cof(i,j,3) = c3
	  cof(i,j,4) = c4
	  cof(i,j,5) = c5
	end do
      end do
c
c     adjust discretization for mixed derivative b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
        !dir$ assume_aligned cof:64
        !dir$ vector aligned
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
         !dir$ assume_aligned cof:64
        !dir$ vector aligned
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
         !dir$ assume_aligned cof:64
        !dir$ vector aligned
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
         !dir$ assume_aligned cof:64
        !dir$ vector aligned
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
	if (cemax.eq.0.0.and.alfmax.eq.0.0) then
	  if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	    if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
c
c     flag singular pde fatal error
c
	      ier = 14
	      return
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
c     set coefficient for specified boundaries
c
      if (nxa.eq.1) then
         i = 1
          !dir$ assume_aligned cof:64
        !dir$ vector aligned
	do j=1,ny
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (nxb.eq.1) then
         i = nx
          !dir$ assume_aligned cof:64
        !dir$ vector aligned
	do j=1,ny
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (nyc.eq.1) then
         j = 1
          !dir$ assume_aligned cof:64
        !dir$ vector aligned
	do i=1,nx
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (nyd.eq.1) then
         j = ny
          !dir$ assume_aligned cof:64
        !dir$ vector aligned
	do i=1,nx
	  do l=1,5
	    cof(i,j,l) = 0.0
	  end do
	  cof(i,j,5) = 1.0
	end do
      end if
      if (klevel .eq. 1) then
c
c     set block tri-diagonal coefficient matrix and do lu decomposition
c     for direct method at coarsest grid level
c
	nx = ixp+1
	ny = jyq+1
	if (nyc .ne. 0) then
c     factor non-periodic block matrix
	  call lud2(nx,ny,cof,wk(ibeta),wk(ialfa),iwk,nxa)
	  return
	else
c     factor periodic block matrix

	  do j =1,ny-1
	    call setbeta(nx,ny,cof,wk(ibeta),j,nxa)
	    call setalfa(nx,ny,cof,wk(ialfa),j)
	  end do
	  call lud2p(nx,ny,cof,wk(ibeta),wk(ialfa),wk(izmat),wk(idmat),
     +               iwk,nxa)
	  return
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
            !dir$ assume_aligned cof:64,tx:64
            !dir$ vector aligned
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
        !dir$ assume_aligned cof:64,tx:64
        !dir$ vector aligned
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
        !dir$ assume_aligned cof:64,ty:64
        !dir$ vector aligned
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
              !dir$ assume_aligned cof:64,ty:64
              !dir$ vector aligned
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

      subroutine lud2(nx,ny,cof,beta,alfa,index,nxa)
          !dir$ optimize:3 
          !dir$ attributes code_align : 32 :: lud2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: lud2
          !dir$ attributes optimization_parameter: 'g2s=on' :: lud2
c
c     decompose nonperiodic block coefficient matrix
c
      use omp_lib
      implicit none
      integer nx,ny,nxa
      real cof(nx,ny,6),beta(nx,nx,*),alfa(nx,nx,*)
      integer index(nx,ny)
      integer iz,i1,jcur,jm1,l,i
      i1 = 1
c
c     set and factor umat(1) in beta(1)
c
      jcur = 1
      call setbeta(nx,ny,cof,beta,jcur,nxa)
      call sgfa(beta,nx,nx,index,iz)
      do jcur=2,ny
	call setalfa(nx,ny,cof,alfa,jcur)
	call transp(nx,alfa(1,1,jcur))
	jm1 = jcur-1
	do l=1,nx
	  call sgsl(beta(1,1,jm1),nx,nx,index(1,jm1),alfa(1,l,jcur),i1)
	end do
	call transp(nx,alfa(1,1,jcur))
	call setbeta(nx,ny,cof,beta,jcur,nxa)
        !dir$ assume_aligned beta:64,alfa:64,cof:64
	do i=1,nx
        !dir$ vector aligned
        !dir$ ivdep
        !dir$ vector vectorlength(4)
        !$omp simd reduction(-:beta)   
	  do l=1,nx
	    beta(i,l,jcur)=beta(i,l,jcur)-alfa(i,l,jcur)*cof(l,jcur-1,4)
	  end do
	end do
c
c     factor current beta for next pass
c
	call sgfa(beta(1,1,jcur),nx,nx,index(1,jcur),iz)
      end do
      return
      end

      subroutine dir2(nx,ny,phi,cof,beta,alfa,index)
        !dir$ attributes forceinline :: dir2
c
c     direct solve at coarsest grid
c
      implicit none
      integer nx,ny,index(nx,ny)
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6)
      real beta(nx,nx,*),alfa(nx,nx,*)
c     forward sweep
      call for2(nx,ny,phi,cof(1,1,6),alfa)
c     backward sweep
      call bkw2(nx,ny,phi,cof,beta,index)
      return
      end

      subroutine for2(nx,ny,phi,frhs,alfa)
          !dir$ optimize:3
          !dir$ attributes forceinline :: for2
          !dir$ attributes code_align : 32 :: for2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: for2
          !dir$ attributes optimization_parameter: 'g2s=on' :: for2
c
c     forward sweep
c
      use omp_lib
      implicit none
      integer nx,ny,i,j,l
      real phi(0:nx+1,0:ny+1),frhs(nx,ny),alfa(nx,nx,*),sum
      !dir$ assume_aligned phi:64,frhs:64
      do j=1,ny
         !dir$ vector aligned
         !dir$ vector vectorlength(4)
         !dir$ unroll(16)
	do i=1,nx
	  phi(i,j)=frhs(i,j)
	end do
      end do
      !dir$ assume_aligned alpha:64,phi:64
      do j=2,ny
	do i=1,nx
           sum=0.0
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector vectorlength(4)
           !$omp simd reduction(+:sum)
	  do l=1,nx
	    sum=sum+alfa(i,l,j)*phi(l,j-1)
	  end do
	  phi(i,j)=phi(i,j)-sum
	end do
      end do
      return                                                                    
      end                                                                       

      
      subroutine bkw2(nx,ny,phi,cof,beta,index)
          !dir$ optimize:3
          !dir$ attributes forceinline :: bkw2
          !dir$ attributes code_align : 32 :: bkw2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: bkw2
          !dir$ attributes optimization_parameter: 'g2s=on' :: bkw2
c
      use omp_lib
      implicit none
      integer nx,ny,index(nx,ny)
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),beta(nx,nx,*)
      integer iz,j,jb,i
      iz = 0
      call sgsl(beta(1,1,ny),nx,nx,index(1,ny),phi(1,ny),iz)
      !dir$ assume_aligned phi:64,cof:64
      do jb=2,ny
         j = ny-jb+1
         !dir$ ivdep
         !dir$ vector aligned
         !dir$ vector vectorlength(4)
         !$omp simd reduction(-:phi)
	do i=1,nx
	  phi(i,j) = phi(i,j) - cof(i,j,4)*phi(i,j+1)
	end do
	call sgsl(beta(1,1,j),nx ,nx ,index(1,j),phi(1,j),iz)
      end do
      return                                                                    
      end                                                                       

      subroutine lud2p(nx,ny,cof,beta,alfa,zmat,dmat,index,nxa)
          !dir$ optimize:3
          !dir$ attributes code_align : 32 :: lud2p
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: lud2p
          !dir$ attributes optimization_parameter: 'g2s=on' :: lud2p
      implicit none
      integer nx,ny,index(nx,ny),i,iz,j,l,jcur,jm1,i1,nxa
      real cof(nx,ny,6),alfa(nx,nx,*),beta(nx,nx,*)
      real dmat(nx,nx,*),zmat(nx,nx,*),sum
      jcur = 1
c
c     set dmat(1)=alfa(1)
c
      call setalfa(nx,ny,cof,alfa,jcur)
      !dir$ assume_aligned dmat:64,alfa:64
      do i=1,nx
         !dir$ ivdep
         !dir$ vector aligned
         !dir$ vector always
         !dir$ vector vectorlength(4)
	do l=1,nx
	  dmat(i,l,1) = alfa(i,l,jcur)
	end do
      end do
c
c     factor umat(1) in beta(1)
c
      call setbeta(nx,ny,cof,beta,jcur,nxa)
      call sgfa(beta(1,1,1),nx,nx,index(1,1),iz)
      do jcur=2,ny-2
c
c     solve transpose of lmat(jcur)umat(jcur-1)=alfa(jcur) in alfa(jcur)
c
	call setalfa(nx,ny,cof,alfa,jcur)
c
c     transpose alfa
c
	call transp(nx,alfa(1,1,jcur))
c
c     solve transpose equation
c
	jm1 = jcur-1
	i1 = 1
	do l=1,nx
	  call sgsl(beta(1,1,jm1),nx,nx,index(1,jm1),alfa(1,l,jcur),i1)
	end do
c
c     re-transpose solution in alfa
c
	call transp(nx,alfa(1,1,jcur))
c
c
	call setbeta(nx,ny,cof,beta,jcur,nxa)
        !dir$ assume_aligned beta:64,alfa:64,cof:64
	do i=1,nx
        !dir$ ivdep
        !dir$ vector aligned
         !dir$ vector vectorlength(4)
          !$omp simd reduction(-:beta)
	  do l=1,nx
	    beta(i,l,jcur)=beta(i,l,jcur)-alfa(i,l,jcur)*cof(l,jcur-1,4)
	  end do
	end do
c
c     factor current beta(1,1,jcur) for next pass
c
	call sgfa(beta(1,1,jcur),nx,nx,index(1,jcur),iz)
c
c     set dmat(jcur) = -alfa(jcur)*dmat(jcur-1)
c
	do i=1,nx
	  do j=1,nx
	    dmat(i,j,jcur) = 0.0
	    do l=1,nx
	      dmat(i,j,jcur) = dmat(i,j,jcur)-alfa(i,l,jcur)*
     +                         dmat(l,j,jcur-1)
	    end do
	  end do
	end do
	if (jcur .eq. ny-2) then
c
c     adjust dmat(ny-2) = gama(ny-2)-alfa(ny-2)*dmat(ny-2)
c
           !dir$ assume_aligned dmat:64,cof:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector vectorlength(4)
           !dir$ vector always
	  do i=1,nx
	    dmat(i,i,jcur) = cof(i,jcur,4)+dmat(i,i,jcur)
	  end do
	end if
      end do
c
c     final phase with periodic factorization
c
c     solve transpose of zmat(1) beta(1) = gama(ny-1)
c
      do i=1,nx
	do j=1,nx
	  zmat(j,i,1) = 0.0
	end do
	zmat(i,i,1) = cof(i,ny-1,4)
      end do
      do l=1,nx
	call sgsl(beta(1,1,1),nx,nx,index(1,1),zmat(1,l,1),i1)
      end do
      call transp(nx,zmat(1,1,1))
c
c     solve transpose of zmat(j) umat(j) = -zmat(j-1) gama(j-1)
c
      !dir$ assume_aligned zmat:64,cof:64
      do jcur = 2,ny-3
         do i=1,nx
         !dir$ vector aligned
         !dir$ vector vectorlength(4)
         !dir$ ivdep
         !$omp simd reduction(-:zmat)
	  do j=1,nx
	      zmat(j,i,jcur) = -zmat(i,j,jcur-1)*cof(j,jcur-1,4)
	  end do
	end do
	do l=1,nx
	 call sgsl(beta(1,1,jcur),nx,nx,index(1,jcur),zmat(1,l,jcur),i1)
	end do
	call transp(nx,zmat(1,1,jcur))
      end do
c
c     solve transpose of zmat(ny-2)umat(ny-2)=alfa(ny-1)-zmat(ny-3)gama(ny-3)
c
      call setalfa(nx,ny,cof,alfa,ny-1)
      do i=1,nx
	do j=1,nx
	  zmat(j,i,ny-2) = alfa(i,j,ny-1)-zmat(i,j,ny-3)*cof(j,ny-3,4)
	end do
      end do
      do l=1,nx
	call sgsl(beta(1,1,ny-2),nx,nx,index(1,ny-2),zmat(1,l,ny-2),i1)
      end do
      call transp(nx,zmat(1,1,ny-2))
c
c     set umat(ny-1) = beta(ny-1)-(zmat(1)*dmat(1)+...+zmat(ny-2)*dmat(ny-2))
c     in beta(ny-1)
c
      call setbeta(nx,ny,cof,beta,ny-1,nxa)
      do  i=1,nx
	do j=1,nx
	  sum = 0.0
	  do jcur=1,ny-2
	    do l=1,nx
	      sum = sum + zmat(i,l,jcur)*dmat(l,j,jcur)
	    end do
	  end do
	  beta(i,j,ny-1) = beta(i,j,ny-1) - sum
	end do
      end do
c
c     factor bmat(ny-1) for backward sweep
c
      call sgfa(beta(1,1,ny-1),nx,nx,index(1,ny-1),iz)
c
c     lud is now complete
c
      return
      end

      subroutine dir2p(nx,ny,phi,cof,beta,alfa,zmat,dmat,index)
      !dir$ optimize:3
      !dir$ attributes forceinline :: dir2p
         
c
c     direct method for periodic b.c.
c
      use omp_lib
      implicit none
      integer nx,ny,index(nx,ny)
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6)
      real beta(nx,nx,*),alfa(nx,nx,*)
      real zmat(nx,nx,*), dmat(nx,nx,*)
c     forward sweep
      call for2p(nx,ny,phi,cof(1,1,6),alfa,zmat)
c     backward sweep
      call bkw2p(nx,ny,phi,cof,beta,dmat,index)
      return
      end

      subroutine for2p(nx,ny,phi,frhs,alfa,zmat)
          !dir$ optimize:3
          !dir$ attributes forceinline :: for2p
          !dir$ attributes code_align : 32 :: for2p
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: for2p
          !dir$ attributes optimization_parameter: 'g2s=on' :: for2p
       use omp_lib
      implicit none
      integer nx,ny,i,j,jcur,l,k
      real phi(0:nx+1,0:ny+1),frhs(nx,ny),sum
      real alfa(nx,nx,*),zmat(nx,nx,*)
      !dir$ assume_aligned phi:64,frhs:64
      do j=1,ny-1
       !dir$ vector aligned
       !dir$ vector vectorlength(4)
       !dir$ unroll(16)
 	do i=1,nx
	  phi(i,j)=frhs(i,j)
	end do
      end do
      do jcur=2,ny-2
	do i=1,nx
           sum=0.0
          !dir$ vector aligned
          !dir$ vector vectorlength(4)
          !dir$ ivdep
          !$omp simd reduction(+:sum)
	  do l=1,nx
	    sum=sum+alfa(i,l,jcur)*phi(l,jcur-1)
	  end do
	  phi(i,jcur)=phi(i,jcur)-sum
	end do
      end do
c
c     solve:
c     zmat(1)*phi(1)+...+zmat(ny-2)*phi(ny-2) + phi(ny-1) = f(ny-1)
c
      !dir$ assume_aligned zmat:64,phi:64
      do i=1,nx
	sum = 0.0
	do k=1,ny-2
             !dir$ vector aligned
          !dir$ vector vectorlength(4)
          !dir$ ivdep
          !$omp simd reduction(+:sum)
	  do l=1,nx
	    sum = sum + zmat(i,l,k)*phi(l,k)
	  end do
	end do
	phi(i,ny-1) = phi(i,ny-1) - sum
      end do
      return
      end

      subroutine bkw2p(nx,ny,phi,cof,beta,dmat,index)
             !dir$ optimize:3
          !dir$ attributes forceinline :: bkw2p
          !dir$ attributes code_align : 32 :: bkw2p
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: bkw2p
          !dir$ attributes optimization_parameter: 'g2s=on' :: bkw2p
c
          use omp_lib 
      implicit none
      integer nx,ny,index(nx,ny),iz,i,l,kb,k
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),beta(nx,nx,ny),dmat(nx,nx,*)
      real sum
      iz = 0
      call sgsl(beta(1,1,ny-1),nx,nx,index(1,ny-1),phi(1,ny-1),iz)
c
c     solve beta(ny-2)*phi(ny-2) = phi(ny-2)-dmat(ny-2)*phi(ny-1)
c
      !dir$ assume_aligned dmat:64,phi:64
      do i=1,nx
         sum = 0.0
            !dir$ vector aligned
          !dir$ vector vectorlength(4)
          !dir$ ivdep
          !$omp simd reduction(+:sum)
         do l=1,nx
	  sum = sum + dmat(i,l,ny-2)*phi(l,ny-1)
	end do
	phi(i,ny-2) = phi(i,ny-2) - sum
      end do
      call sgsl(beta(1,1,ny-2),nx,nx,index(1,ny-2),phi(1,ny-2),iz)
c
c     solve beta(k)*phi(k) = phi(k) - gama(k)*phi(k+1)-dmat(k)*phi(ny-1)
c     k=ny-3,...,1
c
      do kb=4,ny
	k = ny-kb+1
	do i=1,nx
           sum = 0.0
              !dir$ vector aligned
          !dir$ vector vectorlength(4)
          !dir$ ivdep
          !$omp simd reduction(+:sum)
	  do l=1,nx
	    sum = sum+dmat(i,l,k)*phi(l,ny-1)
	  end do
	  phi(i,k) = phi(i,k) - sum - cof(i,k,4)*phi(i,k+1)
	end do
	call sgsl(beta(1,1,k),nx,nx,index(1,k),phi(1,k),iz)
      end do
c
c     set j=ny by periodicity
c
      do i=1,nx
	phi(i,ny) = phi(i,1)
      end do
      return
      end

      subroutine setbeta(nx,ny,cof,beta,jcur,nxa)
             !dir$ optimize:3
          !dir$ attributes forceinline :: setbeta
          !dir$ attributes code_align : 32 :: setbeta
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: setbeta
          !dir$ attributes optimization_parameter: 'g2s=on' :: setbeta
c
c     set diagonal matrix on block
c
      implicit none
      integer nx,ny,jcur,nxa,i,l,j
      real cof(nx,ny,6),beta(nx,nx,*)
      j = jcur
      do i=1,nx
	do l=1,nx
	  beta(i,l,jcur)=0.0
	end do
      end do
      do i=1,nx
	beta(i,i,jcur) = cof(i,j,5)
      end do
      do i=2,nx
	beta(i,i-1,jcur) = cof(i,j,1)
      end do
      do i=1,nx-1
	beta(i,i+1,jcur) = cof(i,j,2)
      end do
      if (nxa.eq.0) then                                                        
	beta(1,nx-1,jcur) = cof(1,j,1)
	beta(nx,2,jcur) = cof(nx,j,2)
      end if                                                                    
      return                                                                    
      end                                                                       

      subroutine setalfa(nx,ny,cof,alfa,jcur)
      implicit none
      integer nx,ny,jcur,i,l
      real cof(nx,ny,6),alfa(nx,nx,*)
      do i=1,nx
	do l=1,nx
	  alfa(i,l,jcur)=0.0
	end do
	alfa(i,i,jcur)=cof(i,jcur,3)
      end do
      return
      end

      subroutine adjmh2(nx,ny,phi,cof,bndyc,coef)
                !dir$ optimize:3
          !dir$ attributes forceinline :: adjmh2
          !dir$ attributes code_align : 32 :: adjmh2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: adjmh2
          !dir$ attributes optimization_parameter: 'g2s=on' :: adjmh2
c
c     adjust righthand side in cof(i,j,6) for boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6)
      real dlx,dlx2,dlxx,dly,dly2,dlyy
      real x,y,cxx,cyy,cx,cy,ce,c1,c2,c3,c4,alfa,gbdy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +     kcycle,iprer,ipost,intpol,kps
     !dir$ attributes align : 64 :: /imud2/
     
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      !dir$ attributes align : 64 :: /fmud2/
      external bndyc,coef
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlx2 = dlx+dlx
      dly2 = dly+dly
      dlxx = dlx*dlx
      dlyy = dly*dly
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
	  call coef(x,y,cxx,cyy,cx,cy,ce)
	  cxx = amax1(cxx,abs(cx)*dlx*0.5)
	  c1 = cxx/dlxx-cx/dlx2
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
	  call coef(x,y,cxx,cyy,cx,cy,ce)
	  cxx = amax1(cxx,abs(cx)*dlx*0.5)
	  c2 = cxx/dlxx+cx/dlx2
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
	  call coef(x,y,cxx,cyy,cx,cy,ce)
	  cyy = amax1(cyy,abs(cy)*dly*0.5)
	  c3 = cyy/dlyy-cy/dly2
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
	  call coef(x,y,cxx,cyy,cx,cy,ce)
	  cyy = amax1(cyy,abs(cy)*dly*0.5)
	  c4 = cyy/dlyy+cy/dly2
	  cof(i,j,6) = cof(i,j,6)-dly2*c4*gbdy
	end do
      end if
c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
         i = 1
         !dir$ assume_aligned cof:64,phi:64
         !dir$ vector aligned
         !dir$ vector always
	do j=1,ny
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      if (nxb.eq.1) then
         i = nx
         !dir$ assume_aligned cof:64,phi:64
         !dir$ vector aligned
         !dir$ vector always
	do j=1,ny
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      if (nyc.eq.1) then
         j = 1
         !dir$ assume_aligned cof:64,phi:64
         !dir$ vector aligned
         !dir$ vector always
	do i=1,nx
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      if (nyd.eq.1) then
         j = ny
            !dir$ assume_aligned cof:64,phi:64
         !dir$ vector aligned
         !dir$ vector always
	do i=1,nx
	  cof(i,j,6) = phi(i,j)
	end do
      end if
      return
      end

      subroutine resmh2(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
                   !dir$ optimize:3
                   !dir$ attributes forceinline :: resmh2
                   !dir$ attributes code_align : 32 :: resmh2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: resmh2
          !dir$ attributes optimization_parameter: 'g2s=on' :: resmh2
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
     +     kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
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
        !dir$ assume_aligned resf:64,cof:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
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

      subroutine relmh2(nx,ny,phi,cof,tx,ty,sum)
              
                   !dir$ attributes forceinline :: relmh2
                   !dir$ attributes code_align : 32 :: relmh2
          
c
c     relaxation for muh2
c
      implicit none
      integer nx,ny
      real phi(*),cof(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +     kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
      if (method.eq.0) then                ! point relaxation
	call relmh2p(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
	call slxmh2(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slymh2(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxmh2(nx,ny,phi,cof,tx,sum)
	call slymh2(nx,ny,phi,cof,ty,sum)
      end if
      return
      end

      subroutine relmh2p(nx,ny,phi,cof)
          !dir$ optimize:3
          !dir$ attributes code_align : 32 :: relmh2p
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: relmh2p
          !dir$ attributes optimization_parameter: 'g2s=on' :: relmh2p
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
     +     kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
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
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
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
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=1,nx,2
	  do j=2,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
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
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=1,nx,2
	do j=1,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
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
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=1,nx,2
	do j=2,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
c
!$OMP PARALLEL DO DO SCHEDULE(STATIC,8) SHARED(cof,phi,nx,ny) PRIVATE(i,j)
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

      subroutine slxmh2(nx,ny,phi,cof,tx,sum)
          !dir$ optimize:3
          !dir$ attributes code_align : 32 :: slxmh2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: slxmh2
          !dir$ attributes optimization_parameter: 'g2s=on' :: slxmh2
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      use omp_lib
      implicit none
      integer nx,ny,i,ib,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +     kcycle,iprer,ipost,intpol,kps
     !dir$ attributes align : 64 :: /imud2/
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),tx(nx,ny,*),sum(ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call relmh2p(nx,ny,phi,cof)
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
             !dir$ assume_aligned cof:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
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
            !dir$ assume_aligned cof:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
	  do i=1,nx
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
           !dir$ assume_aligned tx:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
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
            !dir$ assume_aligned cof:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
	  do i=1,nx-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
         end do
         !dir$ assume_aligned sum:64,tx:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(+:sum)
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
            !dir$ assume_aligned cof:64,phi:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
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
        !dir$ assume_aligned sum:64,phi:64,tx:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(+:sum)
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

      subroutine slymh2(nx,ny,phi,cof,ty,sum)
        !dir$ optimize:3
          !dir$ attributes code_align : 32 :: slymh2
          !dir$ attributes optimization_parameter: 'target_arch=skylake-avx512' :: slymh2
          !dir$ attributes optimization_parameter: 'g2s=on' :: slymh2
      use omp_lib
      implicit none
      integer nx,ny,i,j,jb
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +     kcycle,iprer,ipost,intpol,kps
      !dir$ attributes align : 64 :: /imud2/
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call relmh2p(nx,ny,phi,cof)
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
         !dir$ assume_aligned phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(-:phi)
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
             !dir$ assume_aligned phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(-:phi)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     forward sweep even x lines
c
!$OMP PARALLEL DO SCHEDULE(STTAIC,8) SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
         end do
             !dir$ assume_aligned phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(-:phi)
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
              !dir$ assume_aligned phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(-:phi)
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
            !dir$ assume_aligned phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(-:phi)
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
         end do
             !dir$ assume_aligned sum:64,phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(+:sum)
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +         ty(ny-2,i,2)
           
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
             !dir$ assume_aligned phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(-:phi)
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
         end do
             !dir$ assume_aligned sum:64,phi:64,ty:64
        !dir$ ivdep
        !dir$ vector aligned
        !dir$ vector vectorlength(4)
        !dir$ vector always
          !$omp simd reduction(+:sum)
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +         ty(ny-2,i,2)
           
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
