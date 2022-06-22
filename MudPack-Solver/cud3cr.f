c
c     file cud3cr.f
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
c     cud3cr attempts to produce a second order finite difference
c     approximation to the complex 3-d nonseparable elliptic
c     partial differential equation with cross derivatives:
c
c       cxx(x,y,z)*pxx + cyy(x,y,z)*pyy + czz(x,y,z)*pzz +
c
c       cxy(x,y,z)*pxy + cxz(x,y,z)*pxz + cyz(x,y,z)*pzz +
c
c       cx(x,y,z)*px + cy(x,y,z)*py + cz(x,y,z)*pz +
c
c       ce(x,y,z)*p(x,y,z) = r(x,y,z)
c
c     Here the coefficients cxx,cyy,czz,cxy,cxz,cyz,cx,cy,cz,ce and
c     the right hand side r and the unknown solution function p are all
c     complex valued functions of the real independent variables x,y,z.
c
c ... see documentation and test files provided in this distribution
c
c ... required mudpack files
c
c     cudcom.f
c
      subroutine cud3cr(iparm,fparm,wk,coef,bnd3cr,rhs,phi,mgopt,
     +     icros,crsxy,crsxz,crsyz,tol,maxit,iouter,rmax,ierror)
           !dir$ attributes code_align : 32 :: cud3cr
           !dir$ optimize : 3
      implicit none
      integer iparm(23),mgopt(4),icros(3),maxit,iouter,ierror
      real fparm(8),tol,rmax(maxit)
      complex wk(*),phi(*),rhs(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +     kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/fcd3cr/xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,krbgn,kxybgn,kxzbgn,kyzbgn,ktxbgn,ktybgn,
     +     ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /fcd3cr/
      common/cd3cr/kpbgn(50),kcbgn(50),krbgn(50),kxybgn(50),kxzbgn(50),
     +kyzbgn(50),ktxbgn(50),ktybgn(50),ktzbgn(50),nxk(50),nyk(50),
     +     nzk(50),ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /cd3cr/
      integer int,m,isx,jsy,ksz,ixy,ixz,iyz
      integer nx,ny,nz,itx,ity,itz,k,kb,ip,ic,ir
      real  dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common / incr3 / dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      !dir$ attributes align : 64 :: /incr3/
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +     kxy,kxz,kyz
      !dir$ attributes align : 64 :: /kcrsxyz/
      external coef,bnd3cr,crsxy,crsxz,crsyz
      data int / 0 /
      save int
      ierror = 1
      intl = iparm(1)    ! set and check intl on ALL calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
	int = 1
	if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
c
c     set input parameters from iparm,fparm internally
c
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      nze = iparm(6)
      nzf = iparm(7)
c
c     set grid size params
c
      ixp = iparm(8)
      jyq = iparm(9)
      kzr = iparm(10)
      iex = iparm(11)
      jey = iparm(12)
      kez = iparm(13)
c
c     set number of subgrids for mg cycling
c
      ngrid = max0(iex,jey,kez)
      nfx = iparm(14)
      nfy = iparm(15)
      nfz = iparm(16)

      iguess = iparm(17)
      maxcy = iparm(18)
      method = iparm(19)
      meth2 = iparm(20)
      nwork = iparm(21)
c
c     set floating point params
c
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      ze = fparm(5)
      zf = fparm(6)
      tolmax = fparm(7)
c
c     set default multigrid option parameters
c
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
c
c     use default settings
c
	kcycle = 2
	iprer = 2
	ipost = 1
	intpol = 3
      else
c
c     flag fatal error
c
	ierror = 12
	return
      end if
      if (intl .eq. 0) then  ! intialization call
c
c     check input arguments
c
	ierror = 2   ! check boundary condition flags
	if (max0(nxa,nxb,nyc,nyd,nze,nzf).gt.2) return
	if (min0(nxa,nxb,nyc,nyd,nze,nzf).lt.0) return
	if (nxa.eq.0.and.nxb.ne.0) return
	if (nxa.ne.0.and.nxb.eq.0) return
	if (nyc.eq.0.and.nyd.ne.0) return
	if (nyc.ne.0.and.nyd.eq.0) return
	if (nze.eq.0.and.nzf.ne.0) return
	if (nze.ne.0.and.nzf.eq.0) return
	ierror = 3   ! check grid sizes
	if (ixp.lt.2) return
	if (jyq.lt.2) return
	if (kzr.lt.2) return
	ierror = 4
	ngrid = max0(iex,jey,kez)
	if (iex.lt.1) return
	if (jey.lt.1) return
	if (kez.lt.1) return
	if (ngrid.gt.50) return
	ierror = 5
	if (nfx.ne.ixp*2**(iex-1)+1) return
	if (nfy.ne.jyq*2**(jey-1)+1) return
	if (nfz.ne.kzr*2**(kez-1)+1) return
	ierror = 6
	if (iguess*(iguess-1).ne.0) return
	ierror = 7
	if (maxcy.lt.1 .or. maxcy.gt.2) return
	ierror = 8
	if (method.lt.0 .or. method.gt.7) return  ! verify method chosen
c
c     set and check cross derivative indicators
c
	kxy = icros(1)
	kxz = icros(2)
	kyz = icros(3)
	ierror = 15
	if (kxy*(kxy-1).ne.0) return
	if (kxz*(kxz-1).ne.0) return
	if (kyz*(kyz-1).ne.0) return
c
c     compute required work space length
c
	ierror = 9
	m = method
	isx = 0
	if ((m-1)*(m-4)*(m-5)*(m-7).eq.0) then
	  isx = 3
	  if (nxa.eq.0) then
	    isx = 5
	  end if
	end if
	jsy = 0
	if ((m-2)*(m-4)*(m-6)*(m-7).eq.0) then
	  jsy = 3
	  if (nyc.eq.0) then
	    jsy = 5
	  end if
	end if
	ksz = 0
	if ((m-3)*(m-5)*(m-6)*(m-7).eq.0) then
	  ksz = 3
	  if (nze.eq.0) then
	    ksz = 5
	  end if
	end if
c
c     set subgrid sizes
c
	do k=1,ngrid
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nzk(k) = kzr*2**(max0(k+kez-ngrid,1)-1)+1
	end do
	kps = 1
	do kb=1,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k)
	  ny = nyk(k)
	  nz = nzk(k)
c
c     set pointers
c
	  kpbgn(k) = kps
	  kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)*(nz+2)
	  krbgn(k) = kcbgn(k)+7*nx*ny*nz
	  kxybgn(k) = krbgn(k) + nx*ny*nz
	  kxzbgn(k) = kxybgn(k)+kxy*nx*ny*nz
	  kyzbgn(k) = kxzbgn(k)+kxz*nx*ny*nz
	  ktxbgn(k) = kyzbgn(k) + kyz*nx*ny*nz
	  ktybgn(k) = ktxbgn(k) + isx*nx*ny*nz
	  ktzbgn(k) = ktybgn(k) + jsy*nx*ny*nz
	  kps = ktzbgn(k) + ksz*nx*ny*nz
	end do
c
c     set pointers for saving cross derivative coefficients on
c     fine grid boundaries
c
	kxyxa = kps
	kxyxb = kxyxa+kxy*nfy*nfz
	kxyyc = kxyxb+kxy*nfy*nfz
	kxyyd = kxyyc+kxy*nfx*nfz
	kxyze = kxyyd+kxy*nfx*nfz
	kxyzf = kxyze+kxy*nfx*nfy
	kxzxa = kxyzf+kxy*nfx*nfy
	kxzxb = kxzxa+kxz*nfy*nfz
	kxzyc = kxzxb+kxz*nfy*nfz
	kxzyd = kxzyc+kxz*nfx*nfz
	kxzze = kxzyd+kxz*nfx*nfz
	kxzzf = kxzze+kxz*nfx*nfy
	kyzxa = kxzzf+kxz*nfx*nfy
	kyzxb = kyzxa+kyz*nfy*nfz
	kyzyc = kyzxb+kyz*nfy*nfz
	kyzyd = kyzyc+kyz*nfx*nfz
	kyzze = kyzyd+kyz*nfx*nfz
	kyzzf = kyzze+kyz*nfx*nfy
	kps  = kyzzf + kyz*nfx*nfy
c
c     set and check minimal work space
c
	iparm(22)=kps+(nfx+2)*(nfy+2)*(nfz+2)
	lwork = iparm(22)
	if (lwork .gt. nwork) return
	ierror = 10   ! check solution region
	if (xb.le.xa .or. yd.le.yc .or. zf.le.ze) return
	ierror = 11
	if (tolmax .lt. 0.0) return
	ierror = 13
	if (tol.le.0.0) return
	ierror = 14
	if (maxit.lt.1) return
	if (ierror .gt. 0) ierror = 0   ! no fatal errors
c
c     discretize pde at each grid level
c
	do k=1,ngrid
	  klevel = k
	  nx = nxk(k)
	  ny = nyk(k)
	  nz = nzk(k)
	  ip = kpbgn(k)
	  ic = kcbgn(k)
	  ixy = kxybgn(k)
	  ixz = kxzbgn(k)
	  iyz = kyzbgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  itz = ktzbgn(k)
	  call discd3cr(nx,ny,nz,wk(ic),wk(ixy),wk(ixz),wk(iyz),wk(itx),
     +    wk(ity),wk(itz),bnd3cr,coef,crsxy,crsxz,crsyz,wk,ierror)
	end do
c
c     set fine grid increments
c
	dx = (xb-xa)/(nfx-1)
	dy = (yd-yc)/(nfy-1)
	dz = (zf-ze)/(nfz-1)
	dx2 = dx+dx
	dy2 = dy+dy
	dz2 = dz+dz
	odxy4 = 1.0/(4.*dx*dy)
	odxz4 = 1.0/(4.*dx*dz)
	odyz4 = 1.0/(4.*dy*dz)
	if (max0(kxy,kxz,kyz) .eq. 1) then
c
c     set nonzero cross coefs on boundaries of fine grid
c     for outer iteration modification of right hand side
c
	  nx = nfx
	  ny = nfy
	  nz = nfz
	  call csexyzb(nx,ny,nz,crsxy,crsxz,crsyz,wk(kxyxa),wk(kxyxb),
     +    wk(kxyyc),wk(kxyyd),wk(kxyze),wk(kxyzf),wk(kxzxa),wk(kxzxb),
     +    wk(kxzyc),wk(kxzyd),wk(kxzze),wk(kxzzf),wk(kyzxa),wk(kyzxb),
     +    wk(kyzyc),wk(kyzyd),wk(kyzze),wk(kyzzf))
	end if
c
c       allow for nonfatal error detection in discretization
c
	if (ierror.gt.0) ierror = 0
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      nz = nfz
      ip = kpbgn(ngrid)
      ir = krbgn(ngrid)
      call cud3cr1(nx,ny,nz,wk,coef,bnd3cr,rhs,phi,wk(ir),wk(ip),
     +crsxy,crsxz,crsyz,tol,rmax,maxit,iouter,ierror)
      return
      end

      subroutine cud3cr1(nx,ny,nz,wk,coef,bnd3cr,rhsf,phif,rhs,phi,
     +     crsxy,crsxz,crsyz,tol,rmax,maxit,iouter,ierror)
        !dir$ attributes code_align : 32 :: cud3cr1
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: cud3cr1
      use omp_lib
      implicit none
      integer nx,ny,nz,maxit,iouter,ierror
      real tol,rmax(maxit)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax,difmax,pmax
      integer ip,ic,ir,icc,irc,ncx,ncy,ncz,ipc
      integer i,j,k,kb,iter
      integer ist,ifn,jst,jfn,kst,kfn
      complex phif(nx,ny,nz),rhsf(nx,ny,nz),wk(*)
      complex phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +     kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      common/fcd3cr/xa,xb,yc,yd,ze,zf,tolmax,relmax
      !dir$ attributes align : 64 :: /fcd3cr/
      integer kpbgn,kcbgn,krbgn,kxybgn,kxzbgn,kyzbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      common/cd3cr/kpbgn(50),kcbgn(50),krbgn(50),kxybgn(50),kxzbgn(50),
     +kyzbgn(50),ktxbgn(50),ktybgn(50),ktzbgn(50),nxk(50),nyk(50),
     +     nzk(50),ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /cd3cr/
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +     kxy,kxz,kyz
      !dir$ attributes align : 64 :: /kcrsxyz/
      external coef,bnd3cr,crsxy,crsxz,crsyz
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = krbgn(ngrid)
c
c     set phif,rhsf in phi,rhs
c
      !dir$ assume_aligned phif:64
      !dir$ assume_aligned phi:64
      !dir$ assume_aligned rhsf:64
      !dir$ assume_aligned rhs:64
      !dir$ code_align(32)
      !$omp parallel do default(none) schedule(static,8) private(k,j,j) &
      !$omp& shared(nfz,nfy,nfx,phif,phi,rhsf,rhs)
      do k=1,nfz
         do j=1,nfy
         !dir$ vector aligned
         !dir$ vector always
         !dir$ unroll(16)   
	  do i=1,nfx
	    phi(i,j,k) = phif(i,j,k)
	    rhs(i,j,k) = rhsf(i,j,k)
	  end do
	end do
      end do
      if (iguess.eq.0) then
c
c     no initial guess so ensure phif=0.0 on interior and
c     nonspecified boundaries
c
      ist = 1
      if (nxa.eq.1) ist = 2
      ifn = nfx
      if (nxb.eq.1) ifn = nfx-1
      jst = 1
      if (nyc.eq.1) jst = 2
      jfn = nfy
      if (nyd.eq.1) jfn = nfy-1
      kst = 1
      if (nze.eq.1) kst = 2
      kfn = nfz
      if (nzf.eq.1) kfn = nfz-1
      !dir$ assume_aligned phif:64
      !dir$ code_align(32)
      do k=kst,kfn
         do j=jst,jfn
         !dir$ vector nontemporal(phif)
         !dir$ vector aligned   
         !dir$ vector always
         !dir$ unroll(16)   
	  do i=ist,ifn
	    phif(i,j,k) = (0.,0.)
	  end do
	end do
      end do
c
c     pass phi,rhs down to all grid levels for full multigrid cycling
c
	do kb=2,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  nz = nzk(k+1)
	  ip = kpbgn(k+1)
	  ic = kcbgn(k+1)
	  ir = krbgn(k+1)
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ncz = nzk(k)
	  ipc = kpbgn(k)
	  icc = kcbgn(k)
	  irc = krbgn(k)
	  call ctrsfc3(nx,ny,nz,wk(ip),wk(ir),ncx,ncy,ncz,
     +                 wk(ipc),wk(irc))
c
c     adjust right hand side at k grid level (all but finest)
c
	  call adjcd3cr(ncx,ncy,ncz,wk(ipc),wk(irc),phif,bnd3cr,coef)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcycd3cr(wk)
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  nz = nzk(k+1)
	  ip = kpbgn(k+1)
	  ipc = kpbgn(k)
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ncz = nzk(k)
c
c     prolong approximation from k to k+1
c
	  call cprolon3(ncx,ncy,ncz,wk(ipc),nx,ny,nz,wk(ip),
     +                  nxa,nxb,nyc,nyd,nze,nzf,intpol)
	end do
c
c     set initial guess p(0) from phi in phif
c
        !dir$ assume_aligned phif:64
        !dir$ assume_aligned phi:64
 1      !dir$ code_align(32)
	do k=1,nfz
          do j=1,nfy
          !dir$ vector aligned
          !dir$ vector always
          !dir$ unroll(16)
	    do i=1,nfx
	      phif(i,j,k) = phi(i,j,k)
	    end do
	  end do
	end do
      end if
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
c
c     begin outer loop from finest grid level
c
      kcur = ngrid
      do iouter= 1,maxit
c
c     adjust rhs by (1) re-incorporating b.c. and (2) subtracting
c     cross derivative estimates from nonspecified boundaries
c
	call adjcd3cr(nx,ny,nz,phi,rhs,phif,bnd3cr,coef)
	if (kxy.eq.1) then
	  call csubxy(nx,ny,nz,phif,rhs,wk(kxyxa),wk(kxyxb),
     +               wk(kxyyc),wk(kxyyd),wk(kxyze),wk(kxyzf))
	end if
	if (kxz.eq.1) then
	  call csubxz(nx,ny,nz,phif,rhs,wk(kxzxa),wk(kxzxb),
     +               wk(kxzyc),wk(kxzyd),wk(kxzze),wk(kxzzf))
	end if
	if (kyz.eq.1) then
	  call csubyz(nx,ny,nz,phif,rhs,wk(kyzxa),wk(kyzxb),
     +               wk(kyzyc),wk(kyzyd),wk(kyzze),wk(kyzzf))
	end if
c
c      execute maxcy w cycles from finest grid level
c
	do iter=1,maxcy
	  call kcycd3cr(wk)
	end do
c
c     phi now contains kth iterate where k = iouter
c
	rmax(iouter) = 0.0
	if (tol.gt.0.0 .or. iouter.eq.maxit) then
c
c     compute maximum relative difference between outer and outer-1
c     iterates in rmax
c
	  pmax = 0.0
	  difmax = 0.0
          
	  do k=1,nfz
	    do j=1,nfy
	      do i=1,nfx
		pmax = amax1(cabs(phi(i,j,k)),pmax)
		difmax = amax1(difmax,cabs(phi(i,j,k)-phif(i,j,k)))
	      end do
	    end do
	  end do
	  if (pmax .gt. 0.0) then
	    rmax(iouter) = difmax/pmax
	  else
c
c     degenerate case
c
	    rmax(iter) = difmax
	  end if
	end if
c
c     update current estimate in phif and restore rhs from rhsf
c
       
        !dir$ assume_aligned phif:64
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned rhsf:64
        !dir$ assume_aligned rhs:64
        do k=1,nfz
	  do j=1,nfy
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
            !dir$ unroll(8)
	    do i=1,nfx
	      phif(i,j,k) = phi(i,j,k)
	      rhs(i,j,k) = rhsf(i,j,k)
	    end do
	  end do
	end do
c
c     check for convergence if flagged
c
	if (tol .gt. 0.0) then
	  if (rmax(iouter) .lt. tol) return
	  if (iouter .eq. maxit) then
c
c     flag convergence failture in maxit outer iterations
c
	    ierror = -10
	    return
	  end if
	end if
      end do
      return
      end

      subroutine kcycd3cr(wk)
       !dir$ attributes code_align : 32 :: kcycd3cr
       !dir$ optimize : 3
c
c     perform multigrid k-cycle at kcur level
c     kcycle = 1 corresponds to v cycles
c     kcycle = 2 corresponds to w cycles
c
      implicit none
      complex wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      integer kpbgn,kcbgn,krbgn,kxybgn,kxzbgn,kyzbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      common/cd3cr/kpbgn(50),kcbgn(50),krbgn(50),kxybgn(50),kxzbgn(50),
     +kyzbgn(50),ktxbgn(50),ktybgn(50),ktzbgn(50),nxk(50),nyk(50),
     +nzk(50),ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /cd3cr/
      integer nx,ny,nz,ncx,ncy,ncz
      integer kount(50),ip,ic,ir,ipc,irc,nrel,l,ixy,ixz,iyz
      klevel = kcur
c
c     pre-relax at current finest grid level
c
      do l = 1,iprer
	call relcd3cr(wk)
      end do
c
c     if at coarsest level post-relax
c
      if (kcur .eq. 1) go to 2
c
c     restrict residual to kcur-1
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      ir = krbgn(klevel)
      ixy = kxybgn(klevel)
      ixz = kxzbgn(klevel)
      iyz = kyzbgn(klevel)
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      ncz = nzk(klevel-1)
      irc = krbgn(klevel-1)
c
c     use full weighting with residual restriction
c
      call rescd3cr(nx,ny,nz,wk(ip),wk(ir),wk(ic),ncx,ncy,ncz,wk(ipc),
     +              wk(irc),wk(kps),wk(ixy),wk(ixz),wk(iyz))
c
c     set counter for grid levels to zero
c
      do l = 1,kcur
	kount(l) = 0
      end do
c
c    set new level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c     kcycle control point
c
    1 continue
c
c     post-relax when kcur revisited
c
      if (klevel .eq. kcur) go to 2
c
c     count "hit" at current level
c
      kount(klevel) = kount(klevel)+1
c
c     relax at current level
c
      do l = 1,nrel
	call relcd3cr(wk)
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle(iprer,ipost) complete at klevel
c     inject correction to finer grid
c
	nx = nxk(klevel+1)
	ny = nyk(klevel+1)
	nz = nzk(klevel+1)
	ip = kpbgn(klevel+1)
	ncx = nxk(klevel)
	ncy = nyk(klevel)
	ncz = nzk(klevel)
	ipc = kpbgn(klevel)
	call ccor3(nx,ny,nz,wk(ip),ncx,ncy,ncz,wk(ipc),
     +             nxa,nxb,nyc,nyd,nze,nzf,intpol,wk(kps))
c
c     reset counter to zero at klevel
c
	kount(klevel) = 0
c
c     ascend to next higher level and set to post-relax there
c
	klevel = klevel+1
	nrel = ipost
	go to 1
      else
c
c     kcycle not complete so descend unless at coarsest
c
	if (klevel .gt. 1) then
	  nx = nxk(klevel)
	  ny = nyk(klevel)
	  nz = nzk(klevel)
	  ip = kpbgn(klevel)
	  ic = kcbgn(klevel)
	  ir = krbgn(klevel)
	  ncx = nxk(klevel-1)
	  ncy = nyk(klevel-1)
	  ncz = nzk(klevel-1)
	  irc = krbgn(klevel-1)
	  ipc = kpbgn(klevel-1)
	  ixy = kxybgn(klevel)
	  ixz = kxzbgn(klevel)
	  iyz = kyzbgn(klevel)
c
c     use full weighting with residual restriction
c
	  call rescd3cr(nx,ny,nz,wk(ip),wk(ir),wk(ic),ncx,ncy,ncz,wk(ipc),
     +                  wk(irc),wk(kps),wk(ixy),wk(ixz),wk(iyz))
c
c     pre-relax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 1
	else
c
c     post-relax at coarsest level (klevel=1)
c
	  do l = 1,ipost
	    call relcd3cr(wk)
	  end do
c
c     inject correction to grid level 2
c
	  ipc = kpbgn(1)
	  ncx = nxk(1)
	  ncy = nyk(1)
	  ncz = nzk(1)
	  ip = kpbgn(2)
	  nx = nxk(2)
	  ny = nyk(2)
	  nz = nzk(2)
	  call ccor3(nx,ny,nz,wk(ip),ncx,ncy,ncz,wk(ipc),
     +               nxa,nxb,nyc,nyd,nze,nzf,intpol,wk(kps))
c
c     set to post-relax at level 2
c
	  nrel = ipost
	  klevel = 2
	  go to 1
	end if
      end if
    2 continue
c
c     post-relax at kcur level
c
      do l = 1,ipost
	call relcd3cr(wk)
      end do
      return
      end

      subroutine rescd3cr(nx,ny,nz,phi,rhs,cof,ncx,ncy,ncz,phic,rhsc,
     +resf,coxy,coxz,coyz)
        !dir$ attributes code_align : 32 :: rescd3cr
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: rescd3cr
c
c     compute fully weighted residual restriction in rhsc
c
      use omp_lib
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      complex phi(0:nx+1,0:ny+1,0:nz+1),phic(0:ncx+1,0:ncy+1,0:ncz+1)
      complex rhsc(ncx,ncy,ncz),resf(nx,ny,nz),cof(nx,ny,nz,7)
      complex rhs(nx,ny,nz),coxy(nx,ny,nz),coxz(nx,ny,nz),coyz(nx,ny,nz)
      integer ic,jc,kc,i,j,k
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
     !dir$ attributes align : 64 :: /kcrsxyz/
c
c     initialize phic to zero
c
      !dir$ assume_aligned phic:64
      !dir$ code_align(32)
      do kc=0,ncz+1
	do jc=0,ncy+1
          !dir$ vector aligned
          !dir$ vector always
          !dir$ vector nontemporal(phic)
          !dir$ unroll(8)
	  do ic=0,ncx+1
	    phic(ic,jc,kc) = (0.0,0.0)
	  end do
	end do
      end do
c
c     compute fine grid residual
c
      !dir$ assume_aligned resf:64
      !dir$ assume_aligned rhs:64
      !dir$ assume_aligned cof:64
      !dir$ assume_aligned phi:64
      !dir$ code_align(32)
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,8)
!$OMP&       PRIVATE(i,j,k) SHARED(rhs,resf,cof,phi,nx,ny,nz)
      do k=1,nz
	do j=1,ny
          !dir$ oce_align(32)
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always
	  do i=1,nx
	    resf(i,j,k) =  rhs(i,j,k)-(
     +                     cof(i,j,k,1)*phi(i-1,j,k)+
     +                     cof(i,j,k,2)*phi(i+1,j,k)+
     +                     cof(i,j,k,3)*phi(i,j-1,k)+
     +                     cof(i,j,k,4)*phi(i,j+1,k)+
     +                     cof(i,j,k,5)*phi(i,j,k-1)+
     +                     cof(i,j,k,6)*phi(i,j,k+1)+
     +                     cof(i,j,k,7)*phi(i,j,k))
	  end do
	end do
      end do
c
c   adjust residual with cross coefs as necessary on interior
c
      if (kxy .eq. 1) then
      !dir$ assume_aligned resf:64
      !dir$ assume_aligned coxy:64
      !dir$ assume_aligned phi:64
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,8)
!$OMP&   PRIVATE(i,j,k), SHARED(resf,coxy,phi,nx,ny,nz)
	do k=2,nz-1
	  do j=2,ny-1
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do i=2,nx-1
	    resf(i,j,k) =  resf(i,j,k) - coxy(i,j,k)*(
     +                     phi(i+1,j+1,k)+phi(i-1,j-1,k) -(
     +                     phi(i+1,j-1,k)+phi(i-1,j+1,k)))
	    end do
	  end do
	end do
      end if
      if (kxz .eq. 1) then
      !dir$ assume_aligned resf:64
      !dir$ assume_aligned coxz:64
      !dir$ assume_aligned phi:64
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,8)
!$OMP&  PRIVATE(i,j,k), SHARED(resf,coxz,phi,nx,ny,nz)
	do k=2,nz-1
	  do j=2,ny-1
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do i=2,nx-1
	    resf(i,j,k) =  resf(i,j,k) - coxz(i,j,k)*(
     +                     phi(i+1,j,k+1)+phi(i-1,j,k-1) -(
     +                     phi(i+1,j,k-1)+phi(i-1,j,k+1)))
	    end do
	  end do
	end do
      end if
      if (kyz .eq. 1) then
      !dir$ assume_aligned resf:64
      !dir$ assume_aligned coyz:64
      !dir$ assume_aligned phi:64
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,8)
!$OMP& PRIVATE(i,j,k), SHARED(resf,coyz,phi,nx,ny,nz)
	do k=2,nz-1
	  do j=2,ny-1
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do i=2,nx-1
	    resf(i,j,k) =  resf(i,j,k) - coyz(i,j,k)*(
     +                     phi(i,j+1,k+1)+phi(i,j-1,k-1) - (
     +                     phi(i,j+1,k-1)+phi(i,j-1,k+1)))
	    end do
	  end do
	end do
      end if
c
c     restrict resf to coarse mesh in rhsc
c     using fully weighted residual restriction in 3-d
c
      call cres3(nx,ny,nz,resf,ncx,ncy,ncz,rhsc,nxa,nxb,nyc,nyd,nze,nzf)
      return
      end

      subroutine discd3cr(nx,ny,nz,cof,coxy,coxz,coyz,tx,ty,tz,
     +                    bnd3cr,coef,crsxy,crsxz,crsyz,wk,ier)
        !dir$ attributes code_align : 32 :: discd3cr
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: discd3cr
c
c     discretize the 3-d elliptic pde
c
      use omp_lib
      implicit none
      integer nx,ny,nz,ier
      complex cof(nx,ny,nz,7)
      complex tx(nx,ny,nz,*),ty(ny,nx,nz,*),tz(nz,nx,ny,*),wk(*)
      complex coxy(nx,ny,nz),coxz(nx,ny,nz),coyz(nx,ny,nz)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      common/fcd3cr/xa,xb,yc,yd,ze,zf,tolmax,relmax
      !dir$ attributes align : 64 :: /fcd3cr/
      integer kpbgn,kcbgn,krbgn,kxybgn,kxzbgn,kyzbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      common/cd3cr/kpbgn(50),kcbgn(50),krbgn(50),kxybgn(50),kxzbgn(50),
     +kyzbgn(50),ktxbgn(50),ktybgn(50),ktzbgn(50),nxk(50),nyk(50),
     +nzk(50),ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /cd3cr/
      real dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz,cmin,cemax,alfmax
      complex cxx,cyy,czz,cx,cy,cz,ce,alfa,x,y,z,c1,c2,c3,c4,c5,c6
      complex cxy,cxz,cyz,a,b,c,g
      integer i,j,k,l,ist,ifn,jst,jfn,kst,kfn,kbdy
      integer nxny,nxnz,nynz,im1,jm1,km1
      real odlxy4,odlxz4,odlyz4
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
     !dir$ attributes align : 64 :: /kcrsxyz/
      real cxxr,cxxi,cyyr,cyyi,czzr,czzi
      external bnd3cr,coef,crsxy,crsxz,crsyz
c
c     set current grid increments
c
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlz = (zf-ze)/(nz-1)
      dlz2 = dlz+dlz
      dlzz = dlz*dlz
      odlxy4 = 0.25/(dlx*dly)
      odlxz4 = 0.25/(dlx*dlz)
      odlyz4 = 0.25/(dly*dlz)
      cmin = 1.0
      cemax = 0.0
c
c     set x,y,z subscript limits to bypass specified boundaries
c     when calling coef or bnd3cr
c
      jst = 1
      jfn = ny
      ist = 1
      ifn = nx
      kst = 1
      kfn = nz
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      if (nze.eq.1) kst = 2
      if (nzf.eq.1) kfn = nz-1
      do k=kst,kfn
	z = ze+(k-1)*dlz
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    cxxr = cxx
	    cxxi = aimag(cxx)
	    cyyr = cyy
	    cyyi = aimag(cyy)
	    czzr = czz
	    czzi = aimag(czz)
	    cmin = amin1(cmin,amax1(cxxr*cyyr*czzr,cxxi*cyyi*czzi))
	    cemax = amax1(cabs(ce),cemax)
c
c     flag hyperbolic pde & adjust coeffs if necessary
c
	    if (cabs(cx)*dlx*0.5 .gt. cabs(cxx)) then
	      if (klevel.eq.ngrid) ier = -4
	      cxx = cmplx(cabs(cx)*dlx*0.5,0.0)
	    end if
	    if (cabs(cy)*dly*0.5 .gt. cabs(cyy)) then
	      if (klevel.eq.ngrid) ier = -4
	      cyy = cmplx(cabs(cy)*dly*0.5,0.0)
	    end if
	    if (cabs(cz)*dlz*0.5 .gt. cabs(czz)) then
	      if (klevel.eq.ngrid) ier = -4
	      czz = cmplx(cabs(cz)*dlz*0.5,0.0)
	    end if
	    c1 = cxx/dlxx-cx/dlx2
	    c2 = cxx/dlxx+cx/dlx2
	    c3 = cyy/dlyy-cy/dly2
	    c4 = cyy/dlyy+cy/dly2
	    c5 = czz/dlzz-cz/dlz2
	    c6 = czz/dlzz+cz/dlz2
	    cof(i,j,k,1) = c1
	    cof(i,j,k,2) = c2
	    cof(i,j,k,3) = c3
	    cof(i,j,k,4) = c4
	    cof(i,j,k,5) = c5
	    cof(i,j,k,6) = c6
	    cof(i,j,k,7) = ce-(c1+c2+c3+c4+c5+c6)
	  end do
	end do
      end do

      if (ier .ne. -4) then
c
c     set nonfatal error flag if ellipticity test fails
c
	if (cmin.le.0.0) ier = -2
      end if
c
c     set cross derivative coefficients as necessary
c
      if (kxy.eq.1) then
        !dir$ assume_aligned coxy:64
!$omp parallel default(none) private(k,j,i,z,y,x) 
!$omp& shared(nz,ny,nx,coxy,ze,dlz,yc,dly,xa,dlx,cxy,odlxy4)
!$omp do schedule(static,8)
	do k=1,nz
	  do j=1,ny
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(coxy)
	    do i=1,nx
	      coxy(i,j,k) = (0.0,0.0)
	    end do
	  end do
	end do
!$omp end do
!$omp do schedule(static,8)
	do k=2,nz-1
	  z = ze+(k-1)*dlz
	  do j=2,ny-1
	    y = yc+(j-1)*dly
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do i=2,nx-1
	      x = xa+(i-1)*dlx
	      call crsxy(x,y,z,cxy)
	      coxy(i,j,k) = cxy*odlxy4
	    end do
	  end do
	end do
!$omp  end do
!$omp end parallel
      end if
      if (kxz.eq.1) then
  !dir$ assume_aligned coxz:64
!$omp parallel default(none) private(k,j,i,z,y,x) 
!$omp& shared(nz,ny,nx,ze,coxz,dlz,yc,dly,xa,dlx,cxz,odlxz4)
!$omp do schedule(static,8)
	do k=1,nz
	  do j=1,ny
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(coxz)
	    do i=1,nx
	      coxz(i,j,k) = (0.0,0.0)
	    end do
	  end do
	end do
!$omp end do
!$omp do schedule(static,8)
	do k=2,nz-1
	  z = ze+(k-1)*dlz
	  do j=2,ny-1
	    y = yc+(j-1)*dly
	    do i=2,nx-1
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	      x = xa+(i-1)*dlx
	      call crsxz(x,y,z,cxz)
	      coxz(i,j,k) = cxz*odlxz4
	    end do
	  end do
	end do
!$omp end do
!$omp end parallel
      end if
      if (kyz.eq.1) then
 !dir$ assume_aligned coyz:64
!$omp parallel default(none) private(k,j,i,z,y,x) 
!$omp& shared(nz,ny,nx,ze,coyz,dlz,yc,dly,xa,dlx,cyz,odlyz4)
!$omp do schedule(static,8)
	do k=1,nz
	  do j=1,ny
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(coyz)
	    do i=1,nx
	      coyz(i,j,k) = (0.0,0.0)
	    end do
	  end do
	end do
!$omp end do
!$omp do schedule(static,8)
	do k=2,nz-1
	  z = ze+(k-1)*dlz
	  do j=2,ny-1
	    y = yc+(j-1)*dly
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do i=2,nx-1
	      x = xa+(i-1)*dlx
	      call crsyz(x,y,z,cyz)
	      coyz(i,j,k) = cyz*odlyz4
	    end do
	  end do
	end do
!$omp end do
!$omp end parallel
      end if
      alfmax = 0.0
c
c     adjust coefficients at mixed b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
        
        do k=kst,kfn
	  z = ze+(k-1)*dlz
          !dir$ assume_aligned cof:64
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always
	  do j=jst,jfn
	    y = yc+(j-1)*dly
	    call bnd3cr(kbdy,y,z,a,b,c,g)
	    alfa = c
	    alfmax = amax1(cabs(alfa),alfmax)
	    c1 = cof(i,j,k,1)
	    cof(i,j,k,1) = (0.0,0.0)
	    cof(i,j,k,2) = cof(i,j,k,2)+c1
	    cof(i,j,k,7) = cof(i,j,k,7)+dlx2*alfa*c1
	  end do
	end do
      end if
      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
	do k=kst,kfn
	  z = ze+(k-1)*dlz
          !dir$ assume_aligned cof:64
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always
	  do j=jst,jfn
	    y = yc+(j-1)*dly
	    call bnd3cr(kbdy,y,z,a,b,c,g)
	    alfa = c
	    alfmax = amax1(cabs(alfa),alfmax)
	    c2 = cof(i,j,k,2)
	    cof(i,j,k,1) = cof(i,j,k,1)+c2
	    cof(i,j,k,2) = (0.0,0.0)
	    cof(i,j,k,7) = cof(i,j,k,7)-dlx2*alfa*c2
	  end do
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
	do k=kst,kfn
	  z = ze+(k-1)*dlz
          !dir$ assume_aligned cof:64
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,z,a,b,c,g)
	    alfa = c
	    alfmax = amax1(cabs(alfa),alfmax)
	    c3 = cof(i,j,k,3)
	    cof(i,j,k,3) = (0.0,0.0)
	    cof(i,j,k,4) = cof(i,j,k,4)+c3
	    cof(i,j,k,7) = cof(i,j,k,7)+dly2*alfa*c3
	  end do
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
	do k=kst,kfn
	z = ze+(k-1)*dlz
          !dir$ assume_aligned cof:64
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always 
	do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,z,a,b,c,g)
	    alfa = c
	    alfmax = amax1(cabs(alfa),alfmax)
	    c4 = cof(i,j,k,4)
	    cof(i,j,k,3) = cof(i,j,k,3)+c4
	    cof(i,j,k,4) = (0.0,0.0)
	    cof(i,j,k,7) = cof(i,j,k,7)-dly2*c4*alfa
	  end do
	end do
      end if
      if (nze.eq.2) then
	kbdy = 5
	z = ze
	k = 1
	do j=jst,jfn
	  y = yc+(j-1)*dly
          !dir$ assume_aligned cof:64
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always 
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,y,a,b,c,g)
	    alfa = c
	    alfmax = amax1(cabs(alfa),alfmax)
	    c5 = cof(i,j,k,5)
	    cof(i,j,k,5) = (0.0,0.0)
	    cof(i,j,k,6) = cof(i,j,k,6)+c5
	    cof(i,j,k,7) = cof(i,j,k,7)+dlz2*c5*alfa
	  end do
	end do
      end if
      if (nzf.eq.2) then
	kbdy = 6
	z = zf
	k = nz
	do j=jst,jfn
	  y = yc+(j-1)*dly
          !dir$ assume_aligned cof:64
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always 
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,y,a,b,c,g)
	    alfa = c
	    alfmax = amax1(cabs(alfa),alfmax)
	    c6 = cof(i,j,k,6)
	    cof(i,j,k,5) = cof(i,j,k,5)+c6
	    cof(i,j,k,6) = (0.0,0.0)
	    cof(i,j,k,7) = cof(i,j,k,7)-dlz2*c6*alfa
	  end do
	end do
      end if
c
c     flag continuous singular elliptic pde if detected
c
      if (ier .ne. -4) then
	if (cemax.eq.0.0.and.alfmax.eq.0.0) then
	  if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	    if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
	      if (nze.eq.0.or.(nze.eq.2.and.nzf.eq.2)) then
		ier = -3
	      end if
	    end if
	  end if
	end if
      end if
c
c     reset cof for specified b.c.
c
      if (nxa.eq.1) then
        !dir$ assume_aligned cof:64
	i = 1
	do j=1,ny
	  do k=1,nz
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(cof)
	    do l=1,7
	      cof(i,j,k,l) = (0.0,0.0)
	    end do
	    cof(i,j,k,7) = (1.0,0.0)
	  end do
	end do
      end if
      if (nxb.eq.1) then
        !dir$ assume_aligned cof:64
	i = nx
	do k=1,nz
	  do j=1,ny
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(cof)
	    do l=1,7
	      cof(i,j,k,l) = (0.0,0.0)
	    end do
	    cof(i,j,k,7) = (1.0,0.0)
	  end do
	end do
      end if
      if (nyc.eq.1) then
        !dir$ assume_aligned cof:64
	j = 1
	do k=1,nz
	  do i=1,nx
             !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(cof)
	    do l=1,7
	      cof(i,j,k,l) = (0.0,0.0)
	    end do
	    cof(i,j,k,7) = (1.0,0.0)
	  end do
	end do
      end if
      if (nyd.eq.1) then
          !dir$ assume_aligned cof:64
	j = ny
	do i=1,nx
	  do k=1,nz
             !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(cof)
	    do l=1,7
	      cof(i,j,k,l) = (0.0,0.0)
	    end do
	    cof(i,j,k,7) = 1.0
	  end do
	end do
      end if
      if (nze.eq.1) then
         !dir$ assume_aligned cof:64
	k = 1
	do j=1,ny
	  do i=1,nx
              !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(cof)
	    do l=1,7
	      cof(i,j,k,l) = (0.0,0.0)
	    end do
	    cof(i,j,k,7) = (1.0,0.0)
	  end do
	end do
      end if
      if (nzf.eq.1) then
         !dir$ assume_aligned cof:64
	k = nz
	do j=1,ny
	  do i=1,nx
             !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector nontemporal(cof)
	    do l=1,7
	      cof(i,j,k,l) = (0.0,0.0)
	    end do
	    cof(i,j,k,7) = (1.0,0.0)
	  end do
	end do
      end if


      if (method.eq. 0) return
c
c     set,factor tridiagonal matrices for line relaxation
c
      if ((method-1)*(method-4)*(method-5)*(method-7).eq.0) then
c
c     line relaxation in x used
c
	if (nxa.ne.0) then
c
c     set non-periodic tridiagonal matrices in tx and factor
c
          !dir$ assume_aligned tx:64
          !dir$ assume_aligned cof:64
!$omp parallel do schedule(static,8) default(none) &
!$omp& private(i,im1,k,j) shared(nx,nz,ny,tx,cof)
	  do i=1,nx
	    im1 = max0(i-1,1)
	    do k=1,nz
	      do j=1,ny
		tx(im1,j,k,1) = cof(i,j,k,1)
		tx(i,j,k,2) = cof(i,j,k,7)
		tx(i,j,k,3) = cof(i,j,k,2)
	      end do
	    end do
	  end do
!$omp end parallel do
	  nynz = ny*nz
	  call cfactri(nynz,nx,tx(1,1,1,1),tx(1,1,1,2),tx(1,1,1,3))
	else
	  if (nx .gt. 3) then
c
c     set "periodic" tridiagonal matrices in tx and factor when nx > 3
c
           !dir$ assume_aligned tx:64
          !dir$ assume_aligned cof:64
!$omp parallel do schedule(static,8) default(none) &
!$omp& private(i,k,j) shared(nx,nz,ny,tx,cof)
	    do k=1,nz
	      do j=1,ny
		do i=1,nx-1
		  tx(i,j,k,1) = cof(i,j,k,1)
		  tx(i,j,k,2) = cof(i,j,k,7)
		  tx(i,j,k,3) = cof(i,j,k,2)
		end do
	      end do
	    end do
!$omp end parallel do
	    nynz = ny*nz
	    call cfactrp(nynz,nx,tx,tx(1,1,1,2),tx(1,1,1,3),tx(1,1,1,4),
     +                   tx(1,1,1,5),wk(kps))
	  end if
	end if
      end if
      if ((method-2)*(method-4)*(method-6)*(method-7).eq.0) then
c
c     line relaxation in y used
c
	if (nyc.ne.0) then
c
c     set non-periodic tridiagonal matrices and factor
c
          !dir$ assume_aligned ty:64
          !dir$ assume_aligned cof:64
!$omp parallel do schedule(static,8) default(none) &
!$omp& private(i,jm1,k,j) shared(nx,nz,ny,ty,cof)
	  do j=1,ny
	    jm1 = max0(j-1,1)
	    do k=1,nz
	      do i=1,nx
		ty(jm1,i,k,1) = cof(i,j,k,3)
		ty(j,i,k,2) = cof(i,j,k,7)
		ty(j,i,k,3) = cof(i,j,k,4)
	      end do
	    end do
	  end do
!$omp end parallel do
	  nxnz = nx*nz
	  call cfactri(nxnz,ny,ty(1,1,1,1),ty(1,1,1,2),ty(1,1,1,3))
	else
	  if (ny .gt. 3) then
c
c     set and factor periodic "tridiagonal" matrices when ny > 3
c
          !dir$ assume_aligned ty:64
          !dir$ assume_aligned cof:64
!$omp parallel do schedule(static,8) default(none) &
!$omp& private(i,k,j) shared(nx,nz,ny,ty,cof)
	    do k=1,nz
	      do i=1,nx
		do j=1,ny-1
		  ty(j,i,k,1) = cof(i,j,k,3)
		  ty(j,i,k,2) = cof(i,j,k,7)
		  ty(j,i,k,3) = cof(i,j,k,4)
		end do
	      end do
	    end do
!$omp end parallel do
	    nxnz = nx*nz
	    call cfactrp(nxnz,ny,ty,ty(1,1,1,2),ty(1,1,1,3),ty(1,1,1,4),
     +                   ty(1,1,1,5),wk(kps))
	  end if
	end if
      end if
      if ((method-3)*(method-5)*(method-6)*(method-7).eq.0) then
c
c     line relaxation in z used
c
	if (nze.ne.0) then
c
c     set and factor non-periodic tridiagonal matrices
c
          !dir$ assume_aligned tz:64
          !dir$ assume_aligned cof:64
!$omp parallel do schedule(static,8) default(none) &
!$omp& private(i,km1,k,j) shared(nx,nz,ny,tz,cof)
	  do k=1,nz
	    km1 = max0(k-1,1)
	    do j=1,ny
	      do i=1,nx
		tz(km1,i,j,1) = cof(i,j,k,5)
		tz(k,i,j,2) = cof(i,j,k,7)
		tz(k,i,j,3) = cof(i,j,k,6)
	      end do
	    end do
	  end do
!$omp end parallel do
	  nxny = nx*ny
	  call cfactri(nxny,nz,tz(1,1,1,1),tz(1,1,1,2),tz(1,1,1,3))
	else
	  if (nz .gt. 3) then
c
c     set and factor periodic "tridiagonal matrices when nz > 3
c
          !dir$ assume_aligned tz:64
          !dir$ assume_aligned cof:64
!$omp parallel do schedule(static,8) default(none) &
!$omp& private(i,k,j) shared(nx,nz,ny,tx,cof)
	    do j=1,ny
	      do i=1,nx
		do k=1,nz-1
		  tz(k,i,j,1) = cof(i,j,k,5)
		  tz(k,i,j,2) = cof(i,j,k,7)
		  tz(k,i,j,3) = cof(i,j,k,6)
		end do
	      end do
	    end do
!$omp end parallel do
	    nxny = nx*ny
	    call cfactrp(nxny,nz,tz(1,1,1,1),tz(1,1,1,2),tz(1,1,1,3),
     +                   tz(1,1,1,4),tz(1,1,1,5),wk(kps))
	  end if
	end if
      end if
      return
      end

      subroutine adjcd3cr(nx,ny,nz,phi,rhs,phif,bnd3cr,coef)
          !dir$ attributes code_align : 32 :: adjcd3cr
          !dir$ optimize : 3
          !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: adjcd3cr
c
c     adjust for mixed b.c. or specified boundaries
c     approximate nonnormal b.c. derivatives with finite differences
c     applied to previous fine grid level guess in phif)
c
      implicit none
      integer nx,ny,nz
      integer i,j,k,kbdy,ifine,jfine,kfine,ist,ifn,jst,jfn,kst,kfn
      complex phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
      complex phif(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/fcd3cr/xa,xb,yc,yd,ze,zf,tolmax,relmax
      !dir$ attributes align : 64 :: /fcd3cr/
      integer kpbgn,kcbgn,krbgn,kxybgn,kxzbgn,kyzbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      common/cd3cr/kpbgn(50),kcbgn(50),krbgn(50),kxybgn(50),kxzbgn(50),
     +kyzbgn(50),ktxbgn(50),ktybgn(50),ktzbgn(50),nxk(50),nyk(50),
     +nzk(50),ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /cd3cr/
      real  dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      common / incr3 / dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      !dir$ attributes align : 64 :: /incr3/
      real dlx,dlx2,dlxx,dly,dly2,dlyy,dlz,dlz2,dlzz
      real odlxy4,odlxz4,odlyz4
      complex c1,c2,c3,c4,c5,c6
      real x,y,z
      complex cxx,cyy,czz,cx,cy,cz,ce
      complex px,py,pz,a,b,c,g,gbdy
      external bnd3cr,coef
c
c     set current grid increments
c
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlz = (zf-ze)/(nz-1)
      dlz2 = dlz+dlz
      dlzz = dlz*dlz
      odlxy4 = 0.25/(dlx*dly)
      odlxz4 = 0.25/(dlx*dlz)
      odlyz4 = 0.25/(dly*dlz)
c
c     set x,y,z subscript limits for calls to coef,bnd3cr
c
      jst = 1
      jfn = ny
      ist = 1
      ifn = nx
      kst = 1
      kfn = nz
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      if (nze.eq.1) kst = 2
      if (nzf.eq.1) kfn = nz-1
c
c     adjust mixed derivaative
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
	do k=kst,kfn
	  z = ze+(k-1)*dlz
	  do j=jst,jfn
	    y = yc+(j-1)*dly
c
c     set fine mesh subscripts for current (x,y,z)
c
	    ifine = 1
	    jfine = int((y-yc)/dy+0.5)+1
	    kfine = int((z-ze)/dz+0.5)+1
	    call bnd3cr(kbdy,y,z,a,b,c,g)
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    c1 = cxx/dlxx-cx/dlx2
c
c     estimate non-normal first derivatives and adjust b.c. rhs
c
	    call cdifxyz(nfx,nfy,nfz,phif,ifine,jfine,kfine,px,py,pz)
	    gbdy = g-a*py-b*pz
	    rhs(i,j,k) = rhs(i,j,k)+dlx2*c1*gbdy
	  end do
	end do
      end if
      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
       	do k=kst,kfn
	  z = ze+(k-1)*dlz
          do j=jst,jfn
	    y = yc+(j-1)*dly
	    call bnd3cr(kbdy,y,z,a,b,c,g)
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    c2 = cxx/dlxx+cx/dlx2
	    ifine = nfx
	    jfine = int((y-yc)/dy+0.5)+1
	    kfine = int((z-ze)/dz+0.5)+1
	    call cdifxyz(nfx,nfy,nfz,phif,ifine,jfine,kfine,px,py,pz)
	    gbdy = g-a*py-b*pz
	    rhs(i,j,k) = rhs(i,j,k)-dlx2*gbdy*c2
	  end do
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
	do k=kst,kfn
	  z = ze+(k-1)*dlz
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,z,a,b,c,g)
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    c3 = cyy/dlyy-cy/dly2
	    ifine = int((x-xa)/dx+0.5)+1
	    jfine = 1
	    kfine = int((z-ze)/dz+0.5)+1
	    call cdifxyz(nfx,nfy,nfz,phif,ifine,jfine,kfine,px,py,pz)
	    gbdy = g-a*px-b*pz
	    rhs(i,j,k) = rhs(i,j,k)+dly2*c3*gbdy
	  end do
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
	do k=kst,kfn
	  z = ze+(k-1)*dlz
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,z,a,b,c,g)
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    c4 = cyy/dlyy+cy/dly2
	    ifine = int((x-xa)/dx+0.5)+1
	    jfine = nfy
	    kfine = int((z-ze)/dz+0.5)+1
	    call cdifxyz(nfx,nfy,nfz,phif,ifine,jfine,kfine,px,py,pz)
	    gbdy = g-a*px-b*pz
	    rhs(i,j,k) = rhs(i,j,k)-dly2*c4*gbdy
	  end do
	end do
      end if
      if (nze.eq.2) then
	kbdy = 5
	z = ze
	k = 1
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,y,a,b,c,g)
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    c5 = czz/dlzz-cz/dlz2
	    ifine = int((x-xa)/dx+0.5)+1
	    jfine = int((y-yc)/dy+0.5)+1
	    kfine = 1
	    call cdifxyz(nfx,nfy,nfz,phif,ifine,jfine,kfine,px,py,pz)
	    gbdy = g-a*px-b*py
	    rhs(i,j,k) = rhs(i,j,k)+dlz2*c5*gbdy
	  end do
	end do
      end if
      if (nzf.eq.2) then
	kbdy = 6
	z = zf
	k = nz
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  do i=ist,ifn
	    x = xa+(i-1)*dlx
	    call bnd3cr(kbdy,x,y,a,b,c,g)
	    call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    c6 = czz/dlzz+cz/dlz2
	    ifine = int((x-xa)/dx+0.5)+1
	    jfine = int((y-yc)/dy+0.5)+1
	    kfine = nfz
	    call cdifxyz(nfx,nfy,nfz,phif,ifine,jfine,kfine,px,py,pz)
	    gbdy = g-a*px-b*py
	    rhs(i,j,k) = rhs(i,j,k)-dlz2*c6*gbdy
	  end do
	end do
      end if
c
c     adjust for specified (dirchlet) boundary conditions
c
      if (nxa.eq.1) then
	i = 1
  	do j=1,ny
    	  do k=1,nz
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nxb.eq.1) then
	i = nx

	do j=1,ny
	  do k=1,nz
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do k=1,nz
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned phi:64
          !dir$ vector aligned
          !dir$ vector always
          !dir$ unroll(16)
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do k=1,nz
           !dir$ assume_aligned rhs:64
          !dir$ assume_aligned phi:64
          !dir$ vector aligned
          !dir$ vector always
          !dir$ unroll(16)
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nze.eq.1) then
	k = 1
	do j=1,ny
           !dir$ assume_aligned rhs:64
          !dir$ assume_aligned phi:64
          !dir$ vector aligned
          !dir$ vector always
          !dir$ unroll(16)
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nzf.eq.1) then
	k = nz
	do j=1,ny
           !dir$ assume_aligned rhs:64
          !dir$ assume_aligned phi:64
          !dir$ vector aligned
          !dir$ vector always
          !dir$ unroll(16)
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      return
      end

      subroutine cdifxyz(nx,ny,nz,p,i,j,k,px,py,pz)
        !dir$ attributes forceinline :: cdifxyz
        !dir$ attributes code_align : 32 :: cdifxyz
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: cdifxyz
c
c     estimate first order partial derivatives at (i,j,k) fine grid point
c     in px,py,pz using second order difference formula applied to p
c     these estimates are used to modify nonnormal derivative b.c.
c
      implicit none
      integer nx,ny,nz,i,j,k
      complex p(nx,ny,nz),px,py,pz
      real  dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      common / incr3 / dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
     
      if (i.gt.1 .and. i.lt.nx) then
	px = (p(i+1,j,k)-p(i-1,j,k))/dx2
      else if (i.eq.1) then
	px = (-3.*p(1,j,k)+4.*p(2,j,k)-p(3,j,k))/dx2
      else
	px = (3.*p(nx,j,k)-4.*p(nx-1,j,k)+p(nx-2,j,k))/dx2
      end if
      if (j.gt.1 .and. j.lt.ny) then
	py = (p(i,j+1,k)-p(i,j-1,k))/dy2
      else if (j.eq.1) then
	py = (-3.*p(i,1,k)+4.*p(i,2,k)-p(i,3,k))/dy2
      else
	py = (3.*p(i,ny,k)-4.*p(i,ny-1,k)+p(i,ny-2,k))/dy2
      end if
      if (k.gt.1 .and. k.lt.nz) then
	pz = (p(i,j,k+1)-p(i,j,k-1))/dz2
      else if (k.eq.1) then
	pz = (-3.*p(i,j,1)+4.*p(i,j,2)-p(i,j,3))/dz2
      else
	pz = (3.*p(i,j,nz)-4.*p(i,j,nz-1)+p(i,j,nz-2))/dz2
      end if
      return
      end

      subroutine relcd3cr(wk)
        !dir$ attributes code_align : 32 :: relcd3cr
        !dir$ optimize : 3
c
c     point or line relaxation in the x and/or y and/or z direction(s)
c
      implicit none
      complex wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      common/fcd3cr/xa,xb,yc,yd,ze,zf,tolmax,relmax
      !dir$ attributes align : 64 :: /fcd3cr/
      integer kpbgn,kcbgn,krbgn,kxybgn,kxzbgn,kyzbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      common/cd3cr/kpbgn(50),kcbgn(50),krbgn(50),kxybgn(50),kxzbgn(50),
     +kyzbgn(50),ktxbgn(50),ktybgn(50),ktzbgn(50),nxk(50),nyk(50),
     +nzk(50),ngrid,klevel,kcur,kps
      !dir$ attributes align : 64 :: /cd3cr/
      integer nx,ny,nz,ip,ic,ir,m,itx,ity,itz,ixy,ixz,iyz
      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      ir = krbgn(klevel)
      ixy = kxybgn(klevel)
      ixz = kxzbgn(klevel)
      iyz = kyzbgn(klevel)
      if (method.eq.0) then
c
c     gauss-seidel pointwise red/black relaxation
c
	call relcd3crp(nx,ny,nz,wk(ip),wk(ic),wk(ixy),wk(ixz),
     +                 wk(iyz),wk(ir))
	return
      end if
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      itz = ktzbgn(klevel)
      m = method
c
c     check for line relaxation(s) (in combinations)
c
      if ((m-1)*(m-4)*(m-5)*(m-7) .eq. 0 ) then
c
c     line - x relaxation
c
	if (nxa .ne. 0 .or. nx .gt. 3) then
	 itx = ktxbgn(klevel)
	 ixy = kxybgn(klevel)
	 ixz = kxzbgn(klevel)
	 iyz = kyzbgn(klevel)
	 call slxcd3cr(nx,ny,nz,wk(ip),wk(ir),wk(ic),wk(itx),wk(kps),
     +                 wk(ixy),wk(ixz),wk(iyz),nxa,nyc,nze)
	else
c
c     replace by point if x-periodic and nx=3
c
	call relcd3crp(nx,ny,nz,wk(ip),wk(ic),wk(ixy),wk(ixz),
     +                 wk(iyz),wk(ir))
	end if
	if (method .eq. 1) return
      end if
      if ((m-2)*(m-4)*(m-6)*(m-7) .eq. 0 ) then
c
c     line - y relaxation
c
	if (nyc .ne. 0 .or. ny .gt. 3) then
	 ity = ktybgn(klevel)
	 ixy = kxybgn(klevel)
	 ixz = kxzbgn(klevel)
	 iyz = kyzbgn(klevel)
	 call slycd3cr(nx,ny,nz,wk(ip),wk(ir),wk(ic),wk(ity),wk(kps),
     +                 wk(ixy),wk(ixz),wk(iyz),nxa,nyc,nze)
	else
c
c     replace by point if y-periodic and ny=3
c
	call relcd3crp(nx,ny,nz,wk(ip),wk(ic),wk(ixy),wk(ixz),
     +                 wk(iyz),wk(ir))
	end if
	if ((m-2)*(m-4) .eq. 0) return
      end if
      if ((m-3)*(m-5)*(m-6)*(m-7) .eq. 0 ) then
c
c     line - z relaxation
c
	if (nze .ne. 0 .or. nz .gt. 3) then
	 itz = ktzbgn(klevel)
	 ixy = kxybgn(klevel)
	 ixz = kxzbgn(klevel)
	 iyz = kyzbgn(klevel)
	 call slzcd3cr(nx,ny,nz,wk(ip),wk(ir),wk(ic),wk(itz),wk(kps),
     +                 wk(ixy),wk(ixz),wk(iyz),nxa,nyc,nze)
	else
c
c     replace by point if z-periodic and nz=3
c
	call relcd3crp(nx,ny,nz,wk(ip),wk(ic),wk(ixy),wk(ixz),
     +                 wk(iyz),wk(ir))
	end if
	return
      end if
      return
      end

      subroutine relcd3crp(nx,ny,nz,phi,cof,coxy,coxz,coyz,rhs)
        !dir$ attributes code_align : 32 :: relcd3crp
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: relcd3crp
c
c     gauss-seidel point relaxation with red/black ordering
c     in three dimensions for nonseparable pde with cross terms
c
      use omp_lib
      implicit none
      integer nx,ny,nz
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,7),rhs(nx,ny,nz)
      complex coxy(nx,ny,nz),coxz(nx,ny,nz),coyz(nx,ny,nz)
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
      !dir$ attributes align : 64 :: /kcrsxyz/
      integer i,j,k,nper
c
c     set periodic b.c. indicator
c
      nper = nxa*nyc*nze
c
c     set periodic boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     relax in order:
c     (1) red (x,y) on odd z planes
c     (2) black (x,y) on even z planes
c     (3) black (x,y) on odd z planes
c     (4) red (x,y) on even z planes
c
!$omp parallel do schedule(static,8) private(i,j,k) shared(phi,cof,rhs,coxy,coxz,coyz)
!$omp& shared(nx,ny,nz,kxy,kxz,kyz)
      do k=1,nz,2
c
c     red (x,y) points on odd z planes
c
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned rhs:64
        !dir$ assume_aligned cof:64
        !dir$ code_align(32)
       	do i=1,nx,2
          !dir$ ivdep
          !dir$ vector aligned
          !dir$ vector always
	  do j=1,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
c
c     adjust for cross derivative coefficients as necessary
c
	if (kxy.eq.1) then
        !dir$ assume_aligned phi:64
        !dir$ assume_aligned coxy:64
        !dir$ assume_aligned cof:64
        !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=2,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
c
c     adjust for cross derivative coefficients as necessary
c
	if (kxy.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
      end do
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c    black (x,y) points on even z planes
c
!$omp parallel do schedule(static,8) private(i,j,k), shared(phi,cof,rhs,coxy,coxz,coyz)
!$omp& shared(nx,ny,nz,kxy,kxz,kyz)
      do k=2,nz,2
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=2,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
	if (kxy.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
             !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=1,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
	if (kxy.eq.1) then
           !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
      end do
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     black (x,y) points on odd z planes
c
!$omp parallel do schedule(static,8) private(i,j,k) shared(phi,cof,rhs,coxy,coxz,coyz)
!$omp& shared(nx,ny,nz,kxy,kxz,kyz)
      do k=1,nz,2
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=2,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
	if (kxy.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
             !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
         !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=1,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
	if (kxy.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
             !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
             !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
      end do
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     red (x,y) points on even z planes
c
!$omp parallel do schedule(static,8) private(i,j,k), shared(phi,cof,rhs,coxy,coxz,coyz)
!$omp& shared(nx,ny,nz,kxy,kxz,kyz)
      do k=2,nz,2
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=1,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
	if (kxy.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
           !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=1,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned rhs:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	  do j=2,ny,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
	  end do
	end do
	if (kxy.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxy:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxy(i,j,k)*(
     +        phi(i-1,j-1,k)+phi(i+1,j+1,k)-(
     +        phi(i-1,j+1,k)+phi(i+1,j-1,k))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kxz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coxz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
             !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coxz(i,j,k)*(
     +        phi(i-1,j,k-1)+phi(i+1,j,k+1)-(
     +        phi(i-1,j,k+1)+phi(i+1,j,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
	if (kyz.eq.1) then
          !dir$ assume_aligned phi:64
          !dir$ assume_aligned coyz:64
          !dir$ assume_aligned cof:64
          !dir$ code_align(32)
	  do i=2,nx,2
            !dir$ ivdep
            !dir$ vector aligned
            !dir$ vector always
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k) - (coyz(i,j,k)*(
     +        phi(i,j-1,k-1)+phi(i,j+1,k+1)-(
     +        phi(i,j-1,k+1)+phi(i,j+1,k-1))))/cof(i,j,k,7)
	    end do
	  end do
	end if
      end do
c
c     final set of periodic virtual boundaries if necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      return
      end

      subroutine csexyzb(nx,ny,nz,crsxy,crsxz,crsyz,cxyxa,cxyxb,
     +cxyyc,cxyyd,cxyze,cxyzf,cxzxa,cxzxb,cxzyc,cxzyd,cxzze,
     +cxzzf,cyzxa,cyzxb,cyzyc,cyzyd,cyzze,cyzzf)
        !dir$ attributes code_align : 32 :: csexyzb
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: csexyzb
c
c     set cross derivative coefficient values along boundaries
c     on fine grid where nonzero cross coefs are flagged
c
      implicit none
      integer nx,ny,nz,i,j,k
      complex cxyxa(ny,nz),cxyxb(ny,nz)
      complex cxyyc(nx,nz),cxyyd(nx,nz)
      complex cxyze(nx,ny),cxyzf(nx,ny)
      complex cxzxa(ny,nz),cxzxb(ny,nz)
      complex cxzyc(nx,nz),cxzyd(nx,nz)
      complex cxzze(nx,ny),cxzzf(nx,ny)
      complex cyzxa(ny,nz),cyzxb(ny,nz)
      complex cyzyc(nx,nz),cyzyd(nx,nz)
      complex cyzze(nx,ny),cyzzf(nx,ny)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/fcd3cr/xa,xb,yc,yd,ze,zf,tolmax,relmax
      !dir$ attributes align : 64 :: /icd3cr/
      real dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      common / incr3 / dx,dy,dz,dx2,dy2,dz2,odxy4,odxz4,odyz4
      !dir$ attributes align : 64 :: /icd3cr/
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
     !dir$ attributes align : 64 :: /kcrsxyz/
      real x,y,z
      complex cxy,cxz,cyz
      external crsxy,crsxz,crsyz
      if (nxa.ne.1) then
	x = xa
	i = 1
	do k=1,nz
	  z = ze+(k-1)*dz
	  do j=1,ny
	    y = yc+(j-1)*dy
	    if (kxy.eq.1) then
	      call crsxy(x,y,z,cxy)
	      cxyxa(j,k) = cxy*odxy4
	    end if
	    if (kxz.eq.1) then
	      call crsxz(x,y,z,cxz)
	      cxzxa(j,k) = cxz*odxz4
	    end if
	    if (kyz.eq.1) then
	      call crsyz(x,y,z,cyz)
	      cyzxa(j,k) = cyz*odyz4
	    end if
	  end do
	end do
      end if
      if (nxb.ne.1) then
	x = xb
	i = nx
	do k=1,nz
	  z = ze+(k-1)*dz
	  do j=1,ny
	    y = yc+(j-1)*dy
	    if (kxy.eq.1) then
	      call crsxy(x,y,z,cxy)
	      cxyxb(j,k) = cxy*odxy4
	    end if
	    if (kxz.eq.1) then
	      call crsxz(x,y,z,cxz)
	      cxzxb(j,k) = cxz*odxz4
	    end if
	    if (kyz.eq.1) then
	      call crsyz(x,y,z,cyz)
	      cyzxb(j,k) = cyz*odyz4
	    end if
	  end do
	end do
      end if
      if (nyc.ne.1) then
	y = yc
	j = 1
	do k=1,nz
	  z = ze+(k-1)*dz
	  do i=1,nx
	    x = xa+(i-1)*dx
	    if (kxy.eq.1) then
	      call crsxy(x,y,z,cxy)
	      cxyyc(i,k) = cxy*odxy4
	    end if
	    if (kxz.eq.1) then
	      call crsxz(x,y,z,cxz)
	      cxzyc(i,k) = cxz*odxz4
	    end if
	    if (kyz.eq.1) then
	      call crsyz(x,y,z,cyz)
	      cyzyc(i,k) = cyz*odyz4
	    end if
	  end do
	end do
      end if
      if (nyd.ne.1) then
	y = yd
	j = ny
	do k=1,nz
	  z = ze+(k-1)*dz
	  do i=1,nx
	    x = xa+(i-1)*dx
	    if (kxy.eq.1) then
	      call crsxy(x,y,z,cxy)
	      cxyyd(i,k) = cxy*odxy4
	    end if
	    if (kxz.eq.1) then
	      call crsxz(x,y,z,cxz)
	      cxzyd(i,k) = cxz*odxz4
	    end if
	    if (kyz.eq.1) then
	      call crsyz(x,y,z,cyz)
	      cyzyd(i,k) = cyz*odyz4
	    end if
	  end do
	end do
      end if
      if (nze.ne.1) then
	z = ze
	k = 1
	do j=1,ny
	  y = yc+(j-1)*dy
	  do i=1,nx
	    x = xa+(i-1)*dx
	    if (kxy.eq.1) then
	      call crsxy(x,y,z,cxy)
	      cxyze(i,j) = cxy*odxy4
	    end if
	    if (kxz.eq.1) then
	      call crsxz(x,y,z,cxz)
	      cxzze(i,j) = cxz*odxz4
	    end if
	    if (kyz.eq.1) then
	      call crsyz(x,y,z,cyz)
	      cyzze(i,j) = cyz*odyz4
	    end if
	  end do
	end do
      end if
      if (nzf.ne.1) then
	z = zf
	k = nz
	do j=1,ny
	  y = yc+(j-1)*dy
	  do i=1,nx
	    x = xa+(i-1)*dx
	    if (kxy.eq.1) then
		call crsxy(x,y,z,cxy)
	    cxyzf(i,j) = cxy*odxy4
	    end if
	    if (kxz.eq.1) then
	      call crsxz(x,y,z,cxz)
	      cxzzf(i,j) = cxz*odxz4
	    end if
	    if (kyz.eq.1) then
	      call crsyz(x,y,z,cyz)
	      cyzzf(i,j) = cyz*odyz4
	    end if
	  end do
	end do
      end if
      return
      end

      subroutine csubxy(nx,ny,nz,p,r,cxyxa,cxyxb,cxyyc,cxyyd,cxyze,cxyzf)
        !dir$ attributes forceinline :: csubxy
        !dir$ attributes code_align : 32 :: csubxy
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: csubxy
c
c     adjust right hand side r by subtracting off second order finite
c     difference approximations to pxy over nonspecified boundaries
c     p holds the current solution estimate to be updated
c
      implicit none
      integer nx,ny,nz,i,j,k,ks,kf,im1,ip1,jm1,jp1
      complex cxyxa(ny,nz),cxyxb(ny,nz)
      complex cxyyc(nx,nz),cxyyd(nx,nz)
      complex cxyze(nx,ny),cxyzf(nx,ny)
      complex p(nx,ny,nz),r(nx,ny,nz)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      complex px(3),pxy
      ks = 1
      kf = nz
      if (nze.eq.1) ks = 2
      if (nzf.eq.1) kf = nz-1
      if (nxa.eq.2) then
	i = 1
	if (nyc.eq.2) then
	  j = 1
	  do k=ks,kf
	    px(1) = -3.*p(1,1,k)+4.*p(2,1,k)-p(3,1,k)
	    px(2) = -3.*p(1,2,k)+4.*p(2,2,k)-p(3,2,k)
	    px(3) = -3.*p(1,3,k)+4.*p(2,3,k)-p(3,3,k)
	    pxy = -3.*px(1) + 4.*px(2)- px(3)
	    r(i,j,k) = r(i,j,k) - cxyxa(j,k)*pxy
	  end do
	end if
	do j=2,ny-1
	  do k=ks,kf
	    px(1) = -3.*p(1,j-1,k)+4.*p(2,j-1,k)-p(3,j-1,k)
	    px(3) = -3.*p(1,j+1,k)+4.*p(2,j+1,k)-p(3,j+1,k)
	    pxy = px(3) - px(1)
	    r(i,j,k) = r(i,j,k) - cxyxa(j,k)*pxy
	  end do
	end do
	if (nyd.eq.2) then
	  j = ny
	  do k=ks,kf
	    px(1) = -3.*p(1,ny-2,k)+4.*p(2,ny-2,k)-p(3,ny-2,k)
	    px(2) = -3.*p(1,ny-1,k)+4.*p(2,ny-1,k)-p(3,ny-1,k)
	    px(3) = -3.*p(1,ny,k)+4.*p(2,ny,k)-p(3,ny,k)
	    pxy = 3.*px(3) - 4.*px(2) + px(1)
	    r(i,j,k) = r(i,j,k) - cxyxa(j,k)*pxy
	  end do
	end if
      end if
      if (nxb.eq.2) then
	i = nx
	if (nyc.eq.2) then
	  j = 1
	  do k=ks,kf
	    px(1) = 3.*p(nx,1,k)-4.*p(nx-1,1,k)+p(nx-2,1,k)
	    px(2) = 3.*p(nx,2,k)-4.*p(nx-1,2,k)+p(nx-2,2,k)
	    px(3) = 3.*p(nx,3,k)-4.*p(nx-1,3,k)+p(nx-2,3,k)
	    pxy = -3.*px(1) + 4.*px(2)- px(3)
	    r(i,j,k) = r(i,j,k) - cxyxb(j,k)*pxy
	  end do
	end if
	do j=2,ny-1
	  do k=ks,kf
	    px(1) = 3.*p(nx,j-1,k)-4.*p(nx-1,j-1,k)+p(nx-2,j-1,k)
	    px(3) = 3.*p(nx,j+1,k)-4.*p(nx-1,j+1,k)+p(nx-2,j+1,k)
	    pxy = px(3) - px(1)
	    r(i,j,k) = r(i,j,k) - cxyxb(j,k)*pxy
	  end do
	end do
	if (nyd.eq.2) then
	  j = ny
	  do k=ks,kf
	    px(1) = 3.*p(nx,ny-2,k)-4.*p(nx-1,ny-2,k)+p(nx-2,ny-2,k)
	    px(2) = 3.*p(nx,ny-1,k)-4.*p(nx-1,ny-1,k)+p(nx-2,ny-1,k)
	    px(3) = 3.*p(nx,ny,k)-4.*p(nx-1,ny,k)+p(nx-2,ny,k)
	    pxy = 3.*px(3) - 4.*px(2) + px(1)
	    r(i,j,k) = r(i,j,k) - cxyxb(j,k)*pxy
	  end do
	end if
      end if
      if (nyc.eq.2) then
	j = 1
	do i=2,nx-1
	  do k=ks,kf
	    px(1) = p(i+1,1,k)-p(i-1,1,k)
	    px(2) = p(i+1,2,k)-p(i-1,2,k)
	    px(3) = p(i+1,3,k)-p(i-1,3,k)
	    pxy = -3.*px(1)+4.*px(2)-px(3)
	    r(i,j,k) = r(i,j,k) - cxyyc(i,k)*pxy
	  end do
	end do
      end if
      if (nyd.eq.2) then
	j = ny
	do i=2,nx-1
	  do k=ks,kf
	    px(1) = p(i+1,ny-2,k)-p(i-1,ny-2,k)
	    px(2) = p(i+1,ny-1,k)-p(i-1,ny-1,k)
	    px(3) = p(i+1,ny,k)-p(i-1,ny,k)
	    pxy = 3.*px(3)-4.*px(2)+px(1)
	    r(i,j,k) = r(i,j,k) - cxyyd(i,k)*pxy
	  end do
	end do
      end if
c
c     handle periodic b.c. using centered difference formula
c
      if (nxa.eq.0) then
	ip1 = 2
	im1 = nx-1
	if (nyc.eq.2) then
	  j = 1
	  do k=ks,kf
	    px(1) = p(ip1,1,k)-p(im1,1,k)
	    px(2) = p(ip1,2,k)-p(im1,2,k)
	    px(3) = p(ip1,3,k)-p(im1,3,k)
	    pxy = -3.*px(1) + 4.*px(2)- px(3)
	    r(1,j,k) = r(1,j,k) - cxyxa(j,k)*pxy
	    r(nx,j,k) = r(nx,j,k) - cxyxb(j,k)*pxy
	  end do
	end if
	do  j=2,ny-1
	  do k=ks,kf
	    pxy = p(ip1,j+1,k)+p(im1,j-1,k)-(p(ip1,j-1,k)+p(im1,j+1,k))
	    r(1,j,k) = r(1,j,k) - cxyxa(j,k)*pxy
	    r(nx,j,k) = r(nx,j,k) - cxyxb(j,k)*pxy
	  end do
	end do
	if (nyd.eq.2) then
	  j = ny
	  do k=ks,kf
	    px(1) = p(ip1,ny-2,k)-p(im1,ny-2,k)
	    px(2) = p(ip1,ny-1,k)-p(im1,ny-1,k)
	    px(3) = p(ip1,ny,k)-p(im1,ny,k)
	    pxy = 3.*px(3) - 4.*px(2) + px(1)
	    r(1,j,k) = r(1,j,k) - cxyxa(j,k)*pxy
	    r(nx,j,k) = r(nx,j,k) - cxyxb(j,k)*pxy
	  end do
	end if
c
c     double periodic
c
	if (nyc.eq.0) then
	  jm1 = ny-1
	  jp1 = 2
	  do k=ks,kf
	    pxy = p(ip1,jp1,k)+p(im1,jm1,k)-(p(ip1,jm1,k)+p(im1,jp1,k))
	    r(1,1,k) = r(1,1,k) - cxyxa(1,k)*pxy
	    r(1,ny,k) = r(1,ny,k) - cxyxa(ny,k)*pxy
	    r(nx,1,k) = r(nx,1,k) - cxyxb(1,k)*pxy
	    r(nx,ny,k) = r(nx,ny,k) - cxyxb(ny,k)*pxy
	  end do
	end if
      end if
      if (nyc.eq.0) then
	jm1 = ny-1
	jp1 = 2
	do i=2,nx-1
	  do k=ks,kf
	    pxy = p(i+1,jp1,k)+p(i-1,jm1,k)-(p(i+1,jm1,k)+p(i-1,jp1,k))
	    r(i,1,k) = r(i,1,k) - cxyyc(i,k)*pxy
	    r(i,ny,k) = r(i,ny,k) - cxyyd(i,k)*pxy
	  end do
	end do
      end if
      if (nze.ne.1) then
	k = 1
	do i=2,nx-1
	  do j=2,ny-1
	    pxy = p(i+1,j+1,k)+p(i-1,j-1,k)-(p(i+1,j-1,k)+p(i-1,j+1,k))
	    r(i,j,k) = r(i,j,k) - cxyze(i,j)*pxy
	  end do
	end do
      end if
      if (nzf.ne.1) then
	k = nz
	do i=2,nx-1
	  do j=2,ny-1
	    pxy = p(i+1,j+1,k)+p(i-1,j-1,k)-(p(i+1,j-1,k)+p(i-1,j+1,k))
	    r(i,j,k) = r(i,j,k) - cxyzf(i,j)*pxy
	  end do
	end do
      end if
      return
      end

      subroutine csubxz(nx,ny,nz,p,r,cxzxa,cxzxb,cxzyc,cxzyd,cxzze,cxzzf)
        !dir$ attributes forceinline :: csubxz
        !dir$ attributes code_align : 32 :: csubxz
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: csubxz
c
c     adjust right hand side r by subtracting off second order finite
c     difference approximations to pxz over nonspecified boundaries
c
      implicit none
      integer nx,ny,nz,i,j,k,js,jf,im1,ip1,km1,kp1
      complex cxzxa(ny,nz),cxzxb(ny,nz)
      complex cxzyc(nx,nz),cxzyd(nx,nz)
      complex cxzze(nx,ny),cxzzf(nx,ny)
      complex p(nx,ny,nz),r(nx,ny,nz)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      complex px(3),pxz
      js = 1
      jf = ny
      if (nyc.eq.1) js = 2
      if (nyd.eq.1) jf = ny-1
      if (nxa.eq.2) then
	i = 1
	if (nze.eq.2) then
	  k = 1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = -3.*p(1,j,1)+4.*p(2,j,1)-p(3,j,1)
	    px(2) = -3.*p(1,j,2)+4.*p(2,j,2)-p(3,j,2)
	    px(3) = -3.*p(1,j,3)+4.*p(2,j,3)-p(3,j,3)
	    pxz = -3.*px(1) + 4.*px(2)- px(3)
	    r(1,j,1) = r(1,j,1) - cxzxa(j,k)*pxz
	  end do
	end if
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
	do k=2,nz-1
            !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = -3.*p(1,j,k-1)+4.*p(2,j,k-1)-p(3,j,k-1)
	    px(3) = -3.*p(1,j,k+1)+4.*p(2,j,k+1)-p(3,j,k+1)
	    pxz = px(3) - px(1)
	    r(1,j,k) = r(1,j,k) - cxzxa(j,k)*pxz
	  end do
	end do
	if (nzf.eq.2) then
	  k = nz
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = -3.*p(1,j,nz-2)+4.*p(2,j,nz-2)-p(3,j,nz-2)
	    px(2) = -3.*p(1,j,nz-1)+4.*p(2,j,nz-1)-p(3,j,nz-1)
	    px(3) = -3.*p(1,j,nz)+4.*p(2,j,nz)-p(3,j,nz)
	    pxz = 3.*px(3) - 4.*px(2) + px(1)
	    r(1,j,nz) = r(1,j,nz) - cxzxa(j,k)*pxz
	  end do
	end if
      end if
      if (nxb.eq.2) then
	i = nx
	if (nze.eq.2) then
	  k = 1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxb:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = 3.*p(nx,j,1)-4.*p(nx-1,j,1)+p(nx-2,j,1)
	    px(2) = 3.*p(nx,j,2)-4.*p(nx-1,j,2)+p(nx-2,j,2)
	    px(3) = 3.*p(nx,j,3)-4.*p(nx-1,j,3)+p(nx-2,j,3)
	    pxz = -3.*px(1) + 4.*px(2)- px(3)
	    r(i,j,k) = r(i,j,k) - cxzxb(j,k)*pxz
	  end do
	end if
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxb:64
	do k=2,nz-1
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = 3.*p(nx,j,k-1)-4.*p(nx-1,j,k-1)+p(nx-2,j,k-1)
	    px(3) = 3.*p(nx,j,k+1)-4.*p(nx-1,j,k+1)+p(nx-2,j,k+1)
	    pxz = px(3) - px(1)
	    r(i,j,k) = r(i,j,k) - cxzxb(j,k)*pxz
	  end do
	end do
	if (nzf.eq.2) then
	  k = nz
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxb:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = 3.*p(nx,j,nz-2)-4.*p(nx-1,j,nz-2)+p(nx-2,j,nz-2)
	    px(2) = 3.*p(nx,j,nz-1)-4.*p(nx-1,j,nz-1)+p(nx-2,j,nz-1)
	    px(3) = 3.*p(nx,j,nz)-4.*p(nx-1,j,nz)+p(nx-2,j,nz)
	    pxz = 3.*px(3) - 4.*px(2) + px(1)
	    r(i,j,k) = r(i,j,k) - cxzxb(j,k)*pxz
	  end do
	end if
      end if
      if (nze.eq.2) then
	k = 1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzze:64
	do i=2,nx-1
            !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = p(i+1,j,1)-p(i-1,j,1)
	    px(2) = p(i+1,j,2)-p(i-1,j,2)
	    px(3) = p(i+1,j,3)-p(i-1,j,3)
	    pxz = -3.*px(1)+4.*px(2)-px(3)
	    r(i,j,k) = r(i,j,k) - cxzze(i,j)*pxz
	  end do
	end do
      end if
      if (nzf.eq.2) then
	k = nz
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzzf:64
	do i=2,nx-1
             !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = p(i+1,j,nz-2)-p(i-1,j,nz-2)
	    px(2) = p(i+1,j,nz-1)-p(i-1,j,nz-1)
	    px(3) = p(i+1,j,nz)-p(i-1,j,nz)
	    pxz = 3.*px(3)-4.*px(2)+px(1)
	    r(i,j,k) = r(i,j,k) - cxzzf(i,j)*pxz
	  end do
	end do
      end if
c
c     handle periodic b.c. using centered difference formula
c
      if (nxa.eq.0) then
	ip1 = 2
	im1 = nx-1
	if (nze.eq.2) then
	  k = 1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
           !dir$ assume_aligned cxzxb:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = p(ip1,j,1)-p(im1,j,1)
	    px(2) = p(ip1,j,2)-p(im1,j,2)
	    px(3) = p(ip1,j,3)-p(im1,j,3)
	    pxz = -3.*px(1) + 4.*px(2)- px(3)
	    r(1,j,1) = r(1,j,1) - cxzxa(j,k)*pxz
	    r(nx,j,1) = r(nx,j,1) - cxzxb(j,k)*pxz
	  end do
	end if
          !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
           !dir$ assume_aligned cxzxb:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	do k=2,nz-1
	  do j=js,jf
	    pxz = p(ip1,j,k+1)+p(im1,j,k-1)-(p(ip1,j,k-1)+p(im1,j,k+1))
	    r(1,j,k) = r(1,j,k) - cxzxa(j,k)*pxz
	    r(nx,j,k) = r(nx,j,k) - cxzxb(j,k)*pxz
	  end do
	end do
	if (nzf.eq.2) then
	  k = nz
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
           !dir$ assume_aligned cxzxb:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    px(1) = p(ip1,j,nz-2)-p(im1,j,nz-2)
	    px(2) = p(ip1,j,nz-1)-p(im1,j,nz-1)
	    px(3) = p(ip1,j,nz)-p(im1,j,nz)
	    pxz = 3.*px(3) - 4.*px(2) + px(1)
	    r(1,j,nz) = r(1,j,nz) - cxzxa(j,nz)*pxz
	    r(nx,j,nz) = r(nx,j,nz) - cxzxb(j,nz)*pxz
	  end do
	end if
c
c     double periodic
c
	if (nze.eq.0) then
	  km1 = nz-1
	  kp1 = 2
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzxa:64
           !dir$ assume_aligned cxzxb:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do j=js,jf
	    pxz = p(ip1,j,kp1)+p(im1,j,km1)-(p(ip1,j,km1)+p(im1,j,kp1))
	    r(1,j,1) = r(1,j,1) - cxzxa(j,1)*pxz
	    r(1,j,nz) = r(1,j,nz) - cxzxa(j,nz)*pxz
	    r(nx,j,1) = r(nx,j,1) - cxzxb(j,1)*pxz
	    r(nx,j,nz) = r(nx,j,nz) - cxzxb(j,nz)*pxz
	  end do
	end if
      end if
      if (nze.eq.0) then
	km1 = nz-1
	kp1 = 2
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzze:64
           !dir$ assume_aligned cxzxf:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	do i=2,nx-1
	  do j=js,jf
	    pxz = p(i+1,j,kp1)+p(i-1,j,km1)-(p(i+1,j,km1)+p(i-1,j,kp1))
	    r(i,j,1) = r(i,j,1) - cxzze(i,j)*pxz
	    r(i,j,nz) = r(i,j,nz) - cxzzf(i,j)*pxz
	  end do
	end do
      end if

      if (nyc.ne.1) then
	j = 1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzyc:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	do i=2,nx-1
	  do k=2,nz-1
	    pxz = p(i+1,j,k+1)+p(i-1,j,k-1)-(p(i+1,j,k-1)+p(i-1,j,k+1))
	    r(i,j,k) = r(i,j,k) - cxzyc(i,k)*pxz
	  end do
	end do
      end if
      if (nyd.ne.1) then
	j = ny
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cxzyd:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	do i=2,nx-1
	  do k=2,nz-1
	    pxz = p(i+1,j,k+1)+p(i-1,j,k-1)-(p(i+1,j,k-1)+p(i-1,j,k+1))
	    r(i,j,k) = r(i,j,k) - cxzyd(i,k)*pxz
	  end do
	end do
      end if
      return
      end

      subroutine csubyz(nx,ny,nz,p,r,cyzxa,cyzxb,cyzyc,cyzyd,cyzze,cyzzf)
          !dir$ attributes code_align : 32 :: csubyz
          !dir$ optimize : 3
          !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: csubyz
c
c     adjust right hand side r by subtracting off second order finite
c     difference approximations to pyz over nonspecified boundaries
c
      implicit none
      integer nx,ny,nz,i,j,k,ist,ifn,jm1,jp1,km1,kp1
      complex cyzxa(ny,nz),cyzxb(ny,nz)
      complex cyzyc(nx,nz),cyzyd(nx,nz)
      complex cyzze(nx,ny),cyzzf(nx,ny)
      complex p(nx,ny,nz),r(nx,ny,nz)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/icd3cr/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      !dir$ attributes align : 64 :: /icd3cr/
      complex py(3),pyz
      ist = 1
      ifn = nx
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.2) then
	j = 1
	if (nze.eq.2) then
	  k = 1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = -3.*p(i,1,1)+4.*p(i,2,1)-p(i,3,1)
	    py(2) = -3.*p(i,1,2)+4.*p(i,2,2)-p(i,3,2)
	    py(3) = -3.*p(i,1,3)+4.*p(i,2,3)-p(i,3,3)
	    pyz = -3.*py(1) + 4.*py(2)- py(3)
	    r(i,j,k) = r(i,j,k) - cyzyc(i,k)*pyz
	  end do
	end if
       
	do k=2,nz-1
           !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = -3.*p(i,1,k-1)+4.*p(i,2,k-1)-p(i,3,k-1)
	    py(3) = -3.*p(i,1,k+1)+4.*p(i,2,k+1)-p(i,3,k+1)
	    pyz = py(3) - py(1)
	    r(i,j,k) = r(i,j,k) - cyzyc(i,k)*pyz
	  end do
	end do
	if (nzf.eq.2) then
	  k = nz
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = -3.*p(i,1,nz-2)+4.*p(i,2,nz-2)-p(i,3,nz-2)
	    py(2) = -3.*p(i,1,nz-1)+4.*p(i,2,nz-1)-p(i,3,nz-1)
	    py(3) = -3.*p(i,1,nz)+4.*p(i,2,nz)-p(i,3,nz)
	    pyz = 3.*py(3) - 4.*py(2) + py(1)
	    r(i,j,k) = r(i,j,k) - cyzyc(i,k)*pyz
	  end do
	end if
      end if
      if (nyd.eq.2) then
	j = ny
	if (nze.eq.2) then
	  k = 1
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyd:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = 3.*p(i,ny,1)-4.*p(i,ny-1,1)+p(i,ny-2,1)
	    py(2) = 3.*p(i,ny,2)-4.*p(i,ny-1,2)+p(i,ny-2,2)
	    py(3) = 3.*p(i,ny,3)-4.*p(i,ny-1,3)+p(i,ny-2,3)
	    pyz = -3.*py(1) + 4.*py(2)- py(3)
	    r(i,j,k) = r(i,j,k) - cyzyd(i,k)*pyz
	  end do
	end if
	do k=2,nz-1
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyd:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = 3.*p(i,ny,k-1)-4.*p(i,ny-1,k-1)+p(i,ny-2,k-1)
	    py(3) = 3.*p(i,ny,k+1)-4.*p(i,ny-1,k+1)+p(i,ny-2,k+1)
	    pyz = py(3) - py(1)
	    r(i,j,k) = r(i,j,k) - cyzyd(i,k)*pyz
	  end do
	end do
	if (nzf.eq.2) then
	  k = nz
             !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyd:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = 3.*p(i,ny,nz-2)-4.*p(i,ny-1,nz-2)+p(i,ny-2,nz-2)
	    py(2) = 3.*p(i,ny,nz-1)-4.*p(i,ny-1,nz-1)+p(i,ny-2,nz-1)
	    py(3) = 3.*p(i,ny,nz)-4.*p(i,ny-1,nz)+p(i,ny-2,nz)
	    pyz = 3.*py(3) - 4.*py(2) + py(1)
	    r(i,j,k) = r(i,j,k) - cyzyd(i,k)*pyz
	  end do
	end if
      end if
      if (nze.eq.2) then
	k = 1
	do j=2,ny-1
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzze:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = p(i,j+1,1)-p(i,j-1,1)
	    py(2) = p(i,j+1,2)-p(i,j-1,2)
	    py(3) = p(i,j+1,3)-p(i,j-1,3)
	    pyz = -3.*py(1)+4.*py(2)-py(3)
	    r(i,j,k) = r(i,j,k) - cyzze(i,j)*pyz
	  end do
	end do
      end if
      if (nzf.eq.2) then
	k = nz
	do j=2,ny-1
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzzf:64
           
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = p(i,j+1,nz-2)-p(i,j-1,nz-2)
	    py(2) = p(i,j+1,nz-1)-p(i,j-1,nz-1)
	    py(3) = p(i,j+1,nz)-p(i,j-1,nz)
	    pyz = 3.*py(3)-4.*py(2)+py(1)
	    r(i,j,k) = r(i,j,k) - cyzzf(i,j)*pyz
	  end do
	end do
      end if
c
c     handle periodic b.c. using centered difference formula
c
      if (nyc.eq.0) then
	jp1 = 2
	jm1 = ny-1
	if (nze.eq.2) then
	  k = 1
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           !dir$ assume_aligned cyzyd:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = p(i,jp1,1)-p(i,jm1,1)
	    py(2) = p(i,jp1,2)-p(i,jm1,2)
	    py(3) = p(i,jp1,3)-p(i,jm1,3)
	    pyz = -3.*py(1) + 4.*py(2)- py(3)
	    r(i,1,k) = r(i,1,k) - cyzyc(i,k)*pyz
	    r(i,ny,k) = r(i,ny,k) - cyzyd(i,k)*pyz
	  end do
	end if
	do k=2,nz-1
            !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           !dir$ assume_aligned cyzyd:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    pyz = p(i,jp1,k+1)+p(i,jm1,k-1)-(p(i,jp1,k-1)+p(i,jm1,k+1))
	    r(i,1,k) = r(i,1,k) - cyzyc(i,k)*pyz
	    r(i,ny,k) = r(i,ny,k) - cyzyd(i,k)*pyz
	  end do
	end do
	if (nzf.eq.2) then
	  k = nz
              !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           !dir$ assume_aligned cyzyd:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    py(1) = p(i,jp1,nz-2)-p(i,jm1,nz-2)
	    py(2) = p(i,jp1,nz-1)-p(i,jm1,nz-1)
	    py(3) = p(i,jp1,nz)-p(i,jm1,nz)
	    pyz = 3.*py(3) - 4.*py(2) + py(1)
	    r(i,1,k) = r(i,1,k) - cyzyc(i,k)*pyz
	    r(i,ny,k) = r(i,ny,k) - cyzyd(i,k)*pyz
	  end do
	end if
c
c     double periodic
c
	if (nze.eq.0) then
	  km1 = nz-1
	  kp1 = 2
             !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzyc:64
           !dir$ assume_aligned cyzyd:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    pyz = p(i,jp1,kp1)+p(i,jm1,km1)-(p(i,jp1,km1)+p(i,jm1,kp1))
	    r(i,1,1) = r(i,1,1) - cyzyc(i,1)*pyz
	    r(i,1,nz) = r(i,1,nz) - cyzyc(i,nz)*pyz
	    r(i,ny,1) = r(i,ny,1) - cyzyd(i,1)*pyz
	    r(i,ny,nz) = r(i,ny,nz) - cyzyd(i,nz)*pyz
	  end do
	end if
      end if
      if (nze.eq.0) then
	km1 = nz-1
	kp1 = 2
	do j=2,ny-1
             !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzze:64
           !dir$ assume_aligned cyzzf:64
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do i=ist,ifn
	    pyz = p(i,j+1,kp1)+p(i,j-1,km1)-(p(i,j+1,km1)+p(i,j-1,kp1))
	    r(i,j,1) = r(i,j,1) - cyzze(i,j)*pyz
	    r(i,j,nz) = r(i,j,nz) - cyzzf(i,j)*pyz
	  end do
	end do
      end if
      if (nxa.ne.1) then
	i = 1
	do j=2,ny-1
             !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzxa:64
         
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do k=2,nz-1
	    pyz = p(i,j+1,k+1)+p(i,j-1,k-1)-(p(i,j+1,k-1)+p(i,j-1,k+1))
	    r(i,j,k) = r(i,j,k) - cyzxa(j,k)*pyz
	  end do
	end do
      end if
      if (nxb.ne.1) then
	i = nx
	do j=2,ny-1
             !dir$ assume_aligned p:64
           !dir$ assume_aligned r:64
           !dir$ assume_aligned cyzxb:64
         
           !dir$ ivdep
           !dir$ vector aligned
           !dir$ vector always
	  do k=2,nz-1
	    pyz = p(i,j+1,k+1)+p(i,j-1,k-1)-(p(i,j+1,k-1)+p(i,j-1,k+1))
	    r(i,j,k) = r(i,j,k) - cyzxb(j,k)*pyz
	  end do
	end do
      end if
      return
      end

      subroutine slxcd3cr(nx,ny,nz,phi,rhs,cof,tx,sum,coxy,coxz,coyz,
     +                    nxa,nyc,nze)
          !dir$ attributes code_align : 32 :: slxcd3cr
          !dir$ optimize : 3
          !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: slxcd3cr
c
c     x line relaxation thru red and then black points in the
c     (y,z) plane for periodic or nonperiodic x b.c.
c
      use omp_lib
      implicit none
      integer nx,ny,nz,i,ib,j,k
      integer nxa,nyc,nze,nper
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,7),tx(nx,ny,nz,*)
      complex rhs(nx,ny,nz),sum(ny,nz)
      complex coxy(nx,ny,nz),coxz(nx,ny,nz),coyz(nx,ny,nz)
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
     !dir$ attributes align : 64 :: /kcrsxyz/
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nxa.ne.0) then
c
c     x direction not periodic
c     first solve for x lines thru red points in (y,z) plane
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
c
	do k=1,nz,2
	  do j=1,ny,2
	    do i=1,nx
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
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
	end do
c
c     black lines in z plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do j=2,ny,2
	    do i=1,nx
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do i=2,nx
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
c
c      backward
c
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
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do j=2,ny,2
	    do i=1,nx
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do i=2,nx
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
c
c     backward sweep
c
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
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do j=1,ny,2
	    do i=1,nx
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do i=2,nx
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
c
c      backward
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
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     x direction periodic
c
	do k=1,nz
	  do j=1,ny
	    sum(j,k) = (0.0,0.0)
	  end do
	end do
c
c      sweep x lines thru red (y,z) forward and back
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do j=1,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  if (kxy.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	     do j=3,ny-1,2
	       do i=2,nx-1
		 phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +           phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	       end do
	     end do
	  end if
c
c     forward sweep
c
	  do i=2,nx-2
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
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
     +                        tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do j=2,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  if (kxy.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do i=2,nx-2
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
	    do j=2,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
c
c     backward sweep
c
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) =(phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                       tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
!$OMP PARALLEL DO SCHEDULE(static,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do j=2,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=2,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do i=2,nx-2
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
	    do j=2,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
c
c     backward sweep
c
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) =(phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                       tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tx,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do j=1,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do j=3,ny-1,2
	      do i=2,nx-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do i=2,nx-2
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
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
	    phi(nx-2,j,k) =(phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
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
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      end

      subroutine slycd3cr(nx,ny,nz,phi,rhs,cof,ty,sum,coxy,coxz,coyz,
     +                    nxa,nyc,nze)
         !dir$ attributes code_align : 32 :: slycd3cr
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: slycd3cr
c
c     x line relaxation thru red and then black points in the
c     (y,z) plane for periodic or nonperiodic x b.c.
c
      use omp_lib
      implicit none
      integer nx,ny,nz,i,jb,j,k
      integer nxa,nyc,nze,nper
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,7),ty(ny,nx,nz,*)
      complex rhs(nx,ny,nz),sum(nx,nz)
      complex coxy(nx,ny,nz),coxz(nx,ny,nz),coyz(nx,ny,nz)
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
      !dir$ attributes align : 64 :: /kcrsxyz/
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nyc.ne.0) then
c
c     y direction not periodic
c     first solve for y lines thru red points in (x,z) plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do i=1,nx,2
	    do j=1,ny
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do j=2,ny
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
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	  do i=2,nx,2
	    do j=1,ny
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do j=2,ny
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
c
c      backward
c
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
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do i=2,nx,2
	    do j=1,ny
	      phi(i,j,k) = rhs(i,j,k)- (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do j=2,ny
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
c
c     backward sweep
c
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
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do i=1,nx,2
	    do j=1,ny
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do j=2,ny
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
c
c      backward
c
	  do i=1,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     y direction periodic
c
	do k=1,nz
	  do i=1,nx
	    sum(i,k) = (0.0,0.0)
	  end do
	end do
c
c      sweep y lines thru red (x,z) forward and back
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do i=1,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do j=2,ny-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
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
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do i=2,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do j=2,ny-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
	    do i=2,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
c
c     backward sweep
c
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))/
     +                     ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=1,nz,2
	  do i=2,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do j=2,ny-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
	    do i=2,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
c
c     backward sweep
c
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,ty,kxy,kxz,kyz,nx,ny,nz)
	do k=2,nz,2
	  do i=1,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do j=2,ny-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do j=2,ny-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
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
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))/
     +                     ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      return
      end

      subroutine slzcd3cr(nx,ny,nz,phi,rhs,cof,tz,sum,coxy,coxz,coyz,
     +                    nxa,nyc,nze)
        !dir$ attributes code_align : 32 :: slzcd3cr
        !dir$ optimize : 3
        !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: slzcd3cr
c
c     x line relaxation thru red and then black points in the
c     (y,z) plane for periodic or nonperiodic x b.c.
c
      use omp_lib
      implicit none
      integer nx,ny,nz,i,kb,j,k
      integer nxa,nyc,nze,nper
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,7),tz(nz,nx,ny,*)
      complex rhs(nx,ny,nz),sum(nx,ny)
      complex coxy(nx,ny,nz),coxz(nx,ny,nz),coyz(nx,ny,nz)
      integer kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +        kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +        kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +        kxy,kxz,kyz
      common/kcrsxyz/kxyxa,kxyxb,kxyyc,kxyyd,kxyze,kxyzf,
     +               kxzxa,kxzxb,kxzyc,kxzyd,kxzze,kxzzf,
     +               kyzxa,kyzxb,kyzyc,kyzyd,kyzze,kyzzf,
     +               kxy,kxz,kyz
     !dir$ attributes align : 64 :: /kcrsxyz/
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nze.ne.0) then
c
c     z direction not periodic
c     first solve for z lines thru red points in (x,y) plane
c
!$omp parallel do schedule(static,8) private(i,j,k)
!$omp& shared(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=1,ny,2
	  do i=1,nx,2
	    do k=1,nz
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do k=2,nz
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do

!$omp parallel do schedule(static,8) private(i,j,k)
!$omp& shared(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=2,ny,2
	  do i=2,nx,2
	    do k=1,nz
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do k=2,nz
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c      backward
c
	  do i=2,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
!$omp parallel do schedule(static,8) private(i,j,k)
!$omp& shared(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=1,ny,2
	  do i=2,nx,2
	    do k=1,nz
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do k=2,nz
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c     backward sweep
c
	  do i=2,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
!$omp parallel do schedule(static,8) private(i,j,k)
!$omp& shared(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=2,ny,2
	  do i=1,nx,2
	    do k=1,nz
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward
c
	  do k=2,nz
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c      backward
c
	  do i=1,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     z direction periodic
c
	do j=1,ny
	  do i=1,nx
	    sum(i,j) = (0.0,0.0)
	  end do
	end do
c
c      sweep z lines thru red (x,y) forward and back
c
!$omp parallel do schedule(static,8) private(i,j,k)
!$omp& shared(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=1,ny,2
	  do i=1,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do k=2,nz-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=1,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
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
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
!$omp PARALLEL do schedule(static,8) private(i,j,k)
!$omp& shared(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=2,ny,2
	  do i=2,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do k=2,nz-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=2,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
c
c     backward sweep
c
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=1,ny,2
	  do i=2,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=2,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do k=2,nz-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=2,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
c
c     backward sweep
c
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k)
!$OMP& SHARED(phi,cof,rhs,coxy,coxz,coyz,tz,kxy,kxz,kyz,nx,ny,nz)
	do j=2,ny,2
	  do i=1,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = rhs(i,j,k) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     adjust for cross derivatives on interior
c
	  if (kxy.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxy(i,j,k)*(phi(i+1,j+1,k)+
     +          phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))
	      end do
	    end do
	  end if
	  if (kxz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coxz(i,j,k)*(phi(i+1,j,k+1)+
     +          phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))
	      end do
	    end do
	  end if
	  if (kyz.eq.1) then
	    do i=3,nx-1,2
	      do k=2,nz-1
		phi(i,j,k) = phi(i,j,k) - coyz(i,j,k)*(phi(i,j+1,k+1)+
     +          phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))
	      end do
	    end do
	  end if
c
c     forward sweep
c
	  do k=2,nz-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=1,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
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
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      return
      end
