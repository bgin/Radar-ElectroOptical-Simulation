c
c     file cud24cr.f
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
c     cud24cr attempts to improve the second order approximation generated
c     by cud2cr to a fourth order approximation using difference corrections
c     within multigrid cycling
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cud2cr.f, cudcom.f
c
c
      subroutine cud24cr(work,coef,bndyc,phi,ierror)
      implicit none
      complex work(*),phi(*)
      integer ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      common/cud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      external coef,bndyc
      integer nx,ny
      ierror = 0
      nx = nfx
      ny = nfy
      if (min0(nx,ny).lt.6) then
	ierror = 30
	return
      end if
      call cu24cr(nx,ny,phi,coef,bndyc,work)
      end

      subroutine cu24cr(nx,ny,phi,coef,bndyc,wk)
      implicit none
      integer nx,ny,ipf,irf,i,j,ij
      complex phi(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      common/cud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      integer icf,kb,k,ncx,ncy,ipc,irc,jc,ic,icjc,jnx
      external coef,bndyc
c
c     estimate truncation error at finest grid level
c
      icf = kcbgn(ngrid)
      irf = icf + 9*nxk(ngrid)*nyk(ngrid)
      call ctr2cr(nx,ny,phi,wk(icf),wk(irf),coef,bndyc)
c
c     set phif in wk(ipf) to zero
c
      ipf = kpbgn(ngrid)
      do j=0,ny+1
	do i=0,nx+1
	  ij = j*(nx+2)+i
	  wk(ipf+ij) = 0.0
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
	irf = kcbgn(k+1)+9*nx*ny
	irc = kcbgn(k)+9*ncx*ncy
c
c     set phic in wk(ipc) to zero
c
	do jc=0,ncy+1
	  do ic=0,ncx+1
	    icjc = jc*(ncx+2)+ic
	    wk(ipc+icjc) = 0.0
	  end do
	end do
c
c      full weight passing of truncation error from k+1 to k
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
	call kcycd2cr(wk)
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
      call kcycd2cr(wk)
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

      subroutine ctr2cr(nx,ny,phi,cof,frhs,coef,bndyc)
c
c     estimate truncation error using second order approximation in phi
c
      implicit none
      integer nx,ny
      complex phi(nx,ny),cof(nx,ny,10),frhs(nx,ny)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      common/icud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fcud2cr/xa,xb,yc,yd,tolmax,relmax
      real dlx,dly,dyox,dxoy,dlx2,dly2,dlxx,dlxy,dlyy,dlxy2,
     +             dlxy4,dxxxy4,dxyyy4,dxxyy,tdlx3,tdly3,dlx4,dly4,
     +             dlxxx,dlyyy
      common/com2dcr/dyox,dxoy,dlx2,dly2,dlxy,dlxy2,dlxy4,
     +               dxxxy4,dxyyy4,dxxyy,dlxxx,dlyyy
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      real x,y
      complex cxx,cxy,cyy,cx,cy,ce
      complex alfaa,alfab,alfac,alfad,betaa,betab,betac,betad
      complex gamaa,gamab,gamac,gamad
      complex alfim1,alfi,alfip1,betim1,beti,betip1,gamim1,gami,gamip1
      complex alfjm1,alfj,alfjp1,betjm1,betj,betjp1,gamjm1,gamj,gamjp1
      complex gbdim1,gbdi,gbdip1,gbdj,gbdjm1,gbdjp1
      complex gbdya,gbdyb,gbdyc,gbdyd
      integer i,j,isrt,ifnl,jsrt,jfnl,ii,jj,kbdy
      complex txx,txy,tyy,tx,ty,tim1,ti,tip1,tjm1,tj,tjp1,ta,tb,tc,td
      complex c1,c2,c3,c4,c5,c6,c7,c8
      complex px4,px3y,px2y2,pxy3,py4,px3,py3
      complex px3im1,py3jm1,px3jm1,py3im1

      external coef,bndyc
c
c    preset truncation estimate to zero
c
      do j=1,ny
	do  i=1,nx
	  frhs(i,j) = 0.0
	end do
      end do
c
c     set increment terms
c
      dlx=(xb-xa)/(nx-1)
      dly=(yd-yc)/(ny-1)
      dyox=dly/dlx
      dxoy=dlx/dly
      dlx2=dlx+dlx
      dly2=dly+dly
      dlxx=dlx*dlx
      dlxy=dlx*dly
      dlyy=dly*dly
      dlxy2=dlxy+dlxy
      dlxy4=dlxy2+dlxy2
      dxxxy4=4.0*dlx**3*dly
      dxyyy4=4.0*dlx*dly**3
      dxxyy=dlxx*dlyy
      tdlx3=2.0*dlx**3
      tdly3=2.0*dly**3
      dlx4=dlx**4
      dly4=dly**4
      dlxxx=dlx**3
      dlyyy=dly**3
c     set subscript limits
      if (nyc.ne.1) then
	jsrt=1
      else
	jsrt=2
      end if
      if (nxb.ne.1) then
	ifnl=nx
      else
	ifnl=nx-1
      end if
      if (nyd.ne.1) then
	jfnl=ny
      else
	jfnl=ny-1
      end if
      if (nxa.ne.1) then
	isrt=1
      else
	isrt=2
      end if
c
c     set truncation on deep interior
c
      call ctr2crd(nx,ny,phi,cof,frhs)
c
c     set truncation at and near boundaries
c
      do i=isrt,ifnl
	x = xa+(i-1)*dlx
	do j=jsrt,2
	  y = yc+(j-1)*dly
	  call coef(x,y,cxx,cxy,cyy,cx,cy,ce)
	  ii = i
	  jj = j
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txx=cxx*dlxx*px4/12.0
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  tyy=cyy*dlyy*py4/12.0
	  tx=cx*dlxx*px3/6.0
	  ty=cy*dlyy*py3/6.0
	  frhs(i,j) = txx+txy+tyy+tx+ty
	end do
	do j=ny-1,jfnl
	  y = yc+(j-1)*dly
	  call coef(x,y,cxx,cxy,cyy,cx,cy,ce)
	  ii = i
	  jj = j
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txx=cxx*dlxx*px4/12.0
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  tyy=cyy*dlyy*py4/12.0
	  tx=cx*dlxx*px3/6.0
	  ty=cy*dlyy*py3/6.0
	  frhs(i,j) = txx+txy+tyy+tx+ty
	end do
      end do
      do j=3,ny-2
	y = yc+(j-1)*dly
	do i=isrt,2
	  x = xa+(i-1)*dlx
	  call coef(x,y,cxx,cxy,cyy,cx,cy,ce)
	  ii = i
	  jj = j
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txx=cxx*dlxx*px4/12.0
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  tyy=cyy*dlyy*py4/12.0
	  tx=cx*dlxx*px3/6.0
	  ty=cy*dlyy*py3/6.0
	  frhs(i,j) = txx+txy+tyy+tx+ty
	end do
	do i=nx-1,ifnl
	   x = xa+(i-1)*dlx
	   call coef(x,y,cxx,cxy,cyy,cx,cy,ce)
	   ii = i
	   jj = j
	   call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	   call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	   txx=cxx*dlxx*px4/12.0
	   txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	   tyy=cyy*dlyy*py4/12.0
	   tx=cx*dlxx*px3/6.0
	   ty=cy*dlyy*py3/6.0
	   frhs(i,j) = txx+txy+tyy+tx+ty
	end do
      end do
c
c     adjust along boundaries
c
      if (nyc.eq.2) then
	kbdy=3
	call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
	call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
	ii = 1
	jj = 1
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	px3im1=px3
	py3im1=py3
	ii = 2
	jj = 1
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	do i=2,nx-1
	  x=xa+float(i-1)*dlx
	  call bndyc(kbdy,x+dlx,alfip1,betip1,gamip1,gbdip1)
	  call coef(x,yc,cxx,cxy,cyy,cx,cy,ce)
	  c6=cxy/dlxy4
	  c7=cyy/dlyy-cy/dly2
	  c8=-c6
	  tim1=-alfim1*dlxx*px3im1/3.0+betim1*dlyy*py3im1/6.0
	  px3im1=px3
	  py3im1=py3
	  ti=alfi*dlxx*px3/6.0+beti*dlyy*py3/6.0
	  ii = i+1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tip1=-alfip1*dlxx*px3/3.0+betip1*dlyy*py3/6.0
c     add adjustment to standard truncation error
	  frhs(i,1)=frhs(i,1)+dly2*(c6/betim1*tim1+c7/beti*ti+c8/betip1*
     +    tip1)
	  alfim1=alfi
	  alfi=alfip1
	  betim1=beti
	  beti=betip1
	end do
	if (nxa.eq.0) then
c     periodic in x direction
	  kbdy = 3
	  call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
	  call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
	  call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
	  c8=-cxy/dlxy2
	  c7=cyy/dlyy-cy/dly2-c8
	  ii = 1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
c     adjust truncation error from cross derivative at cornor
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(1,1)=frhs(1,1)-txy
	  txy=cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(1,1)=frhs(1,1)+txy
	  ti=alfi*dlxx*px3/6.0+beti*dlyy*py3/6.0
	  ii = 2
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tip1=-alfip1*dlxx*px3/3.0+betip1*dlyy*py3/6.0
	  frhs(1,1)=frhs(1,1)+dly2*(c7/beti*ti+c8/betip1*tip1)
c     at (xb,yc)
	  kbdy = 3
	  call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
	  call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
	  call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
	  ii = nx
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(nx,1)=frhs(nx,1)-txy
	  txy=-cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(nx,1)=frhs(nx,1)+txy
	  ti=alfi*dlxx*px3/6.0+beti*dlyy*py3/6.0
	  ii = nx-1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tim1=-alfim1*dlxx*px3im1/3.0+betim1*dlyy*py3im1/6.0
	  c6=cxy/dlxy2
	  c7=cyy/dlyy-cy/dly2-c6
	  frhs(nx,1)=frhs(nx,1)+dly2*(c6/betim1*tim1+c7/beti*ti)
	end if
	if (nxa.eq.2) then
c     mixed-mixed at (xa,yc)
	  call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
	  c8=-cxy/dlxy2
	  c4=c8
	  kbdy=3
	  call bndyc(kbdy,xa+dlx,alfac,betac,gamac,gbdyc)
	  ii = 2
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tip1=-alfac*dlxx*px3/3.0+betac*dlyy*py3/6.0
	  kbdy=1
	  call bndyc(kbdy,yc+dly,alfaa,betaa,gamaa,gbdya)
	  ii = 1
	  jj = 2
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjp1=alfaa*dlxx*px3/6.0-betaa*dlyy*py3/3.0
	  frhs(1,1)=frhs(1,1)+dly2*c8/betac*tip1+dlx2*c4/alfaa*tjp1
c     adjust cross derivative truncation error at cornor
	  ii = 1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(1,1)=frhs(1,1)-txy
	  txy=cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*py3/6.0)
	  frhs(1,1)=frhs(1,1)+txy
c     phase 2
	  c5=cxx/dlxx-cx/dlx2-c4+c8*(alfac/betac*dyox)
	  c7=cyy/dlyy-cy/dly2-c4+c4*(betaa/alfaa*dxoy)
	  ii = 1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  ta =-dlyyy*py3/3.0
	  tc =-dlxxx*px3/3.0
	  frhs(1,1)=frhs(1,1)-c5*ta-c7*tc
	end if
	if (nxb.eq.2) then
c     mixed-mixed at (xb,yc)
	  call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
	  c2=cxy/dlxy2
	  c6=c2
	  kbdy=3
	  call bndyc(kbdy,xb-dlx,alfac,betac,gamac,gbdyc)
	  ii = nx-1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tim1=-alfac*dlxx*px3/3.0+betac*dlyy*py3/6.0
	  kbdy=2
	  call bndyc(kbdy,yc+dly,alfab,betab,gamab,gbdyb)
	  ii = nx
	  jj = 2
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjp1=alfab*dlxx*px3/6.0-betab*dlyy*py3/3.0
	  frhs(nx,1)=frhs(nx,1)+dly2*c6/betac*tim1-dlx2*c2/alfab*tjp1
c     phase 2
	  ii = nx
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(nx,1)=frhs(nx,1)-txy
	  txy=-cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(nx,1)=frhs(nx,1)+txy
	  kbdy=3
	  call bndyc(kbdy,xb,alfac,betac,gamac,gbdyc)
	  kbdy=2
	  call bndyc(kbdy,yc,alfab,betab,gamab,gbdyb)
	  c1=cxx/dlxx+cx/dlx2+c6*(-alfac/betac*dyox)-c2
	  c7=cyy/dlyy-cy/dly2-c2*(betab/alfab*dxoy)-c2
	  tb=dlyyy*py3/3.0
	  tc=-dlxxx*px3/3.0
	  frhs(nx,1)=frhs(nx,1)-c1*tb-c7*tc
	end if
      end if
      if (nxb.eq.2) then
	kbdy=2
	call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
	ii = nx
	jj = 1
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	px3jm1=px3
	py3jm1=py3
	ii = nx
	jj = 2
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	do j=2,ny-1
	  y=yc+float(j-1)*dly
	  call bndyc(kbdy,y+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	  call coef(xb,y,cxx,cxy,cyy,cx,cy,ce)
	  c2=cxy/dlxy4
	  c1=cxx/dlxx+cx/dlx2
	  c8=-c2
	  tjm1=alfjm1*dlxx*px3jm1/6.0-betjm1*dlyy*py3jm1/3.0
	  px3jm1=px3
	  py3jm1=py3
	  tj=(alfj*dlxx*px3+betj*dlyy*py3)/6.0
	  ii = nx
	  jj = j+1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjp1=alfjp1*dlxx*px3/6.0-betjp1*dlyy*py3/3.0
	  frhs(nx,j)=frhs(nx,j)-dlx2*(c8/alfjm1*tjm1+c1/alfj*tj+c2/
     +    alfjp1*tj)
	  alfjm1=alfj
	  betjm1=betj
	  alfj=alfjp1
	  betj=betjp1
	end do
	if (nyc.eq.0) then
c
c     need code for y periodic along x = xb
c
	  kbdy=2
	  call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
	  call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	  call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
	  c2=cxy/dlxy2
	  c3=cxx/dlxx+cx/dlx2-c2
c     adjust truncation error from cross derivative at cornor
	  ii = nx
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(nx,1)=frhs(nx,1)-txy
	  txy=-cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(nx,1)=frhs(nx,1)+txy
	  tj=alfj*dlxx*px3/6.0+betj*dlyy*py3/6.0
	  ii = nx
	  jj = 2
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjp1=alfjp1*dlxx*px3/3.0-betjp1*dlyy*py3/6.0
	  frhs(nx,1)=frhs(nx,1)-dlx2*(c3/alfj*tj+c2/alfjp1*tjp1)
c     at (xb,yd)
	  call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
	  call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
	  call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
	  ii = nx
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(nx,ny)=frhs(nx,ny)-txy
	  txy=cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(nx,ny)=frhs(nx,ny)+txy
	  tj=alfj*dlxx*px3/6.0+betj*dlyy*py3/6.0
	  ii = nx
	  jj = ny-1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjm1=-alfjm1*dlxx*px3/3.0+betjm1*dlyy*py3/6.0
	  c4=-cxy/dlxy2
	  c3=cxx/dlxx+cx/dlx2-c4
	  frhs(nx,ny)=frhs(nx,ny)-dlx2*(c4/alfjm1*tjm1+c3/alfj*tj)
	end if
	if (nyd.eq.2) then
c     mixed-mixed at (xb,yd)
	  call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
	  c4=-cxy/dlxy2
	  c8=c4
	  kbdy=4
	  call bndyc(kbdy,xb-dlx,alfad,betad,gamad,gbdyd)
	  kbdy=2
	  call bndyc(kbdy,yd-dly,alfab,betab,gamab,gbdyb)
	  ii = nx-1
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tim1=alfad*dlxx*px3/3.0+betad*dlyy*py3/6.0
	  ii = nx
	  jj = ny-1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjm1=alfab*dlxx*px3/6.0+betab*dlyy*py3/3.0
	  frhs(nx,ny)=frhs(nx,ny)-dly2*c4/betad*tim1-dlx2*c8/alfab*tjp1
c     phase 2
	  ii = nx
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(nx,ny)=frhs(nx,ny)-txy
	  txy=cxy*(dlxx*px3/6.0-dlxy*px2y2/4.0+dlyy*py3/6.0)
	  frhs(nx,ny)=frhs(nx,ny)+txy
	  kbdy=2
	  call bndyc(kbdy,yd,alfab,betab,gamab,gbdyb)
	  kbdy=4
	  call bndyc(kbdy,xb,alfad,betad,gamad,gbdyd)
	  c1=cxx/dlxx+cx/dlx2-c4+c4*(alfad/betad*dyox)
	  c3=cyy/dlyy+cy/dly2-c4+c8*(betab/alfab*dxoy)
	  tb=dlyyy*py3/3.0
	  td=dlxxx*px3/3.0
	  frhs(nx,ny)=frhs(nx,ny)-c1*tb-c3*td
	end if
      end if

      if (nyd.eq.2) then
c     adjust at y=yd when mixed
	kbdy=4
	call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
	call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
	ii = 1
	jj = ny
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	px3im1=px3
	py3im1=py3
	ii = 2
	jj = ny
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	do i=2,nx-1
	  x=xa+float(i-1)*dlx
	  call bndyc(kbdy,x+dlx,alfip1,betip1,gamip1,gbdip1)
	  call coef(x,yd,cxx,cxy,cyy,cx,cy,ce)
	  c2=cxy/dlxy4
	  c3=cyy/dlyy+cy/dly2
	  c4=-c2
	  tim1=-alfim1*dlxx*px3im1/3.0+betim1*dlyy*py3im1/6.0
	  px3im1=px3
	  py3im1=py3
	  ti=(alfi*dlxx*px3+beti*dlyy*py3)/6.0
	  ii = i+1
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tip1=-alfip1*dlxx*px3/3.0+betip1*dlyy*py3/6.0
	  frhs(i,ny)=frhs(i,ny)-dly2*(c4/betim1*tim1+c3/beti*ti+c2/
     +    betip1*tip1)
	  alfim1=alfi
	  alfi=alfip1
	  betim1=beti
	  beti=betip1
	end do
	if (nxa.eq.0) then
c     periodic in x direction
	  kbdy=4
	  call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
	  call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
	  call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
	  c2=cxy/dlxy2
	  c3=cyy/dlyy+cy/dly2-c2
c     adjust truncation error from cross derivative at cornor
	  ii = 1
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(1,ny)=frhs(1,ny)-txy
	  txy=-cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(1,ny)=frhs(1,ny)+txy
	  ti=alfi*dlxx*px3/6.0+beti*dlyy*py3/6.0
	  ii = 2
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tip1=-alfip1*dlxx*px3/3.0+betip1*dlyy*py3/6.0
	  frhs(1,ny)=frhs(1,ny)-dly2*(c3/beti*ti+c2/betip1*tip1)
c     at (xb,yd)
	  call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
	  call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
	  call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
	  ii = nx
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(nx,ny)=frhs(nx,ny)-txy
	  txy=cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(nx,ny)=frhs(nx,ny)+txy
	  ti=alfi*dlxx*px3/6.0+beti*dlyy*py3/6.0
	  ii = nx-1
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tim1=-alfim1*dlxx*px3im1/3.0+betim1*dlyy*py3im1/6.0
	  c4=-cxy/dlxy2
	  c3=cyy/dlyy+cy/dly2-c4
	  frhs(nx,ny)=frhs(nx,ny)-dly2*(c4/betim1*tim1+c3/beti*ti)
	end if
	if (nxa.eq.2) then
	  call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
c     phase 1
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
	  ii = 2
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tip1=-alfad*dlxx*px3/3.0+betad*dlxx*dlyy*py3/6.0
	  ii = 1
	  jj = ny-1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjm1=alfaa*dlxx*px3/6.0-betaa*dlyy*py3/3.0
	  frhs(1,ny)=frhs(1,ny)-dly2*c2/betad*tip1+dlx2*c6/alfaa*tjm1
c     adjust cornor truncation error from skewed difference formula
	  ii = 1
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(1,ny)=frhs(1,ny)-txy
	  txy=-cxy*(dlxx*px3/6.0-dlxy*px2y2/4.0+dlyy*py3/6.0)
	  frhs(1,ny)=frhs(1,ny)+txy
c     phase 2
	  ta=-dlyyy*py3/3.0
	  td=dlxxx*px3/3.0
	  frhs(1,ny)=frhs(1,ny)-c5*ta-c3*td
	end if
      end if

      if (nxa.eq.2) then
c     adjust truncation error along x=xa if mixed
	kbdy=1
	call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
	call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
	ii = 1
	jj = 1
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	px3jm1=px3
	py3jm1=py3
	ii = 1
	jj = 2
	call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	do j=2,ny-1
	  y=yc+float(j-1)*dly
	  call bndyc(kbdy,y+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	  call coef(xa,y,cxx,cxy,cyy,cx,cy,ce)
	  c4=-cxy/dlxy4
	  c5=cxx/dlxx-cx/dlx2
	  c6=-c4
	  tjm1=alfjm1*dlxx*px3jm1/6.0-betjm1*dlyy*py3jm1/3.0
	  px3jm1=px3
	  py3jm1=py3
	  tj=(alfj*dlxx*px3+betj*dlyy*py3)/6.0
	  ii = 1
	  jj = j+1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjp1=alfjp1*dlxx*px3/6.0-betjp1*dlyy*py3/3.0
	  frhs(1,j)=frhs(1,j)+dlx2*(c6/alfjm1*tjm1+c5/alfj*tj+c4/
     +    alfjp1*tjp1)
	  alfjm1=alfj
	  betjm1=betj
	  alfj=alfjp1
	  betj=betjp1
	end do
	if (nyc.eq.0) then
c
c     need code for y periodic along x=xa
c
	  kbdy = 1
	  call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
	  call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
	  call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
	  c8=-cxy/dlxy2
	  c7=cxx/dlxx-cx/dlx2-c8
	  ii = 1
	  jj = 1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
c     adjust truncation error from cross derivative at cornor
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(1,1)=frhs(1,1)-txy
	  txy=cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(1,1)=frhs(1,1)+txy
	  tj=alfj*dlxx*px3/6.0+betj*dlyy*py3/6.0
	  ii = 1
	  jj = 2
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjp1=alfjp1*dlxx*px3/3.0-betjp1*dlyy*py3/6.0
	  frhs(1,1)=frhs(1,1)+dlx2*(c7/alfj*tj+c8/alfjp1*tjp1)
c     at (xa,yd)
	  kbdy = 1
	  call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
	  call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
	  call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
	  ii = 1
	  jj = ny
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  frhs(1,ny)=frhs(1,ny)-txy
	  txy=-cxy*(dlxx*px3y/6.0-dlxy*px2y2/4.0+dlyy*pxy3/6.0)
	  frhs(1,ny)=frhs(1,ny)+txy
	  tj=alfj*dlxx*px3/6.0+betj*dlyy*py3/6.0
	  ii = 1
	  jj = ny-1
	  call cpde2(nx,ny,phi,ii,jj,px3,px4,py3,py4,nxa,nyc)
	  call cpde2cr(nx,ny,phi,ii,jj,px3y,pxy3,px2y2)
	  tjm1=-alfjm1*dlxx*px3/3.0+betim1*dlyy*py3/6.0
	  c6=cxy/dlxy2
	  c7=cxx/dlxx-cx/dlx2-c6
	  frhs(1,ny)=frhs(1,ny)+dlx2*(c6/alfjm1*tjm1+c7/alfj*tj)
	end if
      end if
      return
      end

      subroutine ctr2crd(nx,ny,u,cof,frhs)
c
c     estimate truncation on deep interior
c
      implicit none
      integer nx,ny
      complex cof(nx,ny,10),u(nx,ny),frhs(nx,ny)
      real dlx,dly,dyox,dxoy,dlx2,dly2,dlxx,dlxy,dlyy,dlxy2,
     +             dlxy4,dxxxy4,dxyyy4,dxxyy,tdlx3,tdly3,dlx4,dly4,
     +             dlxxx,dlyyy
      common/com2dcr/dyox,dxoy,dlx2,dly2,dlxy,dlxy2,dlxy4,
     +               dxxxy4,dxyyy4,dxxyy,dlxxx,dlyyy
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      complex cxx,cxy,cyy,cx,cy
      integer i,j
      complex txx,txy,tyy,tx,ty
      complex px4,px3y,pxy3,py4,px3,py3
!$OMP PARALLEL SCHEDULE(STATIC,8) DO PRIVATE(i,j,px3,px4,py3,py4,px3y,pxy3) &
!$OMP PRIVATE(cxx,cx,cyy,cy,cxy,tx,ty,txx,tyy,txy) &
!$OMP SHARED (frhs,u,cof,nx,ny) &
!$OMP SHARED (dlx,dly,dlxx,dlyy,dlx4,dly4,tdlx3,tdly3,dxxxy4,dxyyy4)
      do j=3,ny-2
	do i=3,nx-2
c
c     estimate partial derivatives
c
	  px3 = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
	  px4 = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))
     +          /dlx4
	  py3 = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
	  py4 = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))
     +          /dly4
	  px3y = (u(i-2,j-1)-2*u(i-1,j-1)+2*u(i+1,j-1) -u(i+2,j-1)-
     +    u(i-2,j+1)+2*u(i-1,j+1)-2*u(i+1,j+1)+u(i+2,j+1))/dxxxy4
     +
	  pxy3 = (u(i-1,j-2)-u(i+1,j-2)-2*u(i-1,j-1)+2*u(i+1,j-1)+
     +    2*u(i-1,j+1)-2*u(i+1,j+1)-u(i-1,j+2)+u(i+1,j+2))/dxyyy4
c
c     reset pde coefficients from discretized coefs
c
	  cxx = 0.5*dlxx*(cof(i,j,1)+cof(i,j,5))
	  cxy = dlxy4*cof(i,j,2)
	  cyy = 0.5*dlyy*(cof(i,j,3)+cof(i,j,7))
	  cx = dlx*(cof(i,j,1)-cof(i,j,5))
	  cy = dly*(cof(i,j,3)-cof(i,j,7))
	  txx=cxx*dlxx*px4/12.0
	  txy=cxy*(dlxx*px3y+dlyy*pxy3)/6.0
	  tyy=cyy*dlyy*py4/12.0
	  tx=cx*dlxx*px3/6.0
	  ty=cy*dlyy*py3/6.0
	  frhs(i,j) = txx+txy+tyy+tx+ty
	end do
      end do
      return
      end

