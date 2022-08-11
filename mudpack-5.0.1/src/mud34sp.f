c
c     file mud34sp.f
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
c     mud34sp attempts to improve the second order approximation generated
c     by mud3sp to a fourth order approximation using difference corrections
c     within multigrid cycling
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mud3sp.f, mudcom.f
c
c
c
      subroutine mud34sp(work,phi,ierror)
      implicit none
      real phi(*),work(*)
      integer ierror
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      integer nx,ny,nz
      ierror = 0
      nx = nfx
      ny = nfy
      nz = nfz
      if (min0(nx,ny,nz).lt.6) then
      ierror = 30
      return
      end if
      call md34sp(nx,ny,nz,phi,work)
      end

      subroutine md34sp(nx,ny,nz,phi,wk)
      implicit none
      integer nx,ny,nz
      real phi(nx,ny,nz),wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer ipf,irf,i,j,k,ijk,kb,ncx,ncy,ncz,irc,ic,jc,kc
      integer icx,icy,icz,jn,kn,ipc
      ipf = kpbgn(ngrid)
      irf = krbgn(ngrid)
      icx = kcxbgn(ngrid)
      icy = kcybgn(ngrid)
      icz = kczbgn(ngrid)
c
c     estimate truncation error and pass down (weighted) to all grid levels
c
      call tr3sp(nx,ny,nz,phi,wk(irf),wk(icx),wk(icy),wk(icz))
c
c     set phif in wk(ipf) to zero
c
      do k=0,nz+1
	do j=0,ny+1
	  do i=0,nx+1
	    ijk = (nx+2)*(j+k*(ny+2))+i
	    wk(ipf+ijk) = 0.0
	  end do
	end do
      end do
      do kb=2,ngrid
	k = ngrid-kb+1
	nx = nxk(k+1)
	ny = nyk(k+1)
	nz = nzk(k+1)
	ncx = nxk(k)
	ncy = nyk(k)
	ncz = nzk(k)
	ipf = kpbgn(k+1)
	ipc = kpbgn(k)
	irf = krbgn(k+1)
	irc = krbgn(k)
c
c     set phic in wk(ipc) to zero
c
	do kc=0,ncz+1
	  do jc=0,ncy+1
	    do ic=0,ncx+1
	      ijk = (kc*(ncy+2)+jc)*(ncx+2)+ic
	      wk(ipc+ijk) = 0.0
	    end do
	  end do
	end do
c
c      full weighting of truncation error from k+1 to k
c
	call res3(nx,ny,nz,wk(irf),ncx,ncy,ncz,wk(irc),nxa,nxb,
     +            nyc,nyd,nze,nzf)
      end do
c
c     execute one fmw(2,1) cycle for fourth order correction term
c
      kcycle = 2
      iprer = 2
      ipost = 1
      intpol = 3
      do k=1,ngrid-1
	kcur = k
	call kcymd3sp(wk)
	ipf = kpbgn(k+1)
	nx = nxk(k+1)
	ny = nyk(k+1)
	nz = nzk(k+1)
	ipc = kpbgn(k)
	ncx = nxk(k)
	ncy = nyk(k)
	ncz = nzk(k)
	call prolon3(ncx,ncy,ncz,wk(ipc),nx,ny,nz,wk(ipf),nxa,nxb,nyc,
     +               nyd,nze,nzf,intpol)
      end do
      kcur = ngrid
      call kcymd3sp(wk)
c
c     add fourth order correction term to solution
c
      do k=1,nz
	kn = k*(ny+2)*(nx+2)
	do j=1,ny
	  jn = j*(nx+2)
	  do i=1,nx
	    ijk = kn+jn+i+1
	    phi(i,j,k) = phi(i,j,k)+wk(ijk)
	  end do
	end do
      end do
      return
      end

      subroutine tr3sp(nx,ny,nz,u,r,cofx,cofy,cofz)
c
c     estimate truncation error using second order approximation in u
c
      implicit none
      integer nx,ny,nz
      real r(nx,ny,nz),u(nx,ny,nz)
      real cofx(nx,3),cofy(ny,3),cofz(nz,3)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,dlx4,dly4,dlz4
      common/fmud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      integer i,j,k,istart,ifinal,jstart,jfinal,kstart,kfinal
      real tx,ty,tz,cxx,cyy,czz,ux3,ux4,uy3,uy4,uz3,uz4,cx,cy,cz
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlz = (zf-ze)/(nz-1)
      dlxx = dlx*dlx
      dlyy = dly*dly
      dlzz = dlz*dlz
      tdlx3 = 2.*dlxx*dlx
      tdly3 = 2.*dlyy*dly
      tdlz3 = 2.*dlzz*dlz
      dlx4 = dlxx*dlxx
      dly4 = dlyy*dlyy
      dlz4 = dlzz*dlzz
c
c     set subscript limit to avoid specified (dirchlet) boundaries.
c
      istart=1
      ifinal=nx
      jstart=1
      jfinal=ny
      kstart=1
      kfinal=nz
      if(nxa.eq.1) istart=2
      if(nxb.eq.1) ifinal=nx-1
      if(nyc.eq.1) jstart=2
      if(nyd.eq.1) jfinal=ny-1
      if (nze.eq.1) kstart=2
      if (nzf.eq.1) kfinal=nz-1
      do k=1,nz
	do j=1,ny
	  do i=1,nx
	    r(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     compute truncation on deep interior--this should vectorize
c
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,ux3,ux4,uy3,uy4,uz3,uz4) &
!$OMP PRIVATE(cxx,cx,cyy,cy,czz,cz,tx,ty,tz) &
!$OMP SHARED (r,u,cofx,cofy,cofz,nx,ny,nz,dlx,dly,dlz) &
!$OMP SHARED (dlxx,dlyy,dlzz,dlx4,dly4,dlz4,tdlx3,tdly3,tdlz3)
      do k=3,nz-2
	do j=3,ny-2
	  do i=3,nx-2
	    ux3=(-u(i-2,j,k)+2.0*u(i-1,j,k)-2.0*u(i+1,j,k)+u(i+2,j,k))
     +          /tdlx3
	    ux4=(u(i-2,j,k)-4.0*u(i-1,j,k)+6.0*u(i,j,k)-4.0*u(i+1,j,k)+
     +           u(i+2,j,k))/dlx4
	    uy3=(-u(i,j-2,k)+2.0*u(i,j-1,k)-2.0*u(i,j+1,k)+u(i,j+2,k))
     +          /tdly3
	    uy4=(u(i,j-2,k)-4.0*u(i,j-1,k)+6.0*u(i,j,k)-4.0*u(i,j+1,k)+
     +          u(i,j+2,k))/dly4
	    uz3=(-u(i,j,k-2)+2.0*u(i,j,k-1)-2.0*u(i,j,k+1)+u(i,j,k+2))
     +          /tdlz3
	    uz4=(u(i,j,k-2)-4.0*u(i,j,k-1)+6.0*u(i,j,k)-4.0*u(i,j,k+1)+
     +           u(i,j,k+2))/dlz4
c
c     reset continuous pde coefficients from discretized coefs
c
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
	    r(i,j,k) = dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
      end do
c
c     estimate truncation error at and near nonspecified boundaries
c
      do k=kstart,kfinal
	do i=istart,ifinal
	  do j=jstart,2
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
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
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cofz(k,2)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cofz(k,1)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    r(i,j,k)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	  do j=ny-1,jfinal
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
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
	    if (nyc.ne.0 .and. j.eq.ny) then
	      cyy = 0.5*dlyy*cofy(j,1)
	      ty=cyy*(uy4*0.25-uy3/dly)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cofz(k,2)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cofz(k,1)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    r(i,j,k)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
	do j=3,ny-2
	  do i=istart,2
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
c
c     adjust truncation error terms at mixed boundaries
c
	    if(nxa.ne.0 .and. i.eq.1) then
	      cxx = 0.5*dlxx*cofx(i,2)
	      tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cofz(k,2)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cofz(k,1)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    r(i,j,k)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	  do i=nx-1,ifinal
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
c
c     adjust truncation error terms at mixed boundaries
c
	    if(nxa.ne.0 .and. i.eq.nx) then
	      cxx = 0.5*dlxx*cofx(i,1)
	      tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cofz(k,2)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cofz(k,1)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    r(i,j,k)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
      end do
      do i=3,nx-2
	do j=3,ny-2
	  do k=kstart,2
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
	    if(nze.ne.0 .and. k.eq.1) then
	      czz = 0.5*dlzz*cofz(k,2)
	      tz = czz*(uz4*0.25+uz3/dlz)/3.0
	    end if
	    r(i,j,k)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	  do k=nz-1,kfinal
	    cxx = 0.5*dlxx*(cofx(i,1)+cofx(i,2))
	    cx = dlx*(cofx(i,2)-cofx(i,1))
	    cyy = 0.5*dlyy*(cofy(j,1)+cofy(j,2))
	    cy = dly*(cofy(j,2)-cofy(j,1))
	    czz = 0.5*dlzz*(cofz(k,1)+cofz(k,2))
	    cz = dlz*(cofz(k,2)-cofz(k,1))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
	    if(nze.ne.0 .and. k.eq.nz) then
	      czz = 0.5*dlzz*cofz(k,1)
	      tz = czz*(uz4*0.25-uz3/dlz)/3.0
	    end if
	    r(i,j,k)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
      end do
      return
      end


