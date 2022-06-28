c
c     file mud34.f
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
c     mud34 attempts to improve the second order approximation generated
c     by mud3 to a fourth order approximation using difference corrections
c     within multigrid cycling
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mud3.f, mudcom.f, mud3ln.f, mud3pn.f
c
c
c
      subroutine mud34(work,phi,ierror)
        !dir$ attributes code_align : 32 :: mud34
        !dir$ optimize : 3
      implicit none
      real phi(*),work(*)
      integer ierror
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
       !dir$ attributes align : 64 :: /imud3/
      integer nx,ny,nz
      ierror = 0
      nx = nfx
      ny = nfy
      nz = nfz
      if (min0(nx,ny,nz).lt.6) then
      ierror = 30
      return
      end if
      call md34(nx,ny,nz,phi,work)
      end

      subroutine md34(nx,ny,nz,phi,wk)
        !dir$ attributes forceinline :: md34
       !dir$ attributes code_align : 32 :: md34
       !dir$ optimize : 3
      implicit none
      integer nx,ny,nz
      real phi(nx,ny,nz),wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
       !dir$ attributes align : 64 :: /imud3/
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
       !dir$ attributes align : 64 :: /mud3c/
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid
      integer klevel,kcur,kps
      integer ipf,icf,i,j,k,ijk,kb,ncx,ncy,ncz,irf,irc,ic,jc,kc
      integer jn,kn,ipc
      ipf = kpbgn(ngrid)
      icf = kcbgn(ngrid)
c
c     estimate truncation error and pass down (weighted) to all grid levels
c
      call tr3(nx,ny,nz,phi,wk(icf))
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
	irf = kcbgn(k+1)+7*nx*ny*nz
	irc = kcbgn(k)+7*ncx*ncy*ncz
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
	call kcymd3(wk)
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
      call kcymd3(wk)
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

      subroutine tr3(nx,ny,nz,u,cof)
       !dir$ attributes forceinline :: tr3
       !dir$ attributes code_align : 32 :: tr3
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: tr3
       use omp_lib
c
c     estimate truncation error using second order approximation in u
c
      implicit none
      integer nx,ny,nz
      real cof(nx,ny,nz,8),u(nx,ny,nz)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
     !dir$ attributes align : 64 :: /imud3/
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,dlx4,dly4,dlz4
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
       !dir$ attributes align : 64 :: /fmud3/
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
        !dir$ attributes align : 64 :: /pde3com/
      integer i,j,k,istart,ifinal,jstart,jfinal,kstart,kfinal
      real tx,ty,tz,cxx,cyy,czz,ux3,ux4,uy3,uy4,uz3,uz4,cx,cy,cz
      do k=1,nz
	do j=1,ny
	  do i=1,nx
	    cof(i,j,k,8) = 0.0
	  end do
	end do
      end do
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
	    cof(i,j,k,8) = 0.0
	  end do
	end do
      end do
c
c     compute truncation on deep interior--this should vectorize
c
     !dir$ assume_aligned u:64,r:64
!$OMP PARALLEL DO SCHEDULE(STATIC,8) PRIVATE(i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,cxx,cx,cyy,cy) &
!$OMP PRIVATE(czz,cz,tx,ty,tz)  &
!$OMP SHARED (u,cof,nx,ny,nz)   &
!$OMP SHARED (dlx,dly,dlz,dlxx,dlyy,dlzz,dlx4,dly4,dlz4) &
!$OMP SHARED (tdlx3,tdly3,tdlz3)
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
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
	    cof(i,j,k,8) = dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
      end do
!$omp end parallel do
c
c     estimate truncation error at and near nonspecified boundaries
c
       
       !dir$ assume_aligned cof:64
       !dir$ assume_aligned u:64
!$omp parallel parallel do default(none) schedule(static,8) &
!$omp private(k,i,j,cxx,cx,cyy,czz,cz,tx,ty,tz)             &
!$omp private(ux3,ux4,uy3,uy4,uz3,uz4)                      &
!$omp shared(kstart,kfinal,istart,jstart,ifinal,dlxx,cof)    &
!$omp shared(u,dlyy,dlzz,dlz,nx,ny,nz,nxa,nyc,dlx,nze)   &
!$omp shared(r,jfinal,ifinal,dly)
      do k=kstart,kfinal
	do i=istart,ifinal
	  do j=jstart,2
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
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
		cxx = 0.5*dlxx*cof(i,j,k,2)
		tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	      else
		cxx = 0.5*dlxx*cof(i,j,k,1)
		tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	      end if
	    end if
	    if (nyc.ne.0 .and. j.eq.1) then
	      cyy = 0.5*dlyy*cof(i,j,k,4)
	      ty=cyy*(uy4*0.25+uy3/dly)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cof(i,j,k,6)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cof(i,j,k,5)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    cof(i,j,k,8)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	  do j=ny-1,jfinal
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
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
		cxx = 0.5*dlxx*cof(i,j,k,2)
		tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	      else
		cxx = 0.5*dlxx*cof(i,j,k,1)
		tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	      end if
	    end if
	    if (nyc.ne.0 .and. j.eq.ny) then
	      cyy = 0.5*dlyy*cof(i,j,k,3)
	      ty=cyy*(uy4*0.25-uy3/dly)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cof(i,j,k,6)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cof(i,j,k,5)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    cof(i,j,k,8)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
	do j=3,ny-2
	  do i=istart,2
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
c
c     adjust truncation error terms at mixed boundaries
c
	    if(nxa.ne.0 .and. i.eq.1) then
	      cxx = 0.5*dlxx*cof(i,j,k,2)
	      tx = cxx*(ux4*0.25+ux3/dlx)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cof(i,j,k,6)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cof(i,j,k,5)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    cof(i,j,k,8)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	  do i=nx-1,ifinal
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
c
c     adjust truncation error terms at mixed boundaries
c
	    if(nxa.ne.0 .and. i.eq.nx) then
	      cxx = 0.5*dlxx*cof(i,j,k,1)
	      tx = cxx*(ux4*0.25-ux3/dlx)/3.0
	    end if
	    if(nze.ne.0 .and. (k.eq.1 .or. k.eq.nz)) then
	      if (k.eq.1) then
		czz = 0.5*dlzz*cof(i,j,k,6)
		tz = czz*(uz4*0.25+uz3/dlz)/3.0
	      else
		czz = 0.5*dlzz*cof(i,j,k,5)
		tz = czz*(uz4*0.25-uz3/dlz)/3.0
	      end if
	    end if
	    cof(i,j,k,8)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
      end do
!$omp end parallel do
      do i=3,nx-2
	do j=3,ny-2
	  do k=kstart,2
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
	    if(nze.ne.0 .and. k.eq.1) then
	      czz = 0.5*dlzz*cof(i,j,k,6)
	      tz = czz*(uz4*0.25+uz3/dlz)/3.0
	    end if
	    cof(i,j,k,8)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	  do k=nz-1,kfinal
	    cxx = 0.5*dlxx*(cof(i,j,k,1)+cof(i,j,k,2))
	    cx = dlx*(cof(i,j,k,2)-cof(i,j,k,1))
	    cyy = 0.5*dlyy*(cof(i,j,k,3)+cof(i,j,k,4))
	    cy = dly*(cof(i,j,k,4)-cof(i,j,k,3))
	    czz = 0.5*dlzz*(cof(i,j,k,5)+cof(i,j,k,6))
	    cz = dlz*(cof(i,j,k,6)-cof(i,j,k,5))
	    call pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
	    tx = (cxx*ux4*0.5+cx*ux3)/6.0
	    ty = (cyy*uy4*0.5+cy*uy3)/6.0
	    tz = (czz*uz4*0.5+cz*uz3)/6.0
	    if(nze.ne.0 .and. k.eq.nz) then
	      czz = 0.5*dlzz*cof(i,j,k,5)
	      tz = czz*(uz4*0.25-uz3/dlz)/3.0
	    end if
	    cof(i,j,k,8)=dlxx*tx+dlyy*ty+dlzz*tz
	  end do
	end do
      end do
      return
      end


