c
c     file tcud34sp.f
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
c     test program for the MUDPACK solver cud34sp
c
c ... required MUDPACK files
c
c ... see documentation and test files provided in this distribution
c
c     cud34sp.f, cud3sp.f, cudcom.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for cud34sp
c
c **********************************************************
************************************************************
c
c     a sample program/test driver for cud34sp is below. it can be
c     executed as an initial test.  the output is listed for the
c     test case described.
c
c     test cud34sp below by solving the complex separable elliptic pde
c
c        d(cmplx(x,-x)*dp/dx)/dx + d(cmplx(y,-y)*dp/dy)/dy +
c
c        d(cmplx(z,-z)*dp/dz)/dz + cmplx(1.,-1.)*(dp/dx + dp/dy + dp/dz) +
c
c        cmplx(-(x+y+z),x+y+z)*p(x,y,z) = r(x,y,z)
c
c     on the (x,y,z) region
c
c       1/4 < x < 3/4, 1/3 < y < 2/3, 1/5 < z < 4/5
c
c     with dirchlet boundary conditions at the upper x,y,z boundaries
c     and the mixed complex derirvative boundary conditions:
c
c          (1) dp/dx + cmplx(1.,1.)*p(0.5,y,z) = gbdxa(y,z) at x = 1/4
c
c          (2) dp/dy + cmplx(-1.,-1.)*p(x,0.5,z) = gbdyc(x,z) at y = 1/3
c
c          (3) dp/dz + cmplx(1.,-1.)*p(x,y,0.5) = gbdze(x,y) at z = 1/5
c
c     note the complex number multiplying "p" must be constant along
c     each surface.
c
c     use the exact solution
c
c          pe(x,y,z) = cmplx(exp(x+y),exp(y+z))
c
c     for testing and choose a 49 by 33 by 41 grid.  First cud3sp
c     is called with 5 cycles to ensure second-order discretization
c     level error is reached (a requirement for mudpack fourth-order
c     solvers).  Then cud34sp is called for a fourth-order estimate.
c
c ******************************************************
c     output (64 bit floating point arithmetic)
c *******************************************************
c
c     cud3sp test
c
c     input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
c     nze =  2 nzf =  1
c     ixp =  3 jyq =  2 kzr =  3
c     iex =  5 jey =  5 kez =  5
c     nx =  49 ny =  33 nz =  49 iguess =  0 maxcy =  5
c     method =  0 work space length input =  291605
c     xa =  0.25 xb =  0.75
c     yc =  0.33 yd =  0.67
c     ze =  0.20 zf =  0.80
c     tolmax =  0.000E+00
c
c     multigrid options
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     discretization call to cud3sp intl =  0
c     ierror =  0 minimum work space =  291605
c
c     approximation call to cud3sp
c     intl =  1 method =  0 iguess =  0 maxcy =  5
c     ierror =  0
c     maximum error  =   0.287E-04
c
c     cud34sp test  ierror =  0
c     maximum error  =   0.592E-07
c
c ************************************************************
c     end of output
c ************************************************************
c
c
      program tcud34sp
      implicit none
c
c     set grid sizes with parameter statements
c
      integer iixp,jjyq,kkzr,iiex,jjey,kkez,llwork,nnx,nny,nnz
      parameter(iixp=3,jjyq=2,kkzr=3)
      parameter(iiex=5,jjey=5,kkez=5)
      parameter (nnx = iixp*2**(iiex-1)+1)
      parameter (nny = jjyq*2**(jjey-1)+1)
      parameter (nnz = kkzr*2**(kkez-1)+1)
c
c     set minimal work space (see tcud3sp.f)
c
      parameter (llwork = 291605)
c
c     dimension solution,right hand side, and work arrays
c
      complex phi(nnx,nny,nnz),rhs(nnx,nny,nnz),work(llwork)
      integer iprm(23),mgopt(4)
      real fprm(8)
c
c     put integer and floating point arguments names in contiguous
c     storeage for labelling purposes
c
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,nwork,
     +              lwrkqd,itero
      common/itcd3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,nwork,
     +              lwrkqd,itero
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftcud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,x,y,z,errm
      complex cxx,cyy,czz,cx,cy,cz,ce,cex,cey,cez
      complex pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k,ierror
      equivalence(intl,iprm)
      equivalence(xa,fprm)
c
c     declare coefficient and boundary condition input subroutines external
c
      external cfx,cfy,cfz,bndc
c
c     set for initial call
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 2
      nxb = 1
      nyc = 2
      nyd = 1
      nze = 2
      nzf = 1
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      kzr = kkzr
      iex = iiex
      jey = jjey
      kez = kkez
      nx = nnx
      ny = nny
      nz = nnz
c
c     set  5 multigrid cycles
c
      maxcy = 5
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set point relaxation (only choice with cud3sp)
c
      method = 0
c
c     set full multigrid cycling beginning at coarsest grid
c
      iguess = 0
c
c     set end points of solution region in (x,y,z) space
c
      xa = 0.25
      xb = 0.75
      yc = 1.0/3.0
      yd = 2.0/3.0
      ze = 0.20
      zf = 0.80
c
c     set default multigrid options
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     set mesh increments
c
      dlx = (xb-xa)/float(nx-1)
      dly = (yd-yc)/float(ny-1)
      dlz = (zf-ze)/float(nz-1)
c
c     set for no error control
c
      tolmax = 0.0
c
c     set right hand side in rhs and phi to zero
c
      do k=1,nz
	z = ze+(k-1)*dlz
	call cfz(z,czz,cz,cez)
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call cfy(y,cyy,cy,cey)
	  do i=1,nx
	    x = xa+float(i-1)*dlx
	    call cfx(x,cxx,cx,cex)
	    call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	    ce = cex+cey+cez
	    rhs(i,j,k) = cxx*pxx+cyy*pyy+czz*pzz+cx*px+cy*py+cz*pz+ce*pe
	    phi(i,j,k) = (0.,0.)
	  end do
	end do
      end do
c
c     set specified values at upper x,y,z boundaries in phi
c
      x = xb
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(nx,j,k) = pe
	end do
      end do
      y = yd
      do k=1,nz
	z = ze+(k-1)*dlz
	do i=1,nx
	  x = xa+float(i-1)*dlx
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(i,ny,k) = pe
	end do
      end do
      z = zf
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+float(i-1)*dlx
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(i,j,nz) = pe
	end do
      end do
      write(*,50)
   50 format(//' cud3sp test ')
c
c     print input arguments
c
      write(*,100)intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +            nx,ny,nz,iguess,maxcy,method,nwork,xa,xb,yc,yd,ze,zf,
     +            tolmax,mgopt(1),mgopt(2),mgopt(3),mgopt(4)
  100 format(/' input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' nze = ',i2, ' nzf = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' kzr = ',i2,
     +/' iex = ',i2, ' jey = ',i2, ' kez = ',i2,
     +/' nx = ',i3,' ny = ',i3,' nz = ',i3, ' iguess = ',i2,
     +/' maxcy = ',i2,
     +/' method = ',i2, ' work space length input = ',i7,
     +/' xa = ',f5.2,' xb = ',f5.2,
     +/' yc = ',f5.2,' yd = ',f5.2,
     +/' ze = ',f5.2,' zf = ',f5.2,
     +/' tolmax = ' ,e10.3
     +//' multigrid options '
     +/' kcycle = ',i2
     +/' iprer = ',i2
     +/' ipost = ',i2
     +/' intpol = ',i2 )
c
c     discretize pde only
c
      write(*,104) intl
  104 format(/' discretization call to cud3sp', ' intl = ', i2)
      call cud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
      write (*,105) ierror,iprm(21)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     approximate pde
c
      intl = 1
      write(*,106) intl,method,iguess,maxcy
  106 format(/' approximation call to cud3sp ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      call cud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
c
c     compute and print maximum error
c
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
  108 format(' maximum error  =  ',e10.3)
c
c     attempt fourth-order estimate
c
      call cud34sp(work,phi,ierror)
      write(*,109) ierror
  109 format(/ ' cud34sp test ', ' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
      end

      subroutine error(nx,ny,nz,phi,errm)
c
c     compute the error in the estimate in phi
c
      implicit none
      integer nx,ny,nz
      complex phi(nx,ny,nz)
      real xa,xb,yc,yd,ze,zf,tolmax,relmax,errm
      common/ftcud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,x,y,z
      complex pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlz = (zf-ze)/(nz-1)
      errm = 0.0
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  do i=1,nx
	    x = xa+float(i-1)*dlx
	    call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	    errm = amax1(errm,abs(phi(i,j,k)-pe))
	  end do
	end do
      end do
      return
      end
c
c     complex coefficient subroutines
c
      subroutine cfx(x,cxx,cx,cex)
      implicit none
      real x
      complex cxx,cx,cex
      cxx = cmplx(x,-x)
      cx = (1.,-1.)
      cex = cmplx(-x,x)
      return
      end
      subroutine cfy(y,cyy,cy,cey)
      implicit none
      real y
      complex cyy,cy,cey
      cyy= cmplx(y,-y)
      cy = (1.,-1.)
      cey = cmplx(-y,y)
      return
      end
      subroutine cfz(z,czz,cz,cez)
      implicit none
      real z
      complex czz,cz,cez
      czz = cmplx(z,-z)
      cz = (1.,-1.)
      cez = cmplx(-z,z)
      return
      end

      subroutine bndc(kbdy,xory,yorz,cons,gbdy)
c
c     input complex derivative boundary conditions to cud3sp
c
      implicit none
      integer kbdy
      real xory,yorz,x,y,z
      complex cons,gbdy,pe,px,py,pz,pxx,pyy,pzz
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftcud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      if (kbdy.eq.1) then
c
c     x=xa surface
c
	x = xa
	y = xory
	z = yorz
	call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	cons = (1.,1.)
	gbdy = px+cons*pe
	return
      end if
      if (kbdy.eq.3) then
c
c     y=yc surface
c
	y = yc
	x = xory
	z = yorz
	call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	cons = (-1.,-1.)
	gbdy = py+cons*pe
	return
      end if
      if (kbdy.eq.5) then
c
c     z=ze surface
c
	z = ze
	x = xory
	y = yorz
	call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	cons = (1.,-1.)
	gbdy = pz + cons*pe
	return
      end if
      end

      subroutine exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
c
c     this subroutine is used to set an exact solution for testing cud3sp
c     (i.e., setting the rhs, boundary conditions and computing the exact
c     error)
c
      implicit none
      real x,y,z,exy,eyz
      complex pe,px,py,pz,pxx,pyy,pzz
      exy = exp(x+y)
      eyz = exp(y+z)
      pe = cmplx(exy,eyz)
      px = cmplx(exy,0.)
      pxx = px
      py = pe
      pyy = pe
      pz = cmplx(0.,eyz)
      pzz = pz
      return
      end

