c
c     file tmud34sp.f
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
c     test program for the MUDPACK solver mud34sp
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mud34sp.f, mud3sp.f, mudcom.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud34sp
c
c **********************************************************
c **********************************************************
c
c     a sample program/test driver for mud3sp is below. it can be
c     executed as an initial test.  the output is listed for the
c     test case described.
c
c     test the driver below by solving the separable elliptic pde
c
c        d(x**2*dp/dx)/dx + d(y**2*dp/dy)/dy + d(z**2*dp/dz)/dz
c
c        - (x+y+z)*pe(x,y,z) = r(x,y,z)
c
c
c     on the cube [0.5,1.0]**3 with specified boundary conditions at x,y,z=1.0
c     and derivative boundary conditions of the form
c
c        dp/dx - pe(0.5,y,z) = f(0.5,y,z) at x = 0.5
c
c        dp/dy - pe(x,0.5,z) = g(x,0.5,z) at y = 0.5
c
c        dp/dz - pe(x,y,0.5) = h(x,y,0.5) at z = 0.5
c
c     generate the approximation on a  65 x 49 x 41 grid using the exact
c     solution
c
c          pe(x,y,z) = (x*y*z)**4
c
c     for testing.  First mud3sp is called to produce a second-order estimate.
c     Then mud34sp is called for the fourth-order approximation
c
c ******************************************************
c     output (64 bit floating point arithmetic)
c *******************************************************
c
c     mud3sp test
c
c     input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
c     nze =  2 nzf =  1
c     ixp =  2 jyq =  3 kzr =  5
c     iex =  6 jey =  5 kez =  4
c     nx =  65 ny =  49 nz =  41 iguess =  0 maxcy =  3
c     method =  0 work space length input =  514258
c     xa =  0.50 xb =  1.00
c     yc =  0.50 yd =  1.00
c     ze =  0.50 zf =  1.00
c     tolmax =  0.000E+00
c
c     multigrid options
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     discretization call to mud3sp intl =  0
c     ierror =  0 minimum work space =  472557
c
c     approximation call to mud3sp
c     intl =  1 method =  0 iguess =  0 maxcy =  3
c     ierror =  0
c     maximum error  =   0.741E-05
c
c     mud34sp test  ierror =  0
c     maximum error  =   0.207E-06
c
c ************************************************************
c     end of output
c ************************************************************
c
c
      program tmud34sp
c
c     set grid sizes with parameter statements
c
      implicit none
c
c     set grid sizes with parameter statements
c
      integer iixp,jjyq,kkzr,iiex,jjey,kkez,llwork,nnx,nny,nnz
      parameter(iixp=2,jjyq=3,kkzr=5)
      parameter(iiex=6,jjey=5,kkez=4)
      parameter (nnx = iixp*2**(iiex-1)+1)
      parameter (nny = jjyq*2**(jjey-1)+1)
      parameter (nnz = kkzr*2**(kkez-1)+1)
c
c     set work space length approximation (see mud3sp.d)
c
      parameter(llwork = 7*(nnx+2)*(nny+2)*(nnz+2)/2 )
c
c
c     dimension solution,right hand side, and work arrays
c
      real phi(nnx,nny,nnz),rhs(nnx,nny,nnz),work(llwork)
      integer iprm(22),mgopt(4)
      real fprm(8)
c
c     put integer and floating point arguments names in contiguous
c     storeage for labelling
c
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmsp3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,x,y,z,cxx,cyy,czz,cx,cy,cz,cex,cey,cez,ce,errm
      real pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k,ierror
      equivalence(intl,iprm)
      equivalence(xa,fprm)
c
c     declare coefficient and boundary condition input subroutines external
c
      external cfx,cfy,cfz,bndc
c
c
c     set input integer parameters
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
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set method of relaxation--point gauss/seidel only for mud3sp
c
      method = 0
c
c     set default multigrid options
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     set full multigrid cycling by flagging no initial guess
c
      iguess = 0
c
c     set a limit of three w(2,1) cycles from the finest level
c
      maxcy = 3
c
c     set end points of solution cube in (x,y,z) space
c
      xa = 0.5
      xb = 1.0
      yc = 0.5
      yd = 1.0
      ze = 0.5
      zf = 1.0
c
c     set mesh increments
c
      dlx = (xb-xa)/float(nx-1)
      dly = (yd-yc)/float(ny-1)
      dlz = (zf-ze)/float(nz-1)
c
c     set for no error control (this allows phi to be equivlaenced with work)
c
      tolmax = 0.0
c
c     set right hand side in rhs using exact solution
c     and initialize phi to zero
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
	    phi(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     set specified values at upper boundaries
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
	  x = xa+(i-1)*dlx
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
      write(6,50)
   50 format(//' mud3sp test ')
c
c     print input arguments
c
      write(6,100)intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
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
c     discretize pde
c
      write(*,104) intl
  104 format(/' discretization call to mud3sp', ' intl = ', i2)
      call mud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
      write (*,105) ierror,iprm(21)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     approximate pde
c
      intl = 1
      write(*,106) intl,method,iguess,maxcy
  106 format(/' approximation call to mud3sp ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      call mud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
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
c     compute fourth-order approximation
c
      call mud34sp(work,phi,ierror)
      write (*,109) ierror
  109 format (/' mud34sp test ', ' ierror = ',i2)
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
      real phi(nx,ny,nz),errm
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,x,y,z,pxx,pyy,pzz,px,py,pz,pe
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
c     coefficient subroutines
c
      subroutine cfx(x,cxx,cx,cex)
      implicit none
      real x,cxx,cx,cex
      cxx = x*x
      cx = x+x
      cex = -x
      return
      end

      subroutine cfy(y,cyy,cy,cey)
      implicit none
      real y,cyy,cy,cey
      cyy = y*y
      cy = y+y
      cey = -y
      return
      end

      subroutine cfz(z,czz,cz,cez)
      implicit none
      real z,czz,cz,cez
      czz = z*z
      cz = z+z
      cez = -z
      return
      end

      subroutine bndc(kbdy,xory,yorz,cons,gbdy)
c
c     input derivative boundary conditions to mud3sp
c     at lower x,y,z boundaries
c
      implicit none
      integer kbdy
      real xory,yorz,cons,gbdy,x,y,z,pxx,pyy,pzz,px,py,pz,pe
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      if (kbdy.eq.1) then
c
c     x=xa surface
c
	x = xa
	y = xory
	z = yorz
	call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	cons = -1.0
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
	cons = -1.0
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
	cons = -1.0
	gbdy = pz + cons*pe
	return
      end if
      end

      subroutine exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
c
c     this subroutine is used to set an exact solution for testing mud3sp
c     (i.e., setting the rhs, boundary conditions and computing the exact
c     error)
c
      implicit none
      real x,y,z,pxx,pyy,pzz,px,py,pz,pe
      pe = (x*y*z)**4
      px = 4.*(x*y*z)**3*y*z
      py = 4.*(x*y*z)**3*x*z
      pz = 4.*(x*y*z)**3*x*y
      pxx = 12.*(x*y*z)**2*(y*z)**2
      pyy = 12.*(x*y*z)**2*(x*z)**2
      pzz = 12.*(x*y*z)**2*(x*y)**2
      return
      end
