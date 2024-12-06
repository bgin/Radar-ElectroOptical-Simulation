c
c     file tmud34.f
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
c     test program for the MUDPACK solver mud34
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mud34.f, mud3.f, mudcom.f, mud3ln.f, mud3pn.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud34
c
c **********************************************************
c **********************************************************
c
c     a sample program/test driver for mud34 is below. it can be
c     executed as an initial test.  the output is listed for the
c     test case described.
c
c     test the driver below by solving the nonseparable 3-d elliptic pde
c
c          pxx+pyy+pzz-y*z*px-x*z*py-x*y*pz-x*y*z*pe = r(x,y,z)
c
c     on a 33 by 25 by 65 grid superimposed on the (x,y,z) region
c
c          [0.5,1.0] X [1.0,2.0] X [0.25,0.75]
c
c     the solution is specified at x=0.5,1.0 and y=1.0,2.0 and there are
c     mixed derivative conditions of the form
c
c          dp/dz - (x*y)*p = g(x,z) at z=0.25
c
c          dp/dz + (x*y)*p = h(x,y) at z=0.75
c
c     One full multigrid cycle using the default multigrid options
c     and line relaxation in the z direction is first executed (see
c     tmud3.f) with mud3.  Then two additional cycles are executed
c     with iguess=1 to ensure that second-order discretization level
c     error is reached (a requirement for fourth-order mudpack solvers).
c     Finally mud34 is called to produce a fourth-order estimate.
c     the exact solution
c
c          pe(x,y,z) = exp(x*y*z)
c
c     is used to set the right hand side of the pde, boundary conditions
c     and compute error.
c
c
c ******************************************************
c     output (64 bit floating point arithmetic)
c *******************************************************
c
c     mud3 test
c
c     input arguments
c     intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
c     nze =  2 nzf =  2
c     ixp =  2 jyq =  3 kzr =  2
c     iex =  5 jey =  4 kez =  6
c     nx =  33 ny =  25 nz =  65 iguess =  0 maxcy =  1
c     method =  3 work space length input =  824224
c     xa =  0.50 xb =  1.00
c     yc =  1.00 yd =  2.00
c     ze =  0.25 zf =  0.75
c     tolmax =  0.000E+00
c
c     multigrid options
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     discretization call to mud3 intl =  0
c     ierror =  0 minimum work space =  824224
c
c     approximation call to mud3
c     intl =  1 method =  3 iguess =  0 maxcy =  1
c     ierror =  0
c     maximum error  =   0.146E-04
c
c     approximation call to mud3
c     intl =  1 method =  3 iguess =  1 maxcy =  2
c     ierror =  0
c     maximum error  =   0.146E-04
c
c     mud34 test  ierror =  0
c     maximum error  =   0.373E-06
c
c ************************************************************
c     end of output
c ************************************************************
c
      program tmud34
      implicit none
c
c     set grid sizes with parameter statements
c
      integer iixp,jjyq,kkzr,iiex,jjey,kkez,llwork,nnx,nny,nnz
      parameter(iixp=2,jjyq=3,kkzr=2)
      parameter(iiex=5,jjey=4,kkez=6)
      parameter (nnx = iixp*2**(iiex-1)+1)
      parameter (nny = jjyq*2**(jjey-1)+1)
      parameter (nnz = kkzr*2**(kkez-1)+1)
c
c     set exact minimal work space (see tmud3.f)
c
      parameter (llwork = 824224 )
c
c     dimension solution,right hand side, and work arrays
c
      real phi(nnx,nny,nnz),rhs(nnx,nny,nnz),work(llwork)
      integer iprm(23),mgopt(4)
      real fprm(8)
c
c     put integer and floating point arguments names in contiguous
c     storeage for labelling purposes
c
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,meth2,nwork,
     +              lwrkqd,itero
      common/itmud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,meth2,nwork,
     +              lwrkqd,itero
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,x,y,z,cxx,cyy,czz,cx,cy,cz,ce,errm
      real pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k,ierror
      equivalence(intl,iprm)
      equivalence(xa,fprm)
c
c     declare coefficient and boundary condition input subroutines external
c
      external cof,bndc
c
c     set for initial call
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 1
      nxb = 1
      nyc = 1
      nyd = 1
      nze = 2
      nzf = 2
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
c     set for one multigrid cycles
c
      maxcy = 1
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set method of relaxation--line in the z direction
c
      method = 3
c
c     meth2 only used in planar relaxation--but set
c
      meth2 = 0
c
c     set full multigrid cycling by flagging no initial guess at the finest
c     grid level
c
      iguess = 0
c
c     set end points of solution region in (x,y,z) space
c
      xa = 0.5
      xb = 1.0
      yc = 1.0
      yd = 2.0
      ze = 0.25
      zf = 0.75
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
	do j=1,ny
	  y = yc+float(j-1)*dly
	  do i=1,nx
	    x = xa+float(i-1)*dlx
	    call cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	    rhs(i,j,k) = cxx*pxx+cyy*pyy+czz*pzz+cx*px+cy*py+cz*pz+ce*pe
	    phi(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     set specified values at x and y boundaries in phi
c
      x = xa
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(1,j,k) = pe
	end do
      end do
      x = xb
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(nx,j,k) = pe
	end do
      end do
      y = yc
      do k=1,nz
	z = ze+(k-1)*dlz
	do i=1,nx
	  x = xa+float(i-1)*dlx
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(i,1,k) = pe
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
      write(6,50)
   50 format(//' mud3 test ')
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
  104 format(/' discretization call to mud3', ' intl = ', i2)
      call mud3(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,105) ierror,iprm(22)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     approximate pde
c
      intl = 1
      write(*,106) intl,method,iguess,maxcy
  106 format(/' approximation call to mud3 ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      call mud3(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
  108 format(' maximum error  =  ',e10.3)
c
c     execute two more cycles using this estimate to ensure second-order
c
      iguess = 1
      maxcy = 2
      write(*,106) intl,method,iguess,maxcy
      call mud3(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,107) ierror
      if (ierror.gt.0) call exit(0)
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
c
c     attempt fourth-order approximation
c
      call mud34(work,phi,ierror)
      write (*,109) ierror
  109 format(/ ' mud34 test ', ' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
      end

      subroutine error(nx,ny,nz,phi,errm)
c
c     compute the maximum error in the estimate in phi
c
      implicit none
      integer nx,ny,nz
      real phi(nx,ny,nz),errm
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
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

      subroutine cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
c
c     input pde coefficients at grid point (x,y,z) in the solution region
c     to subroutine mud3
c
      implicit none
      real x,y,z,cxx,cyy,czz,cx,cy,cz,ce
      cxx = 1.0
      cyy = 1.0
      czz = 1.0
      cx = -y*z
      cy = -x*z
      cz = -x*y
      ce = -(x*y*z)
      return
      end

      subroutine bndc(kbdy,xory,yorz,alfa,gbdy)
c
c     input mixed derivative conditions at z boundaries to mud3
c
      implicit none
      integer kbdy
      real xory,yorz,alfa,gbdy
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real x,y,z,pxx,pyy,pzz,px,py,pz,pe
      if (kbdy.eq.5) then
c
c     z = ze (lower z boundary)
c
	z = ze
	x = xory
	y = yorz
	call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	alfa = -(x*y)
	gbdy = pz + alfa*pe
	return
      end if
      if (kbdy.eq.6) then
c
c     z=zf (upper z boundary)
c
	z = zf
	x = xory
	y = yorz
	call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	alfa = (x*y)
	gbdy = pz + alfa*pe
	return
      end if
      end

      subroutine exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
c
c     this subroutine sets an exact solution for testing mud34
c
      implicit none
      real x,y,z,pe,px,py,pz,pxx,pyy,pzz
      pe = exp(x*y*z)
      px = y*z*pe
      py = x*z*pe
      pz = x*y*pe
      pxx = (y*z)*px
      pyy = (x*z)*py
      pzz = (x*y)*pz
      return
      end

