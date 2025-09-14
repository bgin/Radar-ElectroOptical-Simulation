      function bcrh(xll,xrr,iz,c,a,bh,f,sgn)

c*********************************************************************72
c
cc BCRH
c
c  Author:
c
c    Paul Swarztrauber, Roland Sweet.
c
      dimension       a(*)       ,c(*)       ,bh(*)

      external f

      common /ccblk/  npp        ,k          ,eps        ,cnv        , &
                      nm         ,ncmplx     ,ik
     save :: ccblk
!$omp threadprivate (ccblk)
      xl = xll
      xr = xrr
      dx = .5*abs(xr-xl)
  101 x = .5*(xl+xr)
      if (sgn*f(x,iz,c,a,bh)) 103,105,102
  102 xr = x
      go to 104
  103 xl = x
  104 dx = .5*dx
      if (dx-cnv) 105,105,101
  105 bcrh = .5*(xl+xr)
      return
      end
      subroutine blktr1 ( n, an, bn, cn, m, am, bm, cm, idimy, y, b, &
         w1, w2, w3, wd, ww, wu, prdct, cprdct )

c*********************************************************************72
c
cc BLKTR1 is a utility routine for BLKTRI.
c
c  Discussion:
c
c    b  contains the roots of all the b polynomials
c
c    w1,w2,w3,wd,ww,wu  are all working arrays
c
c    prdct  is either prodp or prod depending on whether the boundary
c    conditions in the m direction are periodic or not
c
c    cprdct is either cprodp or cprod which are the complex versions
c    of prodp and prod. these are called in the event that some
c    of the roots of the b sub p polynomial are complex
c
c  Author:
c
c    Paul Swarztrauber, Roland Sweet.
c
      integer idimy
      integer m
      integer n

      real am(m)
      real an(n)
      real b(*)
      real bm(m)
      real bn(n)
      real cm(m)
      real cn(n)
      real cnv
      real eps
      integer ik
      integer k
      integer ncmplx
      integer nm
      integer npp
      real w1(1)
      real w2(1)
      real w3(1)
      real wd(1)
      real wu(1)
      real ww(1)
      real y(idimy,n)
c
      external cprdct
      external prdct
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
      save :: cblkt
!$omp threadprivate (cblkt)
c
c
c  Begin reduction phase.
c
      kdo = k-1
      do l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         call indxb (i2,ir,im2,nm2)
         call indxb (i1,irm1,im3,nm3)
         call indxb (i3,irm1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,y(1,i2),w3, &
                    m,am,bm,cm,wd,ww,wu)
         if = 2**k
 
         do 108 i=i4,if,i4
            if (i-nm) 101,101,108
  101       ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            call indxc (i,ir,idxc,nc)
 
            if (i-if) 102,108,108
102         continue
            call indxa (i,ir,idxa,na)
            call indxb (i-i1,irm1,im1,nm1)
            call indxb (ipi2,ir,ip2,np2)
            call indxb (ipi1,irm1,ip1,np1)
            call indxb (ipi3,irm1,ip3,np3)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w3,w1,m,am, &
                       bm,cm,wd,ww,wu)
 
            if (ipi2-nm) 105,105,103
  103       do j=1,m
               w3(j) = 0.
               w2(j) = 0.
            end do
            go to 106
 
105         continue
            call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum, &
                       y(1,ipi2),w3,m,am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w3,w2,m,am, &
                       bm,cm,wd,ww,wu)
 
106         continue
            do j=1,m
               y(j,i) = w1(j)+w2(j)+y(j,i)
            end do
 
  108    continue
      end do
 
      if (npp) 132,110,132
c
c  The periodic case is treated using the capacitance matrix method
c
  110 if = 2**k
      i = if/2
      i1 = i/2
      call indxb (i-i1,k-2,im1,nm1)
      call indxb (i+i1,k-2,ip1,np1)
      call indxb (i,k-1,iz,nz)
      call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,y(1,i),w1,m,am, &
                 bm,cm,wd,ww,wu)
      izr = i
      do j=1,m
         w2(j) = w1(j)
      end do

      do 113 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i = i2
         call indxc (i,ir,idxc,nc)
         call indxb (i,ir,iz,nz)
         call indxb (i-i1,ir-1,im1,nm1)
         call indxb (i+i1,ir-1,ip1,np1)
         call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w1,w1,m,am,bm, &
                    cm,wd,ww,wu)
         do j=1,m
            w1(j) = y(j,i)+w1(j)
         end do

         call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,w1,m,am, &
                    bm,cm,wd,ww,wu)
  113 continue

      do 118 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 117 i=i2,ifd,i4
            if (i-i2-izr) 117,114,117
  114       if (i-nm) 115,115,118
  115       call indxa (i,ir,idxa,na)
            call indxb (i,ir,iz,nz)
            call indxb (i-i1,ir-1,im1,nm1)
            call indxb (i+i1,ir-1,ip1,np1)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w2,w2,m,am, &
                       bm,cm,wd,ww,wu)
            do j=1,m
               w2(j) = y(j,i)+w2(j)
            end do

            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w2,w2,m, &
                       am,bm,cm,wd,ww,wu)
            izr = i
            if (i-nm) 117,119,117
  117    continue
  118 continue

  119 do j=1,m
         y(j,nm+1) = y(j,nm+1)-cn(nm+1)*w1(j)-an(nm+1)*w2(j)
      end do

      call indxb (if/2,k-1,im1,nm1)
      call indxb (if,k-1,ip,np)
      if (ncmplx) 121,122,121
  121 call cprdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1), &
                  y(1,nm+1),m,am,bm,cm,w1,w3,ww)
      go to 123
  122 call prdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1), &
                 y(1,nm+1),m,am,bm,cm,wd,ww,wu)
  123 do j=1,m
         w1(j) = an(1)*y(j,nm+1)
         w2(j) = cn(nm)*y(j,nm+1)
         y(j,1) = y(j,1)-w1(j)
         y(j,nm) = y(j,nm)-w2(j)
      end do

      do 126 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         i1 = i2/2
         i = i4
         call indxa (i,ir,idxa,na)
         call indxb (i-i2,ir,im2,nm2)
         call indxb (i-i2-i1,ir-1,im3,nm3)
         call indxb (i-i1,ir-1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,w1,w1,m,am, &
                    bm,cm,wd,ww,wu) 
         call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w1,w1,m,am,bm,  &
                    cm,wd,ww,wu)
         do j=1,m
            y(j,i) = y(j,i)-w1(j)
         end do

  126 continue

      izr = nm
      do 131 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         do 130 i=i4,if,i4
            ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            if (ipi2-izr) 127,128,127
  127       if (i-izr) 130,131,130
  128       call indxc (i,ir,idxc,nc)
            call indxb (ipi2,ir,ip2,np2)
            call indxb (ipi1,irm1,ip1,np1)
            call indxb (ipi3,irm1,ip3,np3)
            call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,w2,w2,m, &
                       am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w2,w2,m,am, &
                       bm,cm,wd,ww,wu)
            do j=1,m
               y(j,i) = y(j,i)-w2(j)
            end do

            izr = i
            go to 131
  130    continue
  131 continue
c
c  begin back substitution phase
c
  132 do 144 ll=1,k
         l = k-ll+1
         ir = l-1
         irm1 = ir-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 143 i=i2,ifd,i4
            if (i-nm) 133,133,143
  133       imi1 = i-i1
            imi2 = i-i2
            ipi1 = i+i1
            ipi2 = i+i2
            call indxa (i,ir,idxa,na)
            call indxc (i,ir,idxc,nc)
            call indxb (i,ir,iz,nz)
            call indxb (imi1,irm1,im1,nm1)
            call indxb (ipi1,irm1,ip1,np1)
            if (i-i2) 134,134,136
  134       do j=1,m
               w1(j) = 0.
            end do
            go to 137
  136       call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),y(1,imi2), &
                       w1,m,am,bm,cm,wd,ww,wu)
  137       if (ipi2-nm) 140,140,138
  138       do j=1,m
               w2(j) = 0.
            end do
            go to 141
  140       call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),y(1,ipi2), &
                       w2,m,am,bm,cm,wd,ww,wu)
  141       do j=1,m
               w1(j) = y(j,i)+w1(j)+w2(j)
            end do

            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,y(1,i), &
                       m,am,bm,cm,wd,ww,wu)
  143    continue
  144 continue

      return
      end
      subroutine blktri ( iflg, np, n, an, bn, cn, mp, m, am, bm, cm,
       idimy, y, ierror, w )

c*********************************************************************72
c
cc BLKTRI solves a linear system derived from a separable elliptic equation.
c
c  The routine solves a system of linear equations of the form
c
c          an(j)*x(i,j-1) + am(i)*x(i-1,j) + (bn(j)+bm(i))*x(i,j)
c
c          + cn(j)*x(i,j+1) + cm(i)*x(i+1,j) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     i+1 and i-1 are evaluated modulo m and j+1 and j-1 modulo n, i.e.,
c
c          x(i,0) = x(i,n),  x(i,n+1) = x(i,1),
c          x(0,j) = x(m,j),  x(m+1,j) = x(1,j).
c
c  These equations usually result from the discretization of
c  separable elliptic equations.  boundary conditions may be
c  dirichlet, neumann, or periodic.
c
c
c  on input
c
c     iflg
c       = 0  initialization only.  certain quantities that depend on np,
c            n, an, bn, and cn are computed and stored in the work
c            array  w.
c
c       = 1  the quantities that were computed in the initialization are
c            used to obtain the solution x(i,j).
c
c       note   a call with iflg=0 takes approximately one half the time
c              time as a call with iflg = 1.  however, the
c              initialization does not have to be repeated unless np, n,
c              an, bn, or cn change.
c
c     np
c       = 0  if an(1) and cn(n) are not both zero, which corresponds to
c            periodic boundary conditions.
c       = 1  if an(1) = cn(n) = 0.
c
c     n
c       the number of unknowns in the j-direction. n must be greater
c       than 4. the operation count is proportional to mnlog2(n), hence
c       n should be selected less than or equal to m.
c
c     an,bn,cn
c       one-dimensional arrays of length n that specify the coefficients
c       in the linear equations given above.
c
c     mp
c       = 0  if am(1) and cm(m) are not both zero, which corresponds to
c            periodic boundary conditions.
c       = 1  if am(1) = cm(m) = 0  .
c
c     m
c       the number of unknowns in the i-direction. m must be greater
c       than 4.
c
c     am,bm,cm
c       one-dimensional arrays of length m that specify the coefficients
c       in the linear equations given above.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling blktri.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c       a two-dimensional array that specifies the values of the right
c       side of the linear system of equations given above.  y must be
c       dimensioned at least m*n.
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.
c             if np=1 define k=int(log2(n))+1 and set l=2**(k+1) then
c                     w must have dimension (k-2)*l+k+5+max(2n,6m)
c
c             if np=0 define k=int(log2(n-1))+1 and set l=2**(k+1) then
c                     w must have dimension (k-2)*l+k+5+2n+max(2n,6m)
c
c       **important** for purposes of checking, the required dimension
c                     of w is computed by blktri and stored in w(1)
c                     in floating point format.
c
c  on output
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m is less than 5
c       = 2  n is less than 5
c       = 3  idimy is less than m.
c       = 4  blktri failed while computing results that depend on the
c            coefficient arrays an, bn, cn.  check these arrays.
c       = 5  an(j)*cn(j-1) is less than 0 for some j. possible reasons
c            for this condition are
c            1. the arrays an and cn are not correct
c            2. too large a grid spacing was used in the discretization
c               of the elliptic equation
c            3. the linear equations resulted from a partial
c               differential equation which was not elliptic
c
c     w
c       contains intermediate values that must not be destroyed if
c       blktri will be called again with iflg=1. w(1) contains the
c       number of locations required by w in floating point format.
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   an(n),bn(n),cn(n),am(m),bm(m),cm(m),y(idimy,n)
c     arguments      w(see argument list)
c
c     latest         june 1979
c     revision
c
c     required       blktri,blktri,prod,prodp,cprod,cprodp,compb,indxa,
c     subprograms    indxb,indxc,ppadd,psgf,ppsgf,ppspf,bsrh,tevls,
c                    epmach,store
c
c     special        the algorithm may fail if abs(bm(i)+bn(j)) is less
c     conditions     than abs(am(i))+abs(an(j))+abs(cm(i))+abs(cn(j))
c                    for some i and j. the algorithm will also fail if
c                    an(j)*cn(j-1) is less than zero for some j
c                    see the discription of the output parameter ierror.
c
c     common         cblkt,value
c     blocks
c
c     specialist     paul swarztrauber
c
c     history        version 1 september 1973
c                    version 2 april     1976
c                    version 3 june      1979
c
c     algorithm      generalized cyclic reduction (see reference below)
c
c     references     swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c                    swarztrauber p. n.,a direct method for the discrete
c                    solution of separable elliptic equations, 
c                    SIAM Journal on Numerical Analysis,11(1974) pp. 1136-1150.
c
      integer idimy
      integer m
      integer n

      real am(m)
      real an(n)
      real bm(m)
      real bn(n)
      real cm(m)
      real cn(n) 
      real w(*)
      real y(idimy,n)

      external cprod
      external cprodp
      external prod
      external prodp
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
      save :: cblkt
!$omp threadprivate (cblkt)
c
      nm = n
      ierror = 0
c
c  Test input.
c
      if ( m .lt. 5 ) then
        ierror = 1
        write(*,*)' '
        write(*,*)'BLKTRI - Fatal error!'
        write(*,*)'  M must be at least 5.'
        write(*,*)'  Your input value is M = ', m
        return
      end if

      if ( n .lt. 5 ) then
        ierror = 2
        write(*,*)' '
        write(*,*)'BLKTRI - Fatal error!'
        write(*,*)'  N must be at least 5.'
        write(*,*)'  Your input value is N = ', n
        return
      end if

      if ( idimy .lt. m ) then
        ierror = 3
        write(*,*)' '
        write(*,*)'BLKTRI - Fatal error!'
        write(*,*)'  IDIMY must be at least M.'
        write(*,*)'  Your input values are IDIMY = ', idimy
        write(*,*)'  M = ', m
        return
      end if
c
c  Set up array indices.
c
      nh = n
      npp = np
      if ( npp .ne. 0 ) then
        nh = nh+1
      end if

      ik = 2
      k = 1

10    continue

      ik = ik+ik
      k = k+1
      if ( nh .gt. ik ) go to 10

      nl = ik
      ik = ik+ik
      nl = nl-1
      iwah = (k-2)*ik+k+6
c
c  Divide W into working subarrays.
c
      if ( npp .ne. 0 ) then
        iw1 = iwah
        iwbh = iw1+nm
        w(1) = float(iw1-1+max(2*nm,6*m) )
      else
        iwbh = iwah+nm+nm
        iw1 = iwbh
        w(1) = float(iw1-1+max(2*nm,6*m))
        nm = nm-1
      end if
c
c  Compute the roots of the B polynomials.
c
      iw2 = iw1+m
      iw3 = iw2+m
      iwd = iw3+m
      iww = iwd+m
      iwu = iww+m

      if ( iflg .eq. 0 ) then
        call compb ( nl, ierror, an, bn, cn, w(2), w(iwah), w(iwbh) )
        return
      end if
c
c  Solve the linear system.
c
      if ( mp .ne. 0 ) then

        call blktr1 ( nl, an, bn, cn, m, am, bm, cm, idimy, y, w(2), &
         w(iw1), w(iw2), w(iw3), w(iwd), w(iww), w(iwu), prod,  &
        cprod )

      else

        call blktr1 ( nl, an, bn, cn, m, am, bm, cm, idimy, y, w(2),  &
         w(iw1), w(iw2), w(iw3), w(iwd), w(iww), w(iwu), prodp, &
         cprodp )

      end if

      return
      end
      function bsrh (xll,xrr,iz,c,a,bh,f,sgn)

c*********************************************************************72
c
cc BSRH
c
      dimension       a(*)       ,c(*)       ,bh(*)
c
      external f
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
      save :: cblkt
!$omp threadprivate (cblkt)
c
      xl = xll
      xr = xrr
      dx = .5*abs(xr-xl)

  101 continue

      x = .5*(xl+xr)
      if (sgn*f(x,iz,c,a,bh)) 103,105,102
  102 xr = x
      go to 104
  103 xl = x
  104 dx = .5*dx
      if (dx.gt.cnv) go to 101

  105 bsrh = .5*(xl+xr)
      return
      end
      subroutine cblkt1 (n,an,bn,cn,m,am,bm,cm,idimy,y,b,w1,w2,w3,wd, &
                        ww,wu,prdct,cprdct)

c*********************************************************************72
c
cc cblkt1 solves the linear system
c
c b  contains the roots of all the b polynomials
c w1,w2,w3,wd,ww,wu  are all working arrays
c prdct is either procp or proc depending on whether the boundary
c conditions in the m direction are periodic or not
c cprdct is either cprocp or cproc which are called if some of the zeros
c of the b polynomials are complex.
c
c
      dimension      an(*)      ,bn(*)      ,cn(*)      ,am(*)      , &
                     bm(*)      ,cm(*)      ,b(*)       ,w1(*)      , &
                     w2(*)      ,w3(*)      ,wd(*)      ,ww(*)      , &
                     wu(*)      ,y(idimy,*) 
c
      external cprdct
      external prdct
c
      common /ccblk/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
      complex         am         ,bm         ,cm         ,y          , &
                     w1         ,w2         ,w3         ,wd         , &
                    ww         ,wu
      save :: ccblk
!$omp threadprivate (ccblk)
c
c begin reduction phase
c
      kdo = k-1
      do 109 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         call inxcb (i2,ir,im2,nm2)
         call inxcb (i1,irm1,im3,nm3)
         call inxcb (i3,irm1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,y(1,i2),w3, &
                    m,am,bm,cm,wd,ww,wu)
         if = 2**k
         do 108 i=i4,if,i4
            if (i-nm) 101,101,108
  101       ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            call inxcc (i,ir,idxc,nc)
            if (i-if) 102,108,108
  102       call inxca (i,ir,idxa,na)
            call inxcb (i-i1,irm1,im1,nm1)
            call inxcb (ipi2,ir,ip2,np2)
            call inxcb (ipi1,irm1,ip1,np1)
            call inxcb (ipi3,irm1,ip3,np3)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w3,w1,m,am, &
                       bm,cm,wd,ww,wu)
            if (ipi2-nm) 105,105,103
  103       do j=1,m
               w3(j) = (0.,0.)
               w2(j) = (0.,0.)
            end do

            go to 106
  105       call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum, &
                       y(1,ipi2),w3,m,am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w3,w2,m,am, &
                       bm,cm,wd,ww,wu)
  106       do 107 j=1,m
               y(j,i) = w1(j)+w2(j)+y(j,i)
  107       continue
  108    continue
  109 continue
      if (npp) 132,110,132
c
c     the periodic case is treated using the capacitance matrix method
c
  110 if = 2**k
      i = if/2
      i1 = i/2
      call inxcb (i-i1,k-2,im1,nm1)
      call inxcb (i+i1,k-2,ip1,np1)
      call inxcb (i,k-1,iz,nz)
      call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,y(1,i),w1,m,am, &
                 bm,cm,wd,ww,wu)
      izr = i
      do 111 j=1,m
         w2(j) = w1(j)
  111 continue
      do 113 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i = i2
         call inxcc (i,ir,idxc,nc)
         call inxcb (i,ir,iz,nz)
         call inxcb (i-i1,ir-1,im1,nm1)
         call inxcb (i+i1,ir-1,ip1,np1)
         call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w1,w1,m,am,bm, &
                    cm,wd,ww,wu)
         do j=1,m
            w1(j) = y(j,i)+w1(j)
         end do
         call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,w1,m,am, &
                    bm,cm,wd,ww,wu)
  113 continue
      do 118 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 117 i=i2,ifd,i4
            if (i-i2-izr) 117,114,117
  114       if (i-nm) 115,115,118
  115       call inxca (i,ir,idxa,na)
            call inxcb (i,ir,iz,nz)
            call inxcb (i-i1,ir-1,im1,nm1)
            call inxcb (i+i1,ir-1,ip1,np1)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w2,w2,m,am, &
                       bm,cm,wd,ww,wu)
            do 116 j=1,m
               w2(j) = y(j,i)+w2(j)
  116       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w2,w2,m, &
                       am,bm,cm,wd,ww,wu)
            izr = i
            if (i-nm) 117,119,117
  117    continue
  118 continue
  119 do 120 j=1,m
         y(j,nm+1) = y(j,nm+1)-cn(nm+1)*w1(j)-an(nm+1)*w2(j)
  120 continue
      call inxcb (if/2,k-1,im1,nm1)
      call inxcb (if,k-1,ip,np)
      if (ncmplx) 121,122,121
  121 call cprdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1), &
                  y(1,nm+1),m,am,bm,cm,w1,w3,ww)
      go to 123
  122 call prdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1), &
                 y(1,nm+1),m,am,bm,cm,wd,ww,wu)
  123 do 124 j=1,m
         w1(j) = an(1)*y(j,nm+1)
         w2(j) = cn(nm)*y(j,nm+1)
         y(j,1) = y(j,1)-w1(j)
         y(j,nm) = y(j,nm)-w2(j)
  124 continue
      do 126 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         i1 = i2/2
         i = i4
         call inxca (i,ir,idxa,na)
         call inxcb (i-i2,ir,im2,nm2)
         call inxcb (i-i2-i1,ir-1,im3,nm3)
         call inxcb (i-i1,ir-1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,w1,w1,m,am, &
                    bm,cm,wd,ww,wu)
         call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w1,w1,m,am,bm, &
                    cm,wd,ww,wu)
         do 125 j=1,m
            y(j,i) = y(j,i)-w1(j)
  125    continue
  126 continue
c
      izr = nm
      do 131 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         do 130 i=i4,if,i4
            ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            if (ipi2-izr) 127,128,127
  127       if (i-izr) 130,131,130
  128       call inxcc (i,ir,idxc,nc)
            call inxcb (ipi2,ir,ip2,np2)
            call inxcb (ipi1,irm1,ip1,np1)
            call inxcb (ipi3,irm1,ip3,np3)
            call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,w2,w2,m, &
                       am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w2,w2,m,am, &
                      bm,cm,wd,ww,wu)
            do 129 j=1,m
               y(j,i) = y(j,i)-w2(j)
  129       continue
            izr = i
            go to 131
  130    continue
  131 continue
c
c begin back substitution phase
c
  132 do 144 ll=1,k
         l = k-ll+1
         ir = l-1
         irm1 = ir-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 143 i=i2,ifd,i4
            if (i-nm) 133,133,143
  133       imi1 = i-i1
            imi2 = i-i2
            ipi1 = i+i1
            ipi2 = i+i2
            call inxca (i,ir,idxa,na)
            call inxcc (i,ir,idxc,nc)
            call inxcb (i,ir,iz,nz)
            call inxcb (imi1,irm1,im1,nm1)
            call inxcb (ipi1,irm1,ip1,np1)
            if (i-i2) 134,134,136
  134       do 135 j=1,m
               w1(j) = (0.,0.)
  135       continue
            go to 137
  136       call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),y(1,imi2), &
                       w1,m,am,bm,cm,wd,ww,wu)
  137       if (ipi2-nm) 140,140,138
  138       do 139 j=1,m
               w2(j) = (0.,0.)
  139       continue
            go to 141
  140       call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),y(1,ipi2), &
                       w2,m,am,bm,cm,wd,ww,wu)
  141       do 142 j=1,m
               w1(j) = y(j,i)+w1(j)+w2(j)
  142       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,y(1,i), &
                       m,am,bm,cm,wd,ww,wu)
  143    continue
  144 continue
      return
      end
      subroutine cblktr (iflg,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y, &
                        ierror,w)

c*********************************************************************72
c
cc CBLKTR is a complex version of BLKTRI.
c
c     both routines solve a system of linear equations of the form
c
c          an(j)*x(i,j-1) + am(i)*x(i-1,j) + (bn(j)+bm(i))*x(i,j)
c
c          + cn(j)*x(i,j+1) + cm(i)*x(i+1,j) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     i+1 and i-1 are evaluated modulo m and j+1 and j-1 modulo n, i.e.,
c
c          x(i,0) = x(i,n),  x(i,n+1) = x(i,1),
c          x(0,j) = x(m,j),  x(m+1,j) = x(1,j).
c
c     these equations usually result from the discretization of
c     separable elliptic equations.  boundary conditions may be
c     dirichlet, neumann, or periodic.
c
c
c     * * * * * * * * * *     on input     * * * * * * * * * *
c
c     iflg
c       = 0  initialization only.  certain quantities that depend on np,
c            n, an, bn, and cn are computed and stored in the work
c            array  w.
c       = 1  the quantities that were computed in the initialization are
c            used to obtain the solution x(i,j).
c
c       note   a call with iflg=0 takes approximately one half the time
c              time as a call with iflg = 1  .  however, the
c              initialization does not have to be repeated unless np, n,
c              an, bn, or cn change.
c
c     np
c       = 0  if an(1) and cn(n) are not zero, which corresponds to
c            periodic bounary conditions.
c       = 1  if an(1) and cn(n) are zero.
c
c     n
c       the number of unknowns in the j-direction. n must be greater
c       than 4. the operation count is proportional to mnlog2(n), hence
c       n should be selected less than or equal to m.
c
c     an,bn,cn
c     real one-dimensional arrays of length n that specify the
c     coefficients in the linear equations given above.
c
c     mp
c       = 0  if am(1) and cm(m) are not zero, which corresponds to
c            periodic boundary conditions.
c       = 1  if am(1) = cm(m) = 0  .
c
c     m
c       the number of unknowns in the i-direction. m must be greater
c       than 4.
c
c     am,bm,cm
c     complex one-dimensional arrays of length m that specify the
c     coefficients in the linear equations given above.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling blktri.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c     a complex two-dimensional array that specifies the values of the
c     right side of the linear system of equations given above.  y must
c     be dimensioned y(idimy,n) with idimy .ge. m.
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.
c             if np=1 define k=int(log2(n))+1 and set l=2**(k+1) then
c                     w must have dimension (k-2)*l+k+5+max(2n,12m)
c
c             if np=0 define k=int(log2(n-1))+1 and set l=2**(k+1) then
c                     w must have dimension (k-2)*l+k+5+2n+max(2n,12m)
c
c       **important** for purposes of checking, the required dimension
c                     of w is computed by blktri and stored in w(1)
c                     in floating point format.
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m is less than 5
c       = 2  n is less than 5
c       = 3  idimy is less than m.
c       = 4  blktri failed while computing results that depend on the
c            coefficient arrays an, bn, cn.  check these arrays.
c       = 5  an(j)*cn(j-1) is less than 0 for some j. possible reasons
c            for this condition are
c            1. the arrays an and cn are not correct
c            2. too large a grid spacing was used in the discretization
c               of the elliptic equation
c            3. the linear equations resulted from a partial
c               differential equation which was not elliptic
c
c     w
c       contains intermediate values that must not be destroyed if
c       cblktr will be called again with iflg=1. w(1) contains the
c       number of locations required by w in floating point format.
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   an(n),bn(n),cn(n),am(m),bm(m),cm(m),y(idimy,n)
c     arguments      w(see argument list)
c
c     latest         june 1979
c     revision
c
c     required       cblktr,cblkt1,proc,procp,cproc,cprocp,ccmpb,inxca,
c     subprograms    inxcb,inxcc,cpadd,pgsf,ppgsf,pppsf,bcrh,tevlc,
c                    epmach,store
c
c     special        the algorithm may fail if abs(bm(i)+bn(j)) is less
c     conditions     than abs(am(i))+abs(an(j))+abs(cm(i))+abs(cn(j))
c                    for some i and j. the algorithm will also fail if
c                    an(j)*cn(j-1) is less than zero for some j
c                    see the discription of the output parameter ierror.
c
c     common         ccblk
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     paul swarztrauber
c
c     language       fortran
c
c     history        cblktr is a complex version of blktri (version 3)
c
c     algorithm      generalized cyclic reduction (see reference below)
c
c     space
c     required       control data 7600
c
c     portability    american national standards institute fortran.
c                    the approximate machine accuracy is computed in
c                    function epmach
c
c     required       none
c     resident
c     routines
c
c     references     swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c                    swarztrauber p. n.,a direct method for the discrete
c                    solution of separable elliptic equations, 
c                    SIAM Journal on Numerical Analysis,11(1974) pp. 1136-1150.
c
      dimension       an(*)      ,bn(*)      ,cn(*)      ,am(*)      , &
                     bm(*)      ,cm(*)      ,y(idimy,*) ,w(*)
      external        proc       ,procp      ,cproc      ,cprocp
      common /ccblk/  npp        ,k          ,eps        ,cnv        , &
                      nm         ,ncmplx     ,ik
      complex         am         ,bm         ,cm         ,y
      save :: ccblk
!$omp threadprivate (ccblk)
c
c test m and n for the proper form
c
      nm = n
      m2 = m+m
      ierror = 0
      if (m-5) 101,102,102
  101 ierror = 1
      go to 119
  102 if (nm-5) 103,104,104
  103 ierror = 2
      go to 119
  104 if (idimy-m) 105,106,106
  105 ierror = 3
      go to 119
  106 nh = n
      npp = np
      if (npp) 107,108,107
  107 nh = nh+1
  108 ik = 2
      k = 1
  109 ik = ik+ik
      k = k+1
      if (nh-ik) 110,110,109
  110 nl = ik
      ik = ik+ik
      nl = nl-1
      iwah = (k-2)*ik+k+6
      if (npp) 111,112,111
c
c     divide w into working sub arrays
c
  111 iw1 = iwah
      iwbh = iw1+nm
      w(1) = float(iw1-1+max0(2*nm,12*m))
      go to 113
  112 iwbh = iwah+nm+nm
      iw1 = iwbh
      w(1) = float(iw1-1+max0(2*nm,12*m))
      nm = nm-1
c
c  ccmpb computes the roots of the b polynomials
c
113   continue
c 113 if (ierror) 119,114,119

  114 iw2 = iw1+m2
      iw3 = iw2+m2
      iwd = iw3+m2
      iww = iwd+m2
      iwu = iww+m2
      if (iflg) 116,115,116
  115 call ccmpb (nl,ierror,an,bn,cn,w(2),w(iwah),w(iwbh))
      go to 119
  116 if (mp) 117,118,117
c
c  cblkt1 solves the linear system
c
  117 call cblkt1 (nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2), &
                  w(iw3),w(iwd),w(iww),w(iwu),proc,cproc)
      go to 119
  118 call cblkt1 (nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2), &
                 w(iw3),w(iwd),w(iww),w(iwu),procp,cprocp)
  119 continue
      return
      end
      subroutine ccmpb (n,ierror,an,bn,cn,b,ah,bh)

c*********************************************************************72
c
cc CCMPB computes the roots of the b polynomials using routine
c     tevlc which is a modification the eispack program tqlrat.
c     ierror is set to 4 if either tevlc fails or if a(j+1)*c(j) is
c     less than zero for some j.  ah,bh are temporary work arrays.
c
      dimension       an(1)      ,bn(1)      ,cn(1)      ,b(1)       , &
                     ah(1)      ,bh(1)
      common /ccblk/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
      save :: ccblk
!$omp threadprivate (ccblk)
      eps = epmach()
      bnorm = abs(bn(1))
      do 102 j=2,nm
         bnorm = amax1(bnorm,abs(bn(j)))
         arg = an(j)*cn(j-1)
         if (arg) 119,101,101
  101    b(j) = sign(sqrt(arg),an(j))
  102 continue
      cnv = eps*bnorm
      if = 2**k
      kdo = k-1
      do 108 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         ipl = i4-1
         ifd = if-i4
         do 107 i=i4,ifd,i4
            call inxcb (i,l,ib,nb)
            if (nb) 108,108,103
  103       js = i-ipl
            jf = js+nb-1
            ls = 0
            do 104 j=js,jf
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
  104       continue
            call tevlc (nb,bh,ah,ierror)
            if (ierror) 118,105,118
  105       lh = ib-1
            do 106 j=1,nb
               lh = lh+1
               b(lh) = -bh(j)
  106       continue
  107    continue
  108 continue
      do 109 j=1,nm
         b(j) = -bn(j)
  109 continue
      if (npp) 117,110,117
  110 nmp = nm+1
      nb = nm+nmp
      do 112 j=1,nb
         l1 = mod(j-1,nmp)+1
         l2 = mod(j+nm-1,nmp)+1
         arg = an(l1)*cn(l2)
         if (arg) 119,111,111
  111    bh(j) = sign(sqrt(arg),-an(l1))
         ah(j) = -bn(l1)
  112 continue
      call tevlc (nb,ah,bh,ierror)
      if (ierror) 118,113,118
  113 call inxcb (if,k-1,j2,lh)
      call inxcb (if/2,k-1,j1,lh)
      j2 = j2+1
      lh = j2
      n2m2 = j2+nm+nm-2
  114 d1 = abs(b(j1)-b(j2-1))
      d2 = abs(b(j1)-b(j2))
      d3 = abs(b(j1)-b(j2+1))
      if ((d2 .lt. d1) .and. (d2 .lt. d3)) go to 115
      b(lh) = b(j2)
      j2 = j2+1
      lh = lh+1
      if (j2-n2m2) 114,114,116
  115 j2 = j2+1
      j1 = j1+1
      if (j2-n2m2) 114,114,116
  116 b(lh) = b(n2m2+1)
      call inxcb (if,k-1,j1,j2)
      j2 = j1+nmp+nmp
      call cpadd (nm+1,ierror,an,cn,b(j1),b(j1),b(j2))
  117 return
  118 ierror = 4
      return
  119 ierror = 5
      return
      end
      subroutine chkpr4(iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,idmn,
     &ierror)

c*********************************************************************72
c
cc CHKPR4 checks the input parameters.
c
      external cofx
c
c     check definition of solution region
c
      ierror = 1
      if (a.ge.b .or. c.ge.d) return
c
c     check boundary switches
c
      ierror = 2
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) return
      ierror = 3
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) return
c
c     check first dimension in calling routine
c
      ierror = 5
      if (idmn .lt. 7) return
c
c     check m
c
      ierror = 6
      if (m.gt.(idmn-1) .or. m.lt.6) return
c
c     check n
c
      ierror = 7
      if (n .lt. 5) return
c
c     check iorder
c
      ierror = 8
      if (iorder.ne.2 .and. iorder.ne.4) return
c
c     check intl
c
c
c     check that equation is elliptic
c
      dlx = (b-a)/float(m)
      do  30 i=2,m
         xi = a+float(i-1)*dlx
         call cofx (xi,ai,bi,ci)
      if (ai.gt.0.0) go to 10
      ierror=10
      return
   10 continue
   30 continue
c
c     no error found
c
      ierror = 0
      return
      end
      subroutine chkprm (intl,iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx, &
                       cofy,idmn,ierror)

c*********************************************************************72
c
cc CHKPRM checks the input parameters for errors.
c
      external        cofx       ,cofy
c
c     check definition of solution region
c
      ierror = 1
      if (a.ge.b .or. c.ge.d) return
c
c     check boundary switches
c
      ierror = 2
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) return
      ierror = 3
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) return
c
c     check first dimension in calling routine
c
      ierror = 5
      if (idmn .lt. 7) return
c
c     check m
c
      ierror = 6
      if (m.gt.(idmn-1) .or. m.lt.6) return
c
c     check n
c
      ierror = 7
      if (n .lt. 5) return
c
c     check iorder
c
      ierror = 8
      if (iorder.ne.2 .and. iorder.ne.4) return
c
c     check intl
c
      ierror = 9
      if (intl.ne.0 .and. intl.ne.1) return
c
c     check that equation is elliptic
c
      dlx = (b-a)/float(m)
      dly = (d-c)/float(n)
      do  30 i=2,m
         xi = a+float(i-1)*dlx
         call cofx (xi,ai,bi,ci)
         do  20 j=2,n
            yj = c+float(j-1)*dly
            call cofy (yj,dj,ej,fj)
            if (ai*dj .gt. 0.0) go to  10
            ierror = 10
            return
   10       continue
   20    continue
   30 continue
c
c     no error found
c
      ierror = 0
      return
      end
      subroutine chksn4(mbdcnd,nbdcnd,alpha,beta,cofx,singlr)

c*********************************************************************72
c
cc CHKSN4 checks if the PDE that SEPX4 must solve is a singular operator.
c
      common /spl4/   kswx       ,kswy       ,k          ,l          , &
                      ait        ,bit        ,cit        ,dit        , &
                      mit        ,nit        ,is         ,ms         , &
                      js         ,ns         ,dlx        ,dly        , &
                      tdlx3      ,tdly3      ,dlx4       ,dly4
      save :: spl4
!$omp threadprivate (spl4)
      logical         singlr
      external cofx
      singlr = .false.
c
c     check if the boundary conditions are
c     entirely periodic and/or mixed
c
      if ((mbdcnd.ne.0 .and. mbdcnd.ne.3) .or. &
         (nbdcnd.ne.0 .and. nbdcnd.ne.3)) return
c
c     check that mixed conditions are pure neuman
c
      if (mbdcnd .ne. 3) go to  10
      if (alpha.ne.0.0 .or. beta.ne.0.0) return
   10 continue
c
c     check that non-derivative coefficient functions
c     are zero
c
      do  30 i=is,ms
         xi = ait+float(i-1)*dlx
         call cofx (xi,ai,bi,ci)
         if (ci .ne. 0.0) return
   30 continue
c
c     the operator must be singular if this point is reached
c
      singlr = .true.
      return
      end
      subroutine chksng (mbdcnd,nbdcnd,alpha,beta,gama,xnu,cofx,cofy, &
                        singlr)

c*********************************************************************72
c
cc CHKSNG checks if the PDE that SEPELI must solve is a singular operator.
c
      external cofx
      external cofy
c
      common /splp/  kswx       ,kswy       ,k          ,l          , &
                     ait        ,bit        ,cit        ,dit        , &
                     mit        ,nit        ,is         ,ms         , &
                     js         ,ns         ,dlx        ,dly        , &
                     tdlx3      ,tdly3      ,dlx4       ,dly4
      save :: splp
!$omp threadprivate (splp)
      logical         singlr
      singlr = .false.
c
c     check if the boundary conditions are
c     entirely periodic and/or mixed
c
      if ((mbdcnd.ne.0 .and. mbdcnd.ne.3) .or.   &
         (nbdcnd.ne.0 .and. nbdcnd.ne.3)) return
c
c     check that mixed conditions are pure neuman
c
      if (mbdcnd .ne. 3) go to  10
      if (alpha.ne.0.0 .or. beta.ne.0.0) return
   10 if (nbdcnd .ne. 3) go to  20
      if (gama.ne.0.0 .or. xnu.ne.0.0) return
   20 continue
c
c     check that non-derivative coefficient functions
c     are zero
c
      do  30 i=is,ms
         xi = ait+float(i-1)*dlx
         call cofx (xi,ai,bi,ci)
         if (ci .ne. 0.0) return
   30 continue
      do  40 j=js,ns
         yj = cit+float(j-1)*dly
         call cofy (yj,dj,ej,fj)
         if (fj .ne. 0.0) return
   40 continue
c
c     the operator must be singular if this point is reached
c
      singlr = .true.
      return
      end
      subroutine cmgnbn (nperod,n,mperod,m,a,b,c,idimy,y,ierror,w)

c*********************************************************************72
c
cc CMGNBN: complex generalized Buneman algorithm, linear equation solver.
c
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c
c     cmgnbn solves the complex linear system of equations
c
c          a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c
c          + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.,
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     0, x(i,2), or x(i,n) and x(i,n+1) may be equal to 0, x(i,n-1), or
c     x(i,1) depending on an input parameter.
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     nperod
c       indicates the values that x(i,0) and x(i,n+1) are assumed to
c       have.
c
c       = 0  if x(i,0) = x(i,n) and x(i,n+1) = x(i,1).
c       = 1  if x(i,0) = x(i,n+1) = 0  .
c       = 2  if x(i,0) = 0 and x(i,n+1) = x(i,n-1).
c       = 3  if x(i,0) = x(i,2) and x(i,n+1) = x(i,n-1).
c       = 4  if x(i,0) = x(i,2) and x(i,n+1) = 0.
c
c     n
c       the number of unknowns in the j-direction.  n must be greater
c       than 2.
c
c     mperod
c       = 0 if a(1) and c(m) are not zero
c       = 1 if a(1) = c(m) = 0
c
c     m
c       the number of unknowns in the i-direction.  n must be greater
c       than 2.
c
c     a,b,c
c       one-dimensional complex arrays of length m that specify the
c       coefficients in the linear equations given above.  if mperod = 0
c       the array elements must not depend upon the index i, but must be
c       constant.  specifically, the subroutine checks the following
c       condition
c
c             a(i) = c(1)
c             c(i) = c(1)
c             b(i) = b(1)
c
c       for i=1,2,...,m.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling cmgnbn.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c       a two-dimensional complex array that specifies the values of the
c       right side of the linear system of equations given above.  y
c       must be dimensioned at least m*n.
c
c     w
c       a one-dimensional complex array that must be provided by the
c       user for work space.  w may require up to 4*n +
c       (10 + int(log2(n)))*m locations.  the actual number of locations
c       used is computed by cmgnbn and is returned in location w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag which indicates invalid input parameters  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m .le. 2  .
c       = 2  n .le. 2
c       = 3  idimy .lt. m
c       = 4  nperod .lt. 0 or nperod .gt. 4
c       = 5  mperod .lt. 0 or mperod .gt. 1
c       = 6  a(i) .ne. c(1) or c(i) .ne. c(1) or b(i) .ne. b(1) for
c            some i=1,2,...,m.
c       = 7  a(1) .ne. 0 or c(m) .ne. 0 and mperod = 1
c
c     w
c       w(1) contains the required length of w.
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(m),b(m),c(m),y(idimy,n),w(see parameter list)
c     arguments
c
c     latest         june 1979
c     revision
c
c     subprograms    cmgnbn,cmposd,cmposn,cmposp,cmpcsg,cmpmrg,
c     required       cmptrx,cmptr3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c   history          written by roland sweet at ncar in june, 1977
c
c     algorithm      the linear system is solved by a cyclic reduction
c                    algorithm described in the reference.
c
c     space          4944(decimal) = 11520(octal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine cmgnbn is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameter nperod.  some typical values are listed
c                    in the table below.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(i) = c(i) = -0.5*b(i) = 1,       i=1,2,...,m
c
c                    and, when mperod = 1
c
c                       a(1) = c(m) = 0
c                       a(m) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and a right side y was computed.  using this
c                    array y subroutine cmgnbn was called to produce an
c                    approximate solution z.  then the relative error,
c                    defined as
c
c                       e = max(cabs(z(i,j)-x(i,j)))/max(cabs(x(i,j)))
c
c                    where the two maxima are taken over all i=1,2,...,m
c                    and j=1,2,...,n, was computed.  the value of e is
c                    given in the table below for some typical values of
c                    m and n.
c
c
c                       m (=n)    mperod    nperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         31        0         0          77     1.e-12
c                         31        1         1          45     4.e-13
c                         31        1         3          91     2.e-12
c                         32        0         0          59     7.e-14
c                         32        1         1          65     5.e-13
c                         32        1         3          97     2.e-13
c                         33        0         0          80     6.e-13
c                         33        1         1          67     5.e-13
c                         33        1         3          76     3.e-12
c                         63        0         0         350     5.e-12
c                         63        1         1         215     6.e-13
c                         63        1         3         412     1.e-11
c                         64        0         0         264     1.e-13
c                         64        1         1         287     3.e-12
c                         64        1         3         421     3.e-13
c                         65        0         0         338     2.e-12
c                         65        1         1         292     5.e-13
c                         65        1         3         329     1.e-11
c
c     portability    american national standards institue fortran.
c                    all machine dependent constants are located in the
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      sweet, r., 'a cyclic reduction algorithm for
c                    solving block tridiagonal systems of arbitrary
c                    dimensions,' siam j. on numer. anal.,
c                    14(sept., 1977), pp. 706-720.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
      complex         a          ,b          ,c          ,y          , &
                     w          ,a1
      dimension       y(idimy,1)
      dimension       w(1)       ,b(1)       ,a(1)       ,c(1)
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.0 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 102
      do 101 i=2,m
         if (cabs(a(i)-c(1)) .ne. 0.) go to 103
         if (cabs(c(i)-c(1)) .ne. 0.) go to 103
         if (cabs(b(i)-b(1)) .ne. 0.) go to 103
  101 continue
      go to 104
  102 if (cabs(a(1)).ne.0. .and. cabs(c(m)).ne.0.) ierror = 7
      go to 104
  103 ierror = 6
  104 if (ierror .ne. 0) return
      iwba = m+1
      iwbb = iwba+m
      iwbc = iwbb+m
      iwb2 = iwbc+m
      iwb3 = iwb2+m
      iww1 = iwb3+m
      iww2 = iww1+m
      iww3 = iww2+m
      iwd = iww3+m
      iwtcos = iwd+m
      iwp = iwtcos+4*n
      do 106 i=1,m
         k = iwba+i-1
         w(k) = -a(i)
         k = iwbc+i-1
         w(k) = -c(i)
         k = iwbb+i-1
         w(k) = 2.-b(i)
         do 105 j=1,n
            y(i,j) = -y(i,j)
  105    continue
  106 continue
      mp = mperod+1
      np = nperod+1
      go to (114,107),mp
  107 go to (108,109,110,111,123),np
  108 call cmposp (m,n,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),    & 
                  w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),  &
                  w(iwp))
      go to 112
  109 call cmposd (m,n,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iww1), &
                  w(iwd),w(iwtcos),w(iwp))
      go to 112
  110 call cmposn (m,n,1,2,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2), &
                  w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),   &
                  w(iwp))
      go to 112
  111 call cmposn (m,n,1,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2), &
                  w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),   &
                  w(iwp))
  112 ipstor = real(w(iww1))
      irev = 2
      if (nperod .eq. 4) go to 124
  113 go to (127,133),mp
  114 continue
c
c     reorder unknowns when mp =0
c
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 119 j=1,n
         do 115 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  115    continue
         w(mh) = 2.*y(mh,j)
         go to (117,116),modd
  116    w(m) = 2.*y(m,j)
  117    continue
         do 118 i=1,m
            y(i,j) = w(i)
  118    continue
  119 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = (0.,0.)
      w(i) = (0.,0.)
      w(k+1) = 2.*w(k+1)
      go to (120,121),modd
  120 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 122
  121 w(iwbb-1) = w(k+1)
  122 continue
      go to 107
c
c     reverse columns when nperod = 4
c
  123 irev = 1
      nby2 = n/2
  124 do 126 j=1,nby2
         mskip = n+1-j
         do 125 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  125    continue
  126 continue
      go to (110,113),irev
  127 continue
      do 132 j=1,n
         do 128 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5*(y(mhpi,j)-y(i,j))
  128    continue
         w(mh) = .5*y(mh,j)
         go to (130,129),modd
  129    w(m) = .5*y(m,j)
  130    continue
         do 131 i=1,m
            y(i,j) = w(i)
  131    continue
  132 continue
  133 continue
c
c     return storage requirements for w array.
c
      w(1) = cmplx(float(ipstor+iwp-1),0.)
      return
      end
      subroutine cmpcsg (n,ijump,fnum,fden,a)

c*********************************************************************72
c
cc CMPCSG computes required cosine values in ascending order.
c
c  when ijump .gt. 1 the routine computes values
c
c        2*cos(j*pi/l) , j=1,2,...,l and j .ne. 0(mod n/ijump+1)
c
c     where l = ijump*(n/ijump+1).
c
c
c     when ijump = 1 it computes
c
c            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... ,n
c
c     where
c        fnum = 0.5, fden = 0.0,  for regular reduction values
c        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
c        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
c        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
c                                in cmposn only.
c
      complex         a
      dimension       a(*)

      pi = pimach()
      if (n .eq. 0) go to 105
      if (ijump .eq. 1) go to 103
      k3 = n/ijump+1
      k4 = k3-1
      pibyn = pi/float(n+ijump)
      do 102 k=1,ijump
         k1 = (k-1)*k3
         k5 = (k-1)*k4
         do 101 i=1,k4
            x = k1+i
            k2 = k5+i
            a(k2) = cmplx(-2.*cos(x*pibyn),0.)
  101    continue
  102 continue
      go to 105
  103 continue
      np1 = n+1
      y = pi/(float(n)+fden)
      do 104 i=1,n
         x = float(np1-i)-fnum
         a(i) = cmplx(2.*cos(x*y),0.)
  104 continue
  105 continue
      return
      end
      subroutine cmpmrg (tcos,i1,m1,i2,m2,i3)

c*********************************************************************72
c
cc CMPMRG merges two ascending strings of numbers.
c
      complex         tcos       ,x          ,y
      dimension       tcos(*)
c
c
c     this subroutine merges two ascending strings of numbers in the
c     array tcos.  the first string is of length m1 and starts at
c     tcos(i1+1).  the second string is of length m2 and starts at
c     tcos(i2+1).  the merged string goes into tcos(i3+1).
c
c
      j1 = 1
      j2 = 1
      j = i3
      if (m1 .eq. 0) go to 107
      if (m2 .eq. 0) go to 104
  101 j = j+1
      l = j1+i1
      x = tcos(l)
      l = j2+i2
      y = tcos(l)
      if (real(x-y)) 102,102,103
  102 tcos(j) = x
      j1 = j1+1
      if (j1 .gt. m1) go to 106
      go to 101
  103 tcos(j) = y
      j2 = j2+1
      if (j2 .le. m2) go to 101
      if (j1 .gt. m1) go to 109
  104 k = j-j1+1
      do 105 j=j1,m1
         m = k+j
         l = j+i1
         tcos(m) = tcos(l)
  105 continue
      go to 109
  106 continue
      if (j2 .gt. m2) go to 109
  107 k = j-j2+1
      do 108 j=j2,m2
         m = k+j
         l = j+i2
         tcos(m) = tcos(l)
  108 continue
  109 continue
      return
      end
      subroutine cmposd (mr,nr,istag,ba,bb,bc,q,idimq,b,w,d,tcos,p)

c*********************************************************************72
c
cc CMPOSD solves Poisson's equation for Dirichlet boundary conditions.
c
c     istag = 1 if the last diagonal block is the matrix a.
c     istag = 2 if the last diagonal block is the matrix a+i.
c
      complex        ba         ,bb         ,bc         ,q          , &
                     b          ,w          ,d          ,tcos       , &
                     p          ,t
      dimension      q(idimq,1) ,ba(1)      ,bb(1)      ,bc(1)      , &
                     tcos(1)    ,b(1)       ,d(1)       ,w(1)       , &
                     p(1)
      m = mr
      n = nr
      fi = 1./float(istag)
      ip = -m
      ipstor = 0
      jsh = 0
      go to (101,102),istag
  101 kr = 0
      irreg = 1
      if (n .gt. 1) go to 106
      tcos(1) = (0.,0.)
      go to 103
  102 kr = 1
      jstsav = 1
      irreg = 2
      if (n .gt. 1) go to 106
      tcos(1) = cmplx(-1.,0.)
  103 do 104 i=1,m
         b(i) = q(i,1)
  104 continue
      call cmptrx (1,0,m,ba,bb,bc,b,tcos,d,w)
      do 105 i=1,m
         q(i,1) = b(i)
  105 continue
      go to 183
  106 lr = 0
      do 107 i=1,m
         p(i) = cmplx(0.,0.)
  107 continue
      nun = n
      jst = 1
      jsp = n
c
c     irreg = 1 when no irregularities have occurred, otherwise it is 2.
c
  108 l = 2*jst
      nodd = 2-2*((nun+1)/2)+nun
c
c     nodd = 1 when nun is odd, otherwise it is 2.
c
      go to (110,109),nodd
  109 jsp = jsp-l
      go to 111
  110 jsp = jsp-jst
      if (irreg .ne. 1) jsp = jsp-l
  111 continue
c
c     regular reduction
c
      call cmpcsg (jst,1,0.5,0.0,tcos)
      if (l .gt. jsp) go to 118
      do 117 j=l,jsp,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         jm3 = jm2-jsh
         jp3 = jp2+jsh
         if (jst .ne. 1) go to 113
         do 112 i=1,m
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  112    continue
         go to 115
  113    do 114 i=1,m
            t = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = t+q(i,j)-q(i,jm3)-q(i,jp3)
            q(i,j) = t
  114    continue
  115    continue
         call cmptrx (jst,0,m,ba,bb,bc,b,tcos,d,w)
         do 116 i=1,m
            q(i,j) = q(i,j)+b(i)
  116    continue
  117 continue
c
c     reduction for last unknown
c
  118 go to (119,136),nodd
  119 go to (152,120),irreg
c
c     odd number of unknowns
c
  120 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (123,121),istag
  121 continue
      if (jst .ne. 1) go to 123
      do 122 i=1,m
         b(i) = q(i,j)
         q(i,j) = cmplx(0.,0.)
  122 continue
      go to 130
  123 go to (124,126),noddpr
  124 do 125 i=1,m
         ip1 = ip+i
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+p(ip1)+q(i,j)
  125 continue
      go to 128
  126 do 127 i=1,m
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+q(i,jp2)-q(i,jp1)+q(i,j)
  127 continue
  128 do 129 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  129 continue
  130 call cmptrx (jst,0,m,ba,bb,bc,b,tcos,d,w)
      ip = ip+m
      ipstor = max0(ipstor,ip+m)
      do 131 i=1,m
         ip1 = ip+i
         p(ip1) = q(i,j)+b(i)
         b(i) = q(i,jp2)+p(ip1)
  131 continue
      if (lr .ne. 0) go to 133
      do 132 i=1,jst
         krpi = kr+i
         tcos(krpi) = tcos(i)
  132 continue
      go to 134
  133 continue
      call cmpcsg (lr,jstsav,0.,fi,tcos(jst+1))
      call cmpmrg (tcos,0,jst,jst,lr,kr)
  134 continue
      call cmpcsg (kr,jstsav,0.0,fi,tcos)
      call cmptrx (kr,kr,m,ba,bb,bc,b,tcos,d,w)
      do 135 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+b(i)+p(ip1)
  135 continue
      lr = kr
      kr = kr+l
      go to 152
c
c     even number of unknowns
c
  136 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (137,138),irreg
  137 continue
      jstsav = jst
      ideg = jst
      kr = l
      go to 139
  138 call cmpcsg (kr,jstsav,0.0,fi,tcos)
      call cmpcsg (lr,jstsav,0.0,fi,tcos(kr+1))
      ideg = kr
      kr = kr+jst
  139 if (jst .ne. 1) go to 141
      irreg = 2
      do 140 i=1,m
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  140 continue
      go to 150
  141 do 142 i=1,m
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  142 continue
      go to (143,145),irreg
  143 do 144 i=1,m
         q(i,j) = q(i,jm2)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
  144 continue
      irreg = 2
      go to 150
  145 continue
      go to (146,148),noddpr
  146 do 147 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+p(ip1)
  147 continue
      ip = ip-m
      go to 150
  148 do 149 i=1,m
         q(i,j) = q(i,jm2)+q(i,j)-q(i,jm1)
  149 continue
  150 call cmptrx (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      do 151 i=1,m
         q(i,j) = q(i,j)+b(i)
  151 continue
  152 nun = nun/2
      noddpr = nodd
      jsh = jst
      jst = 2*jst
      if (nun .ge. 2) go to 108
c
c     start solution.
c
      j = jsp
      do 153 i=1,m
         b(i) = q(i,j)
  153 continue
      go to (154,155),irreg
  154 continue
      call cmpcsg (jst,1,0.5,0.0,tcos)
      ideg = jst
      go to 156
  155 kr = lr+jst
      call cmpcsg (kr,jstsav,0.0,fi,tcos)
      call cmpcsg (lr,jstsav,0.0,fi,tcos(kr+1))
      ideg = kr
  156 continue
      call cmptrx (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      jm1 = j-jsh
      jp1 = j+jsh
      go to (157,159),irreg
  157 do 158 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  158 continue
      go to 164
  159 go to (160,162),noddpr
  160 do 161 i=1,m
         ip1 = ip+i
         q(i,j) = p(ip1)+b(i)
  161 continue
      ip = ip-m
      go to 164
  162 do 163 i=1,m
         q(i,j) = q(i,j)-q(i,jm1)+b(i)
  163 continue
  164 continue
c
c     start back substitution.
c
      jst = jst/2
      jsh = jst/2
      nun = 2*nun
      if (nun .gt. n) go to 183
      do 182 j=jst,n,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         if (j .gt. jst) go to 166
         do 165 i=1,m
            b(i) = q(i,j)+q(i,jp2)
  165    continue
         go to 170
  166    if (jp2 .le. n) go to 168
         do 167 i=1,m
            b(i) = q(i,j)+q(i,jm2)
  167    continue
         if (jst .lt. jstsav) irreg = 1
         go to (170,171),irreg
  168    do 169 i=1,m
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  169    continue
  170    continue
         call cmpcsg (jst,1,0.5,0.0,tcos)
         ideg = jst
         jdeg = 0
         go to 172
  171    if (j+l .gt. n) lr = lr-jst
         kr = jst+lr
         call cmpcsg (kr,jstsav,0.0,fi,tcos)
         call cmpcsg (lr,jstsav,0.0,fi,tcos(kr+1))
         ideg = kr
         jdeg = lr
  172    continue
         call cmptrx (ideg,jdeg,m,ba,bb,bc,b,tcos,d,w)
         if (jst .gt. 1) go to 174
         do 173 i=1,m
            q(i,j) = b(i)
  173    continue
         go to 182
  174    if (jp2 .gt. n) go to 177
  175    do 176 i=1,m
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  176    continue
         go to 182
  177    go to (175,178),irreg
  178    if (j+jsh .gt. n) go to 180
         do 179 i=1,m
            ip1 = ip+i
            q(i,j) = b(i)+p(ip1)
  179    continue
         ip = ip-m
         go to 182
  180    do 181 i=1,m
            q(i,j) = b(i)+q(i,j)-q(i,jm1)
  181    continue
  182 continue
      l = l/2
      go to 164
  183 continue
c
c     return storage requirements for p vectors.
c
      w(1) = cmplx(float(ipstor),0.)
      return
      end
      subroutine cmposn (m,n,istag,mixbnd,a,bb,c,q,idimq,b,b2,b3,w,w2, &
                        w3,d,tcos,p)

c*********************************************************************72
c
cc CMPOSN solves Poisson's equation with Neumann boundary conditions.
c
c     istag = 1 if the last diagonal block is a.
c     istag = 2 if the last diagonal block is a-i.
c     mixbnd = 1 if have neumann boundary conditions at both boundaries.
c     mixbnd = 2 if have neumann boundary conditions at bottom and
c     dirichlet condition at top.  (for this case, must have istag = 1.)
c
      complex         a          ,bb         ,c          ,q         ,  &
                     b          ,b2         ,b3         ,w          ,  &
                     w2         ,w3         ,d          ,tcos       ,  &
                     p          ,fi         ,t
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) , &
     &                b(*)       ,b2(*)      ,b3(*)      ,w(*)       , &
     &                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    , &
     &                k(4)       ,p(*)
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
      fistag = 3-istag
      fnum = 1./float(istag)
      fden = 0.5*float(istag-1)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      go to (101,103),istag
  101 continue
      do 102 i=1,mr
         q(i,n) = .5*q(i,n)
  102 continue
      go to (103,104),mixbnd
  103 if (n .le. 3) go to 155
  104 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      go to (105,106),mixbnd
  105 jstart = 1
      go to 107
  106 jstart = jr
      nrod = 1-nrod
  107 continue
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      call cmpcsg (i2r,1,0.5,0.0,tcos)
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 108
      j = jr
      go to 116
  108 continue
c
c     regular reduction.
c
      do 115 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 109
         jm1 = jp1
         jm2 = jp2
         jm3 = jp3
  109    continue
         if (i2r .ne. 1) go to 111
         if (j .eq. 1) jm2 = jp2
         do 110 i=1,mr
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  110    continue
         go to 113
  111    continue
         do 112 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  112    continue
  113    continue
         call cmptrx (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 114 i=1,mr
            q(i,j) = q(i,j)+b(i)
  114    continue
c
c     end of reduction for regular unknowns.
c
  115 continue
c
c     begin special reduction for last unknown.
c
      j = jstop+jr
  116 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 128
c
c     odd number of unknowns
c
      if (i2r .ne. 1) go to 118
      do 117 i=1,mr
         b(i) = fistag*q(i,j)
         q(i,j) = q(i,jm2)
  117 continue
      go to 126
  118 do 119 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  119 continue
      if (nrodpr .ne. 0) go to 121
      do 120 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  120 continue
      ip = ip-mr
      go to 123
  121 continue
      do 122 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  122 continue
  123 if (lr .eq. 0) go to 124
      call cmpcsg (lr,1,0.5,fden,tcos(kr+1))
      go to 126
  124 continue
      do 125 i=1,mr
         b(i) = fistag*b(i)
  125 continue
  126 continue
      call cmpcsg (kr,1,0.5,fden,tcos)
      call cmptrx (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 127 i=1,mr
         q(i,j) = q(i,j)+b(i)
  127 continue
      kr = kr+i2r
      go to 151
  128 continue
c
c     even number of unknowns
c
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 135
      do 129 i=1,mr
         b(i) = q(i,j)
  129 continue
      call cmptrx (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      go to (133,130),istag
  130 do 131 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  131 continue
      tcos(1) = cmplx(1.,0.)
      tcos(2) = cmplx(0.,0.)
      call cmptrx (1,1,mr,a,bb,c,b,tcos,d,w)
      do 132 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  132 continue
      go to 150
  133 continue
      do 134 i=1,mr
         p(i) = b(i)
         q(i,j) = q(i,jm2)+2.*q(i,jp2)+3.*b(i)
  134 continue
      go to 150
  135 continue
      do 136 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  136 continue
      if (nrodpr .ne. 0) go to 138
      do 137 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  137 continue
      go to 140
  138 continue
      do 139 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  139 continue
  140 continue
      call cmptrx (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max0(ipstor,ip+mr)
      do 141 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  141 continue
      if (lr .eq. 0) go to 142
      call cmpcsg (lr,1,0.5,fden,tcos(i2r+1))
      call cmpmrg (tcos,0,i2r,i2r,lr,kr)
      go to 144
  142 do 143 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  143 continue
  144 call cmpcsg (kr,1,0.5,fden,tcos)
      if (lr .ne. 0) go to 145
      go to (146,145),istag
  145 continue
      call cmptrx (kr,kr,mr,a,bb,c,b,tcos,d,w)
      go to 148
  146 continue
      do 147 i=1,mr
         b(i) = fistag*b(i)
  147 continue
  148 continue
      do 149 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  149 continue
  150 continue
      lr = kr
      kr = kr+jr
  151 continue
      go to (152,153),mixbnd
  152 nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 155
      go to 154
  153 nr = nlast/jr
      if (nr .le. 1) go to 192
  154 i2r = jr
      nrodpr = nrod
      go to 104
  155 continue
c
c      begin solution
c
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 184
      if (lr .ne. 0) go to 170
      if (n .ne. 3) go to 161
c
c     case n = 3.
c
      go to (156,168),istag
  156 continue
      do 157 i=1,mr
         b(i) = q(i,2)
  157 continue
      tcos(1) = cmplx(0.,0.)
      call cmptrx (1,0,mr,a,bb,c,b,tcos,d,w)
      do 158 i=1,mr
         q(i,2) = b(i)
         b(i) = 4.*b(i)+q(i,1)+2.*q(i,3)
  158 continue
      tcos(1) = cmplx(-2.,0.)
      tcos(2) = cmplx(2.,0.)
      i1 = 2
      i2 = 0
      call cmptrx (i1,i2,mr,a,bb,c,b,tcos,d,w)
      do 159 i=1,mr
         q(i,2) = q(i,2)+b(i)
         b(i) = q(i,1)+2.*q(i,2)
  159 continue
      tcos(1) = (0.,0.)
      call cmptrx (1,0,mr,a,bb,c,b,tcos,d,w)
      do 160 i=1,mr
         q(i,1) = b(i)
  160 continue
      jr = 1
      i2r = 0
      go to 194
c
c     case n = 2**p+1
c
  161 continue
      go to (162,170),istag
  162 continue
      do 163 i=1,mr
         b(i) = q(i,j)+.5*q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  163 continue
      call cmpcsg (jr,1,0.5,0.0,tcos)
      call cmptrx (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 164 i=1,mr
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
         b(i) = q(i,1)+2.*q(i,nlast)+4.*q(i,j)
  164 continue
      jr2 = 2*jr
      call cmpcsg (jr,1,0.0,0.0,tcos)
      do 165 i=1,jr
         i1 = jr+i
         i2 = jr+1-i
         tcos(i1) = -tcos(i2)
  165 continue
      call cmptrx (jr2,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  166 continue
      call cmpcsg (jr,1,0.5,0.0,tcos)
      call cmptrx (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 167 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  167 continue
      go to 194
c
c     case of general n with nr = 3 .
c
  168 do 169 i=1,mr
         b(i) = q(i,2)
         q(i,2) = (0.,0.)
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  169 continue
      jr = 1
      i2r = 0
      j = 2
      go to 177
  170 continue
      do 171 i=1,mr
         b(i) = .5*q(i,1)-q(i,jm1)+q(i,j)
  171 continue
      if (nrod .ne. 0) go to 173
      do 172 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  172 continue
      go to 175
  173 do 174 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  174 continue
  175 continue
      do 176 i=1,mr
         t = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+2.*t
  176 continue
  177 continue
      k1 = kr+2*jr-1
      k2 = kr+jr
      tcos(k1+1) = (-2.,0.)
      k4 = k1+3-istag
      call cmpcsg (k2+istag-2,1,0.0,fnum,tcos(k4))
      k4 = k1+k2+1
      call cmpcsg (jr-1,1,0.0,1.0,tcos(k4))
      call cmpmrg (tcos,k1,k2,k1+k2,jr-1,0)
      k3 = k1+k2+lr
      call cmpcsg (jr,1,0.5,0.0,tcos(k3+1))
      k4 = k3+jr+1
      call cmpcsg (kr,1,0.5,fden,tcos(k4))
      call cmpmrg (tcos,k3,jr,k3+jr,kr,k1)
      if (lr .eq. 0) go to 178
      call cmpcsg (lr,1,0.5,fden,tcos(k4))
      call cmpmrg (tcos,k3,jr,k3+jr,lr,k3-lr)
      call cmpcsg (kr,1,0.5,fden,tcos(k4))
  178 k3 = kr
      k4 = kr
      call cmptr3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 179 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  179 continue
      tcos(1) = (2.,0.)
      call cmptrx (1,0,mr,a,bb,c,b,tcos,d,w)
      do 180 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  180 continue
      call cmpcsg (jr,1,0.5,0.0,tcos)
      call cmptrx (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 182
      do 181 i=1,mr
         q(i,1) = b(i)
  181 continue
      go to 194
  182 continue
      do 183 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  183 continue
      go to 194
  184 continue
      if (n .ne. 2) go to 188
c
c     case  n = 2
c
      do 185 i=1,mr
         b(i) = q(i,1)
  185 continue
      tcos(1) = (0.,0.)
      call cmptrx (1,0,mr,a,bb,c,b,tcos,d,w)
      do 186 i=1,mr
         q(i,1) = b(i)
         b(i) = 2.*(q(i,2)+b(i))*fistag
  186 continue
      tcos(1) = cmplx(-fistag,0.)
      tcos(2) = cmplx(2.,0.)
      call cmptrx (2,0,mr,a,bb,c,b,tcos,d,w)
      do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
      jr = 1
      i2r = 0
      go to 194
  188 continue
c
c     case of general n and nr = 2 .
c
      do 189 i=1,mr
         ii = ip+i
         b3(i) = (0.,0.)
         b(i) = q(i,1)+2.*p(ii)
         q(i,1) = .5*q(i,1)-q(i,jm1)
         b2(i) = 2.*(q(i,1)+q(i,nlast))
  189 continue
      k1 = kr+jr-1
      tcos(k1+1) = (-2.,0.)
      k4 = k1+3-istag
      call cmpcsg (kr+istag-2,1,0.0,fnum,tcos(k4))
      k4 = k1+kr+1
      call cmpcsg (jr-1,1,0.0,1.0,tcos(k4))
      call cmpmrg (tcos,k1,kr,k1+kr,jr-1,0)
      call cmpcsg (kr,1,0.5,fden,tcos(k1+1))
      k2 = kr
      k4 = k1+k2+1
      call cmpcsg (lr,1,0.5,fden,tcos(k4))
      k3 = lr
      k4 = 0
      call cmptr3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 190 i=1,mr
         b(i) = b(i)+b2(i)
  190 continue
      tcos(1) = (2.,0.)
      call cmptrx (1,0,mr,a,bb,c,b,tcos,d,w)
      do 191 i=1,mr
         q(i,1) = q(i,1)+b(i)
  191 continue
      go to 194
  192 do 193 i=1,mr
         b(i) = q(i,nlast)
  193 continue
      go to 196
  194 continue
c
c     start back substitution.
c
      j = nlast-jr
      do 195 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  195 continue
  196 jm2 = nlast-i2r
      if (jr .ne. 1) go to 198
      do 197 i=1,mr
         q(i,nlast) = (0.,0.)
  197 continue
      go to 202
  198 continue
      if (nrod .ne. 0) go to 200
      do 199 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  199 continue
      ip = ip-mr
      go to 202
  200 do 201 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  201 continue
  202 continue
      call cmpcsg (kr,1,0.5,fden,tcos)
      call cmpcsg (lr,1,0.5,fden,tcos(kr+1))
      if (lr .ne. 0) go to 204
      do 203 i=1,mr
         b(i) = fistag*b(i)
  203 continue
  204 continue
      call cmptrx (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 205 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  205 continue
      nlastp = nlast
  206 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 222
      go to (207,208),mixbnd
  207 jstart = 1+jr
      go to 209
  208 jstart = jr
  209 continue
      kr = kr-jr
      if (nlast+jr .gt. n) go to 210
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 211
  210 continue
      jstop = nlast-jr
  211 continue
      lr = kr-jr
      call cmpcsg (jr,1,0.5,0.0,tcos)
      do 221 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 213
         do 212 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  212    continue
         go to 215
  213    continue
         do 214 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  214    continue
  215    continue
         if (jr .ne. 1) go to 217
         do 216 i=1,mr
            q(i,j) = (0.,0.)
  216    continue
         go to 219
  217    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 218 i=1,mr
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  218    continue
  219    continue
         call cmptrx (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 220 i=1,mr
            q(i,j) = q(i,j)+b(i)
  220    continue
  221 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 194
      go to 206
  222 continue
c
c     return storage requirements for p vectors.
c
      w(1) = cmplx(float(ipstor),0.)
      return
      end
      subroutine cmposp (m,n,a,bb,c,q,idimq,b,b2,b3,w,w2,w3,d,tcos,p)

c*********************************************************************72
c
cc CMPOSP solves poisson equation with periodic boundary conditions.
c
      complex        a          ,bb         ,c          ,q          , &
                     b          ,b2         ,b3         ,w          , &
                     w2         ,w3         ,d          ,tcos       , &
                     p          ,s          ,t
      dimension      a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) , &
                     b(*)       ,b2(*)      ,b3(*)      ,w(*)       , &
                     w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    , &
                     p(*)
      mr = m
      nr = (n+1)/2
      nrm1 = nr-1
      if (2*nr .ne. n) go to 107
c
c     even number of unknowns
c
      do 102 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 101 i=1,mr
            s = q(i,nrmj)-q(i,nrpj)
            t = q(i,nrmj)+q(i,nrpj)
            q(i,nrmj) = s
            q(i,nrpj) = t
  101    continue
  102 continue
      do 103 i=1,mr
         q(i,nr) = 2.*q(i,nr)
         q(i,n) = 2.*q(i,n)
  103 continue
      call cmposd (mr,nrm1,1,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = real(w(1))
      call cmposn (mr,nr+1,1,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d, &
                  tcos,p)
      ipstor = max0(ipstor,int(real(w(1))))
      do 105 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 104 i=1,mr
            s = .5*(q(i,nrpj)+q(i,nrmj))
            t = .5*(q(i,nrpj)-q(i,nrmj))
            q(i,nrmj) = s
            q(i,nrpj) = t
  104    continue
  105 continue
      do 106 i=1,mr
         q(i,nr) = .5*q(i,nr)
         q(i,n) = .5*q(i,n)
  106 continue
      go to 118
  107 continue
c
c     odd  number of unknowns
c
      do 109 j=1,nrm1
         nrpj = n+1-j
         do 108 i=1,mr
            s = q(i,j)-q(i,nrpj)
            t = q(i,j)+q(i,nrpj)
            q(i,j) = s
            q(i,nrpj) = t
  108    continue
  109 continue
      do 110 i=1,mr
         q(i,nr) = 2.*q(i,nr)
  110 continue
      lh = nrm1/2
      do 112 j=1,lh
         nrmj = nr-j
         do 111 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  111    continue
  112 continue
      call cmposd (mr,nrm1,2,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = real(w(1))
      call cmposn (mr,nr,2,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d, &
                  tcos,p)
      ipstor = max0(ipstor,int(real(w(1))))
      do 114 j=1,nrm1
         nrpj = nr+j
         do 113 i=1,mr
            s = .5*(q(i,nrpj)+q(i,j))
            t = .5*(q(i,nrpj)-q(i,j))
            q(i,nrpj) = t
            q(i,j) = s
  113    continue
  114 continue
      do 115 i=1,mr
         q(i,nr) = .5*q(i,nr)
  115 continue
      do 117 j=1,lh
         nrmj = nr-j
         do 116 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  116    continue
  117 continue
  118 continue
c
c     return storage requirements for p vectors.
c
      w(1) = cmplx(float(ipstor),0.)
      return
      end
      subroutine cmptr3 (m,a,b,c,k,y1,y2,y3,tcos,d,w1,w2,w3)

c*********************************************************************72
c
cc CMPTR3 solves a tridiagonal system.
c
      complex        a          ,b          ,c          ,y1         , &
                     y2         ,y3         ,tcos       ,d          , &
                     w1         ,w2         ,w3         ,x          , &
                     xx         ,z
      dimension       a(*)       ,b(*)       ,c(*)       ,k(4)       , &
                     tcos(*)    ,y1(*)      ,y2(*)      ,y3(*)       , &
                     d(*)       ,w1(*)      ,w2(*)      ,w3(*)

      mm1 = m-1
      k1 = k(1)
      k2 = k(2)
      k3 = k(3)
      k4 = k(4)
      f1 = k1+1
      f2 = k2+1
      f3 = k3+1
      f4 = k4+1
      k2k3k4 = k2+k3+k4
      if (k2k3k4 .eq. 0) go to 101
      l1 = f1/f2
      l2 = f1/f3
      l3 = f1/f4
      lint1 = 1
      lint2 = 1
      lint3 = 1
      kint1 = k1
      kint2 = kint1+k2
      kint3 = kint2+k3
  101 continue
      do 115 n=1,k1
         x = tcos(n)
         if (k2k3k4 .eq. 0) go to 107
         if (n .ne. l1) go to 103
         do 102 i=1,m
            w1(i) = y1(i)
  102    continue
  103    if (n .ne. l2) go to 105
         do 104 i=1,m
            w2(i) = y2(i)
  104    continue
  105    if (n .ne. l3) go to 107
         do 106 i=1,m
            w3(i) = y3(i)
  106    continue
  107    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y1(1) = y1(1)*z
         y2(1) = y2(1)*z
         y3(1) = y3(1)*z
         do 108 i=2,m
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y1(i) = (y1(i)-a(i)*y1(i-1))*z
            y2(i) = (y2(i)-a(i)*y2(i-1))*z
            y3(i) = (y3(i)-a(i)*y3(i-1))*z
  108    continue
         do 109 ip=1,mm1
            i = m-ip
            y1(i) = y1(i)-d(i)*y1(i+1)
            y2(i) = y2(i)-d(i)*y2(i+1)
            y3(i) = y3(i)-d(i)*y3(i+1)
  109    continue
         if (k2k3k4 .eq. 0) go to 115
         if (n .ne. l1) go to 111
         i = lint1+kint1
         xx = x-tcos(i)
         do 110 i=1,m
            y1(i) = xx*y1(i)+w1(i)
  110    continue
         lint1 = lint1+1
         l1 = (float(lint1)*f1)/f2
  111    if (n .ne. l2) go to 113
         i = lint2+kint2
         xx = x-tcos(i)
         do 112 i=1,m
            y2(i) = xx*y2(i)+w2(i)
  112    continue
         lint2 = lint2+1
         l2 = (float(lint2)*f1)/f3
  113    if (n .ne. l3) go to 115
         i = lint3+kint3
         xx = x-tcos(i)
         do 114 i=1,m
            y3(i) = xx*y3(i)+w3(i)
  114    continue
         lint3 = lint3+1
         l3 = (float(lint3)*f1)/f4
  115 continue
      return
      end
      subroutine cmptrx (idegbr,idegcr,m,a,b,c,y,tcos,d,w)

c*********************************************************************72
c
cc CMPTRX solves a system of linear equations where the
c     coefficient matrix is a rational function in the matrix given by
c     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
c
      complex        a          ,b          ,c          ,y          , &
                     tcos       ,d          ,w          ,x          , &
                     xx         ,z
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)    , &
                     tcos(*)    ,d(*)       ,w(*)
      mm1 = m-1
      fb = idegbr+1
      fc = idegcr+1
      l = fb/fc
      lint = 1
      do 108 k=1,idegbr
         x = tcos(k)
         if (k .ne. l) go to 102
         i = idegbr+lint
         xx = x-tcos(i)
         do 101 i=1,m
            w(i) = y(i)
            y(i) = xx*y(i)
  101    continue
  102    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y(1) = y(1)*z
         do 103 i=2,mm1
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y(i) = (y(i)-a(i)*y(i-1))*z
  103    continue
         z = b(m)-x-a(m)*d(mm1)
         if (cabs(z) .ne. 0.) go to 104
         y(m) = (0.,0.)
         go to 105
  104    y(m) = (y(m)-a(m)*y(mm1))/z
  105    continue
         do 106 ip=1,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)
  106    continue
         if (k .ne. l) go to 108
         do 107 i=1,m
            y(i) = y(i)+w(i)
  107    continue
         lint = lint+1
         l = (float(lint)*fb)/fc
  108 continue
      return
      end
      subroutine cofx (x,af,bf,cf)

c*********************************************************************72
c
cc COFX sets coefficients in the x-direction.
c
      af = (x+1.)**2
      bf = 2.0*(x+1.)
      cf = -x
      return
      end
      subroutine cofx4(x,af,bf,cf)

c*********************************************************************72
c
cc COFX4 sets coefficients in the x-direction.
c
      af = (x+1.)**2
      bf = 2.0*(x+1.)
      cf = -x
      return
      end
      subroutine cofy (y,df,ef,ff)

c*********************************************************************72
c
cc COFY sets coefficients in y direction
c
      df = exp(y)
      ef = 0.0
      ff = -y
      return
      end
      subroutine compb (n,ierror,an,bn,cn,b,ah,bh)

c*********************************************************************72
c
cc COMPB computes the roots of the b polynomials using subroutine
c     tevls which is a modification the eispack program tqlrat.
c     ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
c     less than zero for some j.  ah,bh are temporary work arrays.
c
      dimension       an(*)      ,bn(*)      ,cn(*)      ,b(*)       , &
                     ah(*)      ,bh(*)
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
      save :: cblkt
!$omp threadprivate (cblkt)
      eps = epmach()
      bnorm = abs(bn(1))
      do 102 j=2,nm
         bnorm = amax1(bnorm,abs(bn(j)))
         arg = an(j)*cn(j-1)
         if (arg) 119,101,101
  101    b(j) = sign(sqrt(arg),an(j))
  102 continue
      cnv = eps*bnorm
      if = 2**k
      kdo = k-1
      do 108 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         ipl = i4-1
         ifd = if-i4
         do 107 i=i4,ifd,i4
            call indxb (i,l,ib,nb)
            if (nb) 108,108,103
  103       js = i-ipl
            jf = js+nb-1
            ls = 0
            do 104 j=js,jf
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
  104       continue
            call tevls (nb,bh,ah,ierror)
            if (ierror) 118,105,118
  105       lh = ib-1
            do 106 j=1,nb
               lh = lh+1
               b(lh) = -bh(j)
  106       continue
  107    continue
  108 continue
      do 109 j=1,nm
         b(j) = -bn(j)
  109 continue
      if (npp) 117,110,117
  110 nmp = nm+1
      nb = nm+nmp
      do 112 j=1,nb
         l1 = mod(j-1,nmp)+1
         l2 = mod(j+nm-1,nmp)+1
         arg = an(l1)*cn(l2)
         if (arg) 119,111,111
  111    bh(j) = sign(sqrt(arg),-an(l1))
         ah(j) = -bn(l1)
  112 continue
      call tevls (nb,ah,bh,ierror)
      if (ierror) 118,113,118
  113 call indxb (if,k-1,j2,lh)
      call indxb (if/2,k-1,j1,lh)
      j2 = j2+1
      lh = j2
      n2m2 = j2+nm+nm-2
  114 d1 = abs(b(j1)-b(j2-1))
      d2 = abs(b(j1)-b(j2))
      d3 = abs(b(j1)-b(j2+1))
      if ((d2 .lt. d1) .and. (d2 .lt. d3)) go to 115
      b(lh) = b(j2)
      j2 = j2+1
      lh = lh+1
      if (j2-n2m2) 114,114,116
  115 j2 = j2+1
      j1 = j1+1
      if (j2-n2m2) 114,114,116
  116 b(lh) = b(n2m2+1)
      call indxb (if,k-1,j1,j2)
      j2 = j1+nmp+nmp
      call ppadd (nm+1,ierror,an,cn,b(j1),b(j1),b(j2))
  117 return
  118 ierror = 4
      return
  119 ierror = 5
      return
      end
      subroutine cosgen (n,ijump,fnum,fden,a)

c*********************************************************************72
c
cc COSGEN computes required cosine values in ascending order.
c
c     this subroutine computes required cosine values in ascending
c     order.  when ijump .gt. 1 the routine computes values
c
c        2*cos(j*pi/l) , j=1,2,...,l and j .ne. 0(mod n/ijump+1)
c
c     where l = ijump*(n/ijump+1).
c
c
c     when ijump = 1 it computes
c
c            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... ,n
c
c     where
c        fnum = 0.5, fden = 0.0,  for regular reduction values
c        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
c        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
c        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
c                                in poisn2 only.
c
      dimension       a(*)

      pi = pimach()
      if (n .eq. 0) go to 105
      if (ijump .eq. 1) go to 103
      k3 = n/ijump+1
      k4 = k3-1
      pibyn = pi/float(n+ijump)
      do 102 k=1,ijump
         k1 = (k-1)*k3
         k5 = (k-1)*k4
         do 101 i=1,k4
            x = k1+i
            k2 = k5+i
            a(k2) = -2.*cos(x*pibyn)
  101    continue
  102 continue
      go to 105
  103 continue
      np1 = n+1
      y = pi/(float(n)+fden)
      do 104 i=1,n
         x = float(np1-i)-fnum
         a(i) = 2.*cos(x*y)
  104 continue
  105 continue
      return
      end
      SUBROUTINE COSQB (N,X,WSAVE)

c*********************************************************************72
c
cc COSQB backward cosine quarter wave transform.
c
      DIMENSION       X(*)       ,WSAVE(*)
      DATA TSQRT2 /2.82842712474619/
      IF (N-2) 101,102,103
  101 X(1) = 4.*X(1)
      RETURN
  102 X1 = 4.*(X(1)+X(2))
      X(2) = TSQRT2*(X(1)-X(2))
      X(1) = X1
      RETURN
  103 CALL COSQB1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COSQB1 (N,X,W,XH)

c*********************************************************************72
c
cc COSQB1 is a utility routine for COSQB.
c
      DIMENSION       X(*)       ,W(*)       ,XH(*)
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 I=3,N,2
         XIM1 = X(I-1)+X(I)
         X(I) = X(I)-X(I-1)
         X(I-1) = XIM1
  101 CONTINUE
      X(1) = X(1)+X(1)
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) X(N) = X(N)+X(N)
      CALL RFFTB (N,X,XH)
      DO 102 K=2,NS2
         KC = NP2-K
         XH(K) = W(K-1)*X(KC)+W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K)-W(KC-1)*X(KC)
  102 CONTINUE
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO 103 K=2,NS2
         KC = NP2-K
         X(K) = XH(K)+XH(KC)
         X(KC) = XH(K)-XH(KC)
  103 CONTINUE
      X(1) = X(1)+X(1)
      RETURN
      END
      SUBROUTINE COSQF (N,X,WSAVE)

c*********************************************************************72
c
cc COSQF forward cosine quarter wave transform.
c
      DIMENSION       X(*)       ,WSAVE(*)
      DATA SQRT2 /1.4142135623731/
      IF (N-2) 102,101,103
  101 TSQX = SQRT2*X(2)
      X(2) = X(1)-TSQX
      X(1) = X(1)+TSQX
  102 RETURN
  103 CALL COSQF1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COSQF1 (N,X,W,XH)

c*********************************************************************72
c
cc COSQF1 is a utility routine for COSQF.
c
      DIMENSION       X(*)       ,W(*)       ,XH(*)
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 K=2,NS2
         KC = NP2-K
         XH(K) = X(K)+X(KC)
         XH(KC) = X(K)-X(KC)
  101 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) XH(NS2+1) = X(NS2+1)+X(NS2+1)
      DO 102 K=2,NS2
         KC = NP2-K
         X(K) = W(K-1)*XH(KC)+W(KC-1)*XH(K)
         X(KC) = W(K-1)*XH(K)-W(KC-1)*XH(KC)
  102 CONTINUE
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*XH(NS2+1)
      CALL RFFTF (N,X,XH)
      DO 103 I=3,N,2
         XIM1 = X(I-1)-X(I)
         X(I) = X(I-1)+X(I)
         X(I-1) = XIM1
  103 CONTINUE
      RETURN
      END
      SUBROUTINE COSQI (N,WSAVE)

c*********************************************************************72
c
cc COSQI initializes the cosine quarter wave transform.
c
      DIMENSION       WSAVE(*)
      DATA PIH /1.57079632679491/
      DT = PIH/FLOAT(N)
      FK = 0.
      DO 101 K=1,N
         FK = FK+1.
         WSAVE(K) = COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (N,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COST (N,X,WSAVE)

c*********************************************************************72
c
cc COST cosine transform.
c
      DIMENSION       X(*)       ,WSAVE(*)
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      IF (N-2) 106,101,102
  101 X1H = X(1)+X(2)
      X(2) = X(1)-X(2)
      X(1) = X1H
      RETURN
  102 IF (N .GT. 3) GO TO 103
      X1P3 = X(1)+X(3)
      TX2 = X(2)+X(2)
      X(2) = X(1)-X(3)
      X(1) = X1P3+TX2
      X(3) = X1P3-TX2
      RETURN
  103 C1 = X(1)-X(N)
      X(1) = X(1)+X(N)
      DO 104 K=2,NS2
         KC = NP1-K
         T1 = X(K)+X(KC)
         T2 = X(K)-X(KC)
         C1 = C1+WSAVE(KC)*T2
         T2 = WSAVE(K)*T2
         X(K) = T1-T2
         X(KC) = T1+T2
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) X(NS2+1) = X(NS2+1)+X(NS2+1)
      CALL RFFTF (NM1,X,WSAVE(N+1))
      XIM2 = X(2)
      X(2) = C1
      DO 105 I=4,N,2
         XI = X(I)
         X(I) = X(I-2)-X(I-1)
         X(I-1) = XIM2
         XIM2 = XI
  105 CONTINUE
      IF (MODN .NE. 0) X(N) = XIM2
  106 RETURN
      END
      SUBROUTINE COSTI (N,WSAVE)

c*********************************************************************72
c
cc COSTI initializes the cosine transform.
c
      DIMENSION       WSAVE(*)
      DATA PI /3.14159265358979/
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DT = PI/FLOAT(NM1)
      FK = 0.
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (NM1,WSAVE(N+1))
      RETURN
      END
      subroutine cpadd (n,ierror,a,c,cbp,bp,bh)

c*********************************************************************72
c
cc CPADD computes the eigenvalues of the periodic tridiagonal matrix
c     with coefficients an,bn,cn
c
c n is the order of the bh and bp polynomials
c     on output bp contians the eigenvalues
c cbp is the same as bp except type complex
c bh is used to temporarily store the roots of the b hat polynomial
c which enters through bp
c
      complex        cf         ,cx         ,fsg        ,hsg        , &
                     dd         ,f          ,fp         ,fpp        , &
                     cdis       ,r1         ,r2         ,r3         , &
                     cbp
      dimension       a(*)       ,c(*)       ,bp(*)      ,bh(*)      , &
                     cbp(*)
      common /ccblk/  npp        ,k          ,eps        ,cnv        , &
                     nm         ,ncmplx     ,ik
     save :: ccblk
!$omp threadprivate (ccblk)
      external        pgsf       ,pppsf      ,ppgsf
      scnv = sqrt(cnv)
      iz = n
      izm = iz-1
      izm2 = iz-2
      if (bp(n)-bp(1)) 101,142,103
  101 do 102 j=1,n
         nt = n-j
         bh(j) = bp(nt+1)
  102 continue
      go to 105
  103 do 104 j=1,n
         bh(j) = bp(j)
  104 continue
  105 ncmplx = 0
      modiz = mod(iz,2)
      is = 1
      if (modiz) 106,107,106
  106 if (a(1)) 110,142,107
  107 xl = bh(1)
      db = bh(3)-bh(1)
  108 xl = xl-db
      if (pgsf(xl,iz,c,a,bh)) 108,108,109
  109 sgn = -1.
      cbp(1) = cmplx(bcrh(xl,bh(1),iz,c,a,bh,pgsf,sgn),0.)
      is = 2
  110 if = iz-1
      if (modiz) 111,112,111
  111 if (a(1)) 112,142,115
  112 xr = bh(iz)
      db = bh(iz)-bh(iz-2)
  113 xr = xr+db
      if (pgsf(xr,iz,c,a,bh)) 113,114,114
  114 sgn = 1.
      cbp(iz) = cmplx(bcrh(bh(iz),xr,iz,c,a,bh,pgsf,sgn),0.)
      if = iz-2
  115 do 136 ig=is,if,2
         xl = bh(ig)
         xr = bh(ig+1)
         sgn = -1.
         xm = bcrh(xl,xr,iz,c,a,bh,pppsf,sgn)
         psg = pgsf(xm,iz,c,a,bh)
         if (abs(psg)-eps) 118,118,116
  116    if (psg*ppgsf(xm,iz,c,a,bh)) 117,118,119
c
c     case of a real zero
c
  117    sgn = 1.
         cbp(ig) = cmplx(bcrh(bh(ig),xm,iz,c,a,bh,pgsf,sgn),0.)
         sgn = -1.
         cbp(ig+1) = cmplx(bcrh(xm,bh(ig+1),iz,c,a,bh,pgsf,sgn),0.)
         go to 136
c
c     case of a multiple zero
c
  118    cbp(ig) = cmplx(xm,0.)
         cbp(ig+1) = cmplx(xm,0.)
         go to 136
c
c     case of a complex zero
c
  119    it = 0
         icv = 0
         cx = cmplx(xm,0.)
  120    fsg = (1.,0.)
         hsg = (1.,0.)
         fp = (0.,0.)
         fpp = (0.,0.)
         do 121 j=1,iz
            dd = 1./(cx-bh(j))
            fsg = fsg*a(j)*dd
            hsg = hsg*c(j)*dd
            fp = fp+dd
            fpp = fpp-dd*dd
  121    continue
         if (modiz) 123,122,123
  122    f = (1.,0.)-fsg-hsg
         go to 124
  123    f = (1.,0.)+fsg+hsg
  124    i3 = 0
         if (cabs(fp)) 126,126,125
  125    i3 = 1
         r3 = -f/fp
  126    i2 = 0
         if (cabs(fpp)) 132,132,127
  127    i2 = 1
         cdis = csqrt(fp**2-2.*f*fpp)
         r1 = cdis-fp
         r2 = -fp-cdis
         if (cabs(r1)-cabs(r2)) 129,129,128
  128    r1 = r1/fpp
         go to 130
  129    r1 = r2/fpp
  130    r2 = 2.*f/fpp/r1
         if (cabs(r2) .lt. cabs(r1)) r1 = r2
         if (i3) 133,133,131
  131    if (cabs(r3) .lt. cabs(r1)) r1 = r3
         go to 133
  132    r1 = r3
  133    cx = cx+r1
         it = it+1
         if (it .gt. 50) go to 142
         if (cabs(r1) .gt. scnv) go to 120
         if (icv) 134,134,135
  134    icv = 1
         go to 120
  135    cbp(ig) = cx
         cbp(ig+1) = conjg(cx)
  136 continue
      if (cabs(cbp(n))-cabs(cbp(1))) 137,142,139
  137 nhalf = n/2
      do 138 j=1,nhalf
         nt = n-j
         cx = cbp(j)
         cbp(j) = cbp(nt+1)
         cbp(nt+1) = cx
  138 continue
  139 ncmplx = 1
      do 140 j=2,iz
         if (aimag(cbp(j))) 143,140,143
  140 continue
      ncmplx = 0
      do 141 j=2,iz
         bp(j) = real(cbp(j))
  141 continue
      go to 143
  142 ierror = 4
  143 continue
      return
      end
      subroutine cproc (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,yy)

c*********************************************************************72
c
cc CPROC applies a sequence of matrix operations to the vector x and
c stores the result in y
c aa   array containing scalar multipliers of the vector x
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c na is the length of the array aa
c x,y the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w are work arrays
c isgn  determines whether or not a change in sign is made
c
      complex         y          ,d          ,w          ,bd        , &
                     crt        ,den        ,y1         ,y2         , &
                     x          ,a          ,b          ,c
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,&
                     y(*)       ,d(*)       ,w(*)       ,bd(*)      ,&
                     bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
      do 101 j=1,m
         y(j) = x(j)
  101 continue
      mm = m-1
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 109,109,103
  103 crt = bd(id)
      id = id-1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-crt)
      w(m) = y(m)/(b(m)-crt)
      do 104 j=2,mm
         k = m-j
         den = b(k+1)-crt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  104 continue
      den = b(1)-crt-c(1)*d(2)
      if (cabs(den)) 105,106,105
  105 y(1) = (y(1)-c(1)*w(2))/den
      go to 107
  106 y(1) = (1.,0.)
  107 do 108 j=2,m
         y(j) = w(j)-d(j)*y(j-1)
  108 continue
  109 if (m1) 110,110,112
  110 if (m2) 121,121,111
  111 rt = bm2(m2)
      m2 = m2-1
      go to 117
  112 if (m2) 113,113,114
  113 rt = bm1(m1)
      m1 = m1-1
      go to 117
  114 if (abs(bm1(m1))-abs(bm2(m2))) 116,116,115
  115 rt = bm1(m1)
      m1 = m1-1
      go to 117
  116 rt = bm2(m2)
      m2 = m2-1
  117 y1 = (b(1)-rt)*y(1)+c(1)*y(2)
      if (mm-2) 120,118,118
c
c matrix multiplication
c
  118 do 119 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  119 continue
  120 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)
      y(m-1) = y1
      iflg = 1
      go to 102
  121 if (ia) 124,124,122
  122 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 123 j=1,m
         y(j) = rt*y(j)
  123 continue
  124 if (iflg) 125,125,102
  125 return
      end
      subroutine cprocp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,u,yy)

c*********************************************************************72
c
cc CPROCP applies a sequence of matrix operations to the vector x and
c stores the result in y
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u are work arrays
c isgn  determines whether or not a change in sign is made
c
      complex        y          ,d          ,u          ,v          , &
                     den        ,bh         ,ym         ,am         , &
                     y1         ,y2         ,yh         ,bd         , &
                     crt        ,x          ,a          ,b          ,c
      dimension      a(*)       ,b(*)       ,c(*)       ,x(*)       , &
                     y(*)       ,d(*)       ,u(*)       ,bd(*)      , &
                     bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
      do 101 j=1,m
         y(j) = x(j)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 111,111,103
  103 crt = bd(id)
      id = id-1
      iflg = 1
c
c begin solution to system
c
      bh = b(m)-crt
      ym = y(m)
      den = b(1)-crt
      d(1) = c(1)/den
      u(1) = a(1)/den
      y(1) = y(1)/den
      v = c(m)
      if (mm2-2) 106,104,104
  104 do 105 j=2,mm2
         den = b(j)-crt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         y(j) = (y(j)-a(j)*y(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*y(j-1)
         v = -v*d(j-1)
  105 continue
  106 den = b(m-1)-crt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*y(m-2)
      den = bh-am*d(m-1)
      if (cabs(den)) 107,108,107
  107 y(m) = (ym-am*y(m-1))/den
      go to 109
  108 y(m) = (1.,0.)
  109 y(m-1) = y(m-1)-d(m-1)*y(m)
      do 110 j=2,mm
         k = m-j
         y(k) = y(k)-d(k)*y(k+1)-u(k)*y(m)
  110 continue
  111 if (m1) 112,112,114
  112 if (m2) 123,123,113
  113 rt = bm2(m2)
      m2 = m2-1
      go to 119
  114 if (m2) 115,115,116
  115 rt = bm1(m1)
      m1 = m1-1
      go to 119
  116 if (abs(bm1(m1))-abs(bm2(m2))) 118,118,117
  117 rt = bm1(m1)
      m1 = m1-1
      go to 119
  118 rt = bm2(m2)
      m2 = m2-1
c
c matrix multiplication
c
  119 yh = y(1)
      y1 = (b(1)-rt)*y(1)+c(1)*y(2)+a(1)*y(m)
      if (mm-2) 122,120,120
  120 do 121 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  121 continue
  122 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)+c(m)*yh
      y(m-1) = y1
      iflg = 1
      go to 102
  123 if (ia) 126,126,124
  124 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 125 j=1,m
         y(j) = rt*y(j)
  125 continue
  126 if (iflg) 127,127,102
  127 return
      end
      subroutine cprod (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,yy,m,a,b,c,d,w,y)

c*********************************************************************72
c
cc cprod applies a sequence of matrix operations to the vector x and
c stores the result in yy           (complex case)
c aa   array containing scalar multipliers of the vector x
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c na is the length of the array aa
c x,yy the matrix operations are applied to x and the result is yy
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w,y are working arrays
c isgn  determines whether or not a change in sign is made
c
      integer na
      integer nd
      integer nm1
      integer nm2
c
      real aa(na)
      complex bd(nd)
      real x(m)
      complex y(m)
      complex        d  ,w ,crt ,den ,y1 ,y2
                     
      dimension       a(*)       ,b(*)       ,c(*)   , &
                            d(*)       ,w(*)         , &
                     bm1(*)     ,bm2(*)           ,yy(*)
c
      do 101 j=1,m
         y(j) = cmplx(x(j),0.)
  101 continue
      mm = m-1
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 109,109,103
  103 crt = bd(id)
      id = id-1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-crt)
      w(m) = y(m)/(b(m)-crt)
      do 104 j=2,mm
         k = m-j
         den = b(k+1)-crt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  104 continue
      den = b(1)-crt-c(1)*d(2)
      if (cabs(den)) 105,106,105
  105 y(1) = (y(1)-c(1)*w(2))/den
      go to 107
  106 y(1) = (1.,0.)
  107 do 108 j=2,m
         y(j) = w(j)-d(j)*y(j-1)
  108 continue
  109 if (m1) 110,110,112
  110 if (m2) 121,121,111
  111 rt = bm2(m2)
      m2 = m2-1
      go to 117
  112 if (m2) 113,113,114
  113 rt = bm1(m1)
      m1 = m1-1
      go to 117
  114 if (abs(bm1(m1))-abs(bm2(m2))) 116,116,115
  115 rt = bm1(m1)
      m1 = m1-1
      go to 117
  116 rt = bm2(m2)
      m2 = m2-1
  117 y1 = (b(1)-rt)*y(1)+c(1)*y(2)
      if (mm-2) 120,118,118
c
c matrix multiplication
c
  118 do 119 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  119 continue
  120 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)
      y(m-1) = y1
      iflg = 1
      go to 102
  121 if (ia) 124,124,122
  122 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 123 j=1,m
         y(j) = rt*y(j)
  123 continue
  124 if (iflg) 125,125,102
  125 do 126 j=1,m
         yy(j) = real(y(j))
  126 continue
      return
      end
      subroutine cprodp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,yy,m,a,b,c,d,u,y)

c*********************************************************************72
c
cc cprodp applies a sequence of matrix operations to the vector x and
c stores the result in yy       periodic boundary conditions
c and  complex  case
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,yy the matrix operations are applied to x and the result is yy
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,y are working arrays
c isgn  determines whether or not a change in sign is made
c
      complex        y          ,d          ,u          ,v    , &
                     den        ,bh         ,ym         ,am   , &
                     y1         ,y2         ,yh         ,bd   , &
                     crt
      dimension      a(*)       ,b(*)       ,c(*)       ,x(*)   , &
                     y(*)       ,d(*)       ,u(*)       ,bd(*)  , &
                     bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
      do 101 j=1,m
         y(j) = cmplx(x(j),0.)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 111,111,103
  103 crt = bd(id)
      id = id-1
      iflg = 1
c
c begin solution to system
c
      bh = b(m)-crt
      ym = y(m)
      den = b(1)-crt
      d(1) = c(1)/den
      u(1) = a(1)/den
      y(1) = y(1)/den
      v = cmplx(c(m),0.)
      if (mm2-2) 106,104,104
  104 do 105 j=2,mm2
         den = b(j)-crt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         y(j) = (y(j)-a(j)*y(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*y(j-1)
         v = -v*d(j-1)
  105 continue
  106 den = b(m-1)-crt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*y(m-2)
      den = bh-am*d(m-1)
      if (cabs(den)) 107,108,107
  107 y(m) = (ym-am*y(m-1))/den
      go to 109
  108 y(m) = (1.,0.)
  109 y(m-1) = y(m-1)-d(m-1)*y(m)
      do 110 j=2,mm
         k = m-j
         y(k) = y(k)-d(k)*y(k+1)-u(k)*y(m)
  110 continue
  111 if (m1) 112,112,114
  112 if (m2) 123,123,113
  113 rt = bm2(m2)
      m2 = m2-1
      go to 119
  114 if (m2) 115,115,116
  115 rt = bm1(m1)
      m1 = m1-1
      go to 119
  116 if (abs(bm1(m1))-abs(bm2(m2))) 118,118,117
  117 rt = bm1(m1)
      m1 = m1-1
      go to 119
  118 rt = bm2(m2)
      m2 = m2-1
c
c matrix multiplication
c
  119 yh = y(1)
      y1 = (b(1)-rt)*y(1)+c(1)*y(2)+a(1)*y(m)
      if (mm-2) 122,120,120
  120 do 121 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  121 continue
  122 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)+c(m)*yh
      y(m-1) = y1
      iflg = 1
      go to 102
  123 if (ia) 126,126,124
  124 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 125 j=1,m
         y(j) = rt*y(j)
  125 continue
  126 if (iflg) 127,127,102
  127 do 128 j=1,m
         yy(j) = real(y(j))
  128 continue
      return
      end
      subroutine defe4(cofx,idmn,usol,grhs)

c*********************************************************************72
c
cc DEFE4 first approximates the truncation error given by
c     trun1(x,y)=dlx**2*tx+dly**2*ty where
c     tx=afun(x)*uxxxx/12.0+bfun(x)*uxxx/6.0 on the interior and
c     at the boundaries if periodic(here uxxx,uxxxx are the third
c     and fourth partial derivatives of u with respect to x).
c     tx is of the form afun(x)/3.0*(uxxxx/4.0+uxxx/dlx)
c     at x=a or x=b if the boundary condition there is mixed.
c     tx=0.0 along specified boundaries.  ty has symmetric form
c     in y with x,afun(x),bfun(x) replaced by y,dfun(y),efun(y).
c     the second order solution in usol is used to approximate
c     (via second order finite differencing) the trun1ation error
c     and the result is added to the right hand side in grhs
c     and then transferred to usol to be used as a new right
c     hand side when calling BLKTRI for a fourth order solution.
c
      common /spl4/  kswx       ,kswy       ,k          ,l          , &
                     ait        ,bit        ,cit        ,dit        , &
                     mit        ,nit        ,is         ,ms         , &
                     js         ,ns         ,dlx        ,dly        , &
                     tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: spl4
!$omp threadprivate (spl4)
      dimension       grhs(idmn,*)           ,usol(idmn,*)
      external cofx
c
c
c     compute truncation error approximation over the entire mesh
c
         do  30 i=is,ms
            xi = ait+float(i-1)*dlx
            call cofx (xi,ai,bi,ci)
         do 30 j=js,ns
c
c     compute partial derivative approximations at (xi,yj)
c
            call dx4(usol,idmn,i,j,uxxx,uxxxx)
            call dy4(usol,idmn,i,j,uyyy,uyyyy)
            tx = ai*uxxxx/12.0+bi*uxxx/6.0
             ty=uyyyy/12.0
c
c     reset form of trun1ation if at boundary which is non-periodic
c
            if (kswx.eq.1 .or. (i.gt.1 .and. i.lt.k)) go to  10
            tx = ai/3.0*(uxxxx/4.0+uxxx/dlx)
   10       if (kswy.eq.1 .or. (j.gt.1 .and. j.lt.l)) go to  20
            ty = (uyyyy/4.0+uyyy/dly)/3.0
   20 grhs(i,j)=grhs(i,j)+dly**2*(dlx**2*tx+dly**2*ty)
   30    continue
c
c     reset the right hand side in usol
c
      do  60 i=is,ms
         do  50 j=js,ns
            usol(i,j) = grhs(i,j)
   50    continue
   60 continue
      return
      end
      subroutine defer (cofx,cofy,idmn,usol,grhs)

c*********************************************************************72
c
cc DEFER first approximates the truncation error given by
c     trun1(x,y)=dlx**2*tx+dly**2*ty where
c     tx=afun(x)*uxxxx/12.0+bfun(x)*uxxx/6.0 on the interior and
c     at the boundaries if periodic(here uxxx,uxxxx are the third
c     and fourth partial derivatives of u with respect to x).
c     tx is of the form afun(x)/3.0*(uxxxx/4.0+uxxx/dlx)
c     at x=a or x=b if the boundary condition there is mixed.
c     tx=0.0 along specified boundaries.  ty has symmetric form
c     in y with x,afun(x),bfun(x) replaced by y,dfun(y),efun(y).
c     the second order solution in usol is used to approximate
c     (via second order finite differencing) the trun1ation error
c     and the result is added to the right hand side in grhs
c     and then transferred to usol to be used as a new right
c     hand side when calling BLKTRI for a fourth order solution.
c
      common /splp/   kswx       ,kswy       ,k          ,l          , &
                     ait        ,bit        ,cit        ,dit        , &
                     mit        ,nit        ,is         ,ms         , &
                     js         ,ns         ,dlx        ,dly        , &
                     tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: splp
!$omp threadprivate (splp)
      dimension       grhs(idmn,*)           ,usol(idmn,*)
      external        cofx       ,cofy
c
c     compute trun1ation error approximation over the entire mesh
c
      do  40 j=js,ns
         yj = cit+float(j-1)*dly
         call cofy (yj,dj,ej,fj)
         do  30 i=is,ms
            xi = ait+float(i-1)*dlx
            call cofx (xi,ai,bi,ci)
c
c     compute partial derivative approximations at (xi,yj)
c
            call dx (usol,idmn,i,j,uxxx,uxxxx)
            call dy (usol,idmn,i,j,uyyy,uyyyy)
            tx = ai*uxxxx/12.0+bi*uxxx/6.0
            ty = dj*uyyyy/12.0+ej*uyyy/6.0
c
c     reset form of trun1ation if at boundary which is non-periodic
c
            if (kswx.eq.1 .or. (i.gt.1 .and. i.lt.k)) go to  10
            tx = ai/3.0*(uxxxx/4.0+uxxx/dlx)
   10       if (kswy.eq.1 .or. (j.gt.1 .and. j.lt.l)) go to  20
            ty = dj/3.0*(uyyyy/4.0+uyyy/dly)
   20       grhs(i,j) = grhs(i,j)+dlx**2*tx+dly**2*ty
   30    continue
   40 continue
c
c     reset the right hand side in usol
c
      do  60 i=is,ms
         do  50 j=js,ns
            usol(i,j) = grhs(i,j)
   50    continue
   60 continue
      return
      end
      subroutine dx (u,idmn,i,j,uxxx,uxxxx)

c*********************************************************************72
c
cc DX computes second order finite difference
c     approximations to the third and fourth x
c     partial derivatives of u at the (i,j) mesh point
c
      common /splp/   kswx       ,kswy       ,k          ,l          , &
                     ait        ,bit        ,cit        ,dit        , &
                     mit        ,nit        ,is         ,ms         , &
                     js         ,ns         ,dlx        ,dly        , &
                     tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: splp
!$omp threadprivate (splp)
      dimension       u(idmn,*)

      if (i.gt.2 .and. i.lt.(k-1)) go to  50
      if (i .eq. 1) go to  10
      if (i .eq. 2) go to  30
      if (i .eq. k-1) go to  60
      if (i .eq. k) go to  80
c
c     compute partial derivative approximations at x=a
c
   10 if (kswx .eq. 1) go to  20
      uxxx = (-5.0*u(1,j)+18.0*u(2,j)-24.0*u(3,j)+14.0*u(4,j)- &
                            3.0*u(5,j))/(tdlx3)
      uxxxx = (3.0*u(1,j)-14.0*u(2,j)+26.0*u(3,j)-24.0*u(4,j)+ &
                                           11.0*u(5,j)-2.0*u(6,j))/dlx4
      return
c
c     periodic at x=a
c
   20 uxxx = (-u(k-2,j)+2.0*u(k-1,j)-2.0*u(2,j)+u(3,j))/(tdlx3)
      uxxxx = (u(k-2,j)-4.0*u(k-1,j)+6.0*u(1,j)-4.0*u(2,j)+u(3,j))/dlx4
      return
c
c     compute partial derivative approximations at x=a+dlx
c
   30 if (kswx .eq. 1) go to  40
      uxxx = (-3.0*u(1,j)+10.0*u(2,j)-12.0*u(3,j)+6.0*u(4,j)-u(5,j))/ &
            tdlx3
      uxxxx = (2.0*u(1,j)-9.0*u(2,j)+16.0*u(3,j)-14.0*u(4,j)+6.0*u(5,j)- &
                                                           u(6,j))/dlx4
      return
c
c     periodic at x=a+dlx
c
   40 uxxx = (-u(k-1,j)+2.0*u(1,j)-2.0*u(3,j)+u(4,j))/(tdlx3)
      uxxxx = (u(k-1,j)-4.0*u(1,j)+6.0*u(2,j)-4.0*u(3,j)+u(4,j))/dlx4
      return
c
c     compute partial derivative approximations on the interior
c
   50 continue
      uxxx = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
      uxxxx = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))/ &
             dlx4
      return
c
c     compute partial derivative approximations at x=b-dlx
c
   60 if (kswx .eq. 1) go to  70
      uxxx = (u(k-4,j)-6.0*u(k-3,j)+12.0*u(k-2,j)-10.0*u(k-1,j)+ &
                                                      3.0*u(k,j))/tdlx3
      uxxxx = (-u(k-5,j)+6.0*u(k-4,j)-14.0*u(k-3,j)+16.0*u(k-2,j)- &
                                          9.0*u(k-1,j)+2.0*u(k,j))/dlx4
      return
c
c     periodic at x=b-dlx
c
   70 uxxx = (-u(k-3,j)+2.0*u(k-2,j)-2.0*u(1,j)+u(2,j))/tdlx3
      uxxxx = (u(k-3,j)-4.0*u(k-2,j)+6.0*u(k-1,j)-4.0*u(1,j)+u(2,j))/ &
             dlx4
      return
c
c     compute partial derivative approximations at x=b
c
   80 uxxx = -(3.0*u(k-4,j)-14.0*u(k-3,j)+24.0*u(k-2,j)-18.0*u(k-1,j)+  &
                                                      5.0*u(k,j))/tdlx3
      uxxxx = (-2.0*u(k-5,j)+11.0*u(k-4,j)-24.0*u(k-3,j)+26.0*u(k-2,j)- &
                                         14.0*u(k-1,j)+3.0*u(k,j))/dlx4
      return
      end
      subroutine dx4(u,idmn,i,j,uxxx,uxxxx)

c*********************************************************************72
c
cc DX4 computes second order finite difference
c     approximations to the third and fourth x
c     partial derivatives of u at the (i,j) mesh point
c
      common /spl4/  kswx       ,kswy       ,k          ,l          , &
                     ait        ,bit        ,cit        ,dit        , &
                     mit        ,nit        ,is         ,ms         , &
                     js         ,ns         ,dlx        ,dly        , &
                     tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: spl4
!$omp threadprivate (spl4)
      dimension       u(idmn,*)

      if (i.gt.2 .and. i.lt.(k-1)) go to  50
      if (i .eq. 1) go to  10
      if (i .eq. 2) go to  30
      if (i .eq. k-1) go to  60
      if (i .eq. k) go to  80
c
c     compute partial derivative approximations at x=a
c
   10 if (kswx .eq. 1) go to  20
      uxxx = (-5.0*u(1,j)+18.0*u(2,j)-24.0*u(3,j)+14.0*u(4,j)-
     &                                               3.0*u(5,j))/(tdlx3)
      uxxxx = (3.0*u(1,j)-14.0*u(2,j)+26.0*u(3,j)-24.0*u(4,j)+
     &                                      11.0*u(5,j)-2.0*u(6,j))/dlx4
      return
c
c     periodic at x=a
c
   20 uxxx = (-u(k-2,j)+2.0*u(k-1,j)-2.0*u(2,j)+u(3,j))/(tdlx3)
      uxxxx = (u(k-2,j)-4.0*u(k-1,j)+6.0*u(1,j)-4.0*u(2,j)+u(3,j))/dlx4
      return
c
c     compute partial derivative approximations at x=a+dlx
c
   30 if (kswx .eq. 1) go to  40
      uxxx = (-3.0*u(1,j)+10.0*u(2,j)-12.0*u(3,j)+6.0*u(4,j)-u(5,j))/
     &       tdlx3
      uxxxx = (2.0*u(1,j)-9.0*u(2,j)+16.0*u(3,j)-14.0*u(4,j)+6.0*u(5,j)-
     &                                                      u(6,j))/dlx4
      return
c
c     periodic at x=a+dlx
c
   40 uxxx = (-u(k-1,j)+2.0*u(1,j)-2.0*u(3,j)+u(4,j))/(tdlx3)
      uxxxx = (u(k-1,j)-4.0*u(1,j)+6.0*u(2,j)-4.0*u(3,j)+u(4,j))/dlx4
      return
c
c     compute partial derivative approximations on the interior
c
   50 continue
      uxxx = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
      uxxxx = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))/
     &        dlx4
      return
c
c     compute partial derivative approximations at x=b-dlx
c
   60 if (kswx .eq. 1) go to  70
      uxxx = (u(k-4,j)-6.0*u(k-3,j)+12.0*u(k-2,j)-10.0*u(k-1,j)+
     &                                                 3.0*u(k,j))/tdlx3
      uxxxx = (-u(k-5,j)+6.0*u(k-4,j)-14.0*u(k-3,j)+16.0*u(k-2,j)-
     &                                     9.0*u(k-1,j)+2.0*u(k,j))/dlx4
      return
c
c     periodic at x=b-dlx
c
   70 uxxx = (-u(k-3,j)+2.0*u(k-2,j)-2.0*u(1,j)+u(2,j))/tdlx3
      uxxxx = (u(k-3,j)-4.0*u(k-2,j)+6.0*u(k-1,j)-4.0*u(1,j)+u(2,j))/
     &        dlx4
      return
c
c     compute partial derivative approximations at x=b
c
   80 uxxx = -(3.0*u(k-4,j)-14.0*u(k-3,j)+24.0*u(k-2,j)-18.0*u(k-1,j)+
     &                                                 5.0*u(k,j))/tdlx3
      uxxxx = (-2.0*u(k-5,j)+11.0*u(k-4,j)-24.0*u(k-3,j)+26.0*u(k-2,j)-
     &                                    14.0*u(k-1,j)+3.0*u(k,j))/dlx4
      return
      end
      subroutine dy (u,idmn,i,j,uyyy,uyyyy)

c*********************************************************************72
c
cc DY computes second order finite difference
c     approximations to the third and fourth y
c     partial derivatives of u at the (i,j) mesh point
c
      common /splp/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: splp
!$omp threadprivate (splp)
      dimension       u(idmn,*)

      if (j.gt.2 .and. j.lt.(l-1)) go to  50
      if (j .eq. 1) go to  10
      if (j .eq. 2) go to  30
      if (j .eq. l-1) go to  60
      if (j .eq. l) go to  80
c
c     compute partial derivative approximations at y=c
c
   10 if (kswy .eq. 1) go to  20
      uyyy = (-5.0*u(i,1)+18.0*u(i,2)-24.0*u(i,3)+14.0*u(i,4)-
     &                                                 3.0*u(i,5))/tdly3
      uyyyy = (3.0*u(i,1)-14.0*u(i,2)+26.0*u(i,3)-24.0*u(i,4)+
     &                                      11.0*u(i,5)-2.0*u(i,6))/dly4
      return
c
c     periodic at x=a
c
   20 uyyy = (-u(i,l-2)+2.0*u(i,l-1)-2.0*u(i,2)+u(i,3))/tdly3
      uyyyy = (u(i,l-2)-4.0*u(i,l-1)+6.0*u(i,1)-4.0*u(i,2)+u(i,3))/dly4
      return
c
c     compute partial derivative approximations at y=c+dly
c
   30 if (kswy .eq. 1) go to  40
      uyyy = (-3.0*u(i,1)+10.0*u(i,2)-12.0*u(i,3)+6.0*u(i,4)-u(i,5))/
     &       tdly3
      uyyyy = (2.0*u(i,1)-9.0*u(i,2)+16.0*u(i,3)-14.0*u(i,4)+6.0*u(i,5)-
     &                                                      u(i,6))/dly4
      return
c
c     periodic at y=c+dly
c
   40 uyyy = (-u(i,l-1)+2.0*u(i,1)-2.0*u(i,3)+u(i,4))/tdly3
      uyyyy = (u(i,l-1)-4.0*u(i,1)+6.0*u(i,2)-4.0*u(i,3)+u(i,4))/dly4
      return
c
c     compute partial derivative approximations on the interior
c
   50 continue
      uyyy = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
      uyyyy = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))/
     &        dly4
      return
c
c     compute partial derivative approximations at y=d-dly
c
   60 if (kswy .eq. 1) go to  70
      uyyy = (u(i,l-4)-6.0*u(i,l-3)+12.0*u(i,l-2)-10.0*u(i,l-1)+
     &                                                 3.0*u(i,l))/tdly3
      uyyyy = (-u(i,l-5)+6.0*u(i,l-4)-14.0*u(i,l-3)+16.0*u(i,l-2)-
     &                                     9.0*u(i,l-1)+2.0*u(i,l))/dly4
      return
c
c     periodic at y=d-dly
c
   70 continue
      uyyy = (-u(i,l-3)+2.0*u(i,l-2)-2.0*u(i,1)+u(i,2))/tdly3
      uyyyy = (u(i,l-3)-4.0*u(i,l-2)+6.0*u(i,l-1)-4.0*u(i,1)+u(i,2))/
     &        dly4
      return
c
c     compute partial derivative approximations at y=d
c
   80 uyyy = -(3.0*u(i,l-4)-14.0*u(i,l-3)+24.0*u(i,l-2)-18.0*u(i,l-1)+
     &                                                 5.0*u(i,l))/tdly3
      uyyyy = (-2.0*u(i,l-5)+11.0*u(i,l-4)-24.0*u(i,l-3)+26.0*u(i,l-2)-
     &                                    14.0*u(i,l-1)+3.0*u(i,l))/dly4
      return
      end
      subroutine dy4(u,idmn,i,j,uyyy,uyyyy)

c*********************************************************************72
c
cc DY4 computes second order finite difference
c     approximations to the third and fourth y
c     partial derivatives of u at the (i,j) mesh point
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: spl4
!$omp threadprivate(spl4)
      dimension       u(idmn,*)

      if (j.gt.2 .and. j.lt.(l-1)) go to  50
      if (j .eq. 1) go to  10
      if (j .eq. 2) go to  30
      if (j .eq. l-1) go to  60
      if (j .eq. l) go to  80
c
c     compute partial derivative approximations at y=c
c
   10 if (kswy .eq. 1) go to  20
      uyyy = (-5.0*u(i,1)+18.0*u(i,2)-24.0*u(i,3)+14.0*u(i,4)-
     &                                                 3.0*u(i,5))/tdly3
      uyyyy = (3.0*u(i,1)-14.0*u(i,2)+26.0*u(i,3)-24.0*u(i,4)+
     &                                      11.0*u(i,5)-2.0*u(i,6))/dly4
      return
c
c     periodic at x=a
c
   20 uyyy = (-u(i,l-2)+2.0*u(i,l-1)-2.0*u(i,2)+u(i,3))/tdly3
      uyyyy = (u(i,l-2)-4.0*u(i,l-1)+6.0*u(i,1)-4.0*u(i,2)+u(i,3))/dly4
      return
c
c     compute partial derivative approximations at y=c+dly
c
   30 if (kswy .eq. 1) go to  40
      uyyy = (-3.0*u(i,1)+10.0*u(i,2)-12.0*u(i,3)+6.0*u(i,4)-u(i,5))/
     &       tdly3
      uyyyy = (2.0*u(i,1)-9.0*u(i,2)+16.0*u(i,3)-14.0*u(i,4)+6.0*u(i,5)-
     &                                                      u(i,6))/dly4
      return
c
c     periodic at y=c+dly
c
   40 uyyy = (-u(i,l-1)+2.0*u(i,1)-2.0*u(i,3)+u(i,4))/tdly3
      uyyyy = (u(i,l-1)-4.0*u(i,1)+6.0*u(i,2)-4.0*u(i,3)+u(i,4))/dly4
      return
c
c     compute partial derivative approximations on the interior
c
   50 continue
      uyyy = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
      uyyyy = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))/
     &        dly4
      return
c
c     compute partial derivative approximations at y=d-dly
c
   60 if (kswy .eq. 1) go to  70
      uyyy = (u(i,l-4)-6.0*u(i,l-3)+12.0*u(i,l-2)-10.0*u(i,l-1)+
     &                                                 3.0*u(i,l))/tdly3
      uyyyy = (-u(i,l-5)+6.0*u(i,l-4)-14.0*u(i,l-3)+16.0*u(i,l-2)-
     &                                     9.0*u(i,l-1)+2.0*u(i,l))/dly4
      return
c
c     periodic at y=d-dly
c
   70 continue
      uyyy = (-u(i,l-3)+2.0*u(i,l-2)-2.0*u(i,1)+u(i,2))/tdly3
      uyyyy = (u(i,l-3)-4.0*u(i,l-2)+6.0*u(i,l-1)-4.0*u(i,1)+u(i,2))/
     &        dly4
      return
c
c     compute partial derivative approximations at y=d
c
   80 uyyy = -(3.0*u(i,l-4)-14.0*u(i,l-3)+24.0*u(i,l-2)-18.0*u(i,l-1)+
     &                                                 5.0*u(i,l))/tdly3
      uyyyy = (-2.0*u(i,l-5)+11.0*u(i,l-4)-24.0*u(i,l-3)+26.0*u(i,l-2)-
     &                                    14.0*u(i,l-1)+3.0*u(i,l))/dly4
      return
      end
      function epmach ()

c*********************************************************************72
c
cc EPMACH computes an approximate machiine epsilon (accuracy)
c
      common /value/  v
      save :: v
!$omp threadprivate(value)
      eps = 1.
  101 eps = eps/10.
      call store (eps+1.)
      if (v-1.) 102,102,101
  102 epmach = 100.*eps
      return
      end
      subroutine fdump

c*********************************************************************72
c
cc FDUMP creates an error dump.
c
c***BEGIN PROLOGUE  FDUMP
c***PURPOSE  Symbolic dump (should be locally written).
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3
c***TYPE      ALL (FDUMP-A)
c***KEYWORDS  ERROR, XERMSG
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c        ***Note*** Machine Dependent Routine
c        FDUMP is intended to be replaced by a locally written
c        version which produces a symbolic dump.  Failing this,
c        it should be replaced by a version which prints the
c        subprogram nesting list.  Note that this dump must be
c        printed on each of up to five files, as indicated by the
c        XGETUA routine.  See XSETUA and XGETUA for details.
c
c     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
c
c***REFERENCES  (NONE)
c***ROUTINES CALLED  (NONE)
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c***END PROLOGUE  FDUMP
c***FIRST EXECUTABLE STATEMENT  FDUMP
      return
      end
      subroutine genbun (nperod,n,mperod,m,a,b,c,idimy,y,ierror,w)

c*********************************************************************72
c
cc GENBUN: generalized Buneman algorithm, linear equation solver.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c
c   genbun solves the linear system of equations
c
c          a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c
c          + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.,
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     0, x(i,2), or x(i,n) and x(i,n+1) may be equal to 0, x(i,n-1), or
c     x(i,1) depending on an input parameter.
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     nperod
c       indicates the values that x(i,0) and x(i,n+1) are assumed to
c       have.
c
c       = 0  if x(i,0) = x(i,n) and x(i,n+1) = x(i,1).
c       = 1  if x(i,0) = x(i,n+1) = 0  .
c       = 2  if x(i,0) = 0 and x(i,n+1) = x(i,n-1).
c       = 3  if x(i,0) = x(i,2) and x(i,n+1) = x(i,n-1).
c       = 4  if x(i,0) = x(i,2) and x(i,n+1) = 0.
c
c     n
c       the number of unknowns in the j-direction.  n must be greater
c       than 2.
c
c     mperod
c       = 0 if a(1) and c(m) are not zero
c       = 1 if a(1) = c(m) = 0
c
c     m
c       the number of unknowns in the i-direction.  m must be greater
c       than 2.
c
c     a,b,c
c       one-dimensional arrays of length m that specify the
c       coefficients in the linear equations given above.  if mperod = 0
c       the array elements must not depend upon the index i, but must be
c       constant.  specifically, the subroutine checks the following
c       condition
c
c             a(i) = c(1)
c             c(i) = c(1)
c             b(i) = b(1)
c
c       for i=1,2,...,m.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling genbun.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c       a two-dimensional array that specifies the values of the right
c       side of the linear system of equations given above.  y must be
c       dimensioned at least m*n.
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*n + (10 + int(log2(n)))*m
c       locations.  the actual number of locations used is computed by
c       genbun and is returned in location w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag that indicates invalid input parameters  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m .le. 2  .
c       = 2  n .le. 2
c       = 3  idimy .lt. m
c       = 4  nperod .lt. 0 or nperod .gt. 4
c       = 5  mperod .lt. 0 or mperod .gt. 1
c       = 6  a(i) .ne. c(1) or c(i) .ne. c(1) or b(i) .ne. b(1) for
c            some i=1,2,...,m.
c       = 7  a(1) .ne. 0 or c(m) .ne. 0 and mperod = 1
c
c     w
c       w(1) contains the required length of w.
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(m),b(m),c(m),y(idimy,n),w(see parameter list)
c     arguments
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    genbun,poisd2,poisn2,poisp2,cosgen,merge,trix,tri3,
c     required       pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized april 1, 1973
c                    revised august 20,1973
c                    revised january 1, 1976
c
c     algorithm      the linear system is solved by a cyclic reduction
c                    algorithm described in the reference.
c
c     space          4944(decimal) = 11520(octal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine genbun is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameter nperod.  some typical values are listed
c                    in the table below.  more comprehensive timing
c                    charts may be found in the reference.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(i) = c(i) = -0.5*b(i) = 1,       i=1,2,...,m
c
c                    and, when mperod = 1
c
c                       a(1) = c(m) = 0
c                       a(m) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine genbun was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                       e = max(abs(z(i,j)-x(i,j)))/max(abs(x(i,j)))
c
c                    where the two maxima are taken over all i=1,2,...,m
c                    and j=1,2,...,n, was computed.  the value of e is
c                    given in the table below for some typical values of
c                    m and n.
c
c
c                       m (=n)    mperod    nperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         31        0         0          36     6.e-14
c                         31        1         1          21     4.e-13
c                         31        1         3          41     3.e-13
c                         32        0         0          29     9.e-14
c                         32        1         1          32     3.e-13
c                         32        1         3          48     1.e-13
c                         33        0         0          36     9.e-14
c                         33        1         1          30     4.e-13
c                         33        1         3          34     1.e-13
c                         63        0         0         150     1.e-13
c                         63        1         1          91     1.e-12
c                         63        1         3         173     2.e-13
c                         64        0         0         122     1.e-13
c                         64        1         1         128     1.e-12
c                         64        1         3         199     6.e-13
c                         65        0         0         143     2.e-13
c                         65        1         1         120     1.e-12
c                         65        1         3         138     4.e-13
c
c     portability    american national standards institue fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      sweet, r., 'a cyclic reduction algorithm for
c                    solving block tridiagonal systems of arbitrary
c                    dimensions,' siam j. on numer. anal.,
c                    14(sept., 1977), pp. 706-720.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
      dimension       y(idimy,1)
      dimension       w(1)       ,b(1)       ,a(1)       ,c(1)
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.0 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 102
      do 101 i=2,m
         if (a(i) .ne. c(1)) go to 103
         if (c(i) .ne. c(1)) go to 103
         if (b(i) .ne. b(1)) go to 103
  101 continue
      go to 104
  102 if (a(1).ne.0. .or. c(m).ne.0.) ierror = 7
      go to 104
  103 ierror = 6
  104 if (ierror .ne. 0) return
      mp1 = m+1
      iwba = mp1
      iwbb = iwba+m
      iwbc = iwbb+m
      iwb2 = iwbc+m
      iwb3 = iwb2+m
      iww1 = iwb3+m
      iww2 = iww1+m
      iww3 = iww2+m
      iwd = iww3+m
      iwtcos = iwd+m
      iwp = iwtcos+4*n
      do 106 i=1,m
         k = iwba+i-1
         w(k) = -a(i)
         k = iwbc+i-1
         w(k) = -c(i)
         k = iwbb+i-1
         w(k) = 2.-b(i)
         do 105 j=1,n
            y(i,j) = -y(i,j)
  105    continue
  106 continue
      mp = mperod+1
      np = nperod+1
      go to (114,107),mp
  107 go to (108,109,110,111,123),np
  108 call poisp2 (m,n,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     &             w(iwp))
      go to 112
  109 call poisd2 (m,n,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iww1),
     &             w(iwd),w(iwtcos),w(iwp))
      go to 112
  110 call poisn2 (m,n,1,2,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     &             w(iwp))
      go to 112
  111 call poisn2 (m,n,1,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     &             w(iwp))
  112 ipstor = w(iww1)
      irev = 2
      if (nperod .eq. 4) go to 124
  113 go to (127,133),mp
  114 continue
c
c     reorder unknowns when mp =0
c
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 119 j=1,n
         do 115 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  115    continue
         w(mh) = 2.*y(mh,j)
         go to (117,116),modd
  116    w(m) = 2.*y(m,j)
  117    continue
         do 118 i=1,m
            y(i,j) = w(i)
  118    continue
  119 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = 0.
      w(i) = 0.
      w(k+1) = 2.*w(k+1)
      go to (120,121),modd
  120 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 122
  121 w(iwbb-1) = w(k+1)
  122 continue
      go to 107
c
c     reverse columns when nperod = 4.
c
  123 irev = 1
      nby2 = n/2
  124 do 126 j=1,nby2
         mskip = n+1-j
         do 125 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  125    continue
  126 continue
      go to (110,113),irev
  127 continue
      do 132 j=1,n
         do 128 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5*(y(mhpi,j)-y(i,j))
  128    continue
         w(mh) = .5*y(mh,j)
         go to (130,129),modd
  129    w(m) = .5*y(m,j)
  130    continue
         do 131 i=1,m
            y(i,j) = w(i)
  131    continue
  132 continue
  133 continue
c
c     return storage requirements for w array.
c
      w(1) = ipstor+iwp-1
      return
      end
      subroutine hstcrt (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HSTCRT: 2D staggered grid Cartesian coordinates five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c      hstcrt solves the standard five-point finite difference
c      approximation on a staggered grid to the helmholtz equation in
c      cartesian coordinates
c
c      (d/dx)(du/dx) + (d/dy)(du/dy) + lambda*u = f(x,y)
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c    a,b
c      the range of x, i.e. a .le. x .le. b.  a must be less than b.
c
c    m
c      the number of grid points in the interval (a,b).  the grid points
c      in the x-direction are given by x(i) = a + (i-0.5)dx for
c      i=1,2,...,m where dx =(b-a)/m.  m must be greater than 2.
c
c    mbdcnd
c      indicates the type of boundary conditions at x = a and x = b.
c
c      = 0  if the solution is periodic in x,
c           u(m+i,j) = u(i,j).
c
c      = 1  if the solution is specified at x = a and x = b.
c
c      = 2  if the solution is specified at x = a and the derivative
c           of the solution with respect to x is specified at x = b.
c
c      = 3  if the derivative of the solution with respect to x is
c           specified at x = a  and x = b.
c
c      = 4  if the derivative of the solution with respect to x is
c           specified at x = a  and the solution is specified at x = b.
c
c    bda
c      a one-dimensional array of length n that specifies the boundary
c      values (if any) of the solution at x = a.  when mbdcnd = 1 or 2,
c
c               bda(j) = u(a,y(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 3 or 4,
c
c               bda(j) = (d/dx)u(a,y(j)) ,    j=1,2,...,n.
c
c    bdb
c      a one-dimensional array of length n that specifies the boundary
c      values of the solution at x = b.  when mbdcnd = 1 or 4
c
c               bdb(j) = u(b,y(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 2 or 3
c
c               bdb(j) = (d/dx)u(b,y(j)) ,    j=1,2,...,n.
c
c    c,d
c      the range of y, i.e. c .le. y .le. d.  c must be less
c      than d.
c
c    n
c      the number of unknowns in the interval (c,d).  the unknowns in
c      the y-direction are given by y(j) = c + (j-0.5)dy,
c      j=1,2,...,n, where dy = (d-c)/n.  n must be greater than 2.
c
c    nbdcnd
c      indicates the type of boundary conditions at y = c
c      and y = d.
c
c      = 0  if the solution is periodic in y, i.e.
c           u(i,j) = u(i,n+j).
c
c      = 1  if the solution is specified at y = c and y = d.
c
c      = 2  if the solution is specified at y = c and the derivative
c           of the solution with respect to y is specified at y = d.
c
c      = 3  if the derivative of the solution with respect to y is
c           specified at y = c and y = d.
c
c      = 4  if the derivative of the solution with respect to y is
c           specified at y = c and the solution is specified at y = d.
c
c    bdc
c      a one dimensional array of length m that specifies the boundary
c      values of the solution at y = c.   when nbdcnd = 1 or 2,
c
c               bdc(i) = u(x(i),c) ,              i=1,2,...,m.
c
c      when nbdcnd = 3 or 4,
c
c               bdc(i) = (d/dy)u(x(i),c),     i=1,2,...,m.
c
c      when nbdcnd = 0, bdc is a dummy variable.
c
c    bdd
c      a one-dimensional array of length m that specifies the boundary
c      values of the solution at y = d.  when nbdcnd = 1 or 4,
c
c               bdd(i) = u(x(i),d) ,              i=1,2,...,m.
c
c      when nbdcnd = 2 or 3,
c
c               bdd(i) = (d/dy)u(x(i),d) ,    i=1,2,...,m.
c
c      when nbdcnd = 0, bdd is a dummy variable.
c
c    elmbda
c      the constant lambda in the helmholtz equation.  if lambda is
c      greater than 0, a solution may not exist.  however, hstcrt will
c      attempt to find a solution.
c
c    f
c      a two-dimensional array that specifies the values of the right
c      side of the helmholtz equation.  for i=1,2,...,m and j=1,2,...,n
c
c               f(i,j) = f(x(i),y(j)) .
c
c      f must be dimensioned at least m x n.
c
c    idimf
c      the row (or first) dimension of the array f as it appears in the
c      program calling hstcrt.  this parameter is used to specify the
c      variable dimension of f.  idimf must be at least m.
c
c    w
c      a one-dimensional array that must be provided by the user for
c      work space.  w may require up to 13m + 4n + m*int(log2(n))
c      locations.  the actual number of locations used is computed by
c      hstcrt and is returned in the location w(1).
c
c
c             * * * * * *   on output   * * * * * *
c
c    f
c      contains the solution u(i,j) of the finite difference
c      approximation for the grid point (x(i),y(j)) for
c      i=1,2,...,m, j=1,2,...,n.
c
c    pertrb
c      if a combination of periodic or derivative boundary conditions is
c      specified for a poisson equation (lambda = 0), a solution may not
c      exist.  pertrb is a constant, calculated and subtracted from f,
c      which ensures that a solution exists.  hstcrt then computes this
c      solution, which is a least squares solution to the original
c      approximation.  this solution plus any constant is also a
c      solution; hence, the solution is not unique.  the value of pertrb
c      should be small compared to the right side f.  otherwise, a
c      solution is obtained to an essentially different problem.  this
c      comparison should always be made to insure that a meaningful
c      solution has been obtained.
c
c    ierror
c      an error flag that indicates invalid input parameters.
c      except to numbers 0 and  6, a solution is not attempted.
c
c      =  0  no error
c
c      =  1  a .ge. b
c
c      =  2  mbdcnd .lt. 0 or mbdcnd .gt. 4
c
c      =  3  c .ge. d
c
c      =  4  n .le. 2
c
c      =  5  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c      =  6  lambda .gt. 0
c
c      =  7  idimf .lt. m
c
c      =  8  m .le. 2
c
c      since this is the only means of indicating a possibly
c      incorrect call to hstcrt, the user should test ierror after
c      the call.
c
c    w
c      w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c     arguments      w(see argument list)
c
c     latest         june 1, 1977
c     revision
c
c     subprograms    hstcrt,poistg,postg2,genbun,poisd2,poisn2,poisp2,
c     required       cosgen,merge,trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in january , 1977
c
c     algorithm      this subroutine defines the finite-difference
c                    equations, incorporates boundary data, adjusts the
c                    right side when the system is singular and calls
c                    either poistg or genbun which solves the linear
c                    system of equations.
c
c     space          8131(decimal) = 17703(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstcrt is roughly proportional
c                    to m*n*log2(n).  some typical values are listed in
c                    the table below.
c                       the solution process employed results in a loss
c                    of no more than four significant digits for n and m
c                    as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    poistg which is the routine that
c                    actually solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32       1-4       1-4         56
c                        64       1-4       1-4        230
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      schumann, u. and r. sweet,"a direct method for
c                    the solution of poisson"s equation with neumann
c                    boundary conditions on a staggered grid of
c                    arbitrary size," j. comp. phys. 20(1976),
c                    pp. 171-182.
c
      dimension       f(idimf,1) ,bda(1)     ,bdb(1)     ,bdc(1)     ,
     &                bdd(1)     ,w(1)
c
c     check for invalid parameters.
c
      ierror = 0
      if (a .ge. b) ierror = 1
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 2
      if (c .ge. d) ierror = 3
      if (n .le. 2) ierror = 4
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 5
      if (idimf .lt. m) ierror = 7
      if (m .le. 2) ierror = 8
      if (ierror .ne. 0) return
      nperod = nbdcnd
      mperod = 0
      if (mbdcnd .gt. 0) mperod = 1
      deltax = (b-a)/float(m)
      twdelx = 1./deltax
      delxsq = 2./deltax**2
      deltay = (d-c)/float(n)
      twdely = 1./deltay
      delysq = deltay**2
      twdysq = 2./delysq
      np = nbdcnd+1
      mp = mbdcnd+1
c
c     define the a,b,c coefficients in w-array.
c
      id2 = m
      id3 = id2+m
      id4 = id3+m
      s = (deltay/deltax)**2
      st2 = 2.*s
      do 101 i=1,m
         w(i) = s
         j = id2+i
         w(j) = -st2+elmbda*delysq
         j = id3+i
         w(j) = s
  101 continue
c
c     enter boundary data for x-boundaries.
c
      go to (111,102,102,104,104),mp
  102 do 103 j=1,n
         f(1,j) = f(1,j)-bda(j)*delxsq
  103 continue
      w(id2+1) = w(id2+1)-w(1)
      go to 106
  104 do 105 j=1,n
         f(1,j) = f(1,j)+bda(j)*twdelx
  105 continue
      w(id2+1) = w(id2+1)+w(1)
  106 go to (111,107,109,109,107),mp
  107 do 108 j=1,n
         f(m,j) = f(m,j)-bdb(j)*delxsq
  108 continue
      w(id3) = w(id3)-w(1)
      go to 111
  109 do 110 j=1,n
         f(m,j) = f(m,j)-bdb(j)*twdelx
  110 continue
      w(id3) = w(id3)+w(1)
  111 continue
c
c     enter boundary data for y-boundaries.
c
      go to (121,112,112,114,114),np
  112 do 113 i=1,m
         f(i,1) = f(i,1)-bdc(i)*twdysq
  113 continue
      go to 116
  114 do 115 i=1,m
         f(i,1) = f(i,1)+bdc(i)*twdely
  115 continue
  116 go to (121,117,119,119,117),np
  117 do 118 i=1,m
         f(i,n) = f(i,n)-bdd(i)*twdysq
  118 continue
      go to 121
  119 do 120 i=1,m
         f(i,n) = f(i,n)-bdd(i)*twdely
  120 continue
  121 continue
      do 123 i=1,m
         do 122 j=1,n
            f(i,j) = f(i,j)*delysq
  122    continue
  123 continue
      if (mperod .eq. 0) go to 124
      w(1) = 0.
      w(id4) = 0.
  124 continue
      pertrb = 0.
      if (elmbda) 133,126,125
  125 ierror = 6
      go to 133
  126 go to (127,133,133,127,133),mp
  127 go to (128,133,133,128,133),np
c
c     for singular problems must adjust data to insure that a solution
c     will exist.
c
  128 continue
      s = 0.
      do 130 j=1,n
         do 129 i=1,m
            s = s+f(i,j)
  129    continue
  130 continue
      pertrb = s/float(m*n)
      do 132 j=1,n
         do 131 i=1,m
            f(i,j) = f(i,j)-pertrb
  131    continue
  132 continue
      pertrb = pertrb/delysq
c
c     solve the equation.
c
  133 continue
      if (nperod .eq. 0) go to 134
      call poistg (nperod,n,mperod,m,w(1),w(id2+1),w(id3+1),idimf,f,
     &             ierr1,w(id4+1))
      go to 135
  134 continue
      call genbun (nperod,n,mperod,m,w(1),w(id2+1),w(id3+1),idimf,f,
     &             ierr1,w(id4+1))
  135 continue
      w(1) = w(id4+1)+3.*float(m)
      return
      end
      subroutine hstcs1 (intl,a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,
     &                   bdd,elmbda,f,idimf,pertrb,ierr1,am,bm,cm,an,
     &                   bn,cn,snth,rsq,wrk)

c*********************************************************************72
c
cc HSTCS1 is a utility routine for HSTCSP.
c
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     &                f(idimf,1) ,am(1)      ,bm(1)      ,cm(1)      ,
     &                an(1)      ,bn(1)      ,cn(1)      ,snth(1)    ,
     &                rsq(1)     ,wrk(1)
      dth = (b-a)/float(m)
      dthsq = dth*dth
      do 101 i=1,m
         snth(i) = sin(a+(float(i)-0.5)*dth)
  101 continue
      dr = (d-c)/float(n)
      do 102 j=1,n
         rsq(j) = (c+(float(j)-0.5)*dr)**2
  102 continue
c
c     multiply right side by r(j)**2
c
      do 104 j=1,n
         x = rsq(j)
         do 103 i=1,m
            f(i,j) = x*f(i,j)
  103    continue
  104 continue
c
c      define coefficients am,bm,cm
c
      x = 1./(2.*cos(dth/2.))
      do 105 i=2,m
         am(i) = (snth(i-1)+snth(i))*x
         cm(i-1) = am(i)
  105 continue
      am(1) = sin(a)
      cm(m) = sin(b)
      do 106 i=1,m
         x = 1./snth(i)
         y = x/dthsq
         am(i) = am(i)*y
         cm(i) = cm(i)*y
         bm(i) = elmbda*x*x-am(i)-cm(i)
  106 continue
c
c     define coefficients an,bn,cn
c
      x = c/dr
      do 107 j=1,n
         an(j) = (x+float(j-1))**2
         cn(j) = (x+float(j))**2
         bn(j) = -(an(j)+cn(j))
  107 continue
      isw = 1
      nb = nbdcnd
      if (c.eq.0. .and. nb.eq.2) nb = 6
c
c     enter data on theta boundaries
c
      go to (108,108,110,110,112,112,108,110,112),mbdcnd
  108 bm(1) = bm(1)-am(1)
      x = 2.*am(1)
      do 109 j=1,n
         f(1,j) = f(1,j)-x*bda(j)
  109 continue
      go to 112
  110 bm(1) = bm(1)+am(1)
      x = dth*am(1)
      do 111 j=1,n
         f(1,j) = f(1,j)+x*bda(j)
  111 continue
  112 continue
      go to (113,115,115,113,113,115,117,117,117),mbdcnd
  113 bm(m) = bm(m)-cm(m)
      x = 2.*cm(m)
      do 114 j=1,n
         f(m,j) = f(m,j)-x*bdb(j)
  114 continue
      go to 117
  115 bm(m) = bm(m)+cm(m)
      x = dth*cm(m)
      do 116 j=1,n
         f(m,j) = f(m,j)-x*bdb(j)
  116 continue
  117 continue
c
c     enter data on r boundaries
c
      go to (118,118,120,120,122,122),nb
  118 bn(1) = bn(1)-an(1)
      x = 2.*an(1)
      do 119 i=1,m
         f(i,1) = f(i,1)-x*bdc(i)
  119 continue
      go to 122
  120 bn(1) = bn(1)+an(1)
      x = dr*an(1)
      do 121 i=1,m
         f(i,1) = f(i,1)+x*bdc(i)
  121 continue
  122 continue
      go to (123,125,125,123,123,125),nb
  123 bn(n) = bn(n)-cn(n)
      x = 2.*cn(n)
      do 124 i=1,m
         f(i,n) = f(i,n)-x*bdd(i)
  124 continue
      go to 127
  125 bn(n) = bn(n)+cn(n)
      x = dr*cn(n)
      do 126 i=1,m
         f(i,n) = f(i,n)-x*bdd(i)
  126 continue
  127 continue
c
c     check for singular problem.  if singular, perturb f.
c
      pertrb = 0.
      go to (137,137,128,137,137,128,137,128,128),mbdcnd
  128 go to (137,137,129,137,137,129),nb
  129 if (elmbda) 137,131,130
  130 ierr1 = 10
      go to 137
  131 continue
      isw = 2
      do 133 i=1,m
         x = 0.
         do 132 j=1,n
            x = x+f(i,j)
  132    continue
         pertrb = pertrb+x*snth(i)
  133 continue
      x = 0.
      do 134 j=1,n
         x = x+rsq(j)
  134 continue
      pertrb = 2.*(pertrb*sin(dth/2.))/(x*(cos(a)-cos(b)))
      do 136 j=1,n
         x = rsq(j)*pertrb
         do 135 i=1,m
            f(i,j) = f(i,j)-x
  135    continue
  136 continue
  137 continue
      a2 = 0.
      do 138 i=1,m
         a2 = a2+f(i,1)
  138 continue
      a2 = a2/rsq(1)
c
c  Initialize BLKTRI.
c
      if (intl .ne. 0) go to 139
      call blktri (0,1,n,an,bn,cn,1,m,am,bm,cm,idimf,f,ierr1,wrk)
  139 continue
c
c  call BLKTRI to solve system of equations.
c
      call blktri (1,1,n,an,bn,cn,1,m,am,bm,cm,idimf,f,ierr1,wrk)
      if (isw.ne.2 .or. c.ne.0. .or. nbdcnd.ne.2) go to 143
      a1 = 0.
      a3 = 0.
      do 140 i=1,m
         a1 = a1+snth(i)*f(i,1)
         a3 = a3+snth(i)
  140 continue
      a1 = a1+rsq(1)*a2/2.
      if (mbdcnd .eq. 3)
     &    a1 = a1+(sin(b)*bdb(1)-sin(a)*bda(1))/(2.*(b-a))
      a1 = a1/a3
      a1 = bdc(1)-a1
      do 142 i=1,m
         do 141 j=1,n
            f(i,j) = f(i,j)+a1
  141    continue
  142 continue
  143 continue
      return
      end
      subroutine hstcsp (intl,a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,
     &                   bdd,elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HSTCSP 2D staggered grid axisymmetric spherical coordinates, five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c     hstcsp solves the standard five-point finite difference
c     approximation on a staggered grid to the modified helmholtz equation
c     in spherical coordinates assuming axisymmetry (no dependence on
c     longitude)
c
c                  (1/r**2)(d/dr)(r**2(du/dr)) +
c
c       1/(r**2*sin(theta))(d/dtheta)(sin(theta)(du/dtheta)) +
c
c            (lambda/(r*sin(theta))**2)u  =  f(theta,r)
c
c     where theta is colatitude and r is the radial coordinate.
c     this two-dimensional modified helmholtz equation results from
c     the fourier transform of the three-dimensional poisson equation.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c   intl
c     = 0  on initial entry to hstcsp or if any of the arguments
c          c, d, n, or nbdcnd are changed from a previous call
c
c     = 1  if c, d, n, and nbdcnd are all unchanged from previous
c          call to hstcsp
c
c     note:  a call with intl = 0 takes approximately 1.5 times as much
c            time as a call with intl = 1.  once a call with intl = 0
c            has been made then subsequent solutions corresponding to
c            different f, bda, bdb, bdc, and bdd can be obtained
c            faster with intl = 1 since initialization is not repeated.
c
c   a,b
c     the range of theta (colatitude), i.e. a .le. theta .le. b.  a
c     must be less than b and a must be non-negative.  a and b are in
c     radians.  a = 0 corresponds to the north pole and b = pi
c     corresponds to the south pole.
c
c                  * * *  important  * * *
c
c     if b is equal to pi, then b must be computed using the statement
c
c     b = pimach()
c
c     this insures that b in the user"s program is equal to pi in this
c     program which permits several tests of the input parameters that
c     otherwise would not be possible.
c
c                  * * * * * * * * * * * *
c
c   m
c     the number of grid points in the interval (a,b).  the grid points
c     in the theta-direction are given by theta(i) = a + (i-0.5)dtheta
c     for i=1,2,...,m where dtheta =(b-a)/m.  m must be greater than 4.
c
c   mbdcnd
c     indicates the type of boundary conditions at theta = a and
c     theta = b.
c
c     = 1  if the solution is specified at theta = a and theta = b.
c          (see notes 1, 2 below)
c
c     = 2  if the solution is specified at theta = a and the derivative
c          of the solution with respect to theta is specified at
c          theta = b (see notes 1, 2 below).
c
c     = 3  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1, 2 below) and theta = b.
c
c     = 4  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1, 2 below) and the
c          solution is specified at theta = b.
c
c     = 5  if the solution is unspecified at theta = a = 0 and the
c          solution is specified at theta = b. (see note 2 below)
c
c     = 6  if the solution is unspecified at theta = a = 0 and the
c          derivative of the solution with respect to theta is
c          specified at theta = b (see note 2 below).
c
c     = 7  if the solution is specified at theta = a and the
c          solution is unspecified at theta = b = pi.
c
c     = 8  if the derivative of the solution with respect to
c          theta is specified at theta = a (see note 1 below)
c          and the solution is unspecified at theta = b = pi.
c
c     = 9  if the solution is unspecified at theta = a = 0 and
c          theta = b = pi.
c
c     notes:  1.  if a = 0, do not use mbdcnd = 1,2,3,4,7 or 8,
c                 but instead use mbdcnd = 5, 6, or 9.
c
c             2.  if b = pi, do not use mbdcnd = 1,2,3,4,5 or 6,
c                 but instead use mbdcnd = 7, 8, or 9.
c
c             when a = 0  and/or b = pi the only meaningful boundary
c             condition is du/dtheta = 0.  (see d. greenspan, 'numerical
c             analysis of elliptic boundary value problems,' harper and
c             row, 1965, chapter 5.)
c
c   bda
c     a one-dimensional array of length n that specifies the boundary
c     values (if any) of the solution at theta = a.  when
c     mbdcnd = 1, 2, or 7,
c
c              bda(j) = u(a,r(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 3, 4, or 8,
c
c              bda(j) = (d/dtheta)u(a,r(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bda is a dummy variable.
c
c   bdb
c     a one-dimensional array of length n that specifies the boundary
c     values of the solution at theta = b.  when mbdcnd = 1, 4, or 5,
c
c              bdb(j) = u(b,r(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 2,3, or 6,
c
c              bdb(j) = (d/dtheta)u(b,r(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bdb is a dummy variable.
c
c   c,d
c     the range of r , i.e. c .le. r .le. d.
c     c must be less than d.  c must be non-negative.
c
c   n
c     the number of unknowns in the interval (c,d).  the unknowns in
c     the r-direction are given by r(j) = c + (j-0.5)dr,
c     j=1,2,...,n, where dr = (d-c)/n.  n must be greater than 4.
c
c   nbdcnd
c     indicates the type of boundary conditions at r = c
c     and r = d.
c
c     = 1  if the solution is specified at r = c and r = d.
c
c     = 2  if the solution is specified at r = c and the derivative
c          of the solution with respect to r is specified at
c          r = d. (see note 1 below)
c
c     = 3  if the derivative of the solution with respect to r is
c          specified at r = c and r = d.
c
c     = 4  if the derivative of the solution with respect to r is
c          specified at r = c and the solution is specified at
c          r = d.
c
c     = 5  if the solution is unspecified at r = c = 0 (see note 2
c          below) and the solution is specified at r = d.
c
c     = 6  if the solution is unspecified at r = c = 0 (see note 2
c          below) and the derivative of the solution with respect to r
c          is specified at r = d.
c
c     note 1:  if c = 0 and mbdcnd = 3,6,8 or 9, the system of equations
c              to be solved is singular.  the unique solution is
c              determined by extrapolation to the specification of
c              u(theta(1),c).  but in these cases the right side of the
c              system will be perturbed by the constant pertrb.
c
c     note 2:  nbdcnd = 5 or 6 cannot be used with mbdcnd = 1, 2, 4, 5,
c              or 7 (the former indicates that the solution is
c              unspecified at r = 0; the latter indicates that the
c              solution is specified).  use instead nbdcnd = 1 or 2.
c
c   bdc
c     a one dimensional array of length m that specifies the boundary
c     values of the solution at r = c.   when nbdcnd = 1 or 2,
c
c              bdc(i) = u(theta(i),c) ,              i=1,2,...,m.
c
c     when nbdcnd = 3 or 4,
c
c              bdc(i) = (d/dr)u(theta(i),c),         i=1,2,...,m.
c
c     when nbdcnd has any other value, bdc is a dummy variable.
c
c   bdd
c     a one-dimensional array of length m that specifies the boundary
c     values of the solution at r = d.  when nbdcnd = 1 or 4,
c
c              bdd(i) = u(theta(i),d) ,              i=1,2,...,m.
c
c     when nbdcnd = 2 or 3,
c
c              bdd(i) = (d/dr)u(theta(i),d) ,        i=1,2,...,m.
c
c     when nbdcnd has any other value, bdd is a dummy variable.
c
c   elmbda
c     the constant lambda in the modified helmholtz equation.  if
c     lambda is greater than 0, a solution may not exist.  however,
c     hstcsp will attempt to find a solution.
c
c   f
c     a two-dimensional array that specifies the values of the right
c     side of the modified helmholtz equation.  for i=1,2,...,m and
c     j=1,2,...,n
c
c              f(i,j) = f(theta(i),r(j)) .
c
c     f must be dimensioned at least m x n.
c
c   idimf
c     the row (or first) dimension of the array f as it appears in the
c     program calling hstcsp.  this parameter is used to specify the
c     variable dimension of f.  idimf must be at least m.
c
c   w
c     a one-dimensional array that must be provided by the user for
c     work space.  with k = int(log2(n))+1 and l = 2**(k+1), w may
c     require up to (k-2)*l+k+max(2n,6m)+4(n+m)+5 locations.  the
c     actual number of locations used is computed by hstcsp and is
c     returned in the location w(1).
c
c
c            * * * * * *   on output   * * * * * *
c
c   f
c     contains the solution u(i,j) of the finite difference
c     approximation for the grid point (theta(i),r(j)) for
c     i=1,2,...,m, j=1,2,...,n.
c
c   pertrb
c     if a combination of periodic, derivative, or unspecified
c     boundary conditions is specified for a poisson equation
c     (lambda = 0), a solution may not exist.  pertrb is a con-
c     stant, calculated and subtracted from f, which ensures
c     that a solution exists.  hstcsp then computes this
c     solution, which is a least squares solution to the
c     original approximation.  this solution plus any constant is also
c     a solution; hence, the solution is not unique.  the value of
c     pertrb should be small compared to the right side f.
c     otherwise, a solution is obtained to an essentially different
c     problem.  this comparison should always be made to insure that
c     a meaningful solution has been obtained.
c
c   ierror
c     an error flag that indicates invalid input parameters.
c     except for numbers 0 and 10, a solution is not attempted.
c
c     =  0  no error
c
c     =  1  a .lt. 0 or b .gt. pi
c
c     =  2  a .ge. b
c
c     =  3  mbdcnd .lt. 1 or mbdcnd .gt. 9
c
c     =  4  c .lt. 0
c
c     =  5  c .ge. d
c
c     =  6  nbdcnd .lt. 1 or nbdcnd .gt. 6
c
c     =  7  n .lt. 5
c
c     =  8  nbdcnd = 5 or 6 and mbdcnd = 1, 2, 4, 5, or 7
c
c     =  9  c .gt. 0 and nbdcnd .ge. 5
c
c     = 10  elmbda .gt. 0
c
c     = 11  idimf .lt. m
c
c     = 12  m .lt. 5
c
c     = 13  a = 0 and mbdcnd =1,2,3,4,7 or 8
c
c     = 14  b = pi and mbdcnd .le. 6
c
c     = 15  a .gt. 0 and mbdcnd = 5, 6, or 9
c
c     = 16  b .lt. pi and mbdcnd .ge. 7
c
c     = 17  lambda .ne. 0 and nbdcnd .ge. 5
c
c     since this is the only means of indicating a possibly
c     incorrect call to hstcsp, the user should test ierror after
c     the call.
c
c   w
c     w(1) contains the required length of w.  also  w contains
c     intermediate values that must not be destroyed if hstcsp
c     will be called again with intl = 1.
c
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c    dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c    arguments      w(see argument list)
c
c    latest         june 1979
c    revision
c
c    subprograms    hstcsp,hstcs1,blktri,blktr1,indxa,indxb,indxc,
c    required       prod,prodp,cprod,cprodp,ppadd,psgf,bsrh,ppsgf,
c                   ppspf,compb,tevls,epmach,store
c
c    special        none
c    conditions
c
c    common         cblkt,valu1
c    blocks
c
c    i/o            none
c
c    precision      single
c
c    specialist     roland sweet
c
c    language       fortran
c
c    history        written by roland sweet at ncar in may, 1977
c
c    algorithm      this subroutine defines the finite-difference
c                   equations, incorporates boundary data, adjusts the
c                   right side when the system is singular and calls
c                   blktri which solves the linear system of equations.
c
c    space          5269(decimal) = 12225(octal) locations on the
c    required       ncar control data 7600
c
c    timing and        the execution time t on the ncar control data
c    accuracy       7600 for subroutine hstcsp is roughly proportional
c                   to m*n*log2(n), but depends on the input parameter
c                   intl.  some values are listed in the table below.
c                      the solution process employed results in a loss
c                   of no more than four significant digits for n and m
c                   as large as 64.  more detailed information about
c                   accuracy can be found in the documentation for
c                   blktri which is the routine that
c                   actually solves the finite difference equations.
c
c
c                      m(=n)     intl      mbdcnd(=nbdcnd)     t(msecs)
c                      -----     ----      ---------------     --------
c
c                       32        0              1-6             132
c                       32        1              1-6              88
c                       64        0              1-6             546
c                       64        1              1-6             380
c
c    portability    american national standards institute fortran.
c                   an approximate machine epsilon is computed in
c                   function pimach.
c
c    required       cos,sin,cabs,csqrt
c    resident
c    routines
c
c    reference      swarztrauber, p.n., "a direct method for the
c                   discrete solution of separable elliptic equations,"
c                   SIAM Journal on Numerical Analysis 11(1974), pp. 1136-1150.
c                   arbitrary size," j. comp. phys. 20(1976),
c                   pp. 171-182.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dimension       f(idimf,1) ,bda(1)     ,bdb(1)     ,bdc(1)     ,
     &                bdd(1)     ,w(1)
      pi = pimach()
c
c     check for invalid input parameters
c
      ierror = 0
      if (a.lt.0. .or. b.gt.pi) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.lt.1 .or. mbdcnd.gt.9) ierror = 3
      if (c .lt. 0.) ierror = 4
      if (c .ge. d) ierror = 5
      if (nbdcnd.lt.1 .or. nbdcnd.gt.6) ierror = 6
      if (n .lt. 5) ierror = 7
      if ((nbdcnd.eq.5 .or. nbdcnd.eq.6) .and. (mbdcnd.eq.1 .or.
     &    mbdcnd.eq.2 .or. mbdcnd.eq.4 .or. mbdcnd.eq.5 .or.
     &                                                     mbdcnd.eq.7))
     &    ierror = 8
      if (c.gt.0. .and. nbdcnd.ge.5) ierror = 9
      if (idimf .lt. m) ierror = 11
      if (m .lt. 5) ierror = 12
      if (a.eq.0. .and. mbdcnd.ne.5 .and. mbdcnd.ne.6 .and. mbdcnd.ne.9)
     &    ierror = 13
      if (b.eq.pi .and. mbdcnd.le.6) ierror = 14
      if (a.gt.0. .and. (mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9))
     &    ierror = 15
      if (b.lt.pi .and. mbdcnd.ge.7) ierror = 16
      if (elmbda.ne.0. .and. nbdcnd.ge.5) ierror = 17
      if (ierror .ne. 0) go to 101
      iwbm = m+1
      iwcm = iwbm+m
      iwan = iwcm+m
      iwbn = iwan+n
      iwcn = iwbn+n
      iwsnth = iwcn+n
      iwrsq = iwsnth+m
      iwwrk = iwrsq+n
      ierr1 = 0
      call hstcs1 (intl,a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &             elmbda,f,idimf,pertrb,ierr1,w,w(iwbm),w(iwcm),
     &             w(iwan),w(iwbn),w(iwcn),w(iwsnth),w(iwrsq),w(iwwrk))
      w(1) = w(iwwrk)+float(iwwrk-1)
      ierror = ierr1
  101 continue
      return
      end
      subroutine hstcyl (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HSTCYL 2D staggered grid cylindrical coordinates five point scheme.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c      hstcyl solves the standard five-point finite difference
c      approximation on a staggered grid to the modified helmholtz
c      equation in cylindrical coordinates
c
c          (1/r)(d/dr)(r(du/dr)) + (d/dz)(du/dz)
c
c                      + lambda*(1/r**2)*u = f(r,z)
c
c      this two-dimensional modified helmholtz equation results
c      from the fourier transform of a three-dimensional poisson
c      equation.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c    a,b
c      the range of r, i.e. a .le. r .le. b.  a must be less than b and
c      a must be non-negative.
c
c    m
c      the number of grid points in the interval (a,b).  the grid points
c      in the r-direction are given by r(i) = a + (i-0.5)dr for
c      i=1,2,...,m where dr =(b-a)/m.  m must be greater than 2.
c
c    mbdcnd
c      indicates the type of boundary conditions at r = a and r = b.
c
c      = 1  if the solution is specified at r = a (see note below) and
c           r = b.
c
c      = 2  if the solution is specified at r = a (see note below) and
c           the derivative of the solution with respect to r is
c           specified at r = b.
c
c      = 3  if the derivative of the solution with respect to r is
c           specified at r = a (see note below) and r = b.
c
c      = 4  if the derivative of the solution with respect to r is
c           specified at r = a (see note below) and the solution is
c           specified at r = b.
c
c      = 5  if the solution is unspecified at r = a = 0 and the solution
c           is specified at r = b.
c
c      = 6  if the solution is unspecified at r = a = 0 and the
c           derivative of the solution with respect to r is specified at
c           r = b.
c
c      note:  if a = 0, do not use mbdcnd = 1,2,3, or 4, but instead
c             use mbdcnd = 5 or 6.  the resulting approximation gives
c             the only meaningful boundary condition, i.e. du/dr = 0.
c             (see d. greenspan, 'introductory numerical analysis of
c             elliptic boundary value problems,' harper and row, 1965,
c             chapter 5.)
c
c    bda
c      a one-dimensional array of length n that specifies the boundary
c      values (if any) of the solution at r = a.  when mbdcnd = 1 or 2,
c
c               bda(j) = u(a,z(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 3 or 4,
c
c               bda(j) = (d/dr)u(a,z(j)) ,    j=1,2,...,n.
c
c      when mbdcnd = 5 or 6, bda is a dummy variable.
c
c    bdb
c      a one-dimensional array of length n that specifies the boundary
c      values of the solution at r = b.  when mbdcnd = 1,4, or 5,
c
c               bdb(j) = u(b,z(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 2,3, or 6,
c
c               bdb(j) = (d/dr)u(b,z(j)) ,    j=1,2,...,n.
c
c    c,d
c      the range of z, i.e. c .le. z .le. d.  c must be less
c      than d.
c
c    n
c      the number of unknowns in the interval (c,d).  the unknowns in
c      the z-direction are given by z(j) = c + (j-0.5)dz,
c      j=1,2,...,n, where dz = (d-c)/n.  n must be greater than 2.
c
c    nbdcnd
c      indicates the type of boundary conditions at z = c
c      and z = d.
c
c      = 0  if the solution is periodic in z, i.e.
c           u(i,j) = u(i,n+j).
c
c      = 1  if the solution is specified at z = c and z = d.
c
c      = 2  if the solution is specified at z = c and the derivative
c           of the solution with respect to z is specified at
c           z = d.
c
c      = 3  if the derivative of the solution with respect to z is
c           specified at z = c and z = d.
c
c      = 4  if the derivative of the solution with respect to z is
c           specified at z = c and the solution is specified at
c           z = d.
c
c    bdc
c      a one dimensional array of length m that specifies the boundary
c      values of the solution at z = c.   when nbdcnd = 1 or 2,
c
c               bdc(i) = u(r(i),c) ,              i=1,2,...,m.
c
c      when nbdcnd = 3 or 4,
c
c               bdc(i) = (d/dz)u(r(i),c),         i=1,2,...,m.
c
c      when nbdcnd = 0, bdc is a dummy variable.
c
c    bdd
c      a one-dimensional array of length m that specifies the boundary
c      values of the solution at z = d.  when nbdcnd = 1 or 4,
c
c               bdd(i) = u(r(i),d) ,              i=1,2,...,m.
c
c      when nbdcnd = 2 or 3,
c
c               bdd(i) = (d/dz)u(r(i),d) ,        i=1,2,...,m.
c
c      when nbdcnd = 0, bdd is a dummy variable.
c
c    elmbda
c      the constant lambda in the modified helmholtz equation.  if
c      lambda is greater than 0, a solution may not exist.  however,
c      hstcyl will attempt to find a solution.  lambda must be zero
c      when mbdcnd = 5 or 6.
c
c    f
c      a two-dimensional array that specifies the values of the right
c      side of the modified helmholtz equation.  for i=1,2,...,m
c      and j=1,2,...,n
c
c               f(i,j) = f(r(i),z(j)) .
c
c      f must be dimensioned at least m x n.
c
c    idimf
c      the row (or first) dimension of the array f as it appears in the
c      program calling hstcyl.  this parameter is used to specify the
c      variable dimension of f.  idimf must be at least m.
c
c    w
c      a one-dimensional array that must be provided by the user for
c      work space.  w may require up to 13m + 4n + m*int(log2(n))
c      locations.  the actual number of locations used is computed by
c      hstcyl and is returned in the location w(1).
c
c
c             * * * * * *   on output   * * * * * *
c
c    f
c      contains the solution u(i,j) of the finite difference
c      approximation for the grid point (r(i),z(j)) for
c      i=1,2,...,m, j=1,2,...,n.
c
c    pertrb
c      if a combination of periodic, derivative, or unspecified
c      boundary conditions is specified for a poisson equation
c      (lambda = 0), a solution may not exist.  pertrb is a con-
c      stant, calculated and subtracted from f, which ensures
c      that a solution exists.  hstcyl then computes this
c      solution, which is a least squares solution to the
c      original approximation.  this solution plus any constant is also
c      a solution; hence, the solution is not unique.  the value of
c      pertrb should be small compared to the right side f.
c      otherwise, a solution is obtained to an essentially different
c      problem.  this comparison should always be made to insure that
c      a meaningful solution has been obtained.
c
c    ierror
c      an error flag that indicates invalid input parameters.
c      except to numbers 0 and 11, a solution is not attempted.
c
c      =  0  no error
c
c      =  1  a .lt. 0
c
c      =  2  a .ge. b
c
c      =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6
c
c      =  4  c .ge. d
c
c      =  5  n .le. 2
c
c      =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c      =  7  a = 0 and mbdcnd = 1,2,3, or 4
c
c      =  8  a .gt. 0 and mbdcnd .ge. 5
c
c      =  9  m .le. 2
c
c      = 10  idimf .lt. m
c
c      = 11  lambda .gt. 0
c
c      = 12  a=0, mbdcnd .ge. 5, elmbda .ne. 0
c
c      since this is the only means of indicating a possibly
c      incorrect call to hstcyl, the user should test ierror after
c      the call.
c
c    w
c      w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c     arguments      w(see argument list)
c
c     latest         june 1, 1977
c     revision
c
c     subprograms    hstcyl,poistg,postg2,genbun,poisd2,poisn2,poisp2,
c     required       cosgen,merge,trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in march, 1977
c
c     algorithm      this subroutine defines the finite-difference
c                    equations, incorporates boundary data, adjusts the
c                    right side when the system is singular and calls
c                    either poistg or genbun which solves the linear
c                    system of equations.
c
c     space          8228(decimal) = 20044(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstcyl is roughly proportional
c                    to m*n*log2(n).  some typical values are listed in
c                    the table below.
c                       the solution process employed results in a loss
c                    of no more than four significant digits for n and m
c                    as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    poistg which is the routine that
c                    actually solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32       1-6       1-4         56
c                        64       1-6       1-4        230
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      schumann, u. and r. sweet,"a direct method for
c                    the solution of poisson"s equation with neumann
c                    boundary conditions on a staggered grid of
c                    arbitrary size," j. comp. phys. 20(1976),
c                    pp. 171-182.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dimension       f(idimf,1) ,bda(1)     ,bdb(1)     ,bdc(1)     ,
     &                bdd(1)     ,w(1)
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 2) ierror = 5
      if (nbdcnd.lt.0 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. mbdcnd.ne.5 .and. mbdcnd.ne.6) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (idimf .lt. m) ierror = 10
      if (m .le. 2) ierror = 9
      if (a.eq.0. .and. mbdcnd.ge.5 .and. elmbda.ne.0.) ierror = 12
      if (ierror .ne. 0) return
      deltar = (b-a)/float(m)
      dlrsq = deltar**2
      deltht = (d-c)/float(n)
      dlthsq = deltht**2
      np = nbdcnd+1
c
c     define a,b,c coefficients in w-array.
c
      iwb = m
      iwc = iwb+m
      iwr = iwc+m
      do 101 i=1,m
         j = iwr+i
         w(j) = a+(float(i)-0.5)*deltar
         w(i) = (a+float(i-1)*deltar)/(dlrsq*w(j))
         k = iwc+i
         w(k) = (a+float(i)*deltar)/(dlrsq*w(j))
         k = iwb+i
         w(k) = elmbda/w(j)**2-2./dlrsq
  101 continue
c
c     enter boundary data for r-boundaries.
c
      go to (102,102,104,104,106,106),mbdcnd
  102 a1 = 2.*w(1)
      w(iwb+1) = w(iwb+1)-w(1)
      do 103 j=1,n
         f(1,j) = f(1,j)-a1*bda(j)
  103 continue
      go to 106
  104 a1 = deltar*w(1)
      w(iwb+1) = w(iwb+1)+w(1)
      do 105 j=1,n
         f(1,j) = f(1,j)+a1*bda(j)
  105 continue
  106 continue
      go to (107,109,109,107,107,109),mbdcnd
  107 w(iwc) = w(iwc)-w(iwr)
      a1 = 2.*w(iwr)
      do 108 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  108 continue
      go to 111
  109 w(iwc) = w(iwc)+w(iwr)
      a1 = deltar*w(iwr)
      do 110 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  110 continue
c
c     enter boundary data for theta-boundaries.
c
  111 a1 = 2./dlthsq
      go to (121,112,112,114,114),np
  112 do 113 i=1,m
         f(i,1) = f(i,1)-a1*bdc(i)
  113 continue
      go to 116
  114 a1 = 1./deltht
      do 115 i=1,m
         f(i,1) = f(i,1)+a1*bdc(i)
  115 continue
  116 a1 = 2./dlthsq
      go to (121,117,119,119,117),np
  117 do 118 i=1,m
         f(i,n) = f(i,n)-a1*bdd(i)
  118 continue
      go to 121
  119 a1 = 1./deltht
      do 120 i=1,m
         f(i,n) = f(i,n)-a1*bdd(i)
  120 continue
  121 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 130,123,122
  122 ierror = 11
      go to 130
  123 go to (130,130,124,130,130,124),mbdcnd
  124 go to (125,130,130,125,130),np
  125 continue
      do 127 i=1,m
         a1 = 0.
         do 126 j=1,n
            a1 = a1+f(i,j)
  126    continue
         j = iwr+i
         pertrb = pertrb+a1*w(j)
  127 continue
      pertrb = pertrb/(float(m*n)*0.5*(a+b))
      do 129 i=1,m
         do 128 j=1,n
            f(i,j) = f(i,j)-pertrb
  128    continue
  129 continue
  130 continue
c
c     multiply i-th equation through by  deltht**2
c
      do 132 i=1,m
         w(i) = w(i)*dlthsq
         j = iwc+i
         w(j) = w(j)*dlthsq
         j = iwb+i
         w(j) = w(j)*dlthsq
         do 131 j=1,n
            f(i,j) = f(i,j)*dlthsq
  131    continue
  132 continue
      lp = nbdcnd
      w(1) = 0.
      w(iwr) = 0.
c
c     call genbun to solve the system of equations.
c
      if (nbdcnd .eq. 0) go to 133
      call poistg (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
      go to 134
  133 call genbun (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
  134 continue
      w(1) = w(iwr+1)+3.*float(m)
      return
      end
      subroutine hstplr (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HSTPLR 2D staggered grid polar coordinates five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c      hstplr solves the standard five-point finite difference
c      approximation on a staggered grid to the helmholtz equation in
c      polar coordinates
c
c      (1/r)(d/dr)(r(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta)
c
c                      + lambda*u = f(r,theta)
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c    a,b
c      the range of r, i.e. a .le. r .le. b.  a must be less than b and
c      a must be non-negative.
c
c    m
c      the number of grid points in the interval (a,b).  the grid points
c      in the r-direction are given by r(i) = a + (i-0.5)dr for
c      i=1,2,...,m where dr =(b-a)/m.  m must be greater than 2.
c
c    mbdcnd
c      indicates the type of boundary conditions at r = a and r = b.
c
c      = 1  if the solution is specified at r = a and r = b.
c
c      = 2  if the solution is specified at r = a and the derivative
c           of the solution with respect to r is specified at r = b.
c           (see note 1 below)
c
c      = 3  if the derivative of the solution with respect to r is
c           specified at r = a (see note 2 below) and r = b.
c
c      = 4  if the derivative of the solution with respect to r is
c           specified at r = a (see note 2 below) and the solution is
c           specified at r = b.
c
c      = 5  if the solution is unspecified at r = a = 0 and the solution
c           is specified at r = b.
c
c      = 6  if the solution is unspecified at r = a = 0 and the
c           derivative of the solution with respect to r is specified at
c           r = b.
c
c      note 1:  if a = 0, mbdcnd = 2, and nbdcnd = 0 or 3, the system of
c               equations to be solved is singular.  the unique solution
c               is determined by extrapolation to the specification of
c               u(0,theta(1)).  but in this case the right side of the
c               system will be perturbed by the constant pertrb.
c
c      note 2:  if a = 0, do not use mbdcnd = 3 or 4, but instead use
c               mbdcnd = 1,2,5, or 6.
c
c    bda
c      a one-dimensional array of length n that specifies the boundary
c      values (if any) of the solution at r = a.  when mbdcnd = 1 or 2,
c
c               bda(j) = u(a,theta(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 3 or 4,
c
c               bda(j) = (d/dr)u(a,theta(j)) ,    j=1,2,...,n.
c
c      when mbdcnd = 5 or 6, bda is a dummy variable.
c
c    bdb
c      a one-dimensional array of length n that specifies the boundary
c      values of the solution at r = b.  when mbdcnd = 1,4, or 5,
c
c               bdb(j) = u(b,theta(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 2,3, or 6,
c
c               bdb(j) = (d/dr)u(b,theta(j)) ,    j=1,2,...,n.
c
c    c,d
c      the range of theta, i.e. c .le. theta .le. d.  c must be less
c      than d.
c
c    n
c      the number of unknowns in the interval (c,d).  the unknowns in
c      the theta-direction are given by theta(j) = c + (j-0.5)dt,
c      j=1,2,...,n, where dt = (d-c)/n.  n must be greater than 2.
c
c    nbdcnd
c      indicates the type of boundary conditions at theta = c
c      and theta = d.
c
c      = 0  if the solution is periodic in theta, i.e.
c           u(i,j) = u(i,n+j).
c
c      = 1  if the solution is specified at theta = c and theta = d
c           (see note below).
c
c      = 2  if the solution is specified at theta = c and the derivative
c           of the solution with respect to theta is specified at
c           theta = d (see note below).
c
c      = 3  if the derivative of the solution with respect to theta is
c           specified at theta = c and theta = d.
c
c      = 4  if the derivative of the solution with respect to theta is
c           specified at theta = c and the solution is specified at
c           theta = d (see note below).
c
c      note:  when nbdcnd = 1, 2, or 4, do not use mbdcnd = 5 or 6 (the
c      former indicates that the solution is specified at r =  0; the
c      latter indicates the solution is unspecified at r = 0).  use
c      instead mbdcnd = 1 or 2.
c
c    bdc
c      a one dimensional array of length m that specifies the boundary
c      values of the solution at theta = c.   when nbdcnd = 1 or 2,
c
c               bdc(i) = u(r(i),c) ,              i=1,2,...,m.
c
c      when nbdcnd = 3 or 4,
c
c               bdc(i) = (d/dtheta)u(r(i),c),     i=1,2,...,m.
c
c      when nbdcnd = 0, bdc is a dummy variable.
c
c    bdd
c      a one-dimensional array of length m that specifies the boundary
c      values of the solution at theta = d.  when nbdcnd = 1 or 4,
c
c               bdd(i) = u(r(i),d) ,              i=1,2,...,m.
c
c      when nbdcnd = 2 or 3,
c
c               bdd(i) = (d/dtheta)u(r(i),d) ,    i=1,2,...,m.
c
c      when nbdcnd = 0, bdd is a dummy variable.
c
c    elmbda
c      the constant lambda in the helmholtz equation.  if lambda is
c      greater than 0, a solution may not exist.  however, hstplr will
c      attempt to find a solution.
c
c    f
c      a two-dimensional array that specifies the values of the right
c      side of the helmholtz equation.  for i=1,2,...,m and j=1,2,...,n
c
c               f(i,j) = f(r(i),theta(j)) .
c
c      f must be dimensioned at least m x n.
c
c    idimf
c      the row (or first) dimension of the array f as it appears in the
c      program calling hstplr.  this parameter is used to specify the
c      variable dimension of f.  idimf must be at least m.
c
c    w
c      a one-dimensional array that must be provided by the user for
c      work space.  w may require up to 13m + 4n + m*int(log2(n))
c      locations.  the actual number of locations used is computed by
c      hstplr and is returned in the location w(1).
c
c
c             * * * * * *   on output   * * * * * *
c
c    f
c      contains the solution u(i,j) of the finite difference
c      approximation for the grid point (r(i),theta(j)) for
c      i=1,2,...,m, j=1,2,...,n.
c
c    pertrb
c      if a combination of periodic, derivative, or unspecified
c      boundary conditions is specified for a poisson equation
c      (lambda = 0), a solution may not exist.  pertrb is a con-
c      stant, calculated and subtracted from f, which ensures
c      that a solution exists.  hstplr then computes this
c      solution, which is a least squares solution to the
c      original approximation.  this solution plus any constant is also
c      a solution; hence, the solution is not unique.  the value of
c      pertrb should be small compared to the right side f.
c      otherwise, a solution is obtained to an essentially different
c      problem.  this comparison should always be made to insure that
c      a meaningful solution has been obtained.
c
c    ierror
c      an error flag that indicates invalid input parameters.
c      except to numbers 0 and 11, a solution is not attempted.
c
c      =  0  no error
c
c      =  1  a .lt. 0
c
c      =  2  a .ge. b
c
c      =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6
c
c      =  4  c .ge. d
c
c      =  5  n .le. 2
c
c      =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c      =  7  a = 0 and mbdcnd = 3 or 4
c
c      =  8  a .gt. 0 and mbdcnd .ge. 5
c
c      =  9  mbdcnd .ge. 5 and nbdcnd .ne. 0 or 3
c
c      = 10  idimf .lt. m
c
c      = 11  lambda .gt. 0
c
c      = 12  m .le. 2
c
c      since this is the only means of indicating a possibly
c      incorrect call to hstplr, the user should test ierror after
c      the call.
c
c    w
c      w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c     arguments      w(see argument list)
c
c     latest         june 1, 1977
c     revision
c
c     subprograms    hstplr,poistg,postg2,genbun,poisd2,poisn2,poisp2,
c     required       cosgen,merge,trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in february, 1977
c
c     algorithm      this subroutine defines the finite-difference
c                    equations, incorporates boundary data, adjusts the
c                    right side when the system is singular and calls
c                    either poistg or genbun which solves the linear
c                    system of equations.
c
c     space          8265(decimal) = 20111(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstplr is roughly proportional
c                    to m*n*log2(n).  some typical values are listed in
c                    the table below.
c                       the solution process employed results in a loss
c                    of no more than four significant digits for n and m
c                    as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    poistg which is the routine that
c                    actually solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32       1-6       1-4         56
c                        64       1-6       1-4        230
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      schumann, u. and r. sweet,"a direct method for
c                    the solution of poisson"s equation with neumann
c                    boundary conditions on a staggered grid of
c                    arbitrary size," j. comp. phys. 20(1976),
c                    pp. 171-182.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dimension       f(idimf,1)
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     &                w(1)
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 2) ierror = 5
      if (nbdcnd.lt.0 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4)) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (mbdcnd.ge.5 .and. nbdcnd.ne.0 .and. nbdcnd.ne.3) ierror = 9
      if (idimf .lt. m) ierror = 10
      if (m .le. 2) ierror = 12
      if (ierror .ne. 0) return
      deltar = (b-a)/float(m)
      dlrsq = deltar**2
      deltht = (d-c)/float(n)
      dlthsq = deltht**2
      np = nbdcnd+1
      isw = 1
      mb = mbdcnd
      if (a.eq.0. .and. mbdcnd.eq.2) mb = 6
c
c     define a,b,c coefficients in w-array.
c
      iwb = m
      iwc = iwb+m
      iwr = iwc+m
      do 101 i=1,m
         j = iwr+i
         w(j) = a+(float(i)-0.5)*deltar
         w(i) = (a+float(i-1)*deltar)/dlrsq
         k = iwc+i
         w(k) = (a+float(i)*deltar)/dlrsq
         k = iwb+i
         w(k) = (elmbda-2./dlrsq)*w(j)
  101 continue
      do 103 i=1,m
         j = iwr+i
         a1 = w(j)
         do 102 j=1,n
            f(i,j) = a1*f(i,j)
  102    continue
  103 continue
c
c     enter boundary data for r-boundaries.
c
      go to (104,104,106,106,108,108),mb
  104 a1 = 2.*w(1)
      w(iwb+1) = w(iwb+1)-w(1)
      do 105 j=1,n
         f(1,j) = f(1,j)-a1*bda(j)
  105 continue
      go to 108
  106 a1 = deltar*w(1)
      w(iwb+1) = w(iwb+1)+w(1)
      do 107 j=1,n
         f(1,j) = f(1,j)+a1*bda(j)
  107 continue
  108 go to (109,111,111,109,109,111),mb
  109 a1 = 2.*w(iwr)
      w(iwc) = w(iwc)-w(iwr)
      do 110 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  110 continue
      go to 113
  111 a1 = deltar*w(iwr)
      w(iwc) = w(iwc)+w(iwr)
      do 112 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  112 continue
c
c     enter boundary data for theta-boundaries.
c
  113 a1 = 2./dlthsq
      go to (123,114,114,116,116),np
  114 do 115 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)-a1*bdc(i)/w(j)
  115 continue
      go to 118
  116 a1 = 1./deltht
      do 117 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)+a1*bdc(i)/w(j)
  117 continue
  118 a1 = 2./dlthsq
      go to (123,119,121,121,119),np
  119 do 120 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  120 continue
      go to 123
  121 a1 = 1./deltht
      do 122 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  122 continue
  123 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 133,125,124
  124 ierror = 11
      go to 133
  125 go to (133,133,126,133,133,126),mb
  126 go to (127,133,133,127,133),np
  127 continue
      isw = 2
      do 129 j=1,n
         do 128 i=1,m
            pertrb = pertrb+f(i,j)
  128    continue
  129 continue
      pertrb = pertrb/(float(m*n)*0.5*(a+b))
      do 131 i=1,m
         j = iwr+i
         a1 = pertrb*w(j)
         do 130 j=1,n
            f(i,j) = f(i,j)-a1
  130    continue
  131 continue
      a2 = 0.
      do 132 j=1,n
         a2 = a2+f(1,j)
  132 continue
      a2 = a2/w(iwr+1)
  133 continue
c
c     multiply i-th equation through by  r(i)*deltht**2
c
      do 135 i=1,m
         j = iwr+i
         a1 = dlthsq*w(j)
         w(i) = a1*w(i)
         j = iwc+i
         w(j) = a1*w(j)
         j = iwb+i
         w(j) = a1*w(j)
         do 134 j=1,n
            f(i,j) = a1*f(i,j)
  134    continue
  135 continue
      lp = nbdcnd
      w(1) = 0.
      w(iwr) = 0.
c
c     call poistg or genbun to solve the system of equations.
c
      if (lp .eq. 0) go to 136
      call poistg (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
      go to 137
  136 call genbun (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
  137 continue
      w(1) = w(iwr+1)+3.*float(m)
      if (a.ne.0. .or. mbdcnd.ne.2 .or. isw.ne.2) go to 141
      a1 = 0.
      do 138 j=1,n
         a1 = a1+f(1,j)
  138 continue
      a1 = (a1-dlrsq*a2/16.)/float(n)
      if (nbdcnd .eq. 3) a1 = a1+(bdd(1)-bdc(1))/(d-c)
      a1 = bda(1)-a1
      do 140 i=1,m
         do 139 j=1,n
            f(i,j) = f(i,j)+a1
  139    continue
  140 continue
  141 continue
      return
      end
      subroutine hstssp (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HSTSSP 2D staggered unit sphere coordinates five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c     hstssp solves the standard five-point finite difference
c     approximation on a staggered grid to the helmholtz equation in
c     spherical coordinates and on the surface of the unit sphere
c     (radius of 1)
c
c             (1/sin(theta))(d/dtheta)(sin(theta)(du/dtheta)) +
c
c       (1/sin(theta)**2)(d/dphi)(du/dphi) + lambda*u = f(theta,phi)
c
c     where theta is colatitude and phi is longitude.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c            * * * * * *   on input    * * * * * *
c
c   a,b
c     the range of theta (colatitude), i.e. a .le. theta .le. b.  a
c     must be less than b and a must be non-negative.  a and b are in
c     radians.  a = 0 corresponds to the north pole and b = pi
c     corresponds to the south pole.
c
c
c                  * * *  important  * * *
c
c     if b is equal to pi, then b must be computed using the statement
c
c     b = pimach()
c
c     this insures that b in the user"s program is equal to pi in this
c     program which permits several tests of the input parameters that
c     otherwise would not be possible.
c
c                  * * * * * * * * * * * *
c
c
c
c   m
c     the number of grid points in the interval (a,b).  the grid points
c     in the theta-direction are given by theta(i) = a + (i-0.5)dtheta
c     for i=1,2,...,m where dtheta =(b-a)/m.  m must be greater than 2.
c
c   mbdcnd
c     indicates the type of boundary conditions at theta = a and
c     theta = b.
c
c     = 1  if the solution is specified at theta = a and theta = b.
c          (see note 3 below)
c
c     = 2  if the solution is specified at theta = a and the derivative
c          of the solution with respect to theta is specified at
c          theta = b (see notes 2 and 3 below).
c
c     = 3  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1, 2 below) and theta = b.
c
c     = 4  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1 and 2 below) and the
c          solution is specified at theta = b.
c
c     = 5  if the solution is unspecified at theta = a = 0 and the
c          solution is specified at theta = b.  (see note 3 below)
c
c     = 6  if the solution is unspecified at theta = a = 0 and the
c          derivative of the solution with respect to theta is
c          specified at theta = b (see note 2 below).
c
c     = 7  if the solution is specified at theta = a and the
c          solution is unspecified at theta = b = pi. (see note 3 below)
c
c     = 8  if the derivative of the solution with respect to
c          theta is specified at theta = a (see note 1 below)
c          and the solution is unspecified at theta = b = pi.
c
c     = 9  if the solution is unspecified at theta = a = 0 and
c          theta = b = pi.
c
c     notes:  1.  if a = 0, do not use mbdcnd = 3, 4, or 8,
c                 but instead use mbdcnd = 5, 6, or 9.
c
c             2.  if b = pi, do not use mbdcnd = 2, 3, or 6,
c                 but instead use mbdcnd = 7, 8, or 9.
c
c             3.  when the solution is specified at theta = 0 and/or
c                 theta = pi and the other boundary conditions are
c                 combinations of unspecified, normal derivative, or
c                 periodicity a singular system results.  the unique
c                 solution is determined by extrapolation to the
c                 specification of the solution at either theta = 0 or
c                 theta = pi.  but in these cases the right side of the
c                 system will be perturbed by the constant pertrb.
c
c   bda
c     a one-dimensional array of length n that specifies the boundary
c     values (if any) of the solution at theta = a.  when
c     mbdcnd = 1, 2, or 7,
c
c              bda(j) = u(a,phi(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 3, 4, or 8,
c
c              bda(j) = (d/dtheta)u(a,phi(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bda is a dummy variable.
c
c   bdb
c     a one-dimensional array of length n that specifies the boundary
c     values of the solution at theta = b.  when mbdcnd = 1,4, or 5,
c
c              bdb(j) = u(b,phi(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 2,3, or 6,
c
c              bdb(j) = (d/dtheta)u(b,phi(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bdb is a dummy variable.
c
c   c,d
c     the range of phi (longitude), i.e. c .le. phi .le. d.
c     c must be less than d.  if d-c = 2*pi, periodic boundary
c     conditions are usually prescribed.
c
c   n
c     the number of unknowns in the interval (c,d).  the unknowns in
c     the phi-direction are given by phi(j) = c + (j-0.5)dphi,
c     j=1,2,...,n, where dphi = (d-c)/n.  n must be greater than 2.
c
c   nbdcnd
c     indicates the type of boundary conditions at phi = c
c     and phi = d.
c
c     = 0  if the solution is periodic in phi, i.e.
c          u(i,j) = u(i,n+j).
c
c     = 1  if the solution is specified at phi = c and phi = d
c          (see note below).
c
c     = 2  if the solution is specified at phi = c and the derivative
c          of the solution with respect to phi is specified at
c          phi = d (see note below).
c
c     = 3  if the derivative of the solution with respect to phi is
c          specified at phi = c and phi = d.
c
c     = 4  if the derivative of the solution with respect to phi is
c          specified at phi = c and the solution is specified at
c          phi = d (see note below).
c
c     note:  when nbdcnd = 1, 2, or 4, do not use mbdcnd = 5, 6, 7, 8,
c     or 9 (the former indicates that the solution is specified at
c     a pole; the latter indicates the solution is unspecified).  use
c     instead mbdcnd = 1 or 2.
c
c   bdc
c     a one dimensional array of length m that specifies the boundary
c     values of the solution at phi = c.   when nbdcnd = 1 or 2,
c
c              bdc(i) = u(theta(i),c) ,              i=1,2,...,m.
c
c     when nbdcnd = 3 or 4,
c
c              bdc(i) = (d/dphi)u(theta(i),c),       i=1,2,...,m.
c
c     when nbdcnd = 0, bdc is a dummy variable.
c
c   bdd
c     a one-dimensional array of length m that specifies the boundary
c     values of the solution at phi = d.  when nbdcnd = 1 or 4,
c
c              bdd(i) = u(theta(i),d) ,              i=1,2,...,m.
c
c     when nbdcnd = 2 or 3,
c
c              bdd(i) = (d/dphi)u(theta(i),d) ,      i=1,2,...,m.
c
c     when nbdcnd = 0, bdd is a dummy variable.
c
c   elmbda
c     the constant lambda in the helmholtz equation.  if lambda is
c     greater than 0, a solution may not exist.  however, hstssp will
c     attempt to find a solution.
c
c   f
c     a two-dimensional array that specifies the values of the right
c     side of the helmholtz equation.  for i=1,2,...,m and j=1,2,...,n
c
c              f(i,j) = f(theta(i),phi(j)) .
c
c     f must be dimensioned at least m x n.
c
c   idimf
c     the row (or first) dimension of the array f as it appears in the
c     program calling hstssp.  this parameter is used to specify the
c     variable dimension of f.  idimf must be at least m.
c
c   w
c     a one-dimensional array that must be provided by the user for
c     work space.  w may require up to 13m + 4n + m*int(log2(n))
c     locations.  the actual number of locations used is computed by
c     hstssp and is returned in the location w(1).
c
c
c            * * * * * *   on output   * * * * * *
c
c   f
c     contains the solution u(i,j) of the finite difference
c     approximation for the grid point (theta(i),phi(j)) for
c     i=1,2,...,m, j=1,2,...,n.
c
c   pertrb
c     if a combination of periodic, derivative, or unspecified
c     boundary conditions is specified for a poisson equation
c     (lambda = 0), a solution may not exist.  pertrb is a con-
c     stant, calculated and subtracted from f, which ensures
c     that a solution exists.  hstssp then computes this
c     solution, which is a least squares solution to the
c     original approximation.  this solution plus any constant is also
c     a solution; hence, the solution is not unique.  the value of
c     pertrb should be small compared to the right side f.
c     otherwise, a solution is obtained to an essentially different
c     problem.  this comparison should always be made to insure that
c     a meaningful solution has been obtained.
c
c   ierror
c     an error flag that indicates invalid input parameters.
c     except to numbers 0 and 14, a solution is not attempted.
c
c     =  0  no error
c
c     =  1  a .lt. 0 or b .gt. pi
c
c     =  2  a .ge. b
c
c     =  3  mbdcnd .lt. 1 or mbdcnd .gt. 9
c
c     =  4  c .ge. d
c
c     =  5  n .le. 2
c
c     =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c     =  7  a .gt. 0 and mbdcnd = 5, 6, or 9
c
c     =  8  a = 0 and mbdcnd = 3, 4, or 8
c
c     =  9  b .lt. pi and mbdcnd .ge. 7
c
c     = 10  b = pi and mbdcnd = 2,3, or 6
c
c     = 11  mbdcnd .ge. 5 and ndbcnd = 1, 2, or 4
c
c     = 12  idimf .lt. m
c
c     = 13  m .le. 2
c
c     = 14  lambda .gt. 0
c
c     since this is the only means of indicating a possibly
c     incorrect call to hstssp, the user should test ierror after
c     the call.
c
c   w
c     w(1) contains the required length of w.
c
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c    dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c    arguments      w(see argument list)
c
c    latest         june 1, 1977
c    revision
c
c    subprograms    hstssp,poistg,postg2,genbun,poisd2,poisn2,poisp2,
c    required       cosgen,merge,trix,tri3,pimach
c
c    special        none
c    conditions
c
c    common         none
c    blocks
c
c    i/o            none
c
c    precision      single
c
c    specialist     roland sweet
c
c    language       fortran
c
c    history        written by roland sweet at ncar in april, 1977
c
c    algorithm      this subroutine defines the finite-difference
c                   equations, incorporates boundary data, adjusts the
c                   right side when the system is singular and calls
c                   either poistg or genbun which solves the linear
c                   system of equations.
c
c    space          8427(decimal) = 20353(octal) locations on the
c    required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstssp is roughly proportional
c                    to m*n*log2(n).  some typical values are listed in
c                    the table below.
c                       the solution process employed results in a loss
c                    of no more than four significant digits for n and m
c                    as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    poistg which is the routine that
c                    actually solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32       1-9       1-4         56
c                        64       1-9       1-4        230
c
c    portability    american national standards institute fortran.
c                   all machine dependent constants are located in the
c                   function pimach.
c
c    required       cos
c    resident
c    routines
c
c    reference      schumann, u. and r. sweet,"a direct method for
c                   the solution of poisson"s equation with neumann
c                   boundary conditions on a staggered grid of
c                   arbitrary size," j. comp. phys. 20(1976),
c                   pp. 171-182.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dimension       f(idimf,1) ,bda(1)     ,bdb(1)     ,bdc(1)     ,
     &                bdd(1)     ,w(1)
      ierror = 0
      pi = pimach()
      if (a.lt.0. .or. b.gt.pi) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.gt.9) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 2) ierror = 5
      if (nbdcnd.lt.0 .or. nbdcnd.ge.5) ierror = 6
      if (a.gt.0. .and. (mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9))
     &    ierror = 7
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4 .or. mbdcnd.eq.8))
     &    ierror = 8
      if (b.lt.pi .and. mbdcnd.ge.7) ierror = 9
      if (b.eq.pi .and. (mbdcnd.eq.2 .or. mbdcnd.eq.3 .or. mbdcnd.eq.6))
     &    ierror = 10
      if (mbdcnd.ge.5 .and.
     &    (nbdcnd.eq.1 .or. nbdcnd.eq.2 .or. nbdcnd.eq.4)) ierror = 11
      if (idimf .lt. m) ierror = 12
      if (m .le. 2) ierror = 13
      if (ierror .ne. 0) return
      deltar = (b-a)/float(m)
      dlrsq = deltar**2
      deltht = (d-c)/float(n)
      dlthsq = deltht**2
      np = nbdcnd+1
      isw = 1
      jsw = 1
      mb = mbdcnd
      if (elmbda .ne. 0.) go to 105
      go to (101,102,105,103,101,105,101,105,105),mbdcnd
  101 if (a.ne.0. .or. b.ne.pi) go to 105
      mb = 9
      go to 104
  102 if (a .ne. 0.) go to 105
      mb = 6
      go to 104
  103 if (b .ne. pi) go to 105
      mb = 8
  104 jsw = 2
  105 continue
c
c     define a,b,c coefficients in w-array.
c
      iwb = m
      iwc = iwb+m
      iwr = iwc+m
      iws = iwr+m
      do 106 i=1,m
         j = iwr+i
         w(j) = sin(a+(float(i)-0.5)*deltar)
         w(i) = sin((a+float(i-1)*deltar))/dlrsq
  106 continue
      mm1 = m-1
      do 107 i=1,mm1
         k = iwc+i
         w(k) = w(i+1)
         j = iwr+i
         k = iwb+i
         w(k) = elmbda*w(j)-(w(i)+w(i+1))
  107 continue
      w(iwr) = sin(b)/dlrsq
      w(iwc) = elmbda*w(iws)-(w(m)+w(iwr))
      do 109 i=1,m
         j = iwr+i
         a1 = w(j)
         do 108 j=1,n
            f(i,j) = a1*f(i,j)
  108    continue
  109 continue
c
c     enter boundary data for theta-boundaries.
c
      go to (110,110,112,112,114,114,110,112,114),mb
  110 a1 = 2.*w(1)
      w(iwb+1) = w(iwb+1)-w(1)
      do 111 j=1,n
         f(1,j) = f(1,j)-a1*bda(j)
  111 continue
      go to 114
  112 a1 = deltar*w(1)
      w(iwb+1) = w(iwb+1)+w(1)
      do 113 j=1,n
         f(1,j) = f(1,j)+a1*bda(j)
  113 continue
  114 go to (115,117,117,115,115,117,119,119,119),mb
  115 a1 = 2.*w(iwr)
      w(iwc) = w(iwc)-w(iwr)
      do 116 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  116 continue
      go to 119
  117 a1 = deltar*w(iwr)
      w(iwc) = w(iwc)+w(iwr)
      do 118 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  118 continue
c
c     enter boundary data for phi-boundaries.
c
  119 a1 = 2./dlthsq
      go to (129,120,120,122,122),np
  120 do 121 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)-a1*bdc(i)/w(j)
  121 continue
      go to 124
  122 a1 = 1./deltht
      do 123 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)+a1*bdc(i)/w(j)
  123 continue
  124 a1 = 2./dlthsq
      go to (129,125,127,127,125),np
  125 do 126 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  126 continue
      go to 129
  127 a1 = 1./deltht
      do 128 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  128 continue
  129 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 139,131,130
  130 ierror = 14
      go to 139
  131 go to (139,139,132,139,139,132,139,132,132),mb
  132 go to (133,139,139,133,139),np
  133 continue
      isw = 2
      do 135 j=1,n
         do 134 i=1,m
            pertrb = pertrb+f(i,j)
  134    continue
  135 continue
      a1 = float(n)*(cos(a)-cos(b))/(2.*sin(0.5*deltar))
      pertrb = pertrb/a1
      do 137 i=1,m
         j = iwr+i
         a1 = pertrb*w(j)
         do 136 j=1,n
            f(i,j) = f(i,j)-a1
  136    continue
  137 continue
      a2 = 0.
      a3 = 0.
      do 138 j=1,n
         a2 = a2+f(1,j)
         a3 = a3+f(m,j)
  138 continue
      a2 = a2/w(iwr+1)
      a3 = a3/w(iws)
  139 continue
c
c     multiply i-th equation through by  r(i)*deltht**2
c
      do 141 i=1,m
         j = iwr+i
         a1 = dlthsq*w(j)
         w(i) = a1*w(i)
         j = iwc+i
         w(j) = a1*w(j)
         j = iwb+i
         w(j) = a1*w(j)
         do 140 j=1,n
            f(i,j) = a1*f(i,j)
  140    continue
  141 continue
      lp = nbdcnd
      w(1) = 0.
      w(iwr) = 0.
c
c     call poistg or genbun to solve the system of equations.
c
      if (nbdcnd .eq. 0) go to 142
      call poistg (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
      go to 143
  142 call genbun (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
  143 continue
      w(1) = w(iwr+1)+3.*float(m)
      if (isw.ne.2 .or. jsw.ne.2) go to 150
      if (mb .ne. 8) go to 145
      a1 = 0.
      do 144 j=1,n
         a1 = a1+f(m,j)
  144 continue
      a1 = (a1-dlrsq*a3/16.)/float(n)
      if (nbdcnd .eq. 3) a1 = a1+(bdd(m)-bdc(m))/(d-c)
      a1 = bdb(1)-a1
      go to 147
  145 a1 = 0.
      do 146 j=1,n
         a1 = a1+f(1,j)
  146 continue
      a1 = (a1-dlrsq*a2/16.)/float(n)
      if (nbdcnd .eq. 3) a1 = a1+(bdd(1)-bdc(1))/(d-c)
      a1 = bda(1)-a1
  147 do 149 i=1,m
         do 148 j=1,n
            f(i,j) = f(i,j)+a1
  148    continue
  149 continue
  150 continue
      return
      end
      subroutine hw3crt (xs,xf,l,lbdcnd,bdxs,bdxf,ys,yf,m,mbdcnd,bdys,
     &                   bdyf,zs,zf,n,nbdcnd,bdzs,bdzf,elmbda,ldimf,
     &                   mdimf,f,pertrb,ierror,w)

c*********************************************************************72
c
cc HW3CRT 3D standard grid Cartesian coordinates seven point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c     hw3crt solves the standard seven-point finite
c     difference approximation to the helmholtz equation in cartesian
c     coordinates:
c
c         (d/dx)(du/dx) + (d/dy)(du/dy) + (d/dz)(du/dz)
c
c                    + lambda*u = f(x,y,z) .
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c     xs,xf
c        the range of x, i.e. xs .le. x .le. xf .
c        xs must be less than xf.
c
c     l
c        the number of panels into which the interval (xs,xf) is
c        subdivided.  hence, there will be l+1 grid points in the
c        x-direction given by x(i) = xs+(i-1)dx for i=1,2,...,l+1,
c        where dx = (xf-xs)/l is the panel width.  l must be at
c        least 5 .
c
c     lbdcnd
c        indicates the type of boundary conditions at x = xs and x = xf.
c
c        = 0  if the solution is periodic in x, i.e.
c             u(l+i,j,k) = u(i,j,k).
c        = 1  if the solution is specified at x = xs and x = xf.
c        = 2  if the solution is specified at x = xs and the derivative
c             of the solution with respect to x is specified at x = xf.
c        = 3  if the derivative of the solution with respect to x is
c             specified at x = xs and x = xf.
c        = 4  if the derivative of the solution with respect to x is
c             specified at x = xs and the solution is specified at x=xf.
c
c     bdxs
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to x at x = xs.
c        when lbdcnd = 3 or 4,
c
c             bdxs(j,k) = (d/dx)u(xs,y(j),z(k)), j=1,2,...,m+1,
c                                                k=1,2,...,n+1.
c
c        when lbdcnd has any other value, bdxs is a dummy variable.
c        bdxs must be dimensioned at least (m+1)*(n+1).
c
c     bdxf
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to x at x = xf.
c        when lbdcnd = 2 or 3,
c
c             bdxf(j,k) = (d/dx)u(xf,y(j),z(k)), j=1,2,...,m+1,
c                                                k=1,2,...,n+1.
c
c        when lbdcnd has any other value, bdxf is a dummy variable.
c        bdxf must be dimensioned at least (m+1)*(n+1).
c
c     ys,yf
c        the range of y, i.e. ys .le. y .le. yf.
c        ys must be less than yf.
c
c     m
c        the number of panels into which the interval (ys,yf) is
c        subdivided.  hence, there will be m+1 grid points in the
c        y-direction given by y(j) = ys+(j-1)dy for j=1,2,...,m+1,
c        where dy = (yf-ys)/m is the panel width.  m must be at
c        least 5 .
c
c     mbdcnd
c        indicates the type of boundary conditions at y = ys and y = yf.
c
c        = 0  if the solution is periodic in y, i.e.
c             u(i,m+j,k) = u(i,j,k).
c        = 1  if the solution is specified at y = ys and y = yf.
c        = 2  if the solution is specified at y = ys and the derivative
c             of the solution with respect to y is specified at y = yf.
c        = 3  if the derivative of the solution with respect to y is
c             specified at y = ys and y = yf.
c        = 4  if the derivative of the solution with respect to y is
c             specified at y = ys and the solution is specified at y=yf.
c
c     bdys
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to y at y = ys.
c        when mbdcnd = 3 or 4,
c
c             bdys(i,k) = (d/dy)u(x(i),ys,z(k)), i=1,2,...,l+1,
c                                                k=1,2,...,n+1.
c
c        when mbdcnd has any other value, bdys is a dummy variable.
c        bdys must be dimensioned at least (l+1)*(n+1).
c
c     bdyf
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to y at y = yf.
c        when mbdcnd = 2 or 3,
c
c             bdyf(i,k) = (d/dy)u(x(i),yf,z(k)), i=1,2,...,l+1,
c                                                k=1,2,...,n+1.
c
c        when mbdcnd has any other value, bdyf is a dummy variable.
c        bdyf must be dimensioned at least (l+1)*(n+1).
c
c     zs,zf
c        the range of z, i.e. zs .le. z .le. zf.
c        zs must be less than zf.
c
c     n
c        the number of panels into which the interval (zs,zf) is
c        subdivided.  hence, there will be n+1 grid points in the
c        z-direction given by z(k) = zs+(k-1)dz for k=1,2,...,n+1,
c        where dz = (zf-zs)/n is the panel width.  n must be at least 5.
c
c     nbdcnd
c        indicates the type of boundary conditions at z = zs and z = zf.
c
c        = 0  if the solution is periodic in z, i.e.
c             u(i,j,n+k) = u(i,j,k).
c        = 1  if the solution is specified at z = zs and z = zf.
c        = 2  if the solution is specified at z = zs and the derivative
c             of the solution with respect to z is specified at z = zf.
c        = 3  if the derivative of the solution with respect to z is
c             specified at z = zs and z = zf.
c        = 4  if the derivative of the solution with respect to z is
c             specified at z = zs and the solution is specified at z=zf.
c
c     bdzs
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to z at z = zs.
c        when nbdcnd = 3 or 4,
c
c             bdzs(i,j) = (d/dz)u(x(i),y(j),zs), i=1,2,...,l+1,
c                                                j=1,2,...,m+1.
c
c        when nbdcnd has any other value, bdzs is a dummy variable.
c        bdzs must be dimensioned at least (l+1)*(m+1).
c
c     bdzf
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to z at z = zf.
c        when nbdcnd = 2 or 3,
c
c             bdzf(i,j) = (d/dz)u(x(i),y(j),zf), i=1,2,...,l+1,
c                                                j=1,2,...,m+1.
c
c        when nbdcnd has any other value, bdzf is a dummy variable.
c        bdzf must be dimensioned at least (l+1)*(m+1).
c
c     elmbda
c        the constant lambda in the helmholtz equation. if
c        lambda .gt. 0, a solution may not exist.  however, hw3crt will
c        attempt to find a solution.
c
c     f
c        a three-dimensional array that specifies the values of the
c        right side of the helmholtz equation and boundary values (if
c        any).  for i=2,3,...,l, j=2,3,...,m, and k=2,3,...,n
c
c                   f(i,j,k) = f(x(i),y(j),z(k)).
c
c        on the boundaries f is defined by
c
c        lbdcnd      f(1,j,k)         f(l+1,j,k)
c        ------   ---------------   ---------------
c
c          0      f(xs,y(j),z(k))   f(xs,y(j),z(k))
c          1      u(xs,y(j),z(k))   u(xf,y(j),z(k))
c          2      u(xs,y(j),z(k))   f(xf,y(j),z(k))   j=1,2,...,m+1
c          3      f(xs,y(j),z(k))   f(xf,y(j),z(k))   k=1,2,...,n+1
c          4      f(xs,y(j),z(k))   u(xf,y(j),z(k))
c
c        mbdcnd      f(i,1,k)         f(i,m+1,k)
c        ------   ---------------   ---------------
c
c          0      f(x(i),ys,z(k))   f(x(i),ys,z(k))
c          1      u(x(i),ys,z(k))   u(x(i),yf,z(k))
c          2      u(x(i),ys,z(k))   f(x(i),yf,z(k))   i=1,2,...,l+1
c          3      f(x(i),ys,z(k))   f(x(i),yf,z(k))   k=1,2,...,n+1
c          4      f(x(i),ys,z(k))   u(x(i),yf,z(k))
c
c        nbdcnd      f(i,j,1)         f(i,j,n+1)
c        ------   ---------------   ---------------
c
c          0      f(x(i),y(j),zs)   f(x(i),y(j),zs)
c          1      u(x(i),y(j),zs)   u(x(i),y(j),zf)
c          2      u(x(i),y(j),zs)   f(x(i),y(j),zf)   i=1,2,...,l+1
c          3      f(x(i),y(j),zs)   f(x(i),y(j),zf)   j=1,2,...,m+1
c          4      f(x(i),y(j),zs)   u(x(i),y(j),zf)
c
c        f must be dimensioned at least (l+1)*(m+1)*(n+1).
c
c        note:
c
c        if the table calls for both the solution u and the right side f
c        on a boundary, then the solution must be specified.
c
c     ldimf
c        the row (or first) dimension of the arrays f,bdys,bdyf,bdzs,
c        and bdzf as it appears in the program calling hw3crt. this
c        parameter is used to specify the variable dimension of these
c        arrays.  ldimf must be at least l+1.
c
c     mdimf
c        the column (or second) dimension of the array f and the row (or
c        first) dimension of the arrays bdxs and bdxf as it appears in
c        the program calling hw3crt.  this parameter is used to specify
c        the variable dimension of these arrays.
c        mdimf must be at least m+1.
c
c     w
c        a one-dimensional array that must be provided by the user for
c        work space.  the length of w must be at least 30 + l + m + 5*n
c        + max(l,m,n) + 7*(int((l+1)/2) + int((m+1)/2))
c
c
c            * * * * * *   on output   * * * * * *
c
c     f
c        contains the solution u(i,j,k) of the finite difference
c        approximation for the grid point (x(i),y(j),z(k)) for
c        i=1,2,...,l+1, j=1,2,...,m+1, and k=1,2,...,n+1.
c
c     pertrb
c        if a combination of periodic or derivative boundary conditions
c        is specified for a poisson equation (lambda = 0), a solution
c        may not exist.  pertrb is a constant, calculated and subtracted
c        from f, which ensures that a solution exists.  pwscrt then
c        computes this solution, which is a least squares solution to
c        the original approximation.  this solution is not unique and is
c        unnormalized.  the value of pertrb should be small compared to
c        the right side f.  otherwise, a solution is obtained to an
c        essentially different problem.  this comparison should always
c        be made to insure that a meaningful solution has been obtained.
c
c     ierror
c        an error flag that indicates invalid input parameters.  except
c        for numbers 0 and 12, a solution is not attempted.
c
c        =  0  no error
c        =  1  xs .ge. xf
c        =  2  l .lt. 5
c        =  3  lbdcnd .lt. 0 .or. lbdcnd .gt. 4
c        =  4  ys .ge. yf
c        =  5  m .lt. 5
c        =  6  mbdcnd .lt. 0 .or. mbdcnd .gt. 4
c        =  7  zs .ge. zf
c        =  8  n .lt. 5
c        =  9  nbdcnd .lt. 0 .or. nbdcnd .gt. 4
c        = 10  ldimf .lt. l+1
c        = 11  mdimf .lt. m+1
c        = 12  lambda .gt. 0
c
c        since this is the only means of indicating a possibly incorrect
c        call to hw3crt, the user should test ierror after the call.
c
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bdxs(mdimf,n+1),bdxf(mdimf,n+1),bdys(ldimf,n+1),
c     arguments      bdyf(ldimf,n+1),bdzs(ldimf,m+1),bdzf(ldimf,m+1),
c                    f(ldimf,mdimf,n+1),w(see argument list)
c
c     latest         december 1, 1978
c     revision
c
c     subprograms    hw3crt,pois3d,pos3d1,trid,rffti,rfftf,rfftf1,
c     required       rfftb,rfftb1,costi,cost,sinti,sint,cosqi,cosqf,
c                    cosqf1,cosqb,cosqb1,sinqi,sinqf,sinqb,cffti,
c                    cffti1,cfftb,cfftb1,passb2,passb3,passb4,passb,
c                    cfftf,cfftf1,passf1,passf2,passf3,passf4,passf,
c                    pimach
c
c     special        none
c     conditions
c
c     common         value
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in july,1977
c
c     algorithm      this subroutine defines the finite difference
c                    equations, incorporates boundary data, and
c                    adjusts the right side of singular systems and
c                    then calls pois3d to solve the system.
c
c     space          7862(decimal) = 17300(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hw3crt is roughly proportional
c                    to l*m*n*(log2(l)+log2(m)+5), but also depends on
c                    input parameters lbdcnd and mbdcnd.  some typical
c                    values are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for l,m an
c                    n as large as 32.  more detailed information about
c                    accuracy can be found in the documentation for
c                    pois3d which is the routine that actuall
c                    solves the finite difference equations.
c
c
c                       l(=m=n)     lbdcnd(=mbdcnd=nbdcnd)      t(msecs)
c                       -------     ----------------------      --------
c
c                         16                  0                    300
c                         16                  1                    302
c                         16                  3                    348
c                         32                  0                   1925
c                         32                  1                   1929
c                         32                  3                   2109
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos,sin,atan
c     resident
c     routines
c
c     reference      none
c
c     required         cos,sin,atan
c     resident
c     routines
c
c     reference        none
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dimension       bdxs(mdimf,*)          ,bdxf(mdimf,*)          ,
     &                bdys(ldimf,*)          ,bdyf(ldimf,*)          ,
     &                bdzs(ldimf,*)          ,bdzf(ldimf,*)          ,
     &                f(ldimf,mdimf,*)       ,w(*)
c
c     check for invalid input.
c
      ierror = 0
      if (xf .le. xs) ierror = 1
      if (l .lt. 5) ierror = 2
      if (lbdcnd.lt.0 .or. lbdcnd.gt.4) ierror = 3
      if (yf .le. ys) ierror = 4
      if (m .lt. 5) ierror = 5
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 6
      if (zf .le. zs) ierror = 7
      if (n .lt. 5) ierror = 8
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 9
      if (ldimf .lt. l+1) ierror = 10
      if (mdimf .lt. m+1) ierror = 11
      if (ierror .ne. 0) go to 188
      dy = (yf-ys)/float(m)
      twbydy = 2./dy
      c2 = 1./(dy**2)
      mstart = 1
      mstop = m
      mp1 = m+1
      mp = mbdcnd+1
      go to (104,101,101,102,102),mp
  101 mstart = 2
  102 go to (104,104,103,103,104),mp
  103 mstop = mp1
  104 munk = mstop-mstart+1
      dz = (zf-zs)/float(n)
      twbydz = 2./dz
      np = nbdcnd+1
      c3 = 1./(dz**2)
      np1 = n+1
      nstart = 1
      nstop = n
      go to (108,105,105,106,106),np
  105 nstart = 2
  106 go to (108,108,107,107,108),np
  107 nstop = np1
  108 nunk = nstop-nstart+1
      lp1 = l+1
      dx = (xf-xs)/float(l)
      c1 = 1./(dx**2)
      twbydx = 2./dx
      lp = lbdcnd+1
      lstart = 1
      lstop = l
c
c     enter boundary data for x-boundaries.
c
      go to (122,109,109,112,112),lp
  109 lstart = 2
      do 111 j=mstart,mstop
         do 110 k=nstart,nstop
            f(2,j,k) = f(2,j,k)-c1*f(1,j,k)
  110    continue
  111 continue
      go to 115
  112 do 114 j=mstart,mstop
         do 113 k=nstart,nstop
            f(1,j,k) = f(1,j,k)+twbydx*bdxs(j,k)
  113    continue
  114 continue
  115 go to (122,116,119,119,116),lp
  116 do 118 j=mstart,mstop
         do 117 k=nstart,nstop
            f(l,j,k) = f(l,j,k)-c1*f(lp1,j,k)
  117    continue
  118 continue
      go to 122
  119 lstop = lp1
      do 121 j=mstart,mstop
         do 120 k=nstart,nstop
            f(lp1,j,k) = f(lp1,j,k)-twbydx*bdxf(j,k)
  120    continue
  121 continue
  122 lunk = lstop-lstart+1
c
c     enter boundary data for y-boundaries.
c
      go to (136,123,123,126,126),mp
  123 do 125 i=lstart,lstop
         do 124 k=nstart,nstop
            f(i,2,k) = f(i,2,k)-c2*f(i,1,k)
  124    continue
  125 continue
      go to 129
  126 do 128 i=lstart,lstop
         do 127 k=nstart,nstop
            f(i,1,k) = f(i,1,k)+twbydy*bdys(i,k)
  127    continue
  128 continue
  129 go to (136,130,133,133,130),mp
  130 do 132 i=lstart,lstop
         do 131 k=nstart,nstop
            f(i,m,k) = f(i,m,k)-c2*f(i,mp1,k)
  131    continue
  132 continue
      go to 136
  133 do 135 i=lstart,lstop
         do 134 k=nstart,nstop
            f(i,mp1,k) = f(i,mp1,k)-twbydy*bdyf(i,k)
  134    continue
  135 continue
  136 continue
c
c     enter boundary data for z-boundaries.
c
      go to (150,137,137,140,140),np
  137 do 139 i=lstart,lstop
         do 138 j=mstart,mstop
            f(i,j,2) = f(i,j,2)-c3*f(i,j,1)
  138    continue
  139 continue
      go to 143
  140 do 142 i=lstart,lstop
         do 141 j=mstart,mstop
            f(i,j,1) = f(i,j,1)+twbydz*bdzs(i,j)
  141    continue
  142 continue
  143 go to (150,144,147,147,144),np
  144 do 146 i=lstart,lstop
         do 145 j=mstart,mstop
            f(i,j,n) = f(i,j,n)-c3*f(i,j,np1)
  145    continue
  146 continue
      go to 150
  147 do 149 i=lstart,lstop
         do 148 j=mstart,mstop
            f(i,j,np1) = f(i,j,np1)-twbydz*bdzf(i,j)
  148    continue
  149 continue
c
c     define a,b,c coefficients in w-array.
c
  150 continue
      iwb = nunk+1
      iwc = iwb+nunk
      iww = iwc+nunk
      do 151 k=1,nunk
         i = iwc+k-1
         w(k) = c3
         w(i) = c3
         i = iwb+k-1
         w(i) = -2.*c3+elmbda
  151 continue
      go to (155,155,153,152,152),np
  152 w(iwc) = 2.*c3
  153 go to (155,155,154,154,155),np
  154 w(iwb-1) = 2.*c3
  155 continue
      pertrb = 0.
c
c     for singular problems adjust data to insure a solution will exist.
c
      go to (156,172,172,156,172),lp
  156 go to (157,172,172,157,172),mp
  157 go to (158,172,172,158,172),np
  158 if (elmbda) 172,160,159
  159 ierror = 12
      go to 172
  160 continue
      mstpm1 = mstop-1
      lstpm1 = lstop-1
      nstpm1 = nstop-1
      xlp = (2+lp)/3
      ylp = (2+mp)/3
      zlp = (2+np)/3
      s1 = 0.
      do 164 k=2,nstpm1
         do 162 j=2,mstpm1
            do 161 i=2,lstpm1
               s1 = s1+f(i,j,k)
  161       continue
            s1 = s1+(f(1,j,k)+f(lstop,j,k))/xlp
  162    continue
         s2 = 0.
         do 163 i=2,lstpm1
            s2 = s2+f(i,1,k)+f(i,mstop,k)
  163    continue
         s2 = (s2+(f(1,1,k)+f(1,mstop,k)+f(lstop,1,k)+f(lstop,mstop,k))/
     &                                                          xlp)/ylp
         s1 = s1+s2
  164 continue
      s = (f(1,1,1)+f(lstop,1,1)+f(1,1,nstop)+f(lstop,1,nstop)+
     &    f(1,mstop,1)+f(lstop,mstop,1)+f(1,mstop,nstop)+
     &                                   f(lstop,mstop,nstop))/(xlp*ylp)
      do 166 j=2,mstpm1
         do 165 i=2,lstpm1
            s = s+f(i,j,1)+f(i,j,nstop)
  165    continue
  166 continue
      s2 = 0.
      do 167 i=2,lstpm1
         s2 = s2+f(i,1,1)+f(i,1,nstop)+f(i,mstop,1)+f(i,mstop,nstop)
  167 continue
      s = s2/ylp+s
      s2 = 0.
      do 168 j=2,mstpm1
         s2 = s2+f(1,j,1)+f(1,j,nstop)+f(lstop,j,1)+f(lstop,j,nstop)
  168 continue
      s = s2/xlp+s
      pertrb = (s/zlp+s1)/((float(lunk+1)-xlp)*(float(munk+1)-ylp)*
     &                                              (float(nunk+1)-zlp))
      do 171 i=1,lunk
         do 170 j=1,munk
            do 169 k=1,nunk
               f(i,j,k) = f(i,j,k)-pertrb
  169       continue
  170    continue
  171 continue
  172 continue
      nperod = 0
      if (nbdcnd .eq. 0) go to 173
      nperod = 1
      w(1) = 0.
      w(iww-1) = 0.
  173 continue
      call pois3d (lbdcnd,lunk,c1,mbdcnd,munk,c2,nperod,nunk,w,w(iwb),
     &             w(iwc),ldimf,mdimf,f(lstart,mstart,nstart),ir,w(iww))
c
c     fill in sides for periodic boundary conditions.
c
      if (lp .ne. 1) go to 180
      if (mp .ne. 1) go to 175
      do 174 k=nstart,nstop
         f(1,mp1,k) = f(1,1,k)
  174 continue
      mstop = mp1
  175 if (np .ne. 1) go to 177
      do 176 j=mstart,mstop
         f(1,j,np1) = f(1,j,1)
  176 continue
      nstop = np1
  177 do 179 j=mstart,mstop
         do 178 k=nstart,nstop
            f(lp1,j,k) = f(1,j,k)
  178    continue
  179 continue
  180 continue
      if (mp .ne. 1) go to 185
      if (np .ne. 1) go to 182
      do 181 i=lstart,lstop
         f(i,1,np1) = f(i,1,1)
  181 continue
      nstop = np1
  182 do 184 i=lstart,lstop
         do 183 k=nstart,nstop
            f(i,mp1,k) = f(i,1,k)
  183    continue
  184 continue
  185 continue
      if (np .ne. 1) go to 188
      do 187 i=lstart,lstop
         do 186 j=mstart,mstop
            f(i,j,np1) = f(i,j,1)
  186    continue
  187 continue
  188 continue
      return
      end
      subroutine hwscrt (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HWSCRT 2D standard grid Cartesian coordinates five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c     hwscrt solves the standard five-point finite
c     difference approximation to the helmholtz equation in cartesian
c     coordinates:
c
c          (d/dx)(du/dx) + (d/dy)(du/dy) + lambda*u = f(x,y).
c
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     a,b
c       the range of x, i.e., a .le. x .le. b.  a must be less than b.
c
c     m
c       the number of panels into which the interval (a,b) is
c       subdivided.  hence, there will be m+1 grid points in the
c       x-direction given by x(i) = a+(i-1)dx for i = 1,2,...,m+1,
c       where dx = (b-a)/m is the panel width. m must be greater than 3.
c
c     mbdcnd
c       indicates the type of boundary conditions at x = a and x = b.
c
c       = 0  if the solution is periodic in x, i.e., u(i,j) = u(m+i,j).
c       = 1  if the solution is specified at x = a and x = b.
c       = 2  if the solution is specified at x = a and the derivative of
c            the solution with respect to x is specified at x = b.
c       = 3  if the derivative of the solution with respect to x is
c            specified at x = a and x = b.
c       = 4  if the derivative of the solution with respect to x is
c            specified at x = a and the solution is specified at x = b.
c
c     bda
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to x at x = a.
c       when mbdcnd = 3 or 4,
c
c            bda(j) = (d/dx)u(a,y(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to x at x = b.
c       when mbdcnd = 2 or 3,
c
c            bdb(j) = (d/dx)u(b,y(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value bdb is a dummy variable.
c
c     c,d
c       the range of y, i.e., c .le. y .le. d.  c must be less than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       y-direction given by y(j) = c+(j-1)dy for j = 1,2,...,n+1, where
c       dy = (d-c)/n is the panel width.  n must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at y = c and y = d.
c
c       = 0  if the solution is periodic in y, i.e., u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at y = c and y = d.
c       = 2  if the solution is specified at y = c and the derivative of
c            the solution with respect to y is specified at y = d.
c       = 3  if the derivative of the solution with respect to y is
c            specified at y = c and y = d.
c       = 4  if the derivative of the solution with respect to y is
c            specified at y = c and the solution is specified at y = d.
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to y at y = c.
c       when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dy)u(x(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to y at y = d.
c       when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dy)u(x(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscrt will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array which specifies the values of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(x(i),y(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd     f(1,j)        f(m+1,j)
c            ------     ---------     --------
c
c              0        f(a,y(j))     f(a,y(j))
c              1        u(a,y(j))     u(b,y(j))
c              2        u(a,y(j))     f(b,y(j))     j = 1,2,...,n+1
c              3        f(a,y(j))     f(b,y(j))
c              4        f(a,y(j))     u(b,y(j))
c
c
c            nbdcnd     f(i,1)        f(i,n+1)
c            ------     ---------     --------
c
c              0        f(x(i),c)     f(x(i),c)
c              1        u(x(i),c)     u(x(i),d)
c              2        u(x(i),c)     f(x(i),d)     i = 1,2,...,m+1
c              3        f(x(i),c)     f(x(i),d)
c              4        f(x(i),c)     u(x(i),d)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note
c
c       if the table calls for both the solution u and the right side f
c       at  a corner then the solution must be specified.
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwscrt.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwscrt and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (x(i),y(j)), i = 1,2,...,m+1,
c       j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic or derivative boundary conditions
c       is specified for a poisson equation (lambda = 0), a solution may
c       not exist.  pertrb is a constant, calculated and subtracted from
c       f, which ensures that a solution exists.  hwscrt then computes
c       this solution, which is a least squares solution to the original
c       approximation.  this solution plus any constant is also a
c       solution.  hence, the solution is not unique.  the value of
c       pertrb should be small compared to the right side f.  otherwise,
c       a solution is obtained to an essentially different problem.
c       this comparison should always be made to insure that a
c       meaningful solution has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 6, a solution is not attempted.
c
c       = 0  no error.
c       = 1  a .ge. b.
c       = 2  mbdcnd .lt. 0 or mbdcnd .gt. 4  .
c       = 3  c .ge. d.
c       = 4  n .le. 3
c       = 5  nbdcnd .lt. 0 or nbdcnd .gt. 4  .
c       = 6  lambda .gt. 0  .
c       = 7  idimf .lt. m+1  .
c       = 8  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwscrt, the user should test ierror after the call.
c
c     w
c       w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c
c     dimension of   bda(n+1),bdb(n+1),bdc(m+1),bdd(m+1),f(idimf,n+1),
c     arguments      w(see argument list)
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    hwscrt,genbun,poisd2,poisn2,poisp2,cosgen,merge,
c     required       trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized september 1, 1973
c                    revised april 1, 1976
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space          13110(octal) = 5704(decimal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwscrt is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameters nbdcnd and mbdcnd.  some typical values
c                    are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for n and
c                    m as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    genbun which is the routine that
c                    solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32        0         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        0         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     reference      swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     &                w(*)
c
c     check for invalid parameters.
c
      ierror = 0
      if (a .ge. b) ierror = 1
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 2
      if (c .ge. d) ierror = 3
      if (n .le. 3) ierror = 4
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 5
      if (idimf .lt. m+1) ierror = 7
      if (m .le. 3) ierror = 8
      if (ierror .ne. 0) return
      nperod = nbdcnd
      mperod = 0
      if (mbdcnd .gt. 0) mperod = 1
      deltax = (b-a)/float(m)
      twdelx = 2./deltax
      delxsq = 1./deltax**2
      deltay = (d-c)/float(n)
      twdely = 2./deltay
      delysq = 1./deltay**2
      np = nbdcnd+1
      np1 = n+1
      mp = mbdcnd+1
      mp1 = m+1
      nstart = 1
      nstop = n
      nskip = 1
      go to (104,101,102,103,104),np
  101 nstart = 2
      go to 104
  102 nstart = 2
  103 nstop = np1
      nskip = 2
  104 nunk = nstop-nstart+1
c
c     enter boundary data for x-boundaries.
c
      mstart = 1
      mstop = m
      mskip = 1
      go to (117,105,106,109,110),mp
  105 mstart = 2
      go to 107
  106 mstart = 2
      mstop = mp1
      mskip = 2
  107 do 108 j=nstart,nstop
         f(2,j) = f(2,j)-f(1,j)*delxsq
  108 continue
      go to 112
  109 mstop = mp1
      mskip = 2
  110 do 111 j=nstart,nstop
         f(1,j) = f(1,j)+bda(j)*twdelx
  111 continue
  112 go to (113,115),mskip
  113 do 114 j=nstart,nstop
         f(m,j) = f(m,j)-f(mp1,j)*delxsq
  114 continue
      go to 117
  115 do 116 j=nstart,nstop
         f(mp1,j) = f(mp1,j)-bdb(j)*twdelx
  116 continue
  117 munk = mstop-mstart+1
c
c     enter boundary data for y-boundaries.
c
      go to (127,118,118,120,120),np
  118 do 119 i=mstart,mstop
         f(i,2) = f(i,2)-f(i,1)*delysq
  119 continue
      go to 122
  120 do 121 i=mstart,mstop
         f(i,1) = f(i,1)+bdc(i)*twdely
  121 continue
  122 go to (123,125),nskip
  123 do 124 i=mstart,mstop
         f(i,n) = f(i,n)-f(i,np1)*delysq
  124 continue
      go to 127
  125 do 126 i=mstart,mstop
         f(i,np1) = f(i,np1)-bdd(i)*twdely
  126 continue
c
c    multiply right side by deltay**2.
c
  127 delysq = deltay*deltay
      do 129 i=mstart,mstop
         do 128 j=nstart,nstop
            f(i,j) = f(i,j)*delysq
  128    continue
  129 continue
c
c     define the a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      s = delysq*delxsq
      st2 = 2.*s
      do 130 i=1,munk
         w(i) = s
         j = id2+i
         w(j) = -st2+elmbda*delysq
         j = id3+i
         w(j) = s
  130 continue
      if (mp .eq. 1) go to 131
      w(1) = 0.
      w(id4) = 0.
  131 continue
      go to (135,135,132,133,134),mp
  132 w(id2) = st2
      go to 135
  133 w(id2) = st2
  134 w(id3+1) = st2
  135 continue
      pertrb = 0.
      if (elmbda) 144,137,136
  136 ierror = 6
      go to 144
  137 if ((nbdcnd.eq.0 .or. nbdcnd.eq.3) .and.
     &    (mbdcnd.eq.0 .or. mbdcnd.eq.3)) go to 138
      go to 144
c
c     for singular problems must adjust data to insure that a solution
c     will exist.
c
  138 a1 = 1.
      a2 = 1.
      if (nbdcnd .eq. 3) a2 = 2.
      if (mbdcnd .eq. 3) a1 = 2.
      s1 = 0.
      msp1 = mstart+1
      mstm1 = mstop-1
      nsp1 = nstart+1
      nstm1 = nstop-1
      do 140 j=nsp1,nstm1
         s = 0.
         do 139 i=msp1,mstm1
            s = s+f(i,j)
  139    continue
         s1 = s1+s*a1+f(mstart,j)+f(mstop,j)
  140 continue
      s1 = a2*s1
      s = 0.
      do 141 i=msp1,mstm1
         s = s+f(i,nstart)+f(i,nstop)
  141 continue
      s1 = s1+s*a1+f(mstart,nstart)+f(mstart,nstop)+f(mstop,nstart)+
     &     f(mstop,nstop)
      s = (2.+float(nunk-2)*a2)*(2.+float(munk-2)*a1)
      pertrb = s1/s
      do 143 j=nstart,nstop
         do 142 i=mstart,mstop
            f(i,j) = f(i,j)-pertrb
  142    continue
  143 continue
      pertrb = pertrb/delysq
c
c     solve the equation.
c
  144 call genbun (nperod,nunk,mperod,munk,w(1),w(id2+1),w(id3+1),
     &             idimf,f(mstart,nstart),ierr1,w(id4+1))
      w(1) = w(id4+1)+3.*float(munk)
c
c     fill in identical values when have periodic boundary conditions.
c
      if (nbdcnd .ne. 0) go to 146
      do 145 i=mstart,mstop
         f(i,np1) = f(i,1)
  145 continue
  146 if (mbdcnd .ne. 0) go to 148
      do 147 j=nstart,nstop
         f(mp1,j) = f(1,j)
  147 continue
      if (nbdcnd .eq. 0) f(mp1,np1) = f(1,np1)
  148 continue
      return
      end
      subroutine hwscs1 (intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,
     &                   bdrs,bdrf,elmbda,f,idimf,pertrb,w,s,an,bn,cn,
     &                   r,am,bm,cm,sint,bmh)

c*********************************************************************72
c
cc HWSCS1 is a utility routine for HWSCSP.
c
      dimension       f(idimf,*) ,bdrs(*)    ,bdrf(*)    ,bdts(*)    ,
     &                bdtf(*)    ,am(*)      ,bm(*)      ,cm(*)      ,
     &                an(*)      ,bn(*)      ,cn(*)      ,s(*)       ,
     &                r(*)       ,sint(*)    ,bmh(*)     ,w(*)
      pi = pimach()
      eps = epmach()
      mp1 = m+1
      dth = (tf-ts)/float(m)
      tdt = dth+dth
      hdth = dth/2.
      sdts = 1./(dth*dth)
      do 102 i=1,mp1
         theta = ts+float(i-1)*dth
         sint(i) = sin(theta)
         if (sint(i)) 101,102,101
  101    t1 = sdts/sint(i)
         am(i) = t1*sin(theta-hdth)
         cm(i) = t1*sin(theta+hdth)
         bm(i) = -(am(i)+cm(i))
  102 continue
      np1 = n+1
      dr = (rf-rs)/float(n)
      hdr = dr/2.
      tdr = dr+dr
      dr2 = dr*dr
      czr = 6.*dth/(dr2*(cos(ts)-cos(tf)))
      do 103 j=1,np1
         r(j) = rs+float(j-1)*dr
         an(j) = (r(j)-hdr)**2/dr2
         cn(j) = (r(j)+hdr)**2/dr2
         bn(j) = -(an(j)+cn(j))
  103 continue
      mp = 1
      np = 1
c
c boundary condition at phi=ps
c
      go to (104,104,105,105,106,106,104,105,106),mbdcnd
  104 at = am(2)
      its = 2
      go to 107
  105 at = am(1)
      its = 1
      cm(1) = cm(1)+am(1)
      go to 107
  106 its = 1
      bm(1) = -4.*sdts
      cm(1) = -bm(1)
c
c boundary condition at phi=pf
c
  107 go to (108,109,109,108,108,109,110,110,110),mbdcnd
  108 ct = cm(m)
      itf = m
      go to 111
  109 ct = cm(m+1)
      am(m+1) = am(m+1)+cm(m+1)
      itf = m+1
      go to 111
  110 itf = m+1
      am(m+1) = 4.*sdts
      bm(m+1) = -am(m+1)
  111 wts = sint(its+1)*am(its+1)/cm(its)
      wtf = sint(itf-1)*cm(itf-1)/am(itf)
      itsp = its+1
      itfm = itf-1
c
c boundary condition at r=rs
c
      ictr = 0
      go to (112,112,113,113,114,114),nbdcnd
  112 ar = an(2)
      jrs = 2
      go to 118
  113 ar = an(1)
      jrs = 1
      cn(1) = cn(1)+an(1)
      go to 118
  114 jrs = 2
      ictr = 1
      s(n) = an(n)/bn(n)
      do 115 j=3,n
         l = n-j+2
         s(l) = an(l)/(bn(l)-cn(l)*s(l+1))
  115 continue
      s(2) = -s(2)
      do 116 j=3,n
         s(j) = -s(j)*s(j-1)
  116 continue
      wtnm = wts+wtf
      do 117 i=itsp,itfm
         wtnm = wtnm+sint(i)
  117 continue
      yps = czr*wtnm*(s(2)-1.)
c
c boundary condition at r=rf
c
  118 go to (119,120,120,119,119,120),nbdcnd
  119 cr = cn(n)
      jrf = n
      go to 121
  120 cr = cn(n+1)
      an(n+1) = an(n+1)+cn(n+1)
      jrf = n+1
  121 wrs = an(jrs+1)*r(jrs)**2/cn(jrs)
      wrf = cn(jrf-1)*r(jrf)**2/an(jrf)
      wrz = an(jrs)/czr
      jrsp = jrs+1
      jrfm = jrf-1
      munk = itf-its+1
      nunk = jrf-jrs+1
      do 122 i=its,itf
         bmh(i) = bm(i)
  122 continue
      ising = 0
      go to (132,132,123,132,132,123),nbdcnd
  123 go to (132,132,124,132,132,124,132,124,124),mbdcnd
  124 if (elmbda) 132,125,125
  125 ising = 1
      sum = wts*wrs+wts*wrf+wtf*wrs+wtf*wrf
      if (ictr) 126,127,126
  126 sum = sum+wrz
  127 do 129 j=jrsp,jrfm
         r2 = r(j)**2
         do 128 i=itsp,itfm
            sum = sum+r2*sint(i)
  128    continue
  129 continue
      do 130 j=jrsp,jrfm
         sum = sum+(wts+wtf)*r(j)**2
  130 continue
      do 131 i=itsp,itfm
         sum = sum+(wrs+wrf)*sint(i)
  131 continue
      hne = sum
  132 go to (133,133,133,133,134,134,133,133,134),mbdcnd
  133 bm(its) = bmh(its)+elmbda/sint(its)**2
  134 go to (135,135,135,135,135,135,136,136,136),mbdcnd
  135 bm(itf) = bmh(itf)+elmbda/sint(itf)**2
  136 do 137 i=itsp,itfm
         bm(i) = bmh(i)+elmbda/sint(i)**2
  137 continue
      go to (138,138,140,140,142,142,138,140,142),mbdcnd
  138 do 139 j=jrs,jrf
         f(2,j) = f(2,j)-at*f(1,j)/r(j)**2
  139 continue
      go to 142
  140 do 141 j=jrs,jrf
         f(1,j) = f(1,j)+tdt*bdts(j)*at/r(j)**2
  141 continue
  142 go to (143,145,145,143,143,145,147,147,147),mbdcnd
  143 do 144 j=jrs,jrf
         f(m,j) = f(m,j)-ct*f(m+1,j)/r(j)**2
  144 continue
      go to 147
  145 do 146 j=jrs,jrf
         f(m+1,j) = f(m+1,j)-tdt*bdtf(j)*ct/r(j)**2
  146 continue
  147 go to (151,151,153,153,148,148),nbdcnd
  148 if (mbdcnd-3) 155,149,155
  149 yhld = f(its,1)-czr/tdt*(sin(tf)*bdtf(2)-sin(ts)*bdts(2))
      do 150 i=1,mp1
         f(i,1) = yhld
  150 continue
      go to 155
  151 rs2 = (rs+dr)**2
      do 152 i=its,itf
         f(i,2) = f(i,2)-ar*f(i,1)/rs2
  152 continue
      go to 155
  153 do 154 i=its,itf
         f(i,1) = f(i,1)+tdr*bdrs(i)*ar/rs**2
  154 continue
  155 go to (156,158,158,156,156,158),nbdcnd
  156 rf2 = (rf-dr)**2
      do 157 i=its,itf
         f(i,n) = f(i,n)-cr*f(i,n+1)/rf2
  157 continue
      go to 160
  158 do 159 i=its,itf
         f(i,n+1) = f(i,n+1)-tdr*bdrf(i)*cr/rf**2
  159 continue
  160 continue
      pertrb = 0.
      if (ising) 161,170,161
  161 sum = wts*wrs*f(its,jrs)+wts*wrf*f(its,jrf)+wtf*wrs*f(itf,jrs)+
     &      wtf*wrf*f(itf,jrf)
      if (ictr) 162,163,162
  162 sum = sum+wrz*f(its,1)
  163 do 165 j=jrsp,jrfm
         r2 = r(j)**2
         do 164 i=itsp,itfm
            sum = sum+r2*sint(i)*f(i,j)
  164    continue
  165 continue
      do 166 j=jrsp,jrfm
         sum = sum+r(j)**2*(wts*f(its,j)+wtf*f(itf,j))
  166 continue
      do 167 i=itsp,itfm
         sum = sum+sint(i)*(wrs*f(i,jrs)+wrf*f(i,jrf))
  167 continue
      pertrb = sum/hne
      do 169 j=1,np1
         do 168 i=1,mp1
            f(i,j) = f(i,j)-pertrb
  168    continue
  169 continue
  170 do 172 j=jrs,jrf
         rsq = r(j)**2
         do 171 i=its,itf
            f(i,j) = rsq*f(i,j)
  171    continue
  172 continue
      iflg = intl
  173 call blktri (iflg,np,nunk,an(jrs),bn(jrs),cn(jrs),mp,munk,
     &             am(its),bm(its),cm(its),idimf,f(its,jrs),ierror,w)
      iflg = iflg+1
      if (iflg-1) 174,173,174
  174 if (nbdcnd) 177,175,177
  175 do 176 i=1,mp1
         f(i,jrf+1) = f(i,jrs)
  176 continue
  177 if (mbdcnd) 180,178,180
  178 do 179 j=1,np1
         f(itf+1,j) = f(its,j)
  179 continue
  180 xp = 0.
      if (ictr) 181,188,181
  181 if (ising) 186,182,186
  182 sum = wts*f(its,2)+wtf*f(itf,2)
      do 183 i=itsp,itfm
         sum = sum+sint(i)*f(i,2)
  183 continue
      yph = czr*sum
      xp = (f(its,1)-yph)/yps
      do 185 j=jrs,jrf
         xps = xp*s(j)
         do 184 i=its,itf
            f(i,j) = f(i,j)+xps
  184    continue
  185 continue
  186 do 187 i=1,mp1
         f(i,1) = xp
  187 continue
  188 return
      end
      subroutine hwscsp (intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,
     &                   bdrs,bdrf,elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HWSCSP 2D standard grid axisymmetric spherical coordinates, five point scheme.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       hwscsp solves a finite difference approximation to the
c       modified helmholtz equation in spherical coordinates assuming
c       axisymmetry  (no dependence on longitude)
c
c          (1/r**2)(d/dr)((r**2)(d/dr)u)
c
c             + (1/(r**2)sin(theta))(d/dtheta)(sin(theta)(d/dtheta)u)
c
c             + (lambda/(rsin(theta))**2)u = f(theta,r).
c
c     this two dimensional modified helmholtz equation results from
c     the fourier transform of the three dimensional poisson equation
c
c     * * * * * * * * * *     on input     * * * * * * * * * *
c
c     intl
c       = 0  on initial entry to hwscsp or if any of the arguments
c            rs, rf, n, nbdcnd are changed from a previous call.
c       = 1  if rs, rf, n, nbdcnd are all unchanged from previous call
c            to hwscsp.
c
c       note   a call with intl=0 takes approximately 1.5 times as
c              much time as a call with intl = 1  .  once a call with
c              intl = 0 has been made then subsequent solutions
c              corresponding to different f, bdts, bdtf, bdrs, bdrf can
c              be obtained faster with intl = 1 since initialization is
c              not repeated.
c
c     ts,tf
c       the range of theta (colatitude), i.e., ts .le. theta .le. tf.
c       ts must be less than tf.  ts and tf are in radians.  a ts of
c       zero corresponds to the north pole and a tf of pi corresponds
c       to the south pole.
c
c     * * * * * * * * * * * * * * important * * * * * * * * * * * * * *
c
c     if tf is equal to pi then it must be computed using the statement
c     tf = pimach(). this insures that tf in the users program is
c     equal to pi in this program which permits several tests of the
c     input parameters that otherwise would not be possible.
c
c     m
c       the number of panels into which the interval (ts,tf) is
c       subdivided.  hence, there will be m+1 grid points in the
c       theta-direction given by theta(k) = (i-1)dtheta+ts for
c       i = 1,2,...,m+1, where dtheta = (tf-ts)/m is the panel width.
c
c     mbdcnd
c       indicates the type of boundary condition at theta = ts and
c       theta = tf.
c
c       = 1  if the solution is specified at theta = ts and theta = tf.
c       = 2  if the solution is specified at theta = ts and the
c            derivative of the solution with respect to theta is
c            specified at theta = tf (see note 2 below).
c       = 3  if the derivative of the solution with respect to theta is
c            specified at theta = ts and theta = tf (see notes 1,2
c            below).
c       = 4  if the derivative of the solution with respect to theta is
c            specified at theta = ts (see note 1 below) and the
c            solution is specified at theta = tf.
c       = 5  if the solution is unspecified at theta = ts = 0 and the
c            solution is specified at theta = tf.
c       = 6  if the solution is unspecified at theta = ts = 0 and the
c            derivative of the solution with respect to theta is
c            specified at theta = tf (see note 2 below).
c       = 7  if the solution is specified at theta = ts and the
c            solution is unspecified at theta = tf = pi.
c       = 8  if the derivative of the solution with respect to theta is
c            specified at theta = ts (see note 1 below) and the solution
c            is unspecified at theta = tf = pi.
c       = 9  if the solution is unspecified at theta = ts = 0 and
c            theta = tf = pi.
c
c       notes:  1.  if ts = 0, do not use mbdcnd = 3,4, or 8, but
c                   instead use mbdcnd = 5,6, or 9  .
c               2.  if tf = pi, do not use mbdcnd = 2,3, or 6, but
c                   instead use mbdcnd = 7,8, or 9  .
c
c     bdts
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = ts.  when mbdcnd = 3,4, or 8,
c
c            bdts(j) = (d/dtheta)u(ts,r(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdts is a dummy variable.
c
c     bdtf
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = tf.  when mbdcnd = 2,3, or 6,
c
c            bdtf(j) = (d/dtheta)u(tf,r(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdtf is a dummy variable.
c
c     rs,rf
c       the range of r, i.e., rs .le. r .lt. rf.  rs must be less than
c       rf.  rs must be non-negative.
c
c       n
c       the number of panels into which the interval (rs,rf) is
c       subdivided.  hence, there will be n+1 grid points in the
c       r-direction given by r(j) = (j-1)dr+rs for j = 1,2,...,n+1,
c       where dr = (rf-rs)/n is the panel width.
c       n must be greater than 2
c
c     nbdcnd
c       indicates the type of boundary condition at r = rs and r = rf.
c
c       = 1  if the solution is specified at r = rs and r = rf.
c       = 2  if the solution is specified at r = rs and the derivative
c            of the solution with respect to r is specified at r = rf.
c       = 3  if the derivative of the solution with respect to r is
c            specified at r = rs and r = rf.
c       = 4  if the derivative of the solution with respect to r is
c            specified at rs and the solution is specified at r = rf.
c       = 5  if the solution is unspecified at r = rs = 0 (see note
c            below) and the solution is specified at r = rf.
c       = 6  if the solution is unspecified at r = rs = 0 (see note
c            below) and the derivative of the solution with respect to
c            r is specified at r = rf.
c
c       note:  nbdcnd = 5 or 6 cannot be used with
c              mbdcnd = 1,2,4,5, or 7 (the former indicates that the
c                       solution is unspecified at r = 0, the latter
c                       indicates that the solution is specified).
c                       use instead
c              nbdcnd = 1 or 2  .
c
c     bdrs
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to r at r = rs.
c       when nbdcnd = 3 or 4,
c
c            bdrs(i) = (d/dr)u(theta(i),rs), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdrs is a dummy variable.
c
c     bdrf
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to r at r = rf.
c       when nbdcnd = 2,3, or 6,
c
c            bdrf(i) = (d/dr)u(theta(i),rf), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdrf is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscsp will
c       attempt to find a solution.  if nbdcnd = 5 or 6 or
c       mbdcnd = 5,6,7,8, or 9, elmbda must be zero.
c
c     f
c       a two-dimensional array that specifies the value of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(theta(i),r(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   ----------        ----------
c
c              1      u(ts,r(j))        u(tf,r(j))
c              2      u(ts,r(j))        f(tf,r(j))
c              3      f(ts,r(j))        f(tf,r(j))
c              4      f(ts,r(j))        u(tf,r(j))
c              5      f(0,r(j))         u(tf,r(j))   j = 1,2,...,n+1
c              6      f(0,r(j))         f(tf,r(j))
c              7      u(ts,r(j))        f(pi,r(j))
c              8      f(ts,r(j))        f(pi,r(j))
c              9      f(0,r(j))         f(pi,r(j))
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   --------------    --------------
c
c              1      u(theta(i),rs)    u(theta(i),rf)
c              2      u(theta(i),rs)    f(theta(i),rf)
c              3      f(theta(i),rs)    f(theta(i),rf)
c              4      f(theta(i),rs)    u(theta(i),rf)   i = 1,2,...,m+1
c              5      f(ts,0)           u(theta(i),rf)
c              6      f(ts,0)           f(theta(i),rf)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note
c
c       if the table calls for both the solution u and the right side f
c       at  a corner then the solution must be specified.
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwscsp.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space. its length can be computed from the formula below
c       which depends on the value of nbdcnd
c
c       if nbdcnd=2,4 or 6 define nunk=n
c       if nbdcnd=1 or 5   define nunk=n-1
c       if nbdcnd=3        define nunk=n+1
c
c       now set k=int(log2(nunk))+1 and l=2**(k+1) then w must be
c       dimensioned at least (k-2)*l+k+5*(m+n)+max(2*n,6*m)+23
c
c       **important** for purposes of checking, the required length
c                     of w is computed by hwscsp and stored in w(1)
c                     in floating point format.
c
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (theta(i),r(j)),
c       i = 1,2,...,m+1,   j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic or derivative boundary conditions
c       is specified for a poisson equation (lambda = 0), a solution may
c       not exist.  pertrb is a constant, calculated and subtracted from
c       f, which ensures that a solution exists.  hwscsp then computes
c       this solution, which is a least squares solution to the original
c       approximation. this solution is not unique and is unnormalized.
c       the value of pertrb should be small compared to the right side
c       f. otherwise , a solution is obtained to an essentially
c       different problem. this comparison should always be made to
c       insure that a meaningful solution has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 10, a solution is not attempted.
c
c       = 1  ts.lt.0. or tf.gt.pi
c       = 2  ts.ge.tf
c       = 3  m.lt.5
c       = 4  mbdcnd.lt.1 or mbdcnd.gt.9
c       = 5  rs.lt.0
c       = 6  rs.ge.rf
c       = 7  n.lt.5
c       = 8  nbdcnd.lt.1 or nbdcnd.gt.6
c       = 9  elmbda.gt.0
c       = 10 idimf.lt.m+1
c       = 11 elmbda.ne.0 and mbdcnd.ge.5
c       = 12 elmbda.ne.0 and nbdcnd equals 5 or 6
c       = 13 mbdcnd equals 5,6 or 9 and ts.ne.0
c       = 14 mbdcnd.ge.7 and tf.ne.pi
c       = 15 ts.eq.0 and mbdcnd equals 3,4 or 8
c       = 16 tf.eq.pi and mbdcnd equals 2,3 or 6
c       = 17 nbdcnd.ge.5 and rs.ne.0
c       = 18 nbdcnd.ge.5 and mbdcnd equals 1,2,4,5 or 7
c
c       since this is the only means of indicating a possliby incorrect
c       call to hwscsp, the user should test ierror after a call.
c
c     w
c       contains intermediate values that must not be destroyed if
c       hwscsp will be called again with intl = 1.  w(1) contains the
c       number of locations which w must have
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bdts(n+1),bdtf(n+1),bdrs(m+1),bdrf(m+1),
c     arguments      f(idimf,n+1),w(see argument list)
c
c     latest         june 1979
c     revision
c
c     subprograms    hwscsp,hwscs1,blktri,blktr1,prod,prodp,cprod,cprodp
c     required       ,combp,ppadd,psgf,bsrh,ppsgf,ppspf,tevls,indxa,
c                    ,indxb,indxc,epmach,store
c
c     special
c     conditions
c
c     common         cblkt,value
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     paul n swarztrauber
c
c     language       fortran
c
c     history        version 1 september 1973
c                    version 2 april     1976
c                    version 3 june      1979
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    blktri to solve the system.
c
c     space
c     required
c
c     portability    american national standards institute fortran.
c                    the machine accuracy is computed approximately
c                    in function epmach
c
c     required       none
c     resident
c     routines
c
c     reference      swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
      dimension       f(idimf,*) ,bdts(*)    ,bdtf(*)    ,bdrs(*)    ,
     &                bdrf(*)    ,w(*)
      pi = pimach()
      ierror = 0
      if (ts.lt.0. .or. tf.gt.pi) ierror = 1
      if (ts .ge. tf) ierror = 2
      if (m .lt. 5) ierror = 3
      if (mbdcnd.lt.1 .or. mbdcnd.gt.9) ierror = 4
      if (rs .lt. 0.) ierror = 5
      if (rs .ge. rf) ierror = 6
      if (n .lt. 5) ierror = 7
      if (nbdcnd.lt.1 .or. nbdcnd.gt.6) ierror = 8
      if (elmbda .gt. 0.) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if (elmbda.ne.0. .and. mbdcnd.ge.5) ierror = 11
      if (elmbda.ne.0. .and. (nbdcnd.eq.5 .or. nbdcnd.eq.6)) ierror = 12
      if ((mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9) .and.
     &    ts.ne.0.) ierror = 13
      if (mbdcnd.ge.7 .and. tf.ne.pi) ierror = 14
      if (ts.eq.0. .and.
     &    (mbdcnd.eq.4 .or. mbdcnd.eq.8 .or. mbdcnd.eq.3)) ierror = 15
      if (tf.eq.pi .and.
     &    (mbdcnd.eq.2 .or. mbdcnd.eq.3 .or. mbdcnd.eq.6)) ierror = 16
      if (nbdcnd.ge.5 .and. rs.ne.0.) ierror = 17
      if (nbdcnd.ge.5 .and. (mbdcnd.eq.1 .or. mbdcnd.eq.2 .or.
     &                                    mbdcnd.eq.5 .or. mbdcnd.eq.7))
     &    ierror = 18
      if (ierror.ne.0 .and. ierror.ne.9) return
      nck = n
      go to (101,103,102,103,101,103),nbdcnd
  101 nck = nck-1
      go to 103
  102 nck = nck+1
  103 l = 2
      k = 1
  104 l = l+l
      k = k+1
      if (nck-l) 105,105,104
  105 l = l+l
      np1 = n+1
      mp1 = m+1
      i1 = (k-2)*l+k+max0(2*n,6*m)+13
      i2 = i1+np1
      i3 = i2+np1
      i4 = i3+np1
      i5 = i4+np1
      i6 = i5+np1
      i7 = i6+mp1
      i8 = i7+mp1
      i9 = i8+mp1
      i10 = i9+mp1
      w(1) = float(i10+m)
      call hwscs1 (intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,
     &             bdrf,elmbda,f,idimf,pertrb,w(2),w(i1),w(i2),w(i3),
     &             w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10))
      return
      end
      subroutine hwscyl (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HWSCYL 2D standard grid cylindrical coordinates five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c
c     hwscyl solves a finite difference approximation to the
c     helmholtz equation in cylindrical coordinates:
c
c          (1/r)(d/dr)(r(du/dr)) + (d/dz)(du/dz)
c
c                                + (lambda/r**2)u = f(r,z)
c
c     this modified helmholtz equation results from the fourier
c     transform of the three-dimensional poisson equation.
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     a,b
c       the range of r, i.e., a .le. r .le. b.  a must be less than b
c       and a must be non-negative.
c
c     m
c       the number of panels into which the interval (a,b) is
c       subdivided.  hence, there will be m+1 grid points in the
c       r-direction given by r(i) = a+(i-1)dr, for i = 1,2,...,m+1,
c       where dr = (b-a)/m is the panel width. m must be greater than 3.
c
c     mbdcnd
c       indicates the type of boundary conditions at r = a and r = b.
c
c       = 1  if the solution is specified at r = a and r = b.
c       = 2  if the solution is specified at r = a and the derivative of
c            the solution with respect to r is specified at r = b.
c       = 3  if the derivative of the solution with respect to r is
c            specified at r = a (see note below) and r = b.
c       = 4  if the derivative of the solution with respect to r is
c            specified at r = a (see note below) and the solution is
c            specified at r = b.
c       = 5  if the solution is unspecified at r = a = 0 and the
c            solution is specified at r = b.
c       = 6  if the solution is unspecified at r = a = 0 and the
c            derivative of the solution with respect to r is specified
c            at r = b.
c
c       note:  if a = 0, do not use mbdcnd = 3 or 4, but instead use
c              mbdcnd = 1,2,5, or 6  .
c
c     bda
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = a.
c       when mbdcnd = 3 or 4,
c
c            bda(j) = (d/dr)u(a,z(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = b.
c       when mbdcnd = 2,3, or 6,
c
c            bdb(j) = (d/dr)u(b,z(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdb is a dummy variable.
c
c     c,d
c       the range of z, i.e., c .le. z .le. d.  c must be less than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       z-direction given by z(j) = c+(j-1)dz, for j = 1,2,...,n+1,
c       where dz = (d-c)/n is the panel width. n must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at z = c and z = d.
c
c       = 0  if the solution is periodic in z, i.e., u(i,1) = u(i,n+1).
c       = 1  if the solution is specified at z = c and z = d.
c       = 2  if the solution is specified at z = c and the derivative of
c            the solution with respect to z is specified at z = d.
c       = 3  if the derivative of the solution with respect to z is
c            specified at z = c and z = d.
c       = 4  if the derivative of the solution with respect to z is
c            specified at z = c and the solution is specified at z = d.
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to z at z = c.
c       when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dz)u(r(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to z at z = d.
c       when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dz)u(r(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscyl will
c       attempt to find a solution.  lambda must be zero when
c       mbdcnd = 5 or 6  .
c
c     f
c       a two-dimensional array that specifies the values of the right
c       side of the helmholtz equation and boundary data (if any).  for
c       i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(r(i),z(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   ---------         ---------
c
c              1      u(a,z(j))         u(b,z(j))
c              2      u(a,z(j))         f(b,z(j))
c              3      f(a,z(j))         f(b,z(j))   j = 1,2,...,n+1
c              4      f(a,z(j))         u(b,z(j))
c              5      f(0,z(j))         u(b,z(j))
c              6      f(0,z(j))         f(b,z(j))
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   ---------         ---------
c
c              0      f(r(i),c)         f(r(i),c)
c              1      u(r(i),c)         u(r(i),d)
c              2      u(r(i),c)         f(r(i),d)   i = 1,2,...,m+1
c              3      f(r(i),c)         f(r(i),d)
c              4      f(r(i),c)         u(r(i),d)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note
c
c       if the table calls for both the solution u and the right side f
c       at  a corner then the solution must be specified.
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwscyl.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwscyl and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (r(i),z(j)), i = 1,2,...,m+1,
c       j = 1,2,...,n+1  .
c
c     pertrb
c       if one specifies a combination of periodic, derivative, and
c       unspecified boundary conditions for a poisson equation
c       (lambda = 0), a solution may not exist.  pertrb is a constant,
c       calculated and subtracted from f, which ensures that a solution
c       exists.  hwscyl then computes this solution, which is a least
c       squares solution to the original approximation.  this solution
c       plus any constant is also a solution.  hence, the solution is
c       not unique.  the value of pertrb should be small compared to the
c       right side f.  otherwise, a solution is obtained to an
c       essentially different problem.  this comparison should always
c       be made to insure that a meaningful solution has been obtained.
c
c     ierror
c       an error flag which indicates invalid input parameters.  except
c       for numbers 0 and 11, a solution is not attempted.
c
c       =  0  no error.
c       =  1  a .lt. 0  .
c       =  2  a .ge. b.
c       =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6  .
c       =  4  c .ge. d.
c       =  5  n .le. 3
c       =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4  .
c       =  7  a = 0, mbdcnd = 3 or 4  .
c       =  8  a .gt. 0, mbdcnd .ge. 5  .
c       =  9  a = 0, lambda .ne. 0, mbdcnd .ge. 5  .
c       = 10  idimf .lt. m+1  .
c       = 11  lambda .gt. 0  .
c       = 12  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwscyl, the user should test ierror after the call.
c
c     w
c       w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n+1),bdb(n+1),bdc(m+1),bdd(m+1),f(idimf,n+1),
c     arguments      w(see argument list)
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    hwscyl,genbun,poisd2,poisn2,poisp2,cosgen,merge,
c     required      trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized september 1, 1973
c                    revised april 1, 1976
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space          5818(decimal) = 13272(octal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwscyl is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameters nbdcnd and mbdcnd.  some typical values
c                    are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for n and
c                    m as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    genbun which is the routine that
c                    solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32        1         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        1         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     &                w(*)
c
c     check for invalid parameters.
c
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 3) ierror = 5
      if (nbdcnd.le.-1 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4)) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (a.eq.0. .and. elmbda.ne.0. .and. mbdcnd.ge.5) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if (m .le. 3) ierror = 12
      if (ierror .ne. 0) return
      mp1 = m+1
      deltar = (b-a)/float(m)
      dlrby2 = deltar/2.
      dlrsq = deltar**2
      np1 = n+1
      deltht = (d-c)/float(n)
      dlthsq = deltht**2
      np = nbdcnd+1
c
c     define range of indices i and j for unknowns u(i,j).
c
      mstart = 2
      mstop = m
      go to (104,103,102,101,101,102),mbdcnd
  101 mstart = 1
      go to 104
  102 mstart = 1
  103 mstop = mp1
  104 munk = mstop-mstart+1
      nstart = 1
      nstop = n
      go to (108,105,106,107,108),np
  105 nstart = 2
      go to 108
  106 nstart = 2
  107 nstop = np1
  108 nunk = nstop-nstart+1
c
c     define a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      id5 = id4+munk
      id6 = id5+munk
      istart = 1
      a1 = 2./dlrsq
      ij = 0
      if (mbdcnd.eq.3 .or. mbdcnd.eq.4) ij = 1
      if (mbdcnd .le. 4) go to 109
      w(1) = 0.
      w(id2+1) = -2.*a1
      w(id3+1) = 2.*a1
      istart = 2
      ij = 1
  109 do 110 i=istart,munk
         r = a+float(i-ij)*deltar
         j = id5+i
         w(j) = r
         j = id6+i
         w(j) = 1./r**2
         w(i) = (r-dlrby2)/(r*dlrsq)
         j = id3+i
         w(j) = (r+dlrby2)/(r*dlrsq)
         k = id6+i
         j = id2+i
         w(j) = -a1+elmbda*w(k)
  110 continue
      go to (114,111,112,113,114,112),mbdcnd
  111 w(id2) = a1
      go to 114
  112 w(id2) = a1
  113 w(id3+1) = a1*float(istart)
  114 continue
c
c     enter boundary data for r-boundaries.
c
      go to (115,115,117,117,119,119),mbdcnd
  115 a1 = w(1)
      do 116 j=nstart,nstop
         f(2,j) = f(2,j)-a1*f(1,j)
  116 continue
      go to 119
  117 a1 = 2.*deltar*w(1)
      do 118 j=nstart,nstop
         f(1,j) = f(1,j)+a1*bda(j)
  118 continue
  119 go to (120,122,122,120,120,122),mbdcnd
  120 a1 = w(id4)
      do 121 j=nstart,nstop
         f(m,j) = f(m,j)-a1*f(mp1,j)
  121 continue
      go to 124
  122 a1 = 2.*deltar*w(id4)
      do 123 j=nstart,nstop
         f(mp1,j) = f(mp1,j)-a1*bdb(j)
  123 continue
c
c     enter boundary data for z-boundaries.
c
  124 a1 = 1./dlthsq
      l = id5-mstart+1
      go to (134,125,125,127,127),np
  125 do 126 i=mstart,mstop
         f(i,2) = f(i,2)-a1*f(i,1)
  126 continue
      go to 129
  127 a1 = 2./deltht
      do 128 i=mstart,mstop
         f(i,1) = f(i,1)+a1*bdc(i)
  128 continue
  129 a1 = 1./dlthsq
      go to (134,130,132,132,130),np
  130 do 131 i=mstart,mstop
         f(i,n) = f(i,n)-a1*f(i,np1)
  131 continue
      go to 134
  132 a1 = 2./deltht
      do 133 i=mstart,mstop
         f(i,np1) = f(i,np1)-a1*bdd(i)
  133 continue
  134 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 146,136,135
  135 ierror = 11
      go to 146
  136 w(id5+1) = .5*(w(id5+2)-dlrby2)
      go to (146,146,138,146,146,137),mbdcnd
  137 w(id5+1) = .5*w(id5+1)
  138 go to (140,146,146,139,146),np
  139 a2 = 2.
      go to 141
  140 a2 = 1.
  141 k = id5+munk
      w(k) = .5*(w(k-1)+dlrby2)
      s = 0.
      do 143 i=mstart,mstop
         s1 = 0.
         nsp1 = nstart+1
         nstm1 = nstop-1
         do 142 j=nsp1,nstm1
            s1 = s1+f(i,j)
  142    continue
         k = i+l
         s = s+(a2*s1+f(i,nstart)+f(i,nstop))*w(k)
  143 continue
      s2 = float(m)*a+(.75+float((m-1)*(m+1)))*dlrby2
      if (mbdcnd .eq. 3) s2 = s2+.25*dlrby2
      s1 = (2.+a2*float(nunk-2))*s2
      pertrb = s/s1
      do 145 i=mstart,mstop
         do 144 j=nstart,nstop
            f(i,j) = f(i,j)-pertrb
  144    continue
  145 continue
  146 continue
c
c     multiply i-th equation through by deltht**2 to put equation into
c     correct form for subroutine genbun.
c
      do 148 i=mstart,mstop
         k = i-mstart+1
         w(k) = w(k)*dlthsq
         j = id2+k
         w(j) = w(j)*dlthsq
         j = id3+k
         w(j) = w(j)*dlthsq
         do 147 j=nstart,nstop
            f(i,j) = f(i,j)*dlthsq
  147    continue
  148 continue
      w(1) = 0.
      w(id4) = 0.
c
c     call genbun to solve the system of equations.
c
      call genbun (nbdcnd,nunk,1,munk,w(1),w(id2+1),w(id3+1),idimf,
     &             f(mstart,nstart),ierr1,w(id4+1))
      w(1) = w(id4+1)+3.*float(munk)
      if (nbdcnd .ne. 0) go to 150
      do 149 i=mstart,mstop
         f(i,np1) = f(i,1)
  149 continue
  150 continue
      return
      end
      subroutine hwsplr (a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     &                   elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HWSPLR 2D standard grid polar coordinates five point scheme.
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c
c     hwsplr solves a finite difference approximation to the
c     helmholtz equation in polar coordinates:
c
c          (1/r)(d/dr)(r(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta)
c
c                                + lambda*u = f(r,theta).
c
c
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     a,b
c       the range of r, i.e., a .le. r .le. b.  a must be less than b
c       and a must be non-negative.
c
c     m
c       the number of panels into which the interval (a,b) is
c       subdivided.  hence, there will be m+1 grid points in the
c       r-direction given by r(i) = a+(i-1)dr, for i = 1,2,...,m+1,
c       where dr = (b-a)/m is the panel width. m must be greater than 3.
c
c     mbdcnd
c       indicates the type of boundary condition at r = a and r = b.
c
c       = 1  if the solution is specified at r = a and r = b.
c       = 2  if the solution is specified at r = a and the derivative of
c            the solution with respect to r is specified at r = b.
c       = 3  if the derivative of the solution with respect to r is
c            specified at r = a (see note below) and r = b.
c       = 4  if the derivative of the solution with respect to r is
c            specified at r = a (see note below) and the solution is
c            specified at r = b.
c       = 5  if the solution is unspecified at r = a = 0 and the
c            solution is specified at r = b.
c       = 6  if the solution is unspecified at r = a = 0 and the
c            derivative of the solution with respect to r is specified
c            at r = b.
c
c       note:  if a = 0, do not use mbdcnd = 3 or 4, but instead use
c              mbdcnd = 1,2,5, or 6  .
c
c     bda
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = a.
c       when mbdcnd = 3 or 4,
c
c            bda(j) = (d/dr)u(a,theta(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = b.
c       when mbdcnd = 2,3, or 6,
c
c            bdb(j) = (d/dr)u(b,theta(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdb is a dummy variable.
c
c     c,d
c       the range of theta, i.e., c .le. theta .le. d.  c must be less
c       than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       theta-direction given by theta(j) = c+(j-1)dtheta for
c       j = 1,2,...,n+1, where dtheta = (d-c)/n is the panel width.  n
c       must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at theta = c and
c       and theta = d.
c
c       = 0  if the solution is periodic in theta, i.e.,
c            u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at theta = c and theta = d
c            (see note below).
c       = 2  if the solution is specified at theta = c and the
c            derivative of the solution with respect to theta is
c            specified at theta = d (see note below).
c       = 4  if the derivative of the solution with respect to theta is
c            specified at theta = c and the solution is specified at
c            theta = d (see note below).
c
c       note:  when nbdcnd = 1,2, or 4, do not use mbdcnd = 5 or 6
c              (the former indicates that the solution is specified at
c              r = 0, the latter indicates the solution is unspecified
c              at r = 0).  use instead mbdcnd = 1 or 2  .
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = c.  when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dtheta)u(r(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = d.  when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dtheta)u(r(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .lt. 0, a solution may not exist.  however, hwsplr will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array that specifies the values of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(r(i),theta(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   -------------     -------------
c
c              1      u(a,theta(j))     u(b,theta(j))
c              2      u(a,theta(j))     f(b,theta(j))
c              3      f(a,theta(j))     f(b,theta(j))
c              4      f(a,theta(j))     u(b,theta(j))   j = 1,2,...,n+1
c              5      f(0,0)            u(b,theta(j))
c              6      f(0,0)            f(b,theta(j))
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   ---------         ---------
c
c              0      f(r(i),c)         f(r(i),c)
c              1      u(r(i),c)         u(r(i),d)
c              2      u(r(i),c)         f(r(i),d)   i = 1,2,...,m+1
c              3      f(r(i),c)         f(r(i),d)
c              4      f(r(i),c)         u(r(i),d)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note
c
c       if the table calls for both the solution u and the right side f
c       at  a corner then the solution must be specified.
c
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwsplr.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwsplr and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (r(i),theta(j)),
c       i = 1,2,...,m+1, j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic, derivative, or unspecified
c       boundary conditions is specified for a poisson equation
c       (lambda = 0), a solution may not exist.  pertrb is a constant,
c       calculated and subtracted from f, which ensures that a solution
c       exists.  hwsplr then computes this solution, which is a least
c       squares solution to the original approximation.  this solution
c       plus any constant is also a solution.  hence, the solution is
c       not unique.  pertrb should be small compared to the right side.
c       otherwise, a solution is obtained to an essentially different
c       problem.  this comparison should always be made to insure that a
c       meaningful solutin has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 11, a solution is not attempted.
c
c       =  0  no error.
c       =  1  a .lt. 0  .
c       =  2  a .ge. b.
c       =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6  .
c       =  4  c .ge. d.
c       =  5  n .le. 3
c       =  6  nbdcnd .lt. 0 or .gt. 4  .
c       =  7  a = 0, mbdcnd = 3 or 4  .
c       =  8  a .gt. 0, mbdcnd .ge. 5  .
c       =  9  mbdcnd .ge. 5, nbdcnd .ne. 0 and nbdcnd .ne. 3  .
c       = 10  idimf .lt. m+1  .
c       = 11  lambda .gt. 0  .
c       = 12  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwsplr, the user should test ierror after the call.
c
c     w
c       w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n+1),bdb(n+1),bdc(m+1),bdd(m+1),f(idimf,n+1),
c     arguments      w(see argument list)
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    hwsplr,genbun,poisd2,poisn2,poisp2,cosgen,merge,
c     required       trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized april 1, 1973
c                    revised january 1, 1976
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space          13430(octal) = 5912(decimal)  locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwsplr is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameters nbdcnd and mbdcnd.  some typical values
c                    are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for n and
c                    m as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    genbun which is the routine that
c                    solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32        1         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        1         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located  in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     &                w(*)
c
c     check for invalid parameters.
c
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 3) ierror = 5
      if (nbdcnd.le.-1 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4)) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (mbdcnd.ge.5 .and. nbdcnd.ne.0 .and. nbdcnd.ne.3) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if (m .le. 3) ierror = 12
      if (ierror .ne. 0) return
      mp1 = m+1
      deltar = (b-a)/float(m)
      dlrby2 = deltar/2.
      dlrsq = deltar**2
      np1 = n+1
      deltht = (d-c)/float(n)
      dlthsq = deltht**2
      np = nbdcnd+1
c
c     define range of indices i and j for unknowns u(i,j).
c
      mstart = 2
      mstop = mp1
      go to (101,105,102,103,104,105),mbdcnd
  101 mstop = m
      go to 105
  102 mstart = 1
      go to 105
  103 mstart = 1
  104 mstop = m
  105 munk = mstop-mstart+1
      nstart = 1
      nstop = n
      go to (109,106,107,108,109),np
  106 nstart = 2
      go to 109
  107 nstart = 2
  108 nstop = np1
  109 nunk = nstop-nstart+1
c
c     define a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      id5 = id4+munk
      id6 = id5+munk
      a1 = 2./dlrsq
      ij = 0
      if (mbdcnd.eq.3 .or. mbdcnd.eq.4) ij = 1
      do 110 i=1,munk
         r = a+float(i-ij)*deltar
         j = id5+i
         w(j) = r
         j = id6+i
         w(j) = 1./r**2
         w(i) = (r-dlrby2)/(r*dlrsq)
         j = id3+i
         w(j) = (r+dlrby2)/(r*dlrsq)
         j = id2+i
         w(j) = -a1+elmbda
  110 continue
      go to (114,111,112,113,114,111),mbdcnd
  111 w(id2) = a1
      go to 114
  112 w(id2) = a1
  113 w(id3+1) = a1
  114 continue
c
c     enter boundary data for r-boundaries.
c
      go to (115,115,117,117,119,119),mbdcnd
  115 a1 = w(1)
      do 116 j=nstart,nstop
         f(2,j) = f(2,j)-a1*f(1,j)
  116 continue
      go to 119
  117 a1 = 2.*deltar*w(1)
      do 118 j=nstart,nstop
         f(1,j) = f(1,j)+a1*bda(j)
  118 continue
  119 go to (120,122,122,120,120,122),mbdcnd
  120 a1 = w(id4)
      do 121 j=nstart,nstop
         f(m,j) = f(m,j)-a1*f(mp1,j)
  121 continue
      go to 124
  122 a1 = 2.*deltar*w(id4)
      do 123 j=nstart,nstop
         f(mp1,j) = f(mp1,j)-a1*bdb(j)
  123 continue
c
c     enter boundary data for theta-boundaries.
c
  124 a1 = 1./dlthsq
      l = id5-mstart+1
      lp = id6-mstart+1
      go to (134,125,125,127,127),np
  125 do 126 i=mstart,mstop
         j = i+lp
         f(i,2) = f(i,2)-a1*w(j)*f(i,1)
  126 continue
      go to 129
  127 a1 = 2./deltht
      do 128 i=mstart,mstop
         j = i+lp
         f(i,1) = f(i,1)+a1*w(j)*bdc(i)
  128 continue
  129 a1 = 1./dlthsq
      go to (134,130,132,132,130),np
  130 do 131 i=mstart,mstop
         j = i+lp
         f(i,n) = f(i,n)-a1*w(j)*f(i,np1)
  131 continue
      go to 134
  132 a1 = 2./deltht
      do 133 i=mstart,mstop
         j = i+lp
         f(i,np1) = f(i,np1)-a1*w(j)*bdd(i)
  133 continue
  134 continue
c
c     adjust right side of equation for unknown at pole when have
c     derivative specified boundary conditions.
c
      if (mbdcnd.ge.5 .and. nbdcnd.eq.3)
     &    f(1,1) = f(1,1)-(bdd(2)-bdc(2))*4./(float(n)*deltht*dlrsq)
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 144,136,135
  135 ierror = 11
      go to 144
  136 if (nbdcnd.ne.0 .and. nbdcnd.ne.3) go to 144
      s2 = 0.
      go to (144,144,137,144,144,138),mbdcnd
  137 w(id5+1) = .5*(w(id5+2)-dlrby2)
      s2 = .25*deltar
  138 a2 = 2.
      if (nbdcnd .eq. 0) a2 = 1.
      j = id5+munk
      w(j) = .5*(w(j-1)+dlrby2)
      s = 0.
      do 140 i=mstart,mstop
         s1 = 0.
         ij = nstart+1
         k = nstop-1
         do 139 j=ij,k
            s1 = s1+f(i,j)
  139    continue
         j = i+l
         s = s+(a2*s1+f(i,nstart)+f(i,nstop))*w(j)
  140 continue
      s2 = float(m)*a+deltar*(float((m-1)*(m+1))*.5+.25)+s2
      s1 = (2.+a2*float(nunk-2))*s2
      if (mbdcnd .eq. 3) go to 141
      s2 = float(n)*a2*deltar/8.
      s = s+f(1,1)*s2
      s1 = s1+s2
  141 continue
      pertrb = s/s1
      do 143 i=mstart,mstop
         do 142 j=nstart,nstop
            f(i,j) = f(i,j)-pertrb
  142    continue
  143 continue
  144 continue
c
c     multiply i-th equation through by (r(i)*deltht)**2.
c
      do 146 i=mstart,mstop
         k = i-mstart+1
         j = i+lp
         a1 = dlthsq/w(j)
         w(k) = a1*w(k)
         j = id2+k
         w(j) = a1*w(j)
         j = id3+k
         w(j) = a1*w(j)
         do 145 j=nstart,nstop
            f(i,j) = a1*f(i,j)
  145    continue
  146 continue
      w(1) = 0.
      w(id4) = 0.
c
c     call genbun to solve the system of equations.
c
      call genbun (nbdcnd,nunk,1,munk,w(1),w(id2+1),w(id3+1),idimf,
     &             f(mstart,nstart),ierr1,w(id4+1))
      iwstor = w(id4+1)+3.*float(munk)
      go to (157,157,157,157,148,147),mbdcnd
c
c     adjust the solution as necessary for the problems where a = 0.
c
  147 if (elmbda .ne. 0.) go to 148
      ypole = 0.
      go to 155
  148 continue
      j = id5+munk
      w(j) = w(id2)/w(id3)
      do 149 ip=3,munk
         i = munk-ip+2
         j = id5+i
         lp = id2+i
         k = id3+i
         w(j) = w(i)/(w(lp)-w(k)*w(j+1))
  149 continue
      w(id5+1) = -.5*dlthsq/(w(id2+1)-w(id3+1)*w(id5+2))
      do 150 i=2,munk
         j = id5+i
         w(j) = -w(j)*w(j-1)
  150 continue
      s = 0.
      do 151 j=nstart,nstop
         s = s+f(2,j)
  151 continue
      a2 = nunk
      if (nbdcnd .eq. 0) go to 152
      s = s-.5*(f(2,nstart)+f(2,nstop))
      a2 = a2-1.
  152 ypole = (.25*dlrsq*f(1,1)-s/a2)/(w(id5+1)-1.+elmbda*dlrsq*.25)
      do 154 i=mstart,mstop
         k = l+i
         do 153 j=nstart,nstop
            f(i,j) = f(i,j)+ypole*w(k)
  153    continue
  154 continue
  155 do 156 j=1,np1
         f(1,j) = ypole
  156 continue
  157 continue
      if (nbdcnd .ne. 0) go to 159
      do 158 i=mstart,mstop
         f(i,np1) = f(i,1)
  158 continue
  159 continue
      w(1) = iwstor
      return
      end
      subroutine hwsss1 (ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,
     &                   bdpf,elmbda,f,idimf,pertrb,am,bm,cm,sn,ss,
     &                   sint,d)

c*********************************************************************72
c
cc HWSSS1 is a utility routine for HWSSSP.
c
      dimension       f(idimf,*) ,bdts(*)    ,bdtf(*)    ,bdps(*)    ,
     &                bdpf(*)    ,am(*)      ,bm(*)      ,cm(*)      ,
     &                ss(*)      ,sn(*)      ,d(*)       ,sint(*)
c
      pi = pimach()
      tpi = pi+pi
      hpi = pi/2.
      mp1 = m+1
      np1 = n+1
      fn = n
      fm = m
      dth = (tf-ts)/fm
      hdth = dth/2.
      tdt = dth+dth
      dphi = (pf-ps)/fn
      tdp = dphi+dphi
      dphi2 = dphi*dphi
      edp2 = elmbda*dphi2
      dth2 = dth*dth
      cp = 4./(fn*dth2)
      wp = fn*sin(hdth)/4.
      do 102 i=1,mp1
         fim1 = i-1
         theta = fim1*dth+ts
         sint(i) = sin(theta)
         if (sint(i)) 101,102,101
  101    t1 = 1./(dth2*sint(i))
         am(i) = t1*sin(theta-hdth)
         cm(i) = t1*sin(theta+hdth)
         bm(i) = -am(i)-cm(i)+elmbda
  102 continue
      inp = 0
      isp = 0
c
c boundary condition at theta=ts
c
      mbr = mbdcnd+1
      go to (103,104,104,105,105,106,106,104,105,106),mbr
  103 its = 1
      go to 107
  104 at = am(2)
      its = 2
      go to 107
  105 at = am(1)
      its = 1
      cm(1) = am(1)+cm(1)
      go to 107
  106 at = am(2)
      inp = 1
      its = 2
c
c boundary condition theta=tf
c
  107 go to (108,109,110,110,109,109,110,111,111,111),mbr
  108 itf = m
      go to 112
  109 ct = cm(m)
      itf = m
      go to 112
  110 ct = cm(m+1)
      am(m+1) = am(m+1)+cm(m+1)
      itf = m+1
      go to 112
  111 itf = m
      isp = 1
      ct = cm(m)
c
c compute homogeneous solution with solution at pole equal to one
c
  112 itsp = its+1
      itfm = itf-1
      wts = sint(its+1)*am(its+1)/cm(its)
      wtf = sint(itf-1)*cm(itf-1)/am(itf)
      munk = itf-its+1
      if (isp) 116,116,113
  113 d(its) = cm(its)/bm(its)
      do 114 i=itsp,m
         d(i) = cm(i)/(bm(i)-am(i)*d(i-1))
  114 continue
      ss(m) = -d(m)
      iid = m-its
      do 115 ii=1,iid
         i = m-ii
         ss(i) = -d(i)*ss(i+1)
  115 continue
      ss(m+1) = 1.
  116 if (inp) 120,120,117
  117 sn(1) = 1.
      d(itf) = am(itf)/bm(itf)
      iid = itf-2
      do 118 ii=1,iid
         i = itf-ii
         d(i) = am(i)/(bm(i)-cm(i)*d(i+1))
  118 continue
      sn(2) = -d(2)
      do 119 i=3,itf
         sn(i) = -d(i)*sn(i-1)
  119 continue
c
c boundary conditions at phi=ps
c
  120 nbr = nbdcnd+1
      wps = 1.
      wpf = 1.
      go to (121,122,122,123,123),nbr
  121 jps = 1
      go to 124
  122 jps = 2
      go to 124
  123 jps = 1
      wps = .5
c
c boundary condition at phi=pf
c
  124 go to (125,126,127,127,126),nbr
  125 jpf = n
      go to 128
  126 jpf = n
      go to 128
  127 wpf = .5
      jpf = n+1
  128 jpsp = jps+1
      jpfm = jpf-1
      nunk = jpf-jps+1
      fjj = jpfm-jpsp+1
c
c scale coefficients for subroutine genbun
c
      do 129 i=its,itf
         cf = dphi2*sint(i)*sint(i)
         am(i) = cf*am(i)
         bm(i) = cf*bm(i)
         cm(i) = cf*cm(i)
  129 continue
      am(its) = 0.
      cm(itf) = 0.
      ising = 0
      go to (130,138,138,130,138,138,130,138,130,130),mbr
  130 go to (131,138,138,131,138),nbr
  131 if (elmbda) 138,132,132
  132 ising = 1
      sum = wts*wps+wts*wpf+wtf*wps+wtf*wpf
      if (inp) 134,134,133
  133 sum = sum+wp
  134 if (isp) 136,136,135
  135 sum = sum+wp
  136 sum1 = 0.
      do 137 i=itsp,itfm
         sum1 = sum1+sint(i)
  137 continue
      sum = sum+fjj*(sum1+wts+wtf)
      sum = sum+(wps+wpf)*sum1
      hne = sum
  138 go to (146,142,142,144,144,139,139,142,144,139),mbr
  139 if (nbdcnd-3) 146,140,146
  140 yhld = f(1,jps)-4./(fn*dphi*dth2)*(bdpf(2)-bdps(2))
      do 141 j=1,np1
         f(1,j) = yhld
  141 continue
      go to 146
  142 do 143 j=jps,jpf
         f(2,j) = f(2,j)-at*f(1,j)
  143 continue
      go to 146
  144 do 145 j=jps,jpf
         f(1,j) = f(1,j)+tdt*bdts(j)*at
  145 continue
  146 go to (154,150,152,152,150,150,152,147,147,147),mbr
  147 if (nbdcnd-3) 154,148,154
  148 yhld = f(m+1,jps)-4./(fn*dphi*dth2)*(bdpf(m)-bdps(m))
      do 149 j=1,np1
         f(m+1,j) = yhld
  149 continue
      go to 154
  150 do 151 j=jps,jpf
         f(m,j) = f(m,j)-ct*f(m+1,j)
  151 continue
      go to 154
  152 do 153 j=jps,jpf
         f(m+1,j) = f(m+1,j)-tdt*bdtf(j)*ct
  153 continue
  154 go to (159,155,155,157,157),nbr
  155 do 156 i=its,itf
         f(i,2) = f(i,2)-f(i,1)/(dphi2*sint(i)*sint(i))
  156 continue
      go to 159
  157 do 158 i=its,itf
         f(i,1) = f(i,1)+tdp*bdps(i)/(dphi2*sint(i)*sint(i))
  158 continue
  159 go to (164,160,162,162,160),nbr
  160 do 161 i=its,itf
         f(i,n) = f(i,n)-f(i,n+1)/(dphi2*sint(i)*sint(i))
  161 continue
      go to 164
  162 do 163 i=its,itf
         f(i,n+1) = f(i,n+1)-tdp*bdpf(i)/(dphi2*sint(i)*sint(i))
  163 continue
  164 continue
      pertrb = 0.
      if (ising) 165,176,165
  165 sum = wts*wps*f(its,jps)+wts*wpf*f(its,jpf)+wtf*wps*f(itf,jps)+
     &      wtf*wpf*f(itf,jpf)
      if (inp) 167,167,166
  166 sum = sum+wp*f(1,jps)
  167 if (isp) 169,169,168
  168 sum = sum+wp*f(m+1,jps)
  169 do 171 i=itsp,itfm
         sum1 = 0.
         do 170 j=jpsp,jpfm
            sum1 = sum1+f(i,j)
  170    continue
         sum = sum+sint(i)*sum1
  171 continue
      sum1 = 0.
      sum2 = 0.
      do 172 j=jpsp,jpfm
         sum1 = sum1+f(its,j)
         sum2 = sum2+f(itf,j)
  172 continue
      sum = sum+wts*sum1+wtf*sum2
      sum1 = 0.
      sum2 = 0.
      do 173 i=itsp,itfm
         sum1 = sum1+sint(i)*f(i,jps)
         sum2 = sum2+sint(i)*f(i,jpf)
  173 continue
      sum = sum+wps*sum1+wpf*sum2
      pertrb = sum/hne
      do 175 j=1,np1
         do 174 i=1,mp1
            f(i,j) = f(i,j)-pertrb
  174    continue
  175 continue
c
c scale right side for subroutine genbun
c
  176 do 178 i=its,itf
         cf = dphi2*sint(i)*sint(i)
         do 177 j=jps,jpf
            f(i,j) = cf*f(i,j)
  177    continue
  178 continue
      call genbun (nbdcnd,nunk,1,munk,am(its),bm(its),cm(its),idimf,
     &             f(its,jps),ierror,d)
      if (ising) 186,186,179
  179 if (inp) 183,183,180
  180 if (isp) 181,181,186
  181 do 182 j=1,np1
         f(1,j) = 0.
  182 continue
      go to 209
  183 if (isp) 186,186,184
  184 do 185 j=1,np1
         f(m+1,j) = 0.
  185 continue
      go to 209
  186 if (inp) 193,193,187
  187 sum = wps*f(its,jps)+wpf*f(its,jpf)
      do 188 j=jpsp,jpfm
         sum = sum+f(its,j)
  188 continue
      dfn = cp*sum
      dnn = cp*((wps+wpf+fjj)*(sn(2)-1.))+elmbda
      dsn = cp*(wps+wpf+fjj)*sn(m)
      if (isp) 189,189,194
  189 cnp = (f(1,1)-dfn)/dnn
      do 191 i=its,itf
         hld = cnp*sn(i)
         do 190 j=jps,jpf
            f(i,j) = f(i,j)+hld
  190    continue
  191 continue
      do 192 j=1,np1
         f(1,j) = cnp
  192 continue
      go to 209
  193 if (isp) 209,209,194
  194 sum = wps*f(itf,jps)+wpf*f(itf,jpf)
      do 195 j=jpsp,jpfm
         sum = sum+f(itf,j)
  195 continue
      dfs = cp*sum
      dss = cp*((wps+wpf+fjj)*(ss(m)-1.))+elmbda
      dns = cp*(wps+wpf+fjj)*ss(2)
      if (inp) 196,196,200
  196 csp = (f(m+1,1)-dfs)/dss
      do 198 i=its,itf
         hld = csp*ss(i)
         do 197 j=jps,jpf
            f(i,j) = f(i,j)+hld
  197    continue
  198 continue
      do 199 j=1,np1
         f(m+1,j) = csp
  199 continue
      go to 209
  200 rtn = f(1,1)-dfn
      rts = f(m+1,1)-dfs
      if (ising) 202,202,201
  201 csp = 0.
      cnp = rtn/dnn
      go to 205
  202 if (abs(dnn)-abs(dsn)) 204,204,203
  203 den = dss-dns*dsn/dnn
      rts = rts-rtn*dsn/dnn
      csp = rts/den
      cnp = (rtn-csp*dns)/dnn
      go to 205
  204 den = dns-dss*dnn/dsn
      rtn = rtn-rts*dnn/dsn
      csp = rtn/den
      cnp = (rts-dss*csp)/dsn
  205 do 207 i=its,itf
         hld = cnp*sn(i)+csp*ss(i)
         do 206 j=jps,jpf
            f(i,j) = f(i,j)+hld
  206    continue
  207 continue
      do 208 j=1,np1
         f(1,j) = cnp
         f(m+1,j) = csp
  208 continue
  209 if (nbdcnd) 212,210,212
  210 do 211 i=1,mp1
         f(i,jpf+1) = f(i,jps)
  211 continue
  212 return
      end
      subroutine hwsssp (ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,
     &                   bdpf,elmbda,f,idimf,pertrb,ierror,w)

c*********************************************************************72
c
cc HWSSSP 2D standard unit sphere coordinates five point scheme.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     hwsssp solves a finite difference approximation to the
c     helmholtz equation in spherical coordinates and on the surface of
c     the unit sphere (radius of 1):
c
c          (1/sin(theta))(d/dtheta)(sin(theta)(du/dtheta))
c
c             + (1/sin(theta)**2)(d/dphi)(du/dphi)
c
c             + lambda*u = f(theta,phi)
c
c     where theta is colatitude and phi is longitude.
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     ts,tf
c       the range of theta (colatitude), i.e., ts .le. theta .le. tf.
c       ts must be less than tf.  ts and tf are in radians.  a ts of
c       zero corresponds to the north pole and a tf of pi corresponds to
c       the south pole.
c
c     * * * * * * * * * * * * * * important * * * * * * * * * * * * * *
c
c     if tf is equal to pi then it must be computed using the statement
c     tf = pimach(). this insures that tf in the users program is
c     equal to pi in this program which permits several tests of the
c     input parameters that otherwise would not be possible.
c
c
c     m
c       the number of panels into which the interval (ts,tf) is
c       subdivided.  hence, there will be m+1 grid points in the
c       theta-direction given by theta(i) = (i-1)dtheta+ts for
c       i = 1,2,...,m+1, where dtheta = (tf-ts)/m is the panel width.
c       m must be greater than 5
c
c     mbdcnd
c       indicates the type of boundary condition at theta = ts and
c       theta = tf.
c
c       = 1  if the solution is specified at theta = ts and theta = tf.
c       = 2  if the solution is specified at theta = ts and the
c            derivative of the solution with respect to theta is
c            specified at theta = tf (see note 2 below).
c       = 3  if the derivative of the solution with respect to theta is
c            specified at theta = ts and theta = tf (see notes 1,2
c            below).
c       = 4  if the derivative of the solution with respect to theta is
c            specified at theta = ts (see note 1 below) and the
c            solution is specified at theta = tf.
c       = 5  if the solution is unspecified at theta = ts = 0 and the
c            solution is specified at theta = tf.
c       = 6  if the solution is unspecified at theta = ts = 0 and the
c            derivative of the solution with respect to theta is
c            specified at theta = tf (see note 2 below).
c       = 7  if the solution is specified at theta = ts and the
c            solution is unspecified at theta = tf = pi.
c       = 8  if the derivative of the solution with respect to theta is
c            specified at theta = ts (see note 1 below) and the
c            solution is unspecified at theta = tf = pi.
c       = 9  if the solution is unspecified at theta = ts = 0 and
c            theta = tf = pi.
c
c       notes:  1.  if ts = 0, do not use mbdcnd = 3,4, or 8, but
c                   instead use mbdcnd = 5,6, or 9  .
c               2.  if tf = pi, do not use mbdcnd = 2,3, or 6, but
c                   instead use mbdcnd = 7,8, or 9  .
c
c     bdts
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = ts.  when mbdcnd = 3,4, or 8,
c
c            bdts(j) = (d/dtheta)u(ts,phi(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdts is a dummy variable.
c
c     bdtf
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = tf.  when mbdcnd = 2,3, or 6,
c
c            bdtf(j) = (d/dtheta)u(tf,phi(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdtf is a dummy variable.
c
c     ps,pf
c       the range of phi (longitude), i.e., ps .le. phi .le. pf.  ps
c       must be less than pf.  ps and pf are in radians.  if ps = 0 and
c       pf = 2*pi, periodic boundary conditions are usually prescribed.
c
c     * * * * * * * * * * * * * * important * * * * * * * * * * * * * *
c
c     if pf is equal to 2*pi then it must be computed using the
c     statement pf = 2.*pimach(). this insures that pf in the users
c     program is equal to 2*pi in this program which permits tests of
c     the input parameters that otherwise would not be possible.
c
c
c     n
c       the number of panels into which the interval (ps,pf) is
c       subdivided.  hence, there will be n+1 grid points in the
c       phi-direction given by phi(j) = (j-1)dphi+ps  for
c       j = 1,2,...,n+1, where dphi = (pf-ps)/n is the panel width.
c       n must be greater than 4
c
c     nbdcnd
c       indicates the type of boundary condition at phi = ps and
c       phi = pf.
c
c       = 0  if the solution is periodic in phi, i.e.,
c            u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at phi = ps and phi = pf
c            (see note below).
c       = 2  if the solution is specified at phi = ps (see note below)
c            and the derivative of the solution with respect to phi is
c            specified at phi = pf.
c       = 3  if the derivative of the solution with respect to phi is
c            specified at phi = ps and phi = pf.
c       = 4  if the derivative of the solution with respect to phi is
c            specified at ps and the solution is specified at phi = pf
c            (see note below).
c
c       note:  nbdcnd = 1,2, or 4 cannot be used with
c              mbdcnd = 5,6,7,8, or 9 (the former indicates that the
c                       solution is specified at a pole, the latter
c                       indicates that the solution is unspecified).
c                       use instead
c              mbdcnd = 1 or 2  .
c
c     bdps
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to phi at
c       phi = ps.  when nbdcnd = 3 or 4,
c
c            bdps(i) = (d/dphi)u(theta(i),ps), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdps is a dummy variable.
c
c     bdpf
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to phi at
c       phi = pf.  when nbdcnd = 2 or 3,
c
c            bdpf(i) = (d/dphi)u(theta(i),pf), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdpf is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwsssp will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array that specifies the value of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m  and  j = 2,3,...,n
c
c            f(i,j) = f(theta(i),phi(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   ------------      ------------
c
c              1      u(ts,phi(j))      u(tf,phi(j))
c              2      u(ts,phi(j))      f(tf,phi(j))
c              3      f(ts,phi(j))      f(tf,phi(j))
c              4      f(ts,phi(j))      u(tf,phi(j))
c              5      f(0,ps)           u(tf,phi(j))   j = 1,2,...,n+1
c              6      f(0,ps)           f(tf,phi(j))
c              7      u(ts,phi(j))      f(pi,ps)
c              8      f(ts,phi(j))      f(pi,ps)
c              9      f(0,ps)           f(pi,ps)
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   --------------    --------------
c
c              0      f(theta(i),ps)    f(theta(i),ps)
c              1      u(theta(i),ps)    u(theta(i),pf)
c              2      u(theta(i),ps)    f(theta(i),pf)   i = 1,2,...,m+1
c              3      f(theta(i),ps)    f(theta(i),pf)
c              4      f(theta(i),ps)    u(theta(i),pf)
c
c       f must be dimensioned at least (m+1)*(n+1).
c       note
c
c       if the table calls for both the solution u and the right side f
c       at  a corner then the solution must be specified.
c
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwsssp.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space. w may require up to 4*(n+1)+(16+int(log2(n+1)))(m+1)
c       locations. the actual number of locations used is computed by
c       hwsssp and is output in location w(1). int( ) denotes the
c       fortran integer function.
c
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (theta(i),phi(j)),
c       i = 1,2,...,m+1,   j = 1,2,...,n+1  .
c
c     pertrb
c       if one specifies a combination of periodic, derivative or
c       unspecified boundary conditions for a poisson equation
c       (lambda = 0), a solution may not exist.  pertrb is a constant,
c       calculated and subtracted from f, which ensures that a solution
c       exists.  hwsssp then computes this solution, which is a least
c       squares solution to the original approximation.  this solution
c       is not unique and is unnormalized. the value of pertrb should
c       be small compared to the right side f. otherwise , a solution
c       is obtained to an essentially different problem. this comparison
c       should always be made to insure that a meaningful solution has
c       been obtained
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 8, a solution is not attempted.
c
c       = 0  no error
c       = 1  ts.lt.0 or tf.gt.pi
c       = 2  ts.ge.tf
c       = 3  mbdcnd.lt.1 or mbdcnd.gt.9
c       = 4  ps.lt.0 or ps.gt.pi+pi
c       = 5  ps.ge.pf
c       = 6  n.lt.5
c       = 7  m.lt.5
c       = 8  nbdcnd.lt.0 or nbdcnd.gt.4
c       = 9  elmbda.gt.0
c       = 10 idimf.lt.m+1
c       = 11 nbdcnd equals 1,2 or 4 and mbdcnd.ge.5
c       = 12 ts.eq.0 and mbdcnd equals 3,4 or 8
c       = 13 tf.eq.pi and mbdcnd equals 2,3 or 6
c       = 14 mbdcnd equals 5,6 or 9 and ts.ne.0
c       = 15 mbdcnd.ge.7 and tf.ne.pi
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwsssp, the user should test ierror after a call.
c
c     w
c       contains intermediate values that must not be destroyed if
c       hwsssp will be called again with intl = 1. w(1) contains the
c       required length of w .
c
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bdts(n+1),bdtf(n+1),bdps(m+1),bdpf(m+1),
c     arguments      f(idimf,n+1),w(see argument list)
c
c     latest         january 1978
c     revision
c
c
c     subprograms    hwsssp,hwsss1,genbun,poisd2,poisn2,poisp2,cosgen,me
c     required       trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     paul swarztrauber
c
c     language       fortran
c
c     history        version 1 - september 1973
c                    version 2 - april     1976
c                    version 3 - january   1978
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwsssp is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameters nbdcnd and mbdcnd.  some typical values
c                    are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for n and
c                    m as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    genbun which is the routine that
c                    solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32        0         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        0         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       sin,cos
c     resident
c     routines
c
c     reperences     p. n. swarztrauber,the direct solution of the
c                    discrete poisson equation on the surface of a
c                    sphere, SIAM Journal on Numerical Analysis,15(1974),pp 212-215
c
c                    swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      dimension       f(idimf,1) ,bdts(1)    ,bdtf(1)    ,bdps(1)    ,
     &                bdpf(1)    ,w(1)
      nbr = nbdcnd+1
      pi = pimach()
      tpi = 2.*pi
      ierror = 0
      if (ts.lt.0. .or. tf.gt.pi) ierror = 1
      if (ts .ge. tf) ierror = 2
      if (mbdcnd.lt.1 .or. mbdcnd.gt.9) ierror = 3
      if (ps.lt.0. .or. pf.gt.tpi) ierror = 4
      if (ps .ge. pf) ierror = 5
      if (n .lt. 5) ierror = 6
      if (m .lt. 5) ierror = 7
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 8
      if (elmbda .gt. 0.) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if ((nbdcnd.eq.1 .or. nbdcnd.eq.2 .or. nbdcnd.eq.4) .and.
     &    mbdcnd.ge.5) ierror = 11
      if (ts.eq.0. .and.
     &    (mbdcnd.eq.3 .or. mbdcnd.eq.4 .or. mbdcnd.eq.8)) ierror = 12
      if (tf.eq.pi .and.
     &    (mbdcnd.eq.2 .or. mbdcnd.eq.3 .or. mbdcnd.eq.6)) ierror = 13
      if ((mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9) .and.
     &    ts.ne.0.) ierror = 14
      if (mbdcnd.ge.7 .and. tf.ne.pi) ierror = 15
      if (ierror.ne.0 .and. ierror.ne.9) return
      call hwsss1 (ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,bdpf,
     &             elmbda,f,idimf,pertrb,w,w(m+2),w(2*m+3),w(3*m+4),
     &             w(4*m+5),w(5*m+6),w(6*m+7))
      w(1) = w(6*m+7)+float(6*(m+1))
      return
      end
      subroutine indxa (i,ir,idxa,na)

c*********************************************************************72
c
cc INDXA
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: cblkt
!$omp threadprivate (cblkt)
      na = 2**ir
      idxa = i-na+1
      if (i-nm) 102,102,101
  101 na = 0
  102 return
      end
      subroutine indxb (i,ir,idx,idp)

c*********************************************************************72
c
cc INDXB indexes the first root of the B(I,IR) polynomial.
c
c b(idx) is the location of the first root of the b(i,ir) polynomial
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: cblkt
!$omp threadprivate (cblkt)
      idp = 0
      if (ir) 107,101,103
  101 if (i-nm) 102,102,107
  102 idx = i
      idp = 1
      return
  103 izh = 2**ir
      id = i-izh-izh
      idx = id+id+(ir-1)*ik+ir+(ik-i)/izh+4
      ipl = izh-1
      idp = izh+izh-1
      if (i-ipl-nm) 105,105,104
  104 idp = 0
      return
  105 if (i+ipl-nm) 107,107,106
  106 idp = nm+ipl-i+1
  107 return
      end
      subroutine indxc (i,ir,idxc,nc)

c*********************************************************************72
c
cc INDXC
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
      save :: cblkt
!$omp threadprivate (cblkt)
      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
  101 nc = 0
  102 return
      end
      subroutine inxca (i,ir,idxa,na)

c*********************************************************************72
c
cc INXCA
c
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: ccblk
!$omp threadprivate (ccblk)
      na = 2**ir
      idxa = i-na+1
      if (i-nm) 102,102,101
  101 na = 0
  102 return
      end
      subroutine inxcb (i,ir,idx,idp)

c*********************************************************************72
c
cc INXCB
c
c b(idx) is the location of the first root of the b(i,ir) polynomial
c
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: ccblk
!$omp threadprivate (ccblk)
      idp = 0
      if (ir) 107,101,103
  101 if (i-nm) 102,102,107
  102 idx = i
      idp = 1
      return
  103 izh = 2**ir
      id = i-izh-izh
      idx = id+id+(ir-1)*ik+ir+(ik-i)/izh+4
      ipl = izh-1
      idp = izh+izh-1
      if (i-ipl-nm) 105,105,104
  104 idp = 0
      return
  105 if (i+ipl-nm) 107,107,106
  106 idp = nm+ipl-i+1
  107 return
      end
      subroutine inxcc (i,ir,idxc,nc)

c*********************************************************************72
c
cc INXCC
c
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: ccblk
!$omp threadprivate (ccblk)
      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
  101 nc = 0
  102 return
      end
      function j4save (iwhich, ivalue, iset)

c*********************************************************************72
c
cc J4SAVE sets or gets variables needed by the error handler.
c
c***BEGIN PROLOGUE  J4SAVE
c***SUBSIDIARY
c***PURPOSE  Save or recall global variables needed by error
c            handling routines.
c***LIBRARY   SLATEC (XERROR)
c***TYPE      INTEGER (J4SAVE-I)
c***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        J4SAVE saves and recalls several global variables needed
c        by the library error handling routines.
c
c     Description of Parameters
c      --Input--
c        IWHICH - Index of item desired.
c                = 1 Refers to current error number.
c                = 2 Refers to current error control flag.
c                = 3 Refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                = 4 Refers to the maximum number of times any
c                     message is to be printed (as set by XERMAX).
c                = 5 Refers to the total number of units to which
c                     each error message is to be written.
c                = 6 Refers to the 2nd unit for error messages
c                = 7 Refers to the 3rd unit for error messages
c                = 8 Refers to the 4th unit for error messages
c                = 9 Refers to the 5th unit for error messages
c        IVALUE - The value to be set for the IWHICH-th parameter,
c                 if ISET is .TRUE. .
c        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
c                 given the value, IVALUE.  If ISET=.FALSE., the
c                 IWHICH-th parameter will be unchanged, and IVALUE
c                 is a dummy parameter.
c      --Output--
c        The (old) value of the IWHICH-th parameter will be returned
c        in the function value, J4SAVE.
c
c***SEE ALSO  XERMSG
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  (NONE)
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900205  Minor modifications to prologue.  (WRB)
c   900402  Added TYPE section.  (WRB)
c   910411  Added KEYWORDS section.  (WRB)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  J4SAVE
      logical iset
      integer iparam(9)
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***FIRST EXECUTABLE STATEMENT  J4SAVE
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end
      subroutine merge (tcos,i1,m1,i2,m2,i3)

c*********************************************************************72
c
cc MERGE merges two ascending strings of numbers in the array TCOS.
c
c     the first string is of length m1 and starts at
c     tcos(i1+1).  the second string is of length m2 and starts at
c     tcos(i2+1).  the merged string goes into tcos(i3+1).
c
      dimension tcos(*)

      j1 = 1
      j2 = 1
      j = i3
      if (m1 .eq. 0) go to 107
      if (m2 .eq. 0) go to 104
  101 j = j+1
      l = j1+i1
      x = tcos(l)
      l = j2+i2
      y = tcos(l)
      if (x-y) 102,102,103
  102 tcos(j) = x
      j1 = j1+1
      if (j1 .gt. m1) go to 106
      go to 101
  103 tcos(j) = y
      j2 = j2+1
      if (j2 .le. m2) go to 101
      if (j1 .gt. m1) go to 109
  104 k = j-j1+1
      do 105 j=j1,m1
         m = k+j
         l = j+i1
         tcos(m) = tcos(l)
  105 continue
      go to 109
  106 continue
      if (j2 .gt. m2) go to 109
  107 k = j-j2+1
      do 108 j=j2,m2
         m = k+j
         l = j+i2
         tcos(m) = tcos(l)
  108 continue
  109 continue
      return
      end
      subroutine minso4 (usol,idmn,zn,zm,pertb)

c*********************************************************************72
c
cc MINSO4 orthogonalizes the array usol with respect to
c     the constant array in a weighted least squares norm
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: spl4
!$omp threadprivate (spl4)
      dimension       usol(idmn,*)           ,zn(*)      ,zm(*)
c
c     entry at minso4 occurrs when the final solution is
c     to be minimized with respect to the weighted
c     least squares norm
c
      istr = 1
      ifnl = k
      jstr = 1
      jfnl = l
c
c     compute weighted inner products
c
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
c
c     set perturbation parameter
c
      pertrb = ute/ete
c
c     subtract off constant pertrb
c
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end
      subroutine minsol (usol,idmn,zn,zm,pertb)

c*********************************************************************72
c
cc MINSOL orthogonalizes the array usol with respect to
c     the constant array in a weighted least squares norm
c
      common /splp/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: splp
!$omp threadprivate (splp)
      dimension       usol(idmn,*)           ,zn(*)      ,zm(*)
c
c     entry at minsol occurrs when the final solution is
c     to be minimized with respect to the weighted
c     least squares norm
c
      istr = 1
      ifnl = k
      jstr = 1
      jfnl = l
c
c     compute weighted inner products
c
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
c
c     set perturbation parameter
c
      pertrb = ute/ete
c
c     subtract off constant pertrb
c
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end
      subroutine ortho4 (usol,idmn,zn,zm,pertrb)

c*********************************************************************72
c
cc ORTHO4 orthogonalizes the array usol with respect to
c     the constant array in a weighted least squares norm
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: spl4
!$omp threadprivate (spl4)
      dimension       usol(idmn,*)           ,zn(*)      ,zm(*)
      istr = is
      ifnl = ms
      jstr = js
      jfnl = ns
c
c     compute weighted inner products
c
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
c
c     set perturbation parameter
c
      pertrb = ute/ete
c
c     subtract off constant pertrb
c
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end
      subroutine orthog (usol,idmn,zn,zm,pertrb)

c*********************************************************************72
c
cc ORTHOG orthogonalizes the array usol with respect to
c     the constant array in a weighted least squares norm
c
      common /splp/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: splp
!$omp threadprivate (splp)
      dimension       usol(idmn,*)           ,zn(*)      ,zm(*)
      istr = is
      ifnl = ms
      jstr = js
      jfnl = ns
c
c     compute weighted inner products
c
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
c
c     set perturbation parameter
c
      pertrb = ute/ete
c
c     subtract off constant pertrb
c
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end
      function pgsf (x,iz,c,a,bh)

c*********************************************************************72
c
cc PGSF
c
      dimension       a(*)       ,c(*)       ,bh(*)
      fsg = 1.
      hsg = 1.
      do 101 j=1,iz
         dd = 1./(x-bh(j))
         fsg = fsg*a(j)*dd
         hsg = hsg*c(j)*dd
  101 continue
      if (mod(iz,2)) 103,102,103
  102 pgsf = 1.-fsg-hsg
      return
  103 pgsf = 1.+fsg+hsg
      return
      end
      function pimach ()

c*********************************************************************72
c
cc PIMACH supplies the value of the constant pi correct to
c     machine precision where
c
c     pi=3.1415926535897932384626433832795028841971693993751058209749446
c
      real pimach

      pimach = 3.14159265358979
      return
      end
      subroutine pois3d (lperod,l,c1,mperod,m,c2,nperod,n,a,b,c,ldimf,
     &                   mdimf,f,ierror,w)

c*********************************************************************72
c
cc POIS3D solves a special set of linear equations.
c
c     pois3d solves the linear system of equations
c
c       c1*(x(i-1,j,k)-2.*x(i,j,k)+x(i+1,j,k))
c     + c2*(x(i,j-1,k)-2.*x(i,j,k)+x(i,j+1,k))
c     + a(k)*x(i,j,k-1)+b(k)*x(i,j,k)+c(k)*x(i,j,k+1) = f(i,j,k)
c
c     for  i=1,2,...,l , j=1,2,...,m , and k=1,2,...,n .
c
c     the indices k-1 and k+1 are evaluated modulo n, i.e.
c     x(i,j,0) = x(i,j,n) and x(i,j,n+1) = x(i,j,1). the unknowns
c     x(0,j,k), x(l+1,j,k), x(i,0,k), and x(i,m+1,k) are assumed to take
c     on certain prescribed values described below.
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c     lperod   indicates the values that x(0,j,k) and x(l+1,j,k) are
c              assumed to have.
c
c              = 0  if x(0,j,k) = x(l,j,k) and x(l+1,j,k) = x(1,j,k).
c              = 1  if x(0,j,k) = x(l+1,j,k) = 0.
c              = 2  if x(0,j,k) = 0  and x(l+1,j,k) = x(l-1,j,k).
c              = 3  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = x(l-1,j,k).
c              = 4  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = 0.
c
c     l        the number of unknowns in the i-direction. l must be at
c              least 3.
c
c     c1       the real constant that appears in the above equation.
c
c     mperod   indicates the values that x(i,0,k) and x(i,m+1,k) are
c              assumed to have.
c
c              = 0  if x(i,0,k) = x(i,m,k) and x(i,m+1,k) = x(i,1,k).
c              = 1  if x(i,0,k) = x(i,m+1,k) = 0.
c              = 2  if x(i,0,k) = 0 and x(i,m+1,k) = x(i,m-1,k).
c              = 3  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = x(i,m-1,k).
c              = 4  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = 0.
c
c     m        the number of unknowns in the j-direction. m must be at
c              least 3.
c
c     c2       the real constant which appears in the above equation.
c
c     nperod   = 0  if a(1) and c(n) are not zero.
c              = 1  if a(1) = c(n) = 0.
c
c     n        the number of unknowns in the k-direction. n must be at
c              least 3.
c
c
c     a,b,c    one-dimensional arrays of length n that specify the
c              coefficients in the linear equations given above.
c
c              if nperod = 0 the array elements must not depend upon the
c              index k, but must be constant.  specifically,the
c              routine checks the following condition
c
c                          a(k) = c(1)
c                          c(k) = c(1)
c                          b(k) = b(1)
c
c                  for k=1,2,...,n.
c
c     ldimf    the row (or first) dimension of the three-dimensional
c              array f as it appears in the program calling pois3d.
c              this parameter is used to specify the variable dimension
c              of f.  ldimf must be at least l.
c
c     mdimf    the column (or second) dimension of the three-dimensional
c              array f as it appears in the program calling pois3d.
c              this parameter is used to specify the variable dimension
c              of f.  mdimf must be at least m.
c
c     f        a three-dimensional array that specifies the values of
c              the right side of the linear system of equations given
c              above.  f must be dimensioned at least l x m x n.
c
c     w        a one-dimensional array that must be provided by the
c              user for work space.  the length of w must be at least
c              30 + l + m + 2*n + max(l,m,n) +
c              7*(int((l+1)/2) + int((m+1)/2)).
c
c
c            * * * * * *   on output   * * * * * *
c
c     f        contains the solution x.
c
c     ierror   an error flag that indicates invalid input parameters.
c              except for number zero, a solution is not attempted.
c              = 0  no error
c              = 1  if lperod .lt. 0 or .gt. 4
c              = 2  if l .lt. 3
c              = 3  if mperod .lt. 0 or .gt. 4
c              = 4  if m .lt. 3
c              = 5  if nperod .lt. 0 or .gt. 1
c              = 6  if n .lt. 3
c              = 7  if ldimf .lt. l
c              = 8  if mdimf .lt. m
c              = 9  if a(k) .ne. c(1) or c(k) .ne. c(1) or b(i) .ne.b(1)
c                      for some k=1,2,...,n.
c              = 10 if nperod = 1 and a(1) .ne. 0 or c(n) .ne. 0
c
c              since this is the only means of indicating a possibly
c              incorrect call to pois3d, the user should test ierror
c              after the call.
c
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(n),b(n),c(n),f(ldimf,mdimf,n),
c     arguments      w(see argument list)
c
c     latest         december 1, 1978
c     revision
c
c     subprograms    pois3d,pos3d1,trid,rffti,rfftf,rfftf1,rfftb,
c     required       rfftb1,costi,cost,sinti,sint,cosqi,cosqf,cosqf1
c                    cosqb,cosqb1,sinqi,sinqf,sinqb,cffti,cffti1,
c                    cfftb,cfftb1,passb2,passb3,passb4,passb,cfftf,
c                    cfftf1,passf1,passf2,passf3,passf4,passf,pimach,
c
c     special        none
c     conditions
c
c     common         value
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in july,1977
c
c     algorithm      this subroutine solves three-dimensional block
c                    tridiagonal linear systems arising from finite
c                    difference approximations to three-dimensional
c                    poisson equations using the fourier transform
c                    package sclrfftpak written by paul swarztrauber.
c
c     space          6561(decimal) = 14641(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine pois3d is roughly proportional
c                    to l*m*n*(log2(l)+log2(m)+5), but also depends on
c                    input parameters lperod and mperod.  some typical
c                    values are listed in the table below when nperod=0.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(k) = c(k) = -0.5*b(k) = 1,       k=1,2,...,n
c
c                    and, when nperod = 1
c
c                       a(1) = c(n) = 0
c                       a(n) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine pois was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                    e = max(abs(z(i,j,k)-x(i,j,k)))/max(abs(x(i,j,k)))
c
c                    where the two maxima are taken over i=1,2,...,l,
c                    j=1,2,...,m and k=1,2,...,n, was computed.  the
c                    value of e is given in the table below for some
c                    typical values of l,m and n.
c
c
c                       l(=m=n)   lperod    mperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         16        0         0         272     1.e-13
c                         15        1         1         287     4.e-13
c                         17        3         3         338     2.e-13
c                         32        0         0        1755     2.e-13
c                         31        1         1        1894     2.e-12
c                         33        3         3        2042     7.e-13
c
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos,sin,atan
c     resident
c     routines
c
c     reference      none
c
      dimension       a(*)       ,b(*)       ,c(*)       ,
     &                f(ldimf,mdimf,*)       ,w(*)       ,save(6)
      lp = lperod+1
      mp = mperod+1
      np = nperod+1
c
c     check for invalid input.
c
      ierror = 0
      if (lp.lt.1 .or. lp.gt.5) ierror = 1
      if (l .lt. 3) ierror = 2
      if (mp.lt.1 .or. mp.gt.5) ierror = 3
      if (m .lt. 3) ierror = 4
      if (np.lt.1 .or. np.gt.2) ierror = 5
      if (n .lt. 3) ierror = 6
      if (ldimf .lt. l) ierror = 7
      if (mdimf .lt. m) ierror = 8
      if (np .ne. 1) go to 103
      do 101 k=1,n
         if (a(k) .ne. c(1)) go to 102
         if (c(k) .ne. c(1)) go to 102
         if (b(k) .ne. b(1)) go to 102
  101 continue
      go to 104
  102 ierror = 9
  103 if (nperod.eq.1 .and. (a(1).ne.0. .or. c(n).ne.0.)) ierror = 10
  104 if (ierror .ne. 0) go to 122
      iwyrt = l+1
      iwt = iwyrt+m
      iwd = iwt+max0(l,m,n)+1
      iwbb = iwd+n
      iwx = iwbb+n
      iwy = iwx+7*((l+1)/2)+15
      go to (105,114),np
c
c     reorder unknowns when nperod = 0.
c
  105 nh = (n+1)/2
      nhm1 = nh-1
      nodd = 1
      if (2*nh .eq. n) nodd = 2
      do 111 i=1,l
         do 110 j=1,m
            do 106 k=1,nhm1
               nhpk = nh+k
               nhmk = nh-k
               w(k) = f(i,j,nhmk)-f(i,j,nhpk)
               w(nhpk) = f(i,j,nhmk)+f(i,j,nhpk)
  106       continue
            w(nh) = 2.*f(i,j,nh)
            go to (108,107),nodd
  107       w(n) = 2.*f(i,j,n)
  108       do 109 k=1,n
               f(i,j,k) = w(k)
  109       continue
  110    continue
  111 continue
      save(1) = c(nhm1)
      save(2) = a(nh)
      save(3) = c(nh)
      save(4) = b(nhm1)
      save(5) = b(n)
      save(6) = a(n)
      c(nhm1) = 0.
      a(nh) = 0.
      c(nh) = 2.*c(nh)
      go to (112,113),nodd
  112 b(nhm1) = b(nhm1)-a(nh-1)
      b(n) = b(n)+a(n)
      go to 114
  113 a(n) = c(nh)
  114 continue
      call pos3d1 (lp,l,mp,m,n,a,b,c,ldimf,mdimf,f,w,w(iwyrt),w(iwt),
     &             w(iwd),w(iwx),w(iwy),c1,c2,w(iwbb))
      go to (115,122),np
  115 do 121 i=1,l
         do 120 j=1,m
            do 116 k=1,nhm1
               nhmk = nh-k
               nhpk = nh+k
               w(nhmk) = .5*(f(i,j,nhpk)+f(i,j,k))
               w(nhpk) = .5*(f(i,j,nhpk)-f(i,j,k))
  116       continue
            w(nh) = .5*f(i,j,nh)
            go to (118,117),nodd
  117       w(n) = .5*f(i,j,n)
  118       do 119 k=1,n
               f(i,j,k) = w(k)
  119       continue
  120    continue
  121 continue
      c(nhm1) = save(1)
      a(nh) = save(2)
      c(nh) = save(3)
      b(nhm1) = save(4)
      b(n) = save(5)
      a(n) = save(6)
  122 continue
      return
      end
      subroutine poisd2 (mr,nr,istag,ba,bb,bc,q,idimq,b,w,d,tcos,p)

c*********************************************************************72
c
cc POISD2 solves Poisson's equation for Dirichlet boundary conditions.
c
c     istag = 1 if the last diagonal block is the matrix a.
c     istag = 2 if the last diagonal block is the matrix a+i.
c
      dimension       q(idimq,*) ,ba(*)      ,bb(*)      ,bc(*)      ,
     &                tcos(*)    ,b(*)       ,d(*)       ,w(*)       ,
     &                p(*)
      m = mr
      n = nr
      jsh = 0
      fi = 1./float(istag)
      ip = -m
      ipstor = 0
      go to (101,102),istag
  101 kr = 0
      irreg = 1
      if (n .gt. 1) go to 106
      tcos(1) = 0.
      go to 103
  102 kr = 1
      jstsav = 1
      irreg = 2
      if (n .gt. 1) go to 106
      tcos(1) = -1.
  103 do 104 i=1,m
         b(i) = q(i,1)
  104 continue
      call trix (1,0,m,ba,bb,bc,b,tcos,d,w)
      do 105 i=1,m
         q(i,1) = b(i)
  105 continue
      go to 183
  106 lr = 0
      do 107 i=1,m
         p(i) = 0.
  107 continue
      nun = n
      jst = 1
      jsp = n
c
c     irreg = 1 when no irregularities have occurred, otherwise it is 2.
c
  108 l = 2*jst
      nodd = 2-2*((nun+1)/2)+nun
c
c     nodd = 1 when nun is odd, otherwise it is 2.
c
      go to (110,109),nodd
  109 jsp = jsp-l
      go to 111
  110 jsp = jsp-jst
      if (irreg .ne. 1) jsp = jsp-l
  111 continue
c
c     regular reduction
c
      call cosgen (jst,1,0.5,0.0,tcos)
      if (l .gt. jsp) go to 118
      do 117 j=l,jsp,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         jm3 = jm2-jsh
         jp3 = jp2+jsh
         if (jst .ne. 1) go to 113
         do 112 i=1,m
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  112    continue
         go to 115
  113    do 114 i=1,m
            t = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = t+q(i,j)-q(i,jm3)-q(i,jp3)
            q(i,j) = t
  114    continue
  115    continue
         call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
         do 116 i=1,m
            q(i,j) = q(i,j)+b(i)
  116    continue
  117 continue
c
c     reduction for last unknown
c
  118 go to (119,136),nodd
  119 go to (152,120),irreg
c
c     odd number of unknowns
c
  120 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (123,121),istag
  121 continue
      if (jst .ne. 1) go to 123
      do 122 i=1,m
         b(i) = q(i,j)
         q(i,j) = 0.
  122 continue
      go to 130
  123 go to (124,126),noddpr
  124 do 125 i=1,m
         ip1 = ip+i
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+p(ip1)+q(i,j)
  125 continue
      go to 128
  126 do 127 i=1,m
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+q(i,jp2)-q(i,jp1)+q(i,j)
  127 continue
  128 do 129 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  129 continue
  130 call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
      ip = ip+m
      ipstor = max0(ipstor,ip+m)
      do 131 i=1,m
         ip1 = ip+i
         p(ip1) = q(i,j)+b(i)
         b(i) = q(i,jp2)+p(ip1)
  131 continue
      if (lr .ne. 0) go to 133
      do 132 i=1,jst
         krpi = kr+i
         tcos(krpi) = tcos(i)
  132 continue
      go to 134
  133 continue
      call cosgen (lr,jstsav,0.,fi,tcos(jst+1))
      call merge (tcos,0,jst,jst,lr,kr)
  134 continue
      call cosgen (kr,jstsav,0.0,fi,tcos)
      call trix (kr,kr,m,ba,bb,bc,b,tcos,d,w)
      do 135 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+b(i)+p(ip1)
  135 continue
      lr = kr
      kr = kr+l
      go to 152
c
c     even number of unknowns
c
  136 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (137,138),irreg
  137 continue
      jstsav = jst
      ideg = jst
      kr = l
      go to 139
  138 call cosgen (kr,jstsav,0.0,fi,tcos)
      call cosgen (lr,jstsav,0.0,fi,tcos(kr+1))
      ideg = kr
      kr = kr+jst
  139 if (jst .ne. 1) go to 141
      irreg = 2
      do 140 i=1,m
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  140 continue
      go to 150
  141 do 142 i=1,m
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  142 continue
      go to (143,145),irreg
  143 do 144 i=1,m
         q(i,j) = q(i,jm2)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
  144 continue
      irreg = 2
      go to 150
  145 continue
      go to (146,148),noddpr
  146 do 147 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+p(ip1)
  147 continue
      ip = ip-m
      go to 150
  148 do 149 i=1,m
         q(i,j) = q(i,jm2)+q(i,j)-q(i,jm1)
  149 continue
  150 call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      do 151 i=1,m
         q(i,j) = q(i,j)+b(i)
  151 continue
  152 nun = nun/2
      noddpr = nodd
      jsh = jst
      jst = 2*jst
      if (nun .ge. 2) go to 108
c
c     start solution.
c
      j = jsp
      do 153 i=1,m
         b(i) = q(i,j)
  153 continue
      go to (154,155),irreg
  154 continue
      call cosgen (jst,1,0.5,0.0,tcos)
      ideg = jst
      go to 156
  155 kr = lr+jst
      call cosgen (kr,jstsav,0.0,fi,tcos)
      call cosgen (lr,jstsav,0.0,fi,tcos(kr+1))
      ideg = kr
  156 continue
      call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      jm1 = j-jsh
      jp1 = j+jsh
      go to (157,159),irreg
  157 do 158 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  158 continue
      go to 164
  159 go to (160,162),noddpr
  160 do 161 i=1,m
         ip1 = ip+i
         q(i,j) = p(ip1)+b(i)
  161 continue
      ip = ip-m
      go to 164
  162 do 163 i=1,m
         q(i,j) = q(i,j)-q(i,jm1)+b(i)
  163 continue
  164 continue
c
c     start back substitution.
c
      jst = jst/2
      jsh = jst/2
      nun = 2*nun
      if (nun .gt. n) go to 183
      do 182 j=jst,n,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         if (j .gt. jst) go to 166
         do 165 i=1,m
            b(i) = q(i,j)+q(i,jp2)
  165    continue
         go to 170
  166    if (jp2 .le. n) go to 168
         do 167 i=1,m
            b(i) = q(i,j)+q(i,jm2)
  167    continue
         if (jst .lt. jstsav) irreg = 1
         go to (170,171),irreg
  168    do 169 i=1,m
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  169    continue
  170    continue
         call cosgen (jst,1,0.5,0.0,tcos)
         ideg = jst
         jdeg = 0
         go to 172
  171    if (j+l .gt. n) lr = lr-jst
         kr = jst+lr
         call cosgen (kr,jstsav,0.0,fi,tcos)
         call cosgen (lr,jstsav,0.0,fi,tcos(kr+1))
         ideg = kr
         jdeg = lr
  172    continue
         call trix (ideg,jdeg,m,ba,bb,bc,b,tcos,d,w)
         if (jst .gt. 1) go to 174
         do 173 i=1,m
            q(i,j) = b(i)
  173    continue
         go to 182
  174    if (jp2 .gt. n) go to 177
  175    do 176 i=1,m
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  176    continue
         go to 182
  177    go to (175,178),irreg
  178    if (j+jsh .gt. n) go to 180
         do 179 i=1,m
            ip1 = ip+i
            q(i,j) = b(i)+p(ip1)
  179    continue
         ip = ip-m
         go to 182
  180    do 181 i=1,m
            q(i,j) = b(i)+q(i,j)-q(i,jm1)
  181    continue
  182 continue
      l = l/2
      go to 164
  183 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine poisn2 (m,n,istag,mixbnd,a,bb,c,q,idimq,b,b2,b3,w,w2,
     &                   w3,d,tcos,p)

c*********************************************************************72
c
cc POISN2 solves Poisson's equation with Neumann boundary conditions.
c
c     istag = 1 if the last diagonal block is a.
c     istag = 2 if the last diagonal block is a-i.
c     mixbnd = 1 if have neumann boundary conditions at both boundaries.
c     mixbnd = 2 if have neumann boundary conditions at bottom and
c     dirichlet condition at top.  (for this case, must have istag = 1.)
c
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) ,
     &                b(*)       ,b2(*)      ,b3(*)      ,w(*)       ,
     &                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    ,
     &                k(4)       ,p(*)
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
      fistag = 3-istag
      fnum = 1./float(istag)
      fden = 0.5*float(istag-1)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      go to (101,103),istag
  101 continue
      do 102 i=1,mr
         q(i,n) = .5*q(i,n)
  102 continue
      go to (103,104),mixbnd
  103 if (n .le. 3) go to 155
  104 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      go to (105,106),mixbnd
  105 jstart = 1
      go to 107
  106 jstart = jr
      nrod = 1-nrod
  107 continue
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      call cosgen (i2r,1,0.5,0.0,tcos)
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 108
      j = jr
      go to 116
  108 continue
c
c     regular reduction.
c
      do 115 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 109
         jm1 = jp1
         jm2 = jp2
         jm3 = jp3
  109    continue
         if (i2r .ne. 1) go to 111
         if (j .eq. 1) jm2 = jp2
         do 110 i=1,mr
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  110    continue
         go to 113
  111    continue
         do 112 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  112    continue
  113    continue
         call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 114 i=1,mr
            q(i,j) = q(i,j)+b(i)
  114    continue
c
c     end of reduction for regular unknowns.
c
  115 continue
c
c     begin special reduction for last unknown.
c
      j = jstop+jr
  116 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 128
c
c     odd number of unknowns
c
      if (i2r .ne. 1) go to 118
      do 117 i=1,mr
         b(i) = fistag*q(i,j)
         q(i,j) = q(i,jm2)
  117 continue
      go to 126
  118 do 119 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  119 continue
      if (nrodpr .ne. 0) go to 121
      do 120 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  120 continue
      ip = ip-mr
      go to 123
  121 continue
      do 122 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  122 continue
  123 if (lr .eq. 0) go to 124
      call cosgen (lr,1,0.5,fden,tcos(kr+1))
      go to 126
  124 continue
      do 125 i=1,mr
         b(i) = fistag*b(i)
  125 continue
  126 continue
      call cosgen (kr,1,0.5,fden,tcos)
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 127 i=1,mr
         q(i,j) = q(i,j)+b(i)
  127 continue
      kr = kr+i2r
      go to 151
  128 continue
c
c     even number of unknowns
c
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 135
      do 129 i=1,mr
         b(i) = q(i,j)
  129 continue
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      go to (133,130),istag
  130 do 131 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  131 continue
      tcos(1) = 1.
      tcos(2) = 0.
      call trix (1,1,mr,a,bb,c,b,tcos,d,w)
      do 132 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  132 continue
      go to 150
  133 continue
      do 134 i=1,mr
         p(i) = b(i)
         q(i,j) = q(i,jm2)+2.*q(i,jp2)+3.*b(i)
  134 continue
      go to 150
  135 continue
      do 136 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  136 continue
      if (nrodpr .ne. 0) go to 138
      do 137 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  137 continue
      go to 140
  138 continue
      do 139 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  139 continue
  140 continue
      call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max0(ipstor,ip+mr)
      do 141 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  141 continue
      if (lr .eq. 0) go to 142
      call cosgen (lr,1,0.5,fden,tcos(i2r+1))
      call merge (tcos,0,i2r,i2r,lr,kr)
      go to 144
  142 do 143 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  143 continue
  144 call cosgen (kr,1,0.5,fden,tcos)
      if (lr .ne. 0) go to 145
      go to (146,145),istag
  145 continue
      call trix (kr,kr,mr,a,bb,c,b,tcos,d,w)
      go to 148
  146 continue
      do 147 i=1,mr
         b(i) = fistag*b(i)
  147 continue
  148 continue
      do 149 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  149 continue
  150 continue
      lr = kr
      kr = kr+jr
  151 continue
      go to (152,153),mixbnd
  152 nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 155
      go to 154
  153 nr = nlast/jr
      if (nr .le. 1) go to 192
  154 i2r = jr
      nrodpr = nrod
      go to 104
  155 continue
c
c      begin solution
c
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 184
      if (lr .ne. 0) go to 170
      if (n .ne. 3) go to 161
c
c     case n = 3.
c
      go to (156,168),istag
  156 continue
      do 157 i=1,mr
         b(i) = q(i,2)
  157 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 158 i=1,mr
         q(i,2) = b(i)
         b(i) = 4.*b(i)+q(i,1)+2.*q(i,3)
  158 continue
      tcos(1) = -2.
      tcos(2) = 2.
      i1 = 2
      i2 = 0
      call trix (i1,i2,mr,a,bb,c,b,tcos,d,w)
      do 159 i=1,mr
         q(i,2) = q(i,2)+b(i)
         b(i) = q(i,1)+2.*q(i,2)
  159 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 160 i=1,mr
         q(i,1) = b(i)
  160 continue
      jr = 1
      i2r = 0
      go to 194
c
c     case n = 2**p+1
c
  161 continue
      go to (162,170),istag
  162 continue
      do 163 i=1,mr
         b(i) = q(i,j)+.5*q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  163 continue
      call cosgen (jr,1,0.5,0.0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 164 i=1,mr
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
         b(i) = q(i,1)+2.*q(i,nlast)+4.*q(i,j)
  164 continue
      jr2 = 2*jr
      call cosgen (jr,1,0.0,0.0,tcos)
      do 165 i=1,jr
         i1 = jr+i
         i2 = jr+1-i
         tcos(i1) = -tcos(i2)
  165 continue
      call trix (jr2,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  166 continue
      call cosgen (jr,1,0.5,0.0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 167 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  167 continue
      go to 194
c
c     case of general n with nr = 3 .
c
  168 do 169 i=1,mr
         b(i) = q(i,2)
         q(i,2) = 0.
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  169 continue
      jr = 1
      i2r = 0
      j = 2
      go to 177
  170 continue
      do 171 i=1,mr
         b(i) = .5*q(i,1)-q(i,jm1)+q(i,j)
  171 continue
      if (nrod .ne. 0) go to 173
      do 172 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  172 continue
      go to 175
  173 do 174 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  174 continue
  175 continue
      do 176 i=1,mr
         t = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+2.*t
  176 continue
  177 continue
      k1 = kr+2*jr-1
      k2 = kr+jr
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (k2+istag-2,1,0.0,fnum,tcos(k4))
      k4 = k1+k2+1
      call cosgen (jr-1,1,0.0,1.0,tcos(k4))
      call merge (tcos,k1,k2,k1+k2,jr-1,0)
      k3 = k1+k2+lr
      call cosgen (jr,1,0.5,0.0,tcos(k3+1))
      k4 = k3+jr+1
      call cosgen (kr,1,0.5,fden,tcos(k4))
      call merge (tcos,k3,jr,k3+jr,kr,k1)
      if (lr .eq. 0) go to 178
      call cosgen (lr,1,0.5,fden,tcos(k4))
      call merge (tcos,k3,jr,k3+jr,lr,k3-lr)
      call cosgen (kr,1,0.5,fden,tcos(k4))
  178 k3 = kr
      k4 = kr
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 179 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  179 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 180 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  180 continue
      call cosgen (jr,1,0.5,0.0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 182
      do 181 i=1,mr
         q(i,1) = b(i)
  181 continue
      go to 194
  182 continue
      do 183 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  183 continue
      go to 194
  184 continue
      if (n .ne. 2) go to 188
c
c     case  n = 2
c
      do 185 i=1,mr
         b(i) = q(i,1)
  185 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 186 i=1,mr
         q(i,1) = b(i)
         b(i) = 2.*(q(i,2)+b(i))*fistag
  186 continue
      tcos(1) = -fistag
      tcos(2) = 2.
      call trix (2,0,mr,a,bb,c,b,tcos,d,w)
      do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
      jr = 1
      i2r = 0
      go to 194
  188 continue
c
c     case of general n and nr = 2 .
c
      do 189 i=1,mr
         ii = ip+i
         b3(i) = 0.
         b(i) = q(i,1)+2.*p(ii)
         q(i,1) = .5*q(i,1)-q(i,jm1)
         b2(i) = 2.*(q(i,1)+q(i,nlast))
  189 continue
      k1 = kr+jr-1
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (kr+istag-2,1,0.0,fnum,tcos(k4))
      k4 = k1+kr+1
      call cosgen (jr-1,1,0.0,1.0,tcos(k4))
      call merge (tcos,k1,kr,k1+kr,jr-1,0)
      call cosgen (kr,1,0.5,fden,tcos(k1+1))
      k2 = kr
      k4 = k1+k2+1
      call cosgen (lr,1,0.5,fden,tcos(k4))
      k3 = lr
      k4 = 0
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 190 i=1,mr
         b(i) = b(i)+b2(i)
  190 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 191 i=1,mr
         q(i,1) = q(i,1)+b(i)
  191 continue
      go to 194
  192 do 193 i=1,mr
         b(i) = q(i,nlast)
  193 continue
      go to 196
  194 continue
c
c     start back substitution.
c
      j = nlast-jr
      do 195 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  195 continue
  196 jm2 = nlast-i2r
      if (jr .ne. 1) go to 198
      do 197 i=1,mr
         q(i,nlast) = 0.
  197 continue
      go to 202
  198 continue
      if (nrod .ne. 0) go to 200
      do 199 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  199 continue
      ip = ip-mr
      go to 202
  200 do 201 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  201 continue
  202 continue
      call cosgen (kr,1,0.5,fden,tcos)
      call cosgen (lr,1,0.5,fden,tcos(kr+1))
      if (lr .ne. 0) go to 204
      do 203 i=1,mr
         b(i) = fistag*b(i)
  203 continue
  204 continue
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 205 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  205 continue
      nlastp = nlast
  206 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 222
      go to (207,208),mixbnd
  207 jstart = 1+jr
      go to 209
  208 jstart = jr
  209 continue
      kr = kr-jr
      if (nlast+jr .gt. n) go to 210
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 211
  210 continue
      jstop = nlast-jr
  211 continue
      lr = kr-jr
      call cosgen (jr,1,0.5,0.0,tcos)
      do 221 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 213
         do 212 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  212    continue
         go to 215
  213    continue
         do 214 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  214    continue
  215    continue
         if (jr .ne. 1) go to 217
         do 216 i=1,mr
            q(i,j) = 0.
  216    continue
         go to 219
  217    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 218 i=1,mr
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  218    continue
  219    continue
         call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 220 i=1,mr
            q(i,j) = q(i,j)+b(i)
  220    continue
  221 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 194
      go to 206
  222 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine poisp2 (m,n,a,bb,c,q,idimq,b,b2,b3,w,w2,w3,d,tcos,p)

c*********************************************************************72
c
cc POISP2 solves Poisson's equation with periodic boundary conditions.
c
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) ,
     &                b(*)       ,b2(*)      ,b3(*)      ,w(*)       ,
     &                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    ,
     &                p(*)
      mr = m
      nr = (n+1)/2
      nrm1 = nr-1
      if (2*nr .ne. n) go to 107
c
c     even number of unknowns
c
      do 102 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 101 i=1,mr
            s = q(i,nrmj)-q(i,nrpj)
            t = q(i,nrmj)+q(i,nrpj)
            q(i,nrmj) = s
            q(i,nrpj) = t
  101    continue
  102 continue
      do 103 i=1,mr
         q(i,nr) = 2.*q(i,nr)
         q(i,n) = 2.*q(i,n)
  103 continue
      call poisd2 (mr,nrm1,1,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr+1,1,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d,
     &             tcos,p)
      ipstor = max0(ipstor,int(w(1)))
      do 105 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 104 i=1,mr
            s = .5*(q(i,nrpj)+q(i,nrmj))
            t = .5*(q(i,nrpj)-q(i,nrmj))
            q(i,nrmj) = s
            q(i,nrpj) = t
  104    continue
  105 continue
      do 106 i=1,mr
         q(i,nr) = .5*q(i,nr)
         q(i,n) = .5*q(i,n)
  106 continue
      go to 118
  107 continue
c
c     odd  number of unknowns
c
      do 109 j=1,nrm1
         nrpj = n+1-j
         do 108 i=1,mr
            s = q(i,j)-q(i,nrpj)
            t = q(i,j)+q(i,nrpj)
            q(i,j) = s
            q(i,nrpj) = t
  108    continue
  109 continue
      do 110 i=1,mr
         q(i,nr) = 2.*q(i,nr)
  110 continue
      lh = nrm1/2
      do 112 j=1,lh
         nrmj = nr-j
         do 111 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  111    continue
  112 continue
      call poisd2 (mr,nrm1,2,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr,2,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d,
     &             tcos,p)
      ipstor = max0(ipstor,int(w(1)))
      do 114 j=1,nrm1
         nrpj = nr+j
         do 113 i=1,mr
            s = .5*(q(i,nrpj)+q(i,j))
            t = .5*(q(i,nrpj)-q(i,j))
            q(i,nrpj) = t
            q(i,j) = s
  113    continue
  114 continue
      do 115 i=1,mr
         q(i,nr) = .5*q(i,nr)
  115 continue
      do 117 j=1,lh
         nrmj = nr-j
         do 116 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  116    continue
  117 continue
  118 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine poistg (nperod,n,mperod,m,a,b,c,idimy,y,ierror,w)

c*********************************************************************72
c
cc POISTG solves a special set of linear equations.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c
c    poistg solves the linear system of equations
c
c       a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c       + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c       for i=1,2,...,m and j=1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     x(i,1) or -x(i,1) and x(i,n+1) may be equal to x(i,n) or -x(i,n)
c     depending on an input parameter.
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c   nperod
c     indicates the values which x(i,0) and x(i,n+1) are assumed
c     to have.
c     = 1 if x(i,0) = -x(i,1) and x(i,n+1) = -x(i,n)
c     = 2 if x(i,0) = -x(i,1) and x(i,n+1) =  x(i,n)
c     = 3 if x(i,0) =  x(i,1) and x(i,n+1) =  x(i,n)
c     = 4 if x(i,0) =  x(i,1) and x(i,n+1) = -x(i,n)
c
c   n
c     the number of unknowns in the j-direction.  n must
c     be greater than 2.
c
c   mperod
c     = 0 if a(1) and c(m) are not zero
c     = 1 if a(1) = c(m) = 0
c
c   m
c     the number of unknowns in the i-direction.  m must
c     be greater than 2.
c
c   a,b,c
c     one-dimensional arrays of length m that specify the coefficients
c     in the linear equations given above.  if mperod = 0 the array
c     elements must not depend on the index i, but must be constant.
c     specifically, the subroutine checks the following condition
c
c           a(i) = c(1)
c           b(i) = b(1)
c           c(i) = c(1)
c
c     for i = 1, 2, ..., m.
c
c   idimy
c     the row (or first) dimension of the two-dimensional array y as
c     it appears in the program calling poistg.  this parameter is
c     used to specify the variable dimension of y.  idimy must be at
c     least m.
c
c   y
c     a two-dimensional array that specifies the values of the
c     right side of the linear system of equations given above.
c     y must be dimensioned at least m x n.
c
c   w
c     a one-dimensional work array that must be provided by the user
c     for work space.  w may require up to 9m + 4n + m(int(log2(n)))
c     locations.  the actual number of locations used is computed by
c     poistg and returned in location w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c   y
c     contains the solution x.
c
c   ierror
c     an error flag that indicates invalid input parameters.  except
c     for number zero, a solution is not attempted.
c     = 0  no error
c     = 1  if m .le. 2
c     = 2  if n .le. 2
c     = 3  idimy .lt. m
c     = 4  if nperod .lt. 1 or nperod .gt. 4
c     = 5  if mperod .lt. 0 or mperod .gt. 1
c     = 6  if mperod = 0 and
c          a(i) .ne. c(1) or b(i) .ne. b(1) or c(i) .ne. c(1)
c          for some i = 1, 2, ..., m.
c       = 7 if mperod .eq. 1 .and. (a(1).ne.0 .or. c(m).ne.0)
c
c   w
c     w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(m),b(m),c(m),y(idimy,n),
c     arguments      w(see argument list)
c
c     latest         june 1, 1977
c     revision
c
c     subprograms    poistg,postg2,cosgen,merge,trix,tri3,pimach
c     required
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet in 1973
c                    revised by roland sweet in 1977
c
c     algorithm      this subroutine is an implementation of the
c                    algorithm presented in the reference.
c
c     space          3297(decimal) = 6341(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine poistg is roughly proportional
c                    to m*n*log2(n).  some typical values are listed
c                    in the table below.  more comprehensive timing
c                    charts may be found in the reference.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(i) = c(i) = -0.5*b(i) = 1,       i=1,2,...,m
c
c                    and, when mperod = 1
c
c                       a(1) = c(m) = 0
c                       b(1) = b(m) =-1.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine poistg was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                       e = max(abs(z(i,j)-x(i,j)))/max(abs(x(i,j)))
c
c                    where the two maxima are taken over all i=1,2,...,m
c                    and j=1,2,...,n, was computed.  the value of e is
c                    given in the table below for some typical values of
c                    m and n.
c
c
c                       m (=n)    mperod    nperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         31        0-1       1-4        45     9.e-13
c                         31        1         1          21     4.e-13
c                         31        1         3          41     3.e-13
c                         32        0-1       1-4        51     3.e-12
c                         32        1         1          32     3.e-13
c                         32        1         3          48     1.e-13
c                         33        0-1       1-4        42     1.e-12
c                         33        1         1          30     4.e-13
c                         33        1         3          34     1.e-13
c                         63        0-1       1-4       186     3.e-12
c                         63        1         1          91     1.e-12
c                         63        1         3         173     2.e-13
c                         64        0-1       1-4       209     4.e-12
c                         64        1         1         128     1.e-12
c                         64        1         3         199     6.e-13
c                         65        0-1       1-4       143     2.e-13
c                         65        1         1         160     1.e-11
c                         65        1         3         138     4.e-13
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      schumann, u. and r. sweet,"a direct method for
c                    the solution of poisson"s equation with neumann
c                    boundary conditions on a staggered grid of
c                    arbitrary size," j. comp. phys. 20(1976),
c                    pp. 171-182.
c
      dimension       y(idimy,*)
      dimension       w(*)       ,b(*)       ,a(*)       ,c(*)
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.1 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 103
      do 101 i=1,m
         if (a(i) .ne. c(1)) go to 102
         if (c(i) .ne. c(1)) go to 102
         if (b(i) .ne. b(1)) go to 102
  101 continue
      go to 104
  102 ierror = 6
      return
  103 if (a(1).ne.0. .or. c(m).ne.0.) ierror = 7
  104 if (ierror .ne. 0) return
      iwba = m+1
      iwbb = iwba+m
      iwbc = iwbb+m
      iwb2 = iwbc+m
      iwb3 = iwb2+m
      iww1 = iwb3+m
      iww2 = iww1+m
      iww3 = iww2+m
      iwd = iww3+m
      iwtcos = iwd+m
      iwp = iwtcos+4*n
      do 106 i=1,m
         k = iwba+i-1
         w(k) = -a(i)
         k = iwbc+i-1
         w(k) = -c(i)
         k = iwbb+i-1
         w(k) = 2.-b(i)
         do 105 j=1,n
            y(i,j) = -y(i,j)
  105    continue
  106 continue
      np = nperod
      mp = mperod+1
      go to (110,107),mp
  107 continue
      go to (108,108,108,119),nperod
  108 continue
      call postg2 (np,n,m,w(iwba),w(iwbb),w(iwbc),idimy,y,w,w(iwb2),
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     &             w(iwp))
      ipstor = w(iww1)
      irev = 2
      if (nperod .eq. 4) go to 120
  109 continue
      go to (123,129),mp
  110 continue
c
c     reorder unknowns when mp =0
c
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 115 j=1,n
         do 111 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  111    continue
         w(mh) = 2.*y(mh,j)
         go to (113,112),modd
  112    w(m) = 2.*y(m,j)
  113    continue
         do 114 i=1,m
            y(i,j) = w(i)
  114    continue
  115 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = 0.
      w(i) = 0.
      w(k+1) = 2.*w(k+1)
      go to (116,117),modd
  116 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 118
  117 w(iwbb-1) = w(k+1)
  118 continue
      go to 107
  119 continue
c
c     reverse columns when nperod = 4.
c
      irev = 1
      nby2 = n/2
      np = 2
  120 do 122 j=1,nby2
         mskip = n+1-j
         do 121 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  121    continue
  122 continue
      go to (108,109),irev
  123 continue
      do 128 j=1,n
         do 124 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5*(y(mhpi,j)-y(i,j))
  124    continue
         w(mh) = .5*y(mh,j)
         go to (126,125),modd
  125    w(m) = .5*y(m,j)
  126    continue
         do 127 i=1,m
            y(i,j) = w(i)
  127    continue
  128 continue
  129 continue
c
c     return storage requirements for w array.
c
      w(1) = ipstor+iwp-1
      return
      end
      subroutine pos3d1 (lp,l,mp,m,n,a,b,c,ldimf,mdimf,f,xrt,yrt,t,d,
     &                   wx,wy,c1,c2,bb)

c*********************************************************************72
c
cc POS3D1
c
      dimension       a(*)       ,b(*)       ,c(*)       ,
     &                f(ldimf,mdimf,*)       ,xrt(*)     ,yrt(*)     ,
     &                t(*)       ,d(*)       ,wx(*)      ,wy(*)      ,
     &                bb(*)
      pi = pimach()
      lr = l
      mr = m
      nr = n
c
c     generate transform roots
c
      lrdel = ((lp-1)*(lp-3)*(lp-5))/3
      scalx = lr+lrdel
      dx = pi/(2.*scalx)
      go to (108,103,101,102,101),lp
  101 di = 0.5
      scalx = 2.*scalx
      go to 104
  102 di = 1.0
      go to 104
  103 di = 0.0
  104 do 105 i=1,lr
         xrt(i) = -4.*c1*(sin((float(i)-di)*dx))**2
  105 continue
      scalx = 2.*scalx
      go to (112,106,110,107,111),lp
  106 call sinti (lr,wx)
      go to 112
  107 call costi (lr,wx)
      go to 112
  108 xrt(1) = 0.
      xrt(lr) = -4.*c1
      do 109 i=3,lr,2
         xrt(i-1) = -4.*c1*(sin(float((i-1))*dx))**2
         xrt(i) = xrt(i-1)
  109 continue
      call rffti (lr,wx)
      go to 112
  110 call sinqi (lr,wx)
      go to 112
  111 call cosqi (lr,wx)
  112 continue
      mrdel = ((mp-1)*(mp-3)*(mp-5))/3
      scaly = mr+mrdel
      dy = pi/(2.*scaly)
      go to (120,115,113,114,113),mp
  113 dj = 0.5
      scaly = 2.*scaly
      go to 116
  114 dj = 1.0
      go to 116
  115 dj = 0.0
  116 do 117 j=1,mr
         yrt(j) = -4.*c2*(sin((float(j)-dj)*dy))**2
  117 continue
      scaly = 2.*scaly
      go to (124,118,122,119,123),mp
  118 call sinti (mr,wy)
      go to 124
  119 call costi (mr,wy)
      go to 124
  120 yrt(1) = 0.
      yrt(mr) = -4.*c2
      do 121 j=3,mr,2
         yrt(j-1) = -4.*c2*(sin(float((j-1))*dy))**2
         yrt(j) = yrt(j-1)
  121 continue
      call rffti (mr,wy)
      go to 124
  122 call sinqi (mr,wy)
      go to 124
  123 call cosqi (mr,wy)
  124 continue
      ifwrd = 1
      is = 1
  125 continue
c
c     transform x
c
      do 141 j=1,mr
         do 140 k=1,nr
            do 126 i=1,lr
               t(i) = f(i,j,k)
  126       continue
            go to (127,130,131,134,135),lp
  127       go to (128,129),ifwrd
  128       call rfftf (lr,t,wx)
            go to 138
  129       call rfftb (lr,t,wx)
            go to 138
  130       call sint (lr,t,wx)
            go to 138
  131       go to (132,133),ifwrd
  132       call sinqf (lr,t,wx)
            go to 138
  133       call sinqb (lr,t,wx)
            go to 138
  134       call cost (lr,t,wx)
            go to 138
  135       go to (136,137),ifwrd
  136       call cosqf (lr,t,wx)
            go to 138
  137       call cosqb (lr,t,wx)
  138       continue
            do 139 i=1,lr
               f(i,j,k) = t(i)
  139       continue
  140    continue
  141 continue
      go to (142,164),ifwrd
c
c     transform y
c
  142 continue
      do 158 i=1,lr
         do 157 k=1,nr
            do 143 j=1,mr
               t(j) = f(i,j,k)
  143       continue
            go to (144,147,148,151,152),mp
  144       go to (145,146),ifwrd
  145       call rfftf (mr,t,wy)
            go to 155
  146       call rfftb (mr,t,wy)
            go to 155
  147       call sint (mr,t,wy)
            go to 155
  148       go to (149,150),ifwrd
  149       call sinqf (mr,t,wy)
            go to 155
  150       call sinqb (mr,t,wy)
            go to 155
  151       call cost (mr,t,wy)
            go to 155
  152       go to (153,154),ifwrd
  153       call cosqf (mr,t,wy)
            go to 155
  154       call cosqb (mr,t,wy)
  155       continue
            do 156 j=1,mr
               f(i,j,k) = t(j)
  156       continue
  157    continue
  158 continue
      go to (159,125),ifwrd
  159 continue
c
c     solve tridiagonal systems in z
c
      do 163 i=1,lr
         do 162 j=1,mr
            do 160 k=1,nr
               bb(k) = b(k)+xrt(i)+yrt(j)
               t(k) = f(i,j,k)
  160       continue
            call trid (nr,a,bb,c,t,d)
            do 161 k=1,nr
               f(i,j,k) = t(k)
  161       continue
  162    continue
  163 continue
      ifwrd = 2
      is = -1
      go to 142
  164 continue
      do 167 i=1,lr
         do 166 j=1,mr
            do 165 k=1,nr
               f(i,j,k) = f(i,j,k)/(scalx*scaly)
  165       continue
  166    continue
  167 continue
      return
      end
      subroutine postg2 (nperod,n,m,a,bb,c,idimq,q,b,b2,b3,w,w2,w3,d,
     &                   tcos,p)

c*********************************************************************72
c
cc POSTG2 solves Poisson's equation on a staggered grid.
c
c
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) ,
     &                b(*)       ,b2(*)      ,b3(*)      ,w(*)       ,
     &                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    ,
     &                k(4)       ,p(*)
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
      np = nperod
      fnum = 0.5*float(np/3)
      fnum2 = 0.5*float(np/2)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      if (nr .le. 3) go to 142
  101 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      jstart = 1
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 102
      j = jr
      go to 115
  102 continue
c
c     regular reduction.
c
      ijump = 1
      do 114 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 106
         call cosgen (i2r,1,fnum,0.5,tcos)
         if (i2r .ne. 1) go to 104
         do 103 i=1,mr
            b(i) = q(i,1)
            q(i,1) = q(i,2)
  103    continue
         go to 112
  104    do 105 i=1,mr
            b(i) = q(i,1)+0.5*(q(i,jp2)-q(i,jp1)-q(i,jp3))
            q(i,1) = q(i,jp2)+q(i,1)-q(i,jp1)
  105    continue
         go to 112
  106    continue
         go to (107,108),ijump
  107    continue
         ijump = 2
         call cosgen (i2r,1,0.5,0.0,tcos)
  108    continue
         if (i2r .ne. 1) go to 110
         do 109 i=1,mr
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  109    continue
         go to 112
  110    do 111 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  111    continue
  112    continue
         call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 113 i=1,mr
            q(i,j) = q(i,j)+b(i)
  113    continue
c
c     end of reduction for regular unknowns.
c
  114 continue
c
c     begin special reduction for last unknown.
c
      j = jstop+jr
  115 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 125
c
c     odd number of unknowns
c
      if (i2r .ne. 1) go to 117
      do 116 i=1,mr
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  116 continue
      go to 123
  117 do 118 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  118 continue
      if (nrodpr .ne. 0) go to 120
      do 119 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  119 continue
      ip = ip-mr
      go to 122
  120 continue
      do 121 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  121 continue
  122 if (lr .eq. 0) go to 123
      call cosgen (lr,1,fnum2,0.5,tcos(kr+1))
  123 continue
      call cosgen (kr,1,fnum2,0.5,tcos)
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 124 i=1,mr
         q(i,j) = q(i,j)+b(i)
  124 continue
      kr = kr+i2r
      go to 141
  125 continue
c
c     even number of unknowns
c
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 129
      do 126 i=1,mr
         b(i) = q(i,j)
  126 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      do 127 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  127 continue
      tcos(1) = -1.+2.*float(np/2)
      tcos(2) = 0.
      call trix (1,1,mr,a,bb,c,b,tcos,d,w)
      do 128 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  128 continue
      go to 140
  129 continue
      do 130 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  130 continue
      if (nrodpr .ne. 0) go to 132
      do 131 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  131 continue
      go to 134
  132 continue
      do 133 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  133 continue
  134 continue
      call cosgen (i2r,1,0.5,0.0,tcos)
      call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max0(ipstor,ip+mr)
      do 135 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  135 continue
      if (lr .eq. 0) go to 136
      call cosgen (lr,1,fnum2,0.5,tcos(i2r+1))
      call merge (tcos,0,i2r,i2r,lr,kr)
      go to 138
  136 do 137 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  137 continue
  138 call cosgen (kr,1,fnum2,0.5,tcos)
      call trix (kr,kr,mr,a,bb,c,b,tcos,d,w)
      do 139 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  139 continue
  140 continue
      lr = kr
      kr = kr+jr
  141 continue
      nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 142
      i2r = jr
      nrodpr = nrod
      go to 101
  142 continue
c
c      begin solution
c
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 180
      if (lr .ne. 0) go to 167
      if (n .ne. 3) go to 156
c
c     case n = 3.
c
      go to (143,148,143),np
  143 do 144 i=1,mr
         b(i) = q(i,2)
         b2(i) = q(i,1)+q(i,3)
         b3(i) = 0.
  144 continue
      go to (146,146,145),np
  145 tcos(1) = -1.
      tcos(2) = 1.
      k1 = 1
      go to 147
  146 tcos(1) = -2.
      tcos(2) = 1.
      tcos(3) = -1.
      k1 = 2
  147 k2 = 1
      k3 = 0
      k4 = 0
      go to 150
  148 do 149 i=1,mr
         b(i) = q(i,2)
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  149 continue
      call cosgen (3,1,0.5,0.0,tcos)
      tcos(4) = -1.
      tcos(5) = 1.
      tcos(6) = -1.
      tcos(7) = 1.
      k1 = 3
      k2 = 2
      k3 = 1
      k4 = 1
  150 call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 151 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  151 continue
      go to (153,153,152),np
  152 tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  153 do 154 i=1,mr
         q(i,2) = b(i)
         b(i) = q(i,1)+b(i)
  154 continue
      tcos(1) = -1.+4.*fnum
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 155 i=1,mr
         q(i,1) = b(i)
  155 continue
      jr = 1
      i2r = 0
      go to 188
c
c     case n = 2**p+1
c
  156 continue
      do 157 i=1,mr
         b(i) = q(i,j)+q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  157 continue
      go to (158,160,158),np
  158 do 159 i=1,mr
         b2(i) = q(i,1)+q(i,nlast)+q(i,j)-q(i,jm1)-q(i,jp1)
         b3(i) = 0.
  159 continue
      k1 = nlast-1
      k2 = nlast+jr-1
      call cosgen (jr-1,1,0.0,1.0,tcos(nlast))
      tcos(k2) = 2.*float(np-2)
      call cosgen (jr,1,0.5-fnum,0.5,tcos(k2+1))
      k3 = (3-np)/2
      call merge (tcos,k1,jr-k3,k2-k3,jr+k3,0)
      k1 = k1-1+k3
      call cosgen (jr,1,fnum,0.5,tcos(k1+1))
      k2 = jr
      k3 = 0
      k4 = 0
      go to 162
  160 do 161 i=1,mr
         fi = (q(i,j)-q(i,jm1)-q(i,jp1))/2.
         b2(i) = q(i,1)+fi
         b3(i) = q(i,nlast)+fi
  161 continue
      k1 = nlast+jr-1
      k2 = k1+jr-1
      call cosgen (jr-1,1,0.0,1.0,tcos(k1+1))
      call cosgen (nlast,1,0.5,0.0,tcos(k2+1))
      call merge (tcos,k1,jr-1,k2,nlast,0)
      k3 = k1+nlast-1
      k4 = k3+jr
      call cosgen (jr,1,0.5,0.5,tcos(k3+1))
      call cosgen (jr,1,0.0,0.5,tcos(k4+1))
      call merge (tcos,k3,jr,k4,jr,k1)
      k2 = nlast-1
      k3 = jr
      k4 = jr
  162 call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 163 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  163 continue
      if (np .ne. 3) go to 164
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  164 do 165 i=1,mr
         q(i,j) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = q(i,j)+q(i,1)
  165 continue
      call cosgen (jr,1,fnum,0.5,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,1) = q(i,1)-q(i,jm1)+b(i)
  166 continue
      go to 188
c
c     case of general n with nr = 3 .
c
  167 continue
      do 168 i=1,mr
         b(i) = q(i,1)-q(i,jm1)+q(i,j)
  168 continue
      if (nrod .ne. 0) go to 170
      do 169 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  169 continue
      go to 172
  170 do 171 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  171 continue
  172 continue
      do 173 i=1,mr
         t = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+t
  173 continue
      k1 = kr+2*jr
      call cosgen (jr-1,1,0.0,1.0,tcos(k1+1))
      k2 = k1+jr
      tcos(k2) = 2.*float(np-2)
      k4 = (np-1)*(3-np)
      k3 = k2+1-k4
      call cosgen (kr+jr+k4,1,float(k4)/2.,1.-float(k4),tcos(k3))
      k4 = 1-np/3
      call merge (tcos,k1,jr-k4,k2-k4,kr+jr+k4,0)
      if (np .eq. 3) k1 = k1-1
      k2 = kr+jr
      k4 = k1+k2
      call cosgen (kr,1,fnum2,0.5,tcos(k4+1))
      k3 = k4+kr
      call cosgen (jr,1,fnum,0.5,tcos(k3+1))
      call merge (tcos,k4,kr,k3,jr,k1)
      k4 = k3+jr
      call cosgen (lr,1,fnum2,0.5,tcos(k4+1))
      call merge (tcos,k3,jr,k4,lr,k1+k2)
      call cosgen (kr,1,fnum2,0.5,tcos(k3+1))
      k3 = kr
      k4 = kr
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 174 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  174 continue
      if (np .ne. 3) go to 175
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  175 do 176 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+q(i,j)
  176 continue
      call cosgen (jr,1,fnum,0.5,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 178
      do 177 i=1,mr
         q(i,1) = b(i)
  177 continue
      go to 188
  178 continue
      do 179 i=1,mr
         q(i,1) = q(i,1)-q(i,jm1)+b(i)
  179 continue
      go to 188
  180 continue
c
c     case of general n and nr = 2 .
c
      do 181 i=1,mr
         ii = ip+i
         b3(i) = 0.
         b(i) = q(i,1)+p(ii)
         q(i,1) = q(i,1)-q(i,jm1)
         b2(i) = q(i,1)+q(i,nlast)
  181 continue
      k1 = kr+jr
      k2 = k1+jr
      call cosgen (jr-1,1,0.0,1.0,tcos(k1+1))
      go to (182,183,182),np
  182 tcos(k2) = 2.*float(np-2)
      call cosgen (kr,1,0.0,1.0,tcos(k2+1))
      go to 184
  183 call cosgen (kr+1,1,0.5,0.0,tcos(k2))
  184 k4 = 1-np/3
      call merge (tcos,k1,jr-k4,k2-k4,kr+k4,0)
      if (np .eq. 3) k1 = k1-1
      k2 = kr
      call cosgen (kr,1,fnum2,0.5,tcos(k1+1))
      k4 = k1+kr
      call cosgen (lr,1,fnum2,0.5,tcos(k4+1))
      k3 = lr
      k4 = 0
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 185 i=1,mr
         b(i) = b(i)+b2(i)
  185 continue
      if (np .ne. 3) go to 186
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
  186 do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
  188 continue
c
c     start back substitution.
c
      j = nlast-jr
      do 189 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  189 continue
      jm2 = nlast-i2r
      if (jr .ne. 1) go to 191
      do 190 i=1,mr
         q(i,nlast) = 0.
  190 continue
      go to 195
  191 continue
      if (nrod .ne. 0) go to 193
      do 192 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  192 continue
      ip = ip-mr
      go to 195
  193 do 194 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  194 continue
  195 continue
      call cosgen (kr,1,fnum2,0.5,tcos)
      call cosgen (lr,1,fnum2,0.5,tcos(kr+1))
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 196 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  196 continue
      nlastp = nlast
  197 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 210
      jstart = 1+jr
      kr = kr-jr
      if (nlast+jr .gt. n) go to 198
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 199
  198 continue
      jstop = nlast-jr
  199 continue
      lr = kr-jr
      call cosgen (jr,1,0.5,0.0,tcos)
      do 209 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 201
         do 200 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  200    continue
         go to 203
  201    continue
         do 202 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  202    continue
  203    continue
         if (jr .ne. 1) go to 205
         do 204 i=1,mr
            q(i,j) = 0.
  204    continue
         go to 207
  205    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 206 i=1,mr
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  206    continue
  207    continue
         call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 208 i=1,mr
            q(i,j) = q(i,j)+b(i)
  208    continue
  209 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 188
      go to 197
  210 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine ppadd (n,ierror,a,c,cbp,bp,bh)

c*********************************************************************72
c
cc PPADD computes the eigenvalues of the periodic tridiagonal matrix
c     with coefficients an,bn,cn
c
c n is the order of the bh and bp polynomials
c     on output bp contians the eigenvalues
c cbp is the same as bp except type complex
c bh is used to temporarily store the roots of the b hat polynomial
c which enters through bp
c
      complex         cf         ,cx         ,fsg        ,hsg        ,
     &                dd         ,f          ,fp         ,fpp        ,
     &                cdis       ,r1         ,r2         ,r3         ,
     &                cbp
      dimension       a(*)       ,c(*)       ,bp(*)      ,bh(*)      ,
     &                cbp(*)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: cblkt
!$omp threadprivate (cblkt)
      external        psgf       ,ppspf      ,ppsgf
      scnv = sqrt(cnv)
      iz = n
      izm = iz-1
      izm2 = iz-2
      if (bp(n)-bp(1)) 101,142,103
  101 do 102 j=1,n
         nt = n-j
         bh(j) = bp(nt+1)
  102 continue
      go to 105
  103 do 104 j=1,n
         bh(j) = bp(j)
  104 continue
  105 ncmplx = 0
      modiz = mod(iz,2)
      is = 1
      if (modiz) 106,107,106
  106 if (a(1)) 110,142,107
  107 xl = bh(1)
      db = bh(3)-bh(1)
  108 xl = xl-db
      if (psgf(xl,iz,c,a,bh)) 108,108,109
  109 sgn = -1.
      cbp(1) = cmplx(bsrh(xl,bh(1),iz,c,a,bh,psgf,sgn),0.)
      is = 2
  110 if = iz-1
      if (modiz) 111,112,111
  111 if (a(1)) 112,142,115
  112 xr = bh(iz)
      db = bh(iz)-bh(iz-2)
  113 xr = xr+db
      if (psgf(xr,iz,c,a,bh)) 113,114,114
  114 sgn = 1.
      cbp(iz) = cmplx(bsrh(bh(iz),xr,iz,c,a,bh,psgf,sgn),0.)
      if = iz-2
  115 do 136 ig=is,if,2
         xl = bh(ig)
         xr = bh(ig+1)
         sgn = -1.
         xm = bsrh(xl,xr,iz,c,a,bh,ppspf,sgn)
         psg = psgf(xm,iz,c,a,bh)
         if (abs(psg)-eps) 118,118,116
  116    if (psg*ppsgf(xm,iz,c,a,bh)) 117,118,119
c
c     case of a real zero
c
  117    sgn = 1.
         cbp(ig) = cmplx(bsrh(bh(ig),xm,iz,c,a,bh,psgf,sgn),0.)
         sgn = -1.
         cbp(ig+1) = cmplx(bsrh(xm,bh(ig+1),iz,c,a,bh,psgf,sgn),0.)
         go to 136
c
c     case of a multiple zero
c
  118    cbp(ig) = cmplx(xm,0.)
         cbp(ig+1) = cmplx(xm,0.)
         go to 136
c
c     case of a complex zero
c
  119    it = 0
         icv = 0
         cx = cmplx(xm,0.)
  120    fsg = (1.,0.)
         hsg = (1.,0.)
         fp = (0.,0.)
         fpp = (0.,0.)
         do 121 j=1,iz
            dd = 1./(cx-bh(j))
            fsg = fsg*a(j)*dd
            hsg = hsg*c(j)*dd
            fp = fp+dd
            fpp = fpp-dd*dd
  121    continue
         if (modiz) 123,122,123
  122    f = (1.,0.)-fsg-hsg
         go to 124
  123    f = (1.,0.)+fsg+hsg
  124    i3 = 0
         if (cabs(fp)) 126,126,125
  125    i3 = 1
         r3 = -f/fp
  126    i2 = 0
         if (cabs(fpp)) 132,132,127
  127    i2 = 1
         cdis = csqrt(fp**2-2.*f*fpp)
         r1 = cdis-fp
         r2 = -fp-cdis
         if (cabs(r1)-cabs(r2)) 129,129,128
  128    r1 = r1/fpp
         go to 130
  129    r1 = r2/fpp
  130    r2 = 2.*f/fpp/r1
         if (cabs(r2) .lt. cabs(r1)) r1 = r2
         if (i3) 133,133,131
  131    if (cabs(r3) .lt. cabs(r1)) r1 = r3
         go to 133
  132    r1 = r3
  133    cx = cx+r1
         it = it+1
         if (it .gt. 50) go to 142
         if (cabs(r1) .gt. scnv) go to 120
         if (icv) 134,134,135
  134    icv = 1
         go to 120
  135    cbp(ig) = cx
         cbp(ig+1) = conjg(cx)
  136 continue
      if (cabs(cbp(n))-cabs(cbp(1))) 137,142,139
  137 nhalf = n/2
      do 138 j=1,nhalf
         nt = n-j
         cx = cbp(j)
         cbp(j) = cbp(nt+1)
         cbp(nt+1) = cx
  138 continue
  139 ncmplx = 1
      do 140 j=2,iz
         if (aimag(cbp(j))) 143,140,143
  140 continue
      ncmplx = 0
      do 141 j=2,iz
         bp(j) = real(cbp(j))
  141 continue
      go to 143
  142 ierror = 4
  143 continue
      return
      end
      function ppgsf (x,iz,c,a,bh)

c*********************************************************************72
c
cc PPGSF
c
      dimension       a(1)       ,c(1)       ,bh(1)
      sum = 0.
      do 101 j=1,iz
         sum = sum-1./(x-bh(j))**2
  101 continue
      ppgsf = sum
      return
      end
      function pppsf (x,iz,c,a,bh)

c*********************************************************************72
c
cc PPPSF
c
      dimension       a(1)       ,c(1)       ,bh(1)
      sum = 0.
      do 101 j=1,iz
         sum = sum+1./(x-bh(j))
  101 continue
      pppsf = sum
      return
      end
      function ppsgf (x,iz,c,a,bh)

c*********************************************************************72
c
cc PPSGF
c
      dimension       a(1)       ,c(1)       ,bh(1)
      sum = 0.
      do 101 j=1,iz
         sum = sum-1./(x-bh(j))**2
  101 continue
      ppsgf = sum
      return
      end
      function ppspf (x,iz,c,a,bh)

c*********************************************************************72
c
cc PPSPF
c
      dimension       a(1)       ,c(1)       ,bh(1)
      sum = 0.
      do 101 j=1,iz
         sum = sum+1./(x-bh(j))
  101 continue
      ppspf = sum
      return
      end
      subroutine proc (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,u)

c*********************************************************************72
c
cc PROC applies a sequence of matrix operations to the vector x and
c stores the result in y
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w,u are working arrays
c is  determines whether or not a change in sign is made
c
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,
     &                y(*)       ,d(*)       ,w(*)       ,bd(*)      ,
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,u(*)
      complex         x          ,y          ,a          ,b          ,
     &                c          ,d          ,w          ,u          ,
     &                den
      do 101 j=1,m
         w(j) = x(j)
         y(j) = w(j)
  101 continue
      mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
c
c scalar multiplication
c
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 125,125,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)
      do 107 j=2,mm
         k = m-j
         den = b(k+1)-rt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  107 continue
      den = b(1)-rt-c(1)*d(2)
      w(1) = (1.,0.)
      if (cabs(den)) 108,109,108
  108 w(1) = (y(1)-c(1)*w(2))/den
  109 do 110 j=2,m
         w(j) = w(j)-d(j)*w(j-1)
  110 continue
      if (na) 113,113,102
  111 do 112 j=1,m
         y(j) = w(j)
  112 continue
      ibr = 1
      go to 102
  113 if (m1) 114,114,115
  114 if (m2) 111,111,120
  115 if (m2) 117,117,116
  116 if (abs(bm1(m1))-abs(bm2(m2))) 120,120,117
  117 if (ibr) 118,118,119
  118 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 111,119,119
  119 rt = rt-bm1(m1)
      m1 = m1-1
      go to 123
  120 if (ibr) 121,121,122
  121 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 111,122,122
  122 rt = rt-bm2(m2)
      m2 = m2-1
  123 do 124 j=1,m
         y(j) = y(j)+rt*w(j)
  124 continue
      go to 102
  125 return
      end
      subroutine procp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,u,w)

c*********************************************************************72
c
cc PROCP applies a sequence of matrix operations to the vector x and
c stores the result in y        periodic boundary conditions
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,w are working arrays
c is  determines whether or not a change in sign is made
c
      dimension       a(*)       ,b(*)       ,c(*)       ,x(*)       ,
     &                y(*)       ,d(*)       ,u(*)       ,bd(*)      ,
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,w(*)
      complex         x          ,y          ,a          ,b          ,
     &                c          ,d          ,u          ,w          ,
     &                den        ,ym         ,v          ,bh         ,am

      do 101 j=1,m
         y(j) = x(j)
         w(j) = y(j)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 128,128,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      bh = b(m)-rt
      ym = y(m)
      den = b(1)-rt
      d(1) = c(1)/den
      u(1) = a(1)/den
      w(1) = y(1)/den
      v = c(m)
      if (mm2-2) 109,107,107
  107 do 108 j=2,mm2
         den = b(j)-rt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         w(j) = (y(j)-a(j)*w(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*w(j-1)
         v = -v*d(j-1)
  108 continue
  109 den = b(m-1)-rt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*w(m-2)
      den = bh-am*d(m-1)
      if (cabs(den)) 110,111,110
  110 w(m) = (ym-am*w(m-1))/den
      go to 112
  111 w(m) = (1.,0.)
  112 w(m-1) = w(m-1)-d(m-1)*w(m)
      do 113 j=2,mm
         k = m-j
         w(k) = w(k)-d(k)*w(k+1)-u(k)*w(m)
  113 continue
      if (na) 116,116,102
  114 do 115 j=1,m
         y(j) = w(j)
  115 continue
      ibr = 1
      go to 102
  116 if (m1) 117,117,118
  117 if (m2) 114,114,123
  118 if (m2) 120,120,119
  119 if (abs(bm1(m1))-abs(bm2(m2))) 123,123,120
  120 if (ibr) 121,121,122
  121 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 114,122,122
  122 rt = rt-bm1(m1)
      m1 = m1-1
      go to 126
  123 if (ibr) 124,124,125
  124 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 114,125,125
  125 rt = rt-bm2(m2)
      m2 = m2-1
  126 do 127 j=1,m
         y(j) = y(j)+rt*w(j)
  127 continue
      go to 102
  128 return
      end
      subroutine prod (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,u)

c*********************************************************************72
c
cc PROD applies a sequence of matrix operations to the vector x and
c stores the result in y
c
c  bd,bm1,bm2 are arrays containing roots of certain b polynomials
c  nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c  aa   array containing scalar multipliers of the vector x
c  na is the length of the array aa
c  x,y  the matrix operations are applied to x and the result is y
c  a,b,c  are arrays which contain the tridiagonal matrix
c  m  is the order of the matrix
c  d,w,u are working arrays
c  is  determines whether or not a change in sign is made
c
      real a(*)
      dimension           b(*)       ,c(*)       ,x(*)       ,
     &                y(*)       ,d(*)       ,w(*)       ,bd(*)      ,
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,u(*)
c
 
      do j=1,m
         w(j) = x(j)
         y(j) = w(j)
      end do
 
      mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
c
c scalar multiplication
c
      do j=1,m
         y(j) = rt*w(j)
      end do
 
105   continue

      if ( id .le. 0 ) then
        return
      end if
 
      rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c  Begin solution.
c
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)

      do j = 2, mm

         k = m-j
         den = b(k+1)-rt-c(k+1)*d(k+2)

         if ( den .eq. 0.0 ) then
           write(*,*)' '
           write(*,*)'PROD - Fatal error!'
           write(*,*)'  Zero denominator, J=',j
           stop
         end if

         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den

      end do
 
      den = b(1)-rt-c(1)*d(2)
 
      if (den .ne. 0.0) then
        w(1) = (y(1)-c(1)*w(2))/den
      else
        w(1) = 1.0
      end if
 
      do j=2,m
        w(j) = w(j)-d(j)*w(j-1)
      end do
 
      if (na) 113,113,102
  111 do j=1,m
         y(j) = w(j)
      end do
 
      ibr = 1
      go to 102
 
  113 if (m1) 114,114,115
  114 if (m2) 111,111,120
  115 if (m2) 117,117,116
  116 if (abs(bm1(m1))-abs(bm2(m2))) 120,120,117
  117 if (ibr) 118,118,119
  118 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 111,119,119
  119 rt = rt-bm1(m1)
      m1 = m1-1
      go to 123
 
  120 if (ibr) 121,121,122
  121 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 111,122,122
  122 rt = rt-bm2(m2)
      m2 = m2-1
  123 do j=1,m
         y(j) = y(j)+rt*w(j)
      end do
      go to 102
 
      end
      subroutine prodp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,u,w)

c*********************************************************************72
c
cc PRODP applies a sequence of matrix operations to the vector x
c  and stores the result in y        periodic boundary conditions
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,w are working arrays
c is  determines whether or not a change in sign is made
c
      dimension       a(1)       ,b(1)       ,c(1)       ,x(1)       ,
     &                y(1)       ,d(1)       ,u(1)       ,bd(1)      ,
     &                bm1(1)     ,bm2(1)     ,aa(1)      ,w(1)
      do j=1,m
         y(j) = x(j)
         w(j) = y(j)
      end do
 
      mm = m-1
      mm2 = m-2
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
 
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
      do j=1,m
         y(j) = rt*w(j)
      end do
 
  105 if (id) 128,128,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      bh = b(m)-rt
      ym = y(m)
      den = b(1)-rt
      d(1) = c(1)/den
      u(1) = a(1)/den
      w(1) = y(1)/den
      v = c(m)
      if (mm2-2) 109,107,107
  107 do j=2,mm2
         den = b(j)-rt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         w(j) = (y(j)-a(j)*w(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*w(j-1)
         v = -v*d(j-1)
      end do
  109 den = b(m-1)-rt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*w(m-2)
      den = bh-am*d(m-1)
      if (den) 110,111,110
  110 w(m) = (ym-am*w(m-1))/den
      go to 112
  111 w(m) = 1.
  112 w(m-1) = w(m-1)-d(m-1)*w(m)
 
      do j=2,mm
         k = m-j
         w(k) = w(k)-d(k)*w(k+1)-u(k)*w(m)
      end do
 
      if (na) 116,116,102
  114 do j=1,m
         y(j) = w(j)
      end do
 
      ibr = 1
      go to 102
  116 if (m1) 117,117,118
  117 if (m2) 114,114,123
  118 if (m2) 120,120,119
  119 if (abs(bm1(m1))-abs(bm2(m2))) 123,123,120
  120 if (ibr) 121,121,122
  121 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 114,122,122
  122 rt = rt-bm1(m1)
      m1 = m1-1
      go to 126
  123 if (ibr) 124,124,125
  124 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 114,125,125
  125 rt = rt-bm2(m2)
      m2 = m2-1
  126 do j=1,m
         y(j) = y(j)+rt*w(j)
      end do
      go to 102
  128 return
      end
      function psgf (x,iz,c,a,bh)

c*********************************************************************72
c
cc PSGF
c
      dimension       a(1)       ,c(1)       ,bh(1)
      fsg = 1.
      hsg = 1.
      do 101 j=1,iz
         dd = 1./(x-bh(j))
         fsg = fsg*a(j)*dd
         hsg = hsg*c(j)*dd
  101 continue
      if (mod(iz,2)) 103,102,103
  102 psgf = 1.-fsg-hsg
      return
  103 psgf = 1.+fsg+hsg
      return
      end
      SUBROUTINE RADB2 (IDO,L1,CC,CH,WA1)

c*********************************************************************72
c
cc RADB2 - backward Fourier transform, radix 2.
c
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(1)
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADB3 (IDO,L1,CC,CH,WA1,WA2)

c*********************************************************************72
c
cc RADB3 - backward Fourier transform, radix 3.
c
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(1)     ,WA2(1)
      DATA TAUR,TAUI /-.5,.866025403784439/
      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE RADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)

c*********************************************************************72
c
cc RADB4 - backward Fourier transform, radix 4.
c
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)
      DATA SQRT2 /1.414213562373095/
      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)

c*********************************************************************72
c
cc RADB5 - backward Fourier transform, radix 5.
c
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE RADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

c*********************************************************************72
c
cc RADBG - backward Fourier transform, general radix.
c
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(1)
      DATA TPI/6.28318530717959/
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
      SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)

c*********************************************************************72
c
cc RADF2 - forward Fourier transform, radix 2.
c
      DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           ,
     1                WA1(1)
      DO 101 K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADF3 (IDO,L1,CC,CH,WA1,WA2)

c*********************************************************************72
c
cc RADF3 - forward Fourier transform, radix 3.
c
      DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           ,
     1                WA1(1)     ,WA2(1)
      DATA TAUR,TAUI /-.5,.866025403784439/
      DO 101 K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)

c*********************************************************************72
c
cc RADF4 - forward Fourier transform, radix 4.
c
      DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)
      DATA HSQT2 /.7071067811865475/
      DO 101 K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(IDO,4,K) = TR2-TR1
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1+CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(IDO,K,3)
         CH(1,4,K) = TI1+CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)

c*********************************************************************72
c
cc RADF5 - forward Fourier transform, radix 5.
c
      DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      DO 101 K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

c*********************************************************************72
c
cc RADFG - forward Fourier transform, general radix.
c
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
      DATA TPI/6.28318530717959/
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,K,J) = C1(1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(I,1,K) = CH(I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(I,1,K) = CH(I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE RFFTB (N,R,WSAVE)

c*********************************************************************72
c
cc RFFTB - backward Fourier transform.
c
      DIMENSION       R(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE RFFTB1 (N,C,CH,WA,IFAC)

c*********************************************************************72
c
cc RFFTB1 is a utility routine for RFFTB
c
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL RADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE RFFTF (N,R,WSAVE)

c*********************************************************************72
c
cc RFFTF - forward Fourier transform.
c
      DIMENSION       R(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE RFFTF1 (N,C,CH,WA,IFAC)

c*********************************************************************72
c
cc RFFTF1 is a utility routine for RFFTF.
c
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL RADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END
      SUBROUTINE RFFTI (N,WSAVE)

c*********************************************************************72
c
cc RFFTI - initialized Fourier transform.
c
      DIMENSION       WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE RFFTI1 (N,WA,IFAC)

c*********************************************************************72
c
cc RFFTI1 is a utility routine for RFFTI.
c
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 6.28318530717959
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
      subroutine sepeli (intl,iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,
     &                   d,n,nbdcnd,bdc,gama,bdd,xnu,cofx,cofy,grhs,
     &                   usol,idmn,w,pertrb,ierror)

c*********************************************************************72
c
cc SEPELI 2D general separable elliptic problem, second or fourth order scheme.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c
c dimension of           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
c arguments              usol(idmn,n+1), grhs(idmn,n+1),
c                        w (see argument list)
c
c latest revision        march 1977
c
c purpose                sepeli solves for either the second-order
c                        finite difference approximation or a
c                        fourth-order approximation to a separable
c                        elliptic equation
c
c                                    2    2
c                             af(x)*d u/dx + bf(x)*du/dx  + cf(x)*u +
c                                    2    2
c                             df(y)*d u/dy  + ef(y)*du/dy + ff(y)*u
c
c                             = g(x,y)
c
c                        on a rectangle (x greater than or equal to a
c                        and less than or equal to b; y greater than
c                        or equal to c and less than or equal to d).
c                        any combination of periodic or mixed boundary
c                        conditions is allowed.
c
c purpose                the possible boundary conditions are:
c                        in the x-direction:
c                         (0) periodic, u(x+b-a,y)=u(x,y) for all y,x
c                         (1) u(a,y), u(b,y) are specified for all y
c                         (2) u(a,y), du(b,y)/dx+beta*u(b,y) are
c                             specified for all y
c                         (3) du(a,y)/dx+alpha*u(a,y),du(b,y)/dx+
c                             beta*u(b,y) are specified for all y
c                         (4) du(a,y)/dx+alpha*u(a,y),u(b,y) are
c                             specified for all y
c
c                        in the y-direction:
c                         (0) periodic, u(x,y+d-c)=u(x,y) for all x,y
c                         (1) u(x,c),u(x,d) are specified for all x
c                         (2) u(x,c),du(x,d)/dy+xnu*u(x,d) are specified
c                             for all x
c                         (3) du(x,c)/dy+gama*u(x,c),du(x,d)/dy+
c                             xnu*u(x,d) are specified for all x
c                         (4) du(x,c)/dy+gama*u(x,c),u(x,d) are
c                             specified for all x
c
c arguments
c
c on input               intl
c                          = 0 on initial entry to sepeli or if any of
c                              the arguments c, d, n, nbdcnd, cofy are
c                              changed from a previous call
c                          = 1 if c, d, n, nbdcnd, cofy are unchanged
c                              from the previous call.
c
c                        iorder
c                          = 2 if a second-order approximation is sought
c                          = 4 if a fourth-order approximation is sought
c
c                        a,b
c                          the range of the x-independent variable;
c                          i.e., x is greater than or equal to a and
c                          less than or equal to b.  a must be less than
c                          b.
c
c                        m
c                          the number of panels into which the interval
c                          [a,b] is subdivided.  hence, there will be
c                          m+1 grid points in the x-direction given by
c                          xi=a+(i-1)*dlx for i=1,2,...,m+1 where
c                          dlx=(b-a)/m is the panel width.  m must be
c                          less than idmn and greater than 5.
c
c                        mbdcnd
c                          indicates the type of boundary condition at
c                          x=a and x=b
c                          = 0 if the solution is periodic in x; i.e.,
c                              u(x+b-a,y)=u(x,y) for all y,x
c                          = 1 if the solution is specified at x=a and
c                              x=b; i.e., u(a,y) and u(b,y) are
c                              specified for all y
c                          = 2 if the solution is specified at x=a and
c                              the boundary condition is mixed at x=b;
c                              i.e., u(a,y) and du(b,y)/dx+beta*u(b,y)
c                              are specified for all y
c                          = 3 if the boundary conditions at x=a and x=b
c                              are mixed; i.e., du(a,y)/dx+alpha*u(a,y)
c                              and du(b,y)/dx+beta*u(b,y) are specified
c                              for all y
c                          = 4 if the boundary condition at x=a is mixed
c                              and the solution is specified at x=b;
c                              i.e., du(a,y)/dx+alpha*u(a,y) and u(b,y)
c                              are specified for all y
c
c                        bda
c                          a one-dimensional array of length n+1 that
c                          specifies the values of du(a,y)/dx+
c                          alpha*u(a,y) at x=a, when mbdcnd=3 or 4.
c                               bda(j) = du(a,yj)/dx+alpha*u(a,yj);
c                               j=1,2,...,n+1
c                          when mbdcnd has any other value, bda is a
c                          dummy parameter.
c
c on input               alpha
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at x=a (see
c                          argument bda).  if mbdcnd ' 3,4 then alpha is
c                          a dummy parameter.
c
c                        bdb
c                          a one-dimensional array of length n+1 that
c                          specifies the values of du(b,y)/dx+
c                          beta*u(b,y) at x=b.  when mbdcnd=2 or 3
c                               bdb(j) = du(b,yj)/dx+beta*u(b,yj);
c                               j=1,2,...,n+1
c                          when mbdcnd has any other value, bdb is a
c                          dummy parameter.
c
c                        beta
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at x=b (see
c                          argument bdb).  if mbdcnd'2,3 then beta is a
c                          dummy parameter.
c
c                        c,d
c                          the range of the y-independent variable;
c                          i.e., y is greater than or equal to c and
c                          less than or equal to d.  c must be less than
c                          d.
c
c                        n
c                          the number of panels into which the interval
c                          [c,d] is subdivided.  hence, there will be
c                          n+1 grid points in the y-direction given by
c                          yj=c+(j-1)*dly for j=1,2,...,n+1 where
c                          dly=(d-c)/n is the panel width.  in addition,
c                          n must be greater than 4.
c
c                        nbdcnd
c                          indicates the types of boundary conditions at
c                          y=c and y=d
c                          = 0 if the solution is periodic in y; i.e.,
c                              u(x,y+d-c)=u(x,y) for all x,y
c                          = 1 if the solution is specified at y=c and
c                              y = d, i.e., u(x,c) and u(x,d) are
c                              specified for all x
c                          = 2 if the solution is specified at y=c and
c                              the boundary condition is mixed at y=d;
c                              i.e., u(x,c) and du(x,d)/dy+xnu*u(x,d)
c                              are specified for all x
c                          = 3 if the boundary conditions are mixed at
c                              y=c and y=d; i.e., du(x,d)/dy+gama*u(x,c)
c                              and du(x,d)/dy+xnu*u(x,d) are specified
c                              for all x
c                          = 4 if the boundary condition is mixed at y=c
c                              and the solution is specified at y=d;
c                              i.e. du(x,c)/dy+gama*u(x,c) and u(x,d)
c                              are specified for all x
c
c                        bdc
c                          a one-dimensional array of length m+1 that
c                          specifies the value of du(x,c)/dy+gama*u(x,c)
c                          at y=c.  when nbdcnd=3 or 4
c                             bdc(i) = du(xi,c)/dy + gama*u(xi,c);
c                             i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdc is a
c                          dummy parameter.
c
c                        gama
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at y=c (see
c                          argument bdc).  if nbdcnd'3,4 then gama is a
c                          dummy parameter.
c
c                        bdd
c                          a one-dimensional array of length m+1 that
c                          specifies the value of du(x,d)/dy +
c                          xnu*u(x,d) at y=c.  when nbdcnd=2 or 3
c                            bdd(i) = du(xi,d)/dy + xnu*u(xi,d);
c                            i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdd is a
c                          dummy parameter.
c
c                        xnu
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at y=d (see
c                          argument bdd).  if nbdcnd'2 or 3 then xnu is
c                          a dummy parameter.
c
c                        cofx
c                          a user-supplied subprogram with
c                          parameters x, afun, bfun, cfun which
c                          returns the values of the x-dependent
c                          coefficients af(x), bf(x), cf(x) in
c                          the elliptic equation at x.
c
c                        cofy
c                          a user-supplied subprogram with
c                          parameters y, dfun, efun, ffun which
c                          returns the values of the y-dependent
c                          coefficients df(y), ef(y), ff(y) in
c                          the elliptic equation at y.
c
c                        note:  cofx and cofy must be declared external
c                        in the calling routine.  the values returned in
c                        afun and dfun must satisfy afun*dfun greater
c                        than 0 for a less than x less than b,
c                        c less than y less than d (see ierror=10).
c                        the coefficients provided may lead to a matrix
c                        equation which is not diagonally dominant in
c                        which case solution may fail (see ierror=4).
c
c                        grhs
c                          a two-dimensional array that specifies the
c                          values of the right-hand side of the elliptic
c                          equation; i.e., grhs(i,j)=g(xi,yi), for
c                          i=2,...,m; j=2,...,n.  at the boundaries,
c                          grhs is defined by
c
c                          mbdcnd   grhs(1,j)   grhs(m+1,j)
c                          ------   ---------   -----------
c                            0      g(a,yj)     g(b,yj)
c                            1         *           *
c                            2         *        g(b,yj)  j=1,2,...,n+1
c                            3      g(a,yj)     g(b,yj)
c                            4      g(a,yj)        *
c
c                          nbdcnd   grhs(i,1)   grhs(i,n+1)
c                          ------   ---------   -----------
c                            0      g(xi,c)     g(xi,d)
c                            1         *           *
c                            2         *        g(xi,d)  i=1,2,...,m+1
c                            3      g(xi,c)     g(xi,d)
c                            4      g(xi,c)        *
c
c                          where * means these quantites are not used.
c                          grhs should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        usol
c                          a two-dimensional array that specifies the
c                          values of the solution along the boundaries.
c                          at the boundaries, usol is defined by
c
c                          mbdcnd   usol(1,j)   usol(m+1,j)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(a,yj)     u(b,yj)
c                            2      u(a,yj)        *     j=1,2,...,n+1
c                            3         *           *
c                            4         *        u(b,yj)
c
c                          nbdcnd   usol(i,1)   usol(i,n+1)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(xi,c)     u(xi,d)
c                            2      u(xi,c)        *     i=1,2,...,m+1
c                            3         *           *
c                            4         *        u(xi,d)
c
c                          where * means the quantites are not used in
c                          the solution.
c
c                          if iorder=2, the user may equivalence grhs
c                          and usol to save space.  note that in this
c                          case the tables specifying the boundaries of
c                          the grhs and usol arrays determine the
c                          boundaries uniquely except at the corners.
c                          if the tables call for both g(x,y) and
c                          u(x,y) at a corner then the solution must be
c                          chosen.  for example, if mbdcnd=2 and
c                          nbdcnd=4, then u(a,c), u(a,d), u(b,d) must be
c                          chosen at the corners in addition to g(b,c).
c
c                          if iorder=4, then the two arrays, usol and
c                          grhs, must be distinct.
c
c                          usol should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        idmn
c                          the row (or first) dimension of the arrays
c                          grhs and usol as it appears in the program
c                          calling sepeli.  this parameter is used to
c                          specify the variable dimension of grhs and
c                          usol.  idmn must be at least 7 and greater
c                          than or equal to m+1.
c
c                        w
c                          a one-dimensional array that must be provided
c                          by the user for work space.  let
c                          k=int(log2(n+1))+1 and set  l=2**(k+1).
c                          then (k-2)*l+k+10*n+12*m+27 will suffice
c                          as a length of w.  the actual length of w in
c                          the calling routine must be set in w(1) (see
c                          ierror=11).
c
c on output              usol
c                          contains the approximate solution to the
c                          elliptic equation.  usol(i,j) is the
c                          approximation to u(xi,yj) for i=1,2...,m+1
c                          and j=1,2,...,n+1.  the approximation has
c                          error o(dlx**2+dly**2) if called with
c                          iorder=2 and o(dlx**4+dly**4) if called with
c                          iorder=4.
c
c                        w
c                          contains intermediate values that must not be
c                          destroyed if sepeli is called again with
c                          intl=1.  in addition w(1) contains the exact
c                          minimal length (in floating point) required
c                          for the work space (see ierror=11).
c
c                        pertrb
c                          if a combination of periodic or derivative
c                          boundary conditions (i.e., alpha=beta=0 if
c                          mbdcnd=3; gama=xnu=0 if nbdcnd=3) is
c                          specified and if the coefficients of u(x,y)
c                          in the separable elliptic equation are zero
c                          (i.e., cf(x)=0 for x greater than or equal to
c                          a and less than or equal to b; ff(y)=0 for
c                          y greater than or equal to c and less than
c                          or equal to d) then a solution may not exist.
c                          pertrb is a constant calculated and
c                          subtracted from the right-hand side of the
c                          matrix equations generated by sepeli which
c                          insures that a solution exists.  sepeli then
c                          computes this solution which is a weighted
c                          minimal least squares solution to the
c                          original problem.
c
c                        ierror
c                          an error flag that indicates invalid input
c                          parameters or failure to find a solution
c                          = 0 no error
c                          = 1 if a greater than b or c greater than d
c                          = 2 if mbdcnd less than 0 or mbdcnd greater
c                              than 4
c                          = 3 if nbdcnd less than 0 or nbdcnd greater
c                              than 4
c                          = 4 if attempt to find a solution fails.
c                              (the linear system generated is not
c                              diagonally dominant.)
c                          = 5 if idmn is too small (see discussion of
c                              idmn)
c                          = 6 if m is too small or too large (see
c                              discussion of m)
c                          = 7 if n is too small (see discussion of n)
c                          = 8 if iorder is not 2 or 4
c                          = 9 if intl is not 0 or 1
c                          = 10 if afun*dfun less than or equal to 0 for
c                               some interior mesh point (xi,yj)
c                          = 11 if the work space length input in w(1)
c                               is less than the exact minimal work
c                               space length required output in w(1).
c
c                          note (concerning ierror=4):  for the
c                          coefficients input through cofx, cofy, the
c                          discretization may lead to a block
c                          tridiagonal linear system which is not
c                          diagonally dominant (for example, this
c                          happens if cfun=0 and bfun/(2.*dlx) greater
c                          than afun/dlx**2).  in this case solution may
c                          fail.  this cannot happen in the limit as
c                          dlx, dly approach zero.  hence, the condition
c                          may be remedied by taking larger values for m
c                          or n.
c
c entry points           sepeli, spelip, chkprm, chksng, orthog, minsol,
c                        trisp, defer, dx, dy, blktri, blktr1, indxb,
c                        indxa, indxc, prod, prodp, cprod, cprodp,
c                        ppadd, psgf, bsrh, ppsgf, ppspf, compb,
c                        trun1, stor1, tqlrat
c
c special conditions     none
c
c common blocks          splp, cblkt, valu1
c
c i/o                    none
c
c precision              single
c
c specialist             john c. adams, ncar, boulder, colorado  80307
c
c language               fortran
c
c history                developed at ncar during 1975-76.
c
c algorithm              sepeli automatically discretizes the separable
c                        elliptic equation which is then solved by a
c                        generalized cyclic reduction algorithm in 
c                        blktri.  the fourth-order solution
c                        is obtained using 'deferred corrections' which
c                        is described and referenced in sections,
c                        references and method.
c
c space required         14654 (octal) = 6572 (decimal)
c
c accuracy and timing    the following computational results were
c                        obtained by solving the sample problem at the
c                        end of this write-up on the control data 7600.
c                        the op count is proportional to m*n*log2(n).
c                        in contrast to the other routines in this
c                        chapter, accuracy is tested by computing and
c                        tabulating second- and fourth-order
c                        discretization errors.  below is a table
c                        containing computational results.  the times
c                        given do not include initialization (i.e.,
c                        times are for intl=1).  note that the
c                        fourth-order accuracy is not realized until the
c                        mesh is sufficiently refined.
c
c              second-order    fourth-order   second-order  fourth-order
c    m    n   execution time  execution time    error         error
c               (m sec)         (m sec)
c     6    6         6              14          6.8e-1        1.2e0
c    14   14        23              58          1.4e-1        1.8e-1
c    30   30       100             247          3.2e-2        9.7e-3
c    62   62       445           1,091          7.5e-3        3.0e-4
c   126  126     2,002           4,772          1.8e-3        3.5e-6
c
c portability            there are no mahcine-dependent constants.
c
c required resident      sqrtf, abs, cabs, logf
c routines
c
c references             keller, h.b., numerical methods for two-point
c                          boundary-value problems, blaisdel (1968),
c                          waltham, mass.
c
c                        swarztrauber, p., and r. sweet (1975):
c                          efficient fortran subprograms for the
c                          solution of elliptic partial differential
c                          equations.  ncar technical note
c                          ncar-tn/ia-109, pp. 135-137.
c
c
c
c
      dimension       grhs(idmn,1)           ,usol(idmn,1)
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     &                w(1)
      external        cofx       ,cofy
c
c     check input parameters
c
      call chkprm (intl,iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,cofy,
     &             idmn,ierror)
      if (ierror .ne. 0) return
c
c     compute minimum work space and check work space length input
c
      l = n+1
      if (nbdcnd .eq. 0) l = n
      logb2n = int(alog(float(l)+0.5)/alog(2.0))+1
      ll = 2**(logb2n+1)
      k = m+1
      l = n+1
      length = (logb2n-2)*ll+logb2n+max0(2*l,6*k)+5
      if (nbdcnd .eq. 0) length = length+2*l
      ierror = 11
      linput = int(w(1)+0.5)
      loutpt = length+6*(k+l)+1
      w(1) = float(loutpt)
      if (loutpt .gt. linput) return
      ierror = 0
c
c     set work space indices
c
      i1 = length+2
      i2 = i1+l
      i3 = i2+l
      i4 = i3+l
      i5 = i4+l
      i6 = i5+l
      i7 = i6+l
      i8 = i7+k
      i9 = i8+k
      i10 = i9+k
      i11 = i10+k
      i12 = i11+k
      i13 = 2
      call spelip (intl,iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n,
     &             nbdcnd,bdc,gama,bdd,xnu,cofx,cofy,w(i1),w(i2),w(i3),
     &             w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10),w(i11),
     &             w(i12),grhs,usol,idmn,w(i13),pertrb,ierror)
      return
      end
      subroutine sepx4(iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n,
     &nbdcnd,bdc,bdd,cofx,grhs,usol,idmn,w,pertrb,ierror)

c*********************************************************************72
c
cc SEPX4 2D restricted separable elliptic problem, second or fourth order scheme.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c
c
c dimension of           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
c arguments              usol(idmn,n+1), grhs(idmn,n+1),
c                        w (see argument list)
c
c latest revision        october 1980
c
c
c purpose                sepx4 solves for either the second-order
c                        finite difference approximation or a
c                        a fourth-order approximation  to the
c                        solution of a separable elliptic equation
c                             af(x)*uxx+bf(x)*ux+cf(x)*u+uyy = g(x,y)
c
c                        on a rectangle (x greater than or equal to a
c                        and less than or equal to b; y greater than
c                        or equal to c and less than or equal to d).
c                        any combination of periodic or mixed boundary
c                        conditions is allowed.
c                        if boundary conditions in the x direction
c                        are periodic (see mbdcnd=0 below) then the
c                        coefficients must satisfy
c                        af(x)=c1,bf(x)=0,cf(x)=c2 for all x.
c                        here c1,c2 are constants, c1.gt.0.
c
c                        the possible boundary conditions are
c                        in the x-direction:
c                         (0) periodic, u(x+b-a,y)=u(x,y) for all y,x
c                         (1) u(a,y), u(b,y) are specified for all y
c                         (2) u(a,y), du(b,y)/dx+beta*u(b,y) are
c                             specified for all y
c                         (3) du(a,y)/dx+alpha*u(a,y),du(b,y)/dx+
c                             beta*u(b,y) are specified for all y
c                         (4) du(a,y)/dx+alpha*u(a,y),u(b,y) are
c                             specified for all y
c
c                        in the y-direction:
c                         (0) periodic, u(x,y+d-c)=u(x,y) for all x,y
c                         (1) u(x,c),u(x,d) are specified for all x
c                         (2) u(x,c),du(x,d)/dy are specified for all x
c                         (3) du(x,c)/dy,du(x,d)/dy are specified for
c                            all x
c                        (4) du(x,c)/dy,u(x,d) are specified for all x
c
c usage                  call sepx4(iorder,a,b,m,mbdcnd,bda,alpha,bdb,
c                                  beta,c,d,n,nbdcnd,bdc,bdd,cofx,
c                                  grhs,usol,idmn,w,pertrb,ierror)
c
c arguments
c
c
c                        iorder
c                          = 2 if a second-order approximation is sought
c                          = 4 if a fourth-order approximation is sought
c
c                        a,b
c                          the range of the x-independent variable;
c                          i.e., x is greater than or equal to a and
c                          less than or equal to b.  a must be less than
c                          b.
c
c                        m
c                          the number of panels into which the interval
c                          [a,b] is subdivided.  hence, there will be
c                          m+1 grid points in the x-direction given by
c                          xi=a+(i-1)*dlx for i=1,2,...,m+1 where
c                          dlx=(b-a)/m is the panel width.  m must be
c                          less than idmn and greater than 5.
c
c                        mbdcnd
c                          indicates the type of boundary condition at
c                          x=a and x=b
c                          = 0 if the solution is periodic in x; i.e.,
c                              u(x+b-a,y)=u(x,y) for all y,x
c                          = 1 if the solution is specified at x=a and
c                              x=b; i.e., u(a,y) and u(b,y) are
c                              specified for all y
c                          = 2 if the solution is specified at x=a and
c                              the boundary condition is mixed at x=b;
c                              i.e., u(a,y) and du(b,y)/dx+beta*u(b,y)
c                              are specified for all y
c                          = 3 if the boundary conditions at x=a and x=b
c                              are mixed; i.e., du(a,y)/dx+alpha*u(a,y)
c                              and du(b,y)/dx+beta*u(b,y) are specified
c                              for all y
c                          = 4 if the boundary condition at x=a is mixed
c                              and the solution is specified at x=b;
c                              i.e., du(a,y)/dx+alpha*u(a,y) and u(b,y)
c                              are specified for all y
c
c                        bda
c                          a one-dimensional array of length n+1 that
c                          specifies the values of du(a,y)/dx+
c                          alpha*u(a,y) at x=a, when mbdcnd=3 or 4.
c                               bda(j) = du(a,yj)/dx+alpha*u(a,yj);
c                               j=1,2,...,n+1
c                          when mbdcnd has any other value, bda is a
c                          dummy parameter.
c
c on input               alpha
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at x=a (see
c                          argument bda).  if mbdcnd ' 3,4 then alpha is
c                          a dummy parameter.
c
c                        bdb
c                          a one-dimensional array of length n+1 that
c                          specifies the values of du(b,y)/dx+
c                          beta*u(b,y) at x=b.  when mbdcnd=2 or 3
c                               bdb(j) = du(b,yj)/dx+beta*u(b,yj);
c                               j=1,2,...,n+1
c                          when mbdcnd has any other value, bdb is a
c                          dummy parameter.
c
c                        beta
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at x=b (see
c                          argument bdb).  if mbdcnd'2,3 then beta is a
c                          dummy parameter.
c
c                        c,d
c                          the range of the y-independent variable;
c                          i.e., y is greater than or equal to c and
c                          less than or equal to d.  c must be less than
c                          d.
c
c                        n
c                          the number of panels into which the interval
c                          [c,d] is subdivided.  hence, there will be
c                          n+1 grid points in the y-direction given by
c                          yj=c+(j-1)*dly for j=1,2,...,n+1 where
c                          dly=(d-c)/n is the panel width.  in addition,
c                          n must be greater than 4.
c
c                        nbdcnd
c                          indicates the types of boundary conditions at
c                          y=c and y=d
c                          = 0 if the solution is periodic in y; i.e.,
c                              u(x,y+d-c)=u(x,y) for all x,y
c                          = 1 if the solution is specified at y=c and
c                              y = d, i.e., u(x,c) and u(x,d) are
c                              specified for all x
c                          = 2 if the solution is specified at y=c and
c                              the boundary condition is mixed at y=d;
c                              i.e., du(x,c)/dy and u(x,d)
c                              are specified for all x
c                          = 3 if the boundary conditions are mixed at
c                              y=cand y=d i.e., du(x,d)/dy
c                              and du(x,d)/dy are specified
c                              for all x
c                          = 4 if the boundary condition is mixed at y=c
c                              and the solution is specified at y=d;
c                              i.e. du(x,c)/dy+gama*u(x,c) and u(x,d)
c                              are specified for all x
c
c                        bdc
c                          a one-dimensional array of length m+1 that
c                          specifies the value du(x,c)/dy
c                          at y=c.  when nbdcnd=3 or 4
c                            bdc(i) = du(xi,c)/dy
c                             i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdc is a
c                          dummy parameter.
c
c
c                        bdd
c                          a one-dimensional array of length m+1 that
c                           specified the value of du(x,d)/dy
c                           at y=d.  when nbdcnd=2 or 3
c                             bdd(i)=du(xi,d)/dy
c                            i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdd is a
c                          dummy parameter.
c
c
c                        cofx
c                          a user-supplied subprogram with
c                          parameters x, afun, bfun, cfun which
c                          returns the values of the x-dependent
c                          coefficients af(x), bf(x), cf(x) in
c                          the elliptic equation at x.
c                        if boundary conditions in the x direction
c                        are periodic then the coefficients
c                        must satisfy af(x)=c1,bf(x)=0,cf(x)=c2 for
c                        all x.  here c1.gt.0 and c2 are constants.
c
c                        note that cofx must be declared external
c                        in the calling routine.
c
c                        grhs
c                          a two-dimensional array that specifies the
c                          values of the right-hand side of the elliptic
c                          equation; i.e., grhs(i,j)=g(xi,yi), for
c                          i=2,...,m; j=2,...,n.  at the boundaries,
c                          grhs is defined by
c
c                          mbdcnd   grhs(1,j)   grhs(m+1,j)
c                          ------   ---------   -----------
c                            0      g(a,yj)     g(b,yj)
c                            1         *           *
c                            2         *        g(b,yj)  j=1,2,...,n+1
c                            3      g(a,yj)     g(b,yj)
c                            4      g(a,yj)        *
c
c                          nbdcnd   grhs(i,1)   grhs(i,n+1)
c                          ------   ---------   -----------
c                            0      g(xi,c)     g(xi,d)
c                            1         *           *
c                            2         *        g(xi,d)  i=1,2,...,m+1
c                            3      g(xi,c)     g(xi,d)
c                            4      g(xi,c)        *
c
c                          where * means these quantites are not used.
c                          grhs should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        usol
c                          a two-dimensional array that specifies the
c                          values of the solution along the boundaries.
c                          at the boundaries, usol is defined by
c
c                          mbdcnd   usol(1,j)   usol(m+1,j)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(a,yj)     u(b,yj)
c                            2      u(a,yj)        *     j=1,2,...,n+1
c                            3         *           *
c                            4         *        u(b,yj)
c
c                          nbdcnd   usol(i,1)   usol(i,n+1)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(xi,c)     u(xi,d)
c                            2      u(xi,c)        *     i=1,2,...,m+1
c                            3         *           *
c                            4         *        u(xi,d)
c
c                          where * means the quantites are not used in
c                          the solution.
c
c                          if iorder=2, the user may equivalence grhs
c                          and usol to save space.  note that in this
c                          case the tables specifying the boundaries of
c                          the grhs and usol arrays determine the
c                          boundaries uniquely except at the corners.
c                          if the tables call for both g(x,y) and
c                          u(x,y) at a corner then the solution must be
c                          chosen.  for example, if mbdcnd=2 and
c                          nbdcnd=4, then u(a,c), u(a,d), u(b,d) must be
c                          chosen at the corners in addition to g(b,c).
c
c                          if iorder=4, then the two arrays, usol and
c                          grhs, must be distinct.
c
c                          usol should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        idmn
c                          the row (or first) dimension of the arrays
c                          grhs and usol as it appears in the program
c                          calling sepx4.  this parameter is used to
c                          specify the variable dimension of grhs and
c                          usol.  idmn must be at least 7 and greater
c                          than or equal to m+1.
c
c                        w
c                          a one-dimensional array that must be provided
c                        by the user for work space.
c                        10*n+(16+int(log2(n)))*(m+1)+23 will suffice
c                        as a length for w.  the actual length of
c                        w in the calling routine must be set in w(1)
c                        (see ierror=11).
c on output              usol
c                          contains the approximate solution to the
c                          elliptic equation.  usol(i,j) is the
c                          approximation to u(xi,yj) for i=1,2...,m+1
c                          and j=1,2,...,n+1.  the approximation has
c                          error o(dlx**2+dly**2) if called with
c                          iorder=2 and o(dlx**4+dly**4) if called with
c                          iorder=4.
c
c                        w
c                          contains intermediate values that must not be
c                          destroyed if sepx4 is called again with
c                          intl=1.  in addition w(1) contains the exact
c                          minimal length (in floating point) required
c                          for the work space (see ierror=11).
c
c                        pertrb
c                        if a combination of periodic or derivative
c                        boundary conditions (i.e., alpha=beta=0 if
c                        mbdcnd=3) is specified and if cf(x)=0 for all x
c                        then a solution to the discretized matrix
c                        equation may not exist (reflecting the non-
c                        uniqueness of solutions to the pde).  pertrb
c                        is a constant calculated and subtracted from
c                        the right hand side of the matrix equation
c                        insuring the existence of a solution.
c                        sepx4 computes this solution which is a
c                        weighted minimal least squares solution to
c                        the original problem.  if singularity is
c                        not detected pertrb=0.0 is returned by
c                        sepx4.
c
c                        ierror
c                          an error flag that indicates invalid input
c                          parameters or failure to find a solution
c                          = 0 no error
c                          = 1 if a greater than b or c greater than d
c                          = 2 if mbdcnd less than 0 or mbdcnd greater
c                              than 4
c                          = 3 if nbdcnd less than 0 or nbdcnd greater
c                              than 4
c                          = 4 if attempt to find a solution fails.
c                              (the linear system generated is not
c                              diagonally dominant.)
c                          = 5 if idmn is too small (see discussion of
c                              idmn)
c                          = 6 if m is too small or too large (see
c                              discussion of m)
c                          = 7 if n is too small (see discussion of n)
c                          = 8 if iorder is not 2 or 4
c                          = 9 if intl is not 0 or 1
c                          = 10 if afun is less than or equal to zero
c                          for some interior mesh point xi
c                               some interior mesh point (xi,yj)
c                          = 11 if the work space length input in w(1)
c                               is less than the exact minimal work
c                               space length required output in w(1).
c                        = 12 if mbdcnd=0 and af(x)=cf(x)=constant
c                             or bf(x)=0 for all x is not true.
c
c
c special conditions     none
c
c common blocks          spl4
c
c i/o                    none
c
c precision              single
c
c required library       none
c files
c
c specialist             john c. adams, ncar, boulder, colorado  80307
c
c language               fortran
c
c
c  entry points          sepx4,speli4,chkpr4,chksn4,ortho4,minso4,tris4,
c                        defe4,dx4,dy4
c
c history                sepx4 was developed by modifying the ulib
c                        routine sepeli during october 1978.
c                        it should be used instead of sepeli whenever
c                        possible.  the increase in speed is at least
c                        a factor of three.
c
c algorithm              sepx4 automatically discretizes the separable
c                        elliptic equation which is then solved by a
c                        generalized cyclic reduction algorithm in the
c                        pois.  the fourth order solution
c                        is obtained using the technique of
c                        defferred corrections referenced below.
c
c
c references             keller, h.b., numerical methods for two-point
c                          boundary-value problems, blaisdel (1968),
c                          waltham, mass.
c
c                        swarztrauber, p., and r. sweet (1975):
c                          efficient fortran subprograms for the
c                          solution of elliptic partial differential
c                          equations.  ncar technical note
c                          ncar-tn/ia-109, pp. 135-137.
c
c
c
c
      dimension       grhs(idmn,1)           ,usol(idmn,1)
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     &                w(1)
      external cofx
c
c     check input parameters
c
      call chkpr4(iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,idmn,ierror)
      if (ierror .ne. 0) return
c
c     compute minimum work space and check work space length input
c
      l = n+1
      if (nbdcnd .eq. 0) l = n
      k = m+1
      l = n+1
c     estimate log base 2 of n
      log2n=int(alog(float(n+1))/alog(2.0)+0.5)
      length=4*(n+1)+(10+log2n)*(m+1)
      ierror = 11
      linput = int(w(1)+0.5)
      loutpt = length+6*(k+l)+1
      w(1) = float(loutpt)
      if (loutpt .gt. linput) return
      ierror = 0
c
c     set work space indices
c
      i1 = length+2
      i2 = i1+l
      i3 = i2+l
      i4 = i3+l
      i5 = i4+l
      i6 = i5+l
      i7 = i6+l
      i8 = i7+k
      i9 = i8+k
      i10 = i9+k
      i11 = i10+k
      i12 = i11+k
      i13 = 2
      call speli4(iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n,
     &nbdcnd,bdc,bdd,cofx,w(i1),w(i2),w(i3),
     &             w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10),w(i11),
     &             w(i12),grhs,usol,idmn,w(i13),pertrb,ierror)
      return
      end
      SUBROUTINE SINQB (N,X,WSAVE)

c*********************************************************************72
c
cc SINQB backward sine quarter wave transform.
c
      DIMENSION       X(*)       ,WSAVE(*)
      IF (N .GT. 1) GO TO 101
      X(1) = 4.*X(1)
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      CALL COSQB (N,X,WSAVE)
      DO 103 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  103 CONTINUE
      RETURN
      END
      SUBROUTINE SINQF (N,X,WSAVE)

c*********************************************************************72
c
cc SINQF forward sine quarter wave transform.
c
      DIMENSION       X(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      NS2 = N/2
      DO 101 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  101 CONTINUE
      CALL COSQF (N,X,WSAVE)
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      RETURN
      END
      SUBROUTINE SINQI (N,WSAVE)

c*********************************************************************72
c
cc SINQI initializes the sine quarter wave transform.
c
      DIMENSION       WSAVE(*)
      CALL COSQI (N,WSAVE)
      RETURN
      END
      SUBROUTINE SINT (N,X,WSAVE)

c*********************************************************************72
c
cc SINT: the sine transform.
c
      DIMENSION       X(*)       ,WSAVE(*)
      NP1 = N+1
      IW1 = N/2+1
      IW2 = IW1+NP1
      IW3 = IW2+NP1
      CALL SINT1(N,X,WSAVE,WSAVE(IW1),WSAVE(IW2),WSAVE(IW3))
      RETURN
      END
      SUBROUTINE SINT1(N,WAR,WAS,XH,X,IFAC)
      DIMENSION WAR(*),WAS(*),X(*),XH(*),IFAC(*)
      DATA SQRT3 /1.73205080756888/
      DO 100 I=1,N
      XH(I) = WAR(I)
      WAR(I) = X(I)
  100 CONTINUE
      IF (N-2) 101,102,103
  101 XH(1) = XH(1)+XH(1)
      GO TO 106
  102 XHOLD = SQRT3*(XH(1)+XH(2))
      XH(2) = SQRT3*(XH(1)-XH(2))
      XH(1) = XHOLD
      GO TO 106
  103 NP1 = N+1
      NS2 = N/2
      X(1) = 0.
      DO 104 K=1,NS2
         KC = NP1-K
         T1 = XH(K)-XH(KC)
         T2 = WAS(K)*(XH(K)+XH(KC))
         X(K+1) = T1+T2
         X(KC+1) = T2-T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) X(NS2+2) = 4.*XH(NS2+1)
      CALL RFFTF1 (NP1,X,XH,WAR,IFAC)
      XH(1) = .5*X(1)
      DO 105 I=3,N,2
         XH(I-1) = -X(I)
         XH(I) = XH(I-2)+X(I-1)
  105 CONTINUE
      IF (MODN .NE. 0) GO TO 106
      XH(N) = -X(N+1)
  106 DO 107 I=1,N
      X(I) = WAR(I)
      WAR(I) = XH(I)
  107 CONTINUE
      RETURN
      END
      SUBROUTINE SINTI (N,WSAVE)

c*********************************************************************72
c
cc SINTI initializes the sine transform.
c
      DIMENSION       WSAVE(*)
      DATA PI /3.14159265358979/
      IF (N .LE. 1) RETURN
      NS2 = N/2
      NP1 = N+1
      DT = PI/FLOAT(NP1)
      DO 101 K=1,NS2
         WSAVE(K) = 2.*SIN(K*DT)
  101 CONTINUE
      CALL RFFTI (NP1,WSAVE(NS2+1))
      RETURN
      END
      subroutine speli4(iorder,a,b,m,mbdcnd,
     &bda,alpha,bdb,beta,c,d,n,nbdcnd,bdc,bdd,cofx,an,bn,
     &                   cn,dn,un,zn,am,bm,cm,dm,um,zm,grhs,usol,idmn,
     &                   w,pertrb,ierror)

c*********************************************************************72
c
cc SPELI4 sets up vectors and arrays for input to BLKTRI
c     and computes a second order solution in usol.  a return jump to
c     sepx4 occurrs if iorder=2.  if iorder=4 a fourth order
c     solution is generated in usol.
c
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     &                w(1)
      dimension       grhs(idmn,1)           ,usol(idmn,1)
      dimension       an(1)      ,bn(1)      ,cn(1)      ,dn(1)      ,
     &                un(1)      ,zn(1)
      dimension       am(1)      ,bm(1)      ,cm(1)      ,dm(1)      ,
     &                um(1)      ,zm(1)
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: spl4
!$omp threadprivate (spl4)
      logical         singlr
      external cofx
c
c
c     set parameters internally
c
      kswx = mbdcnd+1
      kswy = nbdcnd+1
      k = m+1
      l = n+1
      ait = a
      bit = b
      cit = c
      dit = d
      dly=(dit-cit)/float(n)
c
c     set right hand side values from grhs in usol on the interior
c     and non-specified boundaries.
c
      do  20 i=2,m
         do  10 j=2,n
      usol(i,j)=dly**2*grhs(i,j)
   10    continue
   20 continue
      if (kswx.eq.2 .or. kswx.eq.3) go to  40
      do  30 j=2,n
      usol(1,j)=dly**2*grhs(1,j)
   30 continue
   40 continue
      if (kswx.eq.2 .or. kswx.eq.5) go to  60
      do  50 j=2,n
      usol(k,j)=dly**2*grhs(k,j)
   50 continue
   60 continue
      if (kswy.eq.2 .or. kswy.eq.3) go to  80
      do  70 i=2,m
      usol(i,1)=dly**2*grhs(i,1)
   70 continue
   80 continue
      if (kswy.eq.2 .or. kswy.eq.5) go to 100
      do  90 i=2,m
      usol(i,l)=dly**2*grhs(i,l)
   90 continue
  100 continue
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.3)
     &usol(1,1)=dly**2*grhs(1,1)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.3)
     &usol(k,1)=dly**2*grhs(k,1)
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.5)
     &usol(1,l)=dly**2*grhs(1,l)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.5)
     &usol(k,l)=dly**2*grhs(k,l)
      i1 = 1
c
c     set switches for periodic or non-periodic boundaries
c
      mp=1
      if(kswx.eq.1) mp=0
      np=nbdcnd
c
c     set dlx,dly and size of block tri-diagonal system generated
c     in nint,mint
c
      dlx = (bit-ait)/float(m)
      mit = k-1
      if (kswx .eq. 2) mit = k-2
      if (kswx .eq. 4) mit = k
      dly = (dit-cit)/float(n)
      nit = l-1
      if (kswy .eq. 2) nit = l-2
      if (kswy .eq. 4) nit = l
      tdlx3 = 2.0*dlx**3
      dlx4 = dlx**4
      tdly3 = 2.0*dly**3
      dly4 = dly**4
c
c     set subscript limits for portion of array to input to BLKTRI.
c
      is = 1
      js = 1
      if (kswx.eq.2 .or. kswx.eq.3) is = 2
      if (kswy.eq.2 .or. kswy.eq.3) js = 2
      ns = nit+js-1
      ms = mit+is-1
c
c     set x - direction
c
      do 110 i=1,mit
         xi = ait+float(is+i-2)*dlx
         call cofx (xi,ai,bi,ci)
         axi = (ai/dlx-0.5*bi)/dlx
         bxi = -2.*ai/dlx**2+ci
         cxi = (ai/dlx+0.5*bi)/dlx
      am(i)=dly**2*axi
      bm(i)=dly**2*bxi
      cm(i)=dly**2*cxi
  110 continue
c
c     set y direction
c
      do 120 j=1,nit
      dyj=1.0
      eyj=-2.0
      fyj=1.0
         an(j) = dyj
         bn(j) = eyj
         cn(j) = fyj
  120 continue
c
c     adjust edges in x direction unless periodic
c
      ax1 = am(1)
      cxm = cm(mit)
      go to (170,130,150,160,140),kswx
c
c     dirichlet-dirichlet in x direction
c
  130 am(1) = 0.0
      cm(mit) = 0.0
      go to 170
c
c     mixed-dirichlet in x direction
c
  140 am(1) = 0.0
      bm(1) = bm(1)+2.*alpha*dlx*ax1
      cm(1) = cm(1)+ax1
      cm(mit) = 0.0
      go to 170
c
c     dirichlet-mixed in x direction
c
  150 am(1) = 0.0
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*beta*dlx*cxm
      cm(mit) = 0.0
      go to 170
c
c     mixed - mixed in x direction
c
  160 continue
      am(1) = 0.0
      bm(1) = bm(1)+2.*dlx*alpha*ax1
      cm(1) = cm(1)+ax1
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*dlx*beta*cxm
      cm(mit) = 0.0
  170 continue
c
c     adjust in y direction unless periodic
c
      dy1 = an(1)
      fyn = cn(nit)
      gama=0.0
      xnu=0.0
      go to (220,180,200,210,190),kswy
c
c     dirichlet-dirichlet in y direction
c
  180 continue
      an(1) = 0.0
      cn(nit) = 0.0
      go to 220
c
c     mixed-dirichlet in y direction
c
  190 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      cn(nit) = 0.0
      go to 220
c
c     dirichlet-mixed in y direction
c
  200 an(1) = 0.0
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.*dly*xnu*fyn
      cn(nit) = 0.0
      go to 220
c
c     mixed - mixed direction in y direction
c
  210 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.0*dly*xnu*fyn
      cn(nit) = 0.0
  220 if (kswx .eq. 1) go to 270
c
c     adjust usol along x edge
c
      do 260 j=js,ns
         if (kswx.ne.2 .and. kswx.ne.3) go to 230
         usol(is,j) = usol(is,j)-ax1*usol(1,j)
         go to 240
  230    usol(is,j) = usol(is,j)+2.0*dlx*ax1*bda(j)
  240    if (kswx.ne.2 .and. kswx.ne.5) go to 250
         usol(ms,j) = usol(ms,j)-cxm*usol(k,j)
         go to 260
  250    usol(ms,j) = usol(ms,j)-2.0*dlx*cxm*bdb(j)
  260 continue
  270 if (kswy .eq. 1) go to 320
c
c     adjust usol along y edge
c
      do 310 i=is,ms
         if (kswy.ne.2 .and. kswy.ne.3) go to 280
         usol(i,js) = usol(i,js)-dy1*usol(i,1)
         go to 290
  280    usol(i,js) = usol(i,js)+2.0*dly*dy1*bdc(i)
  290    if (kswy.ne.2 .and. kswy.ne.5) go to 300
         usol(i,ns) = usol(i,ns)-fyn*usol(i,l)
         go to 310
  300    usol(i,ns) = usol(i,ns)-2.0*dly*fyn*bdd(i)
  310 continue
  320 continue
c
c     save adjusted edges in grhs if iorder=4
c
      if (iorder .ne. 4) go to 350
      do 330 j=js,ns
         grhs(is,j) = usol(is,j)
         grhs(ms,j) = usol(ms,j)
  330 continue
      do 340 i=is,ms
         grhs(i,js) = usol(i,js)
         grhs(i,ns) = usol(i,ns)
  340 continue
  350 continue
      iord = iorder
      pertrb = 0.0
c
c     check if operator is singular
c
      call chksn4(mbdcnd,nbdcnd,alpha,beta,cofx,singlr)
c
c     compute non-zero eigenvector in null space of transpose
c     if singular
c
      if (singlr) call tris4 (mit,am,bm,cm,dm,um,zm)
      if (singlr) call tris4 (nit,an,bn,cn,dn,un,zn)
c
c     adjust right hand side if necessary
c
  360 continue
      if (singlr) call ortho4 (usol,idmn,zn,zm,pertrb)
c
c     compute solution
c
c     save adjusted right hand side in grhs
      do 444 j=js,ns
      do 444 i=is,ms
      grhs(i,j)=usol(i,j)
  444 continue
      call genbun(np,nit,mp,mit,am,bm,cm,idmn,usol(is,js),ieror,w)
c     check if error detected in pois
c     this can only correspond to ierror=12
      if(ieror.eq.0) go to 224
c     set error flag if improper coefficients input to pois
      ierror=12
      return
  224 continue
      if (ierror .ne. 0) return
c
c     set periodic boundaries if necessary
c
      if (kswx .ne. 1) go to 380
      do 370 j=1,l
         usol(k,j) = usol(1,j)
  370 continue
  380 if (kswy .ne. 1) go to 400
      do 390 i=1,k
         usol(i,l) = usol(i,1)
  390 continue
  400 continue
c
c     minimize solution with respect to weighted least squares
c     norm if operator is singular
c
      if (singlr) call minso4 (usol,idmn,zn,zm,prtrb)
c
c     return if defe4red corrections and a fourth order solution are
c     not flagged
c
      if (iord .eq. 2) return
      iord = 2
c
c     compute new right hand side for fourth order solution
c
      call defe4(cofx,idmn,usol,grhs)
      go to 360
      end
      subroutine spelip (intl,iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,
     &                   d,n,nbdcnd,bdc,gama,bdd,xnu,cofx,cofy,an,bn,
     &                   cn,dn,un,zn,am,bm,cm,dm,um,zm,grhs,usol,idmn,
     &                   w,pertrb,ierror)

c*********************************************************************72
c
cc SPELIP sets up vectors and arrays for input to BLKTRI
c     and computes a second order solution in usol.  a return jump to
c     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
c     solution is generated in usol.
c
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     &                w(1)
      dimension       grhs(idmn,1)           ,usol(idmn,1)
      dimension       an(1)      ,bn(1)      ,cn(1)      ,dn(1)      ,
     &                un(1)      ,zn(1)
      dimension       am(1)      ,bm(1)      ,cm(1)      ,dm(1)      ,
     &                um(1)      ,zm(1)
      common /splp/   kswx       ,kswy       ,k          ,l          ,
     &                ait        ,bit        ,cit        ,dit        ,
     &                mit        ,nit        ,is         ,ms         ,
     &                js         ,ns         ,dlx        ,dly        ,
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
     save :: splp
!$omp threadprivate (splp)
      logical         singlr
      external        cofx       ,cofy
c
c     set parameters internally
c
      kswx = mbdcnd+1
      kswy = nbdcnd+1
      k = m+1
      l = n+1
      ait = a
      bit = b
      cit = c
      dit = d
c
c     set right hand side values from grhs in usol on the interior
c     and non-specified boundaries.
c
      do  20 i=2,m
         do  10 j=2,n
            usol(i,j) = grhs(i,j)
   10    continue
   20 continue
      if (kswx.eq.2 .or. kswx.eq.3) go to  40
      do  30 j=2,n
         usol(1,j) = grhs(1,j)
   30 continue
   40 continue
      if (kswx.eq.2 .or. kswx.eq.5) go to  60
      do  50 j=2,n
         usol(k,j) = grhs(k,j)
   50 continue
   60 continue
      if (kswy.eq.2 .or. kswy.eq.3) go to  80
      do  70 i=2,m
         usol(i,1) = grhs(i,1)
   70 continue
   80 continue
      if (kswy.eq.2 .or. kswy.eq.5) go to 100
      do  90 i=2,m
         usol(i,l) = grhs(i,l)
   90 continue
  100 continue
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.3)
     &    usol(1,1) = grhs(1,1)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.3)
     &    usol(k,1) = grhs(k,1)
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.5)
     &    usol(1,l) = grhs(1,l)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.5)
     &    usol(k,l) = grhs(k,l)
      i1 = 1
c
c     set switches for periodic or non-periodic boundaries
c
      mp = 1
      np = 1
      if (kswx .eq. 1) mp = 0
      if (kswy .eq. 1) np = 0
c
c     set dlx,dly and size of block tri-diagonal system generated
c     in nint,mint
c
      dlx = (bit-ait)/float(m)
      mit = k-1
      if (kswx .eq. 2) mit = k-2
      if (kswx .eq. 4) mit = k
      dly = (dit-cit)/float(n)
      nit = l-1
      if (kswy .eq. 2) nit = l-2
      if (kswy .eq. 4) nit = l
      tdlx3 = 2.0*dlx**3
      dlx4 = dlx**4
      tdly3 = 2.0*dly**3
      dly4 = dly**4
c
c     set subscript limits for portion of array to input to BLKTRI.
c
      is = 1
      js = 1
      if (kswx.eq.2 .or. kswx.eq.3) is = 2
      if (kswy.eq.2 .or. kswy.eq.3) js = 2
      ns = nit+js-1
      ms = mit+is-1
c
c     set x - direction
c
      do 110 i=1,mit
         xi = ait+float(is+i-2)*dlx
         call cofx (xi,ai,bi,ci)
         axi = (ai/dlx-0.5*bi)/dlx
         bxi = -2.*ai/dlx**2+ci
         cxi = (ai/dlx+0.5*bi)/dlx
         am(i) = axi
         bm(i) = bxi
         cm(i) = cxi
  110 continue
c
c     set y direction
c
      do 120 j=1,nit
         yj = cit+float(js+j-2)*dly
         call cofy (yj,dj,ej,fj)
         dyj = (dj/dly-0.5*ej)/dly
         eyj = (-2.*dj/dly**2+fj)
         fyj = (dj/dly+0.5*ej)/dly
         an(j) = dyj
         bn(j) = eyj
         cn(j) = fyj
  120 continue
c
c     adjust edges in x direction unless periodic
c
      ax1 = am(1)
      cxm = cm(mit)
      go to (170,130,150,160,140),kswx
c
c     dirichlet-dirichlet in x direction
c
  130 am(1) = 0.0
      cm(mit) = 0.0
      go to 170
c
c     mixed-dirichlet in x direction
c
  140 am(1) = 0.0
      bm(1) = bm(1)+2.*alpha*dlx*ax1
      cm(1) = cm(1)+ax1
      cm(mit) = 0.0
      go to 170
c
c     dirichlet-mixed in x direction
c
  150 am(1) = 0.0
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*beta*dlx*cxm
      cm(mit) = 0.0
      go to 170
c
c     mixed - mixed in x direction
c
  160 continue
      am(1) = 0.0
      bm(1) = bm(1)+2.*dlx*alpha*ax1
      cm(1) = cm(1)+ax1
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*dlx*beta*cxm
      cm(mit) = 0.0
  170 continue
c
c     adjust in y direction unless periodic
c
      dy1 = an(1)
      fyn = cn(nit)
      go to (220,180,200,210,190),kswy
c
c     dirichlet-dirichlet in y direction
c
  180 continue
      an(1) = 0.0
      cn(nit) = 0.0
      go to 220
c
c     mixed-dirichlet in y direction
c
  190 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      cn(nit) = 0.0
      go to 220
c
c     dirichlet-mixed in y direction
c
  200 an(1) = 0.0
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.*dly*xnu*fyn
      cn(nit) = 0.0
      go to 220
c
c     mixed - mixed direction in y direction
c
  210 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.0*dly*xnu*fyn
      cn(nit) = 0.0
  220 if (kswx .eq. 1) go to 270
c
c     adjust usol along x edge
c
      do 260 j=js,ns
         if (kswx.ne.2 .and. kswx.ne.3) go to 230
         usol(is,j) = usol(is,j)-ax1*usol(1,j)
         go to 240
  230    usol(is,j) = usol(is,j)+2.0*dlx*ax1*bda(j)
  240    if (kswx.ne.2 .and. kswx.ne.5) go to 250
         usol(ms,j) = usol(ms,j)-cxm*usol(k,j)
         go to 260
  250    usol(ms,j) = usol(ms,j)-2.0*dlx*cxm*bdb(j)
  260 continue
  270 if (kswy .eq. 1) go to 320
c
c     adjust usol along y edge
c
      do 310 i=is,ms
         if (kswy.ne.2 .and. kswy.ne.3) go to 280
         usol(i,js) = usol(i,js)-dy1*usol(i,1)
         go to 290
  280    usol(i,js) = usol(i,js)+2.0*dly*dy1*bdc(i)
  290    if (kswy.ne.2 .and. kswy.ne.5) go to 300
         usol(i,ns) = usol(i,ns)-fyn*usol(i,l)
         go to 310
  300    usol(i,ns) = usol(i,ns)-2.0*dly*fyn*bdd(i)
  310 continue
  320 continue
c
c     save adjusted edges in grhs if iorder=4
c
      if (iorder .ne. 4) go to 350
      do 330 j=js,ns
         grhs(is,j) = usol(is,j)
         grhs(ms,j) = usol(ms,j)
  330 continue
      do 340 i=is,ms
         grhs(i,js) = usol(i,js)
         grhs(i,ns) = usol(i,ns)
  340 continue
  350 continue
      iord = iorder
      pertrb = 0.0
c
c     check if operator is singular
c
      call chksng (mbdcnd,nbdcnd,alpha,beta,gama,xnu,cofx,cofy,singlr)
c
c     compute non-zero eigenvector in null space of transpose
c     if singular
c
      if (singlr) call trisp (mit,am,bm,cm,dm,um,zm)
      if (singlr) call trisp (nit,an,bn,cn,dn,un,zn)
c
c     make initialization call to BLKTRI.
c
      if (intl .eq. 0)
     &    call blktri (intl,np,nit,an,bn,cn,mp,mit,am,bm,cm,idmn,
     &                 usol(is,js),ierror,w)
      if (ierror .ne. 0) return
c
c     adjust right hand side if necessary
c
  360 continue
      if (singlr) call orthog (usol,idmn,zn,zm,pertrb)
c
c     compute solution
c
      call blktri (i1,np,nit,an,bn,cn,mp,mit,am,bm,cm,idmn,usol(is,js),
     &             ierror,w)
      if (ierror .ne. 0) return
c
c     set periodic boundaries if necessary
c
      if (kswx .ne. 1) go to 380
      do 370 j=1,l
         usol(k,j) = usol(1,j)
  370 continue
  380 if (kswy .ne. 1) go to 400
      do 390 i=1,k
         usol(i,l) = usol(i,1)
  390 continue
  400 continue
c
c     minimize solution with respect to weighted least squares
c     norm if operator is singular
c
      if (singlr) call minsol (usol,idmn,zn,zm,prtrb)
c
c     return if deferred corrections and a fourth order solution are
c     not flagged
c
      if (iord .eq. 2) return
      iord = 2
c
c     compute new right hand side for fourth order solution
c
      call defer (cofx,cofy,idmn,usol,grhs)
      go to 360
      end
      subroutine store (x)

c*********************************************************************72
c
cc STORE forces its argument to be stored.
c
      common /value/  v
      save :: value
!$omp threadprivate (value)
      v = x
      return
      end
      subroutine tevlc (n,d,e2,ierr)

c*********************************************************************72
c
cc TEVLC
c
      integer         i          ,j          ,l          ,m          ,
     &                n          ,ii         ,l1         ,mml        ,
     &                ierr
      real            d(n)       ,e2(n)
      real            b          ,c          ,f          ,g          ,
     &                h          ,p          ,r          ,s          ,
     &                machep
c
c     real sqrt,abs,sign
c
      common /ccblk/  npp        ,k          ,machep     ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: ccblk
!$omp threadprivate (ccblk)
c
c     this subroutine is a modification of the eispack subroutine tqlrat
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e2 contains the                subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        e2 has been destroyed,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
c
      ierr = 0
      if (n .eq. 1) go to 115

      do 101 i=2,n
         e2(i-1) = e2(i)*e2(i)
  101 continue

      f = 0.0
      b = 0.0
      e2(n) = 0.0

      do 112 l=1,n
         j = 0
         h = machep*(abs(d(l))+sqrt(e2(l)))
         if (b .gt. h) go to 102
         b = h
         c = b*b
c
c     ********** look for small squared sub-diagonal element **********
c
  102    do 103 m=l,n
            if (e2(m) .le. c) go to 104
c
c     ********** e2(n) is always zero, so there is no exit
c                through the bottom of the loop **********
c
  103    continue

  104    if (m .eq. l) go to 108
  105    if (j .eq. 30) go to 114
         j = j+1
c
c  form shift
c
         l1 = l+1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1)-g)/(2.0*s)
         r = sqrt(p*p+1.0)
         d(l) = s/(p+sign(r,p))
         h = g-d(l)
c
         do 106 i=l1,n
            d(i) = d(i)-h
  106    continue
c
         f = f+h
c
c     ********** rational ql transformation **********
c
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m-l
c
c     ********** for i=m-1 step -1 until l do -- **********
c
         do 107 ii=1,mml
            i = m-ii
            p = g*h
            r = p+e2(i)
            e2(i+1) = s*r
            s = e2(i)/r
            d(i+1) = h+s*(h+d(i))
            g = d(i)-e2(i)/g
            if (g .eq. 0.0) g = b
            h = g*p/r
  107    continue
c
         e2(l) = s*g
         d(l) = h
c
c     ********** guard against underflowed h **********
c
         if (h .eq. 0.0) go to 108
         if (abs(e2(l)) .le. abs(c/h)) go to 108
         e2(l) = h*e2(l)
         if (e2(l) .ne. 0.0) go to 105
  108    p = d(l)+f
c
c     ********** order eigenvalues **********
c
         if (l .eq. 1) go to 110
c
c     ********** for i=l step -1 until 2 do -- **********
c
         do 109 ii=2,l
            i = l+2-ii
            if (p .ge. d(i-1)) go to 111
            d(i) = d(i-1)
  109    continue
c
  110    i = 1
  111    d(i) = p
  112 continue
c
      if (abs(d(n)) .ge. abs(d(1))) go to 115
      nhalf = n/2
      do 113 i=1,nhalf
         ntop = n-i
         dhold = d(i)
         d(i) = d(ntop+1)
         d(ntop+1) = dhold
  113 continue
      go to 115
c
c     ********** set error -- no convergence to an
c                eigenvalue after 30 iterations **********
c
  114 ierr = l
  115 return

      end
      subroutine tevls (n,d,e2,ierr)

c*********************************************************************72
c
cc TEVLS finds the eigenvalues of a symmetric tridiagonal matrix.
c
      integer         i          ,j          ,l          ,m          ,
     &                n          ,ii         ,l1         ,mml        ,
     &                ierr
      real            d(n)       ,e2(n)
      real            b          ,c          ,f          ,g          ,
     &                h          ,p          ,r          ,s          ,
     &                machep
c
c     real sqrt,abs,sign
c
      common /cblkt/  npp        ,k          ,machep     ,cnv        ,
     &                nm         ,ncmplx     ,ik
     save :: cblkt
!$omp threadprivate (cblkt)
c
c     this subroutine is a modification of the eispack subroutine tqlrat
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e2 contains the                subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        e2 has been destroyed,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
c
      ierr = 0
      if (n .eq. 1) go to 115
c
      do 101 i=2,n
         e2(i-1) = e2(i)*e2(i)
  101 continue
c
      f = 0.0
      b = 0.0
      e2(n) = 0.0
c
      do 112 l=1,n
         j = 0
         h = machep*(abs(d(l))+sqrt(e2(l)))
         if (b .gt. h) go to 102
         b = h
         c = b*b
c
c     ********** look for small squared sub-diagonal element **********
c
  102    do 103 m=l,n
            if (e2(m) .le. c) go to 104
c
c     ********** e2(n) is always zero, so there is no exit
c                through the bottom of the loop
c
  103    continue
c
  104    if (m .eq. l) go to 108
  105    if (j .eq. 30) go to 114
         j = j+1
c
c     ********** form shift
c
         l1 = l+1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1)-g)/(2.0*s)
         r = sqrt(p*p+1.0)
         d(l) = s/(p+sign(r,p))
         h = g-d(l)
c
         do 106 i=l1,n
            d(i) = d(i)-h
  106    continue
c
         f = f+h
c
c     ********** rational ql transformation
c
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m-l
c
c     ********** for i=m-1 step -1 until l do --
c
         do 107 ii=1,mml
            i = m-ii
            p = g*h
            r = p+e2(i)
            e2(i+1) = s*r
            s = e2(i)/r
            d(i+1) = h+s*(h+d(i))
            g = d(i)-e2(i)/g
            if (g .eq. 0.0) g = b
            h = g*p/r
  107    continue
c
         e2(l) = s*g
         d(l) = h
c
c     ********** guard against underflowed h
c
         if (h .eq. 0.0) go to 108
         if (abs(e2(l)) .le. abs(c/h)) go to 108
         e2(l) = h*e2(l)
         if (e2(l) .ne. 0.0) go to 105
  108    p = d(l)+f
c
c     ********** order eigenvalues
c
         if (l .eq. 1) go to 110
c
c     ********** for i=l step -1 until 2 do --
c
         do 109 ii=2,l
            i = l+2-ii
            if (p .ge. d(i-1)) go to 111
            d(i) = d(i-1)
  109    continue
c
  110    i = 1
  111    d(i) = p
  112 continue
c
      if (abs(d(n)) .ge. abs(d(1))) go to 115
      nhalf = n/2
      do 113 i=1,nhalf
         ntop = n-i
         dhold = d(i)
         d(i) = d(ntop+1)
         d(ntop+1) = dhold
  113 continue
      go to 115
c
c  set error -- no convergence to an
c  eigenvalue after 30 iterations
c
  114 ierr = l
  115 return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine tri3 (m,a,b,c,k,y1,y2,y3,tcos,d,w1,w2,w3)

c*********************************************************************72
c
cc TRI3
c
      dimension       a(1)       ,b(1)       ,c(1)       ,k(4)       ,
     &                tcos(1)    ,y1(1)      ,y2(1)      ,y3(1)      ,
     &                d(1)       ,w1(1)      ,w2(1)      ,w3(1)
c
c     routine to solve three linear systems whose common coefficient
c     matrix is a rational function in the matrix given by
c
c                  tridiagonal (...,a(i),b(i),c(i),...)
c
      mm1 = m-1
      k1 = k(1)
      k2 = k(2)
      k3 = k(3)
      k4 = k(4)
      f1 = k1+1
      f2 = k2+1
      f3 = k3+1
      f4 = k4+1
      k2k3k4 = k2+k3+k4
      if (k2k3k4 .eq. 0) go to 101
      l1 = f1/f2
      l2 = f1/f3
      l3 = f1/f4
      lint1 = 1
      lint2 = 1
      lint3 = 1
      kint1 = k1
      kint2 = kint1+k2
      kint3 = kint2+k3
  101 continue
      do 115 n=1,k1
         x = tcos(n)
         if (k2k3k4 .eq. 0) go to 107
         if (n .ne. l1) go to 103
         do 102 i=1,m
            w1(i) = y1(i)
  102    continue
  103    if (n .ne. l2) go to 105
         do 104 i=1,m
            w2(i) = y2(i)
  104    continue
  105    if (n .ne. l3) go to 107
         do 106 i=1,m
            w3(i) = y3(i)
  106    continue
  107    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y1(1) = y1(1)*z
         y2(1) = y2(1)*z
         y3(1) = y3(1)*z
         do 108 i=2,m
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y1(i) = (y1(i)-a(i)*y1(i-1))*z
            y2(i) = (y2(i)-a(i)*y2(i-1))*z
            y3(i) = (y3(i)-a(i)*y3(i-1))*z
  108    continue
         do 109 ip=1,mm1
            i = m-ip
            y1(i) = y1(i)-d(i)*y1(i+1)
            y2(i) = y2(i)-d(i)*y2(i+1)
            y3(i) = y3(i)-d(i)*y3(i+1)
  109    continue
         if (k2k3k4 .eq. 0) go to 115
         if (n .ne. l1) go to 111
         i = lint1+kint1
         xx = x-tcos(i)
         do 110 i=1,m
            y1(i) = xx*y1(i)+w1(i)
  110    continue
         lint1 = lint1+1
         l1 = (float(lint1)*f1)/f2
  111    if (n .ne. l2) go to 113
         i = lint2+kint2
         xx = x-tcos(i)
         do 112 i=1,m
            y2(i) = xx*y2(i)+w2(i)
  112    continue
         lint2 = lint2+1
         l2 = (float(lint2)*f1)/f3
  113    if (n .ne. l3) go to 115
         i = lint3+kint3
         xx = x-tcos(i)
         do 114 i=1,m
            y3(i) = xx*y3(i)+w3(i)
  114    continue
         lint3 = lint3+1
         l3 = (float(lint3)*f1)/f4
  115 continue
      return
      end
      subroutine trid (mr,a,b,c,y,d)

c*********************************************************************72
c
cc TRID
c
      dimension       a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     &                d(1)
      m = mr
      mm1 = m-1
      z = 1./b(1)
      d(1) = c(1)*z
      y(1) = y(1)*z
      do 101 i=2,mm1
         z = 1./(b(i)-a(i)*d(i-1))
         d(i) = c(i)*z
         y(i) = (y(i)-a(i)*y(i-1))*z
  101 continue
      z = b(m)-a(m)*d(mm1)
      if (z .ne. 0.) go to 102
      y(m) = 0.
      go to 103
  102 y(m) = (y(m)-a(m)*y(mm1))/z
  103 continue
      do 104 ip=1,mm1
         i = m-ip
         y(i) = y(i)-d(i)*y(i+1)
  104 continue
      return
      end
      subroutine tris4 (n,a,b,c,d,u,z)

c*********************************************************************72
c
cc TRIS4 solves for a non-zero eigenvector corresponding
c     to the zero eigenvalue of the transpose of the rank
c     deficient one matrix with subdiagonal a, diagonal b, and
c     superdiagonal c , with a(1) in the (1,n) position, with
c     c(n) in the (n,1) position, and all other elements zero.
c
      dimension       a(n)       ,b(n)       ,c(n)       ,d(n)       ,
     &                u(n)       ,z(n)
      bn = b(n)
      d(1) = a(2)/b(1)
      v = a(1)
      u(1) = c(n)/b(1)
      nm2 = n-2
      do  10 j=2,nm2
         den = b(j)-c(j-1)*d(j-1)
         d(j) = a(j+1)/den
         u(j) = -c(j-1)*u(j-1)/den
         bn = bn-v*u(j-1)
         v = -v*d(j-1)
   10 continue
      den = b(n-1)-c(n-2)*d(n-2)
      d(n-1) = (a(n)-c(n-2)*u(n-2))/den
      an = c(n-1)-v*d(n-2)
      bn = bn-v*u(n-2)
      den = bn-an*d(n-1)
c
c     set last component equal to one
c
      z(n) = 1.0
      z(n-1) = -d(n-1)
      nm1 = n-1
      do  20 j=2,nm1
         k = n-j
         z(k) = -d(k)*z(k+1)-u(k)*z(n)
   20 continue
      return
      end
      subroutine trisp (n,a,b,c,d,u,z)

c*********************************************************************72
c
cc TRISP solves for a non-zero eigenvector corresponding
c     to the zero eigenvalue of the transpose of the rank
c     deficient one matrix with subdiagonal a, diagonal b, and
c     superdiagonal c , with a(1) in the (1,n) position, with
c     c(n) in the (n,1) position, and all other elements zero.
c
      dimension       a(n)       ,b(n)       ,c(n)       ,d(n)       ,
     &                u(n)       ,z(n)
      bn = b(n)
      d(1) = a(2)/b(1)
      v = a(1)
      u(1) = c(n)/b(1)
      nm2 = n-2
      do  10 j=2,nm2
         den = b(j)-c(j-1)*d(j-1)
         d(j) = a(j+1)/den
         u(j) = -c(j-1)*u(j-1)/den
         bn = bn-v*u(j-1)
         v = -v*d(j-1)
   10 continue
      den = b(n-1)-c(n-2)*d(n-2)
      d(n-1) = (a(n)-c(n-2)*u(n-2))/den
      an = c(n-1)-v*d(n-2)
      bn = bn-v*u(n-2)
      den = bn-an*d(n-1)
c
c     set last component equal to one
c
      z(n) = 1.0
      z(n-1) = -d(n-1)
      nm1 = n-1
      do  20 j=2,nm1
         k = n-j
         z(k) = -d(k)*z(k+1)-u(k)*z(n)
   20 continue
      return
      end
      subroutine trix (idegbr,idegcr,m,a,b,c,y,tcos,d,w)

c*********************************************************************72
c
cc TRIX solves a system of linear equations where the
c     coefficient matrix is a rational function in the matrix given by
c     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
c
      dimension       a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     &                tcos(1)    ,d(1)       ,w(1)
      mm1 = m-1
      fb = idegbr+1
      fc = idegcr+1
      l = fb/fc
      lint = 1
      do 108 k=1,idegbr
         x = tcos(k)
         if (k .ne. l) go to 102
         i = idegbr+lint
         xx = x-tcos(i)
         do 101 i=1,m
            w(i) = y(i)
            y(i) = xx*y(i)
  101    continue
  102    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y(1) = y(1)*z
         do 103 i=2,mm1
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y(i) = (y(i)-a(i)*y(i-1))*z
  103    continue
         z = b(m)-x-a(m)*d(mm1)
         if (z .ne. 0.) go to 104
         y(m) = 0.
         go to 105
  104    y(m) = (y(m)-a(m)*y(mm1))/z
  105    continue
         do 106 ip=1,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)
  106    continue
         if (k .ne. l) go to 108
         do 107 i=1,m
            y(i) = y(i)+w(i)
  107    continue
         lint = lint+1
         l = (float(lint)*fb)/fc
  108 continue
      return
      end
      subroutine xercnt (librar, subrou, messg, nerr, level, kontrl)

c*********************************************************************72
c
cc XERCNT allows the user to control error handling.
c
c***BEGIN PROLOGUE  XERCNT
c***SUBSIDIARY
c***PURPOSE  Allow user control over handling of errors.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3C
c***TYPE      ALL (XERCNT-A)
c***KEYWORDS  ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        Allows user control over handling of individual errors.
c        Just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to XERCNT.
c        If the user has provided his own version of XERCNT, he
c        can then override the value of KONTROL used in processing
c        this message by redefining its value.
c        KONTRL may be set to any value from -2 to 2.
c        The meanings for KONTRL are the same as in XSETF, except
c        that the value of KONTRL changes only for this message.
c        If KONTRL is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     Description of Parameters
c
c      --Input--
c        LIBRAR - the library that the routine is in.
c        SUBROU - the subroutine that XERMSG is being called from
c        MESSG  - the first 20 characters of the error message.
c        NERR   - same as in the call to XERMSG.
c        LEVEL  - same as in the call to XERMSG.
c        KONTRL - the current value of the control flag as set
c                 by a call to XSETF.
c
c      --Output--
c        KONTRL - the new value of KONTRL.  If KONTRL is not
c                 defined, it will remain at its original value.
c                 This changed value of control affects only
c                 the current occurrence of the current message.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  (NONE)
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900206  Routine changed from user-callable to subsidiary.  (WRB)
c   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
c           names, changed routine name from XERCTL to XERCNT.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XERCNT
      character*(*) librar, subrou, messg
c***FIRST EXECUTABLE STATEMENT  XERCNT
      return
      end
      subroutine xerhlt (messg)

c*********************************************************************72
c
cc XERHLT aborts the program and prints an error message.
c
c***BEGIN PROLOGUE  XERHLT
c***SUBSIDIARY
c***PURPOSE  Abort program execution and print error message.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3C
c***TYPE      ALL (XERHLT-A)
c***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        ***Note*** machine dependent routine
c        XERHLT aborts the execution of the program.
c        The error message causing the abort is given in the calling
c        sequence, in case one needs it for printing on a dayfile,
c        for example.
c
c     Description of Parameters
c        MESSG is as in XERMSG.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  (NONE)
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900206  Routine changed from user-callable to subsidiary.  (WRB)
c   900510  Changed calling sequence to delete length of character
c           and changed routine name from XERABT to XERHLT.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XERHLT
      character*(*) messg
c***FIRST EXECUTABLE STATEMENT  XERHLT
      stop
      end
      subroutine xermax (max)

c*********************************************************************72
c
cc XERMAX sets the maximum number of appearances of an error message.
c
c***BEGIN PROLOGUE  XERMAX
c***PURPOSE  Set maximum number of times any error message is to be
c            printed.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3C
c***TYPE      ALL (XERMAX-A)
c***KEYWORDS  ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        XERMAX sets the maximum number of times any message
c        is to be printed.  That is, non-fatal messages are
c        not to be printed after they have occurred MAX times.
c        Such non-fatal messages may be printed less than
c        MAX times even if they occur MAX times, if error
c        suppression mode (KONTRL=0) is ever in effect.
c
c     Description of Parameter
c      --Input--
c        MAX - the maximum number of times any one message
c              is to be printed.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  J4SAVE
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XERMAX
c***FIRST EXECUTABLE STATEMENT  XERMAX
      junk = j4save(4,max,.true.)
      return
      end
      subroutine xermsg (librar, subrou, messg, nerr, level)

c*********************************************************************72
c
cc XERMSG processes an error message.
c
c***BEGIN PROLOGUE  XERMSG
c***PURPOSE  Process error messages for SLATEC and other libraries.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3C
c***TYPE      ALL (XERMSG-A)
c***KEYWORDS  ERROR MESSAGE, XERROR
c***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
c***DESCRIPTION
c
c   XERMSG processes a diagnostic message in a manner determined by the
c   value of LEVEL and the current value of the library error control
c   flag, KONTRL.  See subroutine XSETF for details.
c
c    LIBRAR   A character constant (or character variable) with the name
c             of the library.  This will be 'SLATEC' for the SLATEC
c             Common Math Library.  The error handling package is
c             general enough to be used by many libraries
c             simultaneously, so it is desirable for the routine that
c             detects and reports an error to identify the library name
c             as well as the routine name.
c
c    SUBROU   A character constant (or character variable) with the name
c             of the routine that detected the error.  Usually it is the
c             name of the routine that is calling XERMSG.  There are
c             some instances where a user callable library routine calls
c             lower level subsidiary routines where the error is
c             detected.  In such cases it may be more informative to
c             supply the name of the routine the user called rather than
c             the name of the subsidiary routine that detected the
c             error.
c
c    MESSG    A character constant (or character variable) with the text
c             of the error or warning message.  In the example below,
c             the message is a character constant that contains a
c             generic message.
c
c                   CALL XERMSG ('SLATEC', 'MMPY',
c                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
c                  *3, 1)
c
c             It is possible (and is sometimes desirable) to generate a
c             specific message--e.g., one that contains actual numeric
c             values.  Specific numeric values can be converted into
c             character strings using formatted WRITE statements into
c             character variables.  This is called standard Fortran
c             internal file I/O and is exemplified in the first three
c             lines of the following example.  You can also catenate
c             substrings of characters to construct the error message.
c             Here is an example showing the use of both writing to
c             an internal file and catenating character strings.
c
c                   CHARACTER*5 CHARN, CHARL
c                   WRITE (CHARN,10) N
c                   WRITE (CHARL,10) LDA
c                10 FORMAT(I5)
c                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
c                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
c                  *   CHARL, 3, 1)
c
c             There are two subtleties worth mentioning.  One is that
c             the // for character catenation is used to construct the
c             error message so that no single character constant is
c             continued to the next line.  This avoids confusion as to
c             whether there are trailing blanks at the end of the line.
c             The second is that by catenating the parts of the message
c             as an actual argument rather than encoding the entire
c             message into one large character variable, we avoid
c             having to know how long the message will be in order to
c             declare an adequate length for that large character
c             variable.  XERMSG calls XERPRN to print the message using
c             multiple lines if necessary.  If the message is very long,
c             XERPRN will break it into pieces of 72 characters (as
c             requested by XERMSG) for printing on multiple lines.
c             Also, XERMSG asks XERPRN to prefix each line with ' *  '
c             so that the total line length could be 76 characters.
c             Note also that XERPRN scans the error message backwards
c             to ignore trailing blanks.  Another feature is that
c             the substring '$$' is treated as a new line sentinel
c             by XERPRN.  If you want to construct a multiline
c             message without having to count out multiples of 72
c             characters, just use '$$' as a separator.  '$$'
c             obviously must occur within 72 characters of the
c             start of each line to have its intended effect since
c             XERPRN is asked to wrap around at 72 characters in
c             addition to looking for '$$'.
c
c    NERR     An integer value that is chosen by the library routine's
c             author.  It must be in the range -99 to 999 (three
c             printable digits).  Each distinct error should have its
c             own error number.  These error numbers should be described
c             in the machine readable documentation for the routine.
c             The error numbers need be unique only within each routine,
c             so it is reasonable for each routine to start enumerating
c             errors from 1 and proceeding to the next integer.
c
c    LEVEL    An integer value in the range 0 to 2 that indicates the
c             level (severity) of the error.  Their meanings are
c
c            -1  A warning message.  This is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.  An attempt is made to only print this
c                message once.
c
c             0  A warning message.  This is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.
c
c             1  A recoverable error.  This is used even if the error is
c                so serious that the routine cannot return any useful
c                answer.  If the user has told the error package to
c                return after recoverable errors, then XERMSG will
c                return to the Library routine which can then return to
c                the user's routine.  The user may also permit the error
c                package to terminate the program upon encountering a
c                recoverable error.
c
c             2  A fatal error.  XERMSG will not return to its caller
c                after it receives a fatal error.  This level should
c                hardly ever be used; it is much better to allow the
c                user a chance to recover.  An example of one of the few
c                cases in which it is permissible to declare a level 2
c                error is a reverse communication Library routine that
c                is likely to be called repeatedly until it integrates
c                across some interval.  If there is a serious error in
c                the input such that another step cannot be taken and
c                the Library routine is called again without the input
c                error having been corrected by the caller, the Library
c                routine will probably be called forever with improper
c                input.  In this case, it is reasonable to declare the
c                error to be fatal.
c
c    Each of the arguments to XERMSG is input; none will be modified by
c    XERMSG.  A routine may make multiple calls to XERMSG with warning
c    level messages; however, after a call to XERMSG with a recoverable
c    error, the routine should return to the user.  Do not try to call
c    XERMSG with a second recoverable error after the first recoverable
c    error because the error package saves the error number.  The user
c    can retrieve this error number by calling another entry point in
c    the error handling package and then clear the error number when
c    recovering from the error.  Calling XERMSG in succession causes the
c    old error number to be overwritten by the latest error number.
c    This is considered harmless for error numbers associated with
c    warning messages but must not be done for error numbers of serious
c    errors.  After a call to XERMSG with a recoverable error, the user
c    must be given a chance to call NUMXER or XERCLR to retrieve or
c    clear the error number.
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
c***REVISION HISTORY  (YYMMDD)
c   880101  DATE WRITTEN
c   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
c           THERE ARE TWO BASIC CHANGES.
c           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
c               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
c               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
c               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
c               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
c               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
c               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
c               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
c           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
c               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
c               OF LOWER CASE.
c   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
c           THE PRINCIPAL CHANGES ARE
c           1.  CLARIFY COMMENTS IN THE PROLOGUES
c           2.  RENAME XRPRNT TO XERPRN
c           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
c               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
c               CHARACTER FOR NEW RECORDS.
c   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
c           CLEAN UP THE CODING.
c   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
c           PREFIX.
c   891013  REVISED TO CORRECT COMMENTS.
c   891214  Prologue converted to Version 4.0 format.  (WRB)
c   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
c           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
c           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
c           XERCTL to XERCNT.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XERMSG
      character*(*) librar, subrou, messg
      character*8 xlibr, xsubr
      character*72  temp
      character*20  lfirst
c***FIRST EXECUTABLE STATEMENT  XERMSG
      lkntrl = j4save (2, 0, .false.)
      maxmes = j4save (4, 0, .false.)
c
c       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
c       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
c          SHOULD BE PRINTED.
c
c       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
c          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
c          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
c
      if (nerr.lt.-9999999 .or. nerr.gt.99999999 .or. nerr.eq.0 .or.
     *   level.lt.-1 .or. level.gt.2) then
         call xerprn (' ***', -1, 'fatal error in...$$ ' //
     *      'xermsg -- invalid error number or level$$ '//
     *      'job abort due to fatal error.', 72)
         call xersve (' ', ' ', ' ', 0, 0, 0, kdummy)
         call xerhlt (' ***xermsg -- invalid input')
         return
      endif
c
c       RECORD THE MESSAGE.
c
      i = j4save (1, nerr, .true.)
      call xersve (librar, subrou, messg, 1, nerr, level, kount)
c
c       HANDLE PRINT-ONCE WARNING MESSAGES.
c
      if (level.eq.-1 .and. kount.gt.1) return
c
c       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
c
      xlibr  = librar
      xsubr  = subrou
      lfirst = messg
      lerr   = nerr
      llevel = level
      call xercnt (xlibr, xsubr, lfirst, lerr, llevel, lkntrl)
c
      lkntrl = max(-2, min(2,lkntrl))
      mkntrl = abs(lkntrl)
c
c       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
c       ZERO AND THE ERROR IS NOT FATAL.
c
      if (level.lt.2 .and. lkntrl.eq.0) go to 30
      if (level.eq.0 .and. kount.gt.maxmes) go to 30
      if (level.eq.1 .and. kount.gt.maxmes .and. mkntrl.eq.1) go to 30
      if (level.eq.2 .and. kount.gt.max(1,maxmes)) go to 30
c
c       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
c       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
c       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
c       IS NOT ZERO.
c
      if (lkntrl .ne. 0) then
         temp(1:21) = 'message from routine '
         i = min(len(subrou), 16)
         temp(22:21+i) = subrou(1:i)
         temp(22+i:33+i) = ' in library '
         ltemp = 33 + i
         i = min(len(librar), 16)
         temp(ltemp+1:ltemp+i) = librar (1:i)
         temp(ltemp+i+1:ltemp+i+1) = '.'
         ltemp = ltemp + i + 1
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
c       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
c       FROM EACH OF THE FOLLOWING THREE OPTIONS.
c       1.  LEVEL OF THE MESSAGE
c              'INFORMATIVE MESSAGE'
c              'POTENTIALLY RECOVERABLE ERROR'
c              'FATAL ERROR'
c       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
c              'PROG CONTINUES'
c              'PROG ABORTED'
c       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
c           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
c           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
c              'TRACEBACK REQUESTED'
c              'TRACEBACK NOT REQUESTED'
c       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
c       EXCEED 74 CHARACTERS.
c       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
c
      if (lkntrl .gt. 0) then
c
c       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
c
         if (level .le. 0) then
            temp(1:20) = 'informative message,'
            ltemp = 20
         elseif (level .eq. 1) then
            temp(1:30) = 'potentially recoverable error,'
            ltemp = 30
         else
            temp(1:12) = 'fatal error,'
            ltemp = 12
         endif
c
c       THEN WHETHER THE PROGRAM WILL CONTINUE.
c
         if ((mkntrl.eq.2 .and. level.ge.1) .or.
     *       (mkntrl.eq.1 .and. level.eq.2)) then
            temp(ltemp+1:ltemp+14) = ' prog aborted,'
            ltemp = ltemp + 14
         else
            temp(ltemp+1:ltemp+16) = ' prog continues,'
            ltemp = ltemp + 16
         endif
c
c       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
c
         if (lkntrl .gt. 0) then
            temp(ltemp+1:ltemp+20) = ' traceback requested'
            ltemp = ltemp + 20
         else
            temp(ltemp+1:ltemp+24) = ' traceback not requested'
            ltemp = ltemp + 24
         endif
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       NOW SEND OUT THE MESSAGE.
c
      call xerprn (' *  ', -1, messg, 72)
c
c       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
c          TRACEBACK.
c
      if (lkntrl .gt. 0) then
         write (temp, '(''error number = '', i8)') nerr
         do 10 i=16,22
            if (temp(i:i) .ne. ' ') go to 20
   10    continue
c
   20    call xerprn (' *  ', -1, temp(1:15) // temp(i:23), 72)
         call fdump
      endif
c
c       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
c
      if (lkntrl .ne. 0) then
         call xerprn (' *  ', -1, ' ', 72)
         call xerprn (' ***', -1, 'end of message', 72)
         call xerprn ('    ',  0, ' ', 72)
      endif
c
c       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
c       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
c
   30 if (level.le.0 .or. (level.eq.1 .and. mkntrl.le.1)) return
c
c       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
c       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
c       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
c
      if (lkntrl.gt.0 .and. kount.lt.max(1,maxmes)) then
         if (level .eq. 1) then
            call xerprn
     *         (' ***', -1, 'job abort due to unrecovered error.', 72)
         else
            call xerprn(' ***', -1, 'job abort due to fatal error.', 72)
         endif
         call xersve (' ', ' ', ' ', -1, 0, 0, kdummy)
         call xerhlt (' ')
      else
         call xerhlt (messg)
      endif
      return
      end
      subroutine xerprn (prefix, npref, messg, nwrap)

c*********************************************************************72
c
cc XERPRN prints an error message.
c
c***BEGIN PROLOGUE  XERPRN
c***SUBSIDIARY
c***PURPOSE  Print error messages processed by XERMSG.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3C
c***TYPE      ALL (XERPRN-A)
c***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
c***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
c***DESCRIPTION
c
c This routine sends one or more lines to each of the (up to five)
c logical units to which error messages are to be sent.  This routine
c is called several times by XERMSG, sometimes with a single line to
c print and sometimes with a (potentially very long) message that may
c wrap around into multiple lines.
c
c PREFIX  Input argument of type CHARACTER.  This argument contains
c         characters to be put at the beginning of each line before
c         the body of the message.  No more than 16 characters of
c         PREFIX will be used.
c
c NPREF   Input argument of type INTEGER.  This argument is the number
c         of characters to use from PREFIX.  If it is negative, the
c         intrinsic function LEN is used to determine its length.  If
c         it is zero, PREFIX is not used.  If it exceeds 16 or if
c         LEN(PREFIX) exceeds 16, only the first 16 characters will be
c         used.  If NPREF is positive and the length of PREFIX is less
c         than NPREF, a copy of PREFIX extended with blanks to length
c         NPREF will be used.
c
c MESSG   Input argument of type CHARACTER.  This is the text of a
c         message to be printed.  If it is a long message, it will be
c         broken into pieces for printing on multiple lines.  Each line
c         will start with the appropriate prefix and be followed by a
c         piece of the message.  NWRAP is the number of characters per
c         piece; that is, after each NWRAP characters, we break and
c         start a new line.  In addition the characters '$$' embedded
c         in MESSG are a sentinel for a new line.  The counting of
c         characters up to NWRAP starts over for each new line.  The
c         value of NWRAP typically used by XERMSG is 72 since many
c         older error messages in the SLATEC Library are laid out to
c         rely on wrap-around every 72 characters.
c
c NWRAP   Input argument of type INTEGER.  This gives the maximum size
c         piece into which to break MESSG for printing on multiple
c         lines.  An embedded '$$' ends a line, and the count restarts
c         at the following character.  If a line break does not occur
c         on a blank (it would split a word) that word is moved to the
c         next line.  Values of NWRAP less than 16 will be treated as
c         16.  Values of NWRAP greater than 132 will be treated as 132.
c         The actual line length will be NPREF + NWRAP after NPREF has
c         been adjusted to fall between 0 and 16 and NWRAP has been
c         adjusted to fall between 16 and 132.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  I1MACH, XGETUA
c***REVISION HISTORY  (YYMMDD)
c   880621  DATE WRITTEN
c   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
c           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
c           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
c           SLASH CHARACTER IN FORMAT STATEMENTS.
c   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
c           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
c           LINES TO BE PRINTED.
c   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
c           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
c   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
c   891214  Prologue converted to Version 4.0 format.  (WRB)
c   900510  Added code to break messages between words.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XERPRN
      character*(*) prefix, messg
      integer npref, nwrap
      character*148 cbuff
      integer iu(5), nunit
      character*2 newlin
      parameter (newlin = '$$')
c***FIRST EXECUTABLE STATEMENT  XERPRN
      call xgetua(iu,nunit)
c
c       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
c       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
c       ERROR MESSAGE UNIT.
c
      n = 6
      do 10 i=1,nunit
         if (iu(i) .eq. 0) iu(i) = n
   10 continue
c
c       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
c       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
c       THE REST OF THIS ROUTINE.
c
      if ( npref .lt. 0 ) then
         lpref = len(prefix)
      else
         lpref = npref
      endif
      lpref = min(16, lpref)
      if (lpref .ne. 0) cbuff(1:lpref) = prefix
c
c       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
c       TIME FROM MESSG TO PRINT ON ONE LINE.
c
      lwrap = max(16, min(132, nwrap))
c
c       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
c
      lenmsg = len(messg)
      n = lenmsg
      do 20 i=1,n
         if (messg(lenmsg:lenmsg) .ne. ' ') go to 30
         lenmsg = lenmsg - 1
   20 continue
   30 continue
c
c       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
c
      if (lenmsg .eq. 0) then
         cbuff(lpref+1:lpref+1) = ' '
         do 40 i=1,nunit
            write(iu(i), '(a)') cbuff(1:lpref+1)
   40    continue
         return
      endif
c
c       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
c       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
c       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
c       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
c
c       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
c       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
c       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
c       OF THE SECOND ARGUMENT.
c
c       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
c       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
c       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
c       POSITION NEXTC.
c
c       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
c                       REMAINDER OF THE CHARACTER STRING.  LPIECE
c                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
c                       WHICHEVER IS LESS.
c
c       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
c                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
c                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
c                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
c                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
c                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
c                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
c                       SHOULD BE INCREMENTED BY 2.
c
c       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
c
c       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
c                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
c                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
c                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
c                       AT THE END OF A LINE.
c
      nextc = 1
   50 lpiece = index(messg(nextc:lenmsg), newlin)
      if (lpiece .eq. 0) then
c
c       THERE WAS NO NEW LINE SENTINEL FOUND.
c
         idelta = 0
         lpiece = min(lwrap, lenmsg+1-nextc)
         if (lpiece .lt. lenmsg+1-nextc) then
            do 52 i=lpiece+1,2,-1
               if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
                  lpiece = i-1
                  idelta = 1
                  goto 54
               endif
   52       continue
         endif
   54    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      elseif (lpiece .eq. 1) then
c
c       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
c       DON'T PRINT A BLANK LINE.
c
         nextc = nextc + 2
         go to 50
      elseif (lpiece .gt. lwrap+1) then
c
c       LPIECE SHOULD BE SET DOWN TO LWRAP.
c
         idelta = 0
         lpiece = lwrap
         do 56 i=lpiece+1,2,-1
            if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
               lpiece = i-1
               idelta = 1
               goto 58
            endif
   56    continue
   58    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      else
c
c       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
c       WE SHOULD DECREMENT LPIECE BY ONE.
c
         lpiece = lpiece - 1
         cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc  = nextc + lpiece + 2
      endif
c
c       PRINT
c
      do 60 i=1,nunit
         write(iu(i), '(a)') cbuff(1:lpref+lpiece)
   60 continue
c
      if (nextc .le. lenmsg) go to 50
      return
      end
      subroutine xersve (librar, subrou, messg, kflag, nerr, level,
     +   icount)

c*********************************************************************72
c
cc XERSVE records that an error has occurred.
c
c***BEGIN PROLOGUE  XERSVE
c***SUBSIDIARY
c***PURPOSE  Record that an error has occurred.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3
c***TYPE      ALL (XERSVE-A)
c***KEYWORDS  ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c *Usage:
c
c        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
c        CHARACTER * (len) LIBRAR, SUBROU, MESSG
c
c        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
c
c *Arguments:
c
c        LIBRAR :IN    is the library that the message is from.
c        SUBROU :IN    is the subroutine that the message is from.
c        MESSG  :IN    is the message to be saved.
c        KFLAG  :IN    indicates the action to be performed.
c                      when KFLAG > 0, the message in MESSG is saved.
c                      when KFLAG=0 the tables will be dumped and
c                      cleared.
c                      when KFLAG < 0, the tables will be dumped and
c                      not cleared.
c        NERR   :IN    is the error number.
c        LEVEL  :IN    is the error severity.
c        ICOUNT :OUT   the number of times this message has been seen,
c                      or zero if the table has overflowed and does not
c                      contain this message specifically.  When KFLAG=0,
c                      ICOUNT will not be altered.
c
c *Description:
c
c   Record that this error occurred and possibly dump and clear the
c   tables.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  I1MACH, XGETUA
c***REVISION HISTORY  (YYMMDD)
c   800319  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900413  Routine modified to remove reference to KFLAG.  (WRB)
c   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
c           sequence, use IF-THEN-ELSE, make number of saved entries
c           easily changeable, changed routine name from XERSAV to
c           XERSVE.  (RWC)
c   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XERSVE
      parameter (lentab=10)
      integer lun(5)
      character*(*) librar, subrou, messg
      character*8  libtab(lentab), subtab(lentab), lib, sub
      character*20 mestab(lentab), mes
      dimension nertab(lentab), levtab(lentab), kount(lentab)
      save libtab, subtab, mestab, nertab, levtab, kount, kountx, nmsg
      data kountx/0/, nmsg/0/
c***FIRST EXECUTABLE STATEMENT  XERSVE
c
      if (kflag.le.0) then
c
c        Dump the table.
c
         if (nmsg.eq.0) return
c
c        Print to each unit.
c
         call xgetua (lun, nunit)
         do 20 kunit = 1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = 6
c
c           Print the table header.
c
            write (iunit,9000)
c
c           Print body of table.
c
            do 10 i = 1,nmsg
               write (iunit,9010) libtab(i), subtab(i), mestab(i),
     *            nertab(i),levtab(i),kount(i)
   10       continue
c
c           Print number of other errors.
c
            if (kountx.ne.0) write (iunit,9020) kountx
            write (iunit,9030)
   20    continue
c
c        Clear the error tables.
c
         if (kflag.eq.0) then
            nmsg = 0
            kountx = 0
         endif
      else
c
c        PROCESS A MESSAGE...
c        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
c        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
c
         lib = librar
         sub = subrou
         mes = messg
         do 30 i = 1,nmsg
            if (lib.eq.libtab(i) .and. sub.eq.subtab(i) .and.
     *         mes.eq.mestab(i) .and. nerr.eq.nertab(i) .and.
     *         level.eq.levtab(i)) then
                  kount(i) = kount(i) + 1
                  icount = kount(i)
                  return
            endif
   30    continue
c
         if (nmsg.lt.lentab) then
c
c           Empty slot found for new message.
c
            nmsg = nmsg + 1
            libtab(i) = lib
            subtab(i) = sub
            mestab(i) = mes
            nertab(i) = nerr
            levtab(i) = level
            kount (i) = 1
            icount    = 1
         else
c
c           Table is full.
c
            kountx = kountx+1
            icount = 0
         endif
      endif
      return
c
c     Formats.
c
 9000 format ('0          error message summary' /
     +   ' library    routine message start             nerr',
     +   '     level     count')
 9010 format (1x,a,3x,a,3x,a,3i10)
 9020 format ('0other errors not individually tabulated = ', i10)
 9030 format (1x)
      end
      subroutine xgetua (iunita, n)

c*********************************************************************72
c
cc XGETUA returns error unit numbers.
c
c***BEGIN PROLOGUE  XGETUA
c***PURPOSE  Return unit number(s) to which error messages are being
c            sent.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3C
c***TYPE      ALL (XGETUA-A)
c***KEYWORDS  ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        XGETUA may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        These unit numbers may have been set by a call to XSETUN,
c        or a call to XSETUA, or may be a default value.
c
c     Description of Parameters
c      --Output--
c        IUNIT - an array of one to five unit numbers, depending
c                on the value of N.  A value of zero refers to the
c                default unit, as defined by the I1MACH machine
c                constant routine.  Only IUNIT(1),...,IUNIT(N) are
c                defined by XGETUA.  The values of IUNIT(N+1),...,
c                IUNIT(5) are not defined (for N .LT. 5) or altered
c                in any way by XGETUA.
c        N     - the number of units to which copies of the
c                error messages are being sent.  N will be in the
c                range from 1 to 5.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  J4SAVE
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XGETUA
      dimension iunita(5)
c***FIRST EXECUTABLE STATEMENT  XGETUA
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunita(i) = j4save(index,0,.false.)
   30 continue
      return
      end
      subroutine xsetf (kontrl)

c*********************************************************************72
c
cc XSETF sets the error control flag.
c
c***BEGIN PROLOGUE  XSETF
c***PURPOSE  Set the error control flag.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3A
c***TYPE      ALL (XSETF-A)
c***KEYWORDS  ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        XSETF sets the error control flag value to KONTRL.
c        (KONTRL is an input parameter only.)
c        The following table shows how each message is treated,
c        depending on the values of KONTRL and LEVEL.  (See XERMSG
c        for description of LEVEL.)
c
c        If KONTRL is zero or negative, no information other than the
c        message itself (including numeric values, if any) will be
c        printed.  If KONTRL is positive, introductory messages,
c        trace-backs, etc., will be printed in addition to the message.
c
c              ABS(KONTRL)
c        LEVEL        0              1              2
c        value
c          2        fatal          fatal          fatal
c
c          1     not printed      printed         fatal
c
c          0     not printed      printed        printed
c
c         -1     not printed      printed        printed
c                                  only           only
c                                  once           once
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  J4SAVE, XERMSG
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890531  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900510  Change call to XERRWV to XERMSG.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XSETF
      character *8 xern1
c***FIRST EXECUTABLE STATEMENT  XSETF
      if (abs(kontrl) .gt. 2) then
         write (xern1, '(i8)') kontrl
         call xermsg ('slatec', 'xsetf',
     *      'invalid argument = ' // xern1, 1, 2)
         return
      endif
c
      junk = j4save(2,kontrl,.true.)
      return
      end
      subroutine xsetun (iunit)

c*********************************************************************72
c
cc XSETUN sets the error message output unit.
c
c***BEGIN PROLOGUE  XSETUN
c***PURPOSE  Set output file to which error messages are to be sent.
c***LIBRARY   SLATEC (XERROR)
c***CATEGORY  R3B
c***TYPE      ALL (XSETUN-A)
c***KEYWORDS  ERROR, XERROR
c***AUTHOR  Jones, R. E., (SNLA)
c***DESCRIPTION
c
c     Abstract
c        XSETUN sets the output file to which error messages are to
c        be sent.  Only one file will be used.  See XSETUA for
c        how to declare more than one file.
c
c     Description of Parameter
c      --Input--
c        IUNIT - an input parameter giving the logical unit number
c                to which error messages are to be sent.
c
c***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
c                 Error-handling Package, SAND82-0800, Sandia
c                 Laboratories, 1982.
c***ROUTINES CALLED  J4SAVE
c***REVISION HISTORY  (YYMMDD)
c   790801  DATE WRITTEN
c   861211  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  XSETUN
c***FIRST EXECUTABLE STATEMENT  XSETUN
      junk = j4save(3,iunit,.true.)
      junk = j4save(5,1,.true.)
      return
      end
