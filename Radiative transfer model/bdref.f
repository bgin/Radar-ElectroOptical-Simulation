c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: BDREF.f,v 2.1 2000/03/27 21:40:51 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )

c      Supplies surface bi-directional reflectivity.
c
c      This is only a "stub" version. The user must replace this
c      by his/her own BDREF function.
c
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c
c   Called by- DREF, SURFAC

c +-------------------------------------------------------------------+
c
c     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..

      WRITE ( *, '(//,7(1X,A,/))' )
     &  'To use a bidirectionally reflecting lower boundary you must',
     &  'replace file BDREF.f with your own file. In that file, you ',
     &  'should supply the bidirectional reflectivity, as a function ',
     &  'of the cosine of angle of reflection, the cosine of angle ',
     &  'of incidence, and the difference of azimuth angles of ',
     &  'incidence and reflection. See DISORT.doc for more information',
     &  'and subroutine BDREF in file DISOTEST.f for an example.'

      CALL ERRMSG( 'BDREF--Please supply a surface BDRF model', .TRUE. )

      BDREF = 0.0

      RETURN
      END
