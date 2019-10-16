
#ifndef __GMS_ODEPACK_F77_H__
#define __GMS_ODEPACK_F77_H__



namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

  const unsigned int gGMS_ODEPACK_F77_MAJOR = gms::common::gVersionInfo.m_VersionMajor;
  const unsigned int gGMS_ODEPACK_F77_MINOR = gms::common::gVersionInfo.m_VersionMinor;
  const unsigned int gGMS_ODEPACK_F77_MICRO = gms::common::gVersionInfo.m_VersionMicro;
  const unsigned int gGMS_ODEPACK_F77_FULLVER =
    1000U*gGMS_ODEPACK_F77_MAJOR+100U*gGMS_ODEPACK_F77_MINOR+10U*gGMS_ODEPACK_F77_MICRO;
  const char * const pgGMS_ODEPACK_F77_CREATE_DATE = "23-12-2017 09:41 +00200 (SAT 23 DEC 2017 GMT+2)";
  const char * const pgGMS_ODEPACK_F77_BUILD_DATE  = "00-00-0000 00:00";
  const char * const pgGMS_ODEPACK_F77_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_ODEPACK_F77_SYNOPSIS    = "C-wrappers to ODEPACK Library.";
}

#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif

#if defined (GMS_CXX_98) || defined (GMS_CXX_11) || defined (GMS_CXX_14)


extern "C" {

typedef void(*F)(int *, double *, double *, double *);
typedef void(*JAC)(int *, double *, double *, int *, int *, double *, int *);
typedef void(*PJAC)(int *, double *, double *, int *, double *, double *, double *, double *, int *, F,JAC);
typedef void(*SLVS)(double *, int *, double *, double *);
typedef void(*G)(int *, double *, double *, int *, double *);
typedef void(*PSOL)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int * ); 


/*
C-----------------------------------------------------------------------
C This routine computes an adjustment, LDIF, to the required
C integer storage space in IWK (sparse matrix work space).
C It is called only if the word length ratio is LRAT = 1.
C This is to account for the possibility that the symbolic LU phase
C may require more storage than the numerical LU and solution phases.
	Parameters:
				INTEGER N, ISP, LDIF
			    DIMENSION ISP(*)
C-----------------------------------------------------------------------
*/

   void ADJLR(int *, int * __restrict , int *);

/*
C-----------------------------------------------------------------------
C This routine counts the number of nonzero elements in the strict
C upper triangle of the matrix M + M(transpose), where the sparsity
C structure of M is given by pointer arrays IA and JA.
C This is needed to compute the storage requirements for the
C sparse matrix reordering operation in ODRV.
  Parameters:
				INTEGER N, IA, JA, NZSUT
			    DIMENSION IA(*), JA(*)
C-----------------------------------------------------------------------
*/
   void CNTNZU(int *, int * __restrict, int * __restrict, int *);

/*
C-----------------------------------------------------------------------
C This routine manages the solution of the linear system arising from
C a chord iteration.  It is called if MITER .ne. 0.
C If MITER is 1 or 2, it calls CDRV to accomplish this.
C If MITER = 3 it updates the coefficient H*EL0 in the diagonal
C matrix, and then computes the solution.
C communication with DSOLSS uses the following variables:
C WK    = real work space containing the inverse diagonal matrix if
C         MITER = 3 and the LU decomposition of the matrix otherwise.
C         Storage of matrix elements starts at WK(3).
C         WK also contains the following matrix-related data:
C         WK(1) = SQRT(UROUND) (not used here),
C         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
C IWK   = integer work space for matrix-related data, assumed to
C         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
C         are assumed to have identical locations.
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C TEM   = vector of work space of length N, not used in this version.
C IERSL = output flag (in Common).
C         IERSL = 0  if no trouble occurred.
C         IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
C                    This should never occur and is considered fatal.
C         IERSL = 1  if a singular matrix arose with MITER = 3.
C This routine also uses other variables in Common.
C-----------------------------------------------------------------------
*/

   void DSOLSS(double * __restrict, int * __restrict, double * __restrict, double * __restrict);

/*
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSS01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 224 or more.
C ISAV = integer array of length 71 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
*/
   void DSRCMS(double *, int *, int *);

/*
	c                                                                 5/2/83
c***********************************************************************
c  odrv -- driver for sparse matrix reordering routines
c***********************************************************************
c
c  description
c
c    odrv finds a minimum degree ordering of the rows and columns
c    of a matrix m stored in (ia,ja,a) format (see below).  for the
c    reordered matrix, the work and storage required to perform
c    gaussian elimination is (usually) significantly less.
c
c    note.. odrv and its subordinate routines have been modified to
c    compute orderings for general matrices, not necessarily having any
c    symmetry.  the miminum degree ordering is computed for the
c    structure of the symmetric matrix  m + m-transpose.
c    modifications to the original odrv module have been made in
c    the coding in subroutine mdi, and in the initial comments in
c    subroutines odrv and md.
c
c    if only the nonzero entries in the upper triangle of m are being
c    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
c    with the diagonal entries placed first in each row.  this is to
c    ensure that if m(i,j) will be in the upper triangle of m with
c    respect to the new ordering, then m(i,j) is stored in row i (and
c    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
c    strict lower triangle of m, then m(j,i) is stored in row j (and
c    thus m(i,j) is not stored).
c
c
c  storage of sparse matrices
c
c    the nonzero entries of the matrix m are stored row-by-row in the
c    array a.  to identify the individual nonzero entries in each row,
c    we need to know in which column each entry lies.  these column
c    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  to identify the individual rows, we need to know where
c    each row starts.  these row pointers are stored in the array ia.
c    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
c    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
c    the first location following the last element in the last row.
c    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
c    the nonzero entries in the i-th row are stored consecutively in
c
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c
c    and the corresponding column indices are stored consecutively in
c
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c
c    when the coefficient matrix is symmetric, only the nonzero entries
c    in the upper triangle need be stored.  for example, the matrix
c
c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )
c
c    could be stored as
c
c            - 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia - 1  4  5  8 12 14
c         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
c
c    or (symmetrically) as
c
c            - 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia - 1  4  5  7  9 10
c         ja - 1  3  4  2  3  4  4  5  5
c          a - 1  2  3  4  5  6  7  8  9          .
c
c
c  parameters
c
c    n    - order of the matrix
c
c    ia   - integer one-dimensional array containing pointers to delimit
c           rows in ja and a.  dimension = n+1
c
c    ja   - integer one-dimensional array containing the column indices
c           corresponding to the elements of a.  dimension = number of
c           nonzero entries in (the upper triangle of) m
c
c    a    - real one-dimensional array containing the nonzero entries in
c           (the upper triangle of) m, stored by rows.  dimension =
c           number of nonzero entries in (the upper triangle of) m
c
c    p    - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    ip   - integer one-dimensional array used to return the inverse of
c           the permutation returned in p.  dimension = n
c
c    nsp  - declared dimension of the one-dimensional array isp.  nsp
c           must be at least  3n+4k,  where k is the number of nonzeroes
c           in the strict upper triangle of m
c
c    isp  - integer one-dimensional array used for working storage.
c           dimension = nsp
c
c    path - integer path specification.  values and their meanings are -
c             1  find minimum degree ordering only
c             2  find minimum degree ordering and reorder symmetrically
c                  stored matrix (used when only the nonzero entries in
c                  the upper triangle of m are being stored)
c             3  reorder symmetrically stored matrix as specified by
c                  input permutation (used when an ordering has already
c                  been determined and only the nonzero entries in the
c                  upper triangle of m are being stored)
c             4  same as 2 but put diagonal entries at start of each row
c             5  same as 3 but put diagonal entries at start of each row
c
c    flag - integer error flag.  values and their meanings are -
c               0    no errors detected
c              9n+k  insufficient storage in md
c             10n+1  insufficient storage in odrv
c             11n+1  illegal path specification
c
c
c  conversion from real to double precision
c
c    change the real declarations in odrv and sro to double precision
c    declarations.
c
c-----------------------------------------------------------------------
c
*/
	void odrv(int *, int * __restrict, int * __restrict, double * __restrict, int * __restrict, int * __restrict,
			  int * , int *, int *, int * );

/*
c***********************************************************************
c  md -- minimum degree algorithm (based on element model)
c***********************************************************************
c
c  description
c
c    md finds a minimum degree ordering of the rows and columns of a
c    general sparse matrix m stored in (ia,ja,a) format.
c    when the structure of m is nonsymmetric, the ordering is that
c    obtained for the symmetric matrix  m + m-transpose.
c
c
c  additional parameters
c
c    max  - declared dimension of the one-dimensional arrays v and l.
c           max must be at least  n+2k,  where k is the number of
c           nonzeroes in the strict upper triangle of m + m-transpose
c
c    v    - integer one-dimensional work array.  dimension = max
c
c    l    - integer one-dimensional work array.  dimension = max
c
c    head - integer one-dimensional work array.  dimension = n
c
c    last - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    next - integer one-dimensional array used to return the inverse of
c           the permutation returned in last.  dimension = n
c
c    mark - integer one-dimensional work array (may be the same as v).
c           dimension = n
c
c    flag - integer error flag.  values and their meanings are -
c             0     no errors detected
c             9n+k  insufficient storage in md
c
c
c  definitions of internal parameters
c
c    ---------+---------------------------------------------------------
c    v(s)     - value field of list entry
c    ---------+---------------------------------------------------------
c    l(s)     - link field of list entry  (0 =) end of list)
c    ---------+---------------------------------------------------------
c    l(vi)    - pointer to element list of uneliminated vertex vi
c    ---------+---------------------------------------------------------
c    l(ej)    - pointer to boundary list of active element ej
c    ---------+---------------------------------------------------------
c    head(d)  - vj =) vj head of d-list d
c             -  0 =) no vertex in d-list d
c
c
c             -                  vi uneliminated vertex
c             -          vi in ek           -       vi not in ek
c    ---------+-----------------------------+---------------------------
c    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
c             -                             -  0 =) vi tail of d-list
c    ---------+-----------------------------+---------------------------
c    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
c             --vk =) compute degree        - vj =) vj last in d-list
c             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
c             -  0 =) do not compute degree -
c    ---------+-----------------------------+---------------------------
c    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
c
c
c             -                   vi eliminated vertex
c             -      ei active element      -           otherwise
c    ---------+-----------------------------+---------------------------
c    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
c             -       to be eliminated      -       to be eliminated
c    ---------+-----------------------------+---------------------------
c    last(vi) -  m =) size of ei = m        - undefined
c    ---------+-----------------------------+---------------------------
c    mark(vi) - -m =) overlap count of ei   - undefined
c             -       with ek = m           -
c             - otherwise nonnegative tag   -
c             -       .lt. mark(vk)         -
c
c-----------------------------------------------------------------------
c
*/
	void md(int *, int * __restrict, int * __restrict, int *, int * __restrict, int * __restrict,
			int * __restrict, int * __restrict, int * __restrict, int * __restrict, int *  );

/*
c***********************************************************************
c  mdi -- initialization
c***********************************************************************
*/
	void mdi(int *, int * __restrict, int * __restrict, int *, int * __restrict, int * __restrict,
			 int * __restrict, int * __restrict, int * __restrict, int * __restrict, int *, int *);

/*
c***********************************************************************
c  mdm -- form element from uneliminated neighbors of vk
c***********************************************************************
*/
	void mdm(int *, int *, int * __restrict, int * __restrict, int * 
			 __restrict, int * __restrict, int * __restrict);

/*
c***********************************************************************
c  mdp -- purge inactive elements and do mass elimination
c***********************************************************************
*/
	void mdp(int *, int *, int *, int *, int *,
			 int *, int *, int *, int *);

/*
c***********************************************************************
c  mdu -- update degrees of uneliminated vertices in ek
c***********************************************************************
*/
	void mdu(int *, int *, int *, int *, int *, int *, int *, int *);

/*
c***********************************************************************
c  sro -- symmetric reordering of sparse symmetric matrix
c***********************************************************************
c
c  description
c
c    the nonzero entries of the matrix m are assumed to be stored
c    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
c    are stored if i ne j).
c
c    sro does not rearrange the order of the rows, but does move
c    nonzeroes from one row to another to ensure that if m(i,j) will be
c    in the upper triangle of m with respect to the new ordering, then
c    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
c    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
c    stored in row j (and thus m(i,j) is not stored).
c
c
c  additional parameters
c
c    q     - integer one-dimensional work array.  dimension = n
c
c    r     - integer one-dimensional work array.  dimension = number of
c            nonzero entries in the upper triangle of m
c
c    dflag - logical variable.  if dflag = .true., then store nonzero
c            diagonal elements at the beginning of the row
c
c-----------------------------------------------------------------------
*/
	void sro(int *, int *, int *, int *, double *, int *, int *, int *);

/*
c*** subroutine cdrv
c*** driver for subroutines for solving sparse nonsymmetric systems of
c       linear equations (compressed pointer storage)
c
c
c    parameters
c    class abbreviations are--
c       n - integer variable
c       f - real variable
c       v - supplies a value to the driver
c       r - returns a result from the driver
c       i - used internally by the driver
c       a - array
c
c class - parameter
c ------+----------
c       -
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c nv    - n     - number of variables/equations.
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c nva   - ia    - pointers to delimit the rows in a.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c fva   - b     - right-hand side b.  b and z can the same array.
c       -           size = n.
c fra   - z     - solution x.  b and z can be the same array.
c       -           size = n.
c
c         the rows and columns of the original matrix m can be
c    reordered (e.g., to reduce fillin or ensure numerical stability)
c    before calling the driver.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
c    in the original order.
c         if the columns have been reordered (i.e.,  c(i).ne.i  for some
c    i), then the driver will call a subroutine (nroc) which rearranges
c    each row of ja and a, leaving the rows in the original order, but
c    placing the elements of each row in increasing order with respect
c    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
c    been called already.
c
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
c       -           ic(c(i)) = i  for i=1,...,n.
c       -           size = n.
c
c         the solution of the system of linear equations is divided into
c    three stages --
c      nsfc -- the matrix m is processed symbolically to determine where
c               fillin will occur during the numeric factorization.
c      nnfc -- the matrix m is factored numerically into the product ldu
c               of a unit lower triangular matrix l, a diagonal matrix
c               d, and a unit upper triangular matrix u, and the system
c               mx = b  is solved.
c      nnsc -- the linear system  mx = b  is solved using the ldu
c  or           factorization from nnfc.
c      nntc -- the transposed linear system  mt x = b  is solved using
c               the ldu factorization from nnf.
c    for several systems whose coefficient matrices have the same
c    nonzero structure, nsfc need be done only once (for the first
c    system).  then nnfc is done once for each additional system.  for
c    several systems with the same coefficient matrix, nsfc and nnfc
c    need be done only once (for the first system).  then nnsc or nntc
c    is done once for each additional right-hand side.
c
c nv    - path  - path specification.  values and their meanings are --
c       -           1  perform nroc, nsfc, and nnfc.
c       -           2  perform nnfc only  (nsfc is assumed to have been
c       -               done in a manner compatible with the storage
c       -               allocation used in the driver).
c       -           3  perform nnsc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           4  perform nntc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           5  perform nroc and nsfc.
c
c         various errors are detected by the driver and the individual
c    subroutines.
c
c nr    - flag  - error flag.  values and their meanings are --
c       -             0     no errors detected
c       -             n+k   null row in a  --  row = k
c       -            2n+k   duplicate entry in a  --  row = k
c       -            3n+k   insufficient storage in nsfc  --  row = k
c       -            4n+1   insufficient storage in nnfc
c       -            5n+k   null pivot  --  row = k
c       -            6n+k   insufficient storage in nsfc  --  row = k
c       -            7n+1   insufficient storage in nnfc
c       -            8n+k   zero pivot  --  row = k
c       -           10n+1   insufficient storage in cdrv
c       -           11n+1   illegal path specification
c
c         working storage is needed for the factored form of the matrix
c    m plus various temporary vectors.  the arrays isp and rsp should be
c    equivalenced.  integer storage is allocated from the beginning of
c    isp and real storage from the end of rsp.
c
c nv    - nsp   - declared dimension of rsp.  nsp generally must
c       -           be larger than  8n+2 + 2k  (where  k = (number of
c       -           nonzero entries in m)).
c nvira - isp   - integer working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = lratio*nsp.
c fvira - rsp   - real working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = nsp.
c nr    - esp   - if sufficient storage was available to perform the
c       -           symbolic factorization (nsfc), then esp is set to
c       -           the amount of excess storage provided (negative if
c       -           insufficient storage was available to perform the
c       -           numeric factorization (nnfc)).
c
c
c  conversion to double precision
c
c    to convert these routines for double precision arrays..
c    (1) use the double precision declarations in place of the real
c    declarations in each subprogram, as given in comment cards.
c    (2) change the data-loaded value of the integer  lratio
c    in subroutine cdrv, as indicated below.
c    (3) change e0 to d0 in the constants in statement number 10
c    in subroutine nnfc and the line following that.
c
*/
	void cdrv(int *, int *, int *, int *, int *, int *, double *,
			  double *, double *,  int *, int *, double *, int *, int *, int *);

/*
c
c       ----------------------------------------------------------------
c
c               yale sparse matrix package - nonsymmetric codes
c                    solving the system of equations mx = b
c
c    i.   calling sequences
c         the coefficient matrix can be processed by an ordering routine
c    (e.g., to reduce fillin or ensure numerical stability) before using
c    the remaining subroutines.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
c    is used, then nroc should be used to reorder the coefficient matrix
c    the calling sequence is --
c        (       (matrix ordering))
c        (nroc   (matrix reordering))
c         nsfc   (symbolic factorization to determine where fillin will
c                  occur during numeric factorization)
c         nnfc   (numeric factorization into product ldu of unit lower
c                  triangular matrix l, diagonal matrix d, and unit
c                  upper triangular matrix u, and solution of linear
c                  system)
c         nnsc   (solution of linear system for additional right-hand
c                  side using ldu factorization from nnfc)
c    (if only one system of equations is to be solved, then the
c    subroutine trk should be used.)
c
c    ii.  storage of sparse matrices
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c         the strict upper (lower) triangular portion of the matrix
c    u (l) is stored in a similar fashion using the arrays  iu, ju, u
c    (il, jl, l)  except that an additional array iju (ijl) is used to
c    compress storage of ju (jl) by allowing some sequences of column
c    (row) indices to used for more than one row (column)  (n.b., l is
c    stored by columns).  iju(k) (ijl(k)) points to the starting
c    location in ju (jl) of entries for the kth row (column).
c    compression in ju (jl) occurs in two ways.  first, if a row
c    (column) i was merged into the current row (column) k, and the
c    number of elements merged in from (the tail portion of) row
c    (column) i is the same as the final length of row (column) k, then
c    the kth row (column) and the tail of row (column) i are identical
c    and iju(k) (ijl(k)) points to the start of the tail.  second, if
c    some tail portion of the (k-1)st row (column) is identical to the
c    head of the kth row (column), then iju(k) (ijl(k)) points to the
c    start of that tail portion.  for example, the nonzero structure of
c    the strict upper triangular part of the matrix
c            d 0 x x x
c            0 d 0 x x
c            0 0 d x 0
c            0 0 0 d x
c            0 0 0 0 d
c    would be represented as
c                - 1 2 3 4 5 6
c            ----+------------
c             iu - 1 4 6 7 8 8
c             ju - 3 4 5 4
c            iju - 1 2 4 3           .
c    the diagonal entries of l and u are assumed to be equal to one and
c    are not stored.  the array d contains the reciprocals of the
c    diagonal entries of the matrix d.
c
c    iii. additional storage savings
c         in nsfc, r and ic can be the same array in the calling
c    sequence if no reordering of the coefficient matrix has been done.
c         in nnfc, r, c, and ic can all be the same array if no
c    reordering has been done.  if only the rows have been reordered,
c    then c and ic can be the same array.  if the row and column
c    orderings are the same, then r and c can be the same array.  z and
c    row can be the same array.
c         in nnsc or nntc, r and c can be the same array if no
c    reordering has been done or if the row and column orderings are the
c    same.  z and b can be the same array.  however, then b will be
c    destroyed.
c
c    iv.  parameters
c         following is a list of parameters to the programs.  names are
c    uniform among the various subroutines.  class abbreviations are --
c       n - integer variable
c       f - real variable
c       v - supplies a value to a subroutine
c       r - returns a result from a subroutine
c       i - used internally by a subroutine
c       a - array
c
c class - parameter
c ------+----------
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c fva   - b     - right-hand side b.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c fvra  - d     - reciprocals of the diagonal entries of the matrix d.
c       -           size = n.
c nr    - flag  - error flag.  values and their meanings are --
c       -            0     no errors detected
c       -            n+k   null row in a  --  row = k
c       -           2n+k   duplicate entry in a  --  row = k
c       -           3n+k   insufficient storage for jl  --  row = k
c       -           4n+1   insufficient storage for l
c       -           5n+k   null pivot  --  row = k
c       -           6n+k   insufficient storage for ju  --  row = k
c       -           7n+1   insufficient storage for u
c       -           8n+k   zero pivot  --  row = k
c nva   - ia    - pointers to delimit the rows of a.
c       -           size = n+1.
c nvra  - ijl   - pointers to the first element in each column in jl,
c       -           used to compress storage in jl.
c       -           size = n.
c nvra  - iju   - pointers to the first element in each row in ju, used
c       -           to compress storage in ju.
c       -           size = n.
c nvra  - il    - pointers to delimit the columns of l.
c       -           size = n+1.
c nvra  - iu    - pointers to delimit the rows of u.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c nvra  - jl    - row numbers corresponding to the elements of l.
c       -           size = jlmax.
c nv    - jlmax - declared dimension of jl.  jlmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin minus compression.
c nvra  - ju    - column numbers corresponding to the elements of u.
c       -           size = jumax.
c nv    - jumax - declared dimension of ju.  jumax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin minus compression.
c fvra  - l     - nonzero entries in the strict lower triangular portion
c       -           of the matrix l, stored by columns.
c       -           size = lmax.
c nv    - lmax  - declared dimension of l.  lmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin  (il(n+1)-1 after nsfc).
c nv    - n     - number of variables/equations.
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c fvra  - u     - nonzero entries in the strict upper triangular portion
c       -           of the matrix u, stored by rows.
c       -           size = umax.
c nv    - umax  - declared dimension of u.  umax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin  (iu(n+1)-1 after nsfc).
c fra   - z     - solution x.
c       -           size = n.
c
c       ----------------------------------------------------------------
c
c*** subroutine nroc
c*** reorders rows of a, leaving row order unchanged
c
c
c       input parameters.. n, ic, ia, ja, a
c       output parameters.. ja, a, flag
c
c       parameters used internally..
c nia   - p     - at the kth step, p is a linked list of the reordered
c       -           column indices of the kth row of a.  p(n+1) points
c       -           to the first entry in the list.
c       -           size = n+1.
c nia   - jar   - at the kth step,jar contains the elements of the
c       -           reordered column indices of a.
c       -           size = n.
c fia   - ar    - at the kth step, ar contains the elements of the
c       -           reordered row of a.
c       -           size = n.
*/
	void nroc(int *, int *, int *, int *, double *, int *, double *, int *, int *);

/*
c*** subroutine nsfc
c*** symbolic ldu-factorization of nonsymmetric sparse matrix
c      (compressed pointer storage)
c
c
c       input variables.. n, r, ic, ia, ja, jlmax, jumax.
c       output variables.. il, jl, ijl, iu, ju, iju, flag.
c
c       parameters used internally..
c nia   - q     - suppose  m*  is the result of reordering  m.  if
c       -           processing of the ith row of  m*  (hence the ith
c       -           row of  u) is being done,  q(j)  is initially
c       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
c       -           values need not be stored, each entry points to the
c       -           next nonzero and  q(n+1)  points to the first.  n+1
c       -           indicates the end of the list.  for example, if n=9
c       -           and the 5th row of  m*  is
c       -              0 x x 0 x 0 0 x 0
c       -           then  q  will initially be
c       -              a a a a 8 a a 10 5           (a - arbitrary).
c       -           as the algorithm proceeds, other elements of  q
c       -           are inserted in the list because of fillin.
c       -           q  is used in an analogous manner to compute the
c       -           ith column of  l.
c       -           size = n+1.
c nia   - ira,  - vectors used to find the columns of  m.  at the kth
c nia   - jra,      step of the factorization,  irac(k)  points to the
c nia   - irac      head of a linked list in  jra  of row indices i
c       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
c       -           indicates the end of the list.  ira(i)  (i.ge.k)
c       -           points to the smallest j such that j .ge. k and
c       -           m(i,j)  is nonzero.
c       -           size of each = n.
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c nia   - iru,  - vectors used in a manner analogous to  irl and jrl
c nia   - jru       to find the columns of  u.
c       -           size of each = n.
c
c  internal variables..
c    jlptr - points to the last position used in  jl.
c    juptr - points to the last position used in  ju.
c    jmin,jmax - are the indices in  a or u  of the first and last
c                elements to be examined in a given row.
c                for example,  jmin=ia(k), jmax=ia(k+1)-1.
c
*/
	void nsfc(int *, int *, int *, int *, int *, int *, int *,
			  int *, int *, int *, int *, int *, int *, int *,
			  int *, int *, int *, int *, int *, int *, int *, int *);

/*
c*** subroutine nnfc
c*** numerical ldu-factorization of sparse nonsymmetric matrix and
c      solution of system of linear equations (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, ic, ia, ja, a, b,
c                          il, jl, ijl, lmax, iu, ju, iju, umax
c       output variables.. z, l, d, u, flag
c
c       parameters used internally..
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c fia   - row   - holds intermediate values in calculation of  u and l.
c       -           size = n.
c fia   - tmp   - holds new right-hand side  b*  for solution of the
c       -           equation ux = b*.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row to
c      be examined.
c    sum - used in calculating  tmp.
*/
	void nnfc(int *, int *, int *, int *, int *, int *, double *, double *, double *,
			  int *, int *, int *, int *, double *, double *, int *, int *, int *, int *,
			  double *, double *, double *, int *, int *, int *);

/*
c*** subroutine nnsc
c*** numerical solution of sparse nonsymmetric system of linear
c      equations given ldu-factorization (compressed pointer storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving  ly = b.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
*/

	void nnsc(int *, int *, int *, int *, int *, int *, double *, double *,
			  int *, int *, int *, double *, double *, double *, double *);

/*
c*** subroutine nntc
c*** numeric solution of the transpose of a sparse nonsymmetric system
c      of linear equations given lu-factorization (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving ut y = b
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
*/
	void nntc(int *, int *, int *, int *, int *, int *, double *, double *,
			  int *, int *, int *, double *, double *, double *, double * );

/*
C-----------------------------------------------------------------------
C DSTODA performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note: DSTODA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODA is done with the following variables:
C
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C          and P = I - H*EL0*Jac, if a chord method is being used.
C          It also returns an estimate of norm(Jac) in PDNORM.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH   = current method.
C          METH = 1 means Adams method (nonstiff)
C          METH = 2 means BDF method (stiff)
C          METH may be reset by DSTODA.
C MITER  = corrector iteration method.
C          MITER = 0 means functional iteration.
C          MITER = JT .gt. 0 means a chord iteration corresponding
C          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
C          communicated here as JTYP, but is not used in DSTODA
C          except to load MITER following a method switch.)
C          MITER may be reset by DSTODA.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
*/
#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	// If this does not work properly then fall back to alternative.
	void DSTODA(int *, double *, double *, int *, double *, double *, double *,
				double *, double *, int *,F,JAC,PJAC,SLVS);
#else

	void DSTODA(int *, double *, double *, int *, double *, double *, double *,
		double *, double *, int *, void (*)(int *, double *, double *, double *),
		void (*)(int *, double *, double *, int *, int *, double *, int *),
		void(*)(int *, double *, double *, int *, double *, double *, double *,
		double *, int *, void(*)(int *, double *, double *, double *),
		void(*)(int *, double *, double *, int *, int *, double *, int *)),
		void(*)(double *, int *, double *, double *));

#endif

/*
	DPRJA
	C-----------------------------------------------------------------------
	C DPRJA is called by DSTODA to compute and process the matrix
	C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
	C Here J is computed by the user-supplied routine JAC if
	C MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
	C J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
	C matrix norm consistent with the weighted max-norm on vectors given
	C by DMNORM) is computed, and J is overwritten by P.  P is then
	C subjected to LU decomposition in preparation for later solution
	C of linear systems with P as coefficient matrix.  This is done
	C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
	C
	C In addition to variables described previously, communication
	C with DPRJA uses the following:
	C Y     = array containing predicted values on entry.
	C FTEM  = work array of length N (ACOR in DSTODA).
	C SAVF  = array containing f evaluated at predicted y.
	C WM    = real work space for matrices.  On output it contains the
	C         LU decomposition of P.
	C         Storage of matrix elements starts at WM(3).
	C         WM also contains the following matrix-related data:
	C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
	C IWM   = integer work space containing pivot information, starting at
	C         IWM(21).   IWM also contains the band parameters
	C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
	C EL0   = EL(1) (input).
	C PDNORM= norm of Jacobian matrix. (Output).
	C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
	C         P matrix found to be singular.
	C JCUR  = output flag = 1 to indicate that the Jacobian matrix
	C         (or approximation) is now current.
	C This routine also uses the Common variables EL0, H, TN, UROUND,
	C MITER, N, NFE, and NJE.
	C-----------------------------------------------------------------------
*/
#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	void DPRJA(int *, double *, double *, int *, double *, double *, double *,
			   double *, int *, F,JAC);

#else
	
	void DPRJA(int *, double *, double *, int *, double *, double *, double *,
		       double *, int *, void (*)(int *, double *, double *, double *),
			   void (*)(int *, double *, double *, int *, int *, double *, int *) );

#endif

/*
C-----------------------------------------------------------------------
C This function routine computes the weighted max-norm
C of the vector of length N contained in the array V, with weights
C contained in the array w of length N:
C   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
C-----------------------------------------------------------------------
*/
	double DMNORM(int *, double *, double *);

/*
C-----------------------------------------------------------------------
C This function computes the norm of a full N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W:
C   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
*/
	double DFNORM(int *, double *, double *);
	
/*
C-----------------------------------------------------------------------
C This function computes the norm of a banded N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W.
C ML and MU are the lower and upper half-bandwidths of the matrix.
C NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
C In terms of the matrix elements a(i,j), the norm is given by:
C   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
*/
	double DBNORM(int *, double *, int *, int *, int *, double *);

/*
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 240 or more.
C ISAV = integer array of length 46 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
*/
	void DSRCMA(double *, int *, int *);

/*
C-----------------------------------------------------------------------
C This routine checks for the presence of a root in the vicinity of
C the current T, in a manner depending on the input flag JOB.  It calls
C Subroutine DROOTS to locate the root as precisely as possible.
C
C In addition to variables described previously, DRCHEK
C uses the following for communication:
C JOB    = integer flag indicating type of call:
C          JOB = 1 means the problem is being initialized, and DRCHEK
C                  is to look for a root at or very near the initial T.
C          JOB = 2 means a continuation call to the solver was just
C                  made, and DRCHEK is to check for a root in the
C                  relevant part of the step last taken.
C          JOB = 3 means a successful step was just taken, and DRCHEK
C                  is to look for a root in the interval of the step.
C G0     = array of length NG, containing the value of g at T = T0.
C          G0 is input for JOB .ge. 2, and output in all cases.
C G1,GX  = arrays of length NG for work space.
C IRT    = completion flag:
C          IRT = 0  means no root was found.
C          IRT = -1 means JOB = 1 and a root was found too near to T.
C          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
C                   On return, T0 is the root location, and Y is the
C                   corresponding solution vector.
C T0     = value of T at one endpoint of interval of interest.  Only
C          roots beyond T0 in the direction of integration are sought.
C          T0 is input if JOB .ge. 2, and output in all cases.
C          T0 is updated by DRCHEK, whether a root is found or not.
C TLAST  = last value of T returned by the solver (input only).
C TOUTC  = copy of TOUT (input only).
C IRFND  = input flag showing whether the last step taken had a root.
C          IRFND = 1 if it did, = 0 if not.
C ITASKC = copy of ITASK (input only).
C NGC    = copy of NG (input only).
C-----------------------------------------------------------------------
*/
#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1

	void DRCHEK(int *, G, int *, double *, double *, int *, double *, double *, double *,
				int *, int *);

#else
	
	void DRCHEK(int *, void(*)(int *, double *, double *, int *, double *),
				int *, double *, double *, int *, double *, double *, double *,
				int *, int * );

#endif

/*
C-----------------------------------------------------------------------
C This subroutine finds the leftmost root of a set of arbitrary
C functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
C of odd multiplicity (i.e. changes of sign of the gi) are found.
C Here the sign of X1 - X0 is arbitrary, but is constant for a given
C problem, and -leftmost- means nearest to X0.
C The values of the vector-valued function g(x) = (gi, i=1...NG)
C are communicated through the call sequence of DROOTS.
C The method used is the Illinois algorithm.
C
C Reference:
C Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
C Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
C February 1980.
C
C Description of parameters.
C
C NG     = number of functions gi, or the number of components of
C          the vector valued function g(x).  Input only.
C
C HMIN   = resolution parameter in X.  Input only.  When a root is
C          found, it is located only to within an error of HMIN in X.
C          Typically, HMIN should be set to something on the order of
C               100 * UROUND * MAX(ABS(X0),ABS(X1)),
C          where UROUND is the unit roundoff of the machine.
C
C JFLAG  = integer flag for input and output communication.
C
C          On input, set JFLAG = 0 on the first call for the problem,
C          and leave it unchanged until the problem is completed.
C          (The problem is completed when JFLAG .ge. 2 on return.)
C
C          On output, JFLAG has the following values and meanings:
C          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
C                    and call DROOTS again.
C          JFLAG = 2 means a root has been found.  The root is
C                    at X, and GX contains g(X).  (Actually, X is the
C                    rightmost approximation to the root on an interval
C                    (X0,X1) of size HMIN or less.)
C          JFLAG = 3 means X = X1 is a root, with one or more of the gi
C                    being zero at X1 and no sign changes in (X0,X1).
C                    GX contains g(X) on output.
C          JFLAG = 4 means no roots (of odd multiplicity) were
C                    found in (X0,X1) (no sign changes).
C
C X0,X1  = endpoints of the interval where roots are sought.
C          X1 and X0 are input when JFLAG = 0 (first call), and
C          must be left unchanged between calls until the problem is
C          completed.  X0 and X1 must be distinct, but X1 - X0 may be
C          of either sign.  However, the notion of -left- and -right-
C          will be used to mean nearer to X0 or X1, respectively.
C          When JFLAG .ge. 2 on return, X0 and X1 are output, and
C          are the endpoints of the relevant interval.
C
C G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
C          respectively.  When JFLAG = 0, G0 and G1 are input and
C          none of the G0(i) should be zero.
C          When JFLAG .ge. 2 on return, G0 and G1 are output.
C
C GX     = array of length NG containing g(X).  GX is input
C          when JFLAG = 1, and output when JFLAG .ge. 2.
C
C X      = independent variable value.  Output only.
C          When JFLAG = 1 on output, X is the point at which g(x)
C          is to be evaluated and loaded into GX.
C          When JFLAG = 2 or 3, X is the root.
C          When JFLAG = 4, X is the right endpoint of the interval, X1.
C
C JROOT  = integer array of length NG.  Output only.
C          When JFLAG = 2 or 3, JROOT indicates which components
C          of g(x) have a root at X.  JROOT(i) is 1 if the i-th
C          component has a root, and JROOT(i) = 0 otherwise.
C-----------------------------------------------------------------------
*/
	void DROOTS(int *, double *, int *, double *, double *, double *, double *, double *, double *, int *);

/*
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, DLSR01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 245 or more.
C ISAV = integer array of length 55 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
*/
	void DSRCAR(double *, int *, int *);

/*
C-----------------------------------------------------------------------
C DSTODPK performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C-----------------------------------------------------------------------
C The following changes were made to generate Subroutine DSTODPK
C from Subroutine DSTODE:
C 1. The array SAVX was added to the call sequence.
C 2. PJAC and SLVS were replaced by PSOL in the call sequence.
C 3. The Common block /DLPK01/ was added for communication.
C 4. The test constant EPCON is loaded into Common below statement
C    numbers 125 and 155, and used below statement 400.
C 5. The Newton iteration counter MNEWT is set below 220 and 400.
C 6. The call to PJAC was replaced with a call to DPKSET (fixed name),
C    with a longer call sequence, called depending on JACFLG.
C 7. The corrector residual is stored in SAVX (not Y) at 360,
C    and the solution vector is in SAVX in the 380 loop.
C 8. SLVS was renamed DSOLPK and includes NEQ, SAVX, EWT, F, and JAC.
C    SAVX was added because DSOLPK now needs Y and SAVF undisturbed.
C 9. The nonlinear convergence failure count NCFN is set at 430.
C-----------------------------------------------------------------------
C Note: DSTODPK is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODPK is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C          Also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C SAVX   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C CCMAX  = maximum relative change in H*EL0 before DPKSET is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in DPKSET or DSOLPK.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between DPKSET calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
*/

#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	void DSTODPK(int *, double *, double *, int *, double *, double *, double *,
				 double *, double *, double *, int *, F,JAC,PSOL);

#else
	
	void DSTODPK(int *, double *, double *, int *, double *, double *, double *,
		double *, double *, double *, int *, void(*)(int *, double *, double *, double *),
		void (*)(int *, double *, double *, int *, int *, double *, int *),
		void (*) (int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int *) );

#endif

/*
C-----------------------------------------------------------------------
C DPKSET is called by DSTODPK to interface with the user-supplied
C routine JAC, to compute and process relevant parts of
C the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
C as need for preconditioning matrix operations later.
C
C In addition to variables described previously, communication
C with DPKSET uses the following:
C Y     = array containing predicted values on entry.
C YSV   = array containing predicted y, to be saved (YH1 in DSTODPK).
C FTEM  = work array of length N (ACOR in DSTODPK).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.
C         Space for preconditioning data starts at WM(LOCWP).
C IWM   = integer work space.
C         Space for preconditioning data starts at IWM(LOCIWP).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         JAC returned an error flag.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
C-----------------------------------------------------------------------
*/
#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	void DPKSET(int *, double *, double *, double *, double *, double *, double *, int *, F, JAC);

#else

	void DPKSET(int *, double *, double *, double *, double *, double *, double *, int *,
				void(*)(int *, double *, double *, double *),
				void(*)(int *, double *, double *, int *, int *, double *, int *) );

#endif

/*
C-----------------------------------------------------------------------
C This routine interfaces to one of DSPIOM, DSPIGMR, DPCG, DPCGS, or
C DUSOL, for the solution of the linear system arising from a Newton
C iteration.  It is called if MITER .ne. 0.
C In addition to variables described elsewhere,
C communication with DSOLPK uses the following variables:
C WM    = real work space containing data for the algorithm
C         (Krylov basis vectors, Hessenberg matrix, etc.)
C IWM   = integer work space containing data for the algorithm
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C IERSL = output flag (in Common):
C         IERSL =  0 means no trouble occurred.
C         IERSL =  1 means the iterative method failed to converge.
C                    If the preconditioner is out of date, the step
C                    is repeated with a new preconditioner.
C                    Otherwise, the stepsize is reduced (forcing a
C                    new evaluation of the preconditioner) and the
C                    step is repeated.
C         IERSL = -1 means there was a nonrecoverable error in the
C                    iterative solver, and an error exit occurs.
C This routine also uses the Common variables TN, EL0, H, N, MITER,
C DELT, EPCON, SQRTN, RSQRTN, MAXL, KMP, MNEWT, NNI, NLI, NPS, NCFL,
C LOCWP, LOCIWP.
C-----------------------------------------------------------------------
*/
#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	void DSOLPK(int *, double *, double *, double *, double *, double *, int *, F, PSOL);

#else
	
	void DSOLPK(int *, double *, double *, double *, double *, double *, int *, 
				void(*)(int *, double *, double *, double *),
				void(*)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int *) );

#endif

/*
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using a scaled
C preconditioned version of the Incomplete Orthogonalization Method.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C         B    = the right hand side of the system A*x = b.
C                B is also used as work space when computing the
C                final approximation.
C                (B is the same as V(*,MAXL+1) in the call to DSPIOM.)
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C         N    = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, B, WGHT, and X.
C
C         MAXL = the maximum allowable order of the matrix HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array of length N used by DATV and PSOL.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         V    = the N by (LIOM+1) array containing the LIOM
C                orthogonal vectors V(*,1) to V(*,LIOM).
C
C         HES  = the LU factorization of the LIOM by LIOM upper
C                Hessenberg matrix whose entries are the
C                scaled inner products of A*V(*,k) and V(*,i).
C
C         IPVT = an integer array containg pivoting information.
C                It is loaded in DHEFA and used in DHESL.
C
C         LIOM = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LIOM iterations, LIOM.le.MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
*/

#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	void DSPIOM(int *, double *, double *, double *, double *, int *, int *, int *,
			    double *, double *, int *, int *, F, PSOL, int *, double *, double *,
				double *, int *, int *, double *, int *, double *, int * );

#else
	
	void DSPIOM(int *, double *, double *, double *, double *, int *, int *, int *,
		double *, double *, int *, int *, void(*)(int *, double *, double *, double *),
		void (*)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int *),
		int *, double *, double *, double *, int *, int *, double *, int *, double *, int * );

#endif

/*
	C-----------------------------------------------------------------------
C This routine computes the product
C
C   (D-inverse)*(P1-inverse)*(I - hl0*df/dy)*(P2-inverse)*(D*v),
C
C where D is a diagonal scaling matrix, and P1 and P2 are the
C left and right preconditioning matrices, respectively.
C v is assumed to have WRMS norm equal to 1.
C The product is stored in z.  This is computed by a
C difference quotient, a call to F, and two calls to PSOL.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            V = real array of length N (can be the same array as Z).
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the matrix D.
C
C         FTEM = work array of length N.
C
C         VTEM = work array of length N used to store the
C                unscaled version of V.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C
C      On return
C
C            Z = array of length N containing desired scaled
C                matrix-vector product.
C
C          IER = error flag from PSOL.
C
C         NPSL = the number of calls to PSOL.
C
C In addition, this routine uses the Common variables TN, N, NFE.
C-----------------------------------------------------------------------
*/
#if (C_WRAPPER_ODEPACK_FPTR_WORKAROUND) == 1
	
	void DATV(int *, double *, double *, double *, double *, F, PSOL, double *,
			  double *, double *, int *, double *, int *, int *, int *);

#else

	void DATV(int *, double *, double *, double *, double *,
			  void(*) (int *, double *, double *, double *),
			  void(*)(int *, double *, double *, double *, double *, double *, double *, int *, double *, int *, int *),
			  double *, double *, double *, int *, double *, int *, int *, int * );
	
#endif

/*
	C-----------------------------------------------------------------------
C This routine orthogonalizes the vector VNEW against the previous
C KMP vectors in the V array.  It uses a modified Gram-Schmidt
C orthogonalization procedure with conditional reorthogonalization.
C This is the version of 28 may 1986.
C-----------------------------------------------------------------------
C
C      On entry
C
C         VNEW = the vector of length N containing a scaled product
C                of the Jacobian and the vector V(*,LL).
C
C         V    = the N x l array containing the previous LL
C                orthogonal vectors v(*,1) to v(*,LL).
C
C         HES  = an LL x LL upper Hessenberg matrix containing,
C                in HES(i,k), k.lt.LL, scaled inner products of
C                A*V(*,k) and V(*,i).
C
C        LDHES = the leading dimension of the HES array.
C
C         N    = the order of the matrix A, and the length of VNEW.
C
C         LL   = the current order of the matrix HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to (KMP .le. MAXL).
C
C
C      On return
C
C         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
C                where i0 = MAX(1, LL-KMP+1).
C
C         HES  = upper Hessenberg matrix with column LL filled in with
C                scaled inner products of A*V(*,LL) and V(*,i).
C
C       SNORMW = L-2 norm of VNEW.
C
C-----------------------------------------------------------------------
*/
	void DORTHOG(double *, double *, double *, int *, int *, int *, int *, double *);

/*
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using a scaled
C preconditioned version of the Generalized Minimal Residual method.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C                B is also used as work space when computing
C                the final approximation.
C                (B is the same as V(*,MAXL+1) in the call to DSPIGMR.)
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, B and X.
C
C         MAXL = the maximum allowable order of the matrix HES.
C
C       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATV and PSOL.
C
C           DL = real work array used for calculation of the residual
C                norm RHO when the method is incomplete (KMP .lt. MAXL).
C                Not needed or referenced in complete case (KMP = MAXL).
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LGMR = the number of iterations performed and
C                the current order of the upper Hessenberg
C                matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C         V    = the N by (LGMR+1) array containing the LGMR
C                orthogonal vectors V(*,1) to V(*,LGMR).
C
C         HES  = the upper triangular factor of the QR decomposition
C                of the (LGMR+1) by lgmr upper Hessenberg matrix whose
C                entries are the scaled inner-products of A*V(*,i)
C                and V(*,k).
C
C         Q    = real array of length 2*MAXL containing the components
C                of the Givens rotations used in the QR decomposition
C                of HES.  It is loaded in DHEQR and used in DHELS.
C
C        IFLAG = integer error flag:
C                0 means convergence in LGMR iterations, LGMR .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
*/


} // End of extern "C"
#endif


#endif /*__GMS_ODEPACK_F77_H__*/
