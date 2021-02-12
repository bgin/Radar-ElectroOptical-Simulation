      LOGICAL FUNCTION CBIGP ( IP, DIM1, DIM2 )
*     .. Scalar Arguments ..
      INTEGER                IP, DIM1, DIM2
*     ..
*
*  Purpose
*  =======
*
*  CBIGP determines which of two alternative code sections in a GEMM-
*  Based Level 3 BLAS routine that will be the fastest for a particular
*  problem. If the problem is considered large enough CBIGP returns
*  .TRUE., otherwise .FALSE. is returned. The input parameter IP
*  specifies the calling routine and a break point for alternative code
*  sections. The input parameters DIM1 and DIM2 are matrix dimensions.
*  The returned value is a function of the input parameters and the
*  performance characteristics of the two alternative code sections.
*
*  In this simple implementation, the returned values are determined by
*  looking at only one of the two dimensions DIM1 and DIM2. It may be
*  rewarding to rewrite the logical expressions in CBIGP so that both
*  dimensions are involved. The returned values should effectively
*  reflect the performance characteristics of the underlying BLAS
*  routines.
*
*
*  Input
*  =====
*
*  IP     - INTEGER
*           On entry, IP specifies which routine and which alternative
*           code sections that the decision is intended for.
*           Unchanged on exit.
*
*  DIM1   - INTEGER.
*           On entry, DIM1 specifies the first dimension in the calling
*           sequence of the Level 3 routine specified by IP.
*           Unchanged on exit.
*
*  DIM2   - INTEGER.
*           On entry, DIM2 specifies the second dimension in the
*           calling sequence of the Level 3 routine specified by IP.
*           Unchanged on exit.
*
*
*  -- Written in May-1994.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. User specified parameters for CBIGP ..
      INTEGER            CIP41, CIP42,
     $                   CIP51, CIP52,
     $                   CIP81, CIP82, CIP83,
     $                   CIP91, CIP92, CIP93
      PARAMETER        ( CIP41 = 4, CIP42 = 3,
     $                   CIP51 = 4, CIP52 = 3,
     $                   CIP81 = 4, CIP82 = 3, CIP83 = 4,
     $                   CIP91 = 4, CIP92 = 3, CIP93 = 4 )
*     ..
*     .. Executable Statements ..
      IF( IP.EQ.41 )THEN
         CBIGP = DIM1.GE.CIP41
      ELSE IF( IP.EQ.42 )THEN
         CBIGP = DIM2.GE.CIP42
      ELSE IF( IP.EQ.51 )THEN
         CBIGP = DIM1.GE.CIP51
      ELSE IF( IP.EQ.52 )THEN
         CBIGP = DIM2.GE.CIP52
      ELSE IF( IP.EQ.81 )THEN
         CBIGP = DIM2.GE.CIP81
      ELSE IF( IP.EQ.82 )THEN
         CBIGP = DIM2.GE.CIP82
      ELSE IF( IP.EQ.83 )THEN
         CBIGP = DIM1.GE.CIP83
      ELSE IF( IP.EQ.91 )THEN
         CBIGP = DIM2.GE.CIP91
      ELSE IF( IP.EQ.92 )THEN
         CBIGP = DIM2.GE.CIP92
      ELSE IF( IP.EQ.93 )THEN
         CBIGP = DIM1.GE.CIP93
      ELSE
         CBIGP = .FALSE.
      END IF
*
      RETURN
*
*     End of CBIGP.
*
      END
