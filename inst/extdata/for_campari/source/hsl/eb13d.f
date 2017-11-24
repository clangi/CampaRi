! COPYRIGHT (c) 1993 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 1.1.1
! See ChangeLog for version history
!
      SUBROUTINE EB13ID(ICNTL,CNTL)
C
C  This subroutine gives the control parameters default values
C
C *ICNTL    _ Integer array of length 11.
C  ICNTL(1) - Stream number for errors. Default 6.
C  ICNTL(2) - Stream number for warnings. Default 6.
C  ICNTL(3) - Max. degree of iteration polynomial. Default 800.
C  ICNTL(4) - Controls choice of degree of iteration polynomial.
C             Default 0. If ICNTL(4).GE.0, the degree of the
C             iteration polynomial is taken to be ICNTL(4).
C  ICNTL(5) - Max. number of matrix-vector products allowed is
C             ICNTL(5)*NUMEIG. Default 20000.
C  ICNTL(6) - Controls amount of printing.
C             Default 2.
C             ICNTL(6) = 0 no printing.
C                      = 1 errors only.
C                      = 2 errors and warnings only.
C                      = 3 parameters on input and INFO and RINFO
C                          on final successful exit.
C                      = 4 converged eigenvalues on final exit.
C                      = 5 INFO and RINFO and the
C                          computed residual of first unconverged
C                          eigenvalue printed on each iteration.
C                      = 6 computed eigenvalues printed on
C                          each iteration.
C  ICNTL(7) - Default 0.  In this case, the user must set
C             KEEP(1) = ||A|| prior to the first call to
C             EB13A/AD and the stopping criteria will use ||A||
C             as supplied by the user in the stopping criteria.
C             If ICNTL(7) = 1, the Frobenius
C             norm of A will be computed by EB13A/AD
C             (requires N matrix-vector products) and this will
C             be used in the denominator of the stopping criteria.
C             If ICNTL(7) = 2, the stopping criteria will use
C             ||Ax|| in the denominator (norm of A not needed).
C             If ICNTL(7) = 3, the stopping criteria will use
C             the norm of the residual only (norm A not needed).
C  ICNTL(8) _ Default 1.
C             ICNTL(8) = 1 Ho ellipse
C                      = 2 Braconnier ellipse
C                        3 Saad ellipse
C  ICNTL(9) _ Default 2.
C             ICNTL(9) = 1 Algorithm A1/AB1
C                         (algorithm A1 if blocksize is 1
C                          and algorithm AB1 otherwise)
C                      = 2 Algorithm A2/AB2
C                         (algorithm A2 if blocksize is 1
C                          and algorithm AB2 otherwise)
C                      = 3 Algorithm A3/AB3
C                         (algorithm A3 if blocksize is 1
C                          and algorithm A31 otherwise)
C
C  ICNTL(10) _ Default 0. Serves the dual function of allowing the user
C              to supply an initial estimate of the basis vectors and
C              allowing the user to supply the initial seed for the
C              random number generator FA14.  The following values are
C              handled:
C              ICNTL(10) = 0 No basis vectors or seed supplied,
C                        = 1 Only the basis vectors are supplied,
C                        = 2 Only the seed is supplied,
C                        = 3 Both the basis vectors and the seed are
C                            supplied.
C
C  ICNTL(11) - Default 500. Maximum number of Arnoldi iterations
C              allowed is ICNTL(11)*NUMEIG.
C
C  CNTL(1)  - Convergence tolerance. Default u*10**3, u = machine
C             precision.
C
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(1)
      INTEGER ICNTL(11)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 80
      ICNTL(4) = 0
      ICNTL(5) = 20000
      ICNTL(6) = 2
      ICNTL(7) = 0
      ICNTL(8) = 1
      ICNTL(9) = 2
      ICNTL(10) = 0
      ICNTL(11) = 500
C
      CNTL(1) = FD15AD('E')*10**3
C
      RETURN
      END
C**************************************************************
      SUBROUTINE EB13AD(IND,N,NUMEIG,NBLOCK,NSTEPS,ER,EI,LN,X,U,W,IKEEP,
     +                  RKEEP,ICNTL,CNTL,INFO,RINFO)
C
C     This routine computes the right-most eigenvalues
C     (i.e. those with largest real parts)
C     or the eigenvalues of largest imaginary part and a basis for
C     the associated invariant subspace of a real
C     unsymmetric matrix A(). A Chebychev preconditoned Arnoldi
C     method or a block Arnoldi method with Chebychev
C     acceleration is used.
C
C     Arguments (an * indicates that the arguement is altered by the
C     routine).
C
C     IND       Integer variable. Must be set to 0 if the eigenvalues
C               of largest absolute value are required, to
C               1 if the right-most eigenvalues are required and
C               to 2 if the eigenvalues of largest imaginary parts
C               are required.
C     N         Integer variable. Must be set to the order of the
C               matrix (N.gt.1).
C     NUMEIG    Integer variable. On entry,  must be set to
C               the number of eigenvalues the user requires.
C     NBLOCK    Integer variable. On entry,  must be set to
C               the block size.
C     NSTEPS    Integer variable which must be set to
C               the number of Arnoldi steps per iteration
C     LN        An integer which holds the first dimension of the arrays
C               X, U, and W. LN must be at least N.
C   * ER(NSTEPS*NBLOCK) Real (DP) array.
C               On exit with positive IKEEP(1),
C               ER holds successive approximations to the real
C               parts of the sought-after eigenvalues of A.
C   * EI(NSTEPS*NBLOCK) Real (DP) array.
C               On exit with positive IKEEP(1),
C               EI holds successive approximations to the imaginary
C               parts of the sought-after eigenvalues of A.
C               Complex conjugate pairs of eigenvalues appear
C               consecutively with the eigenvalue with positive
C               imaginary part appearing first.
C   * X(LN,NSTEPS*NBLOCK) Real (DP) array.
C               Need not be set on first call.
C               On each exit, with positive IKEEP(1),
C               it holds approximations to the basis vectors
C               which span the space corresponding
C               to the right most (or left most) eigenvalues.
C   * U(LN,NSTEPS*NBLOCK) Real (DP) array used by the routine as a
C               work array.The user forms matrix-matrix multiplications
C               in this array.
C   * W(LN,NSTEPS*NBLOCK)   Real (DP) array. On each exit, the array
C               is used by the user for performing matrix-matrix
C               multiplications.
C   * IKEEP     Integer array of length at least 11+NSTEPS*NBLOCK
C     IKEEP(1)  Holds a pointer IPOS (determines point in code at
C               which control was returned to the user).
C               On first call, user must set IKEEP(1) = 0.
C               If IKEEP(1) = 0 on return, computation has terminated.
C               Otherwise, user must perform matrix-vector products.
C    IKEEP(2)   On each exit with IKEEP(1) greater then zero,
C               IKEEP(2) indicates the first column of W which must be
C               premultiplied by A.
C    IKEEP(3)   On each exit with IKEEP(1) greater then zero,
C               IKEEP(3) indicates the last column of W which must be
C               premultiplied by A.
C               The user should form the matrix product
C               U(IKEEP(2):IKEEP(3)) = A() * W(IKEEP(2):IKEEP(3))
C               (i.e. the product only involves columns IKEEP(2)
C               to IKEEP(3))
C     IKEEP(4)-IKEEP(8)  Used in dividing up the real (DP) workspace.
C     IKEEP(9)  Holds FA14's generator word (used in EB13C/CD).
C     IKEEP(10)  Holds number of points defining complex hull.
C     IKEEP(11)  Used to hold a pointer for iterations.
C     IKEEP(12)-IKEEP(11+NSTEPS*NBLOCK) used as workspace.
C   * RKEEP     Real (DP) array of length
C               NSTEPS*NBLOCK*(10+2*NSTEPS*NBLOCK)+13.
C               If ICNTL(7)=0, RKEEP(1) must be set by the user to
C               hold the norm of A.
C               If ICNTL(7)=1, RKEEP(1) need not be set
C               by the user and must not be altered between calls.
C               If ICNTL(7)=2, RKEEP(1) is set to 0.
C               If ICNTL(7)=3, RKEEP(1) is set to 1.
C     RKEEP(2)-RKEEP(6) Holds residuals for first unconverged eigenvalue
C               for current and previous 4 iterations. Initialised
C               to FD15A/AD('H')
C     RKEEP(7)  Holds FD15A/AD('E')
C     RKEEP(8)-RKEEP(12) Used to hold values for the Chebychev
C               three-term recurrence.
C     RKEEP(13)-RKEEP(12+(NSTEPS*NBLOCK)**2) used to hold Schur form.
C               The remainder of RKEEP is used as workspace.
C     ICNTL     Integer array length 11. Holds control parameters
C               (see EB13I/ID).
C     CNTL      Real (DP) array length 1. Holds control parameters
C               (see EB13I/ID).
C   * INFO      Integer array not set on first entry. On each exit
C               provides information.
C     INFO(1)   Negative value indicates fatal error;
C               positive value (>0) indicates warning message.
C     INFO(2)   Holds number of matrix-vector products formed so far.
C     INFO(3)   Holds current degree of iteration polynomial.
C     INFO(4)   Holds number of eigenvalues which have converged
C               so far.
C     INFO(5)   Holds number of iterations so far.
C     INFO(6)   Holds highest degree of polynomial used so far.
C   * RINFO     Real (DP) array not set on first entry. On each
C               exit provides information.
C     RINFO(1), RINFO(2), RINFO(3) Hold the current ellipse parameters
C               D,C2,A2.
C     RINFO(4)  CNTL(1) is used to determine if the computed eigenvalues
C               are sufficiently accurate. If CNTL(1) is out
C               of range, the default value u*10**3 (u = machine
C               precision) is used and is held in RINFO(4).
C               On exit, RINFO(4) is set to the value of the convergence
C               parameter actually satisfied by the computed
C               eigenvalues.
C     RINFO(5)  Used to hold Frobenius norm of the matrix.
C
C  Local variables
C  I    - DO LOOP variable.
C  IPOS - IKEEP(1). Used in reverse communication with the
C         user, to indicate the current status.
C  LP   - ICNTL(1). Stream number for errors.
C  MP   - ICNTL(2). Stream number for warnings/diagnostics.
C  NLOW - IKEEP(2). On each exit with IPOS greater then zero,
C         NLOW indicates the first column of W which must be
C         premultiplied by A.
C  NUP  - IKEEP(3). On each exit with IPOS greater then zero,
C         NUP indicates the last column of W which must be
C         premultiplied by A. The user should form the matrix
C         product U(NLOW:NUP) = A() * W(NLOW:NUP)
C         (i.e. the product only involves columns NLOW
C         to NUP)
C NUMCOL - Integer variable. Used to hold NSTEPS*NBLOCK
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IND,LN,N,NBLOCK,NSTEPS,NUMEIG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(1),EI(NSTEPS*NBLOCK),ER(NSTEPS*NBLOCK),
     +                 RKEEP(NSTEPS*NBLOCK* (10+2*NSTEPS*NBLOCK)+13),
     +                 RINFO(5),U(LN,NSTEPS*NBLOCK),W(LN,NSTEPS*NBLOCK),
     +                 X(LN,NSTEPS*NBLOCK)
      INTEGER ICNTL(11),IKEEP(11+NSTEPS*NBLOCK),INFO(6)
C     ..
C     .. Local Scalars ..
      INTEGER I,IPOS,LP,MP,NLOW,NUMCOL,NUP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL EB13CD
C     ..
      LP = ICNTL(1)
      MP = ICNTL(2)
      NUMCOL = NSTEPS*NBLOCK
C
      IF (IKEEP(1).EQ.0) THEN
C Initial entry.
C Optional printing of input parameters.
         IF (MP.GT.0 .AND. ICNTL(6).GE.3) THEN
            WRITE (MP,FMT='(/A/A,I2,A,I6,A,I5,A,I5,A,I5,A,I7)')
     +        ' Entering EB13A/AD with',' IND =',IND,' N =',N,
     +        '  NUMEIG =',NUMEIG,'  NBLOCK =',NBLOCK,'  NSTEPS =',
     +        NSTEPS,'  LN =',LN
            WRITE (MP,FMT='(/A,2I3,2I4,I7,5I4,I5/A,1PD12.4)')
     +        ' ICNTL (1:11) =', (ICNTL(I),I=1,11),' CNTL = ',CNTL
         END IF
C Check the validity of the input data.
         IF (N.LT.3) GO TO 50
         IF (NBLOCK.LT.1) GO TO 60
         IF (NUMEIG.LT.1 .OR. NUMEIG.GT.N-2) GO TO 70
         IF (NUMCOL.LT.NUMEIG+2 .OR. NUMCOL.GT.N .OR.
     +       NSTEPS.LE.1) GO TO 80
         IF (N.GT.LN) GO TO 90
C
C Initialise.
         DO 10 I = 1,6
            INFO(I) = 0
   10    CONTINUE
         INFO(5) = 1
         IF (IND.NE.0 .AND. ICNTL(9).NE.1) THEN
            IF (ICNTL(4).GT.0) INFO(3) = ICNTL(4)
         END IF
         DO 20 I = 1,5
            RINFO(I) = ZERO
   20    CONTINUE
C Check value of CNTL(1)
         RINFO(4) = CNTL(1)
         IF (CNTL(1).LT.FD15AD('E') .OR. CNTL(1).GT.ONE) THEN
            INFO(1) = 1
            RINFO(4) = FD15AD('E')*10**3
            IF (MP.GT.0 .AND. ICNTL(6).GE.2) THEN
               WRITE (MP,FMT=9070) INFO(1)
               WRITE (MP,FMT=9080) CNTL(1),RINFO(4)
            END IF
         END IF
         IF (ICNTL(7).EQ.1 .OR. ICNTL(7).EQ.2) RKEEP(1) = ZERO
         IF (ICNTL(7).EQ.3) THEN
            RKEEP(1) = ONE
            RINFO(5) = ONE
         END IF
         DO 30 I = 2,6
            RKEEP(I) = FD15AD('H')
   30    CONTINUE
         RKEEP(7) = FD15AD('E')
         IKEEP(2) = 0
         IKEEP(3) = 0
C Divide up rest of RKEEP.
         IKEEP(4) = 13
         IKEEP(5) = IKEEP(4) + NUMCOL*NUMCOL
         IKEEP(6) = IKEEP(5) + NUMCOL*NUMCOL
         IKEEP(7) = IKEEP(6) + 4*NUMCOL
         IKEEP(8) = IKEEP(7) + 4*NUMCOL
         IKEEP(10) = 0
         IKEEP(11) = 1
C
         DO 40 I = 1,NUMCOL
            ER(I) = ZERO
            EI(I) = ZERO
   40    CONTINUE
C
C Initialise FA14's generator word
         IF (ICNTL(10).NE.2.AND.ICNTL(10).NE.3) THEN
           CALL FA14ID(IKEEP(9))
         END IF
C
      END IF
C
      IPOS = IKEEP(1)
      NLOW = IKEEP(2)
      NUP = IKEEP(3)
      CALL EB13CD(IPOS,IND,N,NUMEIG,NBLOCK,NUMCOL,X,U,W,LN,NLOW,NUP,
     +            RKEEP(IKEEP(4)),RKEEP(IKEEP(5)),ER,EI,IKEEP(12),
     +            RKEEP(IKEEP(6)),RKEEP(IKEEP(7)),RKEEP(IKEEP(8)),
     +            RKEEP(IKEEP(8)+NUMCOL),ICNTL,INFO,RINFO,IKEEP,
     +            RKEEP,CNTL)
      IKEEP(1) = IPOS
      IKEEP(2) = NLOW
      IKEEP(3) = NUP
      IF (INFO(1).EQ.0) GO TO 140
      IF (INFO(1).EQ.-6) GO TO 100
      IF (INFO(1).EQ.-7) GO TO 110
      IF (INFO(1).EQ.-8) GO TO 120
      IF (INFO(1).EQ.-9) GO TO 130
      GO TO 140
C
C  Error returns
   50 INFO(1) = -1
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) N
      END IF
      GO TO 140
   60 INFO(1) = -2
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9020) NBLOCK
      END IF
      GO TO 140
   70 INFO(1) = -3
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9030) NUMEIG
      END IF
      GO TO 140
   80 INFO(1) = -4
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9050) NSTEPS
      END IF
      GO TO 140
   90 INFO(1) = -5
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9090) N
      END IF
      GO TO 140
  100 INFO(1) = -6
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9100) RKEEP(1)
      END IF
      GO TO 140
  110 INFO(1) = -7
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9040)
      END IF
      GO TO 140
  120 INFO(1) = -8
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9060) NUMEIG,INFO(4)
      END IF
      GO TO 140
  130 INFO(1) = -9
      IF (LP.GT.0 .AND. ICNTL(6).GE.1) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9110)
      END IF
C
  140 RETURN
 9000 FORMAT (/'  Error return from EB13A/AD. INFO(1) = ',I2)
 9010 FORMAT ('  Value of N out of range. N = ',I6)
 9020 FORMAT ('  Value of NBLOCK out of range. NBLOCK = ',I6)
 9030 FORMAT ('  Value of NUMEIG out of range. NUMEIG = ',I6)
 9040 FORMAT ('  Number of matrix-vector multiplications required is ',
     +       /'  more than was specified by ICNTL(5)')
 9050 FORMAT ('  Value of NSTEPS out of range. NSTEPS = ',I6)
 9060 FORMAT ('  The algorithm has failed to find all ',I10,/'  eigen',
     +       'values requested by the user. ',/I10,'   eigenvalues ',
     +       'have been found')
 9070 FORMAT (/'  Warning from EB13A/AD.  INFO(1) = ',I2)
 9080 FORMAT (/' Control parameter CNTL(1) = ',D12.5,' out of ',
     +       'range.',/' Code will use the default convergence ',
     +       'parameter ',1PD12.5)
 9090 FORMAT ('  Increase LN to at least = ',I6)
 9100 FORMAT ('  Norm of A is = ',1PD12.5)
 9110 FORMAT ('  Number of Arnoldi iterations required is more',/,
     +        '  than was specified by ICNTL(11)')
      END
C****************************************************************
      SUBROUTINE EB13BD(N,NBLOCK,NSTEPS,ER,EI,LN,X,U,W,NEV,Y,RES,IKEEP,
     +                  RKEEP)
C
C     This routine computes the (right) eigenvectors corresponding
C     to the eigenvalues computed by EB13A/AD
C
C     Arguments (an * indicates that the arguement is altered by the
C     routine).
C
C     N,  NBLOCK, NSTEPS Must be unchanged since EB13A/AD.
C     NEV        Number of eigenvalues computed by EB13A/AD
C                (=INFO(4) )
C   * Y(LN,NEV) Real (DP) array. Not set on entry. On exit,
C                with positive IPOS, this array
C                holds the computed eigenvectors. If EI(J)=0.0, the Jth
C                column of Y contains the corresponding eigenvector.
C                If EI(J).GT.0.0, then the Jth and (J+1)th eigenvalues
C                are a complex conjugate pair, and the Jth and (J+1)th
C                columns of Y contain, respectively, the real and
C                imaginary parts of the corresponding complex conjugate
C                eigenvectors.
C     LN         An integer which holds the first dimension of the
C                arrays X, U, and W. LN must be at least N.
C    X(LN,NSTEPS*NBLOCK) Real (DP) array. This array must be unchanged
C                since the last call to EB13A/AD.
C     ER(NSTEPS*NBLOCK) Real (DP) array. This array must be unchanged
C                since the last call to EB13A/AD.
C     EI(NSTEPS*NBLOCK) Real (DP) array. This array must be unchanged
C                since the last call to EB13A/AD.
C   * U(LN,NSTEPS*NBLOCK) Real (DP) array used by the routine as a work
C                array. The user forms matrix-matrix multiplications in
C                this array.
C   * W(LN,NSTEPS*NBLOCK)   Real (DP) array. On exit with IPOS=1,
C                the array is used by the
C                user for performing matrix-matrix multiplications.
C   * RES(NEV)   Real (DP) array. Not set on entry.
C                On exit with IKEEP(1)=0,
C                holds the residuals of the computed eigensolutions.
C   * IKEEP(4)   Integer array of length 4.
C     IKEEP(1)   Used in reverse communication with the
C                user, to indicate the current status. On initial call
C                to this routine it should be set to 0. On return from
C                this routine it has the following meanings
C               = 0       routine has terminated normally
C               = 1       the eigenvectors have been computed but if the
C                         user wants the residuals to be computed
C                         the user should form the matrix product
C                         U() = A() * W(), cols IKEEP(2),...,IKEEP(3)
C                         and this routine should be
C                         re-called without further alteration to any of
C                         its arguments.
C                On exit with IKEEP(1) = 0, IF ICNTL(7)=2, IKEEP(4)
C                holds the number of eigenvectors for which
C                ||(AY) sub i|| < u, u=machine precision
C   * RKEEP      Real (DP) array.
C                The first 12+(NSTEPS*NBLOCK)*(NSTEPS*NBLOCK)
C                entries of the array RKEEP  must be unchanged
C                since the last call to EB13A/AD.
C
C     .. Scalar Arguments ..
      INTEGER LN,N,NBLOCK,NEV,NSTEPS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EI(NSTEPS*NBLOCK),ER(NSTEPS*NBLOCK),
     +                 RKEEP(12+ (NSTEPS*NBLOCK)**2),RES(NEV),
     +                 U(LN,NSTEPS*NBLOCK),W(LN,NSTEPS*NBLOCK),
     +                 X(LN,NSTEPS*NBLOCK),Y(LN,NEV)
      INTEGER IKEEP(4)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ANORM
      INTEGER NUMCOL
C     ..
C     .. External Subroutines ..
      EXTERNAL EB13QD,EB13RD
C     ..
      IKEEP(4) = 0
      NUMCOL = NSTEPS*NBLOCK
      IF (IKEEP(1).EQ.0) THEN
         CALL EB13RD(N,NUMCOL,NEV,ER,EI,Y,X,W,LN,RKEEP(13))
         IKEEP(1) = 1
         IKEEP(2) = 1
         IKEEP(3) = NEV
      ELSE IF (IKEEP(1).EQ.1) THEN
         ANORM = RKEEP(1)
         CALL EB13QD(N,NEV,LN,Y,U,RES,ER,EI,IKEEP(4),ANORM)
         IKEEP(1) = 0
      END IF
C
      RETURN
      END
C*************************************************************
      SUBROUTINE EB13CD(IPOS,IND,N,NUMEIG,NBLOCK,NUMCOL,X,U,W,LN,NLOW,
     +                  NUP,QTAQ,WORK,ER,EI,ITYPE,CH,EV,EROLD,EIOLD,
     +                  ICNTL,INFO,RINFO,IKEEP,KEEP,CNTL)
C
C     Arguments (an * indicates that the arguement is altered by the
C     routine).
C
C   * IPOS      Integer variable used in reverse communication with the
C               user, to indicate the current status.
C               On return from
C               this routine it has the following meanings
C               = 0       routine has terminated normally
C               > 0       the user should form the matrix product
C                         U(NLOW:NUP) = A() * W(NLOW:NUP)
C                         and this routine
C                         re-called without further alteration to any of
C                         its arguments.
C     IND       Integer variable. Must be set to 0 if the eigenvalues
C               of largest absolute value are required, to
C               1 if the right-most eigenvalues are required and
C               to 2 if the left-most eigenvalues are required.
C     NUMEIG    Integer variable. On entry, NUMEIG must be set to
C               the number of eigensolutions the user requires.
C     NUMCOL    Integer variable which must be set to
C               the (number of Arnoldi steps per iteration*NBLOCK).
C   * X(LN,NUMCOL) Real (DP) array. On each exit, with positive IPOS,
C               the first NUMEIG columns hold approximations to the
C               NUMEIG basis vectors which span the space corresponding
C               to the right most  eigenvalues.
C   * U(LN,NUMCOL) Real (DP) array used by the routine as a work array.
C               The user forms matrix-matrix multiplications in this
C               array.
C   * W(LN,NUMCOL)   Real (DP) array. On each exit, the array is used
C               by the user for performing matrix-matrix
C               multiplications.
C     LN        An integer which holds the first dimension of the arrays
C               X, U, and W. LN must be at least N.
C   * NLOW      Integer variable. Not set on entry. On each exit
C               with IPOS greater then zero,
C               NLOW indicates the first column of W which must be
C               premultiplied by A.
C   * NUP       Integer variable. Not set on entry. On each exit
C               with IPOS greater then zero,
C               NUP indicates the last column of W which must be
C               premultiplied by A.
C               The user should form the matrix product
C               U(NLOW:NUP) = A() * W(NLOW:NUP)
C               (i.e. the product only involves columns NLOW
C               to NUP)
C   * WORK (NUMCOL,NUMCOL) Real (DP) array used as workspace.
C   * QTAQ(NUMCOL,NUMCOL)  Real (DP) array.
C               On each exit, holds X(T) * A * X.
C   * ER(NUMCOL) Real (DP) array. On exit with positive IPOS,
C               ER holds successive approximations to the real
C               parts of the sought-after eigenvalues of A.
C   * EI(NUMCOL) Real (DP) array. On exit with positive IPOS,
C               EI holds successive approximations to the imaginary
C               parts of the sought-after eigenvalues of A.
C   * ITYPE(NUMCOL) Integer array. On exit with positive IPOS, the
C               i-th entry of ITYPE is
C                 0  if the i-th eigenvalue is real
C                 1  if the i-th eigenvalue is complex
C                    with positive real part
C                 2  if the i-th eigenvalue is complex
C                    with negative real part
C                -1  if the i-th eigenvalue was not calculated
C                    successfully.
C   * CH(2,2*NUMCOL) Real (DP) array used as workspace
C   * EV(2,2*NUMCOL) Real (DP) array used as workspace
C   * EROLD(NUMCOL) Real (DP) array. On exit with positive IPOS,
C               EROLD holds approximations to the real
C               parts of the sought-after eigenvalues of A
C               found on previous iteration.
C   * EIOLD(NUMCOL) Real (DP) array. On exit with positive IPOS,
C               EIOLD holds approximations to the imaginary
C               parts of the sought-after eigenvalues of A
C               found on previous iteration.
C     ICNTL     Integer array length 11. Holds control parameters
C               (see EB13I/ID).
C   * INFO      Integer array not set on first entry. On each exit
C               provides information.
C     INFO(1)   Negative value indicates fatal error;
C               positive value (>0) indicates warning message.
C     INFO(2)   Holds number of matrix-vector products formed.
C     INFO(3)   Holds current degree of iteration polynomial.
C     INFO(4)   Holds number of eigenvalues which have converged
C               so far.
C     INFO(5)   Holds number of iterations so far.
C     INFO(6)   Holds highest degree of polynomial used so far.
C   * RINFO     Real (DP) array not set on first entry. On each
C               exit provides information.
C     RINFO(1), RINFO(2), RINFO(3) Hold the current ellipse parameters
C               D,C,A.
C     RINFO(4)  CNTL(1) is used to determine if the computed eigenvalues
C               are sufficiently accurate. If CNTL(1) is out
C               of range, the default value SQRT(u) (u = machine
C               precision) is used and is held in RINFO(4).
C               On exit, RINFO(4) is set to the value of the convergence
C               parameter actually satisfied by the computed
C               eigenvalues.
C     RINFO(5)  Used to hold Frobenius norm of the matrix.
C
C Local variables
C  IPNTR   Used in iterations to keep count.(IKEEP(11))
C  JFOUND  Holds the number of blocks of eigenvalues that
C          have converged (=number of converged eigenvalues
C          for blocksize 1).
C  MHULL   Holds number of points which define complex hull
C
C     .. Parameters ..
      DOUBLE PRECISION HALF,ONE,ZERO,TWO,TOL
      PARAMETER (HALF=0.5D0,ONE=1.0D0,ZERO=0.0D0,TWO=2.0D0,TOL=0.01D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IND,IPOS,LN,N,NBLOCK,NLOW,NUMCOL,NUMEIG,NUP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CH(2,2*NUMCOL),EI(NUMCOL),EIOLD(NUMCOL),
     +                 ER(NUMCOL),EROLD(NUMCOL),EV(2,2*NUMCOL),KEEP(12),
     +                 QTAQ(NUMCOL,NUMCOL),WORK(NUMCOL,NUMCOL),RINFO(5),
     +                 U(LN,NUMCOL),W(LN,NUMCOL),X(LN,NUMCOL),CNTL(1)
      INTEGER ICNTL(11),IKEEP(11),INFO(6),ITYPE(NUMCOL)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2,C2,D,RATIO,T,T1,T2,T3,WRT1,XI,XR,NETA
      INTEGER I,I1,IEV,IFOUND,IPNTR,J,J1,J2,JFOUND,JJ,JPNTR,JPOWER,K,
     +        ICNTL9,KSTOP,MHULL,MP,NEV,ITS
C     ..
C
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,FA14AD
      EXTERNAL DDOT,DNRM2,FA14AD
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMM,DGEMV,DSCAL,EB13DD,EB13ED,EB13FD,
     +         EB13GD,EB13HD,EB13JD,EB14AD,EB14BD,EB14CD,EB14DD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,LOG10,MAX,MIN,MOD,DBLE,SQRT
C     ..
   10 CONTINUE
      IFOUND = INFO(4)
      JFOUND = IFOUND/NBLOCK
      MP = ICNTL(2)
      ICNTL9 = ICNTL(9)
      IF (IND.EQ.0) ICNTL9 = 1
      IF (IPOS.EQ.0) THEN
C This is the initial call to this routine. If ICNTL(7)=1
C compute the Frobenius norm of A. This will be required when checking
C for convergence. It requires the user to form matrix vector products.
C The norm of A is held in KEEP(1) (and RINFO(5)).
C
         IF (ICNTL(7).EQ.1) THEN
            DO 20 I = 1,N
               W(I,1) = ZERO
   20       CONTINUE
            W(1,1) = ONE
            NLOW = 1
            NUP = 1
            INFO(2) = INFO(2) + 1
            IPOS = 1
            GO TO 290
         ELSE
            IPOS = 1
            GO TO 10
         END IF
C
      ELSE IF (IPOS.EQ.1) THEN
         IF (ICNTL(7).EQ.1) THEN
            DO 30 I = 1,N
               KEEP(1) = KEEP(1) + U(I,1)**2
               W(I,1) = ZERO
   30       CONTINUE
            IKEEP(11) = IKEEP(11) + 1
            IF (IKEEP(11).LE.N) THEN
               W(IKEEP(11),1) = ONE
               NLOW = 1
               NUP = 1
               INFO(2) = INFO(2) + 1
               GO TO 290
            ELSE
               KEEP(1) = SQRT(KEEP(1))
            END IF
         END IF
         RINFO(5) = KEEP(1)
C Check that the norm of A (in KEEP(1)) is large enough
         IF (ICNTL(7).NE.2 .AND. KEEP(1).LT.100.0D0*KEEP(7)) THEN
            INFO(1) = -6
            GO TO 290
         END IF
         IKEEP(11) = 0

         IF (ICNTL(10).NE.1.AND.ICNTL(10).NE.3) THEN
C Generate a random array X()
            DO 50 J = 1,NBLOCK
               DO 40 I = 1,N
                  X(I,J) = FA14AD(IKEEP(9),-1)
   40          CONTINUE
   50       CONTINUE
         END IF
         IF (NBLOCK.GT.1) THEN
            CALL EB13JD(X,LN,N,NBLOCK,0,U)
            CALL EB13ED(N,X,LN,W,LN,NBLOCK)
         ELSE
            T = DNRM2(N,X,1)
            CALL DSCAL(N,ONE/T,X,1)
            CALL DCOPY(N,X,1,W,1)
         END IF
         NLOW = 1
         NUP = NBLOCK
         INFO(2) = INFO(2) + NBLOCK
         IPOS = 2
         GO TO 290
C
      ELSE IF (IPOS.EQ.2) THEN
         IPNTR = IKEEP(11)
         IF (IPNTR.EQ.JFOUND*NBLOCK .AND. NBLOCK.GT.1) THEN
C Ensure we have lower part of matrix set to zero
C for block Hessenberg form
            DO 70 J = JFOUND*NBLOCK + 1,NUMCOL
               DO 60 I = JFOUND*NBLOCK + 1,NUMCOL
                  QTAQ(I,J) = ZERO
   60          CONTINUE
   70       CONTINUE
         END IF
         IPNTR = IPNTR + NBLOCK
         J1 = IPNTR + 1
         J = J1 - NBLOCK
         IF (NBLOCK.GT.1) THEN
            CALL EB13ED(N,U(1,J),LN,X(1,J1),LN,NBLOCK)
            NETA = ONE/ (SQRT(TWO))

            DO 80 I = 1,NBLOCK
               W(I,1) = DNRM2(N,X(1,J1+I-1),1)
   80       CONTINUE
            ITS = 0

   90       CALL DCOPY(NBLOCK,W(1,1),1,W(1,2),1)
            ITS = ITS + 1

            CALL DGEMM('T','N',IPNTR,NBLOCK,N,ONE,X,LN,X(1,J1),LN,ZERO,
     +                 U(1,J),LN)
            CALL DGEMM('N','N',N,NBLOCK,IPNTR,-ONE,X,LN,U(1,J),LN,ONE,
     +                 X(1,J1),LN)
            DO 100 I = 1,NBLOCK
               IF (ITS.EQ.1) CALL DCOPY(IPNTR,U(1,J+I-1),1,
     +                                  QTAQ(1,J+I-1),1)
               IF (ITS.GT.1) CALL DAXPY(IPNTR,ONE,U(1,J+I-1),1,
     +                                  QTAQ(1,J+I-1),1)
  100       CONTINUE

C Test to see if we need to iterate GS process
            DO 110 I = 1,NBLOCK
               W(I,1) = DNRM2(N,X(1,J1+I-1),1)
  110       CONTINUE
            DO 120 I = 1,NBLOCK
               IF (W(I,1).LT.NETA*W(I,2) .AND. ITS.LE.3) THEN
                  GO TO 90
               END IF
  120       CONTINUE

c            CALL DGEMM('T','N',IPNTR,NBLOCK,N,ONE,X,LN,X(1,J1),LN,
c     +                 ZERO,QTAQ(1,J),NUMCOL)
c            CALL DGEMM('N','N',N,NBLOCK,IPNTR,-ONE,X,LN,QTAQ(1,J),
c     +                 NUMCOL,ONE,X(1,J1),LN)
c            CALL DGEMM('T','N',IPNTR,NBLOCK,N,ONE,X,LN,X(1,J1),LN,
c     +                 ZERO,U(1,J),LN)
c            CALL DGEMM('N','N',N,NBLOCK,IPNTR,-ONE,X,LN,U(1,J),LN,
c     +                 ONE,X(1,J1),LN)
c            DO 80 I = 1,NBLOCK
c               CALL DAXPY(IPNTR,ONE,U(1,J+I-1),1,QTAQ(1,J+I-1),1)
c   80       CONTINUE
            DO 140 I = 1,NBLOCK
               T = DNRM2(N,X(1,J1+I-1),1)
               IF (T.LT.KEEP(7)) THEN
                  IF (MP.GT.0 .AND. ICNTL(6).GE.2 .AND.
     +                INFO(1).NE.3) WRITE (MP,FMT=9000)
                  INFO(1) = 3
                  T = KEEP(7)
               END IF
               QTAQ(J1+I-1,J+I-1) = T
               CALL DSCAL(N,ONE/T,X(1,J1+I-1),1)
               DO 130 I1 = I + 1,NBLOCK
                  T = DDOT(N,X(1,J1+I1-1),1,X(1,J1+I-1),1)
                  QTAQ(J1+I-1,J+I1-1) = T
                  CALL DAXPY(N,-T,X(1,J1+I-1),1,X(1,J1+I1-1),1)
  130          CONTINUE
  140       CONTINUE
            CALL EB13ED(N,X(1,J1),LN,W(1,J1),LN,NBLOCK)
            NLOW = J1
            NUP = J1 + NBLOCK - 1
            INFO(2) = INFO(2) + NBLOCK
            IKEEP(11) = IPNTR
            IF (IPNTR+NBLOCK.EQ.NUMCOL) THEN
               IPOS = 6
               IKEEP(11) = 0
               CALL EB13JD(X,LN,N,NUMCOL,JFOUND*NBLOCK,U)
            END IF
            GO TO 290

         ELSE

            CALL DCOPY(N,U(1,J),1,X(1,J1),1)
            NETA = ONE/ (SQRT(TWO))
            T = DNRM2(N,X(1,J1),1)
            ITS = 0
  150       T1 = T
            ITS = ITS + 1
C Use classical Gram Schmidt
            CALL DGEMV('T',N,J,ONE,X,LN,X(1,J1),1,ZERO,U(1,J),1)
            CALL DGEMV('N',N,J,-ONE,X,LN,U(1,J),1,ONE,X(1,J1),1)
            IF (ITS.EQ.1) CALL DCOPY(J,U(1,J),1,QTAQ(1,J),1)
            IF (ITS.GT.1) CALL DAXPY(J,ONE,U(1,J),1,QTAQ(1,J),1)
            T = DNRM2(N,X(1,J1),1)
C Test to see if we need to iterate GS process
            IF (T.LT.NETA*T1 .AND. ITS.LE.3) THEN
               GO TO 150
            END IF

            IF (T.LT.KEEP(7)) THEN
               IF (MP.GT.0 .AND. ICNTL(6).GE.2 .AND.
     +             INFO(1).NE.3) WRITE (MP,FMT=9000)
               INFO(1) = 3
               T = KEEP(7)
            END IF

            QTAQ(J1,J) = T

            CALL DSCAL(N,ONE/T,X(1,J1),1)

C Return for matrix-vector product
            CALL DCOPY(N,X(1,J1),1,W(1,J1),1)
            NLOW = J1
            NUP = J1
            INFO(2) = INFO(2) + 1
            IKEEP(11) = IPNTR
            IF (IPNTR+1.EQ.NUMCOL) THEN
               IPOS = 6
               IKEEP(11) = 0
            END IF
            GO TO 290
         END IF
C
      ELSE IF (IPOS.EQ.3) THEN
C
         IF (ICNTL9.EQ.2) THEN
            J1 = (JFOUND+1)*NBLOCK - IFOUND
C (If NBLOCK=1, JFOUND=IFOUND and J1=1)
C Arnoldi iteration with Chebychev acceleration
            IF (IKEEP(11).EQ.0) THEN
               IF (NBLOCK.GT.1) CALL EB13ED(N,X(1,IFOUND+1),LN,
     +                               W(1,IFOUND+1),LN,J1)
               IF (NBLOCK.EQ.1) CALL DCOPY(N,X(1,IFOUND+1),1,
     +                                     W(1,IFOUND+1),1)
               KEEP(9) = TWO/KEEP(8)
               KEEP(10) = ZERO
               T1 = TWO*KEEP(9)/KEEP(12)
            ELSE
               T1 = KEEP(9)/KEEP(12)
            END IF
            T2 = T1*RINFO(1)
            T3 = KEEP(9)*KEEP(10)*KEEP(11)
C
C Recurrence relation
            CALL EB13HD(N,T1,T2,T3,U(1,IFOUND+1),LN,W(1,IFOUND+1),LN,
     +                  X(1,IFOUND+1),LN,J1)
            IKEEP(11) = IKEEP(11) + 1
            IF (IKEEP(11).LT.INFO(3)) THEN
               IF (NBLOCK.GT.1) THEN
                  CALL EB13ED(N,W(1,IFOUND+1),LN,X(1,IFOUND+1),LN,J1)
                  CALL EB13ED(N,U(1,IFOUND+1),LN,W(1,IFOUND+1),LN,J1)
               ELSE
                  CALL DCOPY(N,W(1,IFOUND+1),1,X(1,IFOUND+1),1)
                  CALL DCOPY(N,U(1,IFOUND+1),1,W(1,IFOUND+1),1)
               END IF
               NLOW = IFOUND + 1
               NUP = (JFOUND+1)*NBLOCK
               INFO(2) = INFO(2) + NUP - IFOUND
               KEEP(10) = KEEP(9)
               KEEP(9) = ONE/ (KEEP(8)-KEEP(9)*KEEP(11))
               GO TO 290
            ELSE
               IF (NBLOCK.GT.1) THEN
                  CALL EB13ED(N,U(1,IFOUND+1),LN,X(1,IFOUND+1),LN,J1)
                  CALL EB13JD(X,LN,N, (JFOUND+1)*NBLOCK,IFOUND,U)
                  CALL EB13JD(X,LN,N, (JFOUND+1)*NBLOCK,IFOUND,U)
                  CALL EB13ED(N,X(1,JFOUND*NBLOCK+1),LN,
     +                        W(1,JFOUND*NBLOCK+1),LN,NBLOCK)
               ELSE
                  CALL DCOPY(N,U(1,IFOUND+1),1,X(1,IFOUND+1),1)
                  CALL EB13JD(X,LN,N,IFOUND+1,IFOUND,U)
                  CALL EB13JD(X,LN,N,IFOUND+1,IFOUND,U)
                  CALL DCOPY(N,X(1,IFOUND+1),1,W(1,IFOUND+1),1)
               END IF
               NLOW = JFOUND*NBLOCK + 1
               NUP = (JFOUND+1)*NBLOCK
               INFO(2) = INFO(2) + NBLOCK
               IKEEP(11) = JFOUND*NBLOCK
               IPOS = 2
               GO TO 290
            END IF
         ELSE
C Arnoldi without Chebychev
            IKEEP(11) = JFOUND*NBLOCK
            IPOS = 2
            GO TO 10
         END IF
C
      ELSE IF (IPOS.EQ.4) THEN
C
C  Preconditioned Arnoldi iteration with Chebychev
         IPNTR = IKEEP(11)
         I = MOD(IPNTR,INFO(3)) + 1
C JPNTR counts the number of blocks that have been processed
         JPNTR = (IPNTR/INFO(3)) + JFOUND
         J = JPNTR*NBLOCK + 1
         JJ = J + NBLOCK - 1
C  JJ = J if NBLOCK=1)
         IF (I.EQ.1) THEN
            IF (NBLOCK.GT.1) THEN
               CALL EB13ED(N,X(1,J),LN,W(1,JJ+1),LN,NBLOCK)
               CALL EB13ED(N,U(1,J),LN,U(1,JJ+1),LN,NBLOCK)
            ELSE
               CALL DCOPY(N,X(1,J),1,W(1,JJ+1),1)
               CALL DCOPY(N,U(1,J),1,U(1,JJ+1),1)
            END IF
            KEEP(9) = TWO/KEEP(8)
            KEEP(10) = ZERO
            T1 = TWO*KEEP(9)/KEEP(12)
         ELSE
            T1 = KEEP(9)/KEEP(12)
         END IF
         T2 = T1*RINFO(1)
         T3 = KEEP(11)*KEEP(9)*KEEP(10)
C
C Recurrence relation
         CALL EB13HD(N,T1,T2,T3,U(1,JJ+1),LN,W(1,JJ+1),LN,X(1,JJ+1),LN,
     +               NBLOCK)
C
         IPNTR = IPNTR + 1
         IF (I.LT.INFO(3)) THEN
            IF (NBLOCK.GT.1) THEN
               CALL EB13ED(N,W(1,JJ+1),LN,X(1,JJ+1),LN,NBLOCK)
               CALL EB13ED(N,U(1,JJ+1),LN,W(1,JJ+1),LN,NBLOCK)
            ELSE
               CALL DCOPY(N,W(1,JJ+1),1,X(1,JJ+1),1)
               CALL DCOPY(N,U(1,JJ+1),1,W(1,JJ+1),1)
            END IF
            NLOW = JJ + 1
            NUP = JJ + NBLOCK
            INFO(2) = INFO(2) + NBLOCK
            KEEP(10) = KEEP(9)
            KEEP(9) = ONE/ (KEEP(8)-KEEP(11)*KEEP(9))
         ELSE
            IF (NBLOCK.GT.1) THEN
               NETA = ONE/ (SQRT(TWO))
               CALL EB13ED(N,U(1,JJ+1),LN,X(1,JJ+1),LN,NBLOCK)
               CALL EB13JD(X,LN,N,JJ+NBLOCK,JJ,U)

               DO 160 I = 1,NBLOCK
                  W(I,1) = DNRM2(N,X(1,JJ+I-1),1)
  160          CONTINUE
               ITS = 0

  170          CALL DCOPY(NBLOCK,W(1,1),1,W(1,2),1)
               ITS = ITS + 1

               CALL DGEMM('T','N',JJ,NBLOCK,N,ONE,X,LN,X(1,JJ+1),LN,
     +                    ZERO,U,LN)
               CALL DGEMM('N','N',N,NBLOCK,JJ,-ONE,X,LN,U,LN,ONE,
     +                    X(1,JJ+1),LN)
C Test to see if we need to iterate GS process
               DO 180 I = 1,NBLOCK
                  W(I,1) = DNRM2(N,X(1,JJ+I),1)
  180          CONTINUE
               DO 190 I = 1,NBLOCK
                  IF (W(I,1).LT.NETA*W(I,2) .AND. ITS.LE.3) THEN
                     GO TO 170
                  END IF
  190          CONTINUE

C QR factorisation.
               DO 210 I = 1,NBLOCK
                  T = DNRM2(N,X(1,JJ+I),1)
                  IF (T.LT.KEEP(7)) THEN
                     IF (MP.GT.0 .AND. ICNTL(6).GE.2 .AND.
     +                   INFO(1).NE.3) WRITE (MP,FMT=9000)
                     INFO(1) = 3
                     T = KEEP(7)
                  END IF
                  CALL DSCAL(N,ONE/T,X(1,JJ+I),1)
                  DO 200 I1 = I + 1,NBLOCK
                     T = DDOT(N,X(1,JJ+I1),1,X(1,JJ+I),1)
                     CALL DAXPY(N,-T,X(1,JJ+I),1,X(1,JJ+I1),1)
  200             CONTINUE
  210          CONTINUE
               CALL EB13ED(N,X(1,JJ+1),LN,W(1,JJ+1),LN,NBLOCK)

            ELSE
               NETA = ONE/ (SQRT(TWO))
               CALL DCOPY(N,U(1,JJ+1),1,X(1,JJ+1),1)
               CALL EB13JD(X,LN,N,JJ+1,JJ,U)
               CALL DSCAL(N,ZERO,W,1)
               T = DNRM2(N,X(1,JJ+1),1)
               ITS = 0

  220          T1 = T
               ITS = ITS + 1
               CALL DGEMV('T',N,JJ,ONE,X,LN,X(1,JJ+1),1,ZERO,U,1)
               CALL DGEMV('N',N,JJ,-ONE,X,LN,U,1,ONE,X(1,JJ+1),1)
               T = DNRM2(N,X(1,JJ+1),1)

               IF (T.LT.NETA*T1 .AND. ITS.LE.3) THEN
                  GO TO 220
               END IF

               IF (T.LT.KEEP(7)) THEN
                  IF (MP.GT.0 .AND. ICNTL(6).GE.2 .AND.
     +                INFO(1).NE.3) WRITE (MP,FMT=9000)
                  INFO(1) = 3
                  T = KEEP(7)
               END IF

               CALL DSCAL(N,ONE/T,X(1,JJ+1),1)
               CALL DCOPY(N,X(1,JJ+1),1,W(1,JJ+1),1)

            END IF
            IF (JJ+NBLOCK.LT.NUMCOL) THEN
               NLOW = JJ + 1
               NUP = JJ + NBLOCK
               INFO(2) = INFO(2) + NBLOCK
            ELSE
               J1 = JFOUND*NBLOCK
               CALL EB13JD(X,LN,N,NUMCOL,J1,U)
               CALL EB13JD(X,LN,N,NUMCOL,J1,U)
               CALL EB13ED(N,X(1,J1+1),LN,W(1,J1+1),LN,NUMCOL-J1)
               NLOW = J1 + 1
               NUP = NUMCOL
               INFO(2) = INFO(2) + NUMCOL - J1
               IPOS = 5
            END IF
         END IF
         IKEEP(11) = IPNTR
         GO TO 290
C
      ELSE IF (IPOS.EQ.5) THEN
C
         J1 = JFOUND*NBLOCK
C
C U() represents an N-BY-NUMCOL matrix which holds the product
C A*X (i.e. A*Q2+ ) we now wish to form the product
C
C           B = QTRANSPOSE * U
C
C        I.E. B = QTRANSPOSE * ( A * Q )
C
         CALL DGEMM('T','N',NUMCOL,NUMCOL-J1,N,ONE,X,LN,U(1,J1+1),LN,
     +              ZERO,QTAQ(1,J1+1),NUMCOL)
C
C Convert B into upper Hessenberg form and
C perform a Ritz acceleration step.
C
         CALL DCOPY(NUMCOL,ER,1,EROLD,1)
         CALL DCOPY(NUMCOL,EI,1,EIOLD,1)
         CALL EB13FD(2,IND,N,NUMCOL,J1,QTAQ,WORK,U,X,LN,KEEP(7),ER,EI,
     +               ITYPE)
C
C Schur vectors held in X.
C Copy into W and return to user to form
C U = A*X (required for convergence test)
C
         NEV = NUMEIG
         IF (EI(NUMEIG).GT.ZERO) NEV = NUMEIG + 1
         J2 = NEV
         IF ((NEV/NBLOCK)*NBLOCK.NE.NEV) J2 = (NEV/NBLOCK+1)*NBLOCK
         CALL EB13ED(N,X(1,J1+1),LN,W(1,J1+1),LN,J2-J1)
         NLOW = J1 + 1
         NUP = J2
         INFO(2) = INFO(2) + J2 - J1
         IPOS = 7
         GO TO 290
C
      ELSE IF (IPOS.EQ.6) THEN
         J1 = JFOUND*NBLOCK
         CALL DCOPY(NUMCOL,ER,1,EROLD,1)
         CALL DCOPY(NUMCOL,EI,1,EIOLD,1)
C Store the matrix computed by the user in the last step.
C Copy it into X.
C
         IF (NBLOCK.GT.1) THEN
            CALL DGEMM('T','N',NUMCOL,NBLOCK,N,ONE,X,LN,
     +                 U(1,NUMCOL-NBLOCK+1),LN,ZERO,
     +                 QTAQ(1,NUMCOL-NBLOCK+1),NUMCOL)
            CALL EB13FD(2,IND,N,NUMCOL,J1,QTAQ,WORK,U,X,LN,KEEP(7),ER,
     +                  EI,ITYPE)
         ELSE
            CALL DGEMV('T',N,NUMCOL,ONE,X,LN,U(1,NUMCOL),1,ZERO,
     +                 QTAQ(1,NUMCOL),1)
C QTAQ is already in Hessenberg form
            CALL EB13FD(1,IND,N,NUMCOL,J1,QTAQ,WORK,U,X,LN,KEEP(7),ER,
     +                  EI,ITYPE)
         END IF
C
C Schur vectors held in X.
C Copy into W and return to user to form
C U = A*X (required for convergence test)
C
         NEV = NUMEIG
         IF (EI(NUMEIG).GT.ZERO) NEV = NUMEIG + 1
         J2 = NEV
         IF ((NEV/NBLOCK)*NBLOCK.NE.NEV) J2 = (NEV/NBLOCK+1)*NBLOCK
         CALL EB13ED(N,X(1,J1+1),LN,W(1,J1+1),LN,J2-J1)
         NLOW = J1 + 1
         NUP = J2
         INFO(2) = INFO(2) + J2 - J1
         IPOS = 7
         GO TO 290
C
      ELSE IF (IPOS.EQ.7) THEN
         J1 = JFOUND*NBLOCK
         NEV = NUMEIG
         IF (EI(NUMEIG).GT.ZERO) NEV = NUMEIG + 1
         J2 = NEV
         IF ((NEV/NBLOCK)*NBLOCK.NE.NEV) J2 = (NEV/NBLOCK+1)*NBLOCK
C We need to keep U = A*X for the next iteration.
         CALL EB13ED(N,U(1,J1+1),LN,W(1,J1+1),LN,J2-J1)
C Test for convergence.
         CALL EB13DD(N,NUMCOL,NBLOCK,NEV,ITYPE,QTAQ,W,X,LN,ER,EI,EROLD,
     +               EIOLD,ICNTL,INFO,RINFO,KEEP(2),EV,CNTL)
         IF (INFO(4).GE.NUMEIG) THEN
            IPOS = 0
C Set the lower traingluar part of QTAQ to zero before return
C (this is nonzero because EB12Z/ZD stores information about
C the reduction to Hessenberg form in the lower triangluar
C part)
            DO 240 I = 3,NUMCOL
               IF (EI(I).GE.ZERO) THEN
                  K = I - 1
               ELSE
                  K = I - 2
               END IF
               DO 230 J = 1,K
                  QTAQ(I,J) = ZERO
  230          CONTINUE
  240       CONTINUE
            GO TO 290
         ELSE
C Check number of matrix-vector products has not been exceeded.
            IF (INFO(2).GT.ICNTL(5)*NUMEIG) THEN
               INFO(1) = -7
               GO TO 290
            END IF
            IPOS = 8
            GO TO 10
         END IF
C
      ELSE IF (IPOS.EQ.8) THEN
C
C Prepare for next iteration.
         IF (ICNTL9.NE.1) THEN
            IF (NBLOCK.GE.NUMEIG) THEN
               KSTOP = NBLOCK + 1
            ELSE
               KSTOP = NUMEIG + 1
            END IF
            XR = ER(KSTOP-1)
            XI = EI(KSTOP-1)
            IF (XI.LT.ZERO) XI = -XI
C Get new convex hull and ellipse
            IEV = 0
            MHULL = IKEEP(10)
C Include points from the previous convex hull
            DO 250 K = 1,MHULL
               IF (IEV.EQ.2*NUMCOL) GO TO 260
               IF (IND.EQ.2) THEN
                  IF (CH(2,K).LT.XI) THEN
                     IEV = IEV + 1
                     EV(1,IEV) = CH(1,K)
                     EV(2,IEV) = CH(2,K)
                  END IF
               ELSE IF (CH(1,K).LT.XR) THEN
                  IEV = IEV + 1
                  EV(1,IEV) = CH(1,K)
                  EV(2,IEV) = CH(2,K)
               END IF
  250       CONTINUE
  260       CONTINUE
C Include new eigenvalue estimates
            DO 270 K = NUMCOL,KSTOP,-1
               IF (EI(K).LT.ZERO) GO TO 270
               IF (IEV.EQ.2*NUMCOL) GO TO 280
               IF (IND.EQ.2) THEN
                  IF (EI(K).LT.XI) THEN
                     IEV = IEV + 1
                     EV(1,IEV) = ER(K)
                     EV(2,IEV) = EI(K)
                  END IF
               ELSE IF (ER(K).LT.XR) THEN
                  IEV = IEV + 1
                  EV(1,IEV) = ER(K)
                  EV(2,IEV) = EI(K)
               END IF
  270       CONTINUE
  280       CONTINUE
            IF (IEV.EQ.0) THEN
C Insufficient points on which to find the complex hull
C Return to the user with an error message
               INFO(1) = -8
               GO TO 290
            END IF
            CALL EB14AD(IEV,MHULL,EV,CH)
            IKEEP(10) = MHULL
C Compute new ellipse
            D = RINFO(1)
            C2 = RINFO(2)
            A2 = RINFO(3)
            IF (ICNTL(8).EQ.2) THEN
C Braconnier ellipse
               CALL EB14CD(MHULL,XR,CH,D,C2,A2,A1)
            ELSE IF (ICNTL(8).EQ.3) THEN
C Saad ellipse
               CALL EB14DD(MHULL,NUMEIG,EV,CH,ER,EI,W,D,C2,A2)
            ELSE
C Ho ellipse
               CALL EB14BD(MHULL,XR,XI,CH,D,C2,A2,A1)
            END IF
            RINFO(1) = D
            RINFO(2) = C2
            RINFO(3) = A2
            CALL EB13GD(NUMEIG,ER,EI,W,RATIO,WRT1,D,C2,A2)
C Parameters for recurrence relation on next iteration
            KEEP(11) = ONE
            IF (C2.LT.ZERO) KEEP(11) = -ONE
            KEEP(12) = SQRT(ABS(C2))/TWO
            IF (KEEP(12).LT.SQRT(KEEP(7))) KEEP(12) = SQRT(KEEP(7))
            IF (ICNTL(8).EQ.2 .OR. ICNTL(8).EQ.3) THEN
               KEEP(8) = TWO*WRT1
            ELSE
               KEEP(8) = A1/KEEP(12)
            END IF
C
C Determine parameters for next iteration step
C
            IF (RATIO.LT.ONE+TOL) RATIO = ONE + TOL
            JPOWER = INT(HALF* (ONE+LOG10(1/KEEP(7))/LOG10(RATIO)))
            JPOWER = MAX(JPOWER,1)
            INFO(3) = INT(DBLE(INFO(3))* (ONE+LOG10(DBLE(INFO(5)))))
            IF (INFO(5).EQ.1) INFO(3) = 40
            INFO(3) = MIN(JPOWER,INFO(3))
            INFO(3) = MIN(INFO(3),ICNTL(3))
            INFO(3) = MAX(INFO(3),1)
            IF (ICNTL(4).GT.0) INFO(3) = ICNTL(4)
            INFO(6) = MAX(INFO(6),INFO(3))
         END IF
         IKEEP(11) = 0
         IPOS = 4
         IF (ICNTL9.EQ.1 .OR. ICNTL9.EQ.2) IPOS = 3
         INFO(5) = INFO(5) + 1
         IF (INFO(5).GT.ICNTL(11)*NUMEIG) THEN
            INFO(1) = -9
            GO TO 290
         END IF
         GO TO 10
C
      END IF
  290 RETURN
C
 9000 FORMAT (/'  Warning from EB13A/AD. INFO(1) = 3',/'  Arnoldi ',
     +       'iteration stopped in fewer than NSTEPS.',/'  Some ',
     +       'spurious zero (or almost zero) eigenvalues may be ',
     +       'returned.')
      END
C*****************************************************************
      SUBROUTINE EB13DD(N,NUMCOL,NBLOCK,NUMEIG,ITYPE,T,AQ,Q,LN,ER,EI,
     +                  EROLD,EIOLD,ICNTL,INFO,RINFO,R,WORK,CNTL)
C
C This subroutine checks for convergence.
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C     N       order of the matrix
C     NUMEIG  number of required eigenvalues
C     ITYPE   array of integer flags set in EB13K/KD. The i-th entry
C             of ITYPE indicates the type of the i-th eigenvalue
C             (refer to comments in EB13K/KD for further details)
C     T       the NUMCOL-by-NUMCOL matrix Q(T)*A*Q which is upper
C             quasi-triangular. the eigenvalues of T are the
C             sought-after eigenvalues of A.
C   * AQ      houses the product A*Q.
C     Q       the LN-by-NUMCOL matrix of computed basis vectors.
C    ER, EI   hold real and imaginary parts of computed eigenvalues
C    EROLD, EIOLD  hold real and imaginary parts of eigenvalues
C             computed on previous iteration
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LN,N,NBLOCK,NUMCOL,NUMEIG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AQ(LN,NUMCOL),EI(NUMCOL),EIOLD(NUMCOL),
     +                 ER(NUMCOL),EROLD(NUMCOL),Q(LN,NUMCOL),R(5),
     +                 RINFO(5),T(NUMCOL,NUMCOL),WORK(NUMCOL),CNTL(1)
      INTEGER ICNTL(11),INFO(6),ITYPE(NUMCOL)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ANORM,EPS,R1,R2,RMACH
      INTEGER I,IFOUND,IPOWER,J,JFOUND,KFOUND,MP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2,FD15AD
      EXTERNAL DNRM2,FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,LOG10,MAX,MIN,SQRT
C     ..
C
      IPOWER = INFO(3)
      EPS = RINFO(4)
      ANORM = RINFO(5)
      MP = ICNTL(2)
      RMACH = SQRT(FD15AD('E'))
C Convergence checking.
C This involves looking at the norms of the columns of AQ - QT.
      IFOUND = INFO(4)
      JFOUND = IFOUND/NBLOCK
      J = JFOUND*NBLOCK
      IF (IFOUND-J.EQ.0) GO TO 30
C Block method.
C Check that the IFOUND-J eigenvalues that had previously converged
C within the block have still converged
      KFOUND = J

   10 CONTINUE
      IF (ITYPE(J+1).GT.0) THEN
C
C Next eigenvalue is complex. Check convergence of  next
C two computed basis vectors, i.e., Q(:,J+1) and Q(:,J+2).
C
         IF (ICNTL(7).EQ.2) THEN
C Stopping criteria is going to use ||(AQ) sub i|| in the
C denominator.
            WORK(J+1) = DNRM2(N,AQ(1,J+1),1)
            WORK(J+2) = DNRM2(N,AQ(1,J+2),1)
         END IF
         CALL DGEMV('N',N,J+2,ONE,Q,LN,T(1,J+1),1,-ONE,AQ(1,J+1),1)
         R1 = DNRM2(N,AQ(1,J+1),1)
         CALL DGEMV('N',N,J+2,ONE,Q,LN,T(1,J+2),1,-ONE,AQ(1,J+2),1)
         R2 = DNRM2(N,AQ(1,J+2),1)
         IF (ICNTL(7).NE.2) THEN
            R(1) = MAX(R1,R2)/ANORM
         ELSE IF (WORK(J+1).LT.RMACH .OR. WORK(J+2).LT.RMACH) THEN
            R(1) = MAX(R1,R2)
         ELSE
            R(1) = MAX(R1/WORK(J+1),R2/WORK(J+2))
         END IF
         IF (R(1).LE.EPS) J = J + 2

      ELSE IF (ITYPE(J+1).EQ.0) THEN

C Next eigenvalue is real so check if Q(:,J+1) has converged.
         IF (ICNTL(7).EQ.2) THEN
C Stopping criteria is going to use ||(AQ) sub i|| in the
C denominator.
            WORK(J+1) = DNRM2(N,AQ(1,J+1),1)
         END IF
         CALL DGEMV('N',N,J+1,ONE,Q,LN,T(1,J+1),1,-ONE,AQ(1,J+1),1)
         R1 = DNRM2(N,AQ(1,J+1),1)
         IF (ICNTL(7).NE.2) THEN
            R(1) = R1/ANORM
         ELSE IF (WORK(J+1).LT.RMACH) THEN
            R(1) = R1
         ELSE
            R(1) = R1/WORK(J+1)
         END IF
         IF (R(1).LE.EPS) J = J + 1
      END IF
C
      IF (J.GT.KFOUND) THEN
         KFOUND = J
         IF (KFOUND.LT.IFOUND) GO TO 10
         GO TO 30
      ELSE
         IFOUND = J
         DO 20 I = 2,5
            R(I) = FD15AD('H')
   20    CONTINUE
         GO TO 60
      END IF
C
   30 CONTINUE
      IF (ITYPE(J+1).GT.0) THEN
C
C Next eigenvalue is complex. Check convergence of  next
C two computed basis vectors, i.e., Q(:,J+1) and Q(:,J+2).
C
         IF (ICNTL(7).EQ.2) THEN
C Stopping criteria is going to use ||(AQ) sub i|| in the
C denominator.
            WORK(J+1) = DNRM2(N,AQ(1,J+1),1)
            WORK(J+2) = DNRM2(N,AQ(1,J+2),1)
         END IF
         CALL DGEMV('N',N,J+2,ONE,Q,LN,T(1,J+1),1,-ONE,AQ(1,J+1),1)
         R1 = DNRM2(N,AQ(1,J+1),1)
         CALL DGEMV('N',N,J+2,ONE,Q,LN,T(1,J+2),1,-ONE,AQ(1,J+2),1)
         R2 = DNRM2(N,AQ(1,J+2),1)
         IF (ICNTL(7).NE.2) THEN
            R(1) = MAX(R1,R2)/ANORM
         ELSE IF (WORK(J+1).LT.RMACH .OR. WORK(J+2).LT.RMACH) THEN
            R(1) = MAX(R1,R2)
         ELSE
            R(1) = MAX(R1/WORK(J+1),R2/WORK(J+2))
         END IF
         R1 = SQRT((ER(J+1)-EROLD(J+1))**2+ (EI(J+1)-EIOLD(J+1))**2)
         R2 = MAX(SQRT(ER(J+1)**2+EI(J+1)**2),
     +        SQRT(EROLD(J+1)**2+EIOLD(J+1)**2))
         IF (R2.GT.ZERO) R1 = R1/R2
         R2 = MAX(R(2),R(3),R(4),R(5))
         IF (R(1).LE.EPS) THEN
            J = J + 2
         ELSE IF (R(1).LT.R(2) .AND. R(2).GT.R(3) .AND.
     +            R(3).LT.R(4) .AND. R(4).GT.R(5)) THEN
            IF (R(1).LT.RMACH .AND. EPS.EQ.CNTL(1)) THEN
               J = J + 2
               EPS = R(1)
               IF (MP.GT.0 .AND. ICNTL(6).GE.2) WRITE (MP,FMT=9000)
               INFO(1) = 2
            END IF
        ELSE IF (R(1).LE.EPS*10**2 .AND. R(1).LT.R(2)) THEN
            IPOWER = MIN(IPOWER,INT(40.0D0* (ONE+LOG10(R(1)/EPS))))
         END IF

      ELSE IF (ITYPE(J+1).EQ.0) THEN

C Next eigenvalue is real so check if Q(:,J+1) has converged.
         IF (ICNTL(7).EQ.2) THEN
C Stopping criteria is going to use ||(AQ) sub i|| in the
C denominator.
            WORK(J+1) = DNRM2(N,AQ(1,J+1),1)
         END IF
         CALL DGEMV('N',N,J+1,ONE,Q,LN,T(1,J+1),1,-ONE,AQ(1,J+1),1)
         R1 = DNRM2(N,AQ(1,J+1),1)
         IF (ICNTL(7).NE.2) THEN
            R(1) = R1/ANORM
         ELSE IF (WORK(J+1).LT.RMACH) THEN
            R(1) = R1
         ELSE
            R(1) = R1/WORK(J+1)
         END IF
         R1 = ABS(ER(J+1)-EROLD(J+1))
         R2 = MAX(ABS(ER(J+1)),ABS(EROLD(J+1)))
         IF (R2.GT.ZERO) R1 = R1/R2
         R2 = MAX(R(2),R(3),R(4),R(5))
         IF (R(1).LE.EPS) THEN
            J = J + 1
         ELSE IF (R(1).LT.R(2) .AND. R(2).GT.R(3) .AND.
     +            R(3).LT.R(4) .AND. R(4).GT.R(5)) THEN
            IF (R(1).LT.RMACH .AND. EPS.EQ.CNTL(1)) THEN
               J = J + 1
               EPS = R(1)
               IF (MP.GT.0 .AND. ICNTL(6).GE.2) WRITE (MP,FMT=9000)
               INFO(1) = 2
            END IF
         ELSE IF (R(1).LE.EPS*10**2 .AND. R(1).LT.R(2)) THEN
            IPOWER = MIN(IPOWER,INT(40.0D0* (ONE+LOG10(R(1)/EPS))))
         END IF

      ELSE

         IF (ICNTL(7).EQ.2) THEN
            WORK(J+1) = DNRM2(N,AQ(1,J+1),1)
         END IF
         CALL DGEMV('N',N,J+1,ONE,Q,LN,T(1,J+1),1,-ONE,AQ(1,J+1),1)
         R1 = DNRM2(N,AQ(1,J+1),1)
         IF (ICNTL(7).NE.2) THEN
            R(1) = R1/ANORM
         ELSE IF (WORK(J+1).LT.RMACH) THEN
            R(1) = R1
         ELSE
            R(1) = R1/WORK(J+1)
         END IF
      END IF
C
C  Check if some new columns converged. if so increase the value
C  of IFOUND and continue the search for converged basis vectors.
C
      IF (J.GT.IFOUND) THEN
         IFOUND = J
         IF (IFOUND.LT.NUMEIG) THEN
            DO 40 I = 2,5
               R(I) = FD15AD('H')
   40       CONTINUE
            GO TO 30
         ELSE
            GO TO 60
         END IF
      ELSE
         DO 50 I = 5,2,-1
            R(I) = R(I-1)
   50    CONTINUE
      END IF
   60 RINFO(4) = EPS
      INFO(4) = IFOUND
C Optional printing at end of iteration.
      IF (MP.GT.0 .AND. ICNTL(6).GE.3) THEN
         IF (IFOUND.GE.NUMEIG) THEN
            WRITE (MP,FMT='(/A/A,6I7/A,1P,5D12.4)')
     +        ' On successful exit: ',' INFO (1:6) = ', (INFO(I),I=1,6),
     +        ' RINFO(1:5) = ', (RINFO(I),I=1,5)
            IF (ICNTL(6).GE.4) THEN
               WRITE (MP,FMT=9010)
               WRITE (MP,FMT=9020) (ER(I),EI(I),I=1,NUMEIG)
            END IF
         ELSE IF (ICNTL(6).GE.5) THEN
            WRITE (MP,FMT='(/A,I3,A/A,6I7/A,1P,5D12.4/A,1PD12.4)')
     +        ' On iteration ',INFO(5),':',' INFO (1:6) = ',
     +        (INFO(I),I=1,6),' RINFO(1:5) = ', (RINFO(I),I=1,5),
     +        ' Residual of first unconverged eigenvalue is ',R(1)
            IF (ICNTL(6).GE.6) THEN
               WRITE (MP,FMT=9010)
               WRITE (MP,FMT=9020) (ER(I),EI(I),I=1,NUMEIG)
            END IF
         END IF
      END IF
      IF (ICNTL(9).NE.1 .AND. ICNTL(9).NE.2) INFO(3) = IPOWER
      RETURN
 9000 FORMAT (/'  Warning from EB13A/AD. INFO(1) = 2',/'  Requested ',
     +       'convergence tolerance not achieved.')
 9010 FORMAT (' Computed eigenvalues:',/'  Real part  Imaginary part  ',
     +       ' Real part  Imaginary part')
 9020 FORMAT (1P,2D12.4,4(' '),1P,2D12.4)
      END
C*****************************************************************
      SUBROUTINE EB13ED(N,W,LDW,V,LDV,M)
C
C Copy the matrix W into the matrix V
C
C     .. Scalar Arguments ..
      INTEGER LDV,LDW,M,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION V(LDV,M),W(LDW,M)
C     ..
C     .. Local Scalars ..
      INTEGER J
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY
C     ..
      DO 10 J = 1,M
         CALL DCOPY(N,W(1,J),1,V(1,J),1)
   10 CONTINUE
      RETURN
      END
C*****************************************************************
      SUBROUTINE EB13FD(JND,IND,N,NUMCOL,IFOUND,QTAQ,Z,QWORK,Q,LN,EPS,
     +                  ER,EI,ITYPE)
C
C     This subroutine performs a Ritz acceleration step
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C     JND        JND = 1 if matrix already in Hessenberg form
C                JND = 2 if matrix must first be reduced to
C                Hessenberg form.
C     IND        IND = 0 eigenvalues of largest modulus are required
C                IND = 1 right-most eigenvalues are required
C                IND = 2 eigenvalues of largest imaginary parts
C                required.
C     N          Order of the matrix A
C   * QTAQ       An NUMCOL-by-NUMCOL array which on entry should hold
C                the product Q(T)*A*Q which is to be reduced to upper
C                Hessenberg form, for subsequent use by EB13K/KD. On
C                return from EB13K/KD it contains T which is a
C                quasi-triangular matrix.
C                The eigenvalues of T are the sought after
C                eigenvalues of A.
C   * Z      An NUMCOL-by-NUMCOL work array which on return from
C                ORTRAN holds the transformation matrix Z, such that
C                         H = Z(T)*B*Z
C                where H is upper Hessenberg (B=Q(T)*A*Q).
C                On return from EB13K/KD it holds the product Z*V
C                where V is the transformation matrix such that
C                         V(T)*H*V = T
C                T being quasi-triangular
C     QWORK      An N-by-NUMCOL array which holds the the matrix
C                Q+ which is from the latest QR factorization.
C   * Q          An LN-by-NUMCOL array which on exit holds the computed
C                Schur vectors. Q(1,1:IFOUND) (IFOUND>=0) houses the
C                converged Schur vectors (deflated system)
C     EPS        A value used by EB13K/KD in a convergence criterion
C   * ER         An NUMCOL-by-1 array which on return holds the real
c                parts of the eigenvalues.
C   * EI         An NUMCOL-by-1 array which on return holds the
C                imaginary parts of the eigenvalues
C   * ITYPE      An NUMCOL-by-1 array which on return holds
C                flags about the type of the i-th eigenvalue.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS
      INTEGER IFOUND,IND,JND,LN,N,NUMCOL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EI(NUMCOL),ER(NUMCOL),Q(LN,NUMCOL),
     +                 QTAQ(NUMCOL,NUMCOL),QWORK(N,NUMCOL),
     +                 Z(NUMCOL,NUMCOL)
      INTEGER ITYPE(NUMCOL)
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,EB13ED,EB13KD,EB13LD
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
      IF (JND.EQ.1) THEN
C  Set up the identity matrix in Z
         DO 20 I = 1,NUMCOL
            DO 10 J = 1,NUMCOL
               Z(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1,NUMCOL
            Z(I,I) = ONE
   30    CONTINUE
      ELSE
C Reduce QTAQ to upper Hessenberg form, then explicitly form the
C transformation matrix. Note that EB13L/LD is based on the
C EISPACK routines ORTHES and ORTRAN.
C
         CALL EB13LD(NUMCOL,IFOUND+1,QTAQ,QWORK(1,IFOUND+1),Z)
C
C Now have H = Z(T)*B*Z . H is stored in QTAQ with the
C transformation matrix stored in Z()
C
      END IF
C
      CALL EB13KD(IND,QTAQ,Z,NUMCOL,IFOUND+1,NUMCOL,EPS,ER,EI,ITYPE,
     +            NUMCOL,NUMCOL)
C
C We now have T = V(T)*H*V. We now wish to form
C                   Q = Q+ * (Z*V)
C
      CALL EB13ED(N,Q(1,IFOUND+1),LN,QWORK(1,IFOUND+1),N,NUMCOL-IFOUND)

C Q+    is stored in QWORK()  (N-by-NUMCOL)
C Z*V    is stored in Z()  (NUMCOL-by-NUMCOL)
C Use BLAS-3 routine DGEMM to form matrix product.

      CALL DGEMM('N','N',N,NUMCOL-IFOUND,NUMCOL,ONE,QWORK,N,
     +           Z(1,IFOUND+1),NUMCOL,ZERO,Q(1,IFOUND+1),LN)
      RETURN
      END
C*******************************************************************
      SUBROUTINE EB13GD(NEV,ER,EI,VEC,RATIO,WRT1,D,C2,A2)
C
C Returns convergence ratios for the eigenvalues ER,EI, relative to
C the ellipse defined by D,C2.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO,HALF
      PARAMETER (ONE=1.0D0,ZERO=0.0D0,HALF=0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A2,C2,D,RATIO,WRT1
      INTEGER NEV
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EI(NEV),ER(NEV),VEC(NEV)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,EE,SMALL,T,X,Y
      INTEGER I
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL KB06AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      SMALL = SQRT(FD15AD('E'))
      EE = SQRT(ABS(C2))
      IF (EE.LT.SMALL) EE = SMALL
      DO 10 I = 1,NEV
         X = (ER(I)-D)/EE
         Y = EI(I)/EE
         IF (C2.LT.ZERO) THEN
            T = X
            X = Y
            Y = T
         END IF
C Calculate major semi-axis of ellipse centre D, focus SQRT(C2)
C passing through (ER(I), EI(I))
         A = HALF* (SQRT((X+ONE)**2+Y**2)+SQRT((X-ONE)**2+Y**2))
         VEC(I) = (A+SQRT(ABS(A**2-ONE)))
   10 CONTINUE
      T = A2/ (EE*EE)
      IF (C2.LT.ZERO) T = ONE + T
      T = SQRT(T) + SQRT(ABS(T-ONE))
C Order ratios (descending order).
      CALL KB06AD(VEC,NEV)
      RATIO = VEC(1)/VEC(NEV)
      WRT1 = (VEC(NEV)+ONE/VEC(NEV))*HALF
      RETURN
      END
C*****************************************************************
      SUBROUTINE EB13HD(N,T1,T2,T3,X,LX,Y,LY,Z,LZ,M)
C
C 3-term Chebychev recurrence relation
C
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION T1,T2,T3
      INTEGER LX,LY,LZ,M,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(LX,M),Y(LY,M),Z(LZ,M)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EMAX,EMIN,T
      INTEGER I,J
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN
C     ..
      T = MAX(T1,T2)
      T = MAX(T,T3)
      T = MAX(T,1.0D0)
C      EMIN = FD05AD(3)* (3.0D0*T)
      EMAX = FD15AD('H')/ (3.0D0*T)
      DO 20 J = 1,M
         DO 10 I = 1,N
C            IF (ABS(X(I,J)).LT.EMIN .AND. X(I,J).NE.ZERO) X(I,
C     +          J) = SIGN(EMIN,X(I,J))
C            IF (ABS(Y(I,J)).LT.EMIN .AND. Y(I,J).NE.ZERO) Y(I,
C     +          J) = SIGN(EMIN,Y(I,J))
C            IF (ABS(Z(I,J)).LT.EMIN .AND. Z(I,J).NE.ZERO) Z(I,
C     +          J) = SIGN(EMIN,Z(I,J))
            IF (ABS(X(I,J)).GT.EMAX) X(I,J) = SIGN(EMAX,X(I,J))
            IF (ABS(Y(I,J)).GT.EMAX) Y(I,J) = SIGN(EMAX,Y(I,J))
            IF (ABS(Z(I,J)).GT.EMAX) Z(I,J) = SIGN(EMAX,Z(I,J))
   10    CONTINUE
   20 CONTINUE
      DO 40 J = 1,M
         DO 30 I = 1,N
            X(I,J) = T1*X(I,J) - T2*Y(I,J) - T3*Z(I,J)
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
C******************************************************************
      SUBROUTINE EB13JD(A,ADIM,N,P,R,WORK)
C
C  Given a full column rank matrix A (NxP), this subroutine computes
C  a matrix Q (NxP) such that the columns of Q are orthonormal and
C  satisfy range(Q(:,1:K)) = range(A(:,1:K)), K = 1:N. In particular,
C  the columns of Q form an orthonormal basis for the column space oF A.
C
C  The subroutine is a level-2 implemention of the modified Gram-Schmidt
C  method. See Golub and van Loan (1989),  "Matrix Computations, 2nd ed.
C  Johns Hopkins University Press, Baltimore, pp. 218ff.
C
C  * A(ADIM,P) Real (DP) array. On entry holds A.
C              On exit holds Q.
C    ADIM      Integer variable. First dimension of the array A.
C    N         Integer variable. Row dimension of the matrix A.
C              (N.LE.ADIM)
C    P         Integer variable. Column dimension of the matrix A.
C    R         Integer variable. It is assumed that A(:,1:R) has
C              orthonormal columns. Set R = 0  if all of A needs to
C              be orthonormalised.
C  * WORK(P)   Real (DP) array.  Workspace.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER ADIM,N,P,R
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(ADIM,P),WORK(P)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RHO
      INTEGER J,K
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMV,DGER,DSCAL
C     ..
      IF (R.GE.1 .AND. R.LT.P) THEN
C  Overwrite A(2) with (I - A(1)A(1)(T))A(2) where A(1) = A(:,1:R)
C  and A(2) = A(:,R+1:N)
C
         DO 10 J = R + 1,P
            CALL DGEMV('T',N,R,ONE,A(1,1),ADIM,A(1,J),1,ZERO,WORK,1)
            CALL DGEMV('N',N,R,-ONE,A(1,1),ADIM,WORK,1,ONE,A(1,J),1)
   10    CONTINUE
      END IF
      K = R + 1
      RHO = DNRM2(N,A(1,K),1)
C
   20 IF (RHO.NE.ZERO .AND. K.LE.P) THEN
C Normalize A(:,K) to get Q(:,K)
         CALL DSCAL(N,ONE/RHO,A(1,K),1)
C
         IF (K.LT.P) THEN
C Compute the K-th row of R, i.e., R(K,K+1:P)^T = A(:,K+1:P)^T Q(:,K)
            CALL DGEMV('T',N,P-K,ONE,A(1,K+1),ADIM,A(1,K),1,ZERO,
     +                 WORK(K+1),1)
C
C Rank-1 update: A(:,K+1:P) = A(:,K+1:P) - Q(:,K)R(K,K+1:P)
            CALL DGER(N,P-K,-ONE,A(1,K),1,WORK(K+1),1,A(1,K+1),ADIM)
            RHO = DNRM2(N,A(1,K+1),1)
         END IF
         K = K + 1
         GO TO 20
      END IF
C
      RETURN
      END
C****************************************************************
      SUBROUTINE EB13KD(IND,A,V,N,NLOW,NUP,EPS,ER,EI,TYPE,NA,NV)
C
C This is a modified version of HQR3.
C (Stewart ACM Trans Math Software 1976, Vol 2, 275-280).
C The modification only concerns the ordering along the diagonal.
C This version is suitable for use when computing
C the eigenvalues of largest (or smallest) real parts or for computing
C the eigenvalues of largest modulli. EB13K/KD reduces the upper
C Hessenberg matrix A to quasi-triangular form by unitary similarity
C transformations. The eigenvalues of A, which are contained in the 1x1
C and 2x2 diagonal blocks of the reduced matrix, are ordered in order
C of largest (or smallest) real part along the diagonal.
C The transformations are accumulated in the array V.
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C
C     IND     IND = 0 eigenvalues of largest modulus are required
C             IND = 1 right-most eigenvalues are required
C             IND = 2 eigenvalues of largest imaginary parts required.
C    *A       An array that initially contains the N X N
C             upper hessenberg matrix to be reduced.  On
C             return a contains the reduced, quasi-triangular matrix.
C    *V       An array that contains a matrix into which the reducing
C             transformations are to be multiplied.
C     N       The order of the matrices A and V.
C     NLOW    A(NLOW,NLOW-1) and A(NUP,+1,NUP) are assumed to be
C     NUP     zero, and only rows NLOW through NUP and columns
C             NLOW through NUP are transformed, resulting in the
C             calculation of eigenvalues NLOW through NUP.
C     EPS     A convergence criterion.
C    *ER      An array that on return contains the real parts
C             of the eigenvalues.
C    *EI      An array that on return contains the imaginary parts
C             of the eigenvalues.
C    *TYPE    An integer array whose I-th entry is
C             0   if the I-th eigenvalue is real,
C             1   if the I-th eigenvalue is complex
C                 with positive imaginary part.
C             2   if the I-th eigenvalue is complex
C                 with negative imaginary part.
C            -1   if the I-th eigenvalue not calculated successfully.
C     NA      The first dimension of the array A.
C     NV      The first dimension of the array V.
C
C The convergence criterion EPS is used to determine when a subdiagonal
C element of A is negligible. Specifically, A(I+1,I) is regarded as
C negligible if
C
C        ABS(A(I+1),I)) .LE. EPS*(ABS(A(I,I))+ABS(A(I+1,I+1))).
C
C This means that the final matrix returned by the program will
C be exactly similar to A + E where E is of order EPS*NORM(A),
C for any reasonably balanced norm such as the row-sum norm.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS
      INTEGER IND,N,NA,NLOW,NUP,NV
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NA,N),EI(N),ER(N),V(NV,N)
      INTEGER TYPE(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION E1,E2,P,Q,R,RMACH,S,T,W,X,Y,Z
      INTEGER I,IT,L,MU,NL,NU
      LOGICAL FAIL
C     ..
C     .. External Subroutines ..
      EXTERNAL EB13MD,EB13ND,EB13OD
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      RMACH = FD15AD('E')
C Initialise TYPE
      DO 10 I = NLOW,NUP
         TYPE(I) = -1
   10 CONTINUE
      T = ZERO
C Main loop. Find and order eigenvalues.
      NU = NUP
   20 CONTINUE
      IF (NU.LT.NLOW) GO TO 160
      IT = 0
   30 CONTINUE
C QR loop.  Find negligible elements and perform QR steps
C Search back for negligible elements.
      L = NU
   40 CONTINUE
C       IF (L.NE.NLOW .AND. ABS(A(L,L-1)).GT.
C     +    EPS* (ABS(A(L-1,L-1))+ABS(A(L,L)))) THEN
      IF (L.NE.NLOW) THEN
        IF (ABS(A(L,L-1)).GT.EPS* (ABS(A(L-1,L-1))+ABS(A(L,L)))) THEN
          L = L - 1
          GO TO 40
        END IF
      END IF
C Test to see if an eigenvalue or a 2x2 block has been found.
      X = A(NU,NU)
      IF (L.EQ.NU) GO TO 80
      Y = A(NU-1,NU-1)
      W = A(NU,NU-1)*A(NU-1,NU)
      IF (L.NE.NU-1) THEN
C Test iteration count. If IT is 30 quit.  If
C IT is 10 or 20 set up an ad-hoc shift.
         IF (IT.EQ.30) GO TO 160
         IF (IT.EQ.10 .OR. IT.EQ.20) THEN
C ad-hoc shift.
            T = T + X
            DO 50 I = NLOW,NU
               A(I,I) = A(I,I) - X
   50       CONTINUE
            S = ABS(A(NU,NU-1)) + ABS(A(NU-1,NU-2))
            X = 0.75D0*S
            Y = X
            W = -0.4375D0*S**2
         END IF
         IT = IT + 1
C Look for 2 consecutive small sub-diagonal elements
         NL = NU - 2
   60    CONTINUE
         Z = A(NL,NL)
         R = X - Z
         S = Y - Z
         P = (R*S-W)/A(NL+1,NL) + A(NL,NL+1)
         Q = A(NL+1,NL+1) - Z - R - S
         R = A(NL+2,NL+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
C         IF (NL.NE.L .AND. ABS(A(NL,NL-1))* (ABS(Q)+ABS(R)).GT.
C     +       EPS*ABS(P)* (ABS(A(NL-1,NL-1))+ABS(Z)+ABS(A(NL+1,
C     +       NL+1)))) THEN
         IF (NL.NE.L) THEN
           IF (ABS(A(NL,NL-1))* (ABS(Q)+ABS(R)).GT.
     +         EPS*ABS(P)* (ABS(A(NL-1,NL-1))+ABS(Z)+ABS(A(NL+1,
     +         NL+1)))) THEN
             NL = NL - 1
             GO TO 60
           END IF
         END IF
C Perform a QR step between NL and NU.
         CALL EB13ND(A,V,P,Q,R,NL,NU,N,NA,NV)
         GO TO 30
      END IF
C 2x2 block found.
      IF (NU.NE.NLOW+1) A(NU-1,NU-2) = ZERO
      A(NU,NU) = A(NU,NU) + T
      A(NU-1,NU-1) = A(NU-1,NU-1) + T
      TYPE(NU) = 0
      TYPE(NU-1) = 0
      MU = NU
   70 CONTINUE
C Loop to position  2x2 block.
      NL = MU - 1
C Attempt  to EB13OD the block into two real eigenvalues
      CALL EB13OD(A,V,N,NL,E1,E2,NA,NV,RMACH)
C If the EB13OD was successful, go and order the
C real eigenvalues.
      IF (A(MU,MU-1).EQ.ZERO) GO TO 90
C Test to see if the block is properly positioned,
C and if not exchange it
      IF (MU.EQ.NUP) GO TO 120
      IF (MU.NE.NUP-1) THEN
         IF (A(MU+2,MU+1).NE.ZERO) THEN
C The next block is 2X2.
            IF (IND.EQ.2) THEN
               IF (ABS((A(MU-1,MU-1)-A(MU,MU))**2+4.0D0*A(MU,
     +             MU-1)*A(MU-1,MU)).GE.ABS((A(MU+1,MU+1)-A(MU+2,
     +             MU+2))**2+4.0D0*A(MU+1,MU+2)*A(MU+2,
     +             MU+1))) GO TO 120
            ELSE IF (IND.EQ.0) THEN
               IF (A(MU-1,MU-1)*A(MU,MU)-A(MU-1,MU)*A(MU,MU-1).GE.
     +             A(MU+1,MU+1)*A(MU+2,MU+2)-
     +             A(MU+1,MU+2)*A(MU+2,MU+1)) GO TO 120
            ELSE IF (A(MU-1,MU-1)+A(MU,MU).GE.
     +               A(MU+1,MU+1)+A(MU+2,MU+2)) THEN
               GO TO 120
            END IF
            CALL EB13MD(A,V,N,NL,2,2,EPS,FAIL,NA,NV,RMACH)
            IF (FAIL) GO TO 130
            MU = MU + 2
            GO TO 70
         END IF
      END IF
C The next block is 1X1.
      IF (IND.EQ.2) GO TO 120
      IF (IND.EQ.0) THEN
         IF (A(MU-1,MU-1)*A(MU,MU)-A(MU-1,MU)*A(MU,MU-1).GE.
     +       A(MU+1,MU+1)**2) GO TO 120
      ELSE IF ((A(MU-1,MU-1)+A(MU,MU))/2.0D0.GE.A(MU+1,MU+1)) THEN
         GO TO 120
      END IF
      CALL EB13MD(A,V,N,NL,2,1,EPS,FAIL,NA,NV,RMACH)
      IF (FAIL) GO TO 140
      MU = MU + 1
      GO TO 70
C Single eigenvalue found.
   80 NL = 0
      A(NU,NU) = A(NU,NU) + T
      IF (NU.NE.NLOW) A(NU,NU-1) = ZERO
      TYPE(NU) = 0
      MU = NU
   90 CONTINUE
C Loop to position one or two real eigenvalues.
  100 CONTINUE
C Position the eigenvalue located at A(NL,NL).
      IF (MU.NE.NUP) THEN
         IF (MU.NE.NUP-1) THEN
            IF (A(MU+2,MU+1).NE.ZERO) THEN
C The next block is 2X2.
               IF (IND.EQ.0) THEN
                  IF (A(MU,MU)**2.GE.A(MU+1,MU+1)*A(MU+2,MU+2)-
     +                A(MU+1,MU+2)*A(MU+2,MU+1)) GO TO 110
               ELSE IF (IND.NE.2) THEN
                  IF (A(MU,MU).GE. (A(MU+1,MU+1)+A(MU+2,MU+2))/
     +                2.0D0) GO TO 110
               END IF
               CALL EB13MD(A,V,N,MU,1,2,EPS,FAIL,NA,NV,RMACH)
               IF (FAIL) GO TO 150
               MU = MU + 2
               GO TO 100
            END IF
         END IF
C The next block is 1X1.
         IF (IND.EQ.2) GO TO 110
         IF (IND.EQ.0) THEN
            IF (ABS(A(MU,MU)).GE.ABS(A(MU+1,MU+1))) GO TO 110
         ELSE IF (A(MU,MU).GE.A(MU+1,MU+1)) THEN
            GO TO 110
         END IF
         CALL EB13MD(A,V,N,MU,1,1,EPS,FAIL,NA,NV,RMACH)
         MU = MU + 1
         GO TO 100
      END IF
  110 MU = NL
      NL = 0
      IF (MU.NE.0) GO TO 90
C Go back and get the next eigenvalue.
  120 CONTINUE
      NU = L - 1
      GO TO 20
  130 TYPE(NL) = -1
      TYPE(NL+1) = -1
      TYPE(NL+2) = -1
      TYPE(NL+3) = -1
      GO TO 160
  140 TYPE(NL) = -1
      TYPE(NL+1) = -1
      TYPE(NL+2) = -1
      GO TO 160
  150 TYPE(MU) = -1
      TYPE(MU+1) = -1
      TYPE(MU+2) = -1
C All the eigenvalues have been found and ordered.
C compute their values and type.
  160 DO 170 I = NLOW,NU
         A(I,I) = A(I,I) + T
  170 CONTINUE
      NU = NUP
  180 CONTINUE
      IF (TYPE(NU).NE.-1) THEN
         IF (NU.NE.NLOW) THEN
            IF (A(NU,NU-1).NE.ZERO) THEN
C 2X2 block.
               CALL EB13OD(A,V,N,NU-1,E1,E2,NA,NV,RMACH)
               IF (A(NU,NU-1).NE.ZERO) THEN
                  ER(NU) = E1
                  EI(NU-1) = E2
                  ER(NU-1) = ER(NU)
                  EI(NU) = -EI(NU-1)
                  TYPE(NU-1) = 1
                  TYPE(NU) = 2
                  NU = NU - 2
                  GO TO 190
               END IF
            END IF
         END IF
C Single root.
         ER(NU) = A(NU,NU)
         EI(NU) = ZERO
      END IF
      NU = NU - 1
  190 IF (NU.GE.NLOW) GO TO 180
      RETURN
      END
C*********************************************************************
      SUBROUTINE EB13LD(N,LOW,A,ORT,Z)
C
C     This subroutine based on the EISPACK routines ORTHES and ORTRAN.
C
C     Given a real general matrix A, this subroutine reduces a submatrix
C     situated in rows and columns LOW to N to upper Hessenberg
C     form by orthogonal similarity transformations and
C     accumulates the orthogonal similarity transformations
C     in Z.
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C
C     N      order of the matrix
C     LOW    rows and columns LOW to N are to be reduced to
C            upper Hessenberg form.
C     *A     on input contains the matrix to be reduced. On exit,
C            holds the upper Hessenberg form.
C     *ORT   not set on entry. Used as workspace.
C            Only entries LOW to N are used.
C     *Z     Not set on entry. On exit,
C            contains the transformation matrix.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LOW,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(N,N),ORT(N),Z(N,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,G,H,SCALE
      INTEGER I,II,J,JJ,M,MM,MP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,SQRT
C     ..
C
C Initialize Z to identity matrix
      DO 20 J = 1,N
         DO 10 I = 1,N
            Z(I,J) = ZERO
   10    CONTINUE
         Z(J,J) = ONE
   20 CONTINUE
C
C Loop over columns
      DO 110 M = LOW + 1,N - 1
         H = ZERO
         ORT(M) = ZERO
         SCALE = ZERO
C Scale column
         DO 30 I = M,N
            SCALE = SCALE + ABS(A(I,M-1))
   30    CONTINUE
C Jump if zero column.
         IF (SCALE.EQ.ZERO) GO TO 110
         MP = M + N
         DO 40 II = M,N
            I = MP - II
            ORT(I) = A(I,M-1)/SCALE
            H = H + ORT(I)*ORT(I)
   40    CONTINUE
         G = -SIGN(SQRT(H),ORT(M))
         H = H - ORT(M)*G
         ORT(M) = ORT(M) - G
C Form (I-(U*UT)/H) * A
         DO 70 J = M,N
            F = ZERO
            DO 50 II = M,N
               I = MP - II
               F = F + ORT(I)*A(I,J)
   50       CONTINUE
            F = F/H
            DO 60 I = M,N
               A(I,J) = A(I,J) - F*ORT(I)
   60       CONTINUE
   70    CONTINUE
C Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 100 I = 1,N
            F = ZERO
            DO 80 JJ = M,N
               J = MP - JJ
               F = F + ORT(J)*A(I,J)
   80       CONTINUE
            F = F/H
            DO 90 J = M,N
               A(I,J) = A(I,J) - F*ORT(J)
   90       CONTINUE
  100    CONTINUE
         ORT(M) = SCALE*ORT(M)
         A(M,M-1) = SCALE*G
  110 CONTINUE
C
C Accumulate transformations.
      DO 160 MM = 1,N - LOW - 1
         MP = N - MM
         IF (A(MP,MP-1).EQ.ZERO) GO TO 160
         DO 120 I = MP + 1,N
            ORT(I) = A(I,MP-1)
  120    CONTINUE
         DO 150 J = MP,N
            G = ZERO
            DO 130 I = MP,N
               G = G + ORT(I)*Z(I,J)
  130       CONTINUE
            G = (G/ORT(MP))/A(MP,MP-1)
            DO 140 I = MP,N
               Z(I,J) = Z(I,J) + G*ORT(I)
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
C
      RETURN
      END
C******************************************************************
      SUBROUTINE EB13MD(A,V,N,L,B1,B2,EPS,FAIL,NA,NV,RMACH)
C
C Given the upper Hessenberg matrix A with consecutive B1xB1 and
C B2xB2 diagonal blocks (B1,B2 .le. 2) starting at A(L,L),
C EB13M/MD produces a unitary similarity transformation that
C exchanges the blocks along with their eigenvalues.
C The transformation is accumulated in V.
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C
C    *A       The matrix whose blocks are to be interchanged.
C    *V       The array into which the transformations
C             are to be accumulated.
C     N       The order of the matrix A.
C     L       The position of the blocks.
C     B1      An integer containing the size of the first block.
C     B2      An integer containing the size of the second block.
C     EPS     A convergence criterion (CF. EB13K/KD).
C    *FAIL    A logical variable which is false on a normal
C             return.  If 30 iterations were performed without
C             convergence, FAIL is set TRUE and the element
C             A(L+B2,L+B2-1) cannot be assumed zero.
C     NA      the first dimension of the array A.
C     NV      the first dimension of the array V.
C     RMACH   machine precision.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS,RMACH
      INTEGER B1,B2,L,N,NA,NV
      LOGICAL FAIL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NA,N),V(NV,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION P,Q,R,S,W,X,Y,Z
      INTEGER I,IT,J,M
C     ..
C     .. External Subroutines ..
      EXTERNAL EB13ND
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SQRT
C     ..
      FAIL = .FALSE.
      IF (B1.EQ.2) THEN
C Interchange 2X2 and B2XB2 blocks.
         M = L + 2
         IF (B2.EQ.2) M = M + 1
         X = A(L+1,L+1)
         Y = A(L,L)
         W = A(L+1,L)*A(L,L+1)
         P = ONE
         Q = ONE
         R = ONE
         CALL EB13ND(A,V,P,Q,R,L,M,N,NA,NV)
         IT = 0
   10    CONTINUE
         IT = IT + 1
         IF (IT.GT.30) THEN
            FAIL = .TRUE.
         ELSE
            Z = A(L,L)
            R = X - Z
            S = Y - Z
            P = (R*S-W)/A(L+1,L) + A(L,L+1)
            Q = A(L+1,L+1) - Z - R - S
            R = A(L+2,L+1)
            S = ABS(P) + ABS(Q) + ABS(R)
            P = P/S
            Q = Q/S
            R = R/S
            CALL EB13ND(A,V,P,Q,R,L,M,N,NA,NV)
            IF (ABS(A(M-1,M-2)).GT.EPS* (ABS(A(M-1,M-1))+ABS(A(M-2,
     +          M-2)))) GO TO 10
         END IF
         A(M-1,M-2) = ZERO
      ELSE IF (B2.EQ.2) THEN
C Interchange 1X1 and 2X2 blocks.
         X = A(L,L)
         P = ONE
         Q = ONE
         R = ONE
         CALL EB13ND(A,V,P,Q,R,L,L+2,N,NA,NV)
         IT = 0
   20    CONTINUE
         IT = IT + 1
         IF (IT.GT.30) THEN
            FAIL = .TRUE.
         ELSE
            P = A(L,L) - X
            Q = A(L+1,L)
            R = ZERO
            CALL EB13ND(A,V,P,Q,R,L,L+2,N,NA,NV)
            IF (ABS(A(L+2,L+1)).GT.EPS* (ABS(A(L+1,L+1))+ABS(A(L+2,
     +          L+2)))) GO TO 20
         END IF
         A(L+2,L+1) = ZERO
      ELSE
C Interchange 1X1 and 1X1 blocks.
         Q = A(L+1,L+1) - A(L,L)
         P = A(L,L+1)
         R = MAX(ABS(P),ABS(Q))
         IF (R.GE.RMACH) THEN
            P = P/R
            Q = Q/R
            R = SQRT(P**2+Q**2)
            P = P/R
            Q = Q/R
            DO 30 J = L,N
               S = P*A(L,J) + Q*A(L+1,J)
               A(L+1,J) = P*A(L+1,J) - Q*A(L,J)
               A(L,J) = S
   30       CONTINUE
            DO 40 I = 1,L + 1
               S = P*A(I,L) + Q*A(I,L+1)
               A(I,L+1) = P*A(I,L+1) - Q*A(I,L)
               A(I,L) = S
   40       CONTINUE
            DO 50 I = 1,N
               S = P*V(I,L) + Q*V(I,L+1)
               V(I,L+1) = P*V(I,L+1) - Q*V(I,L)
               V(I,L) = S
   50       CONTINUE
            A(L+1,L) = ZERO
         END IF
      END IF
      RETURN
      END
C******************************************************************
      SUBROUTINE EB13ND(A,V,P,Q,R,NL,NU,N,NA,NV)
C
C EB13N/ND performs one implicit QR step on the upper Hessenberg
C matrix A.  The shift is determined by the numbers P,Q, and R,
C and the step is applied to rows and columns NL through NU.
C The transformations are accumulated in V.
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C
C    *A       The matrix on which the QR step is to be performed.
C    *V       The array into which the transformations
C             are to be accumulated.
C    *P,*Q,*R Parameters that determine the shift.
C     NL,NU   The lower and upper limits of the step.
C     N       The order of the matrix A.
C     NA      the first dimension of the array A.
C     NV      the first dimension of the array V.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION P,Q,R
      INTEGER N,NA,NL,NU,NV
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NA,N),V(NV,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S,X,Y,Z
      INTEGER I,J,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SQRT
C     ..
      DO 10 I = NL + 2,NU
         A(I,I-2) = ZERO
   10 CONTINUE
      DO 20 I = NL + 3,NU
         A(I,I-3) = ZERO
   20 CONTINUE
      DO 60 K = NL,NU - 2
C Determine the transformation.
         IF (K.NE.NL) THEN
            P = A(K,K-1)
            Q = A(K+1,K-1)
            R = A(K+2,K-1)
            X = ABS(P) + ABS(Q) + ABS(R)
            IF (X.EQ.ZERO) GO TO 60
            P = P/X
            Q = Q/X
            R = R/X
            S = SQRT(P**2+Q**2+R**2)
            IF (P.LT.ZERO) S = -S
            A(K,K-1) = -S*X
         ELSE
            S = SQRT(P**2+Q**2+R**2)
            IF (P.LT.ZERO) S = -S
            IF (NL.NE.1) A(K,K-1) = -A(K,K-1)
         END IF
         P = P + S
         X = P/S
         Y = Q/S
         Z = R/S
         Q = Q/P
         R = R/P
C Premultiply.
         DO 30 J = K,N
            P = A(K,J) + Q*A(K+1,J) + R*A(K+2,J)
            A(K+2,J) = A(K+2,J) - P*Z
            A(K+1,J) = A(K+1,J) - P*Y
            A(K,J) = A(K,J) - P*X
   30    CONTINUE
C Postmultiply.
         J = MIN(K+3,NU)
         DO 40 I = 1,J
            P = X*A(I,K) + Y*A(I,K+1) + Z*A(I,K+2)
            A(I,K+2) = A(I,K+2) - P*R
            A(I,K+1) = A(I,K+1) - P*Q
            A(I,K) = A(I,K) - P
   40    CONTINUE
C Accumulate the transformation in V.
         DO 50 I = 1,N
            P = X*V(I,K) + Y*V(I,K+1) + Z*V(I,K+2)
            V(I,K+2) = V(I,K+2) - P*R
            V(I,K+1) = V(I,K+1) - P*Q
            V(I,K) = V(I,K) - P
   50    CONTINUE
   60 CONTINUE
      K = NU - 1
C Determine the transformation.
      IF (K.NE.NL) THEN
         P = A(K,K-1)
         Q = A(K+1,K-1)
         R = ZERO
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X.EQ.ZERO) GO TO 100
         P = P/X
         Q = Q/X
         S = SQRT(P**2+Q**2+R**2)
         IF (P.LT.ZERO) S = -S
         A(K,K-1) = -S*X
      ELSE
         S = SQRT(P**2+Q**2+R**2)
         IF (P.LT.ZERO) S = -S
         IF (NL.NE.1) A(K,K-1) = -A(K,K-1)
      END IF
      P = P + S
      X = P/S
      Y = Q/S
      Q = Q/P
      R = R/P
C Premultiply.
      DO 70 J = K,N
         P = A(K,J) + Q*A(K+1,J)
         A(K+1,J) = A(K+1,J) - P*Y
         A(K,J) = A(K,J) - P*X
   70 CONTINUE
C Postmultiply.
      DO 80 I = 1,NU
         P = X*A(I,K) + Y*A(I,K+1)
         A(I,K+1) = A(I,K+1) - P*Q
         A(I,K) = A(I,K) - P
   80 CONTINUE
C Accumulate the transformation in V.
      DO 90 I = 1,N
         P = X*V(I,K) + Y*V(I,K+1)
         V(I,K+1) = V(I,K+1) - P*Q
         V(I,K) = V(I,K) - P
   90 CONTINUE
  100 CONTINUE
      RETURN
      END
C*******************************************************************
      SUBROUTINE EB13OD(A,V,N,L,E1,E2,NA,NV,RMACH)
C
C Given the upper Hessenberg matrix A with a 2x2 block starting
C at A(L,L), this routine determines if the corresponding eigenvalues
C are real or complex. If they are real, a rotation is determined
C that reduces the block to upper triangular form with the
C eigenvalue  of largest absolute value appearing first.
C The rotation is accumulated in V.  The eigenvalues (real
C or complex) are returned in E1 and E2.
C
C     Arguments to subroutine  (* indicates that it is overwritten
C                                 in the subroutine)
C
C    *A       The upper Hessenberg matrix.
C    *V       The array into which the transformations
C             are to be accumulated.
C     N       The order of the matrix A.
C     L       The position of the 2X2 block.
C    *E1      On return, if the eigenvalues are complex E1 contains
C    *E2      their common real part and E2 contains the positive
C             imaginary part. If the eigenvalues are real, E1
c             contains the one largest in absolute value and E2
c             contains the other one.
C     NA      the first dimension of the array A.
C     NV      the first dimension of the array V.
C     RMACH   machine precision.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION E1,E2,RMACH
      INTEGER L,N,NA,NV
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NA,N),V(NV,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION P,Q,R,T,U,W,X,Y,Z
      INTEGER I,J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      X = A(L+1,L+1)
      Y = A(L,L)
      W = A(L,L+1)*A(L+1,L)
      P = (Y-X)/2.0D0
      Q = P**2 + W
      IF (Q.LT.ZERO) THEN
C Complex eigenvalue.
         E1 = P + X
         E2 = SQRT(-Q)
      ELSE
C Two real eigenvalues.  Set up transformation.
         Z = SQRT(Q)
         IF (P.GE.ZERO) THEN
            Z = P + Z
         ELSE
            Z = P - Z
         END IF
         R = ZERO
         IF (ABS(Z).GE.RMACH) R = -W/Z
         IF (ABS(X+Z).GE.ABS(X+R)) Z = R
         Y = Y - X - Z
         X = -Z
         T = A(L,L+1)
         U = A(L+1,L)
         IF (ABS(Y)+ABS(U).LE.ABS(T)+ABS(X)) THEN
            Q = X
            P = T
         ELSE
            Q = U
            P = Y
         END IF
         R = SQRT(P**2+Q**2)
         IF (R.LE.RMACH) THEN
            E1 = A(L,L)
            E2 = A(L+1,L+1)
            A(L+1,L) = ZERO
         ELSE
            P = P/R
            Q = Q/R
C Premultiply.
            DO 10 J = L,N
               Z = A(L,J)
               A(L,J) = P*Z + Q*A(L+1,J)
               A(L+1,J) = P*A(L+1,J) - Q*Z
   10       CONTINUE
C Postmultiply.
            DO 20 I = 1,L + 1
               Z = A(I,L)
               A(I,L) = P*Z + Q*A(I,L+1)
               A(I,L+1) = P*A(I,L+1) - Q*Z
   20       CONTINUE
C Accumulate the transformation in V.
            DO 30 I = 1,N
               Z = V(I,L)
               V(I,L) = P*Z + Q*V(I,L+1)
               V(I,L+1) = P*V(I,L+1) - Q*Z
   30       CONTINUE
            A(L+1,L) = ZERO
            E1 = A(L,L)
            E2 = A(L+1,L+1)
         END IF
      END IF
      RETURN
      END
C****************************************************************
      SUBROUTINE EB13PD(A,B,C,D,E,F)
C     Complex division
C
C     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,D,E,F
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (ABS(C).LE.ABS(D)) THEN
         Y = C/D
         X = ONE/ (Y*C+D)
         E = (Y*A+B)*X
         F = (Y*B-A)*X
      ELSE
         Y = D/C
         X = ONE/ (Y*D+C)
         E = (A+Y*B)*X
         F = (B-Y*A)*X
      END IF
      RETURN
      END
C**************************************************************
      SUBROUTINE EB13QD(N,NUMEIG,LN,Y,AY,RES,ER,EI,JNFO,ANORM)
C
C      Output residuals ||AY - LAMBDA*Y||/||A||
C      (or ||AY - LAMBDA*Y||/||AY||)  (||Y||=1)
C
C  Argument list (an * indicates the argument is altered by the routine)
C
C     N         Integer variable. Must be set to the order of the
C               matrix
C      Y(LN,NUMEIG)  Columns of Y are computed eigenvectors of A.
C                   Note: the eigenvectors have already been
C                   normalised.
C    * AY(LN,NUMEIG)  On entry, must hold A*Y. The contents of this
C                   array are destroyed by the routine.
C    * RES(NUMEIG) Not set on entry.
C                  On exit, holds the computed residuals.
C      ER(NUMEIG), EI(NUMEIG) Hold the real and imaginary
C                  parts of the computed eigenvalues
C    * JNFO        On exit holds the number of eigenvectors for
C                  which ||(AY) sub i||< u,
C                  u=machine precision... only set if ICNTL(7)= 1
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ANORM
      INTEGER JNFO,LN,N,NUMEIG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AY(LN,NUMEIG),EI(NUMEIG),ER(NUMEIG),RES(NUMEIG),
     +                 Y(LN,NUMEIG)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RMACH,T,T1,TETA,XMU
      INTEGER J,J1,K
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,FD15AD
      EXTERNAL DDOT,DNRM2,FD15AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      J = 1
      IF (ANORM.EQ.ZERO) THEN
C This is the case when we are computing ||AY - LAMBDA*Y||/||AY||
C i.e. ICNTL(7) = 2.
         RMACH = FD15AD('E')
   10    IF (EI(J).EQ.ZERO) THEN
C  Case of a real eigenvalue
            T1 = DNRM2(N,AY(1,J),1)
            XMU = ER(J)
            DO 20 K = 1,N
               AY(K,J) = AY(K,J) - XMU*Y(K,J)
   20       CONTINUE
            T = DNRM2(N,AY(1,J),1)
C  Test to see if norm is small
            IF (T1.LT.RMACH) THEN
               RES(J) = T
               JNFO = JNFO + 1
            ELSE
               RES(J) = T/T1
            END IF
            J = J + 1
         ELSE
C Case of a complex eigenvalue
            J1 = J + 1
            XMU = ER(J)
            TETA = EI(J)
            T1 = DDOT(N,AY(1,J),1,AY(1,J),1) +
     +           DDOT(N,AY(1,J1),1,AY(1,J1),1)
            T1 = SQRT(T1)
            DO 30 K = 1,N
               AY(K,J) = AY(K,J) - XMU*Y(K,J) + TETA*Y(K,J1)
               AY(K,J1) = AY(K,J1) - XMU*Y(K,J1) - TETA*Y(K,J)
   30       CONTINUE
            T = DDOT(N,AY(1,J),1,AY(1,J),1) +
     +          DDOT(N,AY(1,J1),1,AY(1,J1),1)
            T = SQRT(T)
            IF (T1.LT.RMACH) THEN
               RES(J) = T
               JNFO = JNFO + 2
            ELSE
               RES(J) = T/T1
            END IF
            RES(J1) = RES(J)
            J = J + 2
         END IF
         IF (J.LE.NUMEIG) GO TO 10
      ELSE
C We are computing ||AY - LAMBDA*Y||/||A||
   40    IF (EI(J).EQ.ZERO) THEN
C  Case of a real eigenvalue
            XMU = ER(J)
            DO 50 K = 1,N
               AY(K,J) = AY(K,J) - XMU*Y(K,J)
   50       CONTINUE
            RES(J) = DNRM2(N,AY(1,J),1)/ANORM
            J = J + 1
         ELSE
C Case of a complex eigenvalue
            J1 = J + 1
            XMU = ER(J)
            TETA = EI(J)
            DO 60 K = 1,N
               AY(K,J) = AY(K,J) - XMU*Y(K,J) + TETA*Y(K,J1)
               AY(K,J1) = AY(K,J1) - XMU*Y(K,J1) - TETA*Y(K,J)
   60       CONTINUE
            T = DDOT(N,AY(1,J),1,AY(1,J),1) +
     +          DDOT(N,AY(1,J1),1,AY(1,J1),1)
            T = SQRT(T)
            RES(J) = T/ANORM
            RES(J1) = RES(J)
            J = J + 2
         END IF
         IF (J.LE.NUMEIG) GO TO 40
      END IF
      RETURN
      END
C*******************************************************************
      SUBROUTINE EB13RD(N,NUMCOL,NUMEIG,ER,EI,Y,X,W,LN,TW)
C
C     This routine computes the (right) eigenvectors corresponding
C     to the eigenvalues computed by EB13A/AD
C
C     Arguments: see EB13B/BD .
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LN,N,NUMCOL,NUMEIG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EI(NUMCOL),ER(NUMCOL),TW(NUMCOL,NUMCOL),
     +                 W(LN,NUMEIG),X(LN,NUMCOL),Y(LN,NUMEIG)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION E,F,P,Q,R,RA,RMACH,S,SA,T1,U1,VI,VR,WI,WR,WW,X1,
     +                 XN,Y1,Z
      INTEGER I,II,J,M,NA,NE
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,FD15AD
      EXTERNAL DDOT,DNRM2,FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,DSCAL,EB13ED,EB13PD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      RMACH = FD15AD('E')
C  Calculate norm
      XN = ZERO
      DO 20 I = 1,NUMCOL
         DO 10 J = I,NUMCOL
            XN = XN + ABS(TW(I,J))
   10    CONTINUE
   20 CONTINUE
C
C Backsubstitution
      DO 50 NE = NUMCOL,1,-1
         P = ER(NE)
         Q = EI(NE)
         NA = NE - 1
         IF (Q.LT.ZERO) THEN
C
C Complex vector associated with lambda = P - I * Q
            M = NA
            IF (ABS(TW(NE,NA)).GT.ABS(TW(NA,NE))) THEN
               TW(NA,NA) = - (TW(NE,NE)-P)/TW(NE,NA)
               TW(NA,NE) = -Q/TW(NE,NA)
            ELSE
               CALL EB13PD(-TW(NA,NE),ZERO,TW(NA,NA)-P,Q,E,F)
               TW(NA,NA) = E
               TW(NA,NE) = F
            END IF
            TW(NE,NA) = ONE
            TW(NE,NE) = ZERO
            IF (NA.LT.2) GO TO 50
            DO 30 II = 2,NA
               I = NA - II + 1
               WR = ER(I)
               WI = EI(I)
               IF (WI.EQ.ZERO) THEN
                  U1 = TW(I,I) - P
                  RA = TW(I,NE) + DDOT(NA+1-M,TW(I,M),NUMCOL,TW(M,NA),1)
                  SA = DDOT(NA+1-M,TW(I,M),NUMCOL,TW(M,NE),1)
                  CALL EB13PD(-RA,-SA,U1,Q,E,F)
                  TW(I,NA) = E
                  TW(I,NE) = F
                  M = I
               ELSE IF (WI.GT.ZERO) THEN
                  Z = TW(I+1,I+1) - P
                  U1 = TW(I,I) - P
                  R = TW(I+1,NE) + DDOT(NA+1-M,TW(I+1,M),NUMCOL,
     +                TW(M,NA),1)
                  S = DDOT(NA+1-M,TW(I+1,M),NUMCOL,TW(M,NE),1)
                  RA = TW(I,NE) + DDOT(NA+1-M,TW(I,M),NUMCOL,TW(M,NA),1)
                  SA = DDOT(NA+1-M,TW(I,M),NUMCOL,TW(M,NE),1)
                  X1 = TW(I,I+1)
                  Y1 = TW(I+1,I)
                  VR = (WR-P)**2 + WI*WI - Q*Q
                  VI = (WR-P)*TWO*Q
                  IF (ABS(VR)+ABS(VI).EQ.ZERO) VR = XN*RMACH*
     +                (ABS(U1)+ABS(Q)+ABS(X1)+ABS(Y1)+ABS(Z))
                  CALL EB13PD(X1*R-Z*RA+Q*SA,X1*S-Z*SA-Q*RA,VR,VI,E,F)
                  TW(I,NA) = E
                  TW(I,NE) = F
                  IF (ABS(X1).GT.ABS(Z)+ABS(Q)) THEN
                     TW(I+1,NA) = (-RA-U1*TW(I,NA)+Q*TW(I,NE))/X1
                     TW(I+1,NE) = (-SA-U1*TW(I,NE)-Q*TW(I,NA))/X1
                  ELSE
                     CALL EB13PD(-R-Y1*TW(I,NA),-S-Y1*TW(I,NE),Z,Q,E,F)
                     TW(I+1,NA) = E
                     TW(I+1,NE) = F
                  END IF
                  M = I
               END IF
   30       CONTINUE
C
         ELSE IF (Q.EQ.ZERO) THEN
C Real vector
            M = NE
            TW(NE,NE) = ONE
            IF (NA.LT.1) GO TO 50
            DO 40 II = 1,NA
               I = NA + 1 - II
               WR = ER(I)
               WI = EI(I)
               IF (WI.EQ.ZERO) THEN
                  WW = TW(I,I) - P
                  R = TW(I,NE) + DDOT(NA+1-M,TW(I,M),NUMCOL,TW(M,NE),1)
                  IF (WW.EQ.ZERO) WW = RMACH*XN
                  TW(I,NE) = -R/WW
                  M = I
               ELSE IF (WI.GT.ZERO) THEN
                  U1 = TW(I,I) - P
                  Z = TW(I+1,I+1) - P
                  S = TW(I+1,NE) + DDOT(NA+1-M,TW(I+1,M),NUMCOL,
     +                TW(M,NE),1)
                  R = TW(I,NE) + DDOT(NA+1-M,TW(I,M),NUMCOL,TW(M,NE),1)
                  Y1 = TW(I+1,I)
                  X1 = TW(I,I+1)
                  Q = (WR-P)**2 + WI*WI
                  T1 = (X1*S-Z*R)/Q
                  TW(I,NE) = T1
                  IF (ABS(X1).GT.ABS(Z)) THEN
                     TW(I+1,NE) = (-R-U1*T1)/X1
                  ELSE
                     TW(I+1,NE) = (-S-Y1*T1)/Z
                  END IF
                  M = I
               END IF
   40       CONTINUE
         END IF
   50 CONTINUE
C
C Multiply by transformation matrix to give vectors of original
C full matrix
      CALL DGEMM('N','N',N,NUMEIG,NUMEIG,ONE,X,LN,TW,NUMCOL,ZERO,Y,LN)
C
C Normalise eigenvectors
      I = 1
   60 IF (I.LE.NUMEIG) THEN
         IF (EI(I).EQ.ZERO) THEN
            Z = DNRM2(N,Y(1,I),1)
            IF (Z.GT.RMACH) CALL DSCAL(N,ONE/Z,Y(1,I),1)
            I = I + 1
         ELSE
            Z = DDOT(N,Y(1,I),1,Y(1,I),1) +
     +          DDOT(N,Y(1,I+1),1,Y(1,I+1),1)
            Z = SQRT(Z)
            IF (Z.GT.RMACH) THEN
               CALL DSCAL(N,ONE/Z,Y(1,I),1)
               CALL DSCAL(N,ONE/Z,Y(1,I+1),1)
            END IF
            I = I + 2
         END IF
         GO TO 60
      END IF
C
C Compute residuals. For this it is necessary to ask the user to
C compute A*Y.
      CALL EB13ED(N,Y,LN,W,LN,NUMEIG)
C
      RETURN
      END
C January 2002: PA02 (from archive) included inline
      SUBROUTINE EB14AD(N,M,Z,ZOUT)
C
C Determines the convex hull of a set of points in the complex plane.
C This routine is for the case where Z is a symmetric set with respect
C to the real axis and determines only the upper 1/2-part of the convex
C hull.
C
C    Arguments( an * indicates the arguement is changed on exit)
C
C  N      Integer variable. On entry set to number of points in Z.
C *M      Integer variable. Not set on entry.
C         on exit set to number of points of the convex hull
C Z(2,N)  Real/DP array of length 2*N. On entry contains the
C         coordinates of the points for which we are seeking
C         the convex hull with
C         Z(1,I) = real parts,  Z(2,I) = imag. parts
C *ZOUT(2,N)  Real/DP array of length 2*N. Not set on entry.
C             On exit contains the points of the convex hull
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER M,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION Z(2,N),ZOUT(2,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CM,CS,RMIN
      INTEGER I1,I2,I3,J,K
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ZTMP(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL EB14ED
C     ..
      M = 1
      IF (N.EQ.1) THEN
         ZOUT(1,1) = Z(1,1)
         ZOUT(2,1) = Z(2,1)
         RETURN
      END IF
C Determine the first point
      RMIN = Z(1,1)
      I1 = 1
      DO 10 J = 2,N
         IF (Z(1,J).LE.RMIN .AND. Z(2,J).GE.ZERO) THEN
            RMIN = Z(1,J)
            I1 = J
         END IF
   10 CONTINUE
      ZOUT(1,M) = Z(1,I1)
      ZOUT(2,M) = Z(2,I1)
      IF (N.LE.2) RETURN
C Determine second  point
      ZTMP(1) = Z(1,I1)
      ZTMP(2) = Z(2,I1) - 1.0D0
      CM = 2.0D0
      I2 = I1
      DO 20 K = 1,N
         IF (Z(1,K).GT.RMIN) THEN
            CALL EB14ED(ZTMP,Z(1,I1),Z(1,K),CS)
            IF (CS.LE.CM) THEN
               IF (CS.NE.CM .OR. Z(1,I2).LE.Z(1,K)) THEN
                  CM = CS
                  I2 = K
               END IF
            END IF
         END IF
   20 CONTINUE
C Second  point determined.
      M = M + 1
      ZOUT(1,M) = Z(1,I2)
      ZOUT(2,M) = Z(2,I2)
      RMIN = ZOUT(1,M)
C Determine all other points.
   30 CM = 2.0D0
      I3 = I2
      DO 40 K = 1,N
         IF (Z(1,K).GT.RMIN) THEN
            CALL EB14ED(Z(1,I1),Z(1,I2),Z(1,K),CS)
            IF (CS.LE.CM) THEN
               IF (CS.NE.CM .OR. Z(1,I3).LE.Z(1,K)) THEN
                  CM = CS
                  I3 = K
               END IF
            END IF
         END IF
   40 CONTINUE
      IF (I3.NE.I2) THEN
         M = M + 1
         ZOUT(1,M) = Z(1,I3)
         ZOUT(2,M) = Z(2,I3)
         RMIN = ZOUT(1,M)
         I1 = I2
         I2 = I3
         GO TO 30
      END IF
      RETURN
      END
C********************************************************************
      SUBROUTINE EB14BD(MHULL,XR,XI,CH,D,C2,A2,A1)
C
C Find the best ''ellipse'' which contains the unwanted eigenvalues.
C This is Ho's algorithm.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A1,A2,C2,D,XI,XR
      INTEGER MHULL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CH(2,MHULL)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,A1S,AS,B,B1,B1S,B2,BS,C2S,DS,EPSIL,R,RS,S1,S2,
     +                 T1,T2,TEST,TEST1,X1I,X1R,X2I,X2R,X3I,X3R,Y,Z
      INTEGER I,IER,IG,IOPT1,IOPT2,IOPT3,J,JJ
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EB14ID
      EXTERNAL EB14ID
C     ..
C     .. External Subroutines ..
      EXTERNAL EB14JD,EB14KD
C     ..
      EPSIL = 1.D-4
C
      X2R = CH(1,MHULL)
      X2I = CH(2,MHULL)
      X1R = CH(1,1)
      X1I = CH(2,1)
C
      IF (MHULL.LE.2) THEN
         X3R = (X1R+X2R)/2.0D0
         X3I = (X1I+X2I)/2.0D0
      ELSE
         Z = -1.0D0
         DO 10 I = 2,MHULL - 1
            IF (CH(2,I).GT.Z) THEN
               Z = CH(2,I)
               IG = I
            END IF
   10    CONTINUE
         X3R = CH(1,IG)
         X3I = CH(2,IG)
      END IF
C
C Fit 3-point ellipse through (1R,X1I), (X2R,X2I), (X3R,X3I).
   20 CALL EB14JD(X1R,X1I,X2R,X2I,X3R,X3I,XR,XI,D,A,B,A1,B1,R,C2,IER)
C
      IF (IER.NE.1) THEN
C 3-point ellipse found. Test to see if it is feasible.
         Y = EPSIL
         DO 30 J = 1,MHULL
            S1 = CH(1,J)
            T1 = CH(2,J)
C Skip over points defining the ellipse.
            IF ((S1.EQ.X1R) .AND. (T1.EQ.X1I)) GO TO 30
            IF ((S1.EQ.X2R) .AND. (T1.EQ.X2I)) GO TO 30
            IF ((S1.EQ.X3R) .AND. (T1.EQ.X3I)) GO TO 30
            TEST = EB14ID(S1,T1,D,A,B)
            IF (TEST.GT.Y) THEN
C 3-point ellipse not feasible. Record the point which lies furtherest
C from the ellipse.
               Y = TEST
               JJ = J
            END IF
   30    CONTINUE
C
C If 3-point ellipse is not feasible, replace X3R,X3I.
         IF (Y.NE.EPSIL) THEN
            X3R = CH(1,JJ)
            X3I = CH(2,JJ)
            GO TO 20
         END IF
C Feasible 3-point ellipse found.
      END IF
C
C Compute the optimal ellipse which passes through (X1R,X1I)
C and (X2R,X2I) relative to (XR,XI)
   40 IOPT1 = 0
      IOPT2 = 0
      IOPT3 = 0
C      BS = B
      CALL EB14KD(X1R,X1I,X2R,X2I,XR,XI,RS,AS,BS,A1S,B1S,DS,C2S)
      TEST1 = EB14ID(X3R,X3I,DS,AS,BS)
C
      IF (TEST1.LE.EPSIL) THEN
C 2-point ellipse encloses the point (X3R,X3I)
         IOPT1 = 1
      ELSE
C Compute the other 2 point ellipse
C Choose (S1,T1) according to which side of the ellipse centre
C the point (X3R,X3I) lies.
         IF (X3R.GT.D) THEN
C (X3R,X3I) lies on same side as (X1R,X1I)
            S1 = X1R
            T1 = X1I
            S2 = X2R
            T2 = X2I
         ELSE
C (X3R,X3I) lies on same side as (X2R,X2I)
            S1 = X2R
            T1 = X2I
            S2 = X1R
            T2 = X1I
         END IF
C         BS = B
C Compute the optimal ellipse which passes through (S1,T1)
C and (X3R,X3I) relative to (XR,XI)
         CALL EB14KD(S1,T1,X3R,X3I,XR,XI,RS,AS,BS,A1S,B1S,DS,C2S)
         TEST1 = EB14ID(S2,T2,DS,AS,BS)
         IF (TEST1.LE.EPSIL) THEN
C 2-point ellipse encloses the third point (S2,T2)
            IOPT2 = 1
         ELSE
            IOPT3 = 1
         END IF
      END IF
C
      IF (IOPT1.EQ.1 .OR. IOPT2.EQ.1) THEN
C One of the 2-point ellipses encloses the third point.
C We now have to check whether this ellipse is feasible
C (whether it encloses all the other points on the hull).
C If this 2-point ellipse is feasible then it is optimal.
         Y = EPSIL
         DO 50 J = 1,MHULL
            S2 = CH(1,J)
            T2 = CH(2,J)
C Skip the points defining the 3-point ellipse
            IF ((S2.EQ.X1R) .AND. (T2.EQ.X1I)) GO TO 50
            IF ((S2.EQ.X2R) .AND. (T2.EQ.X2I)) GO TO 50
            IF ((S2.EQ.X3R) .AND. (T2.EQ.X3I)) GO TO 50
            TEST = EB14ID(S2,T2,DS,AS,BS)
            IF (TEST.GT.Y) THEN
C 2-point ellipse not feasible. Record the point which lies furtherest
C from the ellipse.
               Y = TEST
               JJ = J
            END IF
   50    CONTINUE
C
         IF (Y.NE.EPSIL) THEN
C 2-point ellipse not feasible.
            IF (IOPT1.EQ.1) THEN
C Swop the point (X3R,X3I) for the point on the hull
C which lies furtherest from the ellipse.
               X3R = CH(1,JJ)
               X3I = CH(2,JJ)
            ELSE
               IF (X1R.EQ.S1) THEN
C Swop the point (X2R,X2I) for the point on the hull
C which lies furtherest from the ellipse.
                  X2R = CH(1,JJ)
                  X2I = CH(2,JJ)
                  IF (X3R.GT.X2R) THEN
C Interchange X2 and X3 (so that X1R < X2R < X3R)
                     S1 = X2R
                     X2R = X3R
                     X3R = S1
                     S1 = X2I
                     X2I = X3I
                     X3I = S1
                  END IF
               ELSE
C Swop the point (X1R,X1I) for the point on the hull
C which lies furtherest from the ellipse.
                  X1R = CH(1,JJ)
                  X1I = CH(2,JJ)
                  IF (X3R.LT.X1R) THEN
C Interchange X1 and X3 (so that X1R < X2R < X3R)
                     S1 = X1R
                     X1R = X3R
                     X3R = S1
                     S1 = X1I
                     X1I = X3I
                     X3I = S1
                  END IF
               END IF
            END IF
C Compute new (feasible) 3-point ellipse
            CALL EB14JD(X1R,X1I,X2R,X2I,X3R,X3I,XR,XI,D,A,B,A1,B1,R,C2,
     +                  IER)
            GO TO 40
         END IF
      END IF
C
      IF (IOPT3.NE.1) THEN
C Neither 2-point ellipse is feasible and they do not enclose
C the third point so the 3-point ellipse is optimal.
         D = DS
         A = AS
         B = BS
         R = RS
         A1 = A1S
         C2 = C2S
         B1 = B1S
      END IF
      A2 = A*A
      B2 = B*B
      IF (B2.GT.A2) C2 = -C2
C
      RETURN
      END
C
C****************************************************************
      SUBROUTINE EB14CD(IVP,XR,EV,D,C2,A2,A1)
C
C This is Braconnier's algorithm for determining a good elipse
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A1,A2,C2,D,XR
      INTEGER IVP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EV(2,IVP)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,BB,DB,EPS,PP,RR,TB,TEST
      INTEGER I
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EB14ID
      EXTERNAL EB14ID
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      EPS = 0.00001D0
      PP = -1.D0
      D = (EV(1,1)+EV(1,IVP))/2.D0
      A = EV(1,IVP) - D
   10 PP = PP + 1.D0
      IF (PP.GT.EPS) A = A + ABS((EV(1,IVP)- (PP*XR))/ (PP+1.D0))
      B = 0.D0
C  Determine B
      DO 20 I = 1,IVP
         RR = EV(1,I) - D
         TB = EV(2,I)*A
         DB = SQRT(ABS(A*A-RR*RR))
         IF (DB.GT.1.D-4) THEN
            BB = ABS(TB/DB)
            IF (BB.GT.B) B = BB
         END IF
   20 CONTINUE
      DO 30 I = 1,IVP
         TEST = EB14ID(EV(1,I),EV(2,I),D,A,B)
         IF (TEST.GT.EPS) GO TO 10
   30 CONTINUE
      C2 = ABS(A*A-B*B)
      IF (B*B.GT.A*A) C2 = -C2
      A2 = A*A
      A1 = XR - D
      IF (A1.LT.A) A1 = A
      RETURN
      END
C********************************************************************
      SUBROUTINE EB14DD(MHULL,NUMEIG,EV,CH,ER,EI,WORK,D,C2,A2)
C
C This is Saad's algorithm for finding an ellipse
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A2,C2,D
      INTEGER MHULL,NUMEIG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CH(2,MHULL),EI(NUMEIG),ER(NUMEIG),EV(2,MHULL),
     +                 WORK(NUMEIG)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A2OLD,C2OLD,DOLD,EPSL,RATIO,RC,TA2,TC2,TD,TEST,
     +                 TRC,WRT1,XREF
      INTEGER I,IFLAG,J,K,L
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL EB14FD,EB14GD,EB14HD
C     ..
C     .. Intrinsic Functions ..
C
      INTRINSIC ABS,MAX,SQRT
C     ..
      XREF = ER(NUMEIG)
      IF (C2.NE.ZERO .AND. EI(NUMEIG).NE.ZERO) THEN
         CALL EB14FD(NUMEIG,ER,EI,WORK,RATIO,WRT1,D,C2,A2)
         XREF = D + WRT1*SQRT(ABS(C2))
      END IF
C Shift convex hull so that ellipse will not contain the origin
      K = MHULL
      DO 10 J = 1,MHULL
         EV(1,J) = XREF - CH(1,K)
         EV(2,J) = CH(2,K)
         K = K - 1
   10 CONTINUE
      A2OLD = A2
      C2OLD = C2
      DOLD = D
C
      IF (MHULL.EQ.1) THEN
         A2 = ZERO
         C2 = -EV(2,1)**2
         D = EV(1,1)
         RETURN
      END IF
      IF (MHULL.EQ.2 .AND. EV(1,1).EQ.EV(1,2)) THEN
C Complex hull is straight line parallel to imaginary axis
         A2 = ZERO
         C2 = MAX(EV(2,1),EV(2,2))
         C2 = -C2**2
         D = EV(1,1)
         GO TO 90
      END IF
C
      TRC = FD15AD('H')
      EPSL = SQRT(FD15AD('E'))
      TD = ZERO
      TC2 = ZERO
      TA2 = ZERO
C
C First see if any pairwise best is the optimum
C
      DO 40 I = 1,MHULL
         DO 30 J = 1,MHULL
            IF (I.EQ.J) GO TO 30
            CALL EB14GD(MHULL,EV,I,J,D,C2,A2,IFLAG)
            IF (IFLAG.EQ.0) THEN
               DO 20 K = 1,MHULL
                  IF (K.NE.I .AND. K.NE.J) THEN
C TEST > 0 means point EV(*,K) is outside ellipse ---> pairwise best
C is not optimal. --> check next pairwise best.
                     IF (EV(2,K).NE.ZERO) THEN
                        TEST = (D-EV(1,K))* (D-EV(1,K))* (A2-C2) +
     +                         EV(2,K)*EV(2,K)*A2 - A2* (A2-C2)
                        IF (TEST.GT.EPSL*A2* (A2-C2)) GO TO 30
                     ELSE
                        TEST = (D-EV(1,K))* (D-EV(1,K)) - A2
                        IF (TEST.GT.EPSL*A2) GO TO 30
                     END IF
                  END IF
   20          CONTINUE
C We have found a pairwise best - ready to return
               GO TO 90
            END IF
   30    CONTINUE
   40 CONTINUE
C
C No pairwise best has been found - look for three point best.
C Try all triples J,K,L.
C
      DO 80 J = 1,MHULL
         DO 70 K = 1,MHULL
            IF (K.EQ.J) GO TO 70
            DO 60 L = 1,MHULL
               IF (L.EQ.K .OR. L.EQ.J) GO TO 60
C  Find the three point ellipse passing by points J,K,L
               CALL EB14HD(MHULL,EV,J,K,L,D,C2,A2,IFLAG,RC)
               IF (IFLAG.EQ.0) THEN
C Test if this three-point  is a candidate
                  DO 50 I = 1,MHULL
                     IF (I.NE.J .AND. I.NE.K .AND. I.NE.L) THEN
                        IF (EV(2,I).NE.ZERO) THEN
                           TEST = (D-EV(1,I))* (D-EV(1,I))* (A2-C2) +
     +                            EV(2,I)*EV(2,I)*A2 - A2* (A2-C2)
                           IF (TEST.GT.EPSL*A2* (A2-C2)) GO TO 60
                        ELSE
                           TEST = (D-EV(1,I))* (D-EV(1,I)) - A2
                           IF (TEST.GT.EPSL*A2) GO TO 60
                        END IF
                     END IF
   50             CONTINUE
C This ellipse is admissible. now compare to other AR'S
                  IF (RC.LT.TRC) THEN
C Note: TR is the conv. coef.  with current ellipse.
C -->  this three-point is best so far
                     TRC = RC
                     TD = D
                     TC2 = C2
                     TA2 = A2
                  END IF
               END IF
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
C
C The best three-point fit has been found
      D = TD
      C2 = TC2
      A2 = TA2
   90 IF (D.EQ.ZERO .AND. A2.EQ.ZERO .AND. C2.EQ.ZERO) THEN
         A2 = A2OLD
         C2 = C2OLD
         D = DOLD
      ELSE
         D = XREF - D
      END IF
      RETURN
      END
C******************************************************************
      SUBROUTINE EB14ED(Z1,Z2,Z3,CS)
C
C Computes the cosine (CS) of the angle between
C the elements Z1, Z2 and Z3:
C                                      * Z2
C                                    .      .
C                                  .          .
C Z2 is the middle point...      * Z1           * Z3
C     .. Scalar Arguments ..
      DOUBLE PRECISION CS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION Z1(2),Z2(2),Z3(2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DEN,DEN21,DEN23,RMACH,X21,X23,Y21,Y23
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      RMACH = FD15AD('E')
      X21 = Z1(1) - Z2(1)
      X23 = Z3(1) - Z2(1)
      Y21 = Z1(2) - Z2(2)
      Y23 = Z3(2) - Z2(2)
      DEN21 = SQRT(X21**2+Y21**2)
      DEN23 = SQRT(X23**2+Y23**2)
      DEN = DEN21*DEN23
      IF (DEN.LE.RMACH) THEN
         CS = -1.0D0
      ELSE
         CS = (X21*X23+Y21*Y23)/DEN
      END IF
      RETURN
      END
C******************************************************************
      SUBROUTINE EB14FD(NEV,ER,EI,VEC,RATIO,WRT1,D,C2,A2)
C
C     Returns convergence ratios for the eigenvalues ER,EI, relative to
C     the ellipse defined by D,C2.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO,HALF
      PARAMETER (ONE=1.0D0,ZERO=0.0D0,HALF=0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A2,C2,D,RATIO,WRT1
      INTEGER NEV
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EI(NEV),ER(NEV),VEC(NEV)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,EE,SMALL,T,X,Y
      INTEGER I
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL KB06AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      SMALL = SQRT(FD15AD('E'))
      EE = SQRT(ABS(C2))
      IF (EE.LT.SMALL) EE = SMALL
      DO 10 I = 1,NEV
         X = (ER(I)-D)/EE
         Y = EI(I)/EE
         IF (C2.LT.ZERO) THEN
            T = X
            X = Y
            Y = T
         END IF
C Calculate major semi-axis of ellipse centre D, focus SQRT(C2)
C passing through (ER(I), EI(I))
         A = HALF* (SQRT((X+ONE)**2+Y**2)+SQRT((X-ONE)**2+Y**2))
         VEC(I) = (A+SQRT(ABS(A**2-ONE)))
   10 CONTINUE
      T = A2/ (EE*EE)
      IF (C2.LT.ZERO) T = ONE + T
      T = SQRT(T) + SQRT(ABS(T-ONE))
C Order ratios (descending order).
      CALL KB06AD(VEC,NEV)
      RATIO = VEC(1)/VEC(NEV)
      WRT1 = (VEC(NEV)+ONE/VEC(NEV))*HALF
      RETURN
      END
C***************************************************************
      SUBROUTINE EB14GD(M,CH,J,K,D,C2,A2,IFLAG)
C
C This subroutine finds the optimum iteration parameters
C for the two points CH(-,J) and CH(-,K)
C
C  Arguments (an * indicates the variable is changed by the routine)
C  M      - Integer variable. on entry set to the second dimension of CH
C  CH     - Real (double precision) array of length 2*M.
C           on entry must contain the real and imaginary parts of
C           the eigenvalues.
C J,K     - Integer variables. on entry set to the numbers of the two
C           eigenvalues through which we want to find the pair-wise
C           best ellipse
C * IFLAG - Integer variable. not set on entry.
C           On exit a value IFLAG=0
C           indicates a successful return, and IFLAG=1 indicates the
C           pair-wise best ellipse for J,K does not exist.
C
C    On exit, D, C2, A2 are set to the iteration parameters
C
C     .. Parameters ..
      DOUBLE PRECISION TWO,ZERO
      PARAMETER (TWO=2.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A2,C2,D
      INTEGER IFLAG,J,K,M
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CH(2,M)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,AA,B,BB,EPS,S,T,TT,X,Y,Y1,Y2,Z
C     ..
C     .. External Subroutines ..
      EXTERNAL EB14MD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,SQRT
C     ..
      EPS = 0.001D0
C
C Initialise IFLAG
      IFLAG = 0
C
C Set the constants
      A = (CH(1,K)-CH(1,J))/TWO
      IF (A.LE.ZERO) THEN
         IFLAG = 1
         RETURN
      END IF
      B = (CH(1,K)+CH(1,J))/TWO
      S = (CH(2,K)-CH(2,J))/TWO
      T = (CH(2,K)+CH(2,J))/TWO
C
C Eliminate degenerate cases
      IF (T/B.LT.EPS) THEN
C     Case i: T small
         D = B
         C2 = A*A
         A2 = C2
      ELSE IF (A/T.LT.EPS) THEN
C Case ii: A small
         D = B
         Y = T + ABS(S)
         C2 = -Y*Y
         A2 = A* (A+SQRT(A*A+ (TWO*Y)**2))/TWO
      ELSE IF (ABS(S)/A.LT.EPS) THEN
C Case iii: S small
C Finds the best iteration parameters
C when both eigenvalues have the same imaginary part.
C A cubic is solved; we are looking for the only real root
C in (A**2, B**2).
C See Abramowitz and Stegun p17 (the real root of a third degree poly.)
         D = B
         TT = T*T
         BB = B*B
         AA = A*A
         X = (BB*AA**2*TT)* ((BB+TT)**2+AA* (TT-BB))/ (TWO* (BB+TT)**3)
         Z = (BB*AA**2*TT)/ ((BB+TT)**2)
         Y1 = X + SQRT(X**2+Z**3)
         Y2 = X - SQRT(X**2+Z**3)
         Y = SIGN(ABS(Y1)** (1.0D0/3.0D0),Y1) +
     +       SIGN(ABS(Y2)** (1.0D0/3.0D0),Y2)
         A2 = Y + (BB*AA)/ (BB+TT)
         C2 = (A2* (A2-TT-AA))/ (A2-AA)
      ELSE
C Case iv: general case
         CALL EB14MD(A,B,S,T,D,C2,A2,IFLAG)
      END IF
      RETURN
      END
C*****************************************************************
      SUBROUTINE EB14HD(M,CH,J,K,L,D,C2,A2,IFLAG,RC)
C
C This subroutine finds the iteration parameters
C associated with the unique ellipse through the three
C eigenvalues: CH(-,J),CH(-,K), and CH(-,L)
C
C  Arguments (an * indicates the variable is changed by the subroutine)
C
C  M       - integer variable. on entry set to the second
C            dimension of CH
C  CH      - real (double precision) array of length 2*M.
C            on entry must contain the real and imaginary parts
C            of the eigenvalues.
C J,K,L    - integer variables. on entry set to the numbers of the three
C            eigenvalues through which we want a three-way ellipse
C * IFLAG  - integer variable. Not set on entry. on exit a value IFLAG=0
C            indicates a successful return, and IFLAG=1 indicates the
C            three-way ellipse for J,K,L does not exist.
C * RC     - real (double precision) variable. not set on entry.
C           on exit, RC is set to the rate of convergence
C    On exit, D, C2, A2 are set to the ellipse parameters
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A2,C2,D,RC
      INTEGER IFLAG,J,K,L,M
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CH(2,M)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DET,DIFF,DKJ,DLJ,DLK,EPS,W,Z
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN,SQRT
C     ..
C     Initialise IFLAG
      IFLAG = 0
      EPS = FD15AD('E')*10**3
C
      DKJ = CH(1,K) - CH(1,J)
      DLK = CH(1,L) - CH(1,K)
      DLJ = CH(1,L) - CH(1,J)
      DIFF = MIN(DKJ,DLK,DLJ)
C If DIFF.LE.EPS then unstable. Set IFLAG and return
      IF (DIFF.LE.EPS) GO TO 10
C
C Compute D
      Z = CH(2,J)*CH(2,J)* (CH(1,K)-CH(1,L)) +
     +    CH(2,L)*CH(2,L)* (CH(1,J)-CH(1,K)) +
     +    CH(2,K)*CH(2,K)* (CH(1,L)-CH(1,J))
C If no ellipse fits these three points then return
      IF (Z.LE.EPS) GO TO 10

      W = CH(2,J)*CH(2,J)* (CH(1,K)*CH(1,K)-CH(1,L)*CH(1,L)) +
     +    CH(2,L)*CH(2,L)* (CH(1,J)*CH(1,J)-CH(1,K)*CH(1,K)) +
     +    CH(2,K)*CH(2,K)* (CH(1,L)*CH(1,L)-CH(1,J)*CH(1,J))
C
      D = W/ (TWO*Z)
C If ellipse contains the origin then return
      IF (D.LE.EPS) GO TO 10
C
C Compute other parameters
      A2 = (CH(2,J)*CH(2,J)*CH(1,L)*CH(1,K)* (CH(1,L)-CH(1,K))+
     +     CH(2,L)*CH(2,L)*CH(1,J)*CH(1,K)* (CH(1,K)-CH(1,J))+
     +     CH(2,K)*CH(2,K)*CH(1,J)*CH(1,L)* (CH(1,J)-CH(1,L)))/Z + D*D
C
C If ellipse contains the origin then return
      IF ((A2-D*D).GE.EPS) GO TO 10
      IF (A2.LE.EPS) THEN
         A2 = ZERO
         C2 = -MAX(CH(2,J),CH(2,K),CH(2,K))
      ELSE
         DET = (CH(1,J)-CH(1,K))* (CH(1,K)-CH(1,L))* (CH(1,L)-CH(1,J))
         IF (DET.LE.EPS) GO TO 10
         C2 = A2* (ONE-Z/DET)
      END IF
C
      RC = (SQRT(A2)+SQRT(A2-C2))/ (D+SQRT(D*D-C2))
      GO TO 20

   10 IFLAG = 1
   20 RETURN
      END
C*************************************************************
      DOUBLE PRECISION FUNCTION EB14ID(X,Y,D,A,B)
C
C     .. Parameters ..
      DOUBLE PRECISION EPS,ZERO
      PARAMETER (EPS=1.0D-4,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,D,X,Y
C     ..
      IF (A.EQ.ZERO) A = EPS
      IF (B.EQ.ZERO) B = EPS
      EB14ID = ((X-D)* (X-D)/ (A*A)) + (Y*Y/ (B*B)) - 1.0D0
C
      RETURN
      END
C**************************************************************
      SUBROUTINE EB14JD(X1,Y1,X2,Y2,X3,Y3,XR,YR,D,A,B,A1,B1,R,C2,IER)
C
C Compute the ellipse that passes through 3 points (X1,Y1), (X2,Y2)
C and (X3,Y3) with X1<X3<X2.
C R    : Rate of convergence relative to referenced point (XR,YR)
C D    : Centre of the ellipse.
C A,B  : Semi-major and minor axis of the calculated ellipse.
C A1,B1: Semi-major and minor axis of the referenced ellipse.
C IER=0: For no error
C IER=1: No possible ellipse, exit
C IER=2: Ellipse contains reference point. A1,B1,R undefined.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,A1,B,B1,C2,D,R,X1,X2,X3,XR,Y1,Y2,Y3,YR
      INTEGER IER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A12,A2,AAA,ASB2,B2,BETA,X12D,X12M,X13D,X13M,X1D,
     +                 XD,Y212,Y312
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
      IER = 0
      X12M = (X1+X2)/2.0D0
      X13M = (X1+X3)/2.0D0
      Y212 = (Y2-Y1)* (Y2+Y1)/2.0D0
      Y312 = (Y3-Y1)* (Y3+Y1)/2.0D0
      X12D = X1 - X2
      X13D = X1 - X3
      IF ((X12D.EQ.ZERO) .OR. (X13D.EQ.ZERO)) GO TO 10
      AAA = (Y212/X12D) - (Y312/X13D)
      IF (AAA.LE.ZERO) GO TO 10
      ASB2 = (X12M-X13M)/AAA
      D = X12M - (ASB2*Y212/X12D)
      X1D = X1 - D
      B2 = ((X1D*X1D)/ASB2) + (Y1*Y1)
      A2 = ASB2*B2
      A = SQRT(A2)
      B = SQRT(B2)
      XD = XR - D
      IF (A.GE.B) THEN
         C2 = A2 - B2
         BETA = (C2+ (XD*XD)+ (YR*YR))/2.0D0
         A12 = BETA + SQRT(BETA*BETA- (C2*XD*XD))
         IF (A12.LT.A2) GO TO 20
         B1 = SQRT(A12-C2)
      ELSE
         C2 = B2 - A2
         BETA = (-C2+ (XD*XD)+ (YR*YR))/2.0D0
         A12 = BETA + SQRT(BETA*BETA+ (C2*XD*XD))
         B1 = SQRT(A12+C2)
      END IF
      A1 = SQRT(A12)
C
      R = (A+B)/ (A1+B1)
      GO TO 30
   10 IER = 1
      GO TO 30
C
   20 IER = 2
   30 CONTINUE
      RETURN
      END
C
C**************************************************************
      SUBROUTINE EB14KD(X1,Y1,X2,Y2,X,Y,R,A,B,A1,B1,D,C2)
C
C Compute the optimal ellipse which passes through lambda(X1,Y1)
C and lambda(X2,Y2) relative to lambdaref(X,Y)
C R: Rate of convergence
C A, B : Semi-major and minor axis of the ellipse
C A1,B1: Semi-major and minor axis of the referenced  ellipse
C D: Center of the ellipse
C C2=A2-B2=A12-B12
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,A1,B,B1,C2,D,R,X,X1,X2,Y,Y1,Y2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A12,A1K1,A1K2,A1PR,A1SB1,A2,A3,ALFA,AM,APR,ASB2,
     +                 B1D,B1K1,B1K2,B1PR,B2,B3,BETA,BM,DR,X1D,X2D,XD,
     +                 XX,Y1SX1,YY,YY1,Z11A,ZZ,ZZ1
      INTEGER II,INV
C     ..
C     .. External Subroutines ..
      EXTERNAL EB14LD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SQRT
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
      INV = 0
      IF (X2.LT.X1) THEN
         XX = X2
         X2 = X1
         X1 = XX
         XX = Y2
         Y2 = Y1
         Y1 = XX
         INV = 1
      END IF
      IF (X2.EQ.X1) THEN
         D = X1
         B = MAX(Y1,Y2)
         A = ZERO
         C2 = B*B
         XD = X - D
         BETA = (-C2+ (XD*XD)+ (Y*Y))/2.0D0
         A12 = BETA + SQRT(BETA*BETA+ (C2*XD*XD))
         B1 = SQRT(A12+C2)
         A1 = SQRT(A12)
         R = B/ (A1+B1)
         GO TO 20
      END IF
      IF (Y2.EQ.Y1) THEN
         D = (X1+X2)/2.0D0
         IF (Y1.EQ.ZERO) THEN
            A = X2 - D
            B = ZERO
            C2 = A*A
            XD = X - D
            BETA = (C2+ (XD*XD)+ (Y*Y))/2.0D0
            A12 = BETA + SQRT(ABS(BETA*BETA- (C2*XD*XD)))
            B1 = SQRT(A12-C2)
            A1 = SQRT(A12)
            R = A/ (A1+B1)
         ELSE
            CALL EB14LD(X2,Y2,X,Y,D,A,B,A1,B1,R,C2)
         END IF
         GO TO 20
      END IF
      B = 10.0D0*Y1
      IF (Y2.GT.Y1) B = 10.0D0*Y2
      AM = X2 - X1
      BM = (Y2-Y1)* (Y2+Y1)/AM
C
      DO 10 II = 1,400
         XX = (B+Y2)* (B-Y2)
         ZZ = (AM*XX/BM) + ((XX*XX)/ (BM*BM))
         YY1 = - (XX/BM) + SQRT(ZZ)
         YY = - (XX/BM) - SQRT(ZZ)
         D = X2 - YY
         IF ((D.LT.X1) .OR. (D.GT.X2)) THEN
            YY = YY1
            D = X2 - YY
         END IF
         B2 = B*B
         A2 = (B2*YY*YY)/XX
         A = SQRT(A2)
         A3 = A2*A
         B3 = B2*B
         X1D = X1 - D
         X2D = X2 - D
         XD = X - D
         IF (A.GE.B) THEN
            C2 = A2 - B2
            BETA = (C2+ (XD*XD)+ (Y*Y))/2.0D0
            A12 = BETA + SQRT(ABS(BETA*BETA- (C2*XD*XD)))
            IF (A12.LT.A2) A12 = A2
            B1 = SQRT(A12-C2)
         ELSE
            C2 = B2 - A2
            BETA = (-C2+ (XD*XD)+ (Y*Y))/2.0D0
            A12 = BETA + SQRT(BETA*BETA+ (C2*XD*XD))
            IF (A12.LT.A2) A12 = A2
            B1 = SQRT(A12+C2)
         END IF
         A1 = SQRT(A12)
C
         R = (A+B)/ (A1+B1)
         YY1 = (Y1*Y1)/X1D
         YY = (Y2*Y2)/X2D
         ZZ = A3/B3
         Z11A = (XD*A)/ (A1*X1D)
         ALFA = (YY1-YY)/AM
         A1SB1 = A1/B1
         Y1SX1 = Y1/X1D
         ZZ1 = (Y*Y*A1SB1*A1SB1*A1SB1)/ (XD*XD)
         B1K1 = (((Y1SX1*Y1SX1)+ALFA)/Z11A) - A*ALFA/A1
         B1D = ZZ1 + (B1/A1)
         B1K1 = B1K1/B1D
         B1K2 = B/ (A1*B1D)
         B1PR = B1K2 + (ZZ*B1K1)
         APR = ALFA*ZZ
C
         A1K1 = (((Y1SX1*Y1SX1)+ALFA)/Z11A) - (ZZ1*B1K1)
         A1K2 = (ZZ1*B/ (A1*B1D))
         A1PR = (ZZ*A1K1) - A1K2
         DR = (1.0D0+APR)/ (A1PR+B1PR)
         IF (ABS(R-DR).LE.1.D-4) GO TO 20
         ZZ = (1.0D0+R* (A1K2-B1K2))/ (R* (A1K1+B1K1)-ALFA)
         IF (ZZ.LT.ZERO) ZZ = -ZZ
         IF (ZZ.EQ.ZERO) ZZ = 1.D-3
         ASB2 = ZZ** (2.0D0/3.0D0)
         B2 = ((X1D*X1D)/ASB2) + (Y1*Y1)
         IF (Y2.GT.Y1) B2 = ((X2D*X2D)/ASB2) + (Y2*Y2)
         B = SQRT(B2)
C
   10 CONTINUE
   20 CONTINUE
C
      IF (INV.EQ.1) THEN
         XX = X2
         X2 = X1
         X1 = XX
         XX = Y2
         Y2 = Y1
         Y1 = XX
      END IF
      RETURN
      END
C**************************************************************
      SUBROUTINE EB14LD(X,Y,X1,Y1,D,A,B,A1,B1,R,C2)
C
C Compute the optimal ellipse that passes through lambda(X,Y)
C relative to lambdaref(X1,Y1) given D the center of the ellipse.
C Ouput :
C A,B : Semi-major and minor axis of the ellipse
C A1,B1 : Semi-major and minor axis of the referenced ellipse
C R: Rate of convergence
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,A1,B,B1,C2,D,R,X,X1,Y,Y1
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A12,A2,ASB2,AX,AX1,B2,BETA,DR,R1,R2,RP,XSY2,YSX12
      INTEGER I,NUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
      X1 = X1 - D
      X = X - D
      YSX12 = (Y1*Y1)/ (X1*X1)
      XSY2 = (X*X)/ (Y*Y)
      ASB2 = XSY2** (2.0D0/3.0D0)
      B = SQRT(Y*Y+ (X*X/ASB2))
      AX = 1.0D0
      NUM = 40
      RP = 10.0D0
C
      DO 10 I = 1,NUM
         ASB2 = (AX*XSY2)** (2.0D0/3.0D0)
         B2 = Y*Y + (X*X/ASB2)
         A2 = B2*ASB2
         A = SQRT(A2)
         B = SQRT(B2)
         IF (A.GE.B) THEN
            C2 = A2 - B2
            BETA = (C2+ (X1*X1)+ (Y1*Y1))/2.0D0
            A12 = BETA + SQRT(BETA*BETA- (C2*X1*X1))
            IF (A12.LT.A2) A12 = A2
            B1 = SQRT(A12-C2)
         ELSE
            C2 = B2 - A2
            BETA = (-C2+ (X1*X1)+ (Y1*Y1))/2.0D0
            A12 = BETA + SQRT(BETA*BETA+ (C2*X1*X1))
            B1 = SQRT(A12+C2)
         END IF
C
         A1 = SQRT(A12)
         AX = ((A/B)**3.0D0)/XSY2
         AX1 = ((A1/B1)**3.0D0)*YSX12
         IF ((AX1.EQ.1.0D0) .AND. (AX.EQ.1.0D0)) THEN
            R = (A+B)/ (A1+B1)
            GO TO 20
         END IF
         R1 = (1.0D0-AX)/ (A*AX+B)
         R2 = (A1*AX1+B1)/ (1.0D0-AX1)
         DR = R1*R2
         R = (A+B)/ (A1+B1)
         IF ((ABS(R-DR).LE.1.D-3) .OR. (ABS(R-RP).LE.1.D-3)) GO TO 20
         RP = R
         R = R/R2
         AX = (1.0D0-B*R)/ (1.0D0+R*A)
         IF (AX.LE.ZERO) AX = ABS(AX)
   10 CONTINUE
   20 CONTINUE
      X1 = X1 + D
      X = X + D
      RETURN
      END
C********************************************************************
      SUBROUTINE EB14MD(A,B,S,T,D,C2,A2,IFLAG)
C
C  A,B,S,T REAL (DOUBLE PRECISION) variables which define the
C          fifth degree polynomial. Unchanged on exit.
C  IFLAG   INTEGER variable. Not set on entry. On exit, IFLAG=0
C          denotes a successful return; IFLAG=1 denotes an
C          unsuccessful return from the routine ZEROIN called
C          by the routine to find the real root of the fifth degree
C          polynomial.
C
C     This subroutine finds the optimal iteration parameters
C     for two complex eigenvalues. The best value of D is found
C     as the real root of a fifth degree polynomial with coefficients
C     P(1), P(2), P(3), P(4), P(5), P(6).
C     This routine is by Manteuffel. For details see thesis:
C     'An Iterative Method for Solving Nonsymmetric Linear Systems
C     with Dynamic Estimation of Parameters'. T A Manteuffel.
C     University of Illinois Report UIUCDCS-R-75-758 (1975).
C
C    On exit, D, C2, A2 are set to the iteration parameters
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,TWO,THREE,FOUR,FIVE
      PARAMETER (ZERO=0.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0,FIVE=5.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,A2,B,C2,D,S,T
      INTEGER IFLAG
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION E1,E2,PSCAL,R1,R2,R3,R4,Y,Z
      INTEGER N,NA,NB,NR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION E(5),P(6),PP(6),ROOT(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL PA02BD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN,SQRT
C     ..
      IFLAG = 0
C     Compute the fifth degree polynomial
      R1 = -A*T/S
      R2 = -A*S/T
      R3 = S*T/A
      R4 = -B
C
      PSCAL = ZERO
      P(1) = R1* (R1+TWO*R2+R3-FOUR*R4) + R2* (R2+R3-FOUR*R4) +
     +       R4* (FOUR*R4-TWO*R3)
      IF (P(1).GT.PSCAL) PSCAL = P(1)
      P(2) = R1* (FOUR*R4* (R3-R4)+R2* (THREE*FOUR*R4-FIVE*R3-FOUR*R2)+
     +       R1* (TWO*R4-R3-FOUR*R2)) + R2*
     +       (FOUR*R4* (R3-R4)+R2* (TWO*R4-R3)) - R3*R4*R4
      IF (P(2).GT.PSCAL) PSCAL = P(2)
      P(3) = R1* (R2* (R4* (TWO*R4-FOUR*R3)+R2* (THREE*R3-FOUR*R4))+
     +       R1* (R4*R4+R2* (FOUR*R2+THREE*R3-FOUR*R4)-R1*R3)) +
     +       R2*R2* (R4*R4-R2*R3)
      IF (P(3).GT.PSCAL) PSCAL = P(3)
      P(4) = R1*R2*R3* (R1* (THREE*R1-FOUR*R4)+R2* (THREE*R2-FOUR*R4)+
     +       TWO*R4*R4)
      IF (P(4).GT.PSCAL) PSCAL = P(4)
      P(5) = THREE*R1*R1*R2*R2*R3* (TWO*R4-R1-R2)
      IF (P(5).GT.PSCAL) PSCAL = P(5)
      P(6) = R1*R1*R2*R2*R3* (R1*R2-R4*R4)
      IF (P(6).GT.PSCAL) PSCAL = P(6)
C    Scale
      P(1) = P(1)/PSCAL
      P(2) = P(2)/PSCAL
      P(3) = P(3)/PSCAL
      P(4) = P(4)/PSCAL
      P(5) = P(5)/PSCAL
      P(6) = P(6)/PSCAL
C
C     Compute y and z
      Y = (T-S)* (T-S)* (B+A)* (B+A) - (T+S)* (T+S)* (B-A)* (B-A)
      Z = (T-S)* (T-S)* (B+A) - (T+S)* (T+S)* (B-A)
C
      IF (S.LT.ZERO) THEN
         E1 = MAX(-A, (Y/ (TWO*Z)-B))
         E2 = ZERO
      ELSE
         E1 = ZERO
         IF (Y.GT.TWO*B*Z) THEN
            E2 = A
         ELSE
            E2 = MIN(A,Y/ (TWO*Z)-B)
         END IF
      END IF
C
C Find  real root. The root is returned in E1.
C
      N = 5
C Reverse coefficient order for PA02B/BD.
      PP(1) = P(6)
      PP(2) = P(5)
      PP(3) = P(4)
      PP(4) = P(3)
      PP(5) = P(2)
      PP(6) = P(1)
      CALL PA02BD(PP,E1,E2,NB,NR,NA,N,ROOT,E)
      IF (NR.NE.1) THEN
         IFLAG = 1
         RETURN
      ELSE
         E1 = ROOT(1)
C        E(1) is an estimate of the error in the real root
         IF (E(1).GT.SQRT(FD15AD('E'))) THEN
            IFLAG = 1
         ELSE
C Compute iteration parameters
            D = E1 + B
            A2 = (E1-R1)* (E1-R2)
            C2 = A2* (E1-R3)/ (E1)
         END IF
      END IF
      RETURN
      END
      SUBROUTINE PA02BD(A,BX1,BX2,NB,NR,NA,N,ROOT,E)
C     ERROR ESTIMATES IN E( ) ARE RELATIVE.
C     .. Scalar Arguments ..
      DOUBLE PRECISION BX1,BX2
      INTEGER N,NA,NB,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),E(*),ROOT(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C1,C2,D1,DZ,E1,E2,E3,EPS,F1,FP,FZ,G1,P1,P2,P3,
     +                 SGN,STO,WLOG,WMAX,X,X1,X2,XXL,XXU,ZM,ZP,ZZ
      INTEGER I,ICR,IK,IR1,IS1,IS2,IS6,ISC,ISC1,ISW,J,M,M1,N2,NL,NR1,
     +        NRC,NRP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ALPHA(50),BETA(50),W1(50),W2(50)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MAX,EXP,SIGN,DBLE
C     ..
C     .. Executable Statements ..
      EPS = FD15AD('E')*10.0
      X1 = BX1
      X2 = BX2
      IF(N-1.LT.0) GO TO 113
      IF(N-1.GT.0) GO TO 220
      ROOT(1) = -A(1)/A(2)
      E(1) = EPS
      NR = 1
      GO TO 113

  220 M = N + 1
      IS6 = 1
      IF(X1-X2.LT.0) GO TO 30
      IF(X1-X2.EQ.0) GO TO 31
      STO = X1
      X1 = X2
      X2 = STO
C  ****** FORMING THE STURM SEQUENCE ******
   30 M1 = M + 1
      DO 1 J = 1,N
        I = M1 - J
        W1(J) = A(I)
        W2(J) = A(I)*DBLE(I-1)
    1 CONTINUE
      W1(M) = A(1)
  602 IS1 = 0
      IS2 = 0
      IF(W1(1)*W2(1).GT.0) GO TO 10
      IS2 = IS2 + 1
      GO TO 9

   10 IS1 = IS1 + 1
    9 ISW = 1
      N2 = M - 2
      NL = N2
      DO 2 I = 1,N2
        GO TO (3,4) ISW

    3   ALPHA(NL) = W1(1)/W2(1)
        BETA(NL) = (W1(2)-ALPHA(NL)*W2(2))/W2(1)
        W2(NL+2) = 0.0D0
        DO 5 J = 1,NL
          W1(J) = ALPHA(NL)*W2(J+2) + BETA(NL)*W2(J+1) - W1(J+2)
    5   CONTINUE
        IF(W1(1).EQ.0) GO TO 701
        ISW = 2
        GO TO 6

    4   ALPHA(NL) = W2(1)/W1(1)
        BETA(NL) = (W2(2)-ALPHA(NL)*W1(2))/W1(1)
        W1(NL+2) = 0.0D0
        DO 7 J = 1,NL
          W2(J) = ALPHA(NL)*W1(J+2) + BETA(NL)*W1(J+1) - W2(J+2)
    7   CONTINUE
        IF(W2(1).EQ.0) GO TO 701
        ISW = 1
    6   NL = NL - 1
        IF(W1(1)*W2(1).GE.0) GO TO 12
        IS2 = IS2 + 1
        GO TO 2

   12   IS1 = IS1 + 1
    2 CONTINUE
      GO TO (13,144) ISW

  701 NR = -1
      RETURN

   13 C2 = W1(2)
      C1 = W1(1)
      D1 = W2(1)
      GO TO 14

  144 C2 = W2(2)
      C1 = W2(1)
      D1 = W1(1)
   14 X = X1
C  ****** FINDING BRACKETS FOR EACH ROOT ******
   24 ISC = 0
      P1 = D1
      SGN = SIGN(1.D0,P1)
      P2 = C1*X + C2
      IF(SGN*P2.GE.0) GO TO 17
      ISC = ISC + 1
      SGN = -SGN
   17 DO 18 I = 1,N2
        P3 = (ALPHA(I)*X+BETA(I))*P2 - P1
        IF(SGN*P3.GE.0) GO TO 20
        ISC = ISC + 1
        SGN = -SGN
   20   P1 = P2
        P2 = P3
   18 CONTINUE
      GO TO (21,22,23) IS6

   21 ISC1 = ISC
      IS6 = 2
      X = X2
      GO TO 24

   22 IS6 = 3
      NR = ABS(ISC1-ISC)
      NA = ABS(IS2-ISC)
      NB = ABS(IS1-ISC1)
      IF(NR-1.LT.0) GO TO 113
      IF(NR-1.GT.0) GO TO 134
      W1(1) = X1
      W2(1) = X2
      NRC = 1
      GO TO 200

  134 NR1 = NR - 1
      DO 120 I = 1,NR1
        W2(I) = X1
  120 CONTINUE
      W1(1) = X1
      W2(NR) = X2
      XXL = X1
      XXU = X2
      ICR = 1
      NRC = NR
      NRP = NR
      GO TO 103

   23 NR = ABS(ISC1-ISC)
      IF(NR-ICR.GT.0) GO TO 102
      IF(NR-ICR.EQ.0) GO TO 101
      W1(ICR) = X
      XXL = X
      NRP = NR
      GO TO 103

  101 W1(ICR+1) = X
      W2(ICR) = X
      XXL = X
      ICR = ICR + 1
      IR1 = ICR
      DO 104 I = IR1,NRC
        IF(W2(I)-W1(1).LE.0) GO TO 104
        IF(I-ICR.LE.0) GO TO 107
        IK = I
        GO TO 109

  107   ICR = ICR + 1
        XXL = W1(I+1)
  104 CONTINUE
      GO TO 200

  109 XXU = W2(IK)
      NRP = IK
      GO TO 103

  102 IF(NR-NRP.EQ.0) GO TO 111
      W1(NR+1) = X
  111 W2(NR) = X
      NRP = NR
      XXU = X
  103 X = (XXU+XXL)*0.5D0
      GO TO 24
C  ****** LOCATING ROOTS (NEWTON-RAPHSON AND BISECTION) ******
  200 DO 202 I = 1,NRC
        ZZ = W1(I)
        ISW = 1
  304   DZ = D1
        FZ = C1*ZZ + C2
        DO 300 J = 1,N2
          P3 = (ALPHA(J)*ZZ+BETA(J))*FZ - DZ
          DZ = FZ
          FZ = P3
  300   CONTINUE
        GO TO (301,302,319,400) ISW

  301   IF(FZ.EQ.0) GO TO 317
        F1 = FZ
        G1 = DZ
        ZZ = W2(I)
        ISW = 2
        GO TO 304

  302   ISW = 3
        IF(ABS(F1)-ABS(FZ).GT.0) GO TO 306
        ZZ = W1(I)
        DZ = F1/G1
        ZM = W2(I)
        FZ = F1
        GO TO 324

  306   ZZ = W2(I)
        DZ = FZ/DZ
        ZM = W1(I)
  324   FP = FZ
        ZP = ZZ
        X = ZZ - DZ
        X2 = (ZM+ZP)*0.5D0
        IF((X2-X)* (X-ZZ).GE.0) GO TO 310
        ZZ = X2
        GO TO 325

  310   ZZ = X
        GO TO 325

  325   IF((ZP-ZZ)* (ZZ-ZM).LE.0) GO TO 317
        GO TO 304
  319   DZ = FZ/DZ
        IF(FZ*FP.GT.0) GO TO 324
        IF(FZ*FP.EQ.0) GO TO 317
        IF(ABS(FZ)-ABS(FP).LE.0) GO TO 353
        ZM = ZZ
        ZZ = (ZM+ZP)*0.5D0
        GO TO 325

  353   ZM = ZP
        GO TO 324

  317   ROOT(I) = ZZ
        ISW = 4
        GO TO 304

  400   E1 = ABS(D1)
        E2 = MAX(ABS(C1*ZZ),ABS(C2))
        DO 606 J = 1,N2
          E3 = MAX(E2*ABS(ALPHA(J)*ZZ+BETA(J)),E1)
          E1 = E2
          E2 = E3
  606   CONTINUE
        E(I) = E2*EPS/ABS(DZ)
  202 CONTINUE
      NR = NRC
  113 RETURN
C  ****** FINDING UPPER AND LOWER BOUNDS TO ROOTS ******
   31 M1 = M + 1
      WMAX = 0.0D0
      W1(1) = 1.0D0
      W2(1) = DBLE(N)
      DO 211 J = 2,M
        I = M1 - J
        W1(J) = A(I)/A(M)
        W2(J) = W1(J)*DBLE(I-1)
        IF(W1(J).EQ.0) GO TO 211
        WLOG = LOG(ABS(W1(J)))/ (DBLE(J)-1.0D0)
        IF(WMAX-WLOG.GE.0) GO TO 211
        WMAX = WLOG
  211 CONTINUE
      X2 = EXP(WMAX)*2.0
      X1 = -X2
      GO TO 602

      END
