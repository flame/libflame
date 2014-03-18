      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      REAL AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      REAL S,ARS,AIS,BRS,BIS
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
      REAL FUNCTION EPSLON (X)
      REAL X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      REAL A,B,C,EPS
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C
C     THIS VERSION DATED 4/6/83.
C
      A = 4.0E0/3.0E0
   10 B = A - 1.0E0
      C = B + B + B
      EPS = ABS(C-1.0E0)
      IF (EPS .EQ. 0.0E0) GO TO 10
      EPSLON = EPS*ABS(X)
      RETURN
      END
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N)
      REAL P,Q,R,S,T,W,X,Y,ZZ,NORM,TST1,TST2
      LOGICAL NOTLAS
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG
C          FORM BY  ELMHES  OR  ORTHES, IF PERFORMED, IS STORED
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED
C          BEFORE CALLING  HQR  IF SUBSEQUENT CALCULATION AND
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C     MODIFIED ON 11/1/89; ADJUSTING INDICES OF LOOPS
C       200, 210, 230, AND 240 TO INCREASE PERFORMANCE. JACK DONGARRA
C
C     ------------------------------------------------------------------
C
*
      EXTERNAL SLAMCH
      REAL SLAMCH, UNFL,OVFL,ULP,SMLNUM,SMALL
      IF (N.LE.0) RETURN
*
*
*     INITIALIZE
      ITCNT = 0
      OPST = 0
      IERR = 0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
*
*        INCREMENT OPCOUNT FOR COMPUTING MATRIX NORM
         OPS = OPS + (IGH-LOW+1)*(IGH-LOW+2)/2
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
      NORM = 0.0E0
      DO 5 J = LOW, IGH
         S = 0.0E0
         DO 4 I = LOW, MIN(IGH,J+1)
              S = S + ABS(H(I,J))
  4      CONTINUE
         NORM = MAX(NORM, S)
  5   CONTINUE
*
      UNFL = SLAMCH( 'SAFE MINIMUM' )
      OVFL = SLAMCH( 'OVERFLOW' )
      ULP = SLAMCH( 'EPSILON' )*SLAMCH( 'BASE' )
      SMLNUM = MAX( UNFL*( N / ULP ), N / ( ULP*OVFL ) )
      SMALL = MAX( SMLNUM, ULP*NORM )
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
*     REPLACE SPLITTING CRITERION WITH NEW ONE AS IN LAPACK
*
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         IF (ABS(H(L,L-1)) .LE. MAX(ULP*S,SMALL))  GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
         OPS = OPS + 2*(EN-L+1)
      X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
*
*        INCREMENT OP COUNT FOR FORMING EXCEPTIONAL SHIFT
         OPS = OPS + (EN-LOW+6)
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
        ITCNT = 30*N - ITN
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
*     REPLACE SPLITTING CRITERION WITH NEW ONE AS IN LAPACK
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF ( TST2 .LE. MAX(ULP*TST1,SMALL) ) GO TO 150
  140 CONTINUE
C
  150 CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 140
         OPST = OPST + 20*(ENM2-M+1)
      MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
*
*        INCREMENT OPCOUNT FOR LOOP 260
         OPST = OPST + 18*(NA-M+1)
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OPCOUNT
         OPS = OPS + 6*(EN-K+1)
         DO 200 J = K, EN
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT
         OPS = OPS + 6*(J-L+1)
         DO 210 I = L, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
         GO TO 255
  225    CONTINUE
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OPCOUNT
         OPS = OPS + 10*(EN-K+1)
         DO 230 J = K, EN
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT
         OPS = OPS + 10*(J-L+1)
         DO 240 I = L, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
  255    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 WR(EN) = X + T
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      X = X + T
*
*        INCREMENT OP COUNT FOR FINDING TWO ROOTS.
         OPST = OPST + 8
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 CONTINUE
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,TST1,TST2
      LOGICAL NOTLAS
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
*
      EXTERNAL SLAMCH
      REAL SLAMCH, UNFL,OVFL,ULP,SMLNUM,SMALL
      IF (N.LE.0) RETURN
*
*     INITIALIZE
*
      ITCNT = 0
      OPST = 0
C
      IERR = 0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
*
*        INCREMENT OPCOUNT FOR COMPUTING MATRIX NORM
         OPS = OPS + (IGH-LOW+1)*(IGH-LOW+2)/2
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
      NORM = 0.0E0
      DO 5 J = LOW, IGH
         S = 0.0E0
         DO 4 I = LOW, MIN(IGH,J+1)
              S = S + ABS(H(I,J))
  4      CONTINUE
         NORM = MAX(NORM, S)
  5   CONTINUE
C
      UNFL = SLAMCH( 'SAFE MINIMUM' )
      OVFL = SLAMCH( 'OVERFLOW' )
      ULP = SLAMCH( 'EPSILON' )*SLAMCH( 'BASE' )
      SMLNUM = MAX( UNFL*( N / ULP ), N / ( ULP*OVFL ) )
      SMALL = MAX( SMLNUM, ULP*NORM )
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
*     REPLACE SPLITTING CRITERION WITH NEW ONE AS IN LAPACK
*
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         IF ( ABS(H(L,L-1)) .LE. MAX(ULP*S,SMALL) )  GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
         OPS = OPS + 2*(EN-L+1)
      X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
*
*        INCREMENT OP COUNT
         OPS = OPS + (EN-LOW+6)
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
        ITCNT = 30*N - ITN
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF ( TST2 .LE. MAX(ULP*TST1,SMALL) ) GO TO 150
  140 CONTINUE
C
  150 CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 140
         OPST = OPST + 20*(ENM2-M+1)
      MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
*
*        INCREMENT OPCOUNT FOR LOOP 260
         OPST = OPST + 18*(NA-M+1)
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OP COUNT FOR LOOP 200
         OPS = OPS + 6*(N-K+1)
         DO 200 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT FOR LOOP 210
         OPS = OPS + 6*J
         DO 210 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
*
*        INCREMENT OPCOUNT FOR LOOP 220
         OPS = OPS + 6*(IGH-LOW + 1)
         DO 220 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
  220    CONTINUE
         GO TO 255
  225    CONTINUE
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OPCOUNT FOR LOOP 230
         OPS = OPS + 10*(N-K+1)
         DO 230 J = K, N
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT FOR LOOP 240
         OPS = OPS + 10*J
         DO 240 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
*
*        INCREMENT OPCOUNT FOR LOOP 250
         OPS = OPS + 10*(IGH-LOW+1)
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1) + ZZ * Z(I,K+2)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K+2) = Z(I,K+2) - P * R
  250    CONTINUE
  255    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
*
*        INCREMENT OP COUNT FOR FINDING TWO ROOTS.
         OPST = OPST + 18
*
*        INCREMENT OP COUNT FOR MODIFICATION AND ACCUMULATION
*        IN LOOP 290, 300, 310
         OPS = OPS + 6*(N-NA+1) + 6*EN + 6*(IGH-LOW+1)
C     .......... ROW MODIFICATION ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     .......... COLUMN MODIFICATION ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
*
*        INCREMENT OP COUNT FOR FINDING COMPLEX PAIR.
         OPST = OPST + 9
  330 EN = ENM2
      GO TO 60
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 IF (NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     .......... REAL VECTOR ..........
  600    M = EN
         H(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = 0.0E0
C
*
*        INCREMENT OP COUNT FOR LOOP 610
         OPST = OPST + 2*(EN - M+1)
            DO 610 J = M, EN
  610       R = R + H(I,J) * H(J,EN)
C
            IF (WI(I) .GE. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 640
            T = W
            IF (T .NE. 0.0E0) GO TO 635
               TST1 = NORM
               T = TST1
  632          T = 0.01E0 * T
               TST2 = NORM + T
               IF (TST2 .GT. TST1) GO TO 632
  635       H(I,EN) = -R / T
            GO TO 680
C     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
*
*        INCREMENT OP COUNT FOR SOLVING REAL EQUATION.
         OPST = OPST + 13
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 680
  650       H(I+1,EN) = (-S - Y * T) / ZZ
C
C     .......... OVERFLOW CONTROL ..........
  680       T = ABS(H(I,EN))
            IF (T .EQ. 0.0E0) GO TO 700
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 700
*
*        INCREMENT OP COUNT.
         OPST = OPST + (EN-I+1)
            DO 690 J = I, EN
               H(J,EN) = H(J,EN)/T
  690       CONTINUE
C
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
*
*        INCREMENT OP COUNT.
         OPST = OPST + 3
         GO TO 730
  720    CALL CDIV(0.0E0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
*
*        INCREMENT OP COUNT IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN)))
         OPST = OPST + 16
  730    H(EN,NA) = 0.0E0
         H(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 795 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0E0
            SA = 0.0E0
C
*
*        INCREMENT OP COUNT FOR LOOP 760
         OPST = OPST + 4*(EN-M+1)
            DO 760 J = M, EN
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0E0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 795
  770       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
*
*        INCREMENT OP COUNT FOR CDIV
         OPST = OPST + 16
            GO TO 790
C     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0E0 * Q
*
*        INCREMENT OPCOUNT (AVERAGE) FOR SOLVING COMPLEX EQUATIONS
         OPST = OPST + 42
            IF (VR .NE. 0.0E0 .OR. VI .NE. 0.0E0) GO TO 784
               TST1 = NORM * (ABS(W) + ABS(Q) + ABS(X)
     X                      + ABS(Y) + ABS(ZZ))
               VR = TST1
  783          VR = 0.01E0 * VR
               TST2 = TST1 + VR
               IF (TST2 .GT. TST1) GO TO 783
  784       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,
     X                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,
     X                H(I+1,NA),H(I+1,EN))
C
C     .......... OVERFLOW CONTROL ..........
  790       T = AMAX1(ABS(H(I,NA)), ABS(H(I,EN)))
            IF (T .EQ. 0.0E0) GO TO 795
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 795
*
*        INCREMENT OP COUNT.
         OPST = OPST + 2*(EN-I+1)
            DO 792 J = I, EN
               H(J,NA) = H(J,NA)/T
               H(J,EN) = H(J,EN)/T
  792       CONTINUE
C
  795    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
*
*        INCREMENT OP COUNT.
         OPS = OPS + 2*(IGH-LOW+1)*(M-LOW+1)
         DO 880 I = LOW, IGH
            ZZ = 0.0E0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 CONTINUE
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE IMTQL1(N,D,E,IERR)
*
*     EISPACK ROUTINE
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEQR.
*
C
      INTEGER I,J,L,M,N,II,MML,IERR
      REAL D(N),E(N)
      REAL B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      REAL             EPS, TST
      REAL             SLAMCH
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE CONTRIBUTIONS TO OPS FROM
*     FUNCTION PYTHAG.  IT IS PASSED TO AND FROM PYTHAG
*     THROUGH COMMON BLOCK PYTHOP.
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
      COMMON             / PYTHOP / OPST
*
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL1,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 40 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
*
*        INITIALIZE ITERATION COUNT AND OPST
            ITCNT = 0
            OPST = 0
*
*     DETERMINE THE UNIT ROUNDOFF FOR THIS ENVIRONMENT.
*
      EPS = SLAMCH( 'EPSILON' )
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST = ABS( E(M) )
            IF( TST .LE. EPS * ( ABS(D(M)) + ABS(D(M+1)) ) ) GO TO 120
*            TST1 = ABS(D(M)) + ABS(D(M+1))
*            TST2 = TST1 + ABS(E(M))
*            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
*
*        INCREMENT OPCOUNT FOR FINDING SMALL SUBDIAGONAL ELEMENT.
            OPS = OPS + 2*( MIN(M,N-1)-L+1 )
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 40) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0E0 * E(L))
         R = PYTHAG(G,1.0E0)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
            OPS = OPS + 7
         S = 1.0E0
         C = 1.0E0
         P = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. 0.0E0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0E0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
            OPS = OPS + MML*14 + 1
*
*        INCREMENT ITERATION COUNTER
            ITCNT = ITCNT + 1
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP, WHEN UNDERFLOW OCCURS.
            OPS = OPS + 2+(II-1)*14 + 1
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 40 ITERATIONS ..........
 1000 IERR = L
 1001 CONTINUE
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
*
*     EISPACK ROUTINE.  MODIFIED FOR COMPARISON WITH LAPACK.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEQR.
*
C
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      REAL             EPS, TST
      REAL             SLAMCH
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE CONTRIBUTIONS TO OPS FROM
*     FUNCTION PYTHAG.  IT IS PASSED TO AND FROM PYTHAG
*     THROUGH COMMON BLOCK PYTHOP.
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
      COMMON             / PYTHOP / OPST
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 40 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
*
*        INITIALIZE ITERATION COUNT AND OPST
            ITCNT = 0
            OPST = 0
*
*     DETERMINE UNIT ROUNDOFF FOR THIS MACHINE.
      EPS = SLAMCH( 'EPSILON' )
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0E0
C
      DO 240 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
*            TST1 = ABS(D(M)) + ABS(D(M+1))
*            TST2 = TST1 + ABS(E(M))
*            IF (TST2 .EQ. TST1) GO TO 120
            TST = ABS( E(M) )
            IF( TST .LE. EPS * ( ABS(D(M)) + ABS(D(M+1)) ) ) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
*
*        INCREMENT OPCOUNT FOR FINDING SMALL SUBDIAGONAL ELEMENT.
            OPS = OPS + 2*( MIN(M,N)-L+1 )
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 40) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0E0 * E(L))
         R = PYTHAG(G,1.0E0)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
            OPS = OPS + 7
         S = 1.0E0
         C = 1.0E0
         P = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. 0.0E0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0E0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
C
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
            OPS = OPS + MML*( 14+6*N ) + 1
*
*        INCREMENT ITERATION COUNTER
            ITCNT = ITCNT + 1
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP, WHEN UNDERFLOW OCCURS.
            OPS = OPS + 2+(II-1)*(14+6*N) + 1
         GO TO 105
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 40 ITERATIONS ..........
 1000 IERR = L
 1001 CONTINUE
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE INVIT(NM,N,A,WR,WI,SELECT,MM,M,Z,IERR,RM1,RV1,RV2)
C
      INTEGER I,J,K,L,M,N,S,II,IP,MM,MP,NM,NS,N1,UK,IP1,ITS,KM1,IERR
      REAL A(NM,N),WR(N),WI(N),Z(NM,MM),RM1(N,N),
     X       RV1(N),RV2(N)
      REAL T,W,X,Y,EPS3,NORM,NORMV,GROWTO,ILAMBD,
     X       PYTHAG,RLAMBD,UKROOT
      LOGICAL SELECT(N)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE INVIT
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A REAL UPPER
C     HESSENBERG MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE HESSENBERG MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVALUES OF THE MATRIX.  THE EIGENVALUES MUST BE
C          STORED IN A MANNER IDENTICAL TO THAT OF SUBROUTINE  HQR,
C          WHICH RECOGNIZES POSSIBLE SPLITTING OF THE MATRIX.
C
C        SELECT SPECIFIES THE EIGENVECTORS TO BE FOUND. THE
C          EIGENVECTOR CORRESPONDING TO THE J-TH EIGENVALUE IS
C          SPECIFIED BY SETTING SELECT(J) TO .TRUE..
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          COLUMNS REQUIRED TO STORE THE EIGENVECTORS TO BE FOUND.
C          NOTE THAT TWO COLUMNS ARE REQUIRED TO STORE THE
C          EIGENVECTOR CORRESPONDING TO A COMPLEX EIGENVALUE.
C
C     ON OUTPUT
C
C        A AND WI ARE UNALTERED.
C
C        WR MAY HAVE BEEN ALTERED SINCE CLOSE EIGENVALUES ARE PERTURBED
C          SLIGHTLY IN SEARCHING FOR INDEPENDENT EIGENVECTORS.
C
C        SELECT MAY HAVE BEEN ALTERED.  IF THE ELEMENTS CORRESPONDING
C          TO A PAIR OF CONJUGATE COMPLEX EIGENVALUES WERE EACH
C          INITIALLY SET TO .TRUE., THE PROGRAM RESETS THE SECOND OF
C          THE TWO ELEMENTS TO .FALSE..
C
C        M IS THE NUMBER OF COLUMNS ACTUALLY USED TO STORE
C          THE EIGENVECTORS.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE NEXT SELECTED EIGENVALUE IS REAL, THE NEXT COLUMN
C          OF Z CONTAINS ITS EIGENVECTOR.  IF THE EIGENVALUE IS
C          COMPLEX, THE NEXT TWO COLUMNS OF Z CONTAIN THE REAL AND
C          IMAGINARY PARTS OF ITS EIGENVECTOR.  THE EIGENVECTORS ARE
C          NORMALIZED SO THAT THE COMPONENT OF LARGEST MAGNITUDE IS 1.
C          ANY VECTOR WHICH FAILS THE ACCEPTANCE TEST IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -(2*N+1)   IF MORE THAN MM COLUMNS OF Z ARE NECESSARY
C                     TO STORE THE EIGENVECTORS CORRESPONDING TO
C                     THE SPECIFIED EIGENVALUES.
C          -K         IF THE ITERATION CORRESPONDING TO THE K-TH
C                     VALUE FAILS,
C          -(N+K)     IF BOTH ERROR SITUATIONS OCCUR.
C
C        RM1, RV1, AND RV2 ARE TEMPORARY STORAGE ARRAYS.  NOTE THAT RM1
C          IS SQUARE OF DIMENSION N BY N AND, AUGMENTED BY TWO COLUMNS
C          OF Z, IS THE TRANSPOSE OF THE CORRESPONDING ALGOL B ARRAY.
C
C     THE ALGOL PROCEDURE GUESSVEC APPEARS IN INVIT IN LINE.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
*
*     GET ULP FROM SLAMCH FOR NEW SMALL PERTURBATION AS IN LAPACK
      EXTERNAL SLAMCH
      REAL SLAMCH, ULP
      IF (N.LE.0) RETURN
      ULP = SLAMCH( 'EPSILON' )
C
*
*     INITIALIZE
      OPST = 0
      IERR = 0
      UK = 0
      S = 1
C     .......... IP = 0, REAL EIGENVALUE
C                     1, FIRST OF CONJUGATE COMPLEX PAIR
C                    -1, SECOND OF CONJUGATE COMPLEX PAIR ..........
      IP = 0
      N1 = N - 1
C
      DO 980 K = 1, N
         IF (WI(K) .EQ. 0.0E0 .OR. IP .LT. 0) GO TO 100
         IP = 1
         IF (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
  100    IF (.NOT. SELECT(K)) GO TO 960
         IF (WI(K) .NE. 0.0E0) S = S + 1
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (A(UK+1,UK) .EQ. 0.0E0) GO TO 140
  120    CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
  140    NORM = 0.0E0
         MP = 1
C
*
*        INCREMENT OPCOUNT FOR COMPUTING MATRIX NORM
         OPS = OPS + UK*(UK-1)/2
         DO 180 I = 1, UK
            X = 0.0E0
C
            DO 160 J = MP, UK
  160       X = X + ABS(A(I,J))
C
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
         IF (NORM .EQ. 0.0E0) NORM = 1.0E0
*        EPS3 = EPSLON(NORM)
*
*        INCREMENT OPCOUNT
         OPST = OPST + 3
         EPS3 = NORM*ULP
C     .......... GROWTO IS THE CRITERION FOR THE GROWTH ..........
         UKROOT = UK
         UKROOT = SQRT(UKROOT)
         GROWTO = 0.1E0 / UKROOT
  200    RLAMBD = WR(K)
         ILAMBD = WI(K)
         IF (K .EQ. 1) GO TO 280
         KM1 = K - 1
         GO TO 240
C     .......... PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE ..........
  220    RLAMBD = RLAMBD + EPS3
C     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
  240    DO 260 II = 1, KM1
            I = K - II
            IF (SELECT(I) .AND. ABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     X         ABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 260 (ASSUME THAT ALL EIGENVALUES
*        ARE DIFFERENT)
         OPST = OPST + 2*(K-1)
C
         WR(K) = RLAMBD
C     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH ..........
         IP1 = K + IP
         WR(IP1) = RLAMBD
C     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
C                AND INITIAL REAL VECTOR ..........
  280    MP = 1
C
*
*        INCREMENT OP COUNT FOR LOOP 320
         OPS = OPS + UK
         DO 320 I = 1, UK
C
            DO 300 J = MP, UK
  300       RM1(J,I) = A(I,J)
C
            RM1(I,I) = RM1(I,I) - RLAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C
         ITS = 0
         IF (ILAMBD .NE. 0.0E0) GO TO 520
C     .......... REAL EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
         IF (UK .EQ. 1) GO TO 420
C
*
*        INCREMENT OPCOUNT LU DECOMPOSITION
         OPS = OPS + (UK-1)*(UK+2)
         DO 400 I = 2, UK
            MP = I - 1
            IF (ABS(RM1(MP,I)) .LE. ABS(RM1(MP,MP))) GO TO 360
C
            DO 340 J = MP, UK
               Y = RM1(J,I)
               RM1(J,I) = RM1(J,MP)
               RM1(J,MP) = Y
  340       CONTINUE
C
  360       IF (RM1(MP,MP) .EQ. 0.0E0) RM1(MP,MP) = EPS3
            X = RM1(MP,I) / RM1(MP,MP)
            IF (X .EQ. 0.0E0) GO TO 400
C
            DO 380 J = I, UK
  380       RM1(J,I) = RM1(J,I) - X * RM1(J,MP)
C
  400    CONTINUE
C
  420    IF (RM1(UK,UK) .EQ. 0.0E0) RM1(UK,UK) = EPS3
C     .......... BACK SUBSTITUTION FOR REAL VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  440    DO 500 II = 1, UK
            I = UK + 1 - II
            Y = RV1(I)
            IF (I .EQ. UK) GO TO 480
            IP1 = I + 1
C
            DO 460 J = IP1, UK
  460       Y = Y - RM1(J,I) * RV1(J)
C
  480       RV1(I) = Y / RM1(I,I)
  500    CONTINUE
*
*        INCREMENT OP COUNT FOR BACK SUBSTITUTION LOOP 500
         OPS = OPS + UK*(UK+1)
C
         GO TO 740
C     .......... COMPLEX EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
C                PARTS IN UPPER TRIANGLE STARTING AT (1,3) ..........
  520    NS = N - S
         Z(1,S-1) = -ILAMBD
         Z(1,S) = 0.0E0
         IF (N .EQ. 2) GO TO 550
         RM1(1,3) = -ILAMBD
         Z(1,S-1) = 0.0E0
         IF (N .EQ. 3) GO TO 550
C
         DO 540 I = 4, N
  540    RM1(1,I) = 0.0E0
C
  550    DO 640 I = 2, UK
            MP = I - 1
            W = RM1(MP,I)
            IF (I .LT. N) T = RM1(MP,I+1)
            IF (I .EQ. N) T = Z(MP,S-1)
            X = RM1(MP,MP) * RM1(MP,MP) + T * T
            IF (W * W .LE. X) GO TO 580
            X = RM1(MP,MP) / W
            Y = T / W
            RM1(MP,MP) = W
            IF (I .LT. N) RM1(MP,I+1) = 0.0E0
            IF (I .EQ. N) Z(MP,S-1) = 0.0E0
C
*
*        INCREMENT OPCOUNT FOR LOOP 560
         OPS = OPS + 4*(UK-I+1)
            DO 560 J = I, UK
               W = RM1(J,I)
               RM1(J,I) = RM1(J,MP) - X * W
               RM1(J,MP) = W
               IF (J .LT. N1) GO TO 555
               L = J - NS
               Z(I,L) = Z(MP,L) - Y * W
               Z(MP,L) = 0.0E0
               GO TO 560
  555          RM1(I,J+2) = RM1(MP,J+2) - Y * W
               RM1(MP,J+2) = 0.0E0
  560       CONTINUE
C
            RM1(I,I) = RM1(I,I) - Y * ILAMBD
            IF (I .LT. N1) GO TO 570
            L = I - NS
            Z(MP,L) = -ILAMBD
            Z(I,L) = Z(I,L) + X * ILAMBD
            GO TO 640
  570       RM1(MP,I+2) = -ILAMBD
            RM1(I,I+2) = RM1(I,I+2) + X * ILAMBD
            GO TO 640
  580       IF (X .NE. 0.0E0) GO TO 600
            RM1(MP,MP) = EPS3
            IF (I .LT. N) RM1(MP,I+1) = 0.0E0
            IF (I .EQ. N) Z(MP,S-1) = 0.0E0
            T = 0.0E0
            X = EPS3 * EPS3
  600       W = W / X
            X = RM1(MP,MP) * W
            Y = -T * W
C
*
*        INCREMENT OPCOUNT FOR LOOP 620
         OPS = OPS + 6*(UK-I+1)
            DO 620 J = I, UK
               IF (J .LT. N1) GO TO 610
               L = J - NS
               T = Z(MP,L)
               Z(I,L) = -X * T - Y * RM1(J,MP)
               GO TO 615
  610          T = RM1(MP,J+2)
               RM1(I,J+2) = -X * T - Y * RM1(J,MP)
  615          RM1(J,I) = RM1(J,I) - X * RM1(J,MP) + Y * T
  620       CONTINUE
C
            IF (I .LT. N1) GO TO 630
            L = I - NS
            Z(I,L) = Z(I,L) - ILAMBD
            GO TO 640
  630       RM1(I,I+2) = RM1(I,I+2) - ILAMBD
  640    CONTINUE
*
*        INCREMENT OP COUNT (AVERAGE) FOR COMPUTING
*        THE SCALARS IN LOOP 640
         OPS = OPS + 10*(UK -1)
C
         IF (UK .LT. N1) GO TO 650
         L = UK - NS
         T = Z(UK,L)
         GO TO 655
  650    T = RM1(UK,UK+2)
  655    IF (RM1(UK,UK) .EQ. 0.0E0 .AND. T .EQ. 0.0E0) RM1(UK,UK) = EPS3
C     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0E0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C
            DO 680 J = IP1, UK
               IF (J .LT. N1) GO TO 670
               L = J - NS
               T = Z(I,L)
               GO TO 675
  670          T = RM1(I,J+2)
  675          X = X - RM1(J,I) * RV1(J) + T * RV2(J)
               Y = Y - RM1(J,I) * RV2(J) - T * RV1(J)
  680       CONTINUE
C
  700       IF (I .LT. N1) GO TO 710
            L = I - NS
            T = Z(I,L)
            GO TO 715
  710       T = RM1(I,I+2)
  715       CALL CDIV(X,Y,RM1(I,I),T,RV1(I),RV2(I))
  720    CONTINUE
*
*        INCREMENT OP COUNT FOR LOOP 720.
         OPS = OPS + 4*UK*(UK+3)
C     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX
C                EIGENVECTOR AND NORMALIZATION ..........
  740    ITS = ITS + 1
         NORM = 0.0E0
         NORMV = 0.0E0
C
         DO 780 I = 1, UK
            IF (ILAMBD .EQ. 0.0E0) X = ABS(RV1(I))
            IF (ILAMBD .NE. 0.0E0) X = PYTHAG(RV1(I),RV2(I))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
*
*        INCREMENT OP COUNT ACCEPTANCE TEST
         IF (ILAMBD .EQ. 0.0E0) OPS = OPS + UK
         IF (ILAMBD .NE. 0.0E0) OPS = OPS + 16*UK
C
         IF (NORM .LT. GROWTO) GO TO 840
C     .......... ACCEPT VECTOR ..........
         X = RV1(J)
         IF (ILAMBD .EQ. 0.0E0) X = 1.0E0 / X
         IF (ILAMBD .NE. 0.0E0) Y = RV2(J)
C
*
*        INCREMENT OPCOUNT FOR LOOP 820
         IF (ILAMBD .EQ. 0.0E0) OPS = OPS + UK
         IF (ILAMBD .NE. 0.0E0) OPS = OPS + 16*UK
         DO 820 I = 1, UK
            IF (ILAMBD .NE. 0.0E0) GO TO 800
            Z(I,S) = RV1(I) * X
            GO TO 820
  800       CALL CDIV(RV1(I),RV2(I),X,Y,Z(I,S-1),Z(I,S))
  820    CONTINUE
C
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0E0)
         RV1(1) = EPS3
C
         DO 860 I = 2, UK
  860    RV1(I) = Y
C
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         IF (ILAMBD .EQ. 0.0E0) GO TO 440
         GO TO 660
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
         IERR = -K
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
            Z(I,S) = 0.0E0
            IF (ILAMBD .NE. 0.0E0) Z(I,S-1) = 0.0E0
  920    CONTINUE
C
  940    S = S + 1
  960    IF (IP .EQ. (-1)) IP = 0
         IF (IP .EQ. 1) IP = -1
  980 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1 - IABS(IP)
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL A(NM,N),ORT(IGH)
      REAL F,G,H,SCALE
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N.LE.0) RETURN
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
*
*     INCREMENT OP COUNR FOR COMPUTING G,H,ORT(M),.. IN LOOP 180
      OPS = OPS + 6*(LA - KP1 + 1)
      DO 180 M = KP1, LA
         H = 0.0E0
         ORT(M) = 0.0E0
         SCALE = 0.0E0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
*
*     INCREMENT OP COUNT FOR LOOP 90
      OPS = OPS + (IGH-M +1)
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(A(I,M-1))
C
         IF (SCALE .EQ. 0.0E0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
*
*     INCREMENT OP COUNT FOR LOOP 100
      OPS = OPS + 3*(IGH-M+1)
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
C
         G = -SIGN(SQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
C     .......... FORM (I-(U*UT)/H) * A ..........
*
*     INCREMENT OP COUNT FOR LOOP 130 AND 160
      OPS = OPS + (N-M+1+IGH)*(4*(IGH-M+1) + 1)
         DO 130 J = M, N
            F = 0.0E0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
C
            F = F / H
C
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            F = 0.0E0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
C
            F = F / H
C
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
C
  160    CONTINUE
C
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
C
  200 RETURN
      END
      REAL FUNCTION PYTHAG(A,B)
      REAL A,B
C
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT
*     OPST IS ONLY INCREMENTED HERE
*     .. COMMON BLOCKS ..
      COMMON             / PYTHOP / OPST
*     ..
*     .. SCALARS IN COMMON
      REAL               OPST
*     ..
      REAL P,R,S,T,U
      P = AMAX1(ABS(A),ABS(B))
      IF (P .EQ. 0.0E0) GO TO 20
      R = (AMIN1(ABS(A),ABS(B))/P)**2
*
*     INCREMENT OPST
      OPST = OPST + 2
   10 CONTINUE
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         U = 1.0E0 + 2.0E0*S
         P = U*P
         R = (S/U)**2 * R
*
*        INCREMENT OPST
            OPST = OPST + 8
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      SUBROUTINE TQLRAT(N,D,E2,IERR)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEQR.
*
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL D(N),E2(N)
      REAL B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
      REAL             EPS, TST
      REAL             SLAMCH
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE CONTRIBUTIONS TO OPS FROM
*     FUNCTION PYTHAG.  IT IS PASSED TO AND FROM PYTHAG
*     THROUGH COMMON BLOCK PYTHOP.
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
      COMMON             / PYTHOP / OPST
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E2 HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
*
*        INITIALIZE ITERATION COUNT AND OPST
            ITCNT = 0
            OPST = 0
*
*     DETERMINE THE UNIT ROUNDOFF FOR THIS ENVIRONMENT.
*
      EPS = SLAMCH( 'EPSILON' )
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0E0
      T = 0.0E0
      E2(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
         H = ABS(D(L)) + SQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
*
*     INCREMENT OPCOUNT FOR THIS SECTION.
*     (FUNCTION EPSLON IS COUNTED AS 6 FLOPS.  THIS IS THE MINIMUM
*     NUMBER REQUIRED, BUT COUNTING THEM EXACTLY WOULD AFFECT
*     THE TIMING.)
         OPS = OPS + 9
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF( M .EQ. N ) GO TO 120
            TST = SQRT( ABS( E2(M) ) )
            IF( TST .LE. EPS * ( ABS(D(M)) + ABS(D(M+1)) ) ) GO TO 120
*            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    CONTINUE
*
*        INCREMENT OPCOUNT FOR FINDING SMALL SUBDIAGONAL ELEMENT.
            OPS = OPS + 3*( MIN(M,N-1)-L+1 )
         IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * S)
         R = PYTHAG(P,1.0E0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT AND SUBTRACTING.
            OPS = OPS + 8 + (I-L1+1)
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0E0) G = B
         H = G
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0E0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
            OPS = OPS + MML*11 + 1
*
*        INCREMENT ITERATION COUNTER
            ITCNT = ITCNT + 1
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0E0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0E0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 CONTINUE
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),E2(N)
      REAL F,G,H,SCALE
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT.
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED.
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*
      OPS = OPS + MAX( 0.0E0, (4.0E0/3.0E0)*REAL(N)**3 +
     $                              12.0E0*REAL(N)**2 +
     $                      (11.0E0/3.0E0)*N - 22 )
*
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(D(K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
C
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0E0
  125    CONTINUE
C
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 300
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0E0
C
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0E0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         H = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
C
  280    CONTINUE
C
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEBZ.
*
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM
      REAL D(N),E(N),E2(N),W(MM),RV4(N),RV5(N)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(MM)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE BISECTION TECHNIQUE
C     IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL,
C     USING BISECTION.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES.
C          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND.
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVALUES IN THE INTERVAL.  WARNING. IF MORE THAN
C          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL,
C          AN ERROR RETURN IS MADE WITH NO EIGENVALUES FOUND.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB).
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF M EXCEEDS MM.
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM
C     APPEARS IN BISECT IN-LINE.
C
C     NOTE THAT SUBROUTINE TQL1 OR IMTQL1 IS GENERALLY FASTER THAN
C     BISECT, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      REAL             ONE
      PARAMETER        ( ONE = 1.0E0 )
      REAL             RELFAC
      PARAMETER        ( RELFAC = 2.0E0 )
      REAL ATOLI, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP
      REAL SLAMCH, PIVMIN
      EXTERNAL SLAMCH
*        INITIALIZE ITERATION COUNT.
            ITCNT = 0
      SAFEMN = SLAMCH( 'S' )
      ULP = SLAMCH( 'E' )*SLAMCH( 'B' )
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
CCC         TST1 = ABS(D(I)) + ABS(D(I-1))
CCC         TST2 = TST1 + ABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
         TMP1 = E( I )**2
         IF( ABS( D(I)*D(I-1) )*ULP**2+SAFEMN.LE.TMP1 )
     $      GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
*           INCREMENT OPCOUNT FOR DETERMINING IF MATRIX SPLITS.
               OPS = OPS + 5*( N-1 )
C
C                COMPUTE QUANTITIES NEEDED FOR CONVERGENCE TEST.
      TMP1 = D( 1 ) - ABS( E( 2 ) )
      TMP2 = D( 1 ) + ABS( E( 2 ) )
      PIVMIN = ONE
      DO 41 I = 2, N - 1
         TMP1 = MIN( TMP1, D( I )-ABS( E( I ) )-ABS( E( I+1 ) ) )
         TMP2 = MAX( TMP2, D( I )+ABS( E( I ) )+ABS( E( I+1 ) ) )
         PIVMIN = MAX( PIVMIN, E( I )**2 )
   41 CONTINUE
      TMP1 = MIN( TMP1, D( N )-ABS( E( N ) ) )
      TMP2 = MAX( TMP2, D( N )+ABS( E( N ) ) )
      PIVMIN = MAX( PIVMIN, E( N )**2 )
      PIVMIN = PIVMIN*SAFEMN
      TNORM = MAX( ABS(TMP1), ABS(TMP2) )
      ATOLI = ULP*TNORM
*        INCREMENT OPCOUNT FOR COMPUTING THESE QUANTITIES.
            OPS = OPS + 4*( N-1 )
C
C     .......... DETERMINE THE NUMBER OF EIGENVALUES
C                IN THE INTERVAL ..........
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 320
   60 M = S
      X1 = LB
      ISTURM = 2
      GO TO 320
   80 M = M - S
      IF (M .GT. MM) GO TO 980
      Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = AMIN1(D(Q)-(X1+U),XU)
         X0 = AMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
            OPS = OPS + ( N-P+1 )*2
C
  140 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
CCC         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
CCC         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
         TMP1 = ABS( X0 - XU )
         TMP2 = MAX( ABS( X0 ), ABS( XU ) )
         IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) )
     $      GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / EPSLON(1.0E0)
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
*           INCREMENT OPCOUNT FOR STURM SEQUENCE.
               OPS = OPS + ( Q-P+1 )*3
*           INCREMENT ITERATION COUNTER.
               ITCNT = ITCNT + 1
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
C                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
*
*     EISPACK ROUTINE.
*
*     CONVERGENCE TEST WAS NOT MODIFIED, SINCE IT SHOULD GIVE
*     APPROXIMATELY THE SAME LEVEL OF ACCURACY AS LAPACK ROUTINE,
*     ALTHOUGH THE EIGENVECTORS MAY NOT BE AS CLOSE TO ORTHOGONAL.
*
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      REAL D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,EPSLON,
     X       PYTHAG
      INTEGER IND(M)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
      COMMON             / PYTHOP / OPST
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0E0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0E0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES.
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT
C
C        ALL INPUT ARRAYS ARE UNALTERED.
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*        INITIALIZE ITERATION COUNT.
            ITCNT = 0
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0E0 - E2(1)
      Q = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     .......... CHECK FOR ISOLATED ROOT ..........
         XU = 1.0E0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0E0
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = AMAX1(NORM, ABS(D(I))+ABS(E(I)))
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
         EPS2 = 1.0E-3 * NORM
         EPS3 = EPSLON(NORM)
         UK = Q - P + 1
         EPS4 = UK * EPS3
         UK = EPS4 / SQRT(UK)
*           INCREMENT OPCOUNT FOR COMPUTING CRITERIA.
               OPS = OPS + ( Q-IP+4 )
         S = P
  505    GROUP = 0
         GO TO 520
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0E0) X1 = X0 + ORDER * EPS3
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
C     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0E0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
*           INCREMENT OPCOUNT FOR ELIMINATION.
               OPS = OPS + ( Q-P+1 )*5
C
         IF (U .EQ. 0.0E0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0E0
         RV3(Q) = 0.0E0
C     .......... BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
*           INCREMENT OPCOUNT FOR BACK SUBSTITUTION.
               OPS = OPS + ( Q-P+1 )*5
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0E0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
*              INCREMENT OPCOUNT FOR ORTHOGONALIZING.
                  OPS = OPS + ( Q-P+1 )*4
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
*           INCREMENT OPCOUNT FOR COMPUTING NORM.
               OPS = OPS + ( Q-P+1 )
C
         IF (NORM .GE. 1.0E0) GO TO 840
C     .......... FORWARD SUBSTITUTION ..........
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0E0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ..........
  780    DO 820 I = IP, Q
            U = RV6(I)
C     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ..........
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
*           INCREMENT OPCOUNT FOR FORWARD SUBSTITUTION.
               OPS = OPS + ( Q-P+1 ) + ( Q-IP+1 )*2
C
         ITS = ITS + 1
*           INCREMENT ITERATION COUNTER.
               ITCNT = ITCNT + 1
         GO TO 600
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
         XU = 0.0E0
         GO TO 870
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = P, Q
  860    U = PYTHAG(U,RV6(I))
C
         XU = 1.0E0 / U
C
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
*           INCREMENT OPCOUNT FOR NORMALIZING.
               OPS = OPS + ( Q-P+1 )
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
*        INCREMENT OPCOUNT FOR USE OF FUNCTION PYTHAG.
            OPS = OPS + OPST
 1001 RETURN
      END
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEBZ.
*
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      REAL D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(M)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES.
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES.
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE.
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      REAL             ONE
      PARAMETER        ( ONE = 1.0E0 )
      REAL             RELFAC
      PARAMETER        ( RELFAC = 2.0E0 )
      REAL ATOLI, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP
      REAL SLAMCH, PIVMIN
      EXTERNAL SLAMCH
*        INITIALIZE ITERATION COUNT.
            ITCNT = 0
      SAFEMN = SLAMCH( 'S' )
      ULP = SLAMCH( 'E' )*SLAMCH( 'B' )
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
      PIVMIN = ONE
      DO 40 I = 1, N
         X1 = U
         U = 0.0E0
         IF (I .NE. N) U = ABS(E(I+1))
         XU = AMIN1(D(I)-(X1+U),XU)
         X0 = AMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
CCC         TST1 = ABS(D(I)) + ABS(D(I-1))
CCC         TST2 = TST1 + ABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
         TMP1 = E( I )**2
         IF( ABS( D(I)*D(I-1) )*ULP**2+SAFEMN.LE.TMP1 ) THEN
            PIVMIN = MAX( PIVMIN, TMP1 )
            GO TO 40
         END IF
   20    E2(I) = 0.0E0
   40 CONTINUE
      PIVMIN = PIVMIN*SAFEMN
      TNORM = MAX( ABS( XU ), ABS( X0 ) )
      ATOLI = ULP*TNORM
*        INCREMENT OPCOUNT FOR DETERMINING IF MATRIX SPLITS.
            OPS = OPS + 9*( N-1 )
C
      X1 = N
      X1 = X1 * EPSLON(AMAX1(ABS(XU),ABS(X0)))
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     .......... DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES ..........
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5E0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = AMIN1(D(Q)-(X1+U),XU)
         X0 = AMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
            OPS = OPS + ( N-P+1 )*2
C
  140 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
CCC         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
CCC         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
         TMP1 = ABS( X0 - XU )
         TMP2 = MAX( ABS( X0 ), ABS( XU ) )
         IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) )
     $      GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / EPSLON(1.0E0)
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
*           INCREMENT OPCOUNT FOR STURM SEQUENCE.
               OPS = OPS + ( Q-P+1 )*3
*           INCREMENT ITERATION COUNTER.
               ITCNT = ITCNT + 1
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES ..........
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE SSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      REAL X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, IOPS IS ONLY INCREMENTED
*     IOPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO IOPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ IOPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL IOPS, ITCNT, IOPST
*     ..
C
C
C     SSVDC IS A SUBROUTINE TO REDUCE A REAL NXP MATRIX X BY
C     ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         REAL(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY SSVDC.
C
C         LDX       INTEGER.
C                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C         N         INTEGER.
C                   N IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C         P         INTEGER.
C                   P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C         LDU       INTEGER.
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
C                   (SEE BELOW).
C
C         WORK      REAL(N).
C                   WORK IS A SCRATCH ARRAY.
C
C         JOB       INTEGER.
C                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
C                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
C                   WITH THE FOLLOWING MEANING
C
C                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
C                                  VECTORS.
C                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
C                                  IN U.
C                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
C                                  VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         REAL(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         REAL(P).
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         REAL(LDU,K), WHERE LDU.GE.N.  IF JOBA.EQ.1 THEN
C                                   K.EQ.N, IF JOBA.GE.2 THEN
C                                   K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF LEFT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         REAL(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WITH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
C                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
C                   VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C              CORRECTION TO SHIFT CALCULATION MADE 2/85.
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     ***** USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL SROT
C     BLAS SAXPY,SDOT,SSCAL,SSWAP,SNRM2,SROTG
C     FORTRAN ABS,AMAX1,MAX0,MIN0,MOD,SQRT
C
C     INTERNAL VARIABLES
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     *        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      REAL SDOT,T
      REAL B,C,CS,EL,EMM1,F,G,SNRM2,SCALE,SHIFT,SL,SM,SN,SMM1,T1,TEST
*     REAL ZTEST,R
      LOGICAL WANTU,WANTV
*
*     GET EPS FROM SLAMCH FOR NEW STOPPING CRITERION
      EXTERNAL SLAMCH
      REAL SLAMCH, EPS
      IF (N.LE.0 .OR. P.LE.0) RETURN
      EPS = SLAMCH( 'EPSILON' )
*
C
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
      MAXIT = 50
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
*
*     INITIALIZE OP COUNT
      IOPST = 0
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
*
*           INCREMENT OP COUNT
            IOPS = IOPS + (2*(N-L+1)+1)
            S(L) = SNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0E0) GO TO 10
               IF (X(L,L) .NE. 0.0E0) S(L) = SIGN(S(L),X(L,L))
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (N-L+3)
               CALL SSCAL(N-L+1,1.0E0/S(L),X(L,L),1)
               X(L,L) = 1.0E0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0E0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (4*(N-L)+5)
               T = -SDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL SAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
*
*           INCREMENT OP COUNT
            IOPS = IOPS + (2*(P-L)+1)
            E(L) = SNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0E0) GO TO 80
               IF (E(LP1) .NE. 0.0E0) E(L) = SIGN(E(L),E(LP1))
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (P-L+2)
               CALL SSCAL(P-L,1.0E0/E(L),E(LP1),1)
               E(LP1) = 1.0E0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0E0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = 0.0E0
   90          CONTINUE
*
*              INCREMENT OP COUNT
               IOPS = IOPS + FLOAT(4*(N-L)+1)*(P-L)
               DO 100 J = LP1, P
                  CALL SAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL SAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0E0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0E0
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0E0
  180       CONTINUE
            U(J,J) = 1.0E0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0E0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (FLOAT(4*(N-L)+5)*(NCU-L)+(N-L+2))
               DO 210 J = LP1, NCU
                  T = -SDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL SAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL SSCAL(N-L+1,-1.0E0,U(L,L),1)
               U(L,L) = 1.0E0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0E0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0E0
  260          CONTINUE
               U(L,L) = 1.0E0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0E0) GO TO 320
*
*              INCREMENT OP COUNT
               IOPS = IOPS + FLOAT(4*(P-L)+1)*(P-L)
               DO 310 J = LP1, P
                  T = -SDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL SAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0E0
  330       CONTINUE
            V(L,L) = 1.0E0
  340    CONTINUE
  350 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
*
*     INITIALIZE ITERATION COUNTER
      ITCNT = 0
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
*
*        UPDATE ITERATION COUNTER
         ITCNT = ITER
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
C     ......EXIT
            GO TO 620
  370    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 390 LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) GO TO 400
*
*           INCREMENT OP COUNT
            IOPST = IOPST + 2
            TEST = ABS(S(L)) + ABS(S(L+1))
*
*           REPLACE STOPPING CRITERION WITH NEW ONE AS IN LAPACK
*
*           ZTEST = TEST + ABS(E(L))
*           IF (ZTEST .NE. TEST) GO TO 380
            IF (ABS(E(L)) .GT. EPS * TEST) GOTO 380
*
               E(L) = 0.0E0
C        ......EXIT
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0E0
*
*              INCREMENT OP COUNT
               IOPST = IOPST + 3
               IF (LS .NE. M) TEST = TEST + ABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + ABS(E(LS-1))
*
*              REPLACE STOPPING CRITERION WITH NEW ONE AS IN LAPACK
*
*              ZTEST = TEST + ABS(S(LS))
*              IF (ZTEST .NE. TEST) GO TO 420
               IF (ABS(S(LS)) .GT. EPS * TEST) GOTO 420
*
                  S(LS) = 0.0E0
C           ......EXIT
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (490,520,540,570), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0E0
*
*           INCREMENT OP COUNT
            IOPS = IOPS + ((MM1-L+1)*13 - 2)
            IF (WANTV) IOPS = IOPS + FLOAT(MM1-L+1)*6*P
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL SROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL SROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0E0
*
*           INCREMENT OP COUNT
            IOPS = IOPS + (M-L+1)*13
            IF (WANTU) IOPS = IOPS + FLOAT(M-L+1)*6*N
            DO 530 K = L, M
               T1 = S(K)
               CALL SROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL SROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
*
*           INCREMENT OP COUNT
            IOPST = IOPST + 23
            SCALE = AMAX1(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),ABS(S(L)),
     *                    ABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0E0
            C = (SM*EMM1)**2
            SHIFT = 0.0E0
            IF (B .EQ. 0.0E0 .AND. C .EQ. 0.0E0) GO TO 550
               SHIFT = SQRT(B**2+C)
               IF (B .LT. 0.0E0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
*
*           INCREMENT OP COUNT
            IOPS = IOPS + (MM1-L+1)*38
            IF (WANTV) IOPS = IOPS+FLOAT(MM1-L+1)*6*P
            IF (WANTU) IOPS = IOPS+FLOAT(MAX((MIN(MM1,N-1)-L+1),0))*6*N
            DO 560 K = L, MM1
               CALL SROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL SROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL SROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     *            CALL SROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .GE. 0.0E0) GO TO 580
               S(L) = -S(L)
*
*              INCREMENT OP COUNT
               IF (WANTV) IOPS = IOPS + P
               IF (WANTV) CALL SSCAL(P,-1.0E0,V(1,L),1)
  580       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  590       IF (L .EQ. MM) GO TO 600
C           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     *            CALL SSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     *            CALL SSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
*
*     COMPUTE FINAL OPCOUNT
      IOPS = IOPS + IOPST
      RETURN
      END
      SUBROUTINE QZHES(NM,N,A,B,MATZ,Z)
C
      INTEGER I,J,K,L,N,LB,L1,NM,NK1,NM1,NM2
      REAL A(NM,N),B(NM,N),Z(NM,N)
      REAL R,S,T,U1,U2,V1,V2,RHO
      LOGICAL MATZ
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE FIRST STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM AND THE OTHER
C     TO UPPER TRIANGULAR FORM USING ORTHOGONAL TRANSFORMATIONS.
C     IT IS USUALLY FOLLOWED BY  QZIT,  QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL GENERAL MATRIX.
C
C        B CONTAINS A REAL GENERAL MATRIX.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO.
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO.
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS IF
C          MATZ HAS BEEN SET TO .TRUE.  OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z ..........
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 J = 1, N
C
         DO 2 I = 1, N
            Z(I,J) = 0.0E0
    2    CONTINUE
C
         Z(J,J) = 1.0E0
    3 CONTINUE
C     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0E0
C
         DO 20 I = L1, N
            S = S + ABS(B(I,L))
   20    CONTINUE
C
         IF (S .EQ. 0.0E0) GO TO 100
         S = S + ABS(B(L,L))
         R = 0.0E0
C
         DO 25 I = L, N
            B(I,L) = B(I,L) / S
            R = R + B(I,L)**2
   25    CONTINUE
C
         R = SIGN(SQRT(R),B(L,L))
         B(L,L) = B(L,L) + R
         RHO = R * B(L,L)
C
         DO 50 J = L1, N
            T = 0.0E0
C
            DO 30 I = L, N
               T = T + B(I,L) * B(I,J)
   30       CONTINUE
C
            T = -T / RHO
C
            DO 40 I = L, N
               B(I,J) = B(I,J) + T * B(I,L)
   40       CONTINUE
C
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0E0
C
            DO 60 I = L, N
               T = T + B(I,L) * A(I,J)
   60       CONTINUE
C
            T = -T / RHO
C
            DO 70 I = L, N
               A(I,J) = A(I,J) + T * B(I,L)
   70       CONTINUE
C
   80    CONTINUE
C
         B(L,L) = -S * R
C
         DO 90 I = L1, N
            B(I,L) = 0.0E0
   90    CONTINUE
C
  100 CONTINUE
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + REAL( 8*N**2 + 17*N + 24 )*REAL( N-1 ) / 3.0E0
*     ----------------------- END TIMING CODE --------------------------
*
C     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      IF (N .EQ. 2) GO TO 170
      NM2 = N - 2
C
      DO 160 K = 1, NM2
         NK1 = NM1 - K
C     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     .......... ZERO A(L+1,K) ..........
            S = ABS(A(L,K)) + ABS(A(L1,K))
            IF (S .EQ. 0.0E0) GO TO 150
            U1 = A(L,K) / S
            U2 = A(L1,K) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 110 J = K, N
               T = A(L,J) + U2 * A(L1,J)
               A(L,J) = A(L,J) + T * V1
               A(L1,J) = A(L1,J) + T * V2
  110       CONTINUE
C
            A(L1,K) = 0.0E0
C
            DO 120 J = L, N
               T = B(L,J) + U2 * B(L1,J)
               B(L,J) = B(L,J) + T * V1
               B(L1,J) = B(L1,J) + T * V2
  120       CONTINUE
C     .......... ZERO B(L+1,L) ..........
            S = ABS(B(L1,L1)) + ABS(B(L1,L))
            IF (S .EQ. 0.0E0) GO TO 150
            U1 = B(L1,L1) / S
            U2 = B(L1,L) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 130 I = 1, L1
               T = B(I,L1) + U2 * B(I,L)
               B(I,L1) = B(I,L1) + T * V1
               B(I,L) = B(I,L) + T * V2
  130       CONTINUE
C
            B(L1,L) = 0.0E0
C
            DO 140 I = 1, N
               T = A(I,L1) + U2 * A(I,L)
               A(I,L1) = A(I,L1) + T * V1
               A(I,L) = A(I,L) + T * V2
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               T = Z(I,L1) + U2 * Z(I,L)
               Z(I,L1) = Z(I,L1) + T * V1
               Z(I,L) = Z(I,L) + T * V2
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      IF( MATZ ) THEN
         OPS = OPS + REAL( 11*N + 20 )*REAL( N-1 )*REAL( N-2 )
      ELSE
         OPS = OPS + REAL( 8*N + 20 )*REAL( N-1 )*REAL( N-2 )
      END IF
*     ----------------------- END TIMING CODE --------------------------
*
  170 RETURN
      END
      SUBROUTINE QZIT(NM,N,A,B,EPS1,MATZ,Z,IERR)
C
      INTEGER I,J,K,L,N,EN,K1,K2,LD,LL,L1,NA,NM,ISH,ITN,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      REAL A(NM,N),B(NM,N),Z(NM,N)
      REAL R,S,T,A1,A2,A3,EP,SH,U1,U2,U3,V1,V2,V3,ANI,A11,
     X       A12,A21,A22,A33,A34,A43,A44,BNI,B11,B12,B22,B33,B34,
     X       B44,EPSA,EPSB,EPS1,ANORM,BNORM,EPSLON
      LOGICAL MATZ,NOTLAS
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
      REAL               OPST
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE SECOND STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN D-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM.
C     IT REDUCES THE HESSENBERG MATRIX TO QUASI-TRIANGULAR FORM USING
C     ORTHOGONAL TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX.  IT IS USUALLY PRECEDED BY  QZHES  AND
C     FOLLOWED BY  QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER HESSENBERG MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C        Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  QZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED TO QUASI-TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL ARE STILL ZERO AND NO TWO
C          CONSECUTIVE SUBDIAGONAL ELEMENTS ARE NONZERO.
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  THE LOCATION B(N,1) IS USED TO STORE
C          EPS1 TIMES THE NORM OF B FOR LATER USE BY  QZVAL  AND  QZVEC.
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... COMPUTE EPSA,EPSB ..........
      ANORM = 0.0E0
      BNORM = 0.0E0
C
      DO 30 I = 1, N
         ANI = 0.0E0
         IF (I .NE. 1) ANI = ABS(A(I,I-1))
         BNI = 0.0E0
C
         DO 20 J = I, N
            ANI = ANI + ABS(A(I,J))
            BNI = BNI + ABS(B(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + REAL( N*( N+1 ) )
      OPST = 0.0E0
      ITCNT = 0
*     ----------------------- END TIMING CODE --------------------------
*
C
      IF (ANORM .EQ. 0.0E0) ANORM = 1.0E0
      IF (BNORM .EQ. 0.0E0) BNORM = 1.0E0
      EP = EPS1
      IF (EP .GT. 0.0E0) GO TO 50
C     .......... USE ROUNDOFF LEVEL IF EPS1 IS ZERO ..........
      EP = EPSLON(1.0E0)
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      LOR1 = 1
      ENORN = N
      EN = N
      ITN = 30*N
C     .......... BEGIN QZ STEP ..........
   60 IF (EN .LE. 2) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
   70 ISH = 2
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + OPST
      OPST = 0.0E0
      ITCNT = ITCNT + 1
*     ----------------------- END TIMING CODE --------------------------
*
C     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- ..........
      DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(A(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 A(L,LM1) = 0.0E0
      IF (L .LT. NA) GO TO 95
C     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ..........
      EN = LM1
      GO TO 60
C     .......... CHECK FOR SMALL TOP OF B ..........
   95 LD = L
  100 L1 = L + 1
      B11 = B(L,L)
      IF (ABS(B11) .GT. EPSB) GO TO 120
      B(L,L) = 0.0E0
      S = ABS(A(L,L)) + ABS(A(L1,L))
      U1 = A(L,L) / S
      U2 = A(L1,L) / S
      R = SIGN(SQRT(U1*U1+U2*U2),U1)
      V1 = -(U1 + R) / R
      V2 = -U2 / R
      U2 = V2 / V1
C
      DO 110 J = L, ENORN
         T = A(L,J) + U2 * A(L1,J)
         A(L,J) = A(L,J) + T * V1
         A(L1,J) = A(L1,J) + T * V2
         T = B(L,J) + U2 * B(L1,J)
         B(L,J) = B(L,J) + T * V1
         B(L1,J) = B(L1,J) + T * V2
  110 CONTINUE
C
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + REAL( 12*( ENORN+1-L ) + 11 )
*     ----------------------- END TIMING CODE --------------------------
      IF (L .NE. 1) A(L,LM1) = -A(L,LM1)
      LM1 = L
      L = L1
      GO TO 90
  120 A11 = A(L,L) / B11
      A21 = A(L1,L) / B11
      IF (ISH .EQ. 1) GO TO 140
C     .......... ITERATION STRATEGY ..........
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10) GO TO 155
C     .......... DETERMINE TYPE OF SHIFT ..........
      B22 = B(L1,L1)
      IF (ABS(B22) .LT. EPSB) B22 = EPSB
      B33 = B(NA,NA)
      IF (ABS(B33) .LT. EPSB) B33 = EPSB
      B44 = B(EN,EN)
      IF (ABS(B44) .LT. EPSB) B44 = EPSB
      A33 = A(NA,NA) / B33
      A34 = A(NA,EN) / B44
      A43 = A(EN,NA) / B33
      A44 = A(EN,EN) / B44
      B34 = B(NA,EN) / B44
      T = 0.5E0 * (A43 * B34 - A33 - A44)
      R = T * T + A34 * A43 - A33 * A44
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + REAL( 16 )
*     ----------------------- END TIMING CODE --------------------------
      IF (R .LT. 0.0E0) GO TO 150
C     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ..........
      ISH = 1
      R = SQRT(R)
      SH = -T + R
      S = -T - R
      IF (ABS(S-A44) .LT. ABS(SH-A44)) SH = S
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS OF A.
C                FOR L=EN-2 STEP -1 UNTIL LD DO -- ..........
      DO 130 LL = LD, ENM2
         L = ENM2 + LD - LL
         IF (L .EQ. LD) GO TO 140
         LM1 = L - 1
         L1 = L + 1
         T = A(L,L)
         IF (ABS(B(L,L)) .GT. EPSB) T = T - SH * B(L,L)
*        --------------------- BEGIN TIMING CODE -----------------------
         IF (ABS(A(L,LM1)) .LE. ABS(T/A(L1,L)) * EPSA) THEN
            OPST = OPST + REAL( 5 + 4*( LL+1-LD ) )
            GO TO 100
         END IF
*        ---------------------- END TIMING CODE ------------------------
  130 CONTINUE
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + REAL( 5 + 4*( ENM2+1-LD ) )
*     ----------------------- END TIMING CODE --------------------------
C
  140 A1 = A11 - SH
      A2 = A21
      IF (L .NE. LD) A(L,LM1) = -A(L,LM1)
      GO TO 160
C     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ..........
  150 A12 = A(L,L1) / B22
      A22 = A(L1,L1) / B22
      B12 = B(L,L1) / B22
      A1 = ((A33 - A11) * (A44 - A11) - A34 * A43 + A43 * B34 * A11)
     X     / A21 + A12 - A11 * B12
      A2 = (A22 - A11) - A21 * B12 - (A33 - A11) - (A44 - A11)
     X     + A43 * B34
      A3 = A(L1+1,L1) / B22
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + REAL( 25 )
*     ----------------------- END TIMING CODE --------------------------
      GO TO 160
C     .......... AD HOC SHIFT ..........
  155 A1 = 0.0E0
      A2 = 1.0E0
      A3 = 1.1605E0
  160 ITS = ITS + 1
      ITN = ITN - 1
      IF (.NOT. MATZ) LOR1 = LD
C     .......... MAIN LOOP ..........
      DO 260 K = L, NA
         NOTLAS = K .NE. NA .AND. ISH .EQ. 2
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
         LL = MIN0(EN,K1+ISH)
         IF (NOTLAS) GO TO 190
C     .......... ZERO A(K+1,K-1) ..........
         IF (K .EQ. L) GO TO 170
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
  170    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0E0) GO TO 70
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 180 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            T = B(K,J) + U2 * B(K1,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
  180    CONTINUE
C
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + REAL( 11 + 12*( ENORN+1-KM1 ) )
*        ---------------------- END TIMING CODE ------------------------
         IF (K .NE. L) A(K1,KM1) = 0.0E0
         GO TO 240
C     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) ..........
  190    IF (K .EQ. L) GO TO 200
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
         A3 = A(K2,KM1)
  200    S = ABS(A1) + ABS(A2) + ABS(A3)
         IF (S .EQ. 0.0E0) GO TO 260
         U1 = A1 / S
         U2 = A2 / S
         U3 = A3 / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 210 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J) + U3 * A(K2,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            A(K2,J) = A(K2,J) + T * V3
            T = B(K,J) + U2 * B(K1,J) + U3 * B(K2,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
            B(K2,J) = B(K2,J) + T * V3
  210    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + REAL( 17 + 20*( ENORN+1-KM1 ) )
*        ---------------------- END TIMING CODE ------------------------
C
         IF (K .EQ. L) GO TO 220
         A(K1,KM1) = 0.0E0
         A(K2,KM1) = 0.0E0
C     .......... ZERO B(K+2,K+1) AND B(K+2,K) ..........
  220    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))
         IF (S .EQ. 0.0E0) GO TO 240
         U1 = B(K2,K2) / S
         U2 = B(K2,K1) / S
         U3 = B(K2,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 230 I = LOR1, LL
            T = A(I,K2) + U2 * A(I,K1) + U3 * A(I,K)
            A(I,K2) = A(I,K2) + T * V1
            A(I,K1) = A(I,K1) + T * V2
            A(I,K) = A(I,K) + T * V3
            T = B(I,K2) + U2 * B(I,K1) + U3 * B(I,K)
            B(I,K2) = B(I,K2) + T * V1
            B(I,K1) = B(I,K1) + T * V2
            B(I,K) = B(I,K) + T * V3
  230    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + REAL( 17 + 20*( LL+1-LOR1 ) )
*        ---------------------- END TIMING CODE ------------------------
C
         B(K2,K) = 0.0E0
         B(K2,K1) = 0.0E0
         IF (.NOT. MATZ) GO TO 240
C
         DO 235 I = 1, N
            T = Z(I,K2) + U2 * Z(I,K1) + U3 * Z(I,K)
            Z(I,K2) = Z(I,K2) + T * V1
            Z(I,K1) = Z(I,K1) + T * V2
            Z(I,K) = Z(I,K) + T * V3
  235    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + REAL( 10*N )
*        ---------------------- END TIMING CODE ------------------------
C     .......... ZERO B(K+1,K) ..........
  240    S = ABS(B(K1,K1)) + ABS(B(K1,K))
         IF (S .EQ. 0.0E0) GO TO 260
         U1 = B(K1,K1) / S
         U2 = B(K1,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 250 I = LOR1, LL
            T = A(I,K1) + U2 * A(I,K)
            A(I,K1) = A(I,K1) + T * V1
            A(I,K) = A(I,K) + T * V2
            T = B(I,K1) + U2 * B(I,K)
            B(I,K1) = B(I,K1) + T * V1
            B(I,K) = B(I,K) + T * V2
  250    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + REAL( 11 + 12*( LL+1-LOR1 ) )
*        ---------------------- END TIMING CODE ------------------------
C
         B(K1,K) = 0.0E0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            T = Z(I,K1) + U2 * Z(I,K)
            Z(I,K1) = Z(I,K1) + T * V1
            Z(I,K) = Z(I,K) + T * V2
  255    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + REAL( 6*N )
*        ---------------------- END TIMING CODE ------------------------
C
  260 CONTINUE
C     .......... END QZ STEP ..........
      GO TO 70
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
C     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC ..........
 1001 IF (N .GT. 1) B(N,1) = EPSB
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + OPST
      OPST = 0.0E0
*     ----------------------- END TIMING CODE --------------------------
*
      RETURN
      END
      SUBROUTINE QZVAL(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z)
C
      INTEGER I,J,N,EN,NA,NM,NN,ISW
      REAL A(NM,N),B(NM,N),ALFR(N),ALFI(N),BETA(N),Z(NM,N)
      REAL C,D,E,R,S,T,AN,A1,A2,BN,CQ,CZ,DI,DR,EI,TI,TR,U1,
     X       U2,V1,V2,A1I,A11,A12,A2I,A21,A22,B11,B12,B22,SQI,SQR,
     X       SSI,SSR,SZI,SZR,A11I,A11R,A12I,A12R,A22I,A22R,EPSB
      LOGICAL MATZ
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
      REAL               OPST, OPST2
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE THIRD STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN QUASI-TRIANGULAR FORM AND THE OTHER IN UPPER TRIANGULAR FORM.
C     IT REDUCES THE QUASI-TRIANGULAR MATRIX FURTHER, SO THAT ANY
C     REMAINING 2-BY-2 BLOCKS CORRESPOND TO PAIRS OF COMPLEX
C     EIGENVALUES, AND RETURNS QUANTITIES WHOSE RATIOS GIVE THE
C     GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY  QZHES
C     AND  QZIT  AND MAY BE FOLLOWED BY  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C          LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C          COMPUTED AND SAVED IN  QZIT.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C        Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTIONS BY QZHES
C          AND QZIT, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED FURTHER TO A QUASI-TRIANGULAR MATRIX
C          IN WHICH ALL NONZERO SUBDIAGONAL ELEMENTS CORRESPOND TO
C          PAIRS OF COMPLEX EIGENVALUES.
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  B(N,1) IS UNALTERED.
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULAR MATRIX THAT WOULD BE
C          OBTAINED IF A WERE REDUCED COMPLETELY TO TRIANGULAR FORM
C          BY UNITARY TRANSFORMATIONS.  NON-ZERO VALUES OF ALFI OCCUR
C          IN PAIRS, THE FIRST MEMBER POSITIVE AND THE SECOND NEGATIVE.
C
C        BETA CONTAINS THE DIAGONAL ELEMENTS OF THE CORRESPONDING B,
C          NORMALIZED TO BE REAL AND NON-NEGATIVE.  THE GENERALIZED
C          EIGENVALUES ARE THEN THE RATIOS ((ALFR+I*ALFI)/BETA).
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR ALL THREE STEPS) IF MATZ HAS BEEN SET TO .TRUE.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      EPSB = B(N,1)
      ISW = 1
C     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
C                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = 0.0E0
      OPST2 = 0.0E0
*     ----------------------- END TIMING CODE --------------------------
*
      DO 510 NN = 1, N
*
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST = OPST + OPST2
         OPST2 = 0.0E0
*        ---------------------- END TIMING CODE ------------------------
*
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 505
         IF (EN .EQ. 1) GO TO 410
         IF (A(EN,NA) .NE. 0.0E0) GO TO 420
C     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
  410    ALFR(EN) = A(EN,EN)
         IF (B(EN,EN) .LT. 0.0E0) ALFR(EN) = -ALFR(EN)
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0E0
         GO TO 510
C     .......... 2-BY-2 BLOCK ..........
  420    IF (ABS(B(NA,NA)) .LE. EPSB) GO TO 455
         IF (ABS(B(EN,EN)) .GT. EPSB) GO TO 430
         A1 = A(EN,EN)
         A2 = A(EN,NA)
         BN = 0.0E0
         GO TO 435
  430    AN = ABS(A(NA,NA)) + ABS(A(NA,EN)) + ABS(A(EN,NA))
     X      + ABS(A(EN,EN))
         BN = ABS(B(NA,NA)) + ABS(B(NA,EN)) + ABS(B(EN,EN))
         A11 = A(NA,NA) / AN
         A12 = A(NA,EN) / AN
         A21 = A(EN,NA) / AN
         A22 = A(EN,EN) / AN
         B11 = B(NA,NA) / BN
         B12 = B(NA,EN) / BN
         B22 = B(EN,EN) / BN
         E = A11 / B11
         EI = A22 / B22
         S = A21 / (B11 * B22)
         T = (A22 - E * B22) / B22
         IF (ABS(E) .LE. ABS(EI)) GO TO 431
         E = EI
         T = (A11 - E * B11) / B11
  431    C = 0.5E0 * (T - S * B12)
         D = C * C + S * (A12 - E * B12)
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST2 = OPST2 + REAL( 28 )
*        ---------------------- END TIMING CODE ------------------------
         IF (D .LT. 0.0E0) GO TO 480
C     .......... TWO REAL ROOTS.
C                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
         E = E + (C + SIGN(SQRT(D),C))
         A11 = A11 - E * B11
         A12 = A12 - E * B12
         A22 = A22 - E * B22
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST2 = OPST2 + REAL( 11 )
*        ---------------------- END TIMING CODE ------------------------
         IF (ABS(A11) + ABS(A12) .LT.
     X       ABS(A21) + ABS(A22)) GO TO 432
         A1 = A12
         A2 = A11
         GO TO 435
  432    A1 = A22
         A2 = A21
C     .......... CHOOSE AND APPLY REAL Z ..........
  435    S = ABS(A1) + ABS(A2)
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 440 I = 1, EN
            T = A(I,EN) + U2 * A(I,NA)
            A(I,EN) = A(I,EN) + T * V1
            A(I,NA) = A(I,NA) + T * V2
            T = B(I,EN) + U2 * B(I,NA)
            B(I,EN) = B(I,EN) + T * V1
            B(I,NA) = B(I,NA) + T * V2
  440    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST2 = OPST2 + REAL( 11 + 12*EN )
*        ---------------------- END TIMING CODE ------------------------
C
         IF (.NOT. MATZ) GO TO 450
C
         DO 445 I = 1, N
            T = Z(I,EN) + U2 * Z(I,NA)
            Z(I,EN) = Z(I,EN) + T * V1
            Z(I,NA) = Z(I,NA) + T * V2
  445    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST2 = OPST2 + REAL( 6*N )
*        ---------------------- END TIMING CODE ------------------------
C
  450    IF (BN .EQ. 0.0E0) GO TO 475
         IF (AN .LT. ABS(E) * BN) GO TO 455
         A1 = B(NA,NA)
         A2 = B(EN,NA)
         GO TO 460
  455    A1 = A(NA,NA)
         A2 = A(EN,NA)
C     .......... CHOOSE AND APPLY REAL Q ..........
  460    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0E0) GO TO 475
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 470 J = NA, N
            T = A(NA,J) + U2 * A(EN,J)
            A(NA,J) = A(NA,J) + T * V1
            A(EN,J) = A(EN,J) + T * V2
            T = B(NA,J) + U2 * B(EN,J)
            B(NA,J) = B(NA,J) + T * V1
            B(EN,J) = B(EN,J) + T * V2
  470    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST2 = OPST2 + REAL( 11 + 12*( N+1-NA ) )
*        ---------------------- END TIMING CODE ------------------------
C
  475    A(EN,NA) = 0.0E0
         B(EN,NA) = 0.0E0
         ALFR(NA) = A(NA,NA)
         ALFR(EN) = A(EN,EN)
         IF (B(NA,NA) .LT. 0.0E0) ALFR(NA) = -ALFR(NA)
         IF (B(EN,EN) .LT. 0.0E0) ALFR(EN) = -ALFR(EN)
         BETA(NA) = ABS(B(NA,NA))
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0E0
         ALFI(NA) = 0.0E0
         GO TO 505
C     .......... TWO COMPLEX ROOTS ..........
  480    E = E + C
         EI = SQRT(-D)
         A11R = A11 - E * B11
         A11I = EI * B11
         A12R = A12 - E * B12
         A12I = EI * B12
         A22R = A22 - E * B22
         A22I = EI * B22
         IF (ABS(A11R) + ABS(A11I) + ABS(A12R) + ABS(A12I) .LT.
     X       ABS(A21) + ABS(A22R) + ABS(A22I)) GO TO 482
         A1 = A12R
         A1I = A12I
         A2 = -A11R
         A2I = -A11I
         GO TO 485
  482    A1 = A22R
         A1I = A22I
         A2 = -A21
         A2I = 0.0E0
C     .......... CHOOSE COMPLEX Z ..........
  485    CZ = SQRT(A1*A1+A1I*A1I)
         IF (CZ .EQ. 0.0E0) GO TO 487
         SZR = (A1 * A2 + A1I * A2I) / CZ
         SZI = (A1 * A2I - A1I * A2) / CZ
         R = SQRT(CZ*CZ+SZR*SZR+SZI*SZI)
         CZ = CZ / R
         SZR = SZR / R
         SZI = SZI / R
         GO TO 490
  487    SZR = 1.0E0
         SZI = 0.0E0
  490    IF (AN .LT. (ABS(E) + EI) * BN) GO TO 492
         A1 = CZ * B11 + SZR * B12
         A1I = SZI * B12
         A2 = SZR * B22
         A2I = SZI * B22
         GO TO 495
  492    A1 = CZ * A11 + SZR * A12
         A1I = SZI * A12
         A2 = CZ * A21 + SZR * A22
         A2I = SZI * A22
C     .......... CHOOSE COMPLEX Q ..........
  495    CQ = SQRT(A1*A1+A1I*A1I)
         IF (CQ .EQ. 0.0E0) GO TO 497
         SQR = (A1 * A2 + A1I * A2I) / CQ
         SQI = (A1 * A2I - A1I * A2) / CQ
         R = SQRT(CQ*CQ+SQR*SQR+SQI*SQI)
         CQ = CQ / R
         SQR = SQR / R
         SQI = SQI / R
         GO TO 500
  497    SQR = 1.0E0
         SQI = 0.0E0
C     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
C                IF TRANSFORMATIONS WERE APPLIED ..........
  500    SSR = SQR * SZR + SQI * SZI
         SSI = SQR * SZI - SQI * SZR
         I = 1
         TR = CQ * CZ * A11 + CQ * SZR * A12 + SQR * CZ * A21
     X      + SSR * A22
         TI = CQ * SZI * A12 - SQI * CZ * A21 + SSI * A22
         DR = CQ * CZ * B11 + CQ * SZR * B12 + SSR * B22
         DI = CQ * SZI * B12 + SSI * B22
         GO TO 503
  502    I = 2
         TR = SSR * A11 - SQR * CZ * A12 - CQ * SZR * A21
     X      + CQ * CZ * A22
         TI = -SSI * A11 - SQI * CZ * A12 + CQ * SZI * A21
         DR = SSR * B11 - SQR * CZ * B12 + CQ * CZ * B22
         DI = -SSI * B11 - SQI * CZ * B12
  503    T = TI * DR - TR * DI
         J = NA
         IF (T .LT. 0.0E0) J = EN
         R = SQRT(DR*DR+DI*DI)
         BETA(J) = BN * R
         ALFR(J) = AN * (TR * DR + TI * DI) / R
         ALFI(J) = AN * T / R
         IF (I .EQ. 1) GO TO 502
*        --------------------- BEGIN TIMING CODE -----------------------
         OPST2 = OPST2 + REAL( 151 )
*        ---------------------- END TIMING CODE ------------------------
  505    ISW = 3 - ISW
  510 CONTINUE
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + ( OPST + OPST2 )
*     ----------------------- END TIMING CODE --------------------------
*
      B(N,1) = EPSB
C
      RETURN
      END
      SUBROUTINE QZVEC(NM,N,A,B,ALFR,ALFI,BETA,Z)
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN,ISW,ENM2
      REAL A(NM,N),B(NM,N),ALFR(N),ALFI(N),BETA(N),Z(NM,N)
      REAL D,Q,R,S,T,W,X,Y,DI,DR,RA,RR,SA,TI,TR,T1,T2,W1,X1,
     X       ZZ,Z1,ALFM,ALMI,ALMR,BETM,EPSB
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      REAL               ITCNT, OPS
*     ..
      INTEGER            IN2BY2
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE OPTIONAL FOURTH STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM IN
C     QUASI-TRIANGULAR FORM (IN WHICH EACH 2-BY-2 BLOCK CORRESPONDS TO
C     A PAIR OF COMPLEX EIGENVALUES) AND THE OTHER IN UPPER TRIANGULAR
C     FORM.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM AND
C     TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  QZHES,  QZIT, AND  QZVAL.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C          LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C          COMPUTED AND SAVED IN  QZIT.
C
C        ALFR, ALFI, AND BETA  ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  QZVAL.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  QZHES,  QZIT, AND  QZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        A IS UNALTERED.  ITS SUBDIAGONAL ELEMENTS PROVIDE INFORMATION
C           ABOUT THE STORAGE OF THE COMPLEX EIGENVECTORS.
C
C        B HAS BEEN DESTROYED.
C
C        ALFR, ALFI, AND BETA ARE UNALTERED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF ALFI(I) .EQ. 0.0, THE I-TH EIGENVALUE IS REAL AND
C            THE I-TH COLUMN OF Z CONTAINS ITS EIGENVECTOR.
C          IF ALFI(I) .NE. 0.0, THE I-TH EIGENVALUE IS COMPLEX.
C            IF ALFI(I) .GT. 0.0, THE EIGENVALUE IS THE FIRST OF
C              A COMPLEX PAIR AND THE I-TH AND (I+1)-TH COLUMNS
C              OF Z CONTAIN ITS EIGENVECTOR.
C            IF ALFI(I) .LT. 0.0, THE EIGENVALUE IS THE SECOND OF
C              A COMPLEX PAIR AND THE (I-1)-TH AND I-TH COLUMNS
C              OF Z CONTAIN THE CONJUGATE OF ITS EIGENVECTOR.
C          EACH EIGENVECTOR IS NORMALIZED SO THAT THE MODULUS
C          OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      EPSB = B(N,1)
      ISW = 1
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
*        --------------------- BEGIN TIMING CODE -----------------------
         IN2BY2 = 0
*        ---------------------- END TIMING CODE ------------------------
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 795
         IF (ALFI(EN) .NE. 0.0E0) GO TO 710
C     .......... REAL VECTOR ..........
         M = EN
         B(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
         ALFM = ALFR(M)
         BETM = BETA(M)
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = BETM * A(I,I) - ALFM * B(I,I)
            R = 0.0E0
C
            DO 610 J = M, EN
  610       R = R + (BETM * A(I,J) - ALFM * B(I,J)) * B(J,EN)
C
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 630
            IF (BETM * A(I,I-1) .EQ. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 690
  630       M = I
            IF (ISW .EQ. 2) GO TO 640
C     .......... REAL 1-BY-1 BLOCK ..........
            T = W
            IF (W .EQ. 0.0E0) T = EPSB
            B(I,EN) = -R / T
            GO TO 700
C     .......... REAL 2-BY-2 BLOCK ..........
  640       X = BETM * A(I,I+1) - ALFM * B(I,I+1)
            Y = BETM * A(I+1,I)
            Q = W * ZZ - X * Y
            T = (X * S - ZZ * R) / Q
            B(I,EN) = T
*           ------------------- BEGIN TIMING CODE ----------------------
            IN2BY2 = IN2BY2 + 1
*           -------------------- END TIMING CODE -----------------------
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            B(I+1,EN) = (-R - W * T) / X
            GO TO 690
  650       B(I+1,EN) = (-S - Y * T) / ZZ
  690       ISW = 3 - ISW
  700    CONTINUE
C     .......... END REAL VECTOR ..........
*        --------------------- BEGIN TIMING CODE -----------------------
         OPS = OPS + ( 5.0E0/2.0E0 )*REAL( ( EN+2 )*( EN-1 ) + IN2BY2 )
*        ---------------------- END TIMING CODE ------------------------
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
         ALMR = ALFR(M)
         ALMI = ALFI(M)
         BETM = BETA(M)
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         Y = BETM * A(EN,NA)
         B(NA,NA) = -ALMI * B(EN,EN) / Y
         B(NA,EN) = (ALMR * B(EN,EN) - BETM * A(EN,EN)) / Y
         B(EN,NA) = 0.0E0
         B(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 795
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 790 II = 1, ENM2
            I = NA - II
            W = BETM * A(I,I) - ALMR * B(I,I)
            W1 = -ALMI * B(I,I)
            RA = 0.0E0
            SA = 0.0E0
C
            DO 760 J = M, EN
               X = BETM * A(I,J) - ALMR * B(I,J)
               X1 = -ALMI * B(I,J)
               RA = RA + X * B(J,NA) - X1 * B(J,EN)
               SA = SA + X * B(J,EN) + X1 * B(J,NA)
  760       CONTINUE
C
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 770
            IF (BETM * A(I,I-1) .EQ. 0.0E0) GO TO 770
            ZZ = W
            Z1 = W1
            R = RA
            S = SA
            ISW = 2
            GO TO 790
  770       M = I
            IF (ISW .EQ. 2) GO TO 780
C     .......... COMPLEX 1-BY-1 BLOCK ..........
            TR = -RA
            TI = -SA
  773       DR = W
            DI = W1
C     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ..........
  775       IF (ABS(DI) .GT. ABS(DR)) GO TO 777
            RR = DI / DR
            D = DR + DI * RR
            T1 = (TR + TI * RR) / D
            T2 = (TI - TR * RR) / D
            GO TO (787,782), ISW
  777       RR = DR / DI
            D = DR * RR + DI
            T1 = (TR * RR + TI) / D
            T2 = (TI * RR - TR) / D
            GO TO (787,782), ISW
C     .......... COMPLEX 2-BY-2 BLOCK ..........
  780       X = BETM * A(I,I+1) - ALMR * B(I,I+1)
            X1 = -ALMI * B(I,I+1)
            Y = BETM * A(I+1,I)
            TR = Y * RA - W * R + W1 * S
            TI = Y * SA - W * S - W1 * R
            DR = W * ZZ - W1 * Z1 - X * Y
            DI = W * Z1 + W1 * ZZ - X1 * Y
*           ------------------- BEGIN TIMING CODE ----------------------
            IN2BY2 = IN2BY2 + 1
*           -------------------- END TIMING CODE -----------------------
            IF (DR .EQ. 0.0E0 .AND. DI .EQ. 0.0E0) DR = EPSB
            GO TO 775
  782       B(I+1,NA) = T1
            B(I+1,EN) = T2
            ISW = 1
            IF (ABS(Y) .GT. ABS(W) + ABS(W1)) GO TO 785
            TR = -RA - X * B(I+1,NA) + X1 * B(I+1,EN)
            TI = -SA - X * B(I+1,EN) - X1 * B(I+1,NA)
            GO TO 773
  785       T1 = (-R - ZZ * B(I+1,NA) + Z1 * B(I+1,EN)) / Y
            T2 = (-S - ZZ * B(I+1,EN) - Z1 * B(I+1,NA)) / Y
  787       B(I,NA) = T1
            B(I,EN) = T2
  790    CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
         OPS = OPS + REAL( ( 6*EN-7 )*( EN-2 ) + 31*IN2BY2 )
*        ---------------------- END TIMING CODE ------------------------
C     .......... END COMPLEX VECTOR ..........
  795    ISW = 3 - ISW
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 1 DO -- ..........
      DO 880 JJ = 1, N
         J = N + 1 - JJ
C
         DO 880 I = 1, N
            ZZ = 0.0E0
C
            DO 860 K = 1, J
  860       ZZ = ZZ + Z(I,K) * B(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + REAL( N**2 )*REAL( N+1 )
*     ------------------------ END TIMING CODE -------------------------
C     .......... NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1.
C                (ISW IS 1 INITIALLY FROM BEFORE) ..........
*     ------------------------ BEGIN TIMING CODE -----------------------
      IN2BY2 = 0
*     ------------------------- END TIMING CODE ------------------------
      DO 950 J = 1, N
         D = 0.0E0
         IF (ISW .EQ. 2) GO TO 920
         IF (ALFI(J) .NE. 0.0E0) GO TO 945
C
         DO 890 I = 1, N
            IF (ABS(Z(I,J)) .GT. D) D = ABS(Z(I,J))
  890    CONTINUE
C
         DO 900 I = 1, N
  900    Z(I,J) = Z(I,J) / D
C
         GO TO 950
C
  920    DO 930 I = 1, N
            R = ABS(Z(I,J-1)) + ABS(Z(I,J))
            IF (R .NE. 0.0E0) R = R * SQRT((Z(I,J-1)/R)**2
     X                                     +(Z(I,J)/R)**2)
            IF (R .GT. D) D = R
  930    CONTINUE
C
         DO 940 I = 1, N
            Z(I,J-1) = Z(I,J-1) / D
            Z(I,J) = Z(I,J) / D
  940    CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
         IN2BY2 = IN2BY2 + 1
*        ----------------------- END TIMING CODE -----------------------
C
  945    ISW = 3 - ISW
  950 CONTINUE
*     ------------------------ BEGIN TIMING CODE -----------------------
      OPS = OPS + REAL( N*( N + 5*IN2BY2 ) )
*     ------------------------- END TIMING CODE ------------------------
C
      RETURN
      END
