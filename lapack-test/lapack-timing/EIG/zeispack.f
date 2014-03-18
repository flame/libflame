      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      DOUBLE PRECISION S,ARS,AIS,BRS,BIS
      S = DABS(BR) + DABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
      SUBROUTINE CINVIT(NM,N,AR,AI,WR,WI,SELECT,MM,M,ZR,ZI,
     X                  IERR,RM1,RM2,RV1,RV2)
C
      INTEGER I,J,K,M,N,S,II,MM,MP,NM,UK,IP1,ITS,KM1,IERR
      DOUBLE PRECISION AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,MM),
     X       ZI(NM,MM),RM1(N,N),RM2(N,N),RV1(N),RV2(N)
      DOUBLE PRECISION X,Y,EPS3,NORM,NORMV,EPSLON,GROWTO,ILAMBD,PYTHAG,
     X       RLAMBD,UKROOT
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
      DOUBLE PRECISION OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE CX INVIT
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A COMPLEX UPPER
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
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVALUES OF THE MATRIX.  THE EIGENVALUES MUST BE
C          STORED IN A MANNER IDENTICAL TO THAT OF SUBROUTINE  COMLR,
C          WHICH RECOGNIZES POSSIBLE SPLITTING OF THE MATRIX.
C
C        SELECT SPECIFIES THE EIGENVECTORS TO BE FOUND.  THE
C          EIGENVECTOR CORRESPONDING TO THE J-TH EIGENVALUE IS
C          SPECIFIED BY SETTING SELECT(J) TO .TRUE..
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVECTORS TO BE FOUND.
C
C     ON OUTPUT
C
C        AR, AI, WI, AND SELECT ARE UNALTERED.
C
C        WR MAY HAVE BEEN ALTERED SINCE CLOSE EIGENVALUES ARE PERTURBED
C          SLIGHTLY IN SEARCHING FOR INDEPENDENT EIGENVECTORS.
C
C        M IS THE NUMBER OF EIGENVECTORS ACTUALLY FOUND.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVECTORS.  THE EIGENVECTORS ARE NORMALIZED
C          SO THAT THE COMPONENT OF LARGEST MAGNITUDE IS 1.
C          ANY VECTOR WHICH FAILS THE ACCEPTANCE TEST IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -(2*N+1)   IF MORE THAN MM EIGENVECTORS HAVE BEEN SPECIFIED,
C          -K         IF THE ITERATION CORRESPONDING TO THE K-TH
C                     VALUE FAILS,
C          -(N+K)     IF BOTH ERROR SITUATIONS OCCUR.
C
C        RM1, RM2, RV1, AND RV2 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE GUESSVEC APPEARS IN CINVIT IN LINE.
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
C
*
*     GET ULP FROM DLAMCH FOR NEW SMALL PERTURBATION AS IN LAPACK
      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH, ULP
      IF (N.LE.0) RETURN
      ULP = DLAMCH( 'EPSILON' )
C
*
*     INITIALIZE
      OPST = 0
      IERR = 0
      UK = 0
      S = 1
C
      DO 980 K = 1, N
         IF (.NOT. SELECT(K)) GO TO 980
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (AR(UK+1,UK) .EQ. 0.0D0 .AND. AI(UK+1,UK) .EQ. 0.0D0)
     X         GO TO 140
  120    CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
  140    NORM = 0.0D0
         MP = 1
C
*
*        INCREMENT OPCOUNT FOR LOOP 180
         OPS = OPS + 6*UK*(UK-1)
         DO 180 I = 1, UK
            X = 0.0D0
C
            DO 160 J = MP, UK
  160       X = X + PYTHAG(AR(I,J),AI(I,J))
C
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
         IF (NORM .EQ. 0.0D0) NORM = 1.0D0
*         EPS3 = EPSLON(NORM)
*
*        INCREMENT OPCOUNT FOR EPS3, UKROOT
         OPST = OPST + 3
         EPS3 = NORM*ULP
C     .......... GROWTO IS THE CRITERION FOR GROWTH ..........
         UKROOT = UK
         UKROOT = DSQRT(UKROOT)
         GROWTO = 0.1D0 / UKROOT
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
            IF (SELECT(I) .AND. DABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     X         DABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
C
*
*        INCREMENT OPCOUNT FOR LOOP 260.
         OPST = OPST + 2*(K-1)
         WR(K) = RLAMBD
C     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I
C                AND INITIAL COMPLEX VECTOR ..........
  280    MP = 1
C
*
*        INCREMENT OP COUNT FOR LOOP 320
         OPS = OPS + 2*UK
         DO 320 I = 1, UK
C
            DO 300 J = MP, UK
               RM1(I,J) = AR(I,J)
               RM2(I,J) = AI(I,J)
  300       CONTINUE
C
            RM1(I,I) = RM1(I,I) - RLAMBD
            RM2(I,I) = RM2(I,I) - ILAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
         IF (UK .EQ. 1) GO TO 420
C
*
*        INCREMENT OP COUNT FOR LOOP 400
         OPS = OPS + (52+4*UK)*(UK-1)
         DO 400 I = 2, UK
            MP = I - 1
            IF (PYTHAG(RM1(I,MP),RM2(I,MP)) .LE.
     X          PYTHAG(RM1(MP,MP),RM2(MP,MP))) GO TO 360
C
            DO 340 J = MP, UK
               Y = RM1(I,J)
               RM1(I,J) = RM1(MP,J)
               RM1(MP,J) = Y
               Y = RM2(I,J)
               RM2(I,J) = RM2(MP,J)
               RM2(MP,J) = Y
  340       CONTINUE
C
  360       IF (RM1(MP,MP) .EQ. 0.0D0 .AND. RM2(MP,MP) .EQ. 0.0D0)
     X         RM1(MP,MP) = EPS3
            CALL CDIV(RM1(I,MP),RM2(I,MP),RM1(MP,MP),RM2(MP,MP),X,Y)
            IF (X .EQ. 0.0D0 .AND. Y .EQ. 0.0D0) GO TO 400
C
            DO 380 J = I, UK
               RM1(I,J) = RM1(I,J) - X * RM1(MP,J) + Y * RM2(MP,J)
               RM2(I,J) = RM2(I,J) - X * RM2(MP,J) - Y * RM1(MP,J)
  380       CONTINUE
C
  400    CONTINUE
C
  420    IF (RM1(UK,UK) .EQ. 0.0D0 .AND. RM2(UK,UK) .EQ. 0.0D0)
     X      RM1(UK,UK) = EPS3
         ITS = 0
C     .......... BACK SUBSTITUTION
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0D0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C
            DO 680 J = IP1, UK
               X = X - RM1(I,J) * RV1(J) + RM2(I,J) * RV2(J)
               Y = Y - RM1(I,J) * RV2(J) - RM2(I,J) * RV1(J)
  680       CONTINUE
C
  700       CALL CDIV(X,Y,RM1(I,I),RM2(I,I),RV1(I),RV2(I))
  720    CONTINUE
*
*        INCREMENT OP COUNT FOR BACK SUBSTITUTION LOOP 720
         OPS = OPS + 4*UK*(UK+3)
C     .......... ACCEPTANCE TEST FOR EIGENVECTOR
C                AND NORMALIZATION ..........
         ITS = ITS + 1
         NORM = 0.0D0
         NORMV = 0.0D0
C
*
*        INCREMENT OP COUNT ACCEPTANCE TEST
         OPS = OPS + 19*UK
         DO 780 I = 1, UK
            X = PYTHAG(RV1(I),RV2(I))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
C
         IF (NORM .LT. GROWTO) GO TO 840
C     .......... ACCEPT VECTOR ..........
         X = RV1(J)
         Y = RV2(J)
C
*
*        INCREMENT OP COUNT ACCEPT VECTOR LOOP 820
         OPS = OPS + 16*UK
         DO 820 I = 1, UK
            CALL CDIV(RV1(I),RV2(I),X,Y,ZR(I,S),ZI(I,S))
  820    CONTINUE
C
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0D0)
         RV1(1) = EPS3
C
         DO 860 I = 2, UK
  860    RV1(I) = Y
C
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         GO TO 660
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
         IERR = -K
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
            ZR(I,S) = 0.0D0
            ZI(I,S) = 0.0D0
  920    CONTINUE
C
  940    S = S + 1
  980 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE COMQR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
C
      INTEGER I,J,L,N,EN,LL,NM,IGH,ITN,ITS,LOW,LP1,ENM1,IERR
      DOUBLE PRECISION HR(NM,N),HI(NM,N),WR(N),WI(N)
      DOUBLE PRECISION SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
      COMMON /PYTHOP/ OPST
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR, NUM. MATH. 12, 369-376(1968) BY MARTIN
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
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
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN
C          INFORMATION ABOUT THE UNITARY TRANSFORMATIONS USED IN
C          THE REDUCTION BY  CORTH, IF PERFORMED.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
C          CALLING  COMQR  IF SUBSEQUENT CALCULATION OF
C          EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*
      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH, UNFL,OVFL,ULP,SMLNUM,SMALL
      INTRINSIC MAX, MIN
*
      IF (N.LE.0) RETURN
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
      NORM = 0.0D0
      DO 5 J = LOW, IGH
         SR = 0.0D0
         DO 4 I = LOW, MIN(IGH,J+1)
              SR = SR + PYTHAG(HR(I,J),HI(I,J))
  4      CONTINUE
         NORM = MAX(NORM, SR)
  5   CONTINUE
*
*     GET SMALL FOR NEW CONVERGENCE CRITERION AS IN LAPACK
*
      UNFL = DLAMCH( 'SAFE MINIMUM' )
      OVFL = DLAMCH( 'OVERFLOW' )
      ULP = DLAMCH( 'EPSILON' )*DLAMCH( 'BASE' )
      SMLNUM = MAX( UNFL*( N / ULP ), N / ( ULP*OVFL ) )
      SMALL = MAX( SMLNUM, ULP*NORM )
*
*
*     INITIALIZE
      ITCNT = 0
      OPST = 0
      IERR = 0
      IF (LOW .EQ. IGH) GO TO 180
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
      L = LOW + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 170
         OPS = OPS + (6*(IGH-LOW+1)+32)*(IGH-L+1)
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 155 J = I, IGH
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = LOW, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
*         TST2 = TST1 + ABS(HR(L,L-1))
*         IF (TST2 .EQ. TST1) GO TO 300
         TST2 = ABS(HR(L,L-1))
         IF ( TST2 .LE. MIN(ULP*TST1,SMALL) ) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
         OPS = OPS + 4*(EN-L+1)
      IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
*
*        INCREMENT OPCOUNT FOR FOMING SHIFT
         OPST = OPST + 58
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL CSROOT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 360
         OPS = OPS + 2*EN
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
        ITCNT = 30*N - ITN
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
*
*        INCREMENT OPCOUNT FOR REDUCING TO TRIANGULAR, LOOP 500
         OPS = OPS + (EN-LP1+1)*(61+10*(EN-LP1))
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0D0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0D0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, EN
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
*
*        INCREMENT OPCOUNT
         OPST = OPST +20
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = L, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
  600 CONTINUE
*
*        INCREMENT OPCOUNT FOR INVERSE OPERATION LOOP 600
         OPS = OPS + 10*(EN-LP1+1)*(EN+LP1)
C
      IF (SI .EQ. 0.0D0) GO TO 240
C
*
*        INCREMENT OP COUNT FOR LOOP 630
         OPS = OPS + 6*(EN-L+1)
      DO 630 I = L, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 CONTINUE
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE COMQR2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1,
     X        ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      DOUBLE PRECISION HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       ORTR(IGH),ORTI(IGH)
      DOUBLE PRECISION SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
      COMMON /PYTHOP/ OPST
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
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
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0D0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
*     THE ORIGINAL DO STATEMENTS
*
*         DO 840 I = 1, ENM1
*         DO 820 J = IP1, N
*         DO 880 JJ = LOW, ENM1
*
*     HAVE BEEN CHANGED TO
*
*         DO 840 I = 1, N
*         DO 820 J = I, N
*         DO 880 JJ = LOW, N
*
*     ACCORDING TO BURT GARBOW'S SUGGESTION ON NA-NET.
*     ZHAOJUN BAI, NOV.28, 1989
C     ------------------------------------------------------------------
C
*
      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH, UNFL,OVFL,ULP,SMLNUM,SMALL
      INTRINSIC MAX, MIN
*
      IF (N.LE.0) RETURN
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
      NORM = 0.0D0
      DO 5 J = 1,N
         SR = 0.0D0
         DO 4 I = 1, MIN(N,J+1)
              SR = SR + PYTHAG(HR(I,J),HI(I,J))
  4      CONTINUE
         NORM = MAX(NORM, SR)
  5   CONTINUE
*
*     GET SMALL FOR NEW CONVERGENCE CRITERION AS IN LAPACK
*
      UNFL = DLAMCH( 'SAFE MINIMUM' )
      OVFL = DLAMCH( 'OVERFLOW' )
      ULP = DLAMCH( 'EPSILON' )*DLAMCH( 'BASE' )
      SMLNUM = MAX( UNFL*( N / ULP ), N / ( ULP*OVFL ) )
      SMALL = MAX( SMLNUM, MIN( ( NORM*SMLNUM )*NORM, ULP*NORM ) )
*
*
*     INITIALIZE
      ITCNT = 0
      OPST = 0
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 101 J = 1, N
C
         DO 100 I = 1, N
            ZR(I,J) = 0.0D0
            ZI(I,J) = 0.0D0
  100    CONTINUE
         ZR(J,J) = 1.0D0
  101 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND) 180, 150, 105
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0D0 .AND. ORTI(I) .EQ. 0.0D0) GO TO 140
         IF (HR(I,I-1) .EQ. 0.0D0 .AND. HI(I,I-1) .EQ. 0.0D0) GO TO 140
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 110 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  110    CONTINUE
C
*
*        INCREMENT OP COUNT FOR LOOP 130
         OPS = OPS + (16*(IGH-I+1)+2)*(IGH-I+1)
         DO 130 J = I, IGH
            SR = 0.0D0
            SI = 0.0D0
C
            DO 115 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 120 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
*
*        INCREMENT OP COUNT FOR COMPUTING NORM IN LOOP 140
         OPS = OPS + 3*IEND
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 170
         OPS = OPS + (12*(IGH-LOW+1)+42)*(IGH-L+1)
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 155 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
         DO 165 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  165    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
*         TST2 = TST1 + ABS(HR(L,L-1))
*         IF (TST2 .EQ. TST1) GO TO 300
         TST2 = ABS(HR(L,L-1))
         IF ( TST2 .LE. MIN(ULP*TST1,SMALL) ) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
         OPS = OPS + 4*(EN-L+1)
      IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
*
*        INCREMENT OPCOUNT FOR FOMING SHIFT
         OPST = OPST + 58
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL CSROOT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 360
         OPS = OPS + 2*(EN-LOW+1)
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
        ITCNT = 30*N - ITN
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
*
*        INCREMENT OPCOUNT FOR REDUCING TO TRIANGULAR, LOOP 500
         OPS = OPS + (EN-LP1+1)*(61+10*(EN-LP1))
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0D0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0D0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
*
*        INCREMENT OP COUNT
         OPST = OPST +20
      IF (EN .EQ. N) GO TO 540
      IP1 = EN + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 520
         OPST = OPST + 6*(N-IP1+1)
      DO 520 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
         DO 590 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
C
  600 CONTINUE
*
*        INCREMENT OPCOUNT FOR INVERSE OPERATION LOOP 600
         OPS = OPS + ( 10*(EN+LP1) + 20*(IGH-LOW+1) )*(EN-LP1+1)
C
      IF (SI .EQ. 0.0D0) GO TO 240
C
*
*        INCREMENT OPCOUNT FOR LOOP 630 AND 640
         OPS = OPS + 6*EN + 6*(IGH-LOW+1)
      DO 630 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      DO 640 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0D0
C
*
*        INCREMENT OP COUNT FOR LOOP 720
         OPS = OPS + N*(N+1)/2
      DO 720 I = 1, N
C
         DO 720 J = I, N
            TR = DABS(HR(I,J)) + DABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  720 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0D0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0D0
         HI(EN,EN) = 0.0D0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
*
*        INCREMENT OP COUNT FOR COMPUT YR, .. IN LOOP 780
         OPS = OPS + 22*ENM1
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = 0.0D0
            ZZI = 0.0D0
            IP1 = I + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 740
         OPS = OPS + 7*(EN-IP1+1)
            DO 740 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0D0 .OR. YI .NE. 0.0D0) GO TO 765
               TST1 = NORM
               YR = TST1
  760          YR = 0.01D0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GO TO 760
  765       CONTINUE
            CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
*
*        INCREMENT OP COUNT FOR CDIV
         OPST = OPST + 16
C     .......... OVERFLOW CONTROL ..........
            TR = DABS(HR(I,EN)) + DABS(HI(I,EN))
            IF (TR .EQ. 0.0D0) GO TO 780
            TST1 = TR
            TST2 = TST1 + 1.0D0/TST1
            IF (TST2 .GT. TST1) GO TO 780
*
*        INCREMENT OP COUNT FOR LOOP 770
         OPS = OPS + 2*(EN-I+1)
            DO 770 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  770       CONTINUE
C
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = I, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
*
*        INCREMENT OP COUNT FOR LOOP 880
         OPS = OPS + 8*(M-LOW+1)*(IGH-LOW+1)
         DO 880 I = LOW, IGH
            ZZR = 0.0D0
            ZZI = 0.0D0
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
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
      SUBROUTINE CORTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      DOUBLE PRECISION AR(NM,N),AI(NM,N),ORTR(IGH),ORTI(IGH)
      DOUBLE PRECISION F,G,H,FI,FR,SCALE,PYTHAG
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ OPS, ITCNT
      COMMON /PYTHOP/ OPST
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION OPS, ITCNT, OPST
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
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
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
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
      IF (N.LE.0) RETURN
***
*     INITIALIZE
      OPST = 0
***
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0D0
         ORTR(M) = 0.0D0
         ORTI(M) = 0.0D0
         SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + DABS(AR(I,M-1)) + DABS(AI(I,M-1))
***
*        INCREMENT OPCOUNT FOR LOOP 90
         OPS = OPS + 2*(IGH-M+1)
***
C
         IF (SCALE .EQ. 0.0D0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORTR(I) = AR(I,M-1) / SCALE
            ORTI(I) = AI(I,M-1) / SCALE
            H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I)
  100    CONTINUE
***
*        INCREMENT OP COUNT FOR LOOP 100 AND SQRT
         OPS = OPS + 6*(IGH-M+1) + 1
***
C
         G = DSQRT(H)
         F = PYTHAG(ORTR(M),ORTI(M))
         IF (F .EQ. 0.0D0) GO TO 103
         H = H + F * G
         G = G / F
         ORTR(M) = (1.0D0 + G) * ORTR(M)
         ORTI(M) = (1.0D0 + G) * ORTI(M)
         OPST = OPST + 7
         GO TO 105
C
  103    ORTR(M) = G
         AR(M,M-1) = SCALE
C     .......... FORM (I-(U*UT)/H) * A ..........
  105    DO 130 J = M, N
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J)
               FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J)
  110       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 120 I = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I)
               AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J)
               FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J)
  140       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 150 J = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J)
               AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J)
  150       CONTINUE
C
  160    CONTINUE
***
*        INCREMENT OP COUNT FOR LOOPS 130 AND 160
         OPS = OPS + (IGH+N-M+1)*((IGH-M+1)*16 + 2)
         OPST = OPST + 4
***
C
         ORTR(M) = SCALE * ORTR(M)
         ORTI(M) = SCALE * ORTI(M)
         AR(M,M-1) = -G * AR(M,M-1)
         AI(M,M-1) = -G * AI(M,M-1)
  180 CONTINUE
      OPS = OPS + OPST
C
  200 RETURN
      END
      SUBROUTINE CSROOT(XR,XI,YR,YI)
      DOUBLE PRECISION XR,XI,YR,YI
C
C     (YR,YI) = COMPLEX SQRT(XR,XI)
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C
      DOUBLE PRECISION S,TR,TI,PYTHAG
      TR = XR
      TI = XI
      S = DSQRT(0.5D0*(PYTHAG(TR,TI) + DABS(TR)))
      IF (TR .GE. 0.0D0) YR = S
      IF (TI .LT. 0.0D0) S = -S
      IF (TR .LE. 0.0D0) YI = S
      IF (TR .LT. 0.0D0) YR = 0.5D0*(TI/YI)
      IF (TR .GT. 0.0D0) YI = 0.5D0*(TI/YR)
      RETURN
      END
      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
C
      INTEGER I,J,K,L,M,N,NM
      DOUBLE PRECISION AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      DOUBLE PRECISION H,S,SI
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT.
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED.
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRIDI.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  HTRIDI  IN THEIR
C          FULL LOWER TRIANGLES EXCEPT FOR THE DIAGONAL OF AR.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
*
      OPS = OPS + MAX( 0.0D0, 8*M*DBLE(N)**2 - 2*M*DBLE(N) - 4*M )
*
C     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. ..........
      DO 50 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
C
      IF (N .EQ. 1) GO TO 200
C     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. 0.0D0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0D0
            SI = 0.0D0
C
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
      DOUBLE PRECISION F,G,H,FI,GI,HH,SI,SCALE,PYTHAG
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT.
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED.
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX
C     TO A REAL SYMMETRIC TRIDIAGONAL MATRIX USING
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX HERMITIAN INPUT MATRIX.
C          ONLY THE LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION IN THEIR FULL LOWER
C          TRIANGLES.  THEIR STRICT UPPER TRIANGLES AND THE
C          DIAGONAL OF AR ARE UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
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
*
      OPS = OPS + MAX( 0.0D0, (16.0D0/3.0D0)*DBLE(N)**3 + 3*DBLE(N)**2
     $                      + (56.0D0/3.0D0)*N - 61 )
*
      TAU(1,N) = 1.0D0
      TAU(2,N) = 0.0D0
C
      DO 100 I = 1, N
  100 D(I) = AR(I,I)
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(AR(I,K)) + DABS(AI(I,K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
         TAU(1,L) = 1.0D0
         TAU(2,L) = 0.0D0
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AI(I,K) = AI(I,K) / SCALE
            H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = DSQRT(H)
         E(I) = SCALE * G
         F = PYTHAG(AR(I,L),AI(I,L))
C     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
         IF (F .EQ. 0.0D0) GO TO 160
         TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
         SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0D0 + G / F
         AR(I,L) = G * AR(I,L)
         AI(I,L) = G * AI(I,L)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         AR(I,L) = G
  170    F = 0.0D0
C
         DO 240 J = 1, L
            G = 0.0D0
            GI = 0.0D0
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
               G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
               GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
  180       CONTINUE
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
               GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
  200       CONTINUE
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AI(I,J)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
C
            DO 260 K = 1, J
               AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K)
     X                           + FI * TAU(2,K) + GI * AI(I,K)
               AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K)
     X                           - FI * E(K) - GI * AR(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AI(I,K) = SCALE * AI(I,K)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AI(I,I) = SCALE * DSQRT(H)
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE IMTQL1(N,D,E,IERR)
*
*     EISPACK ROUTINE
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEQR.
*
C
      INTEGER I,J,L,M,N,II,MML,IERR
      DOUBLE PRECISION D(N),E(N)
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      DOUBLE PRECISION EPS, TST
      DOUBLE PRECISION DLAMCH
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
      DOUBLE PRECISION   ITCNT, OPS, OPST
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
      EPS = DLAMCH( 'EPSILON' )
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0D0
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
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = PYTHAG(G,1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
            OPS = OPS + 7
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. 0.0D0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
            OPS = OPS + MML*14 + 1
*
*        INCREMENT ITERATION COUNTER
            ITCNT = ITCNT + 1
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0D0
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
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEQR.
*
C
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      DOUBLE PRECISION EPS, TST
      DOUBLE PRECISION DLAMCH
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
      DOUBLE PRECISION   ITCNT, OPS, OPST
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
      EPS = DLAMCH( 'EPSILON' )
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0D0
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
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = PYTHAG(G,1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
            OPS = OPS + 7
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. 0.0D0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
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
         E(M) = 0.0D0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
            OPS = OPS + MML*( 14+6*N ) + 1
*
*        INCREMENT ITERATION COUNTER
            ITCNT = ITCNT + 1
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0D0
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
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
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
      DOUBLE PRECISION   OPST
*     ..
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
*
*     INCREMENT OPST
      OPST = OPST + 2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
*
*        INCREMENT OPST
            OPST = OPST + 8
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      DOUBLE PRECISION A,B,C,EPS
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
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      RETURN
      END
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEBZ.
*
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM
      DOUBLE PRECISION D(N),E(N),E2(N),W(MM),RV4(N),RV5(N)
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(MM)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
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
      DOUBLE PRECISION ONE
      PARAMETER        ( ONE = 1.0D0 )
      DOUBLE PRECISION RELFAC
      PARAMETER        ( RELFAC = 2.0D0 )
      DOUBLE PRECISION ATOLI, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP
      DOUBLE PRECISION DLAMCH, PIVMIN
      EXTERNAL DLAMCH
*        INITIALIZE ITERATION COUNT.
            ITCNT = 0
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
CCC         TST1 = DABS(D(I)) + DABS(D(I-1))
CCC         TST2 = TST1 + DABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
         TMP1 = E( I )**2
         IF( ABS( D(I)*D(I-1) )*ULP**2+SAFEMN.LE.TMP1 )
     $      GO TO 40
   20    E2(I) = 0.0D0
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
      U = 0.0D0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = DABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
            OPS = OPS + ( N-P+1 )*2
C
  140 X1 = EPSLON(DMAX1(DABS(XU),DABS(X0)))
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
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
  300    X1 = (XU + X0) * 0.5D0
CCC         IF ((X0 - XU) .LE. DABS(EPS1)) GO TO 420
CCC         TST1 = 2.0D0 * (DABS(XU) + DABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
         TMP1 = ABS( X0 - XU )
         TMP2 = MAX( ABS( X0 ), ABS( XU ) )
         IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) )
     $      GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0D0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = DABS(E(I)) / EPSLON(1.0D0)
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
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
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,EPSLON,
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
      DOUBLE PRECISION   ITCNT, OPS, OPST
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
C          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0
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
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
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
      ORDER = 1.0D0 - E2(1)
      Q = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140
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
         XU = 1.0D0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0D0
         GO TO 870
  490    NORM = DABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = DMAX1(NORM, DABS(D(I))+DABS(E(I)))
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
         EPS2 = 1.0D-3 * NORM
         EPS3 = EPSLON(NORM)
         UK = Q - P + 1
         EPS4 = UK * EPS3
         UK = EPS4 / DSQRT(UK)
*           INCREMENT OPCOUNT FOR COMPUTING CRITERIA.
               OPS = OPS + ( Q-IP+4 )
         S = P
  505    GROUP = 0
         GO TO 520
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
  520    V = 0.0D0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
C     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
*           INCREMENT OPCOUNT FOR ELIMINATION.
               OPS = OPS + ( Q-P+1 )*5
C
         IF (U .EQ. 0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
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
            XU = 0.0D0
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
  700    NORM = 0.0D0
C
         DO 720 I = P, Q
  720    NORM = NORM + DABS(RV6(I))
*           INCREMENT OPCOUNT FOR COMPUTING NORM.
               OPS = OPS + ( Q-P+1 )
C
         IF (NORM .GE. 1.0D0) GO TO 840
C     .......... FORWARD SUBSTITUTION ..........
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0D0) GO TO 740
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
         XU = 0.0D0
         GO TO 870
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0D0
C
         DO 860 I = P, Q
  860    U = PYTHAG(U,RV6(I))
C
         XU = 1.0D0 / U
C
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0D0
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
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEBZ.
*
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(M)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
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
      DOUBLE PRECISION ONE
      PARAMETER        ( ONE = 1.0D0 )
      DOUBLE PRECISION RELFAC
      PARAMETER        ( RELFAC = 2.0D0 )
      DOUBLE PRECISION ATOLI, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP
      DOUBLE PRECISION DLAMCH, PIVMIN
      EXTERNAL DLAMCH
*        INITIALIZE ITERATION COUNT.
            ITCNT = 0
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
      PIVMIN = ONE
      DO 40 I = 1, N
         X1 = U
         U = 0.0D0
         IF (I .NE. N) U = DABS(E(I+1))
         XU = DMIN1(D(I)-(X1+U),XU)
         X0 = DMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
CCC         TST1 = DABS(D(I)) + DABS(D(I-1))
CCC         TST2 = TST1 + DABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
         TMP1 = E( I )**2
         IF( ABS( D(I)*D(I-1) )*ULP**2+SAFEMN.LE.TMP1 ) THEN
            PIVMIN = MAX( PIVMIN, TMP1 )
            GO TO 40
         END IF
   20    E2(I) = 0.0D0
   40 CONTINUE
      PIVMIN = PIVMIN*SAFEMN
      TNORM = MAX( ABS( XU ), ABS( X0 ) )
      ATOLI = ULP*TNORM
*        INCREMENT OPCOUNT FOR DETERMINING IF MATRIX SPLITS.
            OPS = OPS + 9*( N-1 )
C
      X1 = N
      X1 = X1 * EPSLON(DMAX1(DABS(XU),DABS(X0)))
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
      X1 = XU + (X0 - XU) * 0.5D0
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
      U = 0.0D0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = DABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
            OPS = OPS + ( N-P+1 )*2
C
  140 X1 = EPSLON(DMAX1(DABS(XU),DABS(X0)))
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
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
  300    X1 = (XU + X0) * 0.5D0
CCC         IF ((X0 - XU) .LE. DABS(EPS1)) GO TO 420
CCC         TST1 = 2.0D0 * (DABS(XU) + DABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
         TMP1 = ABS( X0 - XU )
         TMP2 = MAX( ABS( X0 ), ABS( XU ) )
         IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) )
     $      GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0D0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = DABS(E(I)) / EPSLON(1.0D0)
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
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
      SUBROUTINE ZSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      COMPLEX*16 X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, IOPS IS ONLY INCREMENTED
*     IOPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO IOPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON /LATIME/ IOPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION IOPS, ITCNT, IOPST
*     ..
C
C
C     ZSVDC IS A SUBROUTINE TO REDUCE A COMPLEX*16 NXP MATRIX X BY
C     UNITARY TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         COMPLEX*16(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY ZSVDC.
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
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V
C                   (SEE BELOW).
C
C         WORK      COMPLEX*16(N).
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
C                        A.GE.2    RETURNS THE FIRST MIN(N,P)
C                                  LEFT SINGULAR VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         COMPLEX*16(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         COMPLEX*16(P).
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         COMPLEX*16(LDU,K), WHERE LDU.GE.N.  IF JOBA.EQ.1
C                                   THEN K.EQ.N, IF JOBA.GE.2 THEN
C                                   K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF LEFT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.GT.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         COMPLEX*16(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOBB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WHTH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = CTRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (CTRANS(U)
C                   IS THE CONJUGATE-TRANSPOSE OF U).  THUS THE
C                   SINGULAR VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C              CORRECTION TO SHIFT CALCULATION MADE 2/85.
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     ZSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL ZDROT
C     BLAS ZAXPY,ZDOTC,ZSCAL,ZSWAP,DZNRM2,DROTG
C     FORTRAN DABS,DMAX1,CDABS,DCMPLX
C     FORTRAN DCONJG,MAX0,MIN0,MOD,DSQRT
C
C     INTERNAL VARIABLES
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     *        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      COMPLEX*16 ZDOTC,T,R
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,DZNRM2,SCALE,SHIFT,SL,SM,SN,
     *                 SMM1,T1,TEST
*     DOUBLE PRECISION ZTEST
      LOGICAL WANTU,WANTV
C
      COMPLEX*16 CSIGN,ZDUM,ZDUM1,ZDUM2
      DOUBLE PRECISION CABS1
*
*     DECLARE EPS AND DLAMCH FOR NEW STOPPING CRITERION
      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH, EPS
*
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
      CSIGN(ZDUM1,ZDUM2) = CDABS(ZDUM1)*(ZDUM2/CDABS(ZDUM2))
*
*     GET EPS FROM DLAMCH FOR NEW STOPPING CRITERION
      IF (N.LE.0 .OR. P.LE.0) RETURN
      EPS = DLAMCH( 'EPSILON' )
*
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
            IOPS = IOPS + (4*(N-L+1)+2)
            S(L) = DCMPLX(DZNRM2(N-L+1,X(L,L),1),0.0D0)
            IF (CABS1(S(L)) .EQ. 0.0D0) GO TO 10
               IF (CABS1(X(L,L)) .NE. 0.0D0) S(L) = CSIGN(S(L),X(L,L))
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (6*(N-L+1)+23)
               CALL ZSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
               X(L,L) = (1.0D0,0.0D0) + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (CABS1(S(L)) .EQ. 0.0D0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (16*(N-L)+26)
               T = -ZDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL ZAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = DCONJG(X(L,J))
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
            IOPS = IOPS + (4*(P-L)+3)
            E(L) = DCMPLX(DZNRM2(P-L,E(LP1),1),0.0D0)
            IF (CABS1(E(L)) .EQ. 0.0D0) GO TO 80
               IF (CABS1(E(LP1)) .NE. 0.0D0) E(L) = CSIGN(E(L),E(LP1))
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (6*(P-L)+23)
               CALL ZSCAL(P-L,1.0D0/E(L),E(LP1),1)
               E(LP1) = (1.0D0,0.0D0) + E(LP1)
   80       CONTINUE
            E(L) = -DCONJG(E(L))
            IF (LP1 .GT. N .OR. CABS1(E(L)) .EQ. 0.0D0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = (0.0D0,0.0D0)
   90          CONTINUE
*
*              INCREMENT OP COUNT
               IOPS = IOPS + DBLE(16*(N-L)+9)*(P-L)
               DO 100 J = LP1, P
                  CALL ZAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL ZAXPY(N-L,DCONJG(-E(J)/E(LP1)),WORK(LP1),1,
     *                       X(LP1,J),1)
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
      IF (N .LT. M) S(M) = (0.0D0,0.0D0)
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = (0.0D0,0.0D0)
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = (0.0D0,0.0D0)
  180       CONTINUE
            U(J,J) = (1.0D0,0.0D0)
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (CABS1(S(L)) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (DBLE(16*(N-L)+25)*(NCU-L)+6*(N-L)+9)
               DO 210 J = LP1, NCU
                  T = -ZDOTC(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL ZAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL ZSCAL(N-L+1,(-1.0D0,0.0D0),U(L,L),1)
               U(L,L) = (1.0D0,0.0D0) + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = (0.0D0,0.0D0)
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = (0.0D0,0.0D0)
  260          CONTINUE
               U(L,L) = (1.0D0,0.0D0)
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
            IF (CABS1(E(L)) .EQ. 0.0D0) GO TO 320
*
*              INCREMENT OP COUNT
               IOPS = IOPS + (DBLE(16*(P-L)+9)*(P-L)+1)
               DO 310 J = LP1, P
                  T = -ZDOTC(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL ZAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = (0.0D0,0.0D0)
  330       CONTINUE
            V(L,L) = (1.0D0,0.0D0)
  340    CONTINUE
  350 CONTINUE
C
C     TRANSFORM S AND E SO THAT THEY ARE DOUBLE PRECISION.
C
*
*     INCREMENT OP COUNT
      IOPS = IOPS + (2*M-1)
      DO 380 I = 1, M
         IF (CABS1(S(I)) .EQ. 0.0D0) GO TO 360
*
*           INCREMENT OP COUNT
            IOPS = IOPS + 23
            IF (WANTU) IOPS = IOPS + 6*N
            T = DCMPLX(CDABS(S(I)),0.0D0)
            R = S(I)/T
            S(I) = T
            IF (I .LT. M) E(I) = E(I)/R
            IF (WANTU) CALL ZSCAL(N,R,U(1,I),1)
  360    CONTINUE
C     ...EXIT
         IF (I .EQ. M) GO TO 390
         IF (CABS1(E(I)) .EQ. 0.0D0) GO TO 370
*
*           INCREMENT OP COUNT
            IOPS = IOPS + 20
            IF (WANTV) IOPS = IOPS + 6*P
            T = DCMPLX(CDABS(E(I)),0.0D0)
            R = T/E(I)
            E(I) = T
            S(I+1) = S(I+1)*R
            IF (WANTV) CALL ZSCAL(P,R,V(1,I+1),1)
  370    CONTINUE
  380 CONTINUE
  390 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
*
*     INITIALIZE ITERATION COUNTER
      ITCNT = 0
      ITER = 0
  400 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 660
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
*
*        UPDATE ITERATION COUNTER
         ITCNT = ITER
         IF (ITER .LT. MAXIT) GO TO 410
            INFO = M
C     ......EXIT
            GO TO 660
  410    CONTINUE
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
         DO 430 LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) GO TO 440
*
*           INCREMENT OP COUNT
            IOPST = IOPST + 17
            TEST = CDABS(S(L)) + CDABS(S(L+1))
*
*           REPLACE STOPPING CRITERION WITH NEW ONE
*
*           ZTEST = TEST + CDABS(E(L))
*           IF (ZTEST .NE. TEST) GO TO 420
            IF (CDABS(E(L)) .GT. EPS * TEST) GOTO 420
*
               E(L) = (0.0D0,0.0D0)
C        ......EXIT
               GO TO 440
  420       CONTINUE
  430    CONTINUE
  440    CONTINUE
         IF (L .NE. M - 1) GO TO 450
            KASE = 4
         GO TO 520
  450    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 470 LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) GO TO 480
               TEST = 0.0D0
*
*              INCREMENT OP COUNT
               IOPST = IOPST + 18
               IF (LS .NE. M) TEST = TEST + CDABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + CDABS(E(LS-1))
*
*              REPLACE STOPPING CRITERION WITH NEW ONE AS IN LAPACK
*
*              ZTEST = TEST + CDABS(S(LS))
*              IF (ZTEST .NE. TEST) GO TO 460
               IF (CDABS(S(LS))  .GT. EPS * TEST) GOTO 460
*
                  S(LS) = (0.0D0,0.0D0)
C           ......EXIT
                  GO TO 480
  460          CONTINUE
  470       CONTINUE
  480       CONTINUE
            IF (LS .NE. L) GO TO 490
               KASE = 3
            GO TO 510
  490       CONTINUE
            IF (LS .NE. M) GO TO 500
               KASE = 1
            GO TO 510
  500       CONTINUE
               KASE = 2
               L = LS
  510       CONTINUE
  520    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (530, 560, 580, 610), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  530    CONTINUE
            MM1 = M - 1
            F = DREAL(E(M-1))
            E(M-1) = (0.0D0,0.0D0)
*
*           INCREMENT OP COUNT
            IOPS = IOPS + ((MM1-L+1)*14 - 3)
            IF (WANTV) IOPS = IOPS + DBLE(MM1-L+1)*12*P
            DO 550 KK = L, MM1
               K = MM1 - KK + L
               T1 = DREAL(S(K))
               CALL DROTG(T1,F,CS,SN)
               S(K) = DCMPLX(T1,0.0D0)
               IF (K .EQ. L) GO TO 540
                  F = -SN*DREAL(E(K-1))
                  E(K-1) = CS*E(K-1)
  540          CONTINUE
               IF (WANTV) CALL ZDROT(P,V(1,K),1,V(1,M),1,CS,SN)
  550       CONTINUE
         GO TO 650
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  560    CONTINUE
            F = DREAL(E(L-1))
            E(L-1) = (0.0D0,0.0D0)
*
*           INCREMENT OP COUNT
            IOPS = IOPS + (M-L+1)*14
            IF (WANTU) IOPS = IOPS + DBLE(M-L+1)*12*N
            DO 570 K = L, M
               T1 = DREAL(S(K))
               CALL DROTG(T1,F,CS,SN)
               S(K) = DCMPLX(T1,0.0D0)
               F = -SN*DREAL(E(K))
               E(K) = CS*E(K)
               IF (WANTU) CALL ZDROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  570       CONTINUE
         GO TO 650
C
C        PERFORM ONE QR STEP.
C
  580    CONTINUE
C
C           CALCULATE THE SHIFT.
C
*
*           INCREMENT OP COUNT
            IOPST = IOPST + 48
            SCALE = DMAX1(CDABS(S(M)),CDABS(S(M-1)),CDABS(E(M-1)),
     *                    CDABS(S(L)),CDABS(E(L)))
            SM = DREAL(S(M))/SCALE
            SMM1 = DREAL(S(M-1))/SCALE
            EMM1 = DREAL(E(M-1))/SCALE
            SL = DREAL(S(L))/SCALE
            EL = DREAL(E(L))/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 590
               SHIFT = DSQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  590       CONTINUE
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
*
*           INCREMENT OP COUNT
            IOPS = IOPS + (MM1-L+1)*46
            IF (WANTV) IOPS = IOPS+DBLE(MM1-L+1)*12*P
            IF (WANTU) IOPS = IOPS+DBLE(MAX((MIN(MM1,N-1)-L+1),0))*12*N
            DO 600 K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = DCMPLX(F,0.0D0)
               F = CS*DREAL(S(K)) + SN*DREAL(E(K))
               E(K) = CS*E(K) - SN*S(K)
               G = SN*DREAL(S(K+1))
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL ZDROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = DCMPLX(F,0.0D0)
               F = CS*DREAL(E(K)) + SN*DREAL(S(K+1))
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*DREAL(E(K+1))
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     *            CALL ZDROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  600       CONTINUE
            E(M-1) = DCMPLX(F,0.0D0)
            ITER = ITER + 1
         GO TO 650
C
C        CONVERGENCE.
C
  610    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE
C
            IF (DREAL(S(L)) .GE. 0.0D0) GO TO 620
               S(L) = -S(L)
*
*              INCREMENT OP COUNT
               IF (WANTV) IOPS = IOPS + 6*P
               IF (WANTV) CALL ZSCAL(P,(-1.0D0,0.0D0),V(1,L),1)
  620       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  630       IF (L .EQ. MM) GO TO 640
C           ...EXIT
               IF (DREAL(S(L)) .GE. DREAL(S(L+1))) GO TO 640
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     *            CALL ZSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     *            CALL ZSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 630
  640       CONTINUE
            ITER = 0
            M = M - 1
  650    CONTINUE
      GO TO 400
  660 CONTINUE
*
*     COMPUTE FINAL OPCOUNT
      IOPS = IOPS + IOPST
      RETURN
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)
C
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      DOUBLE PRECISION AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),
     1       ZI(NM,N)
      DOUBLE PRECISION R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
CC      REAL SQRT,CABS,ABS
      LOGICAL MATZ
CC      COMPLEX*16 DCMPLX
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
      DOUBLE PRECISION   OPST
      INTEGER            IOPST
*     ------------------------ END TIMING CODE -------------------------
*
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            ZR(I,J) = 0.0D0
            ZI(I,J) = 0.0D0
    2    CONTINUE
C
         ZR(I,I) = 1.0D0
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
*        ---------------------- BEGIN TIMING CODE ----------------------
         IOPST = 0
*        ----------------------- END TIMING CODE -----------------------
         L1 = L + 1
         S = 0.0D0
C
         DO 20 I = L, N
            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   20    CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
         IOPST = IOPST + 2*( N+1-L )
*        ----------------------- END TIMING CODE -----------------------
C
         IF (S .EQ. 0.0D0) GO TO 100
         RHO = 0.0D0
C
         DO 25 I = L, N
            BR(I,L) = BR(I,L) / S
            BI(I,L) = BI(I,L) / S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   25    CONTINUE
C
         R = SQRT(RHO)
         XR = ABS(DCMPLX(BR(L,L),BI(L,L)))
         IF (XR .EQ. 0.0D0) GO TO 27
*        ---------------------- BEGIN TIMING CODE ----------------------
         IOPST = IOPST + 8
*        ----------------------- END TIMING CODE -----------------------
         RHO = RHO + XR * R
         U1 = -BR(L,L) / XR
         U1I = -BI(L,L) / XR
         YR = R / XR + 1.0D0
         BR(L,L) = YR * BR(L,L)
         BI(L,L) = YR * BI(L,L)
         GO TO 28
C
   27    BR(L,L) = R
         U1 = -1.0D0
         U1I = 0.0D0
C
   28    DO 50 J = L1, N
            T = 0.0D0
            TI = 0.0D0
C
            DO 30 I = L, N
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
   30       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 40 I = L, N
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
   40       CONTINUE
C
            XI = U1 * BI(L,J) - U1I * BR(L,J)
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
            BI(L,J) = XI
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0D0
            TI = 0.0D0
C
            DO 60 I = L, N
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
   60       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 70 I = L, N
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
   70       CONTINUE
C
            XI = U1 * AI(L,J) - U1I * AR(L,J)
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
            AI(L,J) = XI
   80    CONTINUE
C
         BR(L,L) = R * S
         BI(L,L) = 0.0D0
C
         DO 90 I = L1, N
            BR(I,L) = 0.0D0
            BI(I,L) = 0.0D0
   90    CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
         OPS = OPS + ( DBLE( 16*( N-L ) + 16*N + 30 )*DBLE( N-L ) +
     $                 DBLE( 24*N + 13 + IOPST ) )
*        ----------------------- END TIMING CODE -----------------------
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
*        ---------------------- BEGIN TIMING CODE ----------------------
         OPST = 0.0D0
*        ----------------------- END TIMING CODE -----------------------
         K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (AI(N,K) .EQ. 0.0D0) GO TO 105
         R = ABS(DCMPLX(AR(N,K),AI(N,K)))
         U1 = AR(N,K) / R
         U1I = AI(N,K) / R
         AR(N,K) = R
         AI(N,K) = 0.0D0
C
         DO 103 J = K1, N
            XI = U1 * AI(N,J) - U1I * AR(N,J)
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
            AI(N,J) = XI
  103    CONTINUE
C
         XI = U1 * BI(N,N) - U1I * BR(N,N)
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
         BI(N,N) = XI
*        ---------------------- BEGIN TIMING CODE ----------------------
         OPST = OPST + DBLE( 18 + 6*( N-K ) )
*        ----------------------- END TIMING CODE -----------------------
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
            IF (S .EQ. 0.0D0) GO TO 150
*           -------------------- BEGIN TIMING CODE ---------------------
            OPST = OPST + DBLE( 18 + 20*( 2*N-K-L ) )
*           --------------------- END TIMING CODE ----------------------
            U1 = AR(L,K) / S
            U1I = AI(L,K) / S
            U2 = AR(L1,K) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            AR(L,K) = R * S
            AI(L,K) = 0.0D0
            AR(L1,K) = 0.0D0
C
            DO 110 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110       CONTINUE
C
            XR = BR(L,L)
            BR(L,L) = U1 * XR
            BI(L,L) = -U1I * XR
            BR(L1,L) = -U2 * XR
C
            DO 120 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
            IF (S .EQ. 0.0D0) GO TO 150
*           -------------------- BEGIN TIMING CODE ---------------------
            OPST = OPST + DBLE( 13 + 20*( N+L ) )
*           --------------------- END TIMING CODE ----------------------
            U1 = BR(L1,L1) / S
            U1I = BI(L1,L1) / S
            U2 = BR(L1,L) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            BR(L1,L1) = R * S
            BI(L1,L1) = 0.0D0
            BR(L1,L) = 0.0D0
C
            DO 130 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
  130       CONTINUE
C
            DO 140 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
*           -------------------- BEGIN TIMING CODE ---------------------
            OPST = OPST + 20*N
*           --------------------- END TIMING CODE ----------------------
C
            DO 145 I = 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
  145       CONTINUE
C
  150    CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
         OPS = OPS + ( OPST + DBLE( 2*( N-1-K ) ) )
*        ----------------------- END TIMING CODE -----------------------
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF CQZHES **********
      END
      SUBROUTINE CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,
     X                                       MATZ,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      DOUBLE PRECISION AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),
     X       ALFI(N),BETA(N),ZR(NM,N),ZI(NM,N)
      DOUBLE PRECISION R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,ANI,A1I,A33,
     X       A34,A43,A44,BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,
     X       B33I,B44I,EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
CC      REAL SQRT,CSQRT,ABS
      INTEGER MAX0
      LOGICAL MATZ
      DOUBLE COMPLEX Z3
CC      COMPLEX CSQRT,DCMPLX
CC      REAL REAL,AIMAG
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
      DOUBLE PRECISION   OPST
      INTEGER            IOPST
*     ------------------------ END TIMING CODE -------------------------
*
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0D0
      BNORM = 0.0D0
C
      DO 30 I = 1, N
         ANI = 0.0D0
         IF (I .NE. 1) ANI = ABS(AR(I,I-1))
         BNI = 0.0D0
C
         DO 20 J = I, N
            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. 0.0D0) ANORM = 1.0D0
      IF (BNORM .EQ. 0.0D0) BNORM = 1.0D0
      EP = EPS1
      IF (EP .GT. 0.0D0) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0D0
   40 EP = EP / 2.0D0
      IF (1.0D0 + EP .GT. 1.0D0) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COUNT OPS FOR NORMS, BUT NOT FOR CALCULATION OF "EP"
      OPS = OPS + DBLE( 2*N*( N+1 ) + 2 )
      OPST = 0.0D0
      ITCNT = 0.0D0
*     ------------------------ END TIMING CODE -------------------------
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + OPST
      OPST = 0.0D0
*     ------------------------ END TIMING CODE -------------------------
      DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 AR(L,LM1) = 0.0D0
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = ABS(DCMPLX(BR(L,L),BI(L,L)))
      IF (B11     .EQ. 0.0D0) GO TO 98
      U1 = BR(L,L) / B11
      U1I = BI(L,L) / B11
C
      DO 97 J = L, ENORN
         XI = U1 * AI(L,J) - U1I * AR(L,J)
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
         AI(L,J) = XI
         XI = U1 * BI(L,J) - U1I * BR(L,J)
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
         BI(L,J) = XI
   97 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + DBLE( 7 + 12*( ENORN+1-L ) )
*     ------------------------ END TIMING CODE -------------------------
C
      BI(L,L) = 0.0D0
   98 IF (L .NE. EN) GO TO 100
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      BR(L,L) = 0.0D0
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
      U1 = AR(L,L) / S
      U1I = AI(L,L) / S
      U2 = AR(L1,L) / S
      R = SQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1 / R
      U1I = U1I / R
      U2 = U2 / R
      AR(L,L) = R * S
      AI(L,L) = 0.0D0
C
      DO 110 J = L1, ENORN
         XR = AR(L,J)
         XI = AI(L,J)
         YR = AR(L1,J)
         YI = AI(L1,J)
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
         XR = BR(L,J)
         XI = BI(L,J)
         YR = BR(L1,J)
         YI = BI(L1,J)
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + DBLE( 15 + 40*( ENORN-L ) )
*     ------------------------ END TIMING CODE -------------------------
C
      LM1 = L
      L = L1
      GO TO 90
C     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
C     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (ABS(DCMPLX(B33,B33I)) .GE. EPSB) GO TO 122
      B33 = EPSB
      B33I = 0.0D0
  122 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (ABS(DCMPLX(B44,B44I)) .GE. EPSB) GO TO 124
      B44 = EPSB
      B44I = 0.0D0
  124 B3344 = B33 * B44 - B33I * B44I
      B3344I = B33 * B44I + B33I * B44
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I
     X    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33
     X     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
      A43 = AR(EN,NA) * B44
      A43I = AR(EN,NA) * B44I
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34 * A43 - A34I * A43I
      XI = A34 * A43I + A34I * A43
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + DBLE( 54 )
*     ------------------------ END TIMING CODE -------------------------
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 140
      YR = (A33 - SH) / 2.0D0
      YI = (A33I - SHI) / 2.0D0
      Z3 = SQRT(DCMPLX(YR**2-YI**2+XR,2.0D0*YR*YI+XI))
      U1 = DBLE(Z3)
      U1I = DIMAG(Z3)
      IF (YR * U1 + YI * U1I .GE. 0.0D0) GO TO 125
      U1 = -U1
      U1I = -U1I
  125 Z3 = (DCMPLX(SH,SHI) - DCMPLX(XR,XI) / DCMPLX(YR+U1,YI+U1I))
     X   / DCMPLX(B3344,B3344I)
      SH = DBLE(Z3)
      SHI = DIMAG(Z3)
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + DBLE( 66 )
*     ------------------------ END TIMING CODE -------------------------
      GO TO 140
C     ********** AD HOC SHIFT **********
  135 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0D0
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = AR(L,L) / B11 - SH
      A1I = AI(L,L) / B11 - SHI
      A2 = AR(L1,L) / B11
      ITS = ITS + 1
*     ----------------------- BEGIN TIMING CODE ------------------------
      ITCNT = ITCNT + 1.0D0
*     ------------------------ END TIMING CODE -------------------------
      IF (.NOT. MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = AR(K,KM1)
         A1I = AI(K,KM1)
         A2 = AR(K1,KM1)
  170    S = ABS(A1) + ABS(A1I) + ABS(A2)
         U1 = A1 / S
         U1I = A1I / S
         U2 = A2 / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
C
         DO 180 J = KM1, ENORN
            XR = AR(K,J)
            XI = AI(K,J)
            YR = AR(K1,J)
            YI = AI(K1,J)
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(K,J)
            XI = BI(K,J)
            YR = BR(K1,J)
            YI = BI(K1,J)
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
  180    CONTINUE
C
         IF (K .EQ. L) GO TO 240
         AI(K,KM1) = 0.0D0
         AR(K1,KM1) = 0.0D0
         AI(K1,KM1) = 0.0D0
C     ********** ZERO B(K+1,K) **********
  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
         U1 = BR(K1,K1) / S
         U1I = BI(K1,K1) / S
         U2 = BR(K1,K) / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = AR(K2,K1)
         AR(K2,K1) = U1 * XR
         AI(K2,K1) = -U1I * XR
         AR(K2,K) = -U2 * XR
C
  245    DO 250 I = LOR1, K1
            XR = AR(I,K1)
            XI = AI(I,K1)
            YR = AR(I,K)
            YI = AI(I,K)
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(I,K1)
            XI = BI(I,K1)
            YR = BR(I,K)
            YI = BI(I,K)
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
  250    CONTINUE
C
         BI(K1,K1) = 0.0D0
         BR(K1,K) = 0.0D0
         BI(K1,K) = 0.0D0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            XR = ZR(I,K1)
            XI = ZI(I,K1)
            YR = ZR(I,K)
            YI = ZI(I,K)
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
  255    CONTINUE
C
  260 CONTINUE
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COUNT OPS FOR STATEMENTS 140 -- 260
      IOPST = 29 + 40*( ENORN-LOR1+4 )
      IF( MATZ ) IOPST = IOPST + 20*N
      OPST = OPST + ( DBLE( N-L )*DBLE( IOPST ) + 2 )
      IF( L.LE.1 ) OPST = OPST - 40
*     ------------------------ END TIMING CODE -------------------------
*
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA) .EQ. 0.0D0) GO TO 70
      R = ABS(DCMPLX(AR(EN,NA),AI(EN,NA)))
      U1 = AR(EN,NA) / R
      U1I = AI(EN,NA) / R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0D0
C
      DO 270 J = EN, ENORN
         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
         AI(EN,J) = XI
         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
         BI(EN,J) = XI
  270 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + DBLE( 7 + 12*( EN+1-ENORN ) )
*     ------------------------ END TIMING CODE -------------------------
C
      GO TO 70
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) BR(N,1) = EPSB
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + OPST
      OPST = 0.0D0
*     ------------------------ END TIMING CODE -------------------------
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
      SUBROUTINE CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      DOUBLE PRECISION AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),
     X       ALFI(N),BETA(N),ZR(NM,N),ZI(NM,N)
      DOUBLE PRECISION R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
CC      REAL CABS
      DOUBLE COMPLEX Z3
CC      COMPLEX CMPLX
CC      REAL REAL,AIMAG
C
C
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
*     .. COMMON BLOCKS ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. SCALARS IN COMMON ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
*     ------------------------ END TIMING CODE -------------------------
*
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (N .LE. 1) GO TO 1001
      EPSB = BR(N,1)
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALMR = ALFR(EN)
         ALMI = ALFI(EN)
         BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = 0.0D0
            RI = 0.0D0
            M = I + 1
C
            DO 610 J = M, EN
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
               IF (J .EQ. EN) GO TO 605
               XI = T * BI(J,EN) + TI * BR(J,EN)
               T = T * BR(J,EN) - TI * BI(J,EN)
               TI = XI
  605          R = R + T
               RI = RI + TI
  610       CONTINUE
C
            T = ALMR * BETA(I) - BETM * ALFR(I)
            TI = ALMI * BETA(I) - BETM * ALFI(I)
            IF (T .EQ. 0.0D0 .AND. TI .EQ. 0.0D0) T = EPSB
            Z3 = DCMPLX(R,RI) / DCMPLX(T,TI)
            BR(I,EN) = DBLE(Z3)
            BI(I,EN) = DIMAG(Z3)
  700    CONTINUE
C
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
C
         DO 880 I = 1, N
C
            DO 860 K = 1, M
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)
  860       CONTINUE
C
  880 CONTINUE
C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1 **********
      DO 950 J = 1, N
         T = 0.0D0
C
         DO 930 I = 1, N
            R = ABS(DCMPLX(ZR(I,J),ZI(I,J)))
            IF (R .GT. T) T = R
  930    CONTINUE
C
         DO 940 I = 1, N
            ZR(I,J) = ZR(I,J) / T
            ZI(I,J) = ZI(I,J) / T
  940    CONTINUE
C
  950 CONTINUE
C
 1001 CONTINUE
*
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + DBLE( N )*DBLE( 14*N**2 + 15*N - 15 ) / DBLE( 2 )
*     ------------------------ END TIMING CODE -------------------------
*
      RETURN
C     ********** LAST CARD OF CQZVEC **********
      END
