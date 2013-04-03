C===================================================== REGTST.FOR
      SUBROUTINE REGTST(NSITES,LENGTH,XMOM,A,B,SEED,NSIM,NPROB,PROB,
     *  KPRINT,KOUT,RMOM,D,VOBS,VBAR,VSD,H,Z,PARA,IERR)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C*  VERSION 3.03  JUNE 2000                                            *
C*  * CHARACTER variable declarations changed to conform with          *
C*    Fortran-77 standard                                              *
C*                                                                     *
C***********************************************************************
C
C  CALCULATES THREE STATISTICS USEFUL IN REGIONAL FREQUENCY ANALYSIS
C
C  DISCORDANCY MEASURE, D(I), FOR INDIVIDUAL SITES IN A REGION.
C      LARGE VALUES MIGHT BE USED AS A FLAG TO INDICATE POTENTIAL ERRORS
C      IN THE DATA AT THE SITE.  "LARGE" MIGHT BE 3 FOR REGIONS WITH 15
C      OR MORE SITES, BUT LESS (EXACT VALUES IN ARRAY DC1) FOR SMALLER
C      REGIONS.
C
C  HETEROGENEITY MEASURES, H(J), FOR A REGION BASED UPON EITHER:-
C      J=1: THE WEIGHTED S.D. OF THE L-CVS OR
C      J=2: THE AVERAGE DISTANCE FROM THE SITE TO THE REGIONAL AVERAGE
C           ON A GRAPH OF L-CV VS. L-SKEWNESS
C      J=3: THE AVERAGE DISTANCE FROM THE SITE TO THE REGIONAL AVERAGE
C           ON A GRAPH OF L-SKEWNESS VS. L-KURTOSIS
C
C      IN PRACTICE H(1) IS PROBABLY SUFFICIENT.  A VALUE GREATER THAN
C      (SAY) 1.0 SUGGESTS THAT FURTHER SUBDIVISION OF THE REGION SHOULD
C      BE CONSIDERED AS IT MIGHT IMPROVE QUANTILE ESTIMATES.
C
C  GOODNESS-OF-FIT MEASURES, Z(K), FOR 5 CANDIDATE DISTRIBUTIONS:
C      K=1: GENERALIZED LOGISTIC
C      K=2: GENERALIZED EXTREME VALUE
C      K=3: GENERALIZED NORMAL (LOGNORMAL)
C      K=4: PEARSON TYPE III (3-PARAMETER GAMMA)
C      K=5: GENERALIZED PARETO
C
C      PROVIDED THAT THE REGION IS ACCEPTABLY CLOSE TO HOMOGENEOUS,
C      THE FIT MAY BE JUDGED ACCEPTABLE AT 10% SIGNIFICANCE LEVEL
C      IF Z(K) IS LESS THAN 1.645 IN ABSOLUTE VALUE.
C
C  FOR FURTHER DETAILS SEE J.R.M. HOSKING AND J.R. WALLIS (1997),
C  "REGIONAL FREQUENCY ANALYSIS: AN APPROACH BASED ON L-MOMENTS",
C  CAMBRIDGE UNIVERSITY PRESS, CHAPTERS 3-5.
C
C  PARAMETERS OF ROUTINE:
C  NSITES * INPUT* NUMBER OF SITES IN REGION
C  NAMES  * INPUT* CHARACTER*12 ARRAY OF LENGTH NSITES. SITE NAMES.
C  LENGTH    * INPUT* ARRAY OF LENGTH NSITES. RECORD LENGTHS AT EACH SITE.
C  XMOM   * INPUT* ARRAY OF DIMENSION (5,NSITES). ARRAY CONTAINING
C                  THE FIRST 5 SAMPLE L-MOMENTS FOR EACH SITE, IN THE
C                  ORDER MEAN, L-CV, L-SKEWNESS, L-KURTOSIS, T-5, I.E
C                  XMOM(I,J) CONTAINS THE I'TH L-MOMENT FOR SITE J.
C                    N.B. XMOM(2,.) CONTAINS L-CV, NOT THE USUAL L-2!
C  A      * INPUT* ) PARAMETERS OF
C  B      * INPUT* ) PLOTTING POSITION.
C                  NOTE: A AND B SHOULD BE THE SAME AS THE VALUES USED
C                  TO CALCULATE THE MOMENTS IN THE XMOM ARRAY.
C  SEED   * INPUT* SEED FOR RANDOM NUMBER GENERATOR. SHOULD BE A WHOLE
C                  NUMBER IN THE RANGE 2D0 TO 2147483647D0.
C  NSIM   * INPUT* NUMBER OF SIMULATED WORLDS FOR HETEROGENEITY AND
C                  GOODNESS-OF-FIT TESTS.
C                    NOTE: NSIM=0 WILL FORCE RETURN AT COMPLETION OF
C                  OUTLIER TEST.  NSIM=1 WILL SUPPRESS CALCULATION OF
C                  H AND Z STATISTICS, BUT PARAMETER AND QUANTILE
C                  ESTIMATES WILL BE FOUND.
C  NPROB  * INPUT* NUMBER OF QUANTILES TO BE CALCULATED
C  PROB   * INPUT* ARRAY OF LENGTH NPROB.  PROBABILITIES FOR WHICH
C                  QUANTILES ARE TO BE CALCULATED.
C  KPRINT * INPUT* OUTPUT FLAG. SHOULD BE SET TO
C                  0  TO SUPPRESS OUTPUT
C                  1  TO PRINT OUTPUT
C  KOUT   * INPUT* CHANNEL TO WHICH OUTPUT IS DIRECTED
C  RMOM   *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE REGIONAL
C                  WEIGHTED AVERAGE L-MOMENT RATIOS.
C  D      *OUTPUT* ARRAY OF LENGTH NSITES. ON EXIT, CONTAINS THE
C                  DISCORDANCY MEASURE (D STATISTIC) FOR EACH SITE.
C  VOBS   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE REGIONAL
C                  OBSERVED VALUES OF 3 HETEROGENEITY STATISTICS:
C                  (1) WEIGHTED S.D. OF L-CVS;
C                  (2) AVERAGE OF L-CV/L-SKEW DISTANCES;
C                  (3) AVERAGE OF L-SKEW/L-KURTOSIS DISTANCES.
C  VBAR   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE MEAN OF THE
C                  SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
C  VSD    *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE S.D. OF THE
C                  SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
C  H      *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS HETEROGENEITY
C                  MEASURES (H STATISTICS), I.E. H=(VOBS-VBAR)/VSD.
C  Z      *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS GOODNESS-OF-FIT
C                  MEASURES (Z STATISTICS) FOR 5 DISTRIBUTIONS:
C                  (1) GEN. LOGISTIC, (2) GEN. EXTREME VALUE,
C                  (3) GEN. NORMAL, (4) PEARSON TYPE III,
C                  (5) GEN. PARETO.
C  PARA   *OUTPUT* ARRAY OF DIMENSION (5,6). ON EXIT, IF NSIM.GE.1,
C                  CONTAINS PARAMETERS OF GROWTH CURVES FITTED BY THE
C                  ABOVE 5 DISTRIBUTIONS, PLUS WAKEBY.
C
C  OTHER ROUTINES USED: DERF,DIGAMD,DLGAMA,DURAND,GAMIND,PELGEV,PELGLO,
C    PELGNO,PELGPA,PELKAP,PELPE3,PELWAK,QUAGAM,QUAGEV,QUAGLO,QUAGNO,
C    QUAGPA,QUAKAP,QUAPE3,QUASTN,QUAWAK,SAMLMR,SORT
C
C  QUANTITIES DEFINED IN PARAMETER STATEMENT:
C  MAXNS  - MUST BE AT LEAST AS LARGE AS NSITES
C  MAXREC - MUST BE AT LEAST AS LARGE AS EACH ELEMENT OF ARRAY LEN
C  MAXQ   - MUST BE AT LEAST AS LARGE AS NPROB
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXNS=200,MAXQ=30,MAXREC=200)
C
      CHARACTER*1 BLANK,STAR,LOOK1,LOOK2
C      CHARACTER*12 NAMES(NSITES)
      CHARACTER*18 DISTRI(6)
      DOUBLE PRECISION D(NSITES),DC1(14),DC2(18),H(3),PARA(5,6),
     *  PROB(NPROB),Q(MAXQ),RMOM(5),RPARA(4),SMAT(3,3),TMOM(4),T4FIT(5),
     *  VBAR(3),VOBS(3),VSD(3),WORK(MAXNS,3),X(MAXREC),XMOM(5,NSITES),
     *  Z(5)
      INTEGER LENGTH(NSITES)

      INTEGER :: IERR
C CHANGES
C -------
C  * Removed NAMES(NSITES) and all mentions of it since this
C    is non-trivial to wrap with f2py. DH
C  * Added IERR error code for passing to modified version of
C    SAMLMR and store internal regtst errors. DH


C F2PY comments
C -------------
Cf2py integer, dimension(nsites), intent(in) :: length
Cf2py integer, depend(length), intent(hide) :: nsites = len(length)
Cf2py double precision, dimension(5,nsites), intent(in) :: xmom
Cf2py double precision, intent(in) :: a, b, seed
Cf2py integer, intent(in) :: nsim
Cf2py double precision, dimension(nprob), intent(in) :: prob
Cf2py integer, intent(hide), depend(prob) :: nprob = len(prob)
Cf2py integer, intent(in) :: kprint, kout
Cf2py double precision,intent(out) :: rmom(5), D(nsites), 
Cf2py double precision, dimension(3), intent(out) :: vobs, vbar, vsd, H
Cf2py double precision, intent(out) :: Z(5), para(5,6)
Cf2py integer, intent(out) :: ierr 


      DATA BLANK/' '/,STAR/'*'/
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
      DATA DISTRI/
     *  'GEN. LOGISTIC     ','GEN. EXTREME VALUE','GEN. NORMAL       ',
     *  'PEARSON TYPE III  ','GEN. PARETO       ','WAKEBY            '/
C
C         COEFFICIENTS OF POWER-SERIES APPROXIMATIONS OF TAU-4 IN TERMS
C         OF TAU-3, FOR THE FIRST 5 DISTRIBUTIONS IN ARRAY DISTRI
C
      DATA GLOC0,GLOC2/0.16667D0,0.83333D0/
      DATA GEVC0,GEVC1,GEVC2,GEVC3,GEVC4,GEVC5,GEVC6/
     *  0.10701D0, 0.11090D0, 0.84838D0,-0.06669D0,
     *  0.00567D0,-0.04208D0, 0.03763D0/
      DATA GNOC0,GNOC2,GNOC4,GNOC6,GNOC8/
     *  0.12282D0,0.77518D0,0.12279D0,-0.13638D0,0.11368D0/
      DATA PE3C0,PE3C2,PE3C4,PE3C6,PE3C8/
     *  0.12240D0,0.30115D0,0.95812D0,-0.57488D0,0.19383D0/
      DATA GPAC1,GPAC2,GPAC3,GPAC4/
     *  0.20196D0,0.95924D0,-0.20096D0,0.04061D0/
C
C         CRITICAL VALUES FOR D, H AND Z STATISTICS
C
      DATA DC1/4*3D0,1.3330D0,1.6481D0,1.9166D0,2.1401D0,2.3287D0,
     *               2.4906D0,2.6321D0,2.7573D0,2.8694D0,2.9709D0/
      DATA DC2/4*4D0,1.3333D0,1.6648D0,1.9821D0,2.2728D0,2.5337D0,
     *               2.7666D0,2.9748D0,3.1620D0,3.3310D0,3.4844D0,
     *               3.6246D0,3.7532D0,3.8718D0,3.9816D0/
      DATA HCRIT1,HCRIT2/1D0,2D0/
      DATA ZCRIT/1.645D0/
C
C         INITIALIZE ARRAYS
C
      NMAX=0
      SUMLEN=0
      DO 10 I=1,NSITES
      NREC=LENGTH(I)
      IF(NREC.GT.NMAX)NMAX=NREC
      SUMLEN=SUMLEN+NREC
   10 D(I)=ZERO
      DO 20 K=1,3
      VOBS(K)=ZERO
      VBAR(K)=ZERO
      VSD(K)=ZERO
      H(K)=ZERO
   20 CONTINUE
      DO 30 IDIST=1,5
   30 Z(IDIST)=ZERO
      DO 40 IPARA=1,5
      DO 40 IDIST=1,6
   40 PARA(IPARA,IDIST)=ZERO
      IF(NSITES.GT.MAXNS)GOTO 1000
C
C         CALCULATE THE WEIGHTED MEAN OF L-CV, L-SKEW, L-KURTOSIS
C
      DO 60 K=2,5
      RMOM(K)=ZERO
      DO 50 I=1,NSITES
   50 RMOM(K)=RMOM(K)+LENGTH(I)*XMOM(K,I)
   60 RMOM(K)=RMOM(K)/SUMLEN
      RMOM(1)=ONE
C
C         CALCULATE SUM OF SQUARES MATRIX
C
      IF(NSITES.LE.3)GOTO 135
      SUM2=ZERO
      SUM3=ZERO
      SUM4=ZERO
      DO 70 I=1,NSITES
      SUM2=SUM2+XMOM(2,I)
      SUM3=SUM3+XMOM(3,I)
      SUM4=SUM4+XMOM(4,I)
   70 CONTINUE
      SUM2=SUM2/NSITES
      SUM3=SUM3/NSITES
      SUM4=SUM4/NSITES
      DO 80 I=1,NSITES
      WORK(I,1)=XMOM(2,I)-SUM2
      WORK(I,2)=XMOM(3,I)-SUM3
      WORK(I,3)=XMOM(4,I)-SUM4
   80 CONTINUE
      DO 100 J=1,3
      DO 100 K=J,3
      SMAT(J,K)=ZERO
      DO 90 I=1,NSITES
   90 SMAT(J,K)=SMAT(J,K)+WORK(I,J)*WORK(I,K)
  100 CONTINUE
C
C         INVERT SUM OF SQUARES MATRIX
C
      DO 110 K=1,3
      IF(SMAT(1,1).LE.ZERO)GOTO 1030
      TEMP0=ONE/SMAT(1,1)
      TEMP1=-SMAT(1,2)*TEMP0
      TEMP2=-SMAT(1,3)*TEMP0
      IF(K.GT.2)TEMP1=-TEMP1
      IF(K.GT.1)TEMP2=-TEMP2
      SMAT(1,1)=SMAT(2,2)+TEMP1*SMAT(1,2)
      SMAT(1,2)=SMAT(2,3)+TEMP1*SMAT(1,3)
      SMAT(2,2)=SMAT(3,3)+TEMP2*SMAT(1,3)
      SMAT(1,3)=TEMP1
      SMAT(2,3)=TEMP2
      SMAT(3,3)=TEMP0
  110 CONTINUE
      SMAT(2,1)=SMAT(1,2)
      SMAT(3,1)=SMAT(1,3)
      SMAT(3,2)=SMAT(2,3)
C
C         CALCULATE DISCORDANCY MEASURES (D STATISTICS)
C
      FACTOR=NSITES/THREE
      DO 130 I=1,NSITES
      DO 120 J=1,3
      DO 120 K=1,3
  120 D(I)=D(I)+WORK(I,J)*WORK(I,K)*SMAT(J,K)
      D(I)=D(I)*FACTOR
      WORK(I,1)=D(I)
  130 CONTINUE
      CALL SORT(WORK(1,1),NSITES)
      GOTO 140
  135 DO 138 I=1,NSITES
  138 D(I)=ONE
C
C         PRINT DISCORDANCY MEASURES
C
  140 CONTINUE
      IF(KPRINT.LE.0)GOTO 160
      WRITE(KOUT,6000)
      DCRIT1=DC1(1)
      DCRIT2=DC2(1)
      IF(NSITES.LE.14)DCRIT1=DC1(NSITES)
      IF(NSITES.LE.18)DCRIT2=DC2(NSITES)
      KSTART=1
      DO 150 I=1,NSITES
      LOOK1=BLANK
      LOOK2=BLANK
      IF(D(I).GE.DCRIT1)LOOK1=STAR
      IF(D(I).GE.DCRIT2)LOOK2=STAR
      IF(D(I).LT.DCRIT1)KSTART=KSTART+1
      WRITE(KOUT,6010)I,LENGTH(I),(XMOM(K,I),K=2,4),
     *  D(I),LOOK1,LOOK2
  150 CONTINUE
      WRITE(KOUT,6020)(RMOM(K),K=2,4)
      IF(KSTART.LE.NSITES)WRITE(KOUT,6030)(WORK(K,1),K=KSTART,NSITES)
  160 CONTINUE
C
      IF(NSIM.LE.0)RETURN
      IF(NPROB.GT.MAXQ)GOTO 1010
      IF(NSIM.EQ.1)GOTO 270
      IF(NMAX.GT.MAXREC)GOTO 1020
C
C         FIT KAPPA DISTRIBUTION TO REGIONAL L-MOMENTS
C
      CALL PELKAP(RMOM,RPARA,IFAIL)
      IF(IFAIL.EQ.0)GOTO 180
      CALL PELGLO(RMOM,RPARA)
      RPARA(4)=-ONE
  180 IF(KPRINT.GT.0)WRITE(KOUT,6040)(RPARA(K),K=1,4)
C
C         START THE NSIM REPETITIONS
C
      T4BAR=ZERO
      T4SD=ZERO
      DO 220 ISIM=1,NSIM
      SUM2=ZERO
      SUM3=ZERO
      SUM4=ZERO
C
C         START OF LOOP OVER SITES
C
      DO 200 I=1,NSITES
      NREC=LENGTH(I)
C
C         GET VECTOR OF UNIFORM RANDOM NUMBERS
C
      CALL DURAND(SEED,NREC,X)
C
C         TRANSFORM FROM UNIFORM TO KAPPA
C
      DO 190 J=1,NREC
      X(J)=QUAKAP(X(J),RPARA)
  190 CONTINUE
C
C         FIND L-MOMENTS OF SIMULATED DATA
C
      CALL SORT(X,NREC)
      CALL SAMLMR(X,NREC,TMOM,4,A,B,IERR)
      IF ((IERR == -1) .OR. (IERR == -2)) RETURN

      CV=TMOM(2)/TMOM(1)
      WORK(I,1)=CV
      WORK(I,2)=TMOM(3)
      WORK(I,3)=TMOM(4)
      SUM2=SUM2+NREC*CV
      SUM3=SUM3+NREC*TMOM(3)
      SUM4=SUM4+NREC*TMOM(4)
C
C         END OF LOOP OVER SITES
C
  200 CONTINUE
C
      SUM2=SUM2/SUMLEN
      SUM3=SUM3/SUMLEN
      SUM4=SUM4/SUMLEN
      T4BAR=T4BAR+SUM4
      T4SD=T4SD+SUM4**2
C
C         CALCULATE HETEROGENEITY V-STATISTICS FOR SIMULATED DATA
C
      IF(NSITES.EQ.1)GOTO 215
      V1=ZERO
      V2=ZERO
      V3=ZERO
      DO 210 I=1,NSITES
      NREC=LENGTH(I)
      TEMP2=(WORK(I,1)-SUM2)**2
      TEMP3=(WORK(I,2)-SUM3)**2
      TEMP4=(WORK(I,3)-SUM4)**2
      V1=V1+NREC*TEMP2
      V2=V2+NREC*DSQRT(TEMP2+TEMP3)
      V3=V3+NREC*DSQRT(TEMP3+TEMP4)
  210 CONTINUE
      V1=DSQRT(V1/SUMLEN)
      V2=V2/SUMLEN
      V3=V3/SUMLEN
      VBAR(1)=VBAR(1)+V1
      VBAR(2)=VBAR(2)+V2
      VBAR(3)=VBAR(3)+V3
      VSD(1)=VSD(1)+V1**2
      VSD(2)=VSD(2)+V2**2
      VSD(3)=VSD(3)+V3**2
  215 CONTINUE
C
C         END OF THE NSIM REPETITIONS
C
  220 CONTINUE
C
C         CALCULATE HETEROGENEITY V-STATISTICS FOR OBSERVED DATA
C
      IF(NSITES.EQ.1)GOTO 235
      V1=ZERO
      V2=ZERO
      V3=ZERO
      DO 225 I=1,NSITES
      NREC=LENGTH(I)
      TEMP2=(XMOM(2,I)-RMOM(2))**2
      TEMP3=(XMOM(3,I)-RMOM(3))**2
      TEMP4=(XMOM(4,I)-RMOM(4))**2
      V1=V1+NREC*TEMP2
      V2=V2+NREC*DSQRT(TEMP2+TEMP3)
      V3=V3+NREC*DSQRT(TEMP3+TEMP4)
  225 CONTINUE
      VOBS(1)=DSQRT(V1/SUMLEN)
      VOBS(2)=V2/SUMLEN
      VOBS(3)=V3/SUMLEN
C
C         CALCULATE AND PRINT HETEROGENEITY MEASURES (H STATISTICS)
C
      IF(KPRINT.GT.0)WRITE(KOUT,6050)NSIM
      DO 230 J=1,3
      VBAR(J)=VBAR(J)/NSIM
      VSD(J)=DSQRT((VSD(J)-NSIM*VBAR(J)**2)/(NSIM-ONE))
      H(J)=(VOBS(J)-VBAR(J))/VSD(J)
      IF(KPRINT.LE.0)GOTO 230
      LOOK1=BLANK
      LOOK2=BLANK
      IF(H(J).GE.HCRIT1)LOOK1=STAR
      IF(H(J).GE.HCRIT2)LOOK2=STAR
      IF(J.EQ.1)WRITE(KOUT,6060)VOBS(J),VBAR(J),VSD(J),H(J),LOOK1,LOOK2
      IF(J.EQ.2)WRITE(KOUT,6070)VOBS(J),VBAR(J),VSD(J),H(J),LOOK1,LOOK2
      IF(J.EQ.3)WRITE(KOUT,6080)VOBS(J),VBAR(J),VSD(J),H(J),LOOK1,LOOK2
  230 CONTINUE
  235 CONTINUE
C
C         FIND TAU-4 VALUES OF EACH CANDIDATE DISTRIBUTION
C
      S=RMOM(3)
      SS=S*S
      T4FIT(1)=GLOC0+SS*GLOC2
      T4FIT(2)=
     *  GEVC0+S*(GEVC1+S*(GEVC2+S*(GEVC3+S*(GEVC4+S*(GEVC5+S*GEVC6)))))
      T4FIT(3)=GNOC0+SS*(GNOC2+SS*(GNOC4+SS*(GNOC6+SS*GNOC8)))
      T4FIT(4)=PE3C0+SS*(PE3C2+SS*(PE3C4+SS*(PE3C6+SS*PE3C8)))
      T4FIT(5)=S*(GPAC1+S*(GPAC2+S*(GPAC3+S*GPAC4)))
C
C         CALCULATE GOODNESS-OF-FIT MEASURES (Z STATISTICS)
C
      T4BAR=T4BAR/NSIM
      T4SD=DSQRT((T4SD-NSIM*T4BAR**2)/(NSIM-ONE))
      DO 240 IDIST=1,5
      Z(IDIST)=(T4FIT(IDIST)+T4BAR-TWO*RMOM(4))/T4SD
  240 CONTINUE
C
C         PRINT Z STATISTICS
C
      IF(KPRINT.LE.0)GOTO 260
      WRITE(KOUT,6090)NSIM
      DO 250 IDIST=1,5
      LOOK1=BLANK
      IF(DABS(Z(IDIST)).LT.ZCRIT)LOOK1=STAR
  250 WRITE(KOUT,6100)DISTRI(IDIST),T4FIT(IDIST),Z(IDIST),LOOK1
  260 CONTINUE
C
C         FIT DISTRIBUTIONS
C
  270 CONTINUE
      CALL PELGLO(RMOM,PARA(1,1))
      CALL PELGEV(RMOM,PARA(1,2))
      CALL PELGNO(RMOM,PARA(1,3))
      CALL PELPE3(RMOM,PARA(1,4))
      CALL PELGPA(RMOM,PARA(1,5))
      CALL PELWAK(RMOM,PARA(1,6),IFAIL)
C
C         FOR SUCCESSFUL CANDIDATES AND WAKEBY, PRINT PARAMETERS ...
C
      IF(KPRINT.LE.0)GOTO 320
      IF(NSIM.EQ.1)WRITE(KOUT,6110)
      IF(NSIM.GT.1)WRITE(KOUT,6120)
      DO 280 IDIST=1,5
      IF(DABS(Z(IDIST)).LE.ZCRIT)
     *  WRITE(KOUT,6130)DISTRI(IDIST),(PARA(IPARA,IDIST),IPARA=1,3)
  280 CONTINUE
      WRITE(KOUT,6130)DISTRI(6),(PARA(IPARA,6),IPARA=1,5)
C
C         ... AND ESTIMATE AND PRINT QUANTILES
C
      IF(NPROB.EQ.0)GOTO 320
      WRITE(KOUT,6140)PROB
      DO 300 IDIST=1,5
      IF(DABS(Z(IDIST)).GT.ZCRIT)GOTO 300
      DO 290 IQ=1,NPROB
      IF(IDIST.EQ.1)Q(IQ)=QUAGLO(PROB(IQ),PARA(1,1))
      IF(IDIST.EQ.2)Q(IQ)=QUAGEV(PROB(IQ),PARA(1,2))
      IF(IDIST.EQ.3)Q(IQ)=QUAGNO(PROB(IQ),PARA(1,3))
      IF(IDIST.EQ.4)Q(IQ)=QUAPE3(PROB(IQ),PARA(1,4))
      IF(IDIST.EQ.5)Q(IQ)=QUAGPA(PROB(IQ),PARA(1,5))
  290 CONTINUE
      WRITE(KOUT,6150)DISTRI(IDIST),(Q(IQ),IQ=1,NPROB)
  300 CONTINUE
      DO 310 IQ=1,NPROB
  310 Q(IQ)=QUAWAK(PROB(IQ),PARA(1,6))
      WRITE(KOUT,6150)DISTRI(6),(Q(IQ),IQ=1,NPROB)
  320 CONTINUE
C
      RETURN
C
 1000 WRITE(KOUT,7000)'MAXNS'
      IERR = -10
      RETURN
 1010 WRITE(KOUT,7000)'MAXQ'
      IERR = -20
      RETURN
 1020 WRITE(KOUT,7000)'MAXREC'
      IERR = -30
      RETURN
 1030 WRITE(KOUT,7010)
      IERR = -40
      GOTO 140
C
 6000 FORMAT(/' SITE    N     L-CV   L-SKEW  L-KURT   D(I)')
 6010 FORMAT(2I5,2X,3F8.4,F7.2,2X,2A1)
 6020 FORMAT(/5X,'WEIGHTED MEANS',5X,6F8.4)
 6030 FORMAT(/' FLAGGED TEST VALUES'/(15F5.1))
 6040 FORMAT(/' PARAMETERS OF REGIONAL KAPPA DISTRIBUTION ',4F8.4)
 6050 FORMAT(//' ***** HETEROGENEITY MEASURES *****'/
     *  ' (NUMBER OF SIMULATIONS  =',I6,')')
 6060 FORMAT(/' OBSERVED     S.D. OF GROUP L-CV          =',F8.4/
     *        ' SIM. MEAN OF S.D. OF GROUP L-CV          =',F8.4/
     *        ' SIM. S.D. OF S.D. OF GROUP L-CV          =',F8.4/
     *        ' STANDARDIZED TEST VALUE H(1)             =',F6.2,2X,2A1)
 6070 FORMAT(/' OBSERVED AVE.  OF L-CV / L-SKEW DISTANCE =',F8.4/
     *        ' SIM. MEAN OF AVE. L-CV / L-SKEW DISTANCE =',F8.4/
     *        ' SIM. S.D. OF AVE. L-CV / L-SKEW DISTANCE =',F8.4/
     *        ' STANDARDIZED TEST VALUE H(2)             =',F6.2,2X,2A1)
 6080 FORMAT(/' OBSERVED AVE.  OF L-SKEW/L-KURT DISTANCE =',F8.4/
     *        ' SIM. MEAN OF AVE. L-SKEW/L-KURT DISTANCE =',F8.4/
     *        ' SIM. S.D. OF AVE. L-SKEW/L-KURT DISTANCE =',F8.4/
     *        ' STANDARDIZED TEST VALUE H(3)             =',F6.2,2X,2A1)
 6090 FORMAT(//' ***** GOODNESS-OF-FIT MEASURES *****'/
     *  ' (NUMBER OF SIMULATIONS  =',I6,')'/)
 6100 FORMAT(1X,A18,2X,' L-KURTOSIS=',F6.3,2X,' Z VALUE=',F6.2,1X,A1)
 6110 FORMAT(//' PARAMETER ESTIMATES'/)
 6120 FORMAT(//' PARAMETER ESTIMATES FOR DISTRIBUTIONS ACCEPTED AT THE',
     *  ' 90% LEVEL'/)
 6130 FORMAT(1X,A18,1X,5F7.3)
 6140 FORMAT(/' QUANTILE ESTIMATES'/19X,(1X,14F7.3))
 6150 FORMAT(1X,A18,(1X,14F7.3))
C
 7000 FORMAT(' *** ERROR *** ROUTINE REGTST :',
     *  ' INSUFFICIENT WORKSPACE - RECOMPILE WITH LARGER VALUE OF ',A6)
 7010 FORMAT(' *** ERROR *** ROUTINE REGTST : UNABLE TO INVERT',
     *  ' SUM-OF-SQUARES MATRIX.'/31X,'D STATISTICS NOT CALCULATED.')
C
      END
