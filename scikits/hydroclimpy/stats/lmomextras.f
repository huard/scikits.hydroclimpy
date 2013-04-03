C===================================================== QUAGEV.FOR
      DOUBLE PRECISION FUNCTION QUAGEV(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE GENERALIZED EXTREME-VALUE DISTRIBUTION
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO)GOTO 1000
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Y=-DLOG(-DLOG(F))
      IF(G.NE.ZERO)Y=(ONE-DEXP(-G*Y))/G
      QUAGEV=U+A*Y
      RETURN
C
   10 IF(F.EQ.ZERO.AND.G.LT.ZERO)GOTO 20
      IF(F.EQ.ONE .AND.G.GT.ZERO)GOTO 20
      WRITE(6,7000)
      QUAGEV=ZERO
      RETURN
   20 QUAGEV=U+A/G
      RETURN
C
 1000 WRITE(6,7010)
      QUAGEV=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAGEV :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAGEV : PARAMETERS INVALID')
      END

C===================================================== QUAGLO.FOR
      DOUBLE PRECISION FUNCTION QUAGLO(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE GENERALIZED LOGISTIC DISTRIBUTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO)GOTO 1000
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Y=DLOG(F/(ONE-F))
      IF(G.NE.ZERO)Y=(ONE-DEXP(-G*Y))/G
      QUAGLO=U+A*Y
      RETURN
C
   10 IF(F.EQ.ZERO.AND.G.LT.ZERO)GOTO 20
      IF(F.EQ.ONE .AND.G.GT.ZERO)GOTO 20
      WRITE(6,7000)
      QUAGLO=ZERO
      RETURN
   20 QUAGLO=U+A/G
      RETURN
C
 1000 WRITE(6,7010)
      QUAGLO=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAGLO :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAGLO : PARAMETERS INVALID')
      END

C===================================================== QUAGAM.FOR
      DOUBLE PRECISION FUNCTION QUAGAM(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE GAMMA DISTRIBUTION
C
C  OTHER ROUTINES USED: DERF,DLGAMA,GAMIND,QUASTN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(2)
      DATA ZERO/0D0/,P01/0.01D0/,ONE/1D0/,NINE/9D0/
C
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
C
      DATA EPS/1D-10/,MAXIT/30/
C
      QUAGAM=ZERO
      ALPHA=PARA(1)
      BETA=PARA(2)
      IF(ALPHA.LE.ZERO.OR.BETA.LE.ZERO)GOTO 1000
      IF(F.LT.ZERO.OR.F.GE.ONE)GOTO 1010
      IF(F.EQ.ZERO)RETURN
      AM1=ALPHA-ONE
      IF(AM1.NE.ZERO)GOTO 10
C
C         CASE ALPHA.EQ.1 - GAMMA IS EXPONENTIAL
C
      QUAGAM=(-DLOG(ONE-F))*BETA
      RETURN
C
C         INITIAL ESTIMATE OF ROOT OF EQUATION GAMIND(X)=F:
C         - IF ALPHA.GT.1, USE WILSON-HILFERTY APPROXIMATION IF IT'S
C           POSITIVE AND NOT TOO CLOSE TO ZERO;
C         - IF ALPHA.LT.1, OR IF W-H APPROX. ISN'T POSITIVE ENOUGH,
C           USE THE SMALL-X APPROXIMATION OF IGNORING THE EXP(-T) TERM
C           IN THE INTEGRAL DEFINING GAMIND(X)
C
   10 DLOGG=DLGAMA(ALPHA)
      IF(AM1.LE.ZERO)GOTO 20
      ROOT=ALPHA*(ONE-ONE/(NINE*ALPHA)+QUASTN(F)/DSQRT(NINE*ALPHA))**3
      IF(ROOT.GT.P01*ALPHA)GOTO 30
   20 ROOT=DEXP((DLOG(ALPHA*F)+DLOGG)/ALPHA)
   30 CONTINUE
C
C         REFINE INITIAL ESTIMATE BY NEWTON-RAPHSON ITERATION
C
      DO 40 IT=1,MAXIT
      FUNC=GAMIND(ROOT,ALPHA,DLOGG)-F
      RINC=FUNC*DEXP(DLOGG+ROOT-AM1*DLOG(ROOT))
      ROOT=ROOT-RINC
      IF(DABS(FUNC).LE.EPS)GOTO 50
   40 CONTINUE
      WRITE(6,7020)
C
C         SCALE SOLUTION
C
   50 QUAGAM=ROOT*BETA
      RETURN
C
 1000 WRITE(6,7000)
      RETURN
 1010 WRITE(6,7010)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAGAM : PARAMETERS INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAGAM :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7020 FORMAT(' ** WARNING ** ROUTINE QUAGAM :',
     *  ' ITERATION HAS NOT CONVERGED. RESULT MAY BE UNRELIABLE')
      END

C===================================================== QUAGNO.FOR
      DOUBLE PRECISION FUNCTION QUAGNO(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE GENERALIZED NORMAL DISTRIBUTION
C
C  OTHER ROUTINES USED: QUASTN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO)GOTO 1000
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Y=QUASTN(F)
      IF(G.NE.ZERO)Y=(ONE-DEXP(-G*Y))/G
      QUAGNO=U+A*Y
      RETURN
C
   10 IF(F.EQ.ZERO.AND.G.LT.ZERO)GOTO 20
      IF(F.EQ.ONE .AND.G.GT.ZERO)GOTO 20
      WRITE(6,7000)
      QUAGNO=ZERO
      RETURN
   20 QUAGNO=U+A/G
      RETURN
C
 1000 WRITE(6,7010)
      QUAGNO=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAGNO :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAGNO : PARAMETERS INVALID')
      END

C===================================================== QUAGPA.FOR
      DOUBLE PRECISION FUNCTION QUAGPA(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE GENERALIZED PARETO DISTRIBUTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO)GOTO 1000
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Y=-DLOG(ONE-F)
      IF(G.NE.ZERO)Y=(ONE-DEXP(-G*Y))/G
      QUAGPA=U+A*Y
      RETURN
C
   10 IF(F.EQ.ZERO)QUAGPA=U
      IF(F.EQ.ZERO)RETURN
      IF(F.EQ.ONE.AND.G.GT.ZERO)QUAGPA=U+A/G
      IF(F.EQ.ONE.AND.G.GT.ZERO)RETURN
      WRITE(6,7000)
      QUAGPA=ZERO
      RETURN
C
 1000 WRITE(6,7010)
      QUAGPA=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAGPA :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAGPA : PARAMETERS INVALID')
      END

C===================================================== QUAKAP.FOR
      DOUBLE PRECISION FUNCTION QUAKAP(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE KAPPA DISTRIBUTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(4)
      DATA ZERO/0D0/,ONE/1D0/
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      H=PARA(4)
      IF(A.LE.ZERO)GOTO 1000
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Y=-DLOG(F)
      IF(H.NE.ZERO)Y=(ONE-DEXP(-H*Y))/H
      Y=-DLOG(Y)
      IF(G.NE.ZERO)Y=(ONE-DEXP(-G*Y))/G
      QUAKAP=U+A*Y
      RETURN
C
   10 IF(F.EQ.ZERO)GOTO 20
      IF(F.EQ.ONE)GOTO 30
      GOTO 1010
   20 IF(H.LE.ZERO.AND.G.LT.ZERO)QUAKAP=U+A/G
      IF(H.LE.ZERO.AND.G.GE.ZERO)GOTO 1010
      IF(H.GT.ZERO.AND.G.NE.ZERO)QUAKAP=U+A/G*(ONE-H**(-G))
      IF(H.GT.ZERO.AND.G.EQ.ZERO)QUAKAP=U+A*DLOG(H)
      RETURN
   30 IF(G.LE.ZERO)GOTO 1010
      QUAKAP=U+A/G
      RETURN
C
 1000 WRITE(6,7000)
      QUAKAP=ZERO
      RETURN
 1010 WRITE(6,7010)
      QUAKAP=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAKAP : PARAMETERS INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAKAP :',
     *  ' ARGUMENT OF FUNCTION INVALID')
      END

C===================================================== QUAPE3.FOR
      DOUBLE PRECISION FUNCTION QUAPE3(F,PARA)
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
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE PEARSON TYPE 3 DISTRIBUTION
C
C  OTHER ROUTINES USED: DERF,DLGAMA,GAMIND,QUAGAM,QUASTN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3),PAR(2)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,FOUR/4D0/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
      IF(PARA(2).LE.ZERO)GOTO 1000
      GAMMA=PARA(3)
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 20
      IF(DABS(GAMMA).LT.SMALL)GOTO 10
      ALPHA=FOUR/(GAMMA*GAMMA)
      BETA=DABS(HALF*PARA(2)*GAMMA)
      PAR(1)=ALPHA
      PAR(2)=BETA
      IF(GAMMA.GT.ZERO)QUAPE3=PARA(1)-ALPHA*BETA+QUAGAM(F,PAR)
      IF(GAMMA.LT.ZERO)QUAPE3=PARA(1)+ALPHA*BETA-QUAGAM(ONE-F,PAR)
      RETURN
C
C         ZERO SKEWNESS
C
   10 QUAPE3=PARA(1)+PARA(2)*QUASTN(F)
      RETURN
C
   20 IF(F.EQ.ZERO.AND.GAMMA.GT.ZERO)GOTO 30
      IF(F.EQ.ONE .AND.GAMMA.LT.ZERO)GOTO 30
      WRITE(6,7000)
      QUAPE3=ZERO
      RETURN
   30 QUAPE3=PARA(1)-TWO*PARA(2)/GAMMA
      RETURN
C
 1000 WRITE(6,7010)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAPE3 :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAPE3 : PARAMETERS INVALID')
      END

C===================================================== QUASTN.FOR
      DOUBLE PRECISION FUNCTION QUASTN(F)
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
C*  * Fixed: WRITE(6,7000) and FORMAT statement 7000 incompatible      *
C*                                                                     *
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE STANDARD NORMAL DISTRIBUTION
C
C  BASED ON ALGORITHM AS241, APPL. STATIST. (1988) VOL.37 NO.3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/
      DATA SPLIT1/0.425D0/,SPLIT2/5D0/,CONST1/0.180625D0/,CONST2/1.6D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS
C
      DATA A0,A1,A2,A3,A4,A5,A6,A7,B1,B2,B3,B4,B5,B6,B7/
     *                                0.33871 32872 79636 661D  1,
     *  0.13314 16678 91784 377D  3,  0.19715 90950 30655 144D  4,
     *  0.13731 69376 55094 611D  5,  0.45921 95393 15498 715D  5,
     *  0.67265 77092 70087 009D  5,  0.33430 57558 35881 281D  5,
     *  0.25090 80928 73012 267D  4,  0.42313 33070 16009 113D  2,
     *  0.68718 70074 92057 908D  3,  0.53941 96021 42475 111D  4,
     *  0.21213 79430 15865 959D  5,  0.39307 89580 00927 106D  5,
     *  0.28729 08573 57219 427D  5,  0.52264 95278 85285 456D  4/
      DATA C0,C1,C2,C3,C4,C5,C6,C7,D1,D2,D3,D4,D5,D6,D7/
     *                                0.14234 37110 74968 358D  1,
     *  0.46303 37846 15654 530D  1,  0.57694 97221 46069 141D  1,
     *  0.36478 48324 76320 461D  1,  0.12704 58252 45236 838D  1,
     *  0.24178 07251 77450 612D  0,  0.22723 84498 92691 846D -1,
     *  0.77454 50142 78341 408D -3,  0.20531 91626 63775 882D  1,
     *  0.16763 84830 18380 385D  1,  0.68976 73349 85100 005D  0,
     *  0.14810 39764 27480 075D  0,  0.15198 66656 36164 572D -1,
     *  0.54759 38084 99534 495D -3,  0.10507 50071 64441 684D -8/
      DATA E0,E1,E2,E3,E4,E5,E6,E7,F1,F2,F3,F4,F5,F6,F7/
     *                                0.66579 04643 50110 378D  1,
     *  0.54637 84911 16411 437D  1,  0.17848 26539 91729 133D  1,
     *  0.29656 05718 28504 891D  0,  0.26532 18952 65761 230D -1,
     *  0.12426 60947 38807 844D -2,  0.27115 55568 74348 758D -4,
     *  0.20103 34399 29228 813D -6,  0.59983 22065 55887 938D  0,
     *  0.13692 98809 22735 805D  0,  0.14875 36129 08506 149D -1,
     *  0.78686 91311 45613 259D -3,  0.18463 18317 51005 468D -4,
     *  0.14215 11758 31644 589D -6,  0.20442 63103 38993 979D-14/
C
      Q=F-HALF
      IF(DABS(Q).GT.SPLIT1)GOTO 10
      R=CONST1-Q*Q
      QUASTN=Q*(((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+A0)
     *        /(((((((B7*R+B6)*R+B5)*R+B4)*R+B3)*R+B2)*R+B1)*R+ONE)
      RETURN
   10 R=F
      IF(Q.GE.ZERO)R=ONE-F
      IF(R.LE.ZERO)GOTO 1000
      R=DSQRT(-DLOG(R))
      IF(R.GT.SPLIT2)GOTO 20
      R=R-CONST2
      QUASTN=(((((((C7*R+C6)*R+C5)*R+C4)*R+C3)*R+C2)*R+C1)*R+C0)
     *      /(((((((D7*R+D6)*R+D5)*R+D4)*R+D3)*R+D2)*R+D1)*R+ONE)
      GOTO 30
   20 R=R-SPLIT2
      QUASTN=(((((((E7*R+E6)*R+E5)*R+E4)*R+E3)*R+E2)*R+E1)*R+E0)
     *      /(((((((F7*R+F6)*R+F5)*R+F4)*R+F3)*R+F2)*R+F1)*R+ONE)
   30 IF(Q.LT.ZERO)QUASTN=-QUASTN
      RETURN
C
 1000 WRITE(6,7000)
      QUASTN=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUASTN :',
     *  ' ARGUMENT OF FUNCTION INVALID')
      END
C===================================================== PELGEV.FOR
      SUBROUTINE PELGEV(XMOM,PARA)
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
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED EXTREME-VALUE
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
C  OTHER ROUTINES USED: DLGAMA
C
C  METHOD: FOR  -0.8 LE TAU3 LT 1,  K IS APPROXIMATED BY RATIONAL
C  FUNCTIONS AS IN DONALDSON (1996, COMMUN. STATIST. SIMUL. COMPUT.).
C  IF TAU3 IS OUTSIDE THIS RANGE, NEWTON-RAPHSON ITERATION IS USED.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
      DATA P8/0.8D0/,P97/0.97D0/

CF2PY DOUBLE PRECISION INTENT(IN),  DIMENSION(3):: XMOM
CF2PY DOUBLE PRECISION INTENT(OUT), DIMENSION(3) :: PARA

C
C         SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
C
      DATA SMALL/1D-5/,EPS/1D-6/,MAXIT/20/
C
C         EU IS EULER'S CONSTANT
C         DL2 IS LOG(2), DL3 IS LOG(3)
C
      DATA EU/0.57721566D0/,DL2/0.69314718D0/,DL3/1.0986123D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS FOR K
C
      DATA A0,A1,A2/ 0.28377530D0,-1.21096399D0,-2.50728214D0/
      DATA A3,A4   /-1.13455566D0,-0.07138022D0/
      DATA B1,B2,B3/ 2.06189696D0, 1.31912239D0, 0.25077104D0/
      DATA C1,C2,C3/ 1.59921491D0,-0.48832213D0, 0.01573152D0/
      DATA D1,D2   /-0.64363929D0, 0.08985247D0/
C
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(T3).GE.ONE)GOTO 1000
      IF(T3.LE.ZERO)GOTO 10
C
C         RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN 0 AND 1
C
      Z=ONE-T3
      G=(-ONE+Z*(C1+Z*(C2+Z*C3)))/(ONE+Z*(D1+Z*D2))
      IF(DABS(G).LT.SMALL)GOTO 50
      GOTO 40
C
C         RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN -0.8 AND 0
C
   10 G=(A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(ONE+T3*(B1+T3*(B2+T3*B3)))
      IF(T3.GE.-P8)GOTO 40
C
C         NEWTON-RAPHSON ITERATION FOR TAU3 LESS THAN -0.8
C
      IF(T3.LE.-P97)G=ONE-DLOG(ONE+T3)/DL2
      T0=(T3+THREE)*HALF
      DO 20 IT=1,MAXIT
      X2=TWO**(-G)
      X3=THREE**(-G)
      XX2=ONE-X2
      XX3=ONE-X3
      T=XX3/XX2
      DERIV=(XX2*X3*DL3-XX3*X2*DL2)/(XX2*XX2)
      GOLD=G
      G=G-(T-T0)/DERIV
      IF(DABS(G-GOLD).LE.EPS*G)GOTO 30
   20 CONTINUE
      WRITE(6,7010)
   30 CONTINUE
C
C         ESTIMATE ALPHA,XI
C
   40 PARA(3)=G
      GAM=DEXP(DLGAMA(ONE+G))
      PARA(2)=XMOM(2)*G/(GAM*(ONE-TWO**(-G)))
      PARA(1)=XMOM(1)-PARA(2)*(ONE-GAM)/G
      RETURN
C
C         ESTIMATED K EFFECTIVELY ZERO
C
   50 PARA(3)=ZERO
      PARA(2)=XMOM(2)/DL2
      PARA(1)=XMOM(1)-EU*PARA(2)
      RETURN
C
 1000 WRITE(6,7000)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID')
 7010 FORMAT(' ** WARNING ** ROUTINE PELGEV :',
     *  ' ITERATION HAS NOT CONVERGED. RESULTS MAY BE UNRELIABLE.')
      END
C===================================================== PELGLO.FOR
      SUBROUTINE PELGLO(XMOM,PARA)
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
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED LOGISTIC
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      DATA PI/3.141592653589793238D0/
C
C         SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
C         ESTIMATE K
C
      G=-XMOM(3)
      IF(XMOM(2).LE.ZERO.OR.DABS(G).GE.ONE)GOTO 1000
      IF(DABS(G).LE.SMALL)GOTO 10
C
C         ESTIMATE ALPHA, XI
C
      GG=G*PI/DSIN(G*PI)
      A=XMOM(2)/GG
      PARA(1)=XMOM(1)-A*(ONE-GG)/G
      PARA(2)=A
      PARA(3)=G
      RETURN
C
C         ESTIMATED K EFFECTIVELY ZERO
C
   10 PARA(3)=ZERO
      PARA(2)=XMOM(2)
      PARA(1)=XMOM(1)
      RETURN
C
 1000 WRITE(6,7000)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELGLO : L-MOMENTS INVALID')
      END
C===================================================== PELGNO.FOR
      SUBROUTINE PELGNO(XMOM,PARA)
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
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED NORMAL
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3. ABS(TAU3) MAY NOT EXCEED 0.95.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
C  OTHER ROUTINES USED: DERF
C
C  METHOD: RATIONAL-FUNCTION APPROXIMATION OF K IN TERMS OF TAU-3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/
      DATA P95/0.95D0/
      DATA ROOTPI/1.772453850905516027D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATION
C         A0 IS 0.5*SQRT(3/PI)
C
      DATA A0,A1,A2,A3/
     *  0.20466534D+01,-0.36544371D+01,0.18396733D+01,-0.20360244D+00/
      DATA B1,B2,B3/-0.20182173D+01,0.12420401D+01,-0.21741801D+00/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-8/
C
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO.OR.DABS(T3).GE.ONE)GOTO 1000
      IF(DABS(T3).GE.P95)GOTO 1010
      IF(DABS(T3).LE.SMALL)GOTO 30
C
      TT=T3*T3
      G=-T3*(A0+TT*(A1+TT*(A2+TT*A3)))/(ONE+TT*(B1+TT*(B2+TT*B3)))
      E=DEXP(HALF*G*G)
      A=XMOM(2)*G/(E*DERF(HALF*G))
      U=XMOM(1)+A*(E-ONE)/G
      PARA(1)=U
      PARA(2)=A
      PARA(3)=G
      RETURN
C
   30 PARA(1)=XMOM(1)
      PARA(2)=XMOM(2)*ROOTPI
      PARA(3)=ZERO
      RETURN
C
 1000 WRITE(6,7000)
      RETURN
 1010 WRITE(6,7010)
      PARA(1)=ZERO
      PARA(2)=-ONE
      PARA(3)=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELGNO : L-MOMENTS INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE PELGNO :',
     *  ' TAU-3 TOO LARGE FOR ROUTINE')
      END
C===================================================== PELGPA.FOR
      SUBROUTINE PELGPA(XMOM,PARA)
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
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR  THE GENERALIZED PARETO
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
C
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(T3).GE.ONE)GOTO 1000
      G=(ONE-THREE*T3)/(ONE+T3)
      PARA(3)=G
      PARA(2)=(ONE+G)*(TWO+G)*XMOM(2)
      PARA(1)=XMOM(1)-PARA(2)/(ONE+G)
      RETURN
C
 1000 WRITE(6,7000)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELGPA : L-MOMENTS INVALID')
      END
C===================================================== PELPE3.FOR
      SUBROUTINE PELPE3(XMOM,PARA)
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
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE PEARSON TYPE 3 DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2 AND TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER MU, SIGMA, GAMMA (MEAN, S.D., SKEWNESS).
C
C  OTHER ROUTINES USED: DLGAMA
C
C  METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA, THE SHAPE
C  PARAMETER OF THE GAMMA DISTRIBUTION, AS A FUNCTION OF TAU-3.
C  RELATIVE ACCURACY OF THE APPROXIMATION IS BETTER THAN 3E-5.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,THIRD/0.33333333D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
C         CONSTANTS USED IN MINIMAX APPROXIMATIONS
C
      DATA C1,C2,C3/ 0.2906D0,  0.1882D0,  0.0442D0/
      DATA D1,D2,D3/ 0.36067D0,-0.59567D0, 0.25361D0/
      DATA D4,D5,D6/-2.78861D0, 2.56096D0,-0.77045D0/
      DATA PI3,ROOTPI/9.4247780D0,1.7724539D0/
C
      T3=DABS(XMOM(3))
      IF(XMOM(2).LE.ZERO.OR.T3.GE.ONE)GOTO 1000
      IF(T3.LE.SMALL)GOTO 100
      IF(T3.GE.THIRD)GOTO 10
      T=PI3*T3*T3
      ALPHA=(ONE+C1*T)/(T*(ONE+T*(C2+T*C3)))
      GOTO 20
   10 CONTINUE
      T=ONE-T3
      ALPHA=T*(D1+T*(D2+T*D3))/(ONE+T*(D4+T*(D5+T*D6)))
   20 CONTINUE
      RTALPH=DSQRT(ALPHA)
      BETA=ROOTPI*XMOM(2)*DEXP(DLGAMA(ALPHA)-DLGAMA(ALPHA+HALF))
      PARA(1)=XMOM(1)
      PARA(2)=BETA*RTALPH
      PARA(3)=TWO/RTALPH
      IF(XMOM(3).LT.ZERO)PARA(3)=-PARA(3)
      RETURN
C
C         ZERO SKEWNESS
C
  100 CONTINUE
      PARA(1)=XMOM(1)
      PARA(2)=XMOM(2)*ROOTPI
      PARA(3)=ZERO
      RETURN
C
 1000 WRITE(6,7000)
      DO 1010 I=1,3
 1010 PARA(I)=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELPE3 : L-MOMENTS INVALID')
      END
C===================================================== PELWAK.FOR
      SUBROUTINE PELWAK(XMOM,PARA,IFAIL)
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
C*  VERSION 3.04  JULY 2005                                            *
C*  * Minor bug fix in test for validity of L-moments.                 *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE WAKEBY DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 5. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3, TAU-4, TAU-5.
C  PARA   *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, BETA, GAMMA, DELTA.
C  IFAIL  *OUTPUT* FAIL FLAG. ON EXIT, IT IS SET AS FOLLOWS.
C                  0 SUCCESSFUL EXIT
C                  1 ESTIMATES COULD ONLY BE OBTAINED BY SETTING XI=0
C                  2 ESTIMATES COULD ONLY BE OBTAINED BY FITTING A
C                    GENERALIZED PARETO DISTRIBUTION
C                  3 L-MOMENTS INVALID
C
C  PROCEDURE:
C  1. LOOK FOR A SOLUTION WITH XI UNCONSTRAINED;
C  2. IF NONE FOUND, LOOK FOR A SOLUTION WITH XI=0;
C  3. IF NONE FOUND, FIT A GENERALIZED PARETO DISTRIBUTION TO THE
C     FIRST 3 L-MOMENTS.
C  ESTIMATES ARE CALCULATED USING THE FORMULAS GIVEN BY GREENWOOD ET AL.
C  (1979, WATER RESOUR. RES., TABLE 5), BUT EXPRESSED IN TERMS OF
C  L-MOMENTS RATHER THAN PROBABILITY WEIGHTED MOMENTS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(5),PARA(5)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/
      DATA X2/2D0/,X3/3D0/,X4/4D0/,X5/5D0/,X7/7D0/,X8/8D0/,X9/9D0/,
     *  X10/10D0/,X11/11D0/,X16/16D0/,X25/25D0/,X29/29D0/,X32/32D0/,
     *  X35/35D0/,X85/85D0/,X125/125D0/,X203/203D0/
C
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(XMOM(3)).GE.ONE)GOTO 1000
      IF(DABS(XMOM(4)).GE.ONE)GOTO 1000
      IF(DABS(XMOM(5)).GE.ONE)GOTO 1000
      IFAIL=0
C
C         CALCULATE THE L-MOMENTS (LAMBDA'S)
C
      ALAM1=XMOM(1)
      ALAM2=XMOM(2)
      ALAM3=XMOM(3)*ALAM2
      ALAM4=XMOM(4)*ALAM2
      ALAM5=XMOM(5)*ALAM2
C
C         ESTIMATE N1,N2,N3,C1,C2,C3 WHEN XI.NE.0
C
      XN1= X3*ALAM2-X25*ALAM3 +X32*ALAM4
      XN2=-X3*ALAM2 +X5*ALAM3  +X8*ALAM4
      XN3= X3*ALAM2 +X5*ALAM3  +X2*ALAM4
      XC1= X7*ALAM2-X85*ALAM3+X203*ALAM4-X125*ALAM5
      XC2=-X7*ALAM2+X25*ALAM3  +X7*ALAM4 -X25*ALAM5
      XC3= X7*ALAM2 +X5*ALAM3  -X7*ALAM4  -X5*ALAM5
C
C         ESTIMATE B AND D
C
      XA=XN2*XC3-XC2*XN3
      XB=XN1*XC3-XC1*XN3
      XC=XN1*XC2-XC1*XN2
      DISC=XB*XB-FOUR*XA*XC
      IF(DISC.LT.ZERO)GOTO 10
      DISC=DSQRT(DISC)
      ROOT1=HALF*(-XB+DISC)/XA
      ROOT2=HALF*(-XB-DISC)/XA
      B= DMAX1(ROOT1,ROOT2)
      D=-DMIN1(ROOT1,ROOT2)
      IF(D.GE.ONE)GOTO 10
C
C         ESTIMATE A, C AND XI
C
      A=(ONE+B)*(TWO+B)*(THREE+B)/
     *  (FOUR*(B+D))*((ONE+D)*ALAM2-(THREE-D)*ALAM3)
      C=-(ONE-D)*(TWO-D)*(THREE-D)/
     *  (FOUR*(B+D))*((ONE-B)*ALAM2-(THREE+B)*ALAM3)
      XI=ALAM1-A/(ONE+B)-C/(ONE-D)
C
C         CHECK FOR VALID PARAMETERS
C
      IF(C.GE.ZERO.AND.A+C.GE.ZERO)GOTO 30
C
C         CAN'T FIND VALID ESTIMATES FOR XI UNRESTRICTED, SO TRY XI=0
C
C         ESTIMATE B AND D FOR XI=0
C
   10 IFAIL=1
      XI=ZERO
      ZN1=X4*ALAM1-X11*ALAM2+X9*ALAM3
      ZN2=-ALAM2+X3*ALAM3
      ZN3=ALAM2+ALAM3
      ZC1=X10*ALAM1-X29*ALAM2+X35*ALAM3-X16*ALAM4
      ZC2=-ALAM2+X5*ALAM3-X4*ALAM4
      ZC3=ALAM2-ALAM4
      ZA=ZN2*ZC3-ZC2*ZN3
      ZB=ZN1*ZC3-ZC1*ZN3
      ZC=ZN1*ZC2-ZC1*ZN2
      DISC=ZB*ZB-FOUR*ZA*ZC
      IF(DISC.LT.ZERO)GOTO 20
      DISC=DSQRT(DISC)
      ROOT1=HALF*(-ZB+DISC)/ZA
      ROOT2=HALF*(-ZB-DISC)/ZA
      B= DMAX1(ROOT1,ROOT2)
      D=-DMIN1(ROOT1,ROOT2)
      IF(D.GE.ONE)GOTO 20
C
C         ESTIMATE A AND C
C
      A= (ONE+B)*(TWO+B)/(B+D)*(ALAM1-(TWO-D)*ALAM2)
      C=-(ONE-D)*(TWO-D)/(B+D)*(ALAM1-(TWO+B)*ALAM2)
      IF(C.GE.ZERO.AND.A+C.GE.ZERO)GOTO 30
C
C         CAN'T FIND VALID ESTIMATES EVEN WITH XI=0 -
C         FIT GENERALIZED PARETO DISTRIBUTION INSTEAD
C
   20 IFAIL=2
      D=-(ONE-THREE*XMOM(3))/(ONE+XMOM(3))
      C=(ONE-D)*(TWO-D)*XMOM(2)
      B=ZERO
      A=ZERO
      XI=XMOM(1)-C/(ONE-D)
      IF(D.GT.ZERO)GOTO 30
      A=C
      B=-D
      C=ZERO
      D=ZERO
C
C         COPY RESULTS INTO ARRAY PARA
C
   30 PARA(1)=XI
      PARA(2)=A
      PARA(3)=B
      PARA(4)=C
      PARA(5)=D
      RETURN
C
 1000 IFAIL=3
      DO 1010 I=1,5
 1010 PARA(I)=ZERO
      END
C===================================================== GAMIND.FOR
      DOUBLE PRECISION FUNCTION GAMIND(X,ALPHA,G)
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
C***********************************************************************
C
C  THE INCOMPLETE GAMMA INTEGRAL
C
C  BASED ON ALGORITHM AS239, APPL. STATIST. (1988) VOL.37 NO.3
C
C  PARAMETERS OF ROUTINE:
C  X      * INPUT* ARGUMENT OF FUNCTION (UPPER LIMIT OF INTEGRATION)
C  ALPHA  * INPUT* SHAPE PARAMETER
C  G      * INPUT* LOG(GAMMA(ALPHA)). MUST BE SUPPLIED BY THE PROGRAM,
C                  E.G. AS DLGAMA(ALPHA).
C
C  OTHER ROUTINES USED: DERF
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,X13/13D0/,
     *  X36/36D0/,X42/42D0/,X119/119D0/,X1620/1620D0/,X38880/38880D0/,
     *  RTHALF/0.70710 67811 86547 524D0/
C
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF THE SERIES AND
C           CONTINUED-FRACTION EXPANSIONS.
C         OFL IS A LARGE NUMBER, USED TO RESCALE THE CONTINUED FRACTION.
C         UFL IS SUCH THAT EXP(UFL) IS JUST .GT. ZERO.
C         AHILL CONTROLS THE SWITCH TO HILL'S APPROXIMATION.
C
      DATA EPS/1D-12/,MAXIT/100000/,OFL/1D30/,UFL/-180D0/,AHILL/1D4/
      GAMIND=ZERO
      IF(ALPHA.LE.ZERO)GOTO 1000
      IF(X.LT.ZERO)GOTO 1010
      IF(X.EQ.ZERO)RETURN
C
      IF(ALPHA.GT.AHILL)GOTO 100
      IF(X.GT.ONE.AND.X.GE.ALPHA)GOTO 50
C
C         SERIES EXPANSION
C
      SUM=ONE
      TERM=ONE
      A=ALPHA
      DO 10 IT=1,MAXIT
      A=A+ONE
      TERM=TERM*X/A
      SUM=SUM+TERM
      IF(TERM.LE.EPS)GOTO 20
   10 CONTINUE
      WRITE(6,7020)
   20 ARG=ALPHA*DLOG(X)-X-G+DLOG(SUM/ALPHA)
      GAMIND=ZERO
      IF(ARG.GE.UFL)GAMIND=DEXP(ARG)
      RETURN
C
C         CONTINUED-FRACTION EXPANSION
C
   50 CONTINUE
      A=ONE-ALPHA
      B=A+X+ONE
      TERM=ZERO
      PN1=ONE
      PN2=X
      PN3=X+ONE
      PN4=X*B
      RATIO=PN3/PN4
      DO 70 IT=1,MAXIT
      A=A+ONE
      B=B+TWO
      TERM=TERM+ONE
      AN=A*TERM
      PN5=B*PN3-AN*PN1
      PN6=B*PN4-AN*PN2
      IF(PN6.EQ.ZERO)GOTO 60
      RN=PN5/PN6
      DIFF=DABS(RATIO-RN)
      IF(DIFF.LE.EPS.AND.DIFF.LE.EPS*RN)GOTO 80
      RATIO=RN
   60 PN1=PN3
      PN2=PN4
      PN3=PN5
      PN4=PN6
      IF(DABS(PN5).LT.OFL)GOTO 70
      PN1=PN1/OFL
      PN2=PN2/OFL
      PN3=PN3/OFL
      PN4=PN4/OFL
   70 CONTINUE
      WRITE(6,7020)
   80 ARG=ALPHA*DLOG(X)-X-G+DLOG(RATIO)
      GAMIND=ONE
      IF(ARG.GE.UFL)GAMIND=ONE-DEXP(ARG)
      RETURN
C
C         ALPHA IS LARGE: USE HILL'S APPROXIMATION (N.L. JOHNSON AND
C         S. KOTZ, 1970, 'CONTINUOUS UNIVARIATE DISTRIBUTIONS 1', P.180)
C
C         THE 'DO 110' LOOP CALCULATES 2*(X-ALPHA-ALPHA*DLOG(X/ALPHA)),
C         USING POWER-SERIES EXPANSION TO AVOID ROUNDING ERROR
C
  100 CONTINUE
      R=ONE/DSQRT(ALPHA)
      Z=(X-ALPHA)*R
      TERM=Z*Z
      SUM=HALF*TERM
      DO 110 I=1,12
      TERM=-TERM*Z*R
      SUM=SUM+TERM/(I+TWO)
      IF(DABS(TERM).LT.EPS)GOTO 120
  110 CONTINUE
  120 WW=TWO*SUM
      W=DSQRT(WW)
      IF(X.LT.ALPHA)W=-W
      H1=ONE/THREE
      H2=-W/X36
      H3=(-WW+X13)/X1620
      H4=(X42*WW+X119)*W/X38880
      Z=(((H4*R+H3)*R+H2)*R+H1)*R+W
      GAMIND=HALF+HALF*DERF(Z*RTHALF)
      RETURN
C
 1000 WRITE(6,7000)ALPHA
      RETURN
 1010 WRITE(6,7010)X
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE GAMIND :',
     *  ' SHAPE PARAMETER OUT OF RANGE :',D16.8)
 7010 FORMAT(' *** ERROR *** ROUTINE GAMIND :',
     *  ' ARGUMENT OF FUNCTION OUT OF RANGE :',D16.8)
 7020 FORMAT(' ** WARNING ** ROUTINE GAMIND :',
     *  ' ITERATION HAS NOT CONVERGED. RESULT MAY BE UNRELIABLE.')
      END
C===================================================== SORT.FOR
      SUBROUTINE SORT(X,N)
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
C***********************************************************************
C
C  SORTS THE ARRAY X INTO ASCENDING ORDER
C
C  PARAMETERS OF ROUTINE:
C  X      *IN/OUT* ARRAY OF LENGTH N. CONTAINS THE NUMBERS TO BE SORTED.
C                  ON EXIT, CONTAINS THE SORTED NUMBERS.
C  N      * INPUT* NUMBER OF ELEMENTS TO BE SORTED
C
C  METHOD USED IS SHELL SORT WITH SEQUENCE OF INCREMENTS AS IN
C  D.F.KNUTH (1969) 'THE ART OF COMPUTER PROGRAMMING', VOL.3, P.95
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N)
      IF(N.LE.1)RETURN
      J=4
      DO 10 I=1,100
      J=3*J+1
      IF(J.GE.N)GOTO 20
   10 CONTINUE
   20 CONTINUE
      M=(J/3)
      DO 60 MM=1,100
      M=M/3
      IF(M.EQ.0)RETURN
      DO 50 I=M+1,N
      TEST=X(I)
      J=I
      DO 30 JJ=1,100
      J=J-M
      IF(J.LE.0)GOTO 40
      IF(TEST.GE.X(J))GOTO 40
   30 X(J+M)=X(J)
   40 CONTINUE
   50 X(J+M)=TEST
   60 CONTINUE
      END

C===================================================== DURAND.FOR
      SUBROUTINE DURAND(SEED,N,X)
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
C***********************************************************************
C
C  PSEUDO RANDOM NUMBER GENERATOR
C
C  PARAMETERS OF ROUTINE:
C  SEED   *IN/OUT* SEED FOR RANDOM NUMBER GENERATOR. SHOULD BE A WHOLE
C                  NUMBER IN THE RANGE 2D0 TO 2147483647D0.
C  N      * INPUT* NUMBER OF NUMBERS TO BE GENERATED
C  X      *OUTPUT* ARRAY OF LENGTH N. ON EXIT, CONTAINS RANDOM NUMBERS.
C
C  METHOD USED: MULTIPLICATIVE CONGRUENTIAL GENERATOR WITH BASE 2**31-1
C  AND MULTIPLIER 7**5 (P.A.W. LEWIS ET AL., 1969, IBM SYSTEMS JOURNAL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N)
      DATA AMULT/16807D0/
      DATA BASE,RBASE/2147483647D0,4.65661287524579692D-10/

      DO 10 I=1,N
      SEED=DMOD(SEED*AMULT,BASE)
      X(I)=SEED*RBASE
   10 CONTINUE
      RETURN
      END
