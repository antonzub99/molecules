      PROGRAM OEN

c** this code first runs the hermann-skillman routines SLATER
c** for discrete configuration and screening potential
c** initial data for SLATER in file hesk.inp :

c** first card:   5i3, f10.6: iz, ion, ns, nf, iut, eps
c**  iz = nuclear charge,
c** ion = degree of ionization (1 = neutral),
c**  ns = number of subshells,
c**  nf = free parameter,
c** iut = print options,
c** eps = accuracy

c** next ns cards: 3i3, F11.3 : n l k Ry
c** are subshell data:
c**  n = principal quantum number,
c**  l = orbital momentum,
c**  k = number of electrons,
c** Ry = initial estimate of ionization energy in Ry

      open(1,file='hesk3.txt',status='old')
      open(2,file='hesk3a.out')

      CALL SLATER

c** SLATER saves in arrays radial mesh, discrete functions and
c** screening potential
c** CONTIN calculates radial continuum functions and phases for
c** energy mesh from the

c** next card: 3f8.4: lc, emin, del, emax
c** lc   = orbital momentum of electron
c** emin = minimum energy of electron in eV,
c** del  = energy step in eV,
c** emax = maximal energy in eV

       CALL CONTIN

       close(1)
       close(2)

       END

******************************************************

      SUBROUTINE SLATER
C *** CALCULATION OF DISCRETE STATES IN THE H-S MODEL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*6 LTEXT
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,RMAX,Z
      COMMON/KE/K(4,25)/EK/E(25)
      COMMON/R/R(500)
     */V/V(510)
     */RWFC/FC(510,18)
     */LOCAL/RESERV(1020)
     */LTEXT/LTEXT(8)
      LTEXT(1)='S     '
      LTEXT(2)='P     '
      LTEXT(3)='D     '
      LTEXT(4)='F     '
      LTEXT(5)='G     '
      LTEXT(6)='H     '
      LTEXT(7)='I     '
      LTEXT(8)='J     '
      CALL INITIA
      CALL SCPOT
      CALL OUTPUT
      RETURN
      END
C_________________________________________________________________
C
      SUBROUTINE INITIA
C *** INPUT DATA, ARRAY OF COORDINATE, ZERO-ORDER POTENTIAL
      CALL INDAT
      CALL VARI
      CALL POT1
      CALL FUNC1
      RETURN
      END
C________________________________________________________________
C
      SUBROUTINE INDAT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C *** INPUT DATA
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,RMAX,Z
      COMMON/KE/K(4,25)/EK/E(25)
C *** K( , ) - DATA FOR CONFIGURATION : N ,L , Z(NL), TYPE OF IUTPUT
C *** E( ) - ENERGIES OF SUBSHELLS (RY)
      T=20.d0
      BETA=.5d0
      H=.1d0
      NN=510
      N1=500
      RMAX=40.0d0
      WRITE(2,10)
10    format(27H  IZ ION NS NF IUT      EPS)
      READ(1,20) IZ,ION,NS,NF,IUT,EPS
20    FORMAT(5I3,F10.6)
      Z=IZ
      II=0
      WRITE(2,30)
30    FORMAT(18H   N  L  K   E(RY))
            DO 6 J=1,NS
            READ(1,40) (K(I,J),I=1,3), E(J)
40          FORMAT(3I3,F11.3)
            II=II+K(3,J)
6           CONTINUE
      IF (IZ-II+1.NE.ION) GOTO 7
      R1=1.d0/(Z*T)
      ALFA=((N1-1)*H-BETA*DLOG(RMAX/R1))/(RMAX-R1)
      IF (ALFA.GE.0.d0) GO TO 9
      ALFA=0.
      BETA=(N1-1.d0)*H/DLOG(RMAX/R1)
9     CONTINUE
      IF (IUT.GE.1) WRITE (2,50) Z,H,T,R1,RMAX
50    FORMAT(//,3H Z=,F4.0,2X,2HH=,F6.4,4H  T=,F5.1,5H  R1=,F9.6,
     *7H  RMAX=,F10.5)
      IF (IUT.GE.1) WRITE(2,60) ALFA,BETA,NN,N1,EPS
60    FORMAT(7H  ALFA=,F8.6,7H  BETA=,F6.3,4H NN=,I3,4H N1=,I3,
     *5H EPS=,F8.6)
      RETURN
7     WRITE (2,80)
80    FORMAT(19H INPUT DATA MISTAKE)
      STOP
      RETURN
      END
C_____________________________________________________________
C
      SUBROUTINE VARI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C *** ARRAY OF COORDINATE R(RHO); RHO(J)=RHO1+H*(J-1)=
C                                       =ALFA*R + BETA*DLOG(R)

      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,RMAX,Z
     */R/R(500)
      RT=DBLE(RMAX)
      DA=DBLE(ALFA)
      DB=DBLE(BETA)
      DH=DBLE(H)
      DRM=DBLE(RMAX)
      RHO2=DA*DRM+DB*DLOG(DRM)
            DO 1 J=1,N1
            I=N1-J+1
            RHOT=RHO2 -DH*(J-1)
2           A=(DA*RT+DB*DLOG(RT)-RHOT)/(DA*RT+DB)
            RT=RT*(1.0D+0-A)
            IF (ABS(A).GT.1.0D-7) GO TO 2
            R(I)=RT
1           CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE POT1
C *** ZERO-ORDER POTENTIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,RMAX,Z
     */R/R(500)
     */V/V(510)
      A0=1.05d0
      A1=9.0d0
      B0=0.1837d0
      B1=0.04d0
      C=3.1d0
      C1=-2.d0*C*(Z*Z*3.d0/3.1416d0)**(1.d0/3.d0)
      TM=(Z**(1.d0/3.d0))/0.8853d0
            DO 1 J=1,N1
            RJ=R(J)
            AJ=ALFA*RJ
            BJ= AJ+BETA
            X=TM*RJ
            VJ=-2.d0*(ION+(Z-ION)*EXP(-B0*X)/(1.d0+A0*X))/RJ+
     *      C1*EXP(-B1*X)/(1.d0+A1*X)
            UJ=-2.d0*ION/RJ
            IF (UJ.LT.VJ) VJ=UJ
            V(J)=VJ+BETA*(AJ+0.25d0*BETA)/(BJ*BJ*RJ*RJ)
1           CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE FUNC1
C *** PRERARING OF INITIAL DATA FOR CALCULATION OF WAVE FUNCTIONS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
      COMMON/KE/K(4,25)/EK/E(25)
     */LOCAL/ A(510),RESERV(510)
      II=0
            DO 36 NR=1,NS
                  DO 1 I=1,NN
1                 A(I)=0.d0
            A(N1+1)=Z-II
            II=II+K(3,NR)
            A(N1+2)=K(1,NR)
            A(N1+3)=K(2,NR)
            A(N1+4)=K(3,NR)
            A(N1+7)=HNOR(A)
            A(N1+8)=E(NR)
            A(N1+9)=N1
            CALL WD(NR,A)
36          CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE WD(NR,A)
C *** WRITING OF ARRAY A TO THE PLACE NR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510)
      COMMON/RWFC/FC(510,18)
      NN=510
            DO 8 J=1,NN
8           FC(J,NR)=A(J)
      RETURN
      END
C______________________________________
C
      SUBROUTINE RD(NR,A)
C *** READING OF ARRAY FROM PLACE NR TO THE A
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510)
      COMMON/RWFC/FC(510,18)
      NN=510
      DO 10 I=1,NN
10    A(I)=FC(I,NR)
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE SCPOT
C *** CALCULATING SELF-CONSISTENT POTENTIAL AND WAVE FUNCTIONS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*6 LTEXT
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
      COMMON/KE/K(4,25)/EK/E(25)
     */R/R(500)
     */V/V(510)
     */LTEXT/LTEXT(8)
      DIMENSION A(510),F(510),C(510)
      H3=H/3.d0
      M1=N1+1
      M2=N1+2
      M4=N1+4
      M5=N1+5
      M7=N1+7
      M8=N1+8
      M9=N1+9
          DO 19 IT=1,30
          ITER=IT
          IF (IUT.GE.2) WRITE(2,100) ITER
100       FORMAT(22H NUMBER OF ITERATION =,I3)
          S1=0.d0
          S2=0.d0
          KRMAX=0
              DO 1 J=1,N1
1             C(J)=0.d0
              DO 6 NR=1,NS
C *** CALCULATING OF WAVE FUNCTIONS
              CALL RD(NR,A)
              GO TO 5
2             ITER=-2
5             KA9=A(N1+9)+0.01d0
              KR=MIN0(N1,KA9+20)
              A(N1+9)=KR
              AKV=A(N1+3)*(A(N1+3)+1.0d0)
                  DO 3 J=1,KR
                  RJ=R(J)
3                 F(J)=AKV/(RJ*RJ)+V(J)
              CALL EIGENF(A,F)
              IF(A(N1+9).GT.KA9+15) GO TO 2
              ITER=IT
              KR=A(N1+9)+0.01d0
              KRMAX=MAX0(KR,KRMAX)
              A(N1+1)=(3.d0*A(N1+2)*A(N1+2)-A(N1+3)*(A(N1+3)+1.d0))/
     *        (2.d0*RKIN(1,A,A))
              E(NR)=A(N1+8)
              CALL WD(NR,A)
              IF (IUT.LT.2) GO TO 105
              L1=A(N1+3)+1.01d0
              WRITE(2,103) A(M2),LTEXT(L1),A(M4),A(M8),A(M1)
              WRITE(2,203) A(M7),A(M5),A(M9),R(KR)
103           FORMAT(F4.0,A1,2X,F4.0,2X,2HE=,E11.3,2X,3HZE=,F8.3)
203           FORMAT(2X,3HA7=,F9.4,2X,3HJC=,F4.0,2X,3HKR=,F4.0,
     *        2X,6HR(KR)=,F7.2)
C *** CALCULATING OF CHARGE DENSITY
105           ZNL=K(3,NR)
              AL=A(N1+3)
              TAR=A(N1+1)/(AL+1.d0)/(AL+2.d0)
              S1=S1+ZNL*A(1)*A(1)*R(1)/(ALFA+BETA/R(1))*(1+TAR*R(1))/
     *        (2*AL+3)
              S2=S2+ZNL*A(2)*A(2)*R(2)/(ALFA+BETA/R(2))*(1+TAR*R(2))/
     *        (2*AL+3)
                  DO 4 J=1,KR
4                 C(J)=C(J)+ZNL*A(J)*A(J)/(ALFA+BETA/R(J))
6             CONTINUE
C *** CALCULATING OF EFFECTIVE ELECTRON CHARGE Z(R)
              DO 7 J=1,N1
7             F(J)=0.d0
          F(1)=S1
          F(2)=S2
          P1=C(1)/(ALFA+BETA/R(1))
          P2=C(2)/(ALFA+BETA/R(2))
              DO 8 J=3,KRMAX
              P3=C(J)/(ALFA+BETA/R(J))
              F(J)=F(J-2)+H3*(P1+4.d0*P2+P3)
              P1=P2
              P2=P3
8             CONTINUE
C *** CALCULATING OF EFFECTIVE ELECTRON CHARGE Y(R)
          P1=F(KRMAX)/R(KRMAX)/(ALFA*R(KRMAX)+BETA)
          P2=F(KRMAX-1)/R(KRMAX-1)/(ALFA*R(KRMAX-1)+BETA)
          S1=F(KRMAX)/R(KRMAX)
          S2=F(KRMAX-1)/R(KRMAX-1)
              DO 9 I=3,KRMAX
              J=KRMAX+1-I
              RJ=R(J)
              P3=F(J)/RJ/(ALFA*RJ+BETA)
              S3=S1+H3*(P1+4.d0*P2+P3)
              F(J)=S3*RJ
              S1=S2
              S2=S3
              P1=P2
              P2=P3
9             CONTINUE
C *** ACCOUNT FOR NUCLEUS AND EXCHANGE SLATER POTENTIALS
          CH=6.d0*(3.d0/(32.d0*3.141593d0**2))**0.3333333d0
              DO 11 J=1,KRMAX
              RJ=R(J)
              AJ=ALFA*RJ
              BJ=AJ+BETA
              F(J)=-2.d0*(Z-F(J))/RJ-CH*(C(J)/(RJ*RJ))**0.3333333d0
              IF(F(J).GT.-2.d0*ION/RJ) GO TO 13
              F(J)=F(J)+BETA*(AJ+0.25d0*BETA)/(BJ*BJ*RJ*RJ)
11            CONTINUE
          J=KRMAX+1
13        V(N1+1)=J
              DO 14 JJ=J,N1
              RJ=R(JJ)
              AJ=ALFA*RJ
              BJ=AJ+BETA
14            F(JJ)=-2*ION/RJ+BETA*(AJ+0.25d0*BETA)/(BJ*BJ*RJ*RJ)
C *** CHECKING OF SELFCONSISTENCY
16        TEST = 0.d0
              DO 12 J=1,KRMAX
              DV=ABS(V(J)-F(J))*R(J)
              IF(TEST.LT.DV) TEST=DV
12            CONTINUE
          IF(IUT.GE.2) WRITE(2,107) TEST
107       FORMAT(5HTEST=,d12.3,//)
          IF(TEST-EPS) 20,20,15
15        CALL PRATT(V,F)
19        CONTINUE
20    CONTINUE
      IF(KRMAX.EQ.N1) RETURN
      KR1=KRMAX+1
          DO 21 J=KR1,N1
          RJ=R(J)
          AJ=ALFA*RJ
          BJ=AJ+BETA
          V(J)=-2.d0*ION/RJ+BETA*(AJ+0.25d0*BETA)/(BJ*BJ*RJ*RJ)
21        CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE EIGENF(A,F)
C *** CALCULATION ON DISCRETE STATE WAVE FUNCTION FROM F ARRAY
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510),F(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */LOCAL/C(510),RESERV(510)
      KR=A(N1+9)+0.01d0
      T2=H*H/12.d0
      DU=0.1d0*EPS*ABS(A(N1+8))
      IF (ITER-1) 4,1,2
1     CALL ZERF(A,F)
      KR=MIN1(KR+0.01d0,A(N1+9)+180.01d0)
      A(N1+9)=KR
      GO TO 5
2     CALL PERTUR(A)
C *** CORRECTION TO ENERGY BY PERTURBATION THEORY
4     CALL P12(A,F)
5           DO 10 KENERG=1,15
C *** ITERATIONS BY THE HARTREE METHOD
            EK=A(N1+8)
                  DO 6 J=1,KR
                  BJ=ALFA+BETA/R(J)
                  C(J)=1.d0-T2*(F(J)-EK)/(BJ*BJ)
                  IF (C(J) .GT. 1.d0) JC=J
6                 CONTINUE
            IF(JC.GT.KR-3) JC=KR-3
            CALL START(A,F)
            CALL NMRV1(A,DP,C,JC)
            IF(JC.EQ.1001) GO TO 1
            A(N1+5)=JC
            ANORMA=1./DSQRT(RKIN(0,A,A))
            A(N1+7)=ANORMA*A(N1+7)
            M1=N1+1
            M10=N1+10
            IF(IUT.GE.3) WRITE(2, 100) (A(J),J=M1,M10)
100         FORMAT(8HEIGENF: ,F7.3,4F5.0,F8.3,2d11.4,F5.0,F8.3)
            DE=DP*A(JC)*ANORMA*ANORMA/H
            IF(ABS(DE).LT.DU) GO TO 20
            IF (DE.GT.-0.4d0*EK) DE=-0.2d0*EK
            A(N1+8)=A(N1+8)+DE
10          CONTINUE
20    CONTINUE
16          DO 18 J=1,N1
            I=N1+1-J
            IF(A(I).EQ.0) GO TO 18
            IF(ABS(A(I)).GE.EPS) GO TO 19
            A(I)=0.
18          CONTINUE
19    KR=I
      A(N1+9)=KR
      ANORMA=1./DSQRT(RKIN(0,A,A))
      A(N1+7)=ANORMA*A(N1+7)
23          DO 21 J=1,KR
21          A(J)=ANORMA*A(J)
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE P12(A,F)
C *** PARAMETERS FOR CALCULATION IN 2 POINTS NEAR ZERO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510),F(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */ST/ST(4)
      DIMENSION RL(2),U(2)
      L=A(N1+3)+0.01d0
      T1=Z/(L+1.d0)
      T2=2.d0*Z*T1
      T3=1.d0/(4.d0*L+6.d0)
      T4=1.d0/(6.d0*L+12.d0)
      T5=(6.d0*L+8.d0)*T1*T4
         DO 1 J=1,2
         RJ=R(J)
         AJ=ALFA*RJ
         BJ=AJ+BETA
         RL(J)=RJ**(L+1)*DSQRT(ALFA+BETA/RJ)
         U(J)=F(J)+2.d0*Z/RJ-(L*(L+1.d0)+BETA*(AJ+0.25d0*BETA)/(BJ*BJ))/
     *   (RJ*RJ)
1        CONTINUE
      UC=(R(2)*U(1)-R(1)*U(2))/(R(2)-R(1))
      UD=(U(2)-U(1))/(R(2)-R(1))
      T6=-2.d0*Z*T2*T4-UC*T5+UD*T4/T3
           DO 2 J=1,2
           RJ=R(J)
           T=T3*RJ*RJ
           ST(J)=RL(J)*(1.d0-RJ*T1+T*(T2+UC+RJ*T6))
           ST(2+J)=-RL(J)*T*(1.d0-RJ*T5)
2          CONTINUE
      RETURN
      END
C_________________________
C
      SUBROUTINE START(A,F)
C *** CALCULATION OF WAVE FUNCTION IN 2 POINTS NEAR ZERO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510),F(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */ST/ST(4)
           DO 3 J=1,2
           A(J)=A(N1+7)*(ST(J)+A(N1+8)*ST(2+J))
3          CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE ZERF(A,F)
C *** FINDING FOR ZERO APROXIMATION TO WAVE FUNCTION (ARRAY A)
C *** IN POTENTIAL FROM ARRAY F. CONDITIONS FOR WAVE FUNCTION
C *** - NUMBER OF ZEROS AND FINITY
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*6 LTEXT
      DIMENSION A(510),F(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */LOCAL/C(510),RESERV(510)
     */LTEXT/ LTEXT(8)
      PMAX=300.d0
      DEMAX=10.d0*EPS*ABS(A(N1+8))
      KR=A(N1+9)+0.01d0
      T2=H**2/12.d0
      EM=-(ION/A(N1+2))**2
      IF(A(N1+8).GT.EM) A(N1+8)=EM
      EM=-(Z/A(N1+2))**2
      IF (A(N1+8).LT.EM) A(N1+8)=EM
      EMIN=0.d0
      EMAX=0.d0
      IR=A(N1+2)-A(N1+3)-0.99d0
      CALL P12(A,F)
1     CALL START(A,F)
      A(N1+10)=A(N1+10)+1.0d0
      EK=A(N1+8)
            DO 2 J=1,KR
            BJ=ALFA+BETA/R(J)
            C(J)=1.d0-T2*(F(J)-EK)/(BJ*BJ)
2           CONTINUE
7     IT=0
            DO 11 J=3,KR
            A(J)=((12.d0-10.d0*C(J-1))*A(J-1)-C(J-2)*A(J-2))/C(J)
            IF (A(J)*A(J-1)) 9,8,10
8           IF (A(J)*A(J-2)) 9,10,10
9           IT=IT+1
            IF(IT.EQ.IR) JB=J
            IF(IT.GT.IR) GO TO 13
10          IF(ABS(A(J)).GT.PMAX) GO TO 14
11          CONTINUE
      IF(IR-IT) 13,16,14
13    EMAX=A(N1+8)
      IF (EMIN.LE.EMAX)  GO TO 15
      A(N1+8) = EMAX*1.5d0
      GO TO 1
14    EMIN=A(N1+8)
      IF(EMAX-EMIN.LE.DEMAX) GO TO 16
15    A(N1+8)=(EMAX+EMIN)/2.d0
      GO TO 1
16    IF(IR.EQ.0) JB=1
      IF (IT.LT.IR-1) GO TO 27
      JB=JB+1
            DO 17 J=JB,N1
            IF (ABS(A(J)).LT.ABS(A(J-1))) GO TO 20
17          CONTINUE
      GO TO 25
20    JBM=J
            DO 21 J=JBM,N1
            IF (ABS(A(J)).GT.ABS(A(J-1))) GO TO 22
21          CONTINUE
      GO TO 25
22    KR=J
      A(N1+9)=KR
      JB=KR+1
            DO 23 J=JB,N1
23          A(J)=0.d0
      GO TO 26
25    KR=N1
      A(N1+9)=KR
26    CONTINUE
      RRKKI=RKIN(0,A,A)
      ANORMA=1.d0/DSQRT(RRKKI)
      A(N1+7)=A(N1+7)*ANORMA
            DO 24 J=1,KR
24          A(J)=ANORMA*A(J)
      RETURN
27    M3=A(N1+3)+1.01d0
      WRITE(2,100) A(N1+2), LTEXT(M3)
100   FORMAT(10X,I3,A1,' - RMAX SMALL')
      STOP
      END
C_______________________________________________________________
C
      SUBROUTINE PERTUR(A)
C *** ENERGY CORRECTION ACCORDING PERTURBATION THEORY
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */LOCAL/RESERV(510),V1(510)
     */V/V(510)
2     CALL RD(NS+1,V1)
      KR=A(N1+9)+0.01d0
      J=1
      RJ=R(J)
      TJ=ALFA+BETA/RJ
      T2=A(J)*A(J)*(V(J)-V1(J))
      T1=T2/(TJ*TJ)
      ID=1
            DO 1 J=2,KR
            RJ=R(J)
            T2=A(J)*A(J)*(V(J)-V1(J))
            T1=T1+(3+ID)*T2/(ALFA+BETA/RJ)**2
1           ID=-ID
      DE=T1*H/3.d0
      A(N1+8)=A(N1+8)+DE
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE NMRV1(A,DP,C,JC)
C *** SOLVING OF ONE-PARTICLE EQUATION WITH POTENTIAL FROM ARRAY F
C *** WITHOUT EXCHANGE. LEFT FROM JC - NUMEROV METHOD
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510),C(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */LOCAL/EIGEN(510),D(510)
      KR=A(N1+9)+0.01d0
      JC1=JC+1
      DEL=0.1d0*EPS
      IR=A(N1+2)-A(N1+3)-0.99d0
      IT=0
            DO 1 J=3,JC1
            A(J)=((12.d0-10.d0*C(J-1))*A(J-1)-C(J-2)*A(J-2))/C(J)
            IF (A(J)*A(J-1)) 9,8,1
8           IF(A(J)*A(J-2)) 9,1,1
9           IT=IT+1
1           CONTINUE
      IF(IT.NE.IR) GO TO 6
      DP=A(JC1)
      A(JC1)=-12.d0+10.d0*C(JC1)
      D(JC1)=-C(JC)*A(JC)
      JC2=JC+2
            DO 2 J=JC2,KR
            T1=C(J-1)/A(J-1)
            A(J)=C(J)*(10.d0-T1)-12.0d0
            D(J)=-T1*D(J-1)
            IF(ABS(D(J)).LE.DEL) GO TO 3
2           CONTINUE
      J=KR
3     JK=J
            DO 4 J=JK,KR
4           A(J)=0.
      HKAPPA=DSQRT(12.d0*ABS(1.d0-C(JK-1)))
      A(JK-1)=D(JK-1)/(A(JK-1)+C(JK)*EXP(-HKAPPA))
      M=JK-JC1
            DO 5 K=2,M
            J=JK-K
5           A(J)=(D(J)-C(J+1)*A(J+1))/A(J)
      DP=DP-A(JC1)
      RETURN
6     JC=1001
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE PRATT(V,F)
C *** DEMPHING FOR BETTER CONVERGENCY (PRATT METHOD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(510),F(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */LOCAL/V1(510),F1(510)
      IF(ITER-1) 1,1,3
1     CONTINUE
      CALL WD(NS+1,V)
      CALL WD(NS+2,F)
            DO 2 J=1,N1
2           V(J)=0.5d0*(V(J)+F(J))
      RETURN
3     CALL RD(NS+1,V1)
      CALL RD(NS+2,F1)
      CALL WD(NS+1,V)
      CALL WD(NS+2,F)
            DO 9 J=1,N1
            DEN=F(J)+V1(J)-V(J)-F1(J)
            IF(DEN) 6,5,6
5           CJ=0.
            GO TO 8
6           CJ=(F(J)-F1(J))/DEN
            IF(CJ.LE.0) CJ=0.d0
            IF(CJ.GT.0.5d0) CJ=0.5d0
8           V(J)=CJ*V(J)+(1.d0-CJ)*F(J)
9           CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      FUNCTION HNOR(A)
C *** NORMALIZATION FOR HYDROGENLIKE FUNCTION:
C *** FOR DISCRETE STATES (2*Z/N)**(L+1)/N/(2*L+1)!*
C ***               SQRT(Z*(N+L)!/(N-L-1)!)
C *** FOR CONTINUUM (2*K)**(L+1)/K/(2*L+1)!*
C ***   SQRT(Z/(1-EXP(-2*PI*Z/K))*DERIV.(J**2+(Z/K)**2))
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
      ZZ=Z
      Z=A(N1+1)
      L=A(N1+3)+0.01d0
      A2=A(N1+2)
      IF (A2) 3,3,1
1     T1=1.d0
      T2=Z
      T3=A2+L+1.d0
      JM=2*L+1
      DO 2 J=1,JM
      T2=T2*(T3-J)
2     T1=T1*J
      HNOR=(2.d0*Z/A2)**(L+1)*DSQRT(T2)/(A2*T1)
      Z=ZZ
      RETURN
3     CONTINUE
      P=DSQRT(A(N1+8))
      PI=3.14159265d0
      HN=2.d0*DSQRT(Z/(1.d0-EXP(-2.d0*PI*Z/P)))
      IF (L) 4,5,6
4     WRITE(2,100)
100   FORMAT('MISTAKE, L<0')
      STOP
5     HNOR=H
      Z=ZZ
      RETURN
6     T1=1.d0
      T2=2.d0*P
      T3=1.d0
      T4=(Z/P)**2
      DO 7 J=1,L
      T1=T1*T2/(2.d0*J)/(2.d0*J+1.d0)
      T3=T3*(J*J+T4)
7     CONTINUE
      HNOR=HN*T1*DSQRT(T3)
      Z=ZZ
      RETURN
      END
C_______________________________________________________________
C
      SUBROUTINE OUTPUT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*6 LTEXT
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
     */LOCAL/A(510),ZV(510)/kv/kv/DAT1/zz,dra
     */V/V(510)
     */LTEXT/LTEXT(8)
      M1=N1+1
      M2=N1+2
      M4=N1+4
      M5=N1+5
      M7=N1+7
      M8=N1+8
      M9=N1+9
      zz = z
      IF (IUT) 12,12,11
11    WRITE(2,100) (R(J),J=1,KRMAX)
100   FORMAT(/,' COORDINATES R=',/,(7X,5d12.5))
      KV=V(N1+1)+0.01
            DO 3 J=1,KV
            RJ=R(J)
            AJ=ALFA*RJ
            BJ=AJ+BETA
3           ZV(J)=(-V(J)*RJ+BETA*(AJ+0.25d0*BETA)/(RJ*BJ*BJ))/(2.d0*Z)
      WRITE(2,101) (ZV(J),J=1,KV)
101   FORMAT(/,' SCREENING POTENTIAL',/,4X,3HZV=,5d12.5,/,
     *(7X,5d12.5))
12    WRITE(2,102) TEST
102   FORMAT (//,' TEST=',d12.5,6X,'RESULTS FOR SUBSHELLS',/)
            DO 1 NR=1,NS
            CALL RD(NR,A)
            L1=A(N1+3)+1.01d0
            II2=A(M2)+0.01d0
            II4=A(M4)+0.01d0
            II5=A(M5)+0.01d0
            KR=A(N1+9)+0.01d0
            WRITE (2,103) II2,LTEXT(L1),II4,A(M8),A(M1),A(M7),
     *      II5,KR,R(KR)
103         FORMAT(I3,A1,I2,2X,2HE=,F11.5,2X,4HZEF=,F8.4,2X,3HA7=,
     *      F9.4,2X,3HJC=,I3,2X,3HKR=,I3,2X,6HR(KR)=,F5.2)
            IF (IUT.LE.0) GO TO 1
            WRITE(2,104) (A(J),J=1,KR)
            WRITE(2,105)
104         FORMAT(/,5X,2HA=,5d12.5,/,(7X,5d12.5))
105         FORMAT(//)
1           CONTINUE
      RETURN
      END
C_______________________________________________________________
C
      FUNCTION RKIN(K,A,B)
C *** RKIN=INTEGRAL R**K*P(A,R)*P(B,R)*DR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(510),B(510)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/inn/H,ALFA,BETA,EPS,TEST,Z
     */R/R(500)
      AL=A(N1+3)
      BL=B(N1+3)
      J=1
      RJ=R(J)
      TJ=ALFA+BETA/RJ
      T1=A(N1+1)*(AL+BL+2.d0)/(AL+1.d0)/(BL+1.d0)/(AL+BL+K+4.d0)
      T2=A(J)*B(J)*RJ**K
      T3=(1.d0+RJ*T1)/(AL+BL+K+3.d0)
      S1=T2*T3*RJ/TJ
      KR=DMIN1(A(N1+9),B(N1+9))+0.01d0
      T1=T2/(TJ*TJ)
      ID=1
            DO 1 J=2,KR
            RJ=R(J)
            TJ=ALFA+BETA/RJ
            T2=A(J)*B(J)*RJ**K
            T1=T1+(3+ID)*T2/(TJ*TJ)
1           ID=-ID
      S2=T1*H/3.d0
      RKIN=S1+S2
      RETURN
      END
C_________________________

      SUBROUTINE CONTIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/IN/KRMAX,ION,NN,N1,NS,NF,ITER,IUT
      common/r/r(500)/local/a(510),zv(510)/kv/kv
      COMMON/VD/VD(5000)/PE/DW(5000)/RP/RP(5000)
      COMMON/DAT1/Z,DRA/DAT2/NR,NRP
      dimension vv(500),faza(100),pcnt(5000,100),wen(100)

      read(1,10) ll,emin,del,emax
10    format(i3,3f8.4)

      if(z.eq.1.d0) then
         rmaxv=40.d0
         goto 11
      end if

      rmaxv = 2.d0*r(kv) - r(kv-1)

11    continue

      dra = 3.d0/(20.d0*dsqrt(emax))
      if(dra .gt. 0.1d0) dra = 0.1d0
      nr = 40
      r0 = r(1)
      hr = r0/2.d0



c** transformation from herman-skillman to dwave potential
      do  1 i=1,kv
         vv(i) = z*(1.d0 - zv(i))
1     continue
      zasim = dble(ion)

C*** RADIAL MESH AND FINDING THE BOUNDARY POINT
      call meshdw(r0,rmaxv,hr)
      hasim = rp(nrp) - rp(nrp-1)
      krp = nrp
c** making krp odd
      if((krp/2)*2 .eq. krp) krp = krp-1
c** interpolation of the dwave potential into the new mesh

c** for z=1 interpolation is not needed
      if(z.eq.1.d0) then
        do 12 i=1,krp
          vd(i) = 0.d0
12      continue
        goto 13
      end if

      do 2 i=1,krp
          if(rp(i).gt.r(kv)) then
              vd(i) = z - 1.d0
              goto 6
          endif
          call parinv(rp(i),r,vv,kv,vd(i))
6         continue
2      continue

13     continue

c** continuum waves and phases for all electron energies
       erg = emin
       j=1
4      continue
       wen(j) = erg
       call dwave(ll,erg,ph)
       do 3 i=1,krp
          pcnt(i,j) = dw(i)
3      continue
       faza(j) = ph
       erg = erg + del
       if(erg .gt. emax) goto 5
       j = j+1
       goto 4
5      continue
       nj = j

c** printout of phases = fcoulomb + fshort
       write(2,60) ll
60     format(/,'  orbital angular momentum in continuum = ',i2,/)
       write(2,30)
30     format(' Attention! in scat ph (and therefore full phase)',/,
     1      ' there is additional +Pi/2 from Berezhko code)',//,
     2      '   i',3x,'Ee(Ry)',2x,'full ph',2x,'coul ph',2x,'scat ph')
       do 7 i=1,nj
         b = facouz(wen(i),ll,zasim)
         write(2,20) i, wen(i), faza(i), b, faza(i)-b
20       format(i4,4f9.4)
7      continue

c*** printout of continuum functions
       write(2,50) '   r(a.u)   ', (wen(k), k=1,nj)
50     format(/,A12,20f12.4)
       do 8 n=1,krp
         write(2,40) rp(n), (pcnt(n,j), j=1,nj)
40       format(21d12.5)
8      continue

       return
       end

*************************************

      SUBROUTINE MESHDW(R0,RMAXV,HR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DAT1/Z,DRA/DAT2/NR,NRP
     */rp/rp(5000)
      h = hr
      rp(1) = R0
      nrp=1
      A=2.D0
          DO 2 K=1,5000
             DO 1 J=1,NR
             NRP=NRP+1
             rp(NRP)=rp(NRP-1)+H
             IF(rp(NRP).GT.RMAXV) RETURN
             IF(NRP.GT.5000) GO TO 3
1            CONTINUE
          H=A*H
          IF(H.GE.DRA) A=1.D0
2         CONTINUE
3     WRITE(6,100)
100   FORMAT(' LAST POINT LESS THEN RMAXV')
      RETURN
      END

*****************************************

      SUBROUTINE DWAVE(L,E,PHASE)
      implicit double precision (a-h,o-z)
      COMMON/VD/VD(5000)/PE/DW(5000)/RP/RP(5000)
      COMMON/DAT1/Z,DRA/DAT2/NR,NRP
      DIMENSION F1(50),F2(50),G1(50),G2(50)
      DIMENSION FC1(50),FC2(50),GC1(50),GC2(50)
      dimension fc(0:50),fcp(0:50),gc(0:50),gcp(0:50)
      pi=dacos(-1.d0)
      L1=L+1
      PHASE=0.0d0
      ANORM=0.0d0
      SE=DSQRT(E)
      AL=L
      PHI=0.0d0
      M=1
      DO 5 I=1,L1
5     PHI=PHI+DLOG(dble(2*I+1))
      R1=DEXP((PHI-35.0d0)/(AL+1.0d0))/SE
      NRP2=NRP-2
      DO 1 J=1,NRP2
      IF(RP(J).GE.R1) GO TO 2
1     DW(J)=0.0d0
      RETURN
2     J1=(J-1-(J-1)/NR*NR)/(NR-1)
      ALL=AL*(AL+1.0d0)
      NR1=NR+1
      NR2=NR+2
      NR3=NR2+1
      N1=J-J1
      N2=N1+1
      N3=N2+1
      X=RP(N1)*SE
      AN=(VD(N1)-Z)/SE
      AL1=AL+1.0d0
6     B1=1.d0
      XX=X*X
      B2=AN/AL1*X
      FL=B1+B2
      S=0.0d0
7     S=S+1
      BS=B2
      B2=(2.d0*AN*B2*X-B1*XX)/((2.d0*AL1+S)*(S+1.d0))
      B1=BS
      FL=FL+B2
      IF(B1*1.d-15.NE.0..OR.B2*1.d-15.NE.0) GO TO 7
      DO 8 I=1,L1
8     FL=FL*X/(dble(2*I-1))
      GO TO (11,12) , M
11    DW(N1)=FL
      AN=(VD(N2)-Z)/SE
      X=RP(N2)*SE
      M=2
      GO TO 6
12    DW(N2)= FL
      H=RP(N2)-RP(N1)
      HH=H*H
      C1=5.d0/6.d0
      C2=1.d0/12.d0
      DO 3 J=N3,NRP
      J1=J-1
      RJ=RP(J)
      RJ1=RP(J1)
      H1=H
      H=RJ-RJ1
      HH=H*H
      JD=(H1+H)/(2.9d0*H1)
      J2=J-2-JD
      RJ2=RP(J2)
      Y=(2.d0+C1*HH*(ALL/RJ1**2-2.d0*(Z-VD(J1))/RJ1-E))*DW(J1)
     *-(1.d0-C2*HH*(ALL/RJ2**2-2.d0*(Z-VD(J2))/RJ2-E))*DW(J2)
      DW(J)=Y/(1.d0-C2*HH*(ALL/RJ**2-2.d0*(Z-VD(J))/RJ-E))
3     CONTINUE
      DW1=DW(NRP-1)
      DW2=DW(NRP)
      ZA=Z-VD(NRP)
      AN=-ZA/SE
      X1=RP(NRP-1)*SE
      X2=RP(NRP)*SE
      CALL COUL90(X1,AN,0.d0,L,FC,GC,FCP,GCP,0,ifail)
      F1(L1)=FC(L)
      FC1(L1)=FCP(L)
      G1(L1)=GC(L)
      GC1(L1)=GCP(L)
      CALL COUL90(X2,AN,0.d0,L,FC,GC,FCP,GCP,0,ifail)
      F2(L1)=FC(L)
      FC2(L1)=FCP(L)
      G2(L1)=GC(L)
      GC2(L1)=GCP(L)
c      CALL COUL(X1,AN,L,L,F1,FC1,G1,GC1,1.0E-6,100.0)
c      CALL COUL(X2,AN,L,L,F2,FC2,G2,GC2,1.0E-6,100.0)
13    DET=F1(L1)*G2(L1)-F2(L1)*G1(L1)
      IF(DET.EQ.0.d0) GO TO 10
      AC=(DW1*G2(L1)-DW2*G1(L1))/DET
      AS=(DW2*F1(L1)-DW1*F2(L1))/DET
      ANORM=AC*AC+AS*AS
      PHASE=pi/2.d0
      IF(AC.EQ.0.d0) GO TO 9
      PHASE=PHASE+DATAN(AS/AC)+FACOUZ(E,L,ZA)
      IF(AC.LT.0.d0) PHASE=PHASE+pi
9     ANORM=1./DSQRT(ANORM*pi*SE)
10    DO 4 J=1,NRP
4     DW(J)=DW(J)*ANORM
      RETURN
      END

********************************************

      DOUBLE PRECISION FUNCTION FACOUZ(E,L,Z)
      implicit double precision (a-h,o-z)
      GAM=-Z/DSQRT(E)
      M=201-L
      DO 19 K=1,M
      AL=DATAN(GAM/(202-K))
      IF(K.EQ.1) GO TO 18
      FACOUZ=FACOUZ-AL
      GO TO 19
18    BE=DSQRT(GAM*GAM+(202-K)**2)
      FACOUZ=AL*200.5d0+GAM*(DLOG(BE)-1.d0)
     *+(-DSIN(AL)/12.0d0+DSIN(3.d0*AL)/(360.d0*BE*BE))/BE
19    CONTINUE
      RETURN
      END

******************************************

      SUBROUTINE PARINV(X,A,F,N,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),F(N)
      IF(X.LT.A(1)) GO TO 11
      IF(X.GT.A(N)) GO TO 4
      K1=1
      K2=N
2     K3=K2-K1
      IF(K3.LE.1) GO TO 6
      K3=K1+K3/2
      IF(A(K3)-X) 7,8,9
7     K1=K3
      GO TO 2
9     K2=K3
      GO TO 2
8     R=F(K3)
      RETURN
3     B1=A(K1)
      B2=A(K1+1)
      B3=A(K1+2)
      B4=F(K1)
      B5=F(K1+1)
      B6=F(K1+2)
      R=B4*((X-B2)*(X-B3))/((B1-B2)*(B1-B3))+B5*((X-B1)*(X-B3))/
     *((B2-B1)*(B2-B3))+B6*((X-B1)*(X-B2))/((B3-B1)*(B3-B2))
      RETURN
6     IF(K2.NE.N) GO TO 3
      K1=N-2
      GO TO 3
4     C=DABS(X-A(N))
      IF(C.LT.0.1D-7) GO TO 5
      K1=N-2
13    WRITE(*,41) X
41    FORMAT(25H X IS OUT OF THE INTERVAL,3H X=,F13.7)
      GO TO 3
5     R=F(N)
      RETURN
11    C=DABS(X-A(1))
      IF(C.LT.0.1D-7) GO TO 12
      K1=1
      GO TO 13
12    R=F(1)
      RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE COUL90(X, ETA, XLMIN,LRANGE, FC,GC,FCP,GCP, KFN,IFAIL)
C----------------------------------------------------------------------
C
C  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
C
C  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
C   FOR REAL X .GT. 0. ,REAL ETA (INCLUDING 0.), AND REAL XLMIN .GT.-1.
C   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
C   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
C   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
C   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
C   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
C----------------------------------------------------------------------
C   CALLING VARIABLES; ALL REALS ARE DOUBLE PRECISION (REAL*8)
C
C   X       - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0
C             [ X > SQRT(ACCUR) : ACCUR IS TARGET ACCURACY 1.0D-14 ]
C   ETA     - REAL SOMMERFELD PARAMETER, UNRESTRICTED > = < 0.0
C   XLMIN   - REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER),
C             GENERALLY IN RANGE 0.0 - 1.0 AND MOST USUALLY 0.0
C   LRANGE  - INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED
C             FOR L-VALUES XLMIN TO XLMIN + LRANGE INCLUSIVE
C   FC ,GC  - REAL VECTORS F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
C   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
C             THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
C             STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
C   KFN     - INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
C           = 0         REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
C           = 1    SPHERICAL BESSEL      "      "     "        j & y
C           = 2  CYLINDRICAL BESSEL      "      "     "        J & Y
C
C   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
C   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
C   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8 IS
C   THE SMALLEST NUMBER WITH 1.+ACC8.NE.1. FOR OUR WORKING PRECISION.
C   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
C   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
C   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
C   THE VARIABLE PACCQ IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
C----------------------------------------------------------------------
C   ERROR RETURNS                THE USER SHOULD TEST IFAIL ON EXIT
C
C   IFAIL ON INPUT IS SET TO 0                        LIMIT = 20000
C   IFAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
C                   =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C                   =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C                   = -1 : X < 1D-7 = SQRT(ACCUR)
C                   = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES)
C----------------------------------------------------------------------
C  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
C           SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
C           IE 1D-30 FOR IBM REAL*8,    1D-150 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
C  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
C  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
C                 +  RCWFF       IN    CPC 11 (1976) 141-142
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
C  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
C  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188
C  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
C  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
C  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
C  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509     1.4.94
C----------------------------------------------------------------------
C  AUTHOR: A. R. BARNETT           MANCHESTER  MARCH   1981/95
C                                  AUCKLAND    MARCH   1991
C----------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          LRANGE, KFN, IFAIL
      DOUBLE PRECISION X, ETA, XLMIN
      DOUBLE PRECISION FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
C----- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL
      DOUBLE PRECISION ACCUR,ACCH,SMALL, ONE,ZERO,HALF,TWO,TEN2, RT2DPI
      DOUBLE PRECISION XINV,PK,CF1,C,D,PK1,ETAK,RK2,TK,DCF1,DEN,XLM,XLL
      DOUBLE PRECISION EL,XL,RL,SL, F,FCMAXL,FCMINL,GCMINL,OMEGA,WRONSK
      DOUBLE PRECISION WI, A,B, AR,AI,BR,BI,DR,DI,DP,DQ, ALPHA,BETA
      DOUBLE PRECISION E2MM1, FJWKB,GJWKB, P,Q,PACCQ, GAMMA,GAMMAI
      INTEGER          IEXP, NFP, NPQ, L, MINL,MAXL, LIMIT
      LOGICAL          ETANE0, XLTURN
      PARAMETER      ( LIMIT = 30000, SMALL = 1.0D-150 )
      COMMON  /STEED/  PACCQ,NFP,NPQ,IEXP,MINL    !not required in code
      COMMON  /DESET/  CF1,P,Q,F,GAMMA,WRONSK     !information only
C----------------------------------------------------------------------
C     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
C----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
      DATA RT2DPI /0.79788 45608 02865  D0/
CQ    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
C-----THIS CONSTANT IS  DSQRT(TWO / PI):
C-----USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND DOUBLE PRECISION
C----------------CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
                        ACCUR = 1.0D-18
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACCH  = DSQRT(ACCUR)
C-----   TEST RANGE OF X, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
                IF( X .LE. ACCH )                GO TO 100
      IF( KFN.EQ.2 )   THEN
         XLM = XLMIN - HALF
        ELSE
         XLM = XLMIN
        ENDIF
      IF( XLM.LE.-ONE .OR. LRANGE.LT.0 )         GO TO 105
      E2MM1  = XLM * XLM + XLM
      XLTURN = X * (X -  TWO * ETA) .LT. E2MM1
      E2MM1  = E2MM1  +  ETA * ETA
      XLL    = XLM + DFLOAT(LRANGE)
C-----  LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C-----  XLL  IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
C-----  DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
      MINL  = MAX0( IDINT(XLMIN + ACCUR),0 )     ! index from 0
      MAXL  = MINL + LRANGE
C-----   EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
      XINV = ONE / X
      DEN  = ONE                       ! unnormalised F(MAXL,ETA,X)
      PK   = XLL + ONE
      CF1  = ETA / PK  +  PK * XINV
           IF( DABS(CF1).LT.SMALL )    CF1 = SMALL
      RK2  = ONE
         D = ZERO
         C = CF1
C----- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
      DO 10 L =  1 , LIMIT             ! abort if reach LIMIT (20000)
          PK1 = PK + ONE
          IF( ETANE0 ) THEN
                ETAK = ETA / PK
                RK2  = ONE + ETAK * ETAK
                 TK  = (PK + PK1) * (XINV + ETAK / PK1)
             ELSE
                 TK  = (PK + PK1) * XINV
             ENDIF
          D   =  TK - RK2 * D          ! direct  ratio of B convergents
          C   =  TK - RK2 / C          ! inverse ratio of A convergents
            IF( DABS(C).LT.SMALL ) C = SMALL
            IF( DABS(D).LT.SMALL ) D = SMALL
          D   = ONE / D
          DCF1=   D * C
          CF1 = CF1 * DCF1
              IF( D.LT.ZERO )    DEN = -DEN
          PK  = PK1
          IF( DABS(DCF1-ONE).LT.ACCUR )     GO TO  20 ! proper exit
   10 CONTINUE
                                            GO TO 110 ! error exit
   20       NFP = PK - XLL - 1                        ! number of steps
              F = CF1                                 ! need DEN later
C----DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
      IF( LRANGE.GT.0 )       THEN
          FCMAXL    = SMALL  * DEN
          FCP(MAXL) = FCMAXL * CF1
          FC (MAXL) = FCMAXL
                    XL = XLL
                    RL = ONE
          DO 30 L =  MAXL, MINL+1, -1
             IF( ETANE0 )  THEN
                    EL = ETA / XL
                    RL = DSQRT( ONE + EL * EL )
                    SL = XL * XINV  + EL
                    GC (L) = RL                  ! storage
                    GCP(L) = SL
                ELSE
                    SL = XL * XINV
                ENDIF
             FC (L-1)  = ( FC(L)   * SL  +  FCP(L) ) / RL
             FCP(L-1)  =   FC(L-1) * SL  -  FC (L) * RL
             XL    =  XL - ONE                   ! end value is XLM
   30     CONTINUE
         IF( DABS(FC(MINL)).LT.ACCUR*SMALL )  FC(MINL) = ACCUR * SMALL
          F   = FCP(MINL) / FC(MINL)             ! F'/F at min L-value
          DEN = FC (MINL)                        ! normalisation
      ENDIF
C---------------------------------------------------------------------
C-----   NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C-----   EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
C---------------------------------------------------------------------
      IF( XLTURN ) CALL JWKB( X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP )
      IF( IEXP.GT.1 .OR. GJWKB.GT.(ONE / (ACCH*TEN2)) ) THEN
          OMEGA = FJWKB
          GAMMA = GJWKB * OMEGA
          P     = F
          Q     = ONE
        ELSE                                     ! find cf2
          XLTURN = .FALSE.
          PK =  ZERO
          WI =  ETA + ETA
          P  =  ZERO
          Q  =  ONE - ETA * XINV
          AR = -E2MM1
          AI =  ETA
          BR =  TWO * (X - ETA)
          BI =  TWO
          DR =  BR / (BR * BR + BI * BI)
          DI = -BI / (BR * BR + BI * BI)
          DP = -XINV * (AR * DI + AI * DR)
          DQ =  XINV * (AR * DR - AI * DI)
          DO 40 L = 1, LIMIT
             P  = P  + DP
             Q  = Q  + DQ
             PK = PK + TWO
             AR = AR + PK
             AI = AI + WI
             BI = BI + TWO
             D  = AR * DR - AI * DI + BR
             DI = AI * DR + AR * DI + BI
             C  = ONE / (D * D + DI * DI)
             DR =  C * D
             DI = -C * DI
             A  = BR * DR - BI * DI - ONE
             B  = BI * DR + BR * DI
             C  = DP * A  - DQ * B
             DQ = DP * B  + DQ * A
             DP = C
      IF( DABS(DP)+DABS(DQ).LT.(DABS(P)+DABS(Q)) * ACCUR ) GO TO 50
   40     CONTINUE
                                              GO TO 120 ! error exit
   50     NPQ   = PK / TWO                              ! proper exit
          PACCQ = HALF * ACCUR / DMIN1( DABS(Q),ONE )
          IF( DABS(P).GT.DABS(Q) ) PACCQ = PACCQ * DABS(P)
C---------------------------------------------------------------------
C    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
C---------------------------------------------------------------------
          GAMMA   = (F - P) / Q
          GAMMAI  = ONE / GAMMA
          IF( DABS(GAMMA) .LE. ONE )  THEN
                 OMEGA  = DSQRT( ONE  +  GAMMA * GAMMA )
            ELSE
                 OMEGA  = DSQRT( ONE  +  GAMMAI* GAMMAI) * DABS(GAMMA)
            ENDIF
          OMEGA  = ONE / ( OMEGA * DSQRT(Q) )
          WRONSK = OMEGA
        ENDIF
C---------------------------------------------------------------------
C    RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C---------------------------------------------------------------------
      IF( KFN.EQ.1 )       THEN         !   spherical Bessel functions
                 ALPHA = XINV
                 BETA  = XINV
        ELSEIF( KFN.EQ.2 ) THEN         ! cylindrical Bessel functions
                 ALPHA = HALF * XINV
                 BETA  = DSQRT( XINV ) * RT2DPI
        ELSE                            ! kfn = 0,   Coulomb functions
                 ALPHA = ZERO
                 BETA  = ONE
        ENDIF
      FCMINL = DSIGN( OMEGA,DEN ) * BETA
      IF( XLTURN )   THEN
                        GCMINL =   GJWKB * BETA
        ELSE
                        GCMINL =  FCMINL * GAMMA
        ENDIF
      IF( KFN.NE.0 )    GCMINL = -GCMINL         ! Bessel sign differs
      FC (MINL) = FCMINL
      GC (MINL) = GCMINL
      GCP(MINL) = GCMINL * (P - Q * GAMMAI - ALPHA)
      FCP(MINL) = FCMINL * (F - ALPHA)
      IF( LRANGE.EQ.0 )                          RETURN
C---------------------------------------------------------------------
C    UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
C    RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C      XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
C---------------------------------------------------------------------
      OMEGA = BETA * OMEGA / DABS(DEN)
                 XL = XLM
                 RL = ONE
      DO 60  L = MINL+1 , MAXL                   ! indexed from 0
                 XL = XL + ONE
          IF( ETANE0 ) THEN
                 RL = GC (L)
                 SL = GCP(L)
            ELSE
                 SL =  XL * XINV
            ENDIF
          GC (L)  = ( (SL - ALPHA) * GC(L-1) - GCP(L-1) ) / RL
          GCP(L)  =    RL *  GC(L-1)  -  (SL + ALPHA) * GC(L)
          FCP(L)  = OMEGA * ( FCP(L)  -  ALPHA * FC(L) )
          FC (L)  = OMEGA *   FC (L)
   60 CONTINUE
      RETURN
C------------------   ERROR MESSAGES
  100 IFAIL = -1
      WRITE(6,1000) X,ACCH
 1000 FORMAT(' FOR X = ',1PD12.3,'     TRY SMALL-X  SOLUTIONS,',
     *' OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)
                     RETURN
  105 IFAIL = -2
      WRITE (6,1005) LRANGE,XLMIN,XLM
 1005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',
     *I10,1P2D15.6/)
                     RETURN
  110 IFAIL =  1
      WRITE (6,1010) LIMIT, CF1,DCF1, PK,ACCUR
 1010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS',/
     *' CF1,DCF1,PK,ACCUR =  ',1P4D12.3/)
                     RETURN
  120 IFAIL =  2
      WRITE (6,1020) LIMIT,P,Q,DP,DQ,ACCUR
 1020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3/)
                     RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE  JWKB   (X,ETA,XL, FJWKB,GJWKB, IEXP)
      DOUBLE PRECISION    X,ETA,XL, FJWKB,GJWKB, DZERO
C----------------------------------------------------------------------
C-----COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
C-----AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C-----CALCULATED IN SINGLE, RETURNED IN DOUBLE PRECISION VARIABLES
C-----CALLS DMAX1, SQRT, ALOG, EXP, ATAN2, FLOAT, INT
C     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
C----------------------------------------------------------------------
      REAL    ZERO,HALF,ONE,SIX,TEN,RL35,ALOGE
      REAL    GH2,XLL1,HLL,HL,SL,RL2,GH,PHI,PHI10
      INTEGER IEXP, MAXEXP
      PARAMETER  ( MAXEXP = 300 )
      DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
      DATA DZERO,RL35,ALOGE /0.0D0, 35.0E0, 0.43429 45 E0 /
C----------------------------------------------------------------------
CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
      GH2   =  X * (ETA + ETA - X)
      XLL1  = DMAX1( XL * XL + XL, DZERO )
      IF( GH2 + XLL1 .LE. ZERO )                 RETURN
      HLL  = XLL1 + SIX / RL35
      HL   = SQRT(HLL)
      SL   = ETA / HL + HL / X
      RL2  = ONE + ETA * ETA / HLL
      GH   = SQRT(GH2 + HLL) / X
      PHI  = X*GH - HALF*( HL*ALOG((GH + SL)**2 / RL2) - ALOG(GH) )
      IF ( ETA.NE.ZERO ) PHI = PHI - ETA * ATAN2(X*GH,X - ETA)
      PHI10 = -PHI * ALOGE
      IEXP  =  INT(PHI10)
      IF ( IEXP.GT.MAXEXP ) THEN
           GJWKB = TEN**(PHI10 - FLOAT(IEXP))
      ELSE
           GJWKB = EXP(-PHI)
           IEXP  = 0
      ENDIF
      FJWKB = HALF / (GH * GJWKB)
C---------------------------------------------------------------------
C     END OF CONTINUED-FRACTION COULOMB & BESSEL PROGRAM  COUL90
C---------------------------------------------------------------------
      RETURN
      END

