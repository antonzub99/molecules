module procedures
  implicit none

contains
!
subroutine PL10(X, Y, XX, M, Eps, YY, Error)
! Вычисление значения функции с помощью
! интерполяционного процесса Эйткена-Лагранжа.
!declaration
real*8, intent(in):: X(:), Y(:), XX, Eps
integer, intent(in):: M
real*8:: YY
integer:: Error
real*8, parameter:: OneP=1.000001
real*8, allocatable:: Arg(:), Val(:)
integer:: i, j, k, L, R, n, M1
real*8 E, E0
!begin
  n=size(X)
  if (n==1) then
    YY=Y(1); Error=1; return
  end if
  L=1; R=n
  do while(L+1<R)
    k=(L+R)/2
    if (X(k)<=XX) then; L=k
                  else; R=k
    end if
  end do
  if (XX==X(L)) then
    YY=Y(L); Error=0; return
  end if
  if (M==1)then
    YY=Y(L); Error=1; return
  end if
  M1=M
  if (M1>n) M1=n
  allocate(Arg(M1), Val(M1))
  do i=1, M1
    if ((XX-X(L))<=(X(R)-XX)) then
      Arg(i)=X(L); Val(i)=Y(L)
      if (L==1) then
        do j=i+1, M1
          Arg(j)=X(R); Val(j)=Y(R); R=R+1
        end do
        goto 40
      end if
      L=L-1
    else
      Arg(i)=X(R); Val(i)=Y(R)
      if (R==n) then
        do j=i+1, M1
          Arg(j)=X(L); Val(j)=Y(L); L=L-1
        end do
        goto 40
      end if
      R=R+1
    end if
  end do
40 continue
  Arg(:)=Arg(:)-XX
  do k=2, M1
    do i=1, k-1
      Val(k)=(Arg(k)*Val(i)-Arg(i)*Val(k))/(Arg(k)-Arg(i))
    end do
    if (k>=3) then
      E=abs((Val(k)-Val(k-1))/Val(k))
      if (E<Eps) then
        YY=Val(k); deallocate(Arg, Val)
        Error=0; return
      end if
      if (k>=4.AND.E>E0*OneP) then
        YY=Val(k-1); deallocate(Arg, Val)
        Error=2; return
      end if
      E0=E
    end if
  end do
  YY=Val(M1); deallocate(Arg, Val)
  Error=1
  return
end subroutine PL10

double precision function DF63(F, H)
! Интегрирование функции F, заданнной таблично, по формуле Симпсона
double precision, intent(in):: F(:), H
double precision S, F12
integer i, n, n2
!begin
    n=size(F)
    if(n<3) then
        if(n<2) then
            DF63=0.0; goto 10
        endif
        DF63=0.5*H*(F(1)+F(2))
10      return
    endif
    S=0.0D0; n2=1
    if(mod(n,2)==0) then
        F12=(5.0D0*F(1)+15.0D0*F(2)-5.0D0*F(3)+F(4))/16.0D0
        S=(F(1)+4.0D0*F12+F(2))/2.0D0
        n2=2
    endif
    S=S+F(n2)-F(n)
    n2=n2+1
    do i=n2, n, 2
        S=S+4.0D0*F(i)+2.0D0*F(i+1)
    enddo
    DF63=sngl(H*S/3.0D0)
    return
end function DF63

double precision function DF67(X, F)
    !Интегрирование функции F, заданнной таблично
! на неравномерной сетке по формуле трапеций
double precision, intent(in):: X(:), F(:)
double precision S
integer i, n
!begin
    n=size(X)
    S=0.0D0
    do i=1, n-1
        S=S+(F(i)+dble(F(i+1)))*(X(i+1)-dble(X(i)))
    enddo
    DF67=sngl(0.5D0*S)
    return
end function DF67

double precision function Lezh(n,x) result(c)
integer, intent(in):: n
double precision, intent(in):: x
double precision a, b
integer i
!begin
  if(n==0) then
    c=1.0
  else if(n==1) then
    c=x
  else
    a=1.0; b=x
    do i=2,n
      c=b*x+float(i-1)*(b*x-a)/float(i)
      a=b; b=c
    end do
  end if
  return
end function Lezh

end module procedures

program main
    use procedures
!------------------------
!
!  E  - energy in Ry,
!  L  - orbital momentum
!  Z > 0 for attraction - charge
! -------------------------
    implicit double precision (a-h,o-z)
    COMMON/VD/VD(20000)/PE/DW(20000)/RP/RP(20000)
    COMMON/dwx/dwx(20000)/ox/ox(20000)/dip/dip(20000)
    common/AM/AM(20000)/scr/scr(20000)/arregw/arregw(20000)/arregwx/arregwx(20000)
    COMMON/DAT1/Z,DRA/DAT2/NR,NRP
    COMMON/SCAT/SCAT
    dimension points(78),potent(55),hesk_mesh(78),hesk_wf(78)
    DIMENSION FC1(0:50),FC0(0:50),GC1(0:50),GC0(0:50)
    dimension fcp1(0:50),fcp0(0:50),gcp1(0:50),gcp0(0:50)
    dimension H00(0:400), H30(0:400), H45(0:400), H60(0:400), H90(0:400)

    open(unit=1,file="carbon_pot.txt",status="old",iostat=ios)
    open(unit=2,file="carbon_mesh.txt",status="old")

    do k=1,55
        read(1,*) potent(k)
    end do

    do k=1,78
        read(2,*) points(k)
    end do

    close(1)
    close(2)

    z=6.0d0
    r0=0.001d0
    rmaxv=20.d0
    hr=0.001d0
    dra=0.01d0
    nr=50
    m=50
    asimpt=1.0d0
    eps=1.e-5

    call meshdw(r0,rmaxv,hr)

    open(unit=225,file="interpolated_potential.txt",status="replace")
    do k=1,nrp
        xx=rp(k)
        call pl10(points,potent,xx,20,eps,ww,ker)
        scr(k)=Z*ww
        if (rp(k)> points(55)) scr(k)=asimpt
        write(225,*)rp(k)," ",scr(k)
        vd(k)=Z-scr(k)
    end do

    write (*,*) "Enter Energy E, Orbital Momentum l"
    read (*,*) E,L

    call dwave(0,E,deltaphase)
    !deltaphase for later calculations
    call dwave(L,E,PHASE)
    SE=DSQRT(E)

    X00=RP(NRP)*SE
    X11=RP(NRP-1)*SE

    PI=3.141592653589793238462643D0

    ZA=Z-VD(NRP)
    AN1=-ZA/SE
    phcoul=facouz(e,l,za)

    CALL COUL90(X11,AN1,0.d0,L,FC1,GC1,FCP1,GCP1,0,ifail)
    CALL COUL90(X00,AN1,0.d0,L,FC0,GC0,FCP0,GCP0,0,ifail)
    arregw(NRP)=-dtan(phase-phcoul-pi/2.d0)*FC0(L)+GC0(L)
    arregw(NRP-1)=-dtan(phase-phcoul-pi/2.d0)*FC1(L)+GC1(L)

    CALL NUMEROV(E,L,HR)

    !normalizing
    do k=1,nrp
        arregw(k)=arregw(k)/dsqrt(PI*SE)
    end do

    open(664, file="carbon_mesh.txt", status="old")
    open(665, file="1s_carb.txt",status="old")
    do k=1, 78
        read(664,*)hesk_mesh(k)
        read(665,*)hesk_wf(k)
    end do

    open(unit=229,file="interpolated_carbon_hesk.txt",status="replace")
    do k=1,nrp
        xx=rp(k)
        call pl10(hesk_mesh,hesk_wf,xx,20,eps,ww,ker)
        ox(k)=ww*dsqrt(rp(k)/(1.141783 * rp(k) + 0.5))
        dip(k)=ox(k)*dw(k)*rp(k)
        write(229,*)rp(k)," ",ox(k)
    end do

    dmatrixel = DF67(RP, dip)

    open (001, file="angular_distribution_0.txt",status="replace")
    open (301, file="angular_distribution_30.txt",status="replace")
    open (451, file="angular_distribution_45.txt",status="replace")
    open (601, file="angular_distribution_60.txt",status="replace")
    open (901, file="angular_distribution_90.txt",status="replace")

    call angular_distrib(H00, E, PHASE, 0, dmatrixel, deltaphase)
    call angular_distrib(H30, E, PHASE, 30, dmatrixel, deltaphase)
    call angular_distrib(H45, E, PHASE, 45, dmatrixel, deltaphase)
    call angular_distrib(H60, E, PHASE, 60, dmatrixel, deltaphase)
    call angular_distrib(H90, E, PHASE, 90, dmatrixel, deltaphase)

    do i = 0, 360
        write(001, *)i,'  ', H00(i)
        write(301, *)i,'  ', H30(i)
        write(451, *)i,'  ', H45(i)
        write(601, *)i,'  ', H60(i)
        write(901, *)i,'  ', H90(i)
    end do

    close (001)
    close (301)
    close (451)
    close (601)
    close (901)

    RETURN
end

subroutine angular_distrib(Hd, E, PHASE, mu, dmatrixel, deltaphase)
    use procedures

    implicit double precision (a-h, o-z)
    COMMON/VD/VD(20000)/PE/DW(20000)/RP/RP(20000)
    COMMON/dwx/dwx(20000)/ox/ox(20000)/dip/dip(20000)
    common/AM/AM(20000)/scr/scr(20000)/arregw/arregw(20000)/arregwx/arregwx(20000)
    dimension Hd(0:400),faza(0:10),anwf(0:10)
    complex(8) :: imaginary_one
    complex(8) :: angular_factor
    complex(8) :: scatter_fact
    complex(8) :: hi_plus
    complex(8) :: psipr
    complex(8) :: ci,c2i,czero,cone,ctwo,hiplus,sctfac

    SE=DSQRT(E)
    PI=3.141592653589793238462643D0

    ci = dcmplx(0.d0,1.d0)
    c2i = dcmplx(0.d0,2.d0)
    czero = dcmplx(0.d0,0.d0)
    cone = dcmplx(1.d0,0.d0)
    ctwo = dcmplx(2.d0,0.d0)

    psipr = dcmplx(0.d0,0.d0)
    imaginary_one = dcmplx(0.d0,1.d0)
    scatter_fact = (cdexp(2*imaginary_one*deltaphase) - 1)/(2*imaginary_one*SE)
    sctfac = (cdexp(c2i*dcmplx(deltaphase,0.d0)) - cone)/(c2i*dcmplx(SE,0.d0))
    amolecular_coordinate = 2.324364d0

    ! pl10 - subroutine for interpolation

    call pl10(rp, dw, amolecular_coordinate, 20, eps, dwmol, ker)
    call pl10(rp,arregw,amolecular_coordinate,20,eps, arregwmol,ker)


    hi_plus = imaginary_one * dwmol - arregwmol
    hiplus = ci * dcmplx(dwmol-arregwmol,0.d0)

    !anewwf = interpolated dwave regular wavefunction

    do k=0, 10
        call dwave(k, E, phaseo)
        call pl10(rp,dw,amolecular_coordinate,20,eps,anewwf,ker)
        anwf(k) = anewwf
        faza(k) = phaseo
        write(*,*)anewwf,' ', phaseo
    end do

    do itheta = 0, 360
        psipr = (0.d0,0.d0)
        tmptheta = dble(itheta) * 2.d0 * PI / 360.d0
        do k=0, 10
            psipr =(((imaginary_one)**k)*(cdexp(imaginary_one*faza(k)))*anwf(k)*((-1.d0)**k)*Lezh(k,dcos(tmptheta)))+psipr
        end do
        psipr = psipr * dcmplx(4 * PI)
        !ASSUMING RHO - MOLECULAR VECTOR - IS FIXED
        !mu = angle(e, rho) - free parameter
        !theta = angle(rho, k) - angular distribution as function of theta
        !mu + theta = angle (e, k)
        amu1 = dble(mu) * 2.d0 * PI / 360.d0
        scal1 = dcos(amu1+tmptheta)
        scal2 = dcos(amu1)
        angular_factor=scal1+imaginary_one*cdexp(-imaginary_one*PHASE)*scatter_fact*E*hi_plus*psipr*scal2
        Hd(itheta) = (dmatrixel**2) * (cdabs(angular_factor))**2
    end do
    return
    end





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
      R=B4*((X-B2)*(X-B3))/((B1-B2)*(B1-B3))+B5*((X-B1)*(X-B3))*((B2-B1)*(B2-B3))+B6*((X-B1)*(X-B2))/((B3-B1)*(B3-B2))
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

subroutine NUMEROV(E,L,hr)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    common/vd/vd(20000)/rp/rp(20000)
    COMMON/DAT1/Z,DRA/DAT2/NR,NRP
    common/arregw/arregw(20000)
    common/arregwx/arregwx(20000)
    common/dwx/dwx(20000)
    common/AM/AM(20000)
    COMMON/scr/scr(20000)
    LL=L*(L+1)
    do J=1,NRP
        AM(J)=(-1)*LL/((RP(J))**2.d0)+2.d0*(scr(J))/RP(J)+E
    end do
    do I=NRP-1,2,-1
        A1=1.d0+AM(I+1)*hr*hr/12.d0
        A2=2.d0-(5.d0/6.d0)*hr*hr*AM(I)
        A3=1.d0+AM(I-1)*hr*hr/12.d0
        arregw(I-1)=(A2*arregw(I)-A1*arregw(I+1))/A3
        arregwx(I-1)=(A2*arregwx(I)-A1*arregwx(I+1))/A3
        dwx(I-1)=(A2*dwx(I)-A1*dwx(I+1))/A3
    end do
    return
    end

DOUBLE PRECISION FUNCTION FACOUZ(E,L,Z)
!------------------------
!  Coulomb phase for
!  E  - energy in Ry,
!  L  - orbital momentum
!  Z > 0 for attraction - charge
! -------------------------
      implicit double precision (a-h,o-z)
      GAM=-Z/DSQRT(E)
      M=201-L
      DO 19 K=1,M
      AL=DATAN(GAM/(202-K))
      IF(K.EQ.1) GO TO 18
      FACOUZ=FACOUZ-AL
      GO TO 19
18    BE=DSQRT(GAM*GAM+(202-K)**2)
      FACOUZ=AL*200.5d0+GAM*(DLOG(BE)-1.d0)*+(-DSIN(AL)/12.0d0+DSIN(3.d0*AL)/(360.d0*BE*BE))/BE
19    CONTINUE
      RETURN
      END


       SUBROUTINE MESHDW(R0,RMAXV,HR)
!-----------------------------------

!--------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DAT1/Z,DRA/DAT2/NR,NRP/rp/rp(20000)
      h = hr
      rp(1) = R0
      nrp=1
      A=1.D0
          DO 2 K=1,20000
             DO 1 J=1,NR
             NRP=NRP+1
             rp(NRP)=rp(NRP-1)+H
             IF(rp(NRP).GT.RMAXV) RETURN
             IF(NRP.GT.20000) GO TO 3
1            CONTINUE
          H=A*H
          IF(H.GE.DRA) A=1.D0
2         CONTINUE
3     WRITE(6,100)
100   FORMAT(' LAST POINT LESS THEN RMAXV')
      RETURN
      END



!------------------------------------------
      SUBROUTINE DWAVE(L,E,PHASE)

!-----------------------------------------
      implicit double precision (a-h,o-z)
      COMMON/VD/VD(20000)/PE/DW(20000)/RP/RP(20000)
      COMMON/DAT1/Z,DRA/DAT2/NR,NRP
      COMMON/SCAT/SCAT
      DIMENSION F1(50),F2(50),G1(50),G2(50)
      DIMENSION FC1(50),FC2(50),GC1(50),GC2(50)
      dimension fc(0:50),fcp(0:50),gc(0:50),gcp(0:50)
      pi=dacos(-1.d0)
      L1=L+1
      PHASE=0.0d0
      ANORM=0.0d0
      SCAT=0.d0
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
      Y=(2.d0+C1*HH*(ALL/RJ1**2-2.d0*(Z-VD(J1))/RJ1-E))*DW(J1)-(1.d0-C2*HH*(ALL/RJ2**2-2.d0*(Z-VD(J2))/RJ2-E))*DW(J2)
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
13    DET=F1(L1)*G2(L1)-F2(L1)*G1(L1)
      IF(DET.EQ.0.d0) GO TO 10
      AC=(DW1*G2(L1)-DW2*G1(L1))/DET
      AS=(DW2*F1(L1)-DW1*F2(L1))/DET
      ANORM=AC*AC+AS*AS
      PHASE=pi/2.d0
      IF(AC.EQ.0.d0) GO TO 9
      SCAT=datan(as/ac)
      PHASE=PHASE+DATAN(AS/AC)+FACOUZ(E,L,ZA)
      IF(AC.LT.0.d0) PHASE=PHASE+pi
9     ANORM=1./DSQRT(ANORM*pi*SE)
10    DO 4 J=1,NRP
4     DW(J)=DW(J)*ANORM
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE COUL90(X, ETA, XLMIN,LRANGE, FC,GC,FCP,GCP, KFN,IFAIL)

!----------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          LRANGE, KFN, IFAIL
      DOUBLE PRECISION X, ETA, XLMIN
      DOUBLE PRECISION FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
!----- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL
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
!----------------------------------------------------------------------
!     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
!----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
      DATA RT2DPI /0.797884560802865D0/
      DATA RT2DPI /0.79788456080286535587989211986876373D0/
!-----THIS CONSTANT IS  DSQRT(TWO / PI):
!-----USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND DOUBLE PRECISION
!----------------CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
                        ACCUR = 1.0D-18
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACCH  = DSQRT(ACCUR)
!-----   TEST RANGE OF X, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
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
!-----  LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
!-----  XLL  IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
!-----  DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
      MINL  = MAX0( IDINT(XLMIN + ACCUR),0 )     ! index from 0
      MAXL  = MINL + LRANGE
!-----   EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
      XINV = ONE / X
      DEN  = ONE                       ! unnormalised F(MAXL,ETA,X)
      PK   = XLL + ONE
      CF1  = ETA / PK  +  PK * XINV
           IF( DABS(CF1).LT.SMALL )    CF1 = SMALL
      RK2  = ONE
         D = ZERO
         C = CF1
!----- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
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
!----DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
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
!---------------------------------------------------------------------
!-----   NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
!-----   EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
!---------------------------------------------------------------------
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
!---------------------------------------------------------------------
!    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
!---------------------------------------------------------------------
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
!---------------------------------------------------------------------
!    RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
!---------------------------------------------------------------------
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
!---------------------------------------------------------------------
!    UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
!    RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
!      XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
!---------------------------------------------------------------------
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
!------------------   ERROR MESSAGES
  100 IFAIL = -1
      WRITE(6,1000) X,ACCH
 1000 FORMAT(' FOR X = ',1PD12.3,'     TRY SMALL-X  SOLUTIONS,',/' OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)
                     RETURN
  105 IFAIL = -2
      WRITE (6,1005) LRANGE,XLMIN,XLM
 1005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',I10,1P2D15.6/)
                     RETURN
  110 IFAIL =  1
      WRITE (6,1010) LIMIT, CF1,DCF1, PK,ACCUR
 1010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS',/' CF1,DCF1,PK,ACCUR =  ',1P4D12.3/)
                     RETURN
  120 IFAIL =  2
      WRITE (6,1020) LIMIT,P,Q,DP,DQ,ACCUR
 1020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS',/' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3/)
                     RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE  JWKB   (X,ETA,XL, FJWKB,GJWKB, IEXP)
      DOUBLE PRECISION    X,ETA,XL, FJWKB,GJWKB, DZERO
!----------------------------------------------------------------------
!-----COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
!-----AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
!-----CALCULATED IN SINGLE, RETURNED IN DOUBLE PRECISION VARIABLES
!-----CALLS DMAX1, SQRT, ALOG, EXP, ATAN2, FLOAT, INT
!     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
!----------------------------------------------------------------------
      REAL    ZERO,HALF,ONE,SIX,TEN,RL35,ALOGE
      REAL    GH2,XLL1,HLL,HL,SL,RL2,GH,PHI,PHI10
      INTEGER IEXP, MAXEXP
      PARAMETER  ( MAXEXP = 300 )
      DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
      DATA DZERO,RL35,ALOGE /0.0D0, 35.0E0, 0.4342945E0 /
!----------------------------------------------------------------------
!CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR DOUBLE PRECISION
!----------------------------------------------------------------------
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
!---------------------------------------------------------------------
!     END OF CONTINUED-FRACTION COULOMB & BESSEL PROGRAM  COUL90
!---------------------------------------------------------------------
      RETURN
      END SUBROUTINE
