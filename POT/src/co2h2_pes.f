c gfortran *.f -o pes.dll -lgfortran -shared -Ofast
c=============================================
C Our paper used these potential, Feb.3,2010 
c---------------------------------------------
c     IMPLICIT REAL*8 (A-H,O-Z)
c     OPEN(11,FILE='h2-co2-pes.chk')
c     Do 300 XPHI=0.0D0, 90.0D0, 30.D0
c     Do 300 XTH1=0.0D0, 90.0D0, 15.0D0
c     DO 300 XTH2=0.0D0, 90.0D0, 15.0D0
c     DO 300   RX=2.0D0, 10.0D0, 0.1
c       CALL  CO2H2PES(XPHI,XTH1,XTH2,RX,V,1) 
c     WRITE(11,663)XPHI,XTH1,XTH2,RX,V
c 300  CONTINUE
c 663 FORMAT(1x,3f7.1,f8.2,4f20.4)
c     END
c=============================================

c******************************************************************************
C      REAL*8 FUNCTION CO2H2PES(XR,XTH1,XTH2,XPHI) 
        SUBROUTINE CO2H2PES( XR, XTH1, XTH2, XPHI, POTVALUE ) bind(C)
c******************************************************************************
c** Subroutine to generate values of the vibrationally averaged 4D-MLR analyic
c  potential energy surfaces for complexes formed between H2 and CO2 
c  isotopologues {12}C{16}O2  in vibrational level v3= 0 or 1, as determined by
c   Hui Li, Pierre-Nicholas Roy and Robert J. Le Roy [JCP, (2010, in press)]. 
c* On first call, input values of select 4D-MLR expansion parameters are taken
c  from DATA statements in subroutine PARAREAD;  subsequent calls generate 
C  additional potential function values for that same case.
c-------------------
c** Input variables:
c-------------------
c XR - distance between CO2 and H2 centre of mass in [Angst], pointing from
c      the center of mass of CO2 to the center of mass of H2.
c XTH1 - Jacobi angular coordinate 'theta1' in degrees, which is the angle 
c      between the vector XR pointing from the center of mass of CO2 to the 
c      center of mass of H2 and the vector pointing from O2 atom to O1. 
c XTH2 - Jacobi angular coordinate 'theta2' in degrees, which is the angle 
c      between the vector XR pointing from the center of mass of CO2 to the 
c      center of mass of H2 and the vector pointing from H2 atom to H1. 
c XPHI - Jacobi dihedral angular coordinate 'phi' in degrees,dihedral angle 
c      between the two half planes extending from the vector R to H1 and the 
c      O1 atom.
c---------------------
c** Output:   V [cm-1]  is the calculated interaction energy '\Delta{V}'.
c-----------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MXPARM
      PARAMETER (MXPARM=600,MMAX=16)
      INTEGER iv3,NPARM,NP,NQ,NDE,NRE,NCN,MCM,NS,NL,NPOW
      INTEGER NC6L1,NC6L2,NC6LMAX,NC8L1,NC8L2,NC8LMAX
      INTEGER NBETA(0:50),NBETAL1(0:50),NBETAL2(0:50),NBETALMAX(0:50)
      REAL*8  BETA(0:50),Pn1(0:50,0:50), Pn2(0:50,0:50)
      REAL*8  PV(MXPARM),PD(MXPARM),C6(MXPARM),C8(MXPARM)
      REAL*8  YC,Re,De,Vasy,RREF,AREF,AREFp,AREFq,Rep,CN,VLRe,
     2 dVLRedRe,phiINF,RTPp,RTPq,yp,yq,ype,dype,yPOW,XP,SUM,DSUM,VLR,
     3 XPW,DER,XDE,XRE,YTH1,YTH2,YPHI,CTH1,CTH2,STH1,STH2,CPHI,
     4 C5Sum,C6Sum,C8Sum,T0,DM(MMAX),DMP(MMAX),DMPP(MMAX)
      COMMON /DATABLK/PV,PI,ua2ba,Qa,Qb,C6,C8,RREF,CN,NBETAL1,
     1 NBETAL2,NBETALMAX,NDEL1,NDEL2,NDELMAX,NREL1,NREL2,NRELMAX,NCN,
     2 MCM,NP,NQ,NS,NL,NC6L1,NC6L2,NC6LMAX,NC8L1,NC8L2,NC8LMAX
c-----------------------------------------------------------------------
c      DATA IPAR/0/
c      SAVE IPAR
c-----------------------------------------------------------------------------
c      IF(IPAR.EQ.0) THEN
c          CALL PARAREAD(iv3)
c           IPAR= 1
c           ENDIF
        PI=DACOS(-1.0D0)
        RY=XR
        !YTH1=XTH1*PI/180.D0
        YTH1=XTH1
        !YTH2=XTH2*PI/180.D0
        YTH2=XTH2
       ! YPHI=XPHI*PI/180.D0
        YPHI=XPHI
        CTH1=DCOS(YTH1)
        CTH2=DCOS(YTH2)
        STH1=DSIN(YTH1)
        STH2=DSIN(YTH2)
        CPHI=DCOS(YPHI)
        !CTH1 = YTH1
        !CTH2 = YTH2
        !STH1 = SQRT(1-YTH1*YTH1)
        !STH2 = SQRT(1-YTH2*YTH2)
        !CPHI = YPHI


        CALL plmrb(Pn1,CTH1,2)
        CALL plmrb(Pn2,CTH2,2)
c caculate the derivative of the parameters of De not including
c the coefficient before, only the three angle function A(TH1,TH2,PHI) 
c with associated Legendre expansion.
c     De expantion 
      IP=0
      De=0.0d0
      DO 201 L2=0,NDEL2,2
      DO 201 L1=0,NDEL1,2
       LMIN=IABS(L1-L2)
c        MM=MIN0(L1,L2)
       LLMAX=MIN0(NDELMAX,L1+L2)
        DO 201 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IP=IP+1
          PD(IP)=al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)
          De=De+PD(IP)*PV(IP)
         ENDIF
 201   CONTINUE
c caculate the derivative of the parameters of Re not including
c the coefficient before, only the three angle function A(TH1,TH2,PHI) 
c with associated Legendre expansion.
c     Re expantion 
      NDE=IP
      Re=0.0d0
      DO 202 L2=0,NREL2,2
      DO 202 L1=0,NREL1,2
       LMIN=IABS(L1-L2)
c        MM=MIN0(L1,L2)
       LLMAX=MIN0(NRELMAX,L1+L2)
        DO 202 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IP=IP+1
          PD(IP)=al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)
          Re=Re+PD(IP)*PV(IP)
         ENDIF
  202  CONTINUE
       NRE=IP-NDE
c C6 coefficient form from Kumar and Meath Chem. Phys. 189 (1994) 467-477 function (11)
c which derived from Langhoff, Gordon and Karplus J. Chem. Phys. 55 (1971) 2126 function (9) 
c here we rewrite it using product of Associated Legendre function and exp(im\phi). 
       AREF= RREF*Re
       IF(RREF.LE.0.d0) AREF= Re
       AREFp= AREF**NP
       AREFq= AREF**NQ
       Rep= Re**NP
c C6 and C8 Coeffcients calculated based on the geometric mean averages of the 
c coefficients for CO2-CO2 and H2-H2 complex which are from F. Visser, P.E.S. Wormer 
c and W.P.J.H. Jacobs JCP 82,3753(1985) for H2 and R. Bukowski and K. Szalewicz 
c JCP 110,3785(1999)for CO2. 
c
c due to the instantaneous dipole of CO2, Q3 not equal to 0. the induction term is not
c zero, and can be expreesed as function 57 of "Permanent
c  and Induced Molecular Moments and Long-Range Intermolecular Forces"
c  A.D. Buckingham vol.12, page 107, Advances in Chemical Physics.) 
c
c  here, ua2ba is the vibrational averaged value of <u(Q3)^2> 
c  facter from hatree a0^6 to cm-1 A^6 is 4819.379496 
c  so the input alfbba should be atomic unit
c  alfbba, alfbper, alfbpar are the isotropic average, perpendicular and parallel 
c  polarizability of H2, which are the vibrationally averaged values for gound state H2 
c  reported by Bishop and Cheung JCP 72, 5125(1980) 
       alfbba=5.4140d0
       alfbpar=6.7632d0
       alfbper=4.7393d0
       C6ind=4819.379496d0*ua2ba*(0.5d0*alfbba*(3.0d0*CTH1**2+1.0d0)
     & +(1.0d0/6.0d0)*(alfbpar-alfbper)*(12.0d0*CTH1**2*CTH2**2
     & +3.0D0*STH1**2*STH2**2*CPHI**2-3.0D0*CTH1**2-1.0D0
     & +12.0D0*STH1*CTH1*STH2*CTH2*CPHI))
      fC6ind=C6ind/CN
c   obove function could also written as 
c   C6ind=4819.379496d0*ub2ba*(alfbba*(1.0d0*A^000+sqrt(5)*A^202)
c        +(alfbpar-alfbper)*(sqrt(5)/3*A^022+3/sqrt(5)*A^220-
c         sqrt(10)/sqrt(7)*A^222+sqrt(8)/sqrt(35)*A^224))   
      IK=0
      C6Sum=0.0d0
      DO 331 L2=0,NC6L2,2
      DO 331 L1=0,NC6L1,2
       LMIN=IABS(L1-L2)
c        MM=MIN0(L1,L2)
       LLMAX=MIN0(NC6LMAX,L1+L2)
        DO 331 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IK=IK+1
          C6Sum=C6Sum+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*C6(IK)
         ENDIF
  331  CONTINUE
c  added the contribution from induce force from averaged instantaneous dipole of CO2
       C6Sum=C6Sum+fC6ind
c** For normal inverse-power sum MLR case, with or without damping
c  calculate damping coeffient for H2-CO2 as below
c Ip(H2)=15.45 ev W. Kolos and J. Rychlewski JCP 98,3960(1993) 
c Ip(CO2)=13.79 ev N. Bussieres and P. Marmet Can. J. Phys. 55, 1889(1977) 
c Ip(H)=13.60 ev 
c Pd(H2)=[Ip(H2)/Ip(H)]^(2/3)=1.210833319 
c Pd(CO2)=[Ip(CO2)/Ip(H)]^(2/3)=1.021028904 
c bd(H2,CO2)=2.78*[2*Pd(H2)*Pd(CO2)/(Pd(H2)+Pd(CO2))]=3.079851734 
      VLRe= CN*C6Sum/Re**NCN
      dVLRedRe= -NCN*VLRe/Re
      IK=0
      C8Sum=0.0d0
      DO 332 L2=0,NC8L2,2
      DO 332 L1=0,NC8L1,2
       LMIN=IABS(L1-L2)
c        MM=MIN0(L1,L2)
       LLMAX=MIN0(NC8LMAX,L1+L2)
        DO 332 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IK=IK+1
          C8Sum=C8Sum+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*C8(IK)
         ENDIF 
  332  CONTINUE
c  V_elect accounts for the electrostatic interaction due to the 
c  molecular permanent dipole, quadrupoles etc. For CO2-H2 system
c  the dipole both for CO2 and H2 are zero, while quadrupoles are
c  not zero. Here, we decided to consider only the main quadrupole-
c  -quadrupole contribution. expreesed as function 56 of "Permanent
c  and Induced Molecular Moments and Long-Range Intermolecular Forces"
c  A.D. Buckingham vol.12, page 107, Advances in Chemical Physics. 
c  Quadrupole moment of H2 = 0.481 (ea0^2) from D.B. Lawson and
c  J. F. Harrison J. Phys. Chem. A 101,4781(1997) 
c  Quadrupole moment of CO2 =-3.201 (ea0^2) from A. Haskopoulos and
c  G. Maroulis  Chem. Phys. Lett., 417, 235(2006)  
c      Qb=0.481d0
c      Qa=-3.201d0
      C5Sum=0.75D0*Qa*Qb*9107.307339*(1.0-5.0*CTH1**2-5.0*CTH2**2
     & -15.0*CTH1**2*CTH2**2+2.0*(4.0*CTH1*CTH2-STH1*STH2*CPHI)**2)
       IF(MCM.GT.NCN) THEN
        MMN= MCM - NCN
       IF(NP.LE.MMN) MMN= 0
       IF(MMN.GT.0) THEN
         VLRe= (CN/Re**NCN)*(C6Sum+C8Sum/Re**MMN)
     &                                         -C5Sum/Re**5   ! electrostatic  C5 is positive     
         dVLRedRe= dVLRedRe - MCM*CN*C8Sum/Re**(MCM+1)
     &                        +5.0D0*C5Sum/Re**6
       ENDIF
         phiINF= DLOG(2.d0*De/VLRe)
       ENDIF
       RTPp= RY**NP
       RTPq= RY**NQ
       yp= (RTPp - AREFp)/(RTPp + AREFp)
       yq= (RTPq - AREFq)/(RTPq + AREFq)
       ype= (RTPp - Rep)/(RTPp + Rep)
c caculate the derivative of the parameters of BETA(N) not including
c the coefficient before, only the three angle function A(TH1,TH2,PHI)
c with associated Legendre expansion.
       NPOW= NS
       IF(RY.GE.Re) NPOW= NL
        yPOW= 1.d0 - yp
        SUM=0.0
        DSUM=0.0
        IP=NDE+NRE
        DO 204 J=0,NPOW
           BETA(J)=0.0
           NPS=0
         DO 203 L2=0,NBETAL2(J),2
         DO 203 L1=0,NBETAL1(J),2
            LMIN=IABS(L1-L2)
c            MM=MIN0(L1,L2)
            LLMAX=MIN0(NBETALMAX(J),L1+L2)
         DO 203 L=LMIN,LLMAX
            LTOT=L1+L2+L
           IF(MOD(LTOT,2).EQ.0) THEN
            IP=IP+1
            NPS=NPS+1
            PD(IP)=al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*yq**(J)
            BETA(J)=BETA(J)+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*PV(IP)
           ENDIF
  203  CONTINUE
         NBETA(J)=NPS
        SUM=SUM+BETA(J)*yq**(J)
        IF(RREF.LE.0.D0) DSUM=DSUM+yPOW*BETA(J)*(J)*yq**(J-1)
  204  CONTINUE
c  caculate the derivative of the parameters of Vasy 
        PD(IP+1)=1.0D0
        Vasy=PD(IP+1)*PV(IP+1)
        XP= SUM*yPOW+ phiINF*yp
c with and without damping function 
        VLR= CN*C6Sum/RY**NCN
     &                       -C5Sum/RY**5
        IF(MMN.GT.0) THEN
          VLR= (CN/RY**NCN)*(C6Sum+C8Sum/RY**MMN)
     &                           -C5Sum/RY**5
        ENDIF
         XPW= DEXP(-XP*ype) * VLR/VLRe
         YC= De*(1.d0 - XPW)**2-De+Vasy
         V=YC
C		 CO2H2PES=V
        POTVALUE = V
      RETURN
      END 
c-----------------------------------------------------------------------------
      SUBROUTINE PARAREAD() bind(C) 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MXPARM=600,MMAX=16,MXDE=82, MXRE=51, MXPHI=87,ISP=2)
      PARAMETER (MXC6=6, MXC8=10, isC=12, isO=16)
      INTEGER I,J,isv,iv3,IP,NDEL1,NDEL2,NDELMAX,NREL1,NREL2,NRELMAX
      INTEGER IFXP(MXPARM),NPARM,NP,NQ,NDE,NRE,NCN,MCM,NS,NL,NPOW
      INTEGER NC6L1,NC6L2,NC6LMAX,NC8L1,NC8L2,NC8LMAX
      INTEGER NBETA(0:50),NBETAL1(0:50),NBETAL2(0:50),NBETALMAX(0:50)
      INTEGER NB(0:5) 
      REAL*8  BETA(0:50),Pn1(0:50,0:50),Pn2(0:50,0:50),U2(2),CCN(2)
      REAL*8  PV(MXPARM),C6(MXPARM),C8(MXPARM),QCO2(2)
      REAL*8  DE(ISP,MXDE),RE(ISP,MXRE),PHIPHI(ISP,MXPHI)
      REAL*8  YC,RREF,Qa,Qb,ua2ba
      COMMON /DATABLK/PV,PI,ua2ba,Qa,Qb,C6,C8,RREF,CN,NBETAL1,
     1 NBETAL2,NBETALMAX,NDEL1,NDEL2,NDELMAX,NREL1,NREL2,NRELMAX,NCN,
     2 MCM,NP,NQ,NS,NL,NC6L1,NC6L2,NC6LMAX,NC8L1,NC8L2,NC8LMAX
c-----------------------------------------------------------------------
      iv3=0
      CALL fct(40)
      CALL fill3j(13,13,13)
      DATA  NDEL1/12/,NDEL2/6/,NDELMAX/12/
      DATA NREL1/12/,NREL2/4/,NRELMAX/12/
      DATA NP/4/,NQ/3/,NS/5/,NL/5/,RREF/1.5d0/
      DATA U2(1)/0.016690d0/,U2(2)/0.050617d0/
      DATA QCO2(1)/-3.1969d0/,QCO2(2)/-3.1887d0/,Qb/0.481d0/
      DATA NPOW/5/
      DATA NB(0)/0/, NBETAL1(0)/6/, NBETAL2(0)/6/,NBETALMAX(0)/12/
      DATA NB(1)/1/, NBETAL1(1)/4/, NBETAL2(1)/4/,NBETALMAX(1)/8/
      DATA NB(2)/2/, NBETAL1(2)/2/, NBETAL2(2)/2/,NBETALMAX(2)/4/
      DATA NB(3)/3/, NBETAL1(3)/2/, NBETAL2(3)/2/,NBETALMAX(3)/4/
      DATA NB(4)/4/, NBETAL1(4)/2/, NBETAL2(4)/2/,NBETALMAX(4)/4/
      DATA NB(5)/5/, NBETAL1(5)/2/, NBETAL2(5)/2/,NBETALMAX(5)/4/
      DATA NCN/6/, MCM/8/, CCN(1)/2.088237d5/,CCN(2)/2.089114d5/
      DATA NC6L1/2/,NC6L2/2/,NC6LMAX/4/
      DATA NC8L1/4/,NC8L2/2/,NC8LMAX/6/
      PI=DACOS(-1.0D0)
c     read C6 expantion coefficients      
      DATA (C6(I),I=1, MXC6)/1.00000D+00, 0.242102D+00, 0.112800D+00,
     1                      0.010447D+00, 0.027890D+00, 0.301495D+01/
c     read C8 expantion coefficients      
      DATA (C8(I),I=1, MXC8)/27.852490D+00,39.345277D+00,4.9441800D+00,
     1       8.1635650D+00, 0.5952160D+00, 2.4683140D+00,15.124290D+00, 
     2       0.0465500D+00, 0.0923000D+00, 2.5798610D+00/     
c    isC=12, isO=16, v3=0
      DATA (DE(1,I),I=1,MXDE)/7.7548D+01, -9.1820D+01, 9.4860D+01,
     &-6.4980D+01, 3.9230D+01, -2.2190D+01, 7.7900D+00, -4.0200D+01,
     & 2.3030D+01, 3.8790D+01, 4.6077D+02, -2.6010D+01, 1.0560D+01,
     &-2.3978D+02, 2.1210D+01, -1.0900D+01, 1.1460D+02, -1.4900D+01,
     & 6.7000D+00, -6.3500D+01, 1.4600D+01, -4.5000D+00, 2.6400D+01,
     &-1.3700D+01, 2.2000D+00, 7.7500D+00, -4.5600D+00, 7.4600D+00,
     &-1.3420D+01, 1.6300D+00, -4.5000D-01, 7.4000D+00, 5.3000D+00,
     & 6.0600D+01, -2.5300D+00, 1.5000D+00, -3.8000D+00, 0.0000D+00,
     &-4.7500D+01, 2.1000D+00, -1.1000D+00, 2.2000D+00, -1.7000D+00,
     & 1.8900D+01, -3.0000D+00, 1.6000D+00, -2.2000D+00, 2.4000D+00,
     & 2.9000D+00, -9.0000D-01, 5.0000D-01, -4.8000D-01, 6.1000D-01,
     &-4.0000D-01, 5.0000D-01, -3.1000D-01, 2.0000D-01, 0.0000D+00,
     & 5.0000D-01, -8.0000D-01, 1.1000D-01, -1.0000D-01, 0.0000D+00,
     &-2.0000D-01, 0.0000D+00, -1.6000D+00, -3.7000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/

      DATA (RE(1,I),I=1,MXRE)/3.74646D+00, 1.9456D+00, -6.7800D-01,
     & 2.8740D-01, -1.1050D-01, 2.6300D-02, 0.0000D+00, 2.4160D-01,
     &-3.8400D-02, -3.0100D-02, -1.4715D+00, 5.5900D-02,-7.1800D-02,
     & 5.9990D-01, -4.0900D-02, 2.4600D-02, 7.0000D-03, 3.0100D-02,
     & 0.0000D+00, -6.4700D-02, -4.7000D-03, -3.0000D-03, 4.0600D-02,
     &-1.4400D-02, 6.4000D-03, 5.6000D-03, 8.1000D-03, 2.1900D-02,
     &-2.0900D-02, -2.1000D-03, 5.2000D-03, 5.0000D-03, 4.1100D-02,
     & 3.2690D-01, 1.7000D-03, -2.4000D-03, 2.0000D-03, 6.0000D-03,
     &-1.8420D-01, -3.1000D-03, 0.0000D+00, 0.0000D+00, -7.0000D-03,
     & 4.4000D-02, 7.4000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     &-7.0000D-03, 0.0000D+00, 0.0000D+00/

      DATA (PHIPHI(1,I),I=1,MXPHI)/-2.5310D-01,1.3524D+00,-4.0990D-01,
     & 1.2400D-01, 2.1630D-01, 0.0000D+00, -6.7000D-02, -1.2400D-01,
     & 7.0000D-02, -3.4000D-02, 3.6300D-01, -6.3000D-02, 1.0000D-02,
     &-4.5000D-02, 0.0000D+00, 1.7000D-02, 1.4900D-01, -2.7400D-01,
     & 0.0000D+00, 3.0000D-03, -4.9000D-02, 7.2000D-02, 1.4000D-01,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.2000D-02,
     & 4.0000D-03, -1.0000D-03, 0.0000D+00, 0.0000D+00, 3.0000D-03,
     &-2.0000D-03, 7.0000D-03, -1.9000D-02, 2.9000D-02, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, -2.0000D-03, -2.8000D-02,
     &-6.9000D-02,                                                  ! end of Beta_0
     & 4.2100D-01, -8.4600D-01, -1.1200D-01, 1.5000D-01, 1.4800D-01,
     &-1.7000D-01, -2.4000D-01, 0.0000D+00, 0.0000D+00, 3.1000D-01,
     & 1.8000D-02, 2.0000D-02, 2.5800D-01, -3.7000D-01, 0.0000D+00,
     & 0.0000D+00, -7.0000D-02, 9.0000D-02, 6.1000D-01,             ! end of Beta_1
     &-2.3000D-01, -7.8000D-01, -2.2000D-01, 0.0000D+00, 0.0000D+00,
     & 2.3200D+00,                                                  ! end of Beta_2
     & 1.5000D-01, 0.0000D+00, -2.7800D-01, 0.0000D+00, 0.0000D+00,
     & 5.0000D-01,                                                  ! end of Beta_3
     & 1.3600D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, -9.9100D-01,
     &-7.8970D+00,                                                  ! end of Beta_4
     & 1.9449D+00, -3.8000D-01, 0.0000D+00, -1.4000D-01, -1.0000D+00,
     &-6.8000D+00/                                                  ! end of Beta_5

c    isC=12, isO=16, v3=1
      DATA (DE(2,I),I=1,MXDE)/7.7819D+01, -9.1340D+01, 9.4600D+01,
     &-6.4620D+01, 3.9020D+01, -2.2010D+01, 7.7000D+00, -3.9770D+01,
     & 2.2880D+01, 3.8860D+01, 4.5848D+02, -2.5770D+01, 1.0210D+01,
     &-2.4073D+02, 2.1020D+01, -1.0800D+01, 1.1378D+02, -1.4800D+01,
     & 6.7000D+00, -6.3200D+01, 1.4500D+01, -4.4000D+00, 2.6200D+01,
     &-1.3600D+01, 2.1000D+00, 7.6700D+00, -4.5100D+00, 7.4200D+00,
     &-1.3240D+01, 1.6100D+00, -4.8000D-01, 7.4000D+00, 5.3000D+00,
     & 5.9790D+01, -2.5100D+00, 1.5000D+00, -3.7000D+00, 0.0000D+00,
     &-4.7700D+01, 2.1000D+00, -1.1000D+00, 2.2000D+00, -1.6000D+00,
     & 1.8900D+01, -3.0000D+00, 1.6000D+00, -2.2000D+00, 2.3000D+00,
     & 2.9000D+00, -9.0000D-01, 5.0000D-01, -4.9000D-01, 5.7000D-01,
     &-4.0000D-01, 5.0000D-01, -3.0000D-01, 2.0000D-01, 0.0000D+00,
     & 6.0000D-01, -8.0000D-01, 1.1000D-01, -1.0000D-01, 0.0000D+00,
     &-3.0000D-01, 0.0000D+00, -1.6000D+00, -3.6000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/

      DATA (RE(2,I),I=1,MXRE)/3.74749D+00, 1.94547D+00, -6.7850D-01,
     & 2.8590D-01, -1.1020D-01, 2.6500D-02, 0.0000D+00, 2.4029D-01,
     &-3.8600D-02, -3.0700D-02,-1.4600D+00, 5.5700D-02, -7.0800D-02,
     & 6.1030D-01, -4.0800D-02, 2.4500D-02, 9.0000D-03, 2.9800D-02,
     & 0.0000D+00, -6.4300D-02, -5.2000D-03, -3.0000D-03, 4.0300D-02,
     &-1.4300D-02, 6.0000D-03, 5.5000D-03, 8.0000D-03, 2.1500D-02,
     &-1.9700D-02, -1.9000D-03, 5.7000D-03, 4.0000D-03, 4.0600D-02,
     & 3.2230D-01, 1.9000D-03, -2.3000D-03, 2.0000D-03, 6.0000D-03,
     &-1.8560D-01, -2.8000D-03, 0.0000D+00, 0.0000D+00, -7.0000D-03,
     & 4.4500D-02, 7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     &-7.0000D-03, 0.0000D+00, 0.0000D+00/

      DATA (PHIPHI(2,I),I=1,MXPHI)/-2.5750D-01,1.3569D+00,-4.1260D-01,
     & 1.2300D-01, 2.0690D-01, 0.0000D+00, -7.0000D-02, -6.0000D-02,
     & 7.2000D-02, -3.5000D-02, 3.6600D-01, -6.3000D-02, 1.3000D-02,
     &-4.9000D-02, 0.0000D+00, 1.1000D-02, 1.4400D-01, -2.5800D-01,
     & 0.0000D+00, 0.0000D+00, -5.0000D-02, 7.4000D-02, 9.0000D-02,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.5000D-02,
     & 3.3000D-03, -3.0000D-03, 0.0000D+00, 0.0000D+00, 3.0000D-03,
     &-3.0000D-03, 6.0000D-03, -1.6000D-02, 3.1000D-02, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, -2.0000D-03, -2.8000D-02,
     &-5.9000D-02,                                                  ! end of Beta_0
     & 4.0400D-01, -8.6000D-01, -1.1600D-01, 1.3000D-01, 1.5160D-01,
     &-1.7000D-01, -6.0000D-02, 0.0000D+00, 0.0000D+00, 3.2000D-01,
     & 1.8000D-02, 1.0000D-02, 2.5100D-01, -3.4000D-01, 0.0000D+00,
     & 0.0000D+00, -7.6000D-02, 9.0000D-02, 5.4600D-01,             ! end of Beta_1
     &-2.5000D-01, -8.0000D-01, -2.3000D-01, 0.0000D+00, 0.0000D+00,
     & 2.5700D+00,                                                  ! end of Beta_2
     & 1.2000D-01, 0.0000D+00, -2.8000D-01, 0.0000D+00, 0.0000D+00,
     & 5.0000D-01,                                                  ! end of Beta_3
     & 1.3100D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, -9.1300D-01,
     &-8.5110D+00,                                                  ! end of Beta_4
     & 1.9109D+00, -4.0000D-01, 0.0000D+00, -1.6000D-01, -9.0000D-01,
     &-7.3000D+00/                                                  ! end of Beta_5

      isv= 0
      IF((isC.EQ.12).and.(isO.EQ.16).and.(iv3.EQ.0)) isv= 1
      IF((isC.EQ.12).and.(isO.EQ.16).and.(iv3.EQ.1)) isv= 2
      ua2ba=U2(isv) 
      Qa=QCO2(isv) 
      CN=CCN(isv) 

c     IF(isv.GT.0) THEN
c         WRITE(6,600) isC,isO,iv3

c       ELSE
c         WRITE(6,602) isC,isO,iv3
c         isv= 1
c       ENDIF
      
      DO I=1,MXDE
c... first prepare well depth expansion parameters
          PV(I)=DE(isv,I)
          ENDDO
      IP=MXDE
      DO I=1,MXRE
c... next prepare potential minimum position expansion parameters
          IP=IP+1
          PV(IP)=RE(isv,I)
          ENDDO
      IP=MXDE+MXRE
c... then, prepare exponent coefficient \beta_i expansion parameters
      DO I=1,MXPHI
          IP=IP+1
          PV(IP)=PHIPHI(isv,I)
          ENDDO
c** Finally ... (if desired) read shift for definition of energy zero
      IP=MXDE+MXRE+MXPHI
      NPARM=IP+1
      PV(NPARM)=0.0D0
      RETURN
  600 FORMAT(/' Generate 4D-MLRR potential for {',i2,'}C{',i2,'}O2(v3=',
     1  i1,')-H2')
  602 FORMAT(/' *** 4D Potential H2-CO2 Function selection parameters  i
     1sC=',i2,'   isO=',i3,'   v3=',i2/'   do not match defined cases.'
     2  '  Generate potential for {12}C{16}O2(v=0)-H2 instead.')
      END
c***********************************************************************
c    
c This subroutine fills up the matrix w3j for C6 with values of 3-j coefficient
c     
      subroutine fill3j(l1max,l2max,lmax)
      implicit real*8 (a-h,o-z)
      common /w3jcg/w3j(0:50,0:30,0:30,0:60)
      dimension x(60)

      do l1=0,l1max
       do l2=0,l2max
        lmin=iabs(l1-l2)
        mm=min0(l1,l2)
        llmax=min0(lmax,l1+l2)

        do l=lmin,llmax
         do m=0,mm
          m1=m
          m2=-m
          mmm=0
          call cgc(l1,m1,l2,m2,l,mmm,c,1)
          w3j(m,l1,l2,l)=c
         end do
        end do
       end do
      end do

      return
      end
C ----------------------------------------------------------------------------
c
c Calculate the function C6Al1l2L for a given set of angles...
c It is assumed that the th1 and th2 angles are between the monomer bond
c and the "inner" part of the intermolecular axis.
c
       function al1l2l0(l1,l2,l,th1,th2,phi)
       implicit real*8 (a-h,o-z)
       dimension p1(0:50,0:50), p2(0:50,0:50)
       common /w3jcg/w3j(0:50,0:30,0:30,0:60)
       data izer/0/, ione/1/, pifact/12.566370614359d0/
c
       c1 = dcos(th1)
       !c1 = th1
       c2 = dcos(th2)
       !c2 = th2

       call plmrb(p1,c1,l1)
       call plmrb(p2,c2,l2)
       mmax = min(l1,l2)
       sum = 0.d0
       do m=1,mmax
        value=w3j(m,l1,l2,l)
        sum = sum + (-1)**m*value*p1(l1,m)*p2(l2,m)*dcos(m*phi)
       end do
       value=w3j(0,l1,l2,l)
       sum = 2*sum + value*p1(l1,0)*p2(l2,0)
c
       al1l2l0 = sum*pifact*(-1.d0)**(l1+l2+l)/dsqrt((2.d0*l1+1.d0)*
     1           (2.d0*l2+1.d0))
c
       return
       end

C --------------------------------------------------------------------------
c
c Compute the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials....
c
        subroutine plmrb(p,x,lmax)
        implicit real*8 (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
c inverse of dsqrt(2Pi)
        data twopinv /0.3989422804014d0/
c
c starting value
c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
c
c compute the diagonal elements
c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
c
c compute P_lm along the columns with fixed m

c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do
        end do
c
c Renormalize values...
c
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do
c
        return
        end

C -------------------------------------------------------------------------
c
c compute the matrix of N!
c
        subroutine fct(nmax)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        f(0) = 1.d0
        do i=1,nmax
         f(i) = f(i-1)*i
        end do
        return
        end
C -------------------------------------------------------------------------
c
Calculate the Clebsh-Gordan coefficient (or the 3-j symbol)
c The parameter ind3j.eq.1 indicates that the 3-J symbol is returned
c
        subroutine cgc(j1,m1,j2,m2,j,m,value,ind3j)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        d3jfact = 1.d0
        if(ind3j.eq.1) then
         d3jfact = ((-1.d0)**(j1-j2-m))/dsqrt(dfloat(2*j+1))
         m = -m
        endif

c
c Check the triangle conditions
c
        if(j.gt.(j1+j2)) write(6,*)'triangle violated'
        if(j.lt.abs(j1-j2)) write(6,*)'triangle violated'
        if((m1+m2).ne.m) then
          value = 0.d0
          return
        endif


c Calculation proper... the pre-sum factor....
c
        facn = (2*j+1)*f(j1+j2-j)*f(j1-m1)*f(j2-m2)*f(j+m)*f(j-m)
        facd = f(j1+j2+j+1)*f(j+j1-j2)*f(j+j2-j1)*f(j1+m1)*f(j2+m2)
        fac = dsqrt(facn/facd)

c
c determine the limit of k summation...
c
        kmax = min(j2+j-m1,j-m,j1-m1)
        if(kmax.lt.0) kmax = 0
        kmin = max(-j1-m1,-j2+j-m1,0)

c
c perform the summation (at least one cycle must be completed...
c
        sum = 0.d0
        do k = kmin,kmax
         facn = f(j1+m1+k)*f(j2+j-m1-k)
         facd = f(k)*f(j-m-k)*f(j1-m1-k)*f(j2-j+m1+k)
         sum = sum + (facn/facd)*(-1)**k
        end do
        value = d3jfact*fac*sum*(-1)**(j1-m1)
       return
       end
C -------------------------------------------------------------------------
       SUBROUTINE flush(nunit)
       endfile nunit
       backspace nunit
       end
C -------------------------------------------------------------------------


c***********************************************************************
      SUBROUTINE WGHT(NGP)
c** Subroutine to generate points (XG & X2) and weights (WG & W2) for
c  both regular (XG & WG) and singular  1/sqrt(1-X) (X2 & W2) 
c  Gaussian integration for  NGP=8 and 16 (for both) and for NGP=7 & 32
c  (for regular Gaussian only)
c-----------------------------------------------------------------------
c  Integraton formulas of Gaussian type(for orthogonal polynomials)
c  Gaus' Formula:
c    int f(x)dx -1-->1= sum[wi*f(xi)]+Rn 
c  related orthogonal polynomials: Legendre polynomials Pn(x), Pn(1)=1
c  abscissas: xi is the ith zero of Pn(x)
c  weights: wi=2/(1-xi^2)[P'_n(xi)]^2 
c     (b-a)^(2n+1)(n!)^4
c  Rn=------------------ f^(2n)(ksi)
c     (2n+1)[(2n)!]^3
c-----------------------------------------------------------------------
c  Gaus' Formula, Arbitrary Interval
c   int f(y)dy=[(b-a)/2]*sum[wi*f(yi)]+Rn
c   yi=[(b-a)/2]*xi+[(b+a)/2]   
c  Related orthogonal polynomials: Pn(x),Pn(1)=1
c  abscissas: xi is the ith zero of Pn(x)
c  Weights: wi=2/(1-xi^2)[P'_n(xi)]^2
c     (b-a)^(2n+1)(n!)^4
c  Rn=------------------ f^(2n)(ksi)
c     (2n+1)[(2n)!]^3
c----------------------------------------------------------------------

      INTEGER I,J,NGP,NGPH
      REAL*8  x7(4),w7(4),aa(4),bb(4),A(8),B(8),XX(16),WX(16)
c** Common block for quadrature weights & points
      REAL*8 XG(32),WG(32),X2(16),W2(16)
      COMMON /GWGHT/XG,WG,X2,W2
c
      data x7/0.949107912342759d0, 0.741531185599394d0,
     1        0.405845151377397d0, 0.d0/,
     2     w7/0.129484966168870d0, 0.279705391489277d0,
     3        0.381830550505119d0, 0.417959183673469d0/
      data aa/0.960289856497536d0, 0.796666477413627d0,
     1        0.525532409916329d0, 0.183434642495650d0/,
     2     bb/0.101228536290376d0, 0.222381034453374d0,
     3        0.313706645877887d0, 0.362683783378362d0/
      DATA A/0.989400934991649932596154D0,0.944575023073232576077988D0,
     1       0.865631202387831743880468D0,0.755404408355003033895101D0,
     2       0.617876244402643748446672D0,0.458016777657227386342420D0,
     3       0.281603550779258913230461D0,0.095012509837637440185319D0/,
     4     B/0.02715245941175409485178D0,0.06225352393864789286284D0,
     5       0.09515851168249278480993D0,0.12462897125553387205248D0,
     6       0.14959598881657673208150D0,0.16915651939500253818931D0,
     7       0.18260341504492358886676D0,0.18945061045506849628540D0/
      DATA XX/0.997263861849481563544981D0,0.985611511545268335400175D0,
     1        0.964762255587506430773812D0,0.934906075937739689170919D0,
     2        0.896321155766052123965307D0,0.849367613732569970133693D0,
     3        0.794483795967942406963097D0,0.732182118740289680387427D0,
     4        0.663044266930215200975115D0,0.587715757240762329040746D0,
     5        0.506899908932229390023747D0,0.421351276130635345364120D0,
     6        0.331868602282127649779917D0,0.239287362252137074544603D0,
     7        0.144471961582796493485186D0,0.048307665687738316234813D0/
     8    ,WX/0.00701861000947009660041D0,0.01627439473090567060517D0,
     9        0.02539206530926205945575D0,0.03427386291302143310269D0,
     1        0.04283589802222668065688D0,0.05099805926237617619616D0,
     2        0.05868409347853554714528D0,0.06582222277636184683765D0,
     3        0.07234579410884850622540D0,0.07819389578707030647174D0,
     4        0.08331192422694675522220D0,0.08765209300440381114277D0,
     5        0.09117387869576388471287D0,0.09384439908080456563918D0,
     6        0.09563872007927485941908D0,0.09654008851472780056676D0/
c** For simple Gaussian case
      NGPH= (NGP+1)/2
      J= NGP+1
      if(ngp.eq.7) then
          DO  I= 1,4
              J= J-1
              XG(I)= -x7(I)
              WG(I)= w7(I)
              XG(J)= x7(I)
              WG(J)= w7(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I= 1,NGP)
c*** For  1/SQRT(E-V)  case
c%%   WRITE(6,603) (XX(I),WX(I),I= 1,NGP)
          write(6,606) ngp
          go to 99
          endif
      if(ngp.eq.8) then
          DO  I= 1,4
              J= J-1
              XG(I)= -aa(I)
              WG(I)= bb(I)
              XG(J)= aa(I)
              WG(J)= bb(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I= 1,NGP)
c*** For  1/SQRT(E-V)  case
c%%   WRITE(6,603) (XX(I),WX(I),I= 1,NGP)
          DO    I= 1,8
              X2(I)= 1.d0-A(I)**2
c** Include SQRT(1-X) factor in weight, not integrand.
              W2(I)= 2.d0*B(I)*A(I)
              ENDDO
c%%       WRITE(6,604) NGP,(X2(I),W2(I),I=1,NGP)
          go to 99
          ENDIF
      if(ngp.eq.16) then
          DO  I=1,8
              J=J-1
              XG(I)=-A(I)
              WG(I)=B(I)
              XG(J)=A(I)
              WG(J)=B(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I=1,NGP)
c*** For  1/SQRT(E-V)  case
c%%   WRITE(6,603) (XX(I),WX(I),I=1,NGP)
          DO  I=1,16
              X2(I)=1.d0-XX(I)**2
c** Include SQRT(1-X) factor in weight, not integrand.
              W2(I)=2.d0*WX(I)*XX(I)
              ENDDO
c%%       WRITE(6,604) NGP,(X2(I),W2(I),I=1,NGP)
          go to 99
          endif
      if(ngp.eq.32) then
          DO  I=1,NGPH
              J=J-1
              XG(I)=-XX(I)
              WG(I)=WX(I)
              XG(J)=XX(I)
              WG(J)=WX(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I=1,NGP)
c          write(6,606) ngp
          go to 99
          endif
c      write(6,608) ngp
   99 RETURN
c 601 FORMAT(//2X,'Points and weights for simple Gaussian',i3,
c    1 '-point quadrature'/(1x,2(F25.20,F23.20)))
c 603 FORMAT(/' Simple Gaussian points and weights used for generating G
c    1auss-Mehler points and weights:'/(5X,2(F25.20,F23.20)))
c 604     FORMAT(/'  Points and weights for singular integrand',i3,
c    1     '-point quadratures'/(1x,2(F25.20,F23.20)))
  606 format(/'  NOTE: cannot generate singular integrand Gaussian point
     1s and weights for   NGP =',i3)
  608 format(/' NOTE:  cannot generate Gaussian points & weights for',
     1   '   NGP =',i3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12


      SUBROUTINE TQL(MD,N,Z,D,E)
C       Z(-1) A  Z  =D                                    
C       A = Z
C       EIGENVALUE         D(I)
C       EIGENFUNCTION      Z(J,I),J=1,N
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION   D(MD),E(MD),Z(MD,MD)
      EPS=1D-12
      NITER=50
      CALL TRED2(MD,N,Z,D,E)
      DO 10 I=2,N
  10  E(I-1)=E(I)
      F=0.0D0
      B=0.0D0
      E(N)=0.0D0
      DO 20 L=1,N
      J=0
      H=EPS*(DABS(D(L))+DABS(E(L)))
      LP1=L+1
      IF (B-H) 30,40,40
  30  B=H
  40  DO 50 M=L,N
      IF (DABS(E(M))-B) 60,60,50
  50  CONTINUE
  60  IF (M-L) 70,80,70
  70  IF (J-NITER) 90,100,90
  90  J=J+1
      P=(D(LP1)-D(L))/(2*E(L))
      R=DSQRT(P*P+1)
      IF (P) 110,111,111
  110 H=D(L)-E(L)/(P-R)
      GOTO 130
  111 H=D(L)-E(L)/(P+R)
  130 DO 140 I=L,N
  140 D(I)=D(I)-H
      F=F+H
      P=D(M)
      C=1.0D0
      S=0.0D0
      MM1=M-1
      IF (MM1-L) 270,280,280
  280 DO 120 LMIP=L,MM1
      I=L+MM1-LMIP
      IP1=I+1
      G=C*E(I)
      H=C*P
      IF (DABS(P)-DABS(E(I))) 160,170,170
  170 C=E(I)/P
      R=DSQRT(C*C+1.0D0)
      E(IP1)=S*P*R
      S=C/R
      C=1.0D0/R
      GOTO 180
  160 C=P/E(I)
      R=DSQRT(C*C+1)
      E(IP1)=S*E(I)*R
      S=1/R
      C=C/R
  180 P=C*D(I)-S*G
      D(IP1)=H+S*(C*G+S*D(I))
      DO 190 K=1,N
      H=Z(K,IP1)
      Z(K,IP1)=S*Z(K,I)+C*H
  190 Z(K,I)=C*Z(K,I)-S*H
  120 CONTINUE
  270 E(L)=S*P
      D(L)=C*P
      IF (DABS(E(L))-B) 80,80,70
  80  D(L)=D(L)+F
  20  CONTINUE
      DO 112 I=1,N
      IP1=I+1
      K=I
      P=D(I)
      IF (N-I) 230,230,300
  300 DO 210 J=IP1,N
      IF (D(J)-P) 220,210,210
  220 K=J
      P=D(J)
  210 CONTINUE
  230 IF (K-I) 240,112,240
  240 D(K)=D(I)
      D(I)=P
      DO 260 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
  260 Z(J,K)=P
  112 CONTINUE
      RETURN
  100 STOP '  FAIL'
      END




      SUBROUTINE TRED2(MD,N,Z,D,E)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  D(MD),E(MD),Z(MD,MD)
      BETA=1D-20
      DO 20 NMIP2=2,N
      I=N+2-NMIP2
      IM1=I-1
      IM2=I-2
      L=IM2
      F=Z(I,IM1)
      G=0.0D0
      IF (L) 30,30,40
  40  DO 50 K=1,L
  50  G=G+Z(I,K)*Z(I,K)
  30  H=G+F*F
      IF (G-BETA) 60,60,70
  60  E(I)=F
      H=0.0D0
      GOTO 180
  70  L=L+1
      IF (F) 80,90,90
  90  E(I)=-DSQRT(H)
      G=E(I)
      GOTO 100
  80  E(I)=DSQRT(H)
      G=E(I)
 100  H=H-F*G
      Z(I,IM1)=F-G
      F=0.0D0
      DO 110 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.0D0
      DO 201 K=1,J
  201 G=G+Z(J,K)*Z(I,K)
      JP1=J+1
      IF (JP1-L) 130,130,140
  130 DO 120 K=JP1,L
  120 G=G+Z(K,J)*Z(I,K)
  140 E(J)=G/H
      F=F+G*Z(J,I)
  110 CONTINUE
      HH=F/(H+H)
      DO 160    J=1,L
      F=Z(I,J)
      E(J)=E(J)-HH*F
      G=E(J)
      DO 170 K=1,J
  170 Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  160 CONTINUE
  180 D(I)=H
  20  CONTINUE
      D(1)=0.0D0
      E(1)=0.0D0
      DO 190 I=1,N
      L=I-1
      IF (D(I)) 202,210,202
  202 IF (L) 210,210,220
  220 DO 230 J=1,L
      G=0.0D0
      DO 240 K=1,L
  240 G=G+Z(I,K)*Z(K,J)
      DO 250 K=1,L
  250 Z(K,J)=Z(K,J)-G*Z(K,I)
  230 CONTINUE
  210 D(I)=Z(I,I)
      Z(I,I)=1.0D0
      IF (L) 260,260,270
  270 DO 280 J=1,L
      Z(I,J)=0.0D0
  280 Z(J,I)=0.0D0
  260 CONTINUE
  190 CONTINUE
      return
      END


c***********************************************************************
c      SUBROUTINE DAMPIG(r,b,MMAX,DM,DMP,DMPP)
c** Subroutine to generate values DM(m) and the first and second radial
c  derivatives DMP(m) and DMPP(m) of normalized incomplete gamma 
c  functions of orders m=0 to MMAX, at radial distance 'r', for damping
c  parameter 'b'.  NOTE that  DM(m)= {Tang-Toennies function}(m+1).
c***********************************************************************
c     INTEGER MMAX,I,m
c     REAL*8 b,r,br,XP,SSm,SSm1,SSm2,TK,DM(MMAX),DMP(MMAX),DMPP(MMAX)
c     br= b*r
c     XP= DEXP(-br)
c     SSm= 0.d0
c     SSm1= 0.d0
c     SSm2= 0.d0
c     TK= 1.d0
c     DO  m=1, MMAX
c         SSm2= SSm1
c         SSm1= SSm
c         SSm= SSm+ TK
c         DM(m)= 1.d0 - XP*SSM
c         DMP(m)= b*XP*(SSm - SSm1)
c         DMPP(m)= b**2 *XP*(2.d0*SSm1 - SSm - SSm2)
c         TK= TK*br/DFLOAT(m)
c         ENDDO
c     RETURN
c     END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DAMPIG(r,b,MMAX,DM,DMP,DMPP)
c** Subroutine to generate values DM and the first and second radial
c  derivatives DMP and DMPP of the Tang-Toennies/incomplete gamma 
c  function damping function of all orders up to NMAX at radial distance
c   'r' for damping parameter 'b'
c***********************************************************************
      INTEGER MMAX,m,n
      REAL*8 b,r,br,XP,SSm,SSm1,SSm2,mfact,TK,DM(MMAX),DMP(MMAX),
     1       DMPP(MMAX),SM(-1:MMAX)

      br= b*r
      XP= DEXP(-br)
      SSm= 0.d0
      SSm1= 0.d0
      SSm2= 0.d0
      TK= 1.d0
      DO m=1,MMAX
         SSm2= SSm1
         SSm1= SSm
         SSm= SSm + TK
         DM(m)= 1.d0 - XP*SSm
         DMP(m)= b*XP*(SSm - SSm1)
         DMPP(m)= b**2 *XP* (2.d0*SSm1 - SSm - SSm2)
         TK= TK*br/DFLOAT(m)
      ENDDO
      IF(DM(MMAX).LT.1.0D-13) THEN
         mfact= 1.d0
         DO n=1,MMAX
            mfact= mfact*DFLOAT(n)
         ENDDO
         SSm= 0.d0
         TK= (br)**MMAX/mfact
         DO n=MMAX,MMAX+3
            SSm= SSm+TK
            TK= TK*br/DFLOAT(n+1)
         ENDDO
         SM(MMAX)= SSm
         TK= (br)**MMAX/mfact
         DO n=1,MMAX-1
            TK= TK*DFLOAT(MMAX+1-n)/br
            SM(MMAX-n)= TK+SM(MMAX-n+1)
         ENDDO
         SM(0)=1.d0 + SM(1)
         SM(-1)=SM(0)
         DO m=1,MMAX
            DM(m)= XP*SM(m)
            DMP(m)= -b*XP*(SM(m)-SM(m-1))
            DMPP(m)= b**2 *XP*(SM(m)-2.d0*SM(m-1)+SM(m-2))
         ENDDO
      ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12




