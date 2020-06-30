      SUBROUTINE SONNE(DJ,EQUIN,uvw)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      dimension uvw(3)

C  Program zum Berechnen der rechtwinkligen Sonnenkoordinaten
C  nach der Theorie von NEWCOMB, fuer den zeitpunkt DJ (Ephemeridenzeit).
C  Die gewuenschte Normalepoche ist als Julianisches Datum einzugeben.
C  (JD=2415020.313+365.24219878*(EP-1900)-1.228e-5*(EP-1900)**2).

C     MKR98xxxx ??? zulässig
      IF(EQUIN.EQ.1950.D0) THEN
         EPS=23.44578787D0
         IPREC=1
      ELSEIF(EQUIN.EQ.2000.D0) THEN
         EPS=23.43929111D0
         IPREC=2
      END IF

      T=(DJ-2415020.D0)/36525.D0
      !EPS=23.45229444D0-1.30125D-2*T-1.64D-6*T**2+5.03D-7*T**3
      EPS=EPS*DATAN(1.D0)/45.D0
      CALL IAUSUN(DJ,SR,SL,SB,DL,DB,DR)
      CALL POLREC(SL,SB,SR,U,VE,WE)
      CALL EKLAEQ(VE,WE,EPS,V,W)
      CALL RECPOL(U,V,W,SR0,D1,R1)
      CALL PRAEZ(IPREC,(DJ-2415020.313D0)/365.242192D0+1900.D0,
     &R1,D1,EQUIN,R2,D2)
      CALL POLREC(R2,D2,SR,U,V,W)
      uvw(1)=u
      uvw(2)=v
      uvw(3)=w
      return
      IF(NR.EQ.0)THEN
      DAT1=DATUM(DJ)
      DJ1=DJ
      U1=U
      V1=V
      W1=W
      RETURN
      END IF
      WRITE(11,10)DAT1,DJ1,U1,V1,W1,DATUM(DJ),DJ,U,V,W
 10   FORMAT(1X,2(F11.2,2X,F10.2,3(2X,F10.7),8X))
      RETURN
      END

      SUBROUTINE IAUSUN (JD,SR,SL,SB,DBL,DBB,DBR)
C ************************************************************************
C
C  NEWCOMBs Theorie der Erde
C  Streng entsprechend den 'Tables of the Sun'(A.P.A.E. VI Pt.1),
C  ohne irgendwelche Korrekturen.     Argument JD (Ephemeridenzeit)
C
C  Defintionen:
C  JD    = Eingabeargument in Ephemeridenzeit, julianische Tageszahl.
C  EL( 2)= numerische Exzentrizitaet, ungestoert
C  EL( 5)= Laenge des Perihels, ungestoert
C  EL( 9)= ungestoerte mittlere Laenge, in  rad
C  EL(10)= mittlere Anomalie, einschliesslich langperiodische Stoerungen
C  EL(11)= ungestoerter Radius, in AE
C  EL(12)= Mittelpunktsgleichung, in rad
C  EL(14)= mittlere Laenge des Mondes
C  EL(15)= mittlere Anomalie des Mondes
C  EL(16)= Argument der Breite des Mondes
C  EL(17)= mittlere Laenge des aufsteigenden Mondknotens
C  EL(21)= mittlere Laenge des Mondperigaeums
C
C  SR    = gestoerte Distanz Erde-Sonne, in AE
C  SL    = gestoerte Sonnenlaenge in rad, mittl. Aequinoktium des Datums
C  SB    = gestoerte Breite der Sonne in rad, m. Aequinoktium des Datums
C
C ***********************************************************************
C
      IMPLICIT REAL(A-H,O-Z)
      DOUBLE PRECISION JD,T,TP,R0,GME,GV,GMA,GJ,GS,DL,DR,DBLARG,DB,DG
      DOUBLE PRECISION D,T18,EL(21),ZWEIPI,STR,SR,SL,SB,DBL,DBB,DBR
      DIMENSION X(8,121),X1(160),X2(160),X3(160),X4(160),X5(160),X6(8)
      DIMENSION XA(8,10),XB(8,10),X1A(80),X1B(80),X2A(80),X2B(80)
     &,X3A(80),X3B(80),X4A(80),X4B(80),X5A(80),X5B(80)
      EQUIVALENCE (X(1, 21),X1(1))
      EQUIVALENCE (X(1, 41),X2(1))
      EQUIVALENCE (X(1, 61),X3(1))
      EQUIVALENCE (X(1, 81),X4(1))
      EQUIVALENCE (X(1,101),X5(1))
      EQUIVALENCE (X(1,121),X6(1))
      DATA ZWEIPI,STR/6.283185307179586D0,206264806.2470964D0/

C  Teil 1:  Tafel der Stoerungsterme:
C           Stoerungen durch Merkur
C                 J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA XA/  - 1.,+ 1., -   6.,-  11., +  26.,-  12.,     0.,    0.,
     &          - 1.,+ 2., -   3.,-   3., -   4.,+   5.,     0.,    0.,
     &          - 1.,+ 3., +  15.,-   1., -   1.,-  18.,     0.,    0.,
     &          - 1.,+ 4., +  19.,-  13., -   2.,-   4.,     0.,    0.,
     &          - 1.,  0., +  33.,-  67., -  85.,-  39., +  24.,-  17.,
     &          - 1.,+ 1., +2353.,-4228., -2062.,-1146., -   4.,+   3.,
     &          - 1.,+ 2., -  65.,-  34., +  68.,-  14., +   6.,-  92.,
     &          - 1.,+ 3., -   3.,-   8., +  14.,-   8., +   1.,+   7.,
     &          - 2.,  0.,     0.,    0.,     0.,+   4.,     0.,    0.,
     &          - 2.,+ 1., -  99.,+  60., +  84.,+ 136., +  23.,-   3./
      DATA XB/  - 2.,+ 2., -4702.,+2903., +3593.,+5822., +  10.,-   6.,
     &          - 2.,+ 3., +1795.,-1737., - 596.,- 632., +  37.,-  56.,
     &          - 2.,+ 4., +  30.,-  33., +  40.,+  33., +   5.,-  13.,
     &          - 3.,+ 2., -  13.,    0.,     0.,+  21., +  13.,+   5.,
     &          - 3.,+ 3., - 666.,+  27., +  44.,+1044., +   8.,+   1.,
     &          - 3.,+ 4., +1508.,- 397., - 381.,-1448., + 185.,- 100.,
     &          - 3.,+ 5., + 763.,- 684., + 126.,+ 148., +   6.,-   3.,
     &          - 3.,+ 6., +  12.,-  12., +  14.,+  13., -   2.,+   4.,
     &          - 4.,+ 3.,     0.,    0.,     0.,+   6., +   4.,+   5.,
     &          - 4.,+ 4., - 188.,-  93., - 166.,+ 337.,     0.,    0./
      DATA X1A/ - 4.,+ 5., - 139.,-  38., -  51.,+ 189., -  31.,-   1.,
     &          - 4.,+ 6., + 146.,-  42., -  25.,-  91., +  12.,    0.,
     &          - 4.,+ 7., +   5.,-   4., +   3.,+   5.,     0.,    0.,
     &          - 5.,+ 5., -  47.,-  69., - 134.,+  93.,     0.,    0.,
     &          - 5.,+ 6., -  28.,-  25., -  39.,+  43., -   8.,-   4.,
     &          - 5.,+ 7., - 119.,-  33., -  37.,+ 136., -  18.,-   6.,
     &          - 5.,+ 8., + 154.,    0.,     0.,-  26.,     0.,    0.,
     &          - 6.,+ 5.,     0.,    0.,     0.,    0., -   2.,+   6.,
     &          - 6.,+ 6., -   4.,-  38., -  80.,+   8.,     0.,    0.,
     &          - 6.,+ 7., -   4.,-  13., -  24.,+   7., -   2.,-   3./
      DATA X1B/ - 6.,+ 8., -   6.,-   7., -  10.,+  10., -   2.,-   3.,
     &          - 6.,+ 9., +  14.,+   3., +   3.,-  12.,     0.,    0.,
     &          - 7.,+ 7., +   8.,-  18., -  38.,-  17.,     0.,    0.,
     &          - 7.,+ 8., +   1.,-   6., -  12.,-   3.,     0.,    0.,
     &          - 7.,+ 9., +   1.,-   3., -   4.,+   1.,     0.,    0.,
     &          - 7.,+10., -   2.,-   2., -   3.,+   3.,     0.,    0.,
     &          - 8.,+ 8., +   9.,-   7., -  14.,-  19.,     0.,    0.,
     &          - 8.,+ 9., +   2.,-   3., -   5.,-   4.,     0.,    0.,
     &          - 8.,+12., -  33.,-  54., -  43.,+   8., -   5.,-   9.,
     &          - 8.,+13.,     0.,    0., -   9.,-   8.,     0.,    0./
      DATA X2A/ - 8.,+14.,     0.,    0., -  25.,+  22.,     0.,    0.,
     &          - 9.,+ 9., +   6.,-   1., -   2.,-  13.,     0.,    0.,
     &          - 9.,+10.,     0.,    0., -   1.,-   4.,     0.,    0.,
     &          -10.,+10.,     0.,    0., +   3.,-   7.,     0.,    0.,
     &          + 1.,- 2., -   5.,-   4., -   5.,+   6.,     0.,    0.,
     &          + 1.,- 1., - 216.,- 167., -  92.,+ 119.,     0.,    0.,
     &          + 1.,  0., -   8.,-  47., +  27.,-   6.,     0.,    0.,
     &          + 2.,- 3., +  40.,-  10., -  13.,-  50.,     0.,    0.,
     &          + 2.,- 2., +1963.,- 567., - 573.,-1976.,     0.,-   8.,
     &          + 2.,- 1., -1659.,- 617., +  64.,- 137.,     0.,    0./
      DATA X2B/ + 2.,  0., -  24.,+  15., -  18.,-  25., -   8.,+   2.,
     &          + 3.,- 4., +   1.,-   4., -   6.,    0.,     0.,    0.,
     &          + 3.,- 3., +  53.,- 118., - 154.,-  67.,     0.,    0.,
     &          + 3.,- 2., + 396.,- 153., -  77.,- 201.,     0.,    0.,
     &          + 3.,- 1., +   8.,+   1.,     0.,+   6.,     0.,    0.,
     &          + 4.,- 4., +  11.,+  32., +  46.,-  17.,     0.,    0.,
     &          + 4.,- 3., - 131.,+ 483., + 461.,+ 125., +   7.,+   1.,
     &          + 4.,- 2., + 526.,- 256., +  43.,+  96.,     0.,    0.,
     &          + 4.,- 1., +   7.,-   5., +   6.,+   8.,     0.,    0.,
     &          + 5.,- 5., -   7.,+   1.,     0.,+  12.,     0.,    0./
      DATA X3A/ + 5.,- 4., +  49.,+  69., +  87.,-  62.,     0.,    0.,
     &          + 5.,- 3., -  38.,+ 200., +  87.,+  17.,     0.,    0.,
     &          + 5.,- 2., +   3.,+   1., -   1.,+   3.,     0.,    0.,
     &          + 6.,- 6.,     0.,    0., -   4.,-   3.,     0.,    0.,
     &          + 6.,- 5., -  20.,-   2., -   3.,+  30.,     0.,    0.,
     &          + 6.,- 4., - 104.,- 113., - 102.,+  94.,     0.,    0.,
     &          + 6.,- 3., -  11.,+ 100., -  27.,-   4.,     0.,    0.,
     &          + 7.,- 6., +   3.,-   5., -   9.,-   5.,     0.,    0.,
     &          + 7.,- 5., -  49.,+   3., +   4.,+  60.,     0.,    0.,
     &          + 7.,- 4., -  78.,-  72., -  26.,+  28.,     0.,    0./
      DATA X3B/ + 8.,- 7., +   1.,+   3., +   5.,-   1.,     0.,    0.,
     &          + 8.,- 6., +   6.,-   8., -  12.,-   9.,     0.,    0.,
     &          + 8.,- 5., +  51.,-  10., -   8.,-  44.,     0.,    0.,
     &          + 8.,- 4., -  17.,-  12., +   5.,-   6.,     0.,    0.,
     &          + 9.,- 7., +   2.,+   3., +   5.,-   3.,     0.,    0.,
     &          + 9.,- 6., +  13.,-  25., -  30.,-  16.,     0.,    0.,
     &          + 9.,- 5., +  60.,-  15., -   4.,-  17.,     0.,    0.,
     &          +10.,- 7., +   2.,+   5., +   7.,-   3.,     0.,    0.,
     &          +10.,- 6., -   7.,+  18., +  14.,+   6.,     0.,    0.,
     &          +10.,- 5., +   5.,-   2.,     0.,    0.,     0.,    0./
      DATA X4A/ +11.,- 7., +   9.,+  15., +  17.,-  10.,     0.,    0.,
     &          +11.,- 6., -  12.,+  42., +   8.,+   3.,     0.,    0.,
     &          +12.,- 7., -   4.,-   5., -   4.,+   3.,     0.,    0.,
     &          +13.,- 8., -  13.,-   1., -   1.,+  15.,     0.,    0.,
     &          +13.,- 7., -  30.,-  33., -   4.,+   3.,     0.,    0.,
     &          +15.,- 9., +  13.,-  16., -  17.,-  14.,     0.,    0.,
     &          +15.,- 8., + 200.,-  30., -   1.,-   6.,     0.,    0.,
     &          +17.,-10., -   2.,-   4., -   4.,+   2.,     0.,    0.,
     &          +17.,- 9., -  10.,+  24.,     0.,    0.,     0.,    0.,
     &          + 1.,- 3., -   3.,-   1., -   2.,+   5.,     0.,    0./
      DATA X4B/ + 1.,- 2., - 155.,-  52., -  78.,+ 193., +   7.,    0.,
     &          + 1.,- 1., -7208.,+  59., +  56.,+7067., -   1.,+  17.,
     &          + 1.,  0., - 307.,-2582., + 227.,-  89., +  16.,    0.,
     &          + 1.,+ 1., +   8.,-  73., +  79.,+   9., +   1.,+  23.,
     &          + 2.,- 3., +  11.,+  68., + 102.,-  17.,     0.,    0.,
     &          + 2.,- 2., + 136.,+2728., +4021.,- 203.,     0.,    0.,
     &          + 2.,- 1., - 537.,+1518., +1376.,+ 486., +  13.,+ 166.,
     &          + 2.,  0., -  22.,-  70., -   1.,-   8.,     0.,    0.,
     &          + 3.,- 4., -   5.,+   2., +   3.,+   8.,     0.,    0.,
     &          + 3.,- 3., - 162.,+  27., +  43.,+ 278.,     0.,    0./
      DATA X5A/ + 3.,- 2., +  71.,+ 551., + 796.,- 104., +   6.,-   1.,
     &          + 3.,- 1., -  31.,+ 208., + 172.,+  26., +   1.,+  18.,
     &          + 4.,- 4., -   3.,-  16., -  29.,+   5.,     0.,    0.,
     &          + 4.,- 3., -  43.,+   9., +  13.,+  73.,     0.,    0.,
     &          + 4.,- 2., +  17.,+  78., + 110.,-  24.,     0.,    0.,
     &          + 4.,- 1., -   1.,+  23., +  17.,+   1.,     0.,    0.,
     &          + 5.,- 5.,     0.,    0., -   1.,-   3.,     0.,    0.,
     &          + 5.,- 4., -   1.,-   5., -  10.,+   2.,     0.,    0.,
     &          + 5.,- 3., -   7.,+   2., +   3.,+  12.,     0.,    0.,
     &          + 5.,- 2., +   3.,+   9., +  13.,-   4.,     0.,    0./
      DATA X5B/ + 1.,- 2., -   3.,+  11., +  15.,+   3.,     0.,    0.,
     &          + 1.,- 1., -  77.,+ 412., + 422.,+  79., +   1.,+   6.,
     &          + 1.,  0., -   3.,- 320., +   8.,-   1.,     0.,    0.,
     &          + 1.,+ 1.,     0.,-   8., +   8.,    0., -   1.,+   6.,
     &          + 2.,- 3.,     0.,    0., -   3.,-   1.,     0.,    0.,
     &          + 2.,- 2., +  38.,- 101., - 152.,-  57.,     0.,    0.,
     &          + 2.,- 1., +  45.,- 103., - 103.,-  44.,     0.,    0.,
     &          + 2.,  0., +   2.,-  17.,     0.,    0.,     0.,    0.,
     &          + 3.,- 2., +   7.,-  20., -  30.,-  11.,     0.,    0.,
     &          + 3.,- 1., +   6.,-  16., -  16.,-   6.,     0.,    0./
      DATA X6/  + 4.,- 2., +   1.,-   3., -   4.,-   1.,     0.,    0./
C
C  Feldzuweisungen:
      DO 1988 I=1,80
      J=80+I
      X1(I)=X1A(I)
      X1(J)=X1B(I)
      X2(I)=X2A(I)
      X2(J)=X2B(I)
      X3(I)=X3A(I)
      X3(J)=X3B(I)
      X4(I)=X4A(I)
      X4(J)=X4B(I)
      X5(I)=X5A(I)
      X5(J)=X5B(I)
 1988 CONTINUE
      DO 1987 I=1,10
      K=10+I
      DO 1986 J=1,8
      X(J,I)=XA(J,I)
      X(J,K)=XB(J,I)
 1986 CONTINUE
 1987 CONTINUE


C  Teil 2: Berechnung von Hilfsargumenten.
C  T    = Zeit in Julian. Jahrhunderten von 1900 Januar 0 Mittag Greenwich.
      T=(JD-2415020.D0)/36525.D0

C  TP   = Zeit in Julian. Jahren        von 1850 Januar 0 Mittag Greenwich.
      TP=(JD-2396758.D0)/365.25D0

C  T18  = Zeit in Julian. Jahrhunderten von 1800 Januar 0 Mittag Greenwich.
      T18=(JD-2378496.D0)/36525.D0

C  Teil 3: Berechnung der ungestoerten Groessen.
      EL( 2)=0.1675104D-1-T*(0.4180D-4+0.126D-6*T)
      EL( 5)=0.4908229466868880D+1 + (0.3000526416797348D-1
     &     +(0.7902463002085429D-5 +  0.5817764173314427D-7*T)*T)*T
      EL( 9)=0.4881627934111871D+1 + (0.6283319509909086D+3
     &     + 0.52796209872828420D-5*T)*T
      EL(10)=0.6256583580497099D+1 + (0.6283019457267408D+3
     &     -(0.2617993877991491D-5 +  0.5817764173314427D-7*T)*T)*T
C  Reduktion auf ersten KreisumfangFANG
      EL( 9)=DMOD(EL( 9),ZWEIPI)+ZWEIPI
      EL(10)=DMOD(EL(10),ZWEIPI)
C
C
C  Teil 4: Mondstoerungen
C  Hauptargumente nach HANSENs Mondtheorie, mit NEWCOMBs Korrekturen
C
      EL(14)= 0.5859485105530519D+1 + (0.8399709200339593D+4
     &      +(0.4619304753611657D-4 +  0.6544984694978728D-7*T18)
     &      *T18)*T18
      EL(17)= 0.5807638839584459D0  - (0.33757282382233870D+2
     &      -(0.3964806284113784D-4 +  0.3470781143063166D-7*T18)
     &      *T18)*T18
      EL(21)= 0.3933938487925745D+1 + (0.7101833330913285D+2
     &      -(0.1752456013106637D-3 -  0.1773109075921903D-6*T18)
     &      *T18)*T18
      EL(14)=DMOD(EL(14),ZWEIPI)
      EL(17)=DMOD(EL(17),ZWEIPI)
      EL(21)=DMOD(EL(21),ZWEIPI)
C
C  Abgeleitete Argumente:
      EL(15)=EL(14)-EL(21)
      EL(16)=EL(14)-EL(17)
      EL(15)=DMOD(EL(15),ZWEIPI)
      EL(16)=DMOD(EL(16),ZWEIPI)
C
C  Mittlere Elongation:
      D=EL(14)-EL(9)
C
C  Aufsummierung der Stoerterme durch den Mond:
      D=DMOD(D,ZWEIPI)
      ARG=D
      DL=+6469.*SIN(ARG)+13.*SIN(3.*ARG)
      DR=+13390.*COS(ARG)+30.*COS(3.*ARG)
      DBLARG=D+EL(15)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DL=DL+177.*SIN(ARG)
      DR=DR+370.*COS(ARG)
      DBLARG=D-EL(15)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DL=DL-424.*SIN(ARG)
      DR=DR-1330.*COS(ARG)
      DBLARG=3.D0*D-EL(15)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DL=DL+39.*SIN(ARG)
      DR=DR+80.*COS(ARG)
      DBLARG=D+EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DL=DL-64.*SIN(ARG)
      DR=DR-140.*COS(ARG)
      DBLARG=D-EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DL=DL+172.*SIN(ARG)
      DR=DR+360.*COS(ARG)
      EL(16)=DMOD(EL(16),ZWEIPI)
      ARG=EL(16)
      DB=+576.*SIN(ARG)
      DBLARG=EL(16)+EL(15)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DB=DB+16.*SIN(ARG)
      DBLARG=EL(16)-EL(15)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DB=DB-47.*SIN(ARG)
      DBLARG=EL(16)-2.D0*(EL(9)-EL(17))
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      DB=DB+21.*SIN(ARG)

C  Koordinaten des Baryzentrums Erde-Mond
      DBL=-DL/STR
      DBB=-DB/STR
      DBR=-DR

C  Teil 5:  Berechnung der Mittelpunktsgleichung und Planetenstoerungen
C           Mittlere Anomalien der Planeten:
      GME = 0.43296383D+01+0.2608784648D+02*TP
      GV  = 0.19984020D+01+0.1021322923D+02*TP
      GMA = 0.19173489D+01+0.3340556174D+01*TP
      GJ  = 0.25836283D+01+0.5296346478D+00*TP
      GS  = 0.49692316D+01+0.2132432808D+00*TP
C  Reduktion auf den ersten Kreisumfang
      GME=DMOD(GME,ZWEIPI)
      GV =DMOD(GV ,ZWEIPI)
      GMA=DMOD(GMA,ZWEIPI)
      GJ =DMOD(GJ ,ZWEIPI)
      GS =DMOD(GS ,ZWEIPI)
C  Anbringung der grossen Jupiter-Saturn-Ungleichung
C  an Jupiters mittlerer Anomalie
      GJ =GJ +0.579904067D-2*DSIN(5.D0*GS-2.D0*GJ+1.1719644977D0
     &       -0.397401726D-3*TP)
C  Langperiodische Stoerungen der Mittleren Anomalie der Erde
C  ST ist die Zeit in einfacher Genauigkeit.
      ST=T
C  Argument (4*Mars-7*Erde+3*Venus) in zweiter Ordnung der Massen
C  Argument (3*Jupiter-8*Mars+4*Erde) in zweiter Ordnung der Massen
C  Argument (13*Erde-8*Venus) in erster Ordnung der Massen
C  und fuenfter Ordnung in den Exzentrizitaeten und Neigungen
      DG=266.*SIN(0.555015+2.076942*ST)
     &   +6400.*SIN(4.035027+0.3525565*ST)
     &   +(1882.-16.*T)*SIN(0.9990265+2.622706*ST)
      GJ=DMOD(GJ,ZWEIPI)
C  Berechnung der Mittelpunktsgleichung aus der gestoerten mittleren Anomalie
      EL(10)=DG/STR + EL(10)
      EL(10)=DMOD(EL(10),ZWEIPI)
      EL(12)= DSIN(     EL(10))*(6910057.-(17240.+52.*T)*T)
     &       +DSIN(2.D0*EL(10))*(  72338.-   361.*T)
     &       +DSIN(3.D0*EL(10))*(   1054.-     1.*T)
     &       +DSIN(4.D0*EL(10))*      18.
C  DIE ungestoerte Distanz Erde-Sonne
      R0    =                      30570.-   150.*T
     &       -DCOS(     EL(10))*(7274120.-(18140.+50.*T)*T)
     &       -DCOS(2.D0*EL(10))*(  91380.-   460.*T)
     &       -DCOS(3.D0*EL(10))*(   1450.-    10.*T)
     &       -DCOS(4.D0*EL(10))*      20.
      EL(11)=10.D0**(R0*1.D-9)
C  Stoerungen von Laenge und Distanz durch Merkur
      DO 10 K=1,4
C  Argument J*Merkur + I*Erde
      DBLARG=X(1,K)*GME + X(2,K)*EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      CS=COS(ARG)
      SS=SIN(ARG)
      DL=(X(3,K)*CS + X(4,K)*SS)+DL
      DR=(X(5,K)*CS + X(6,K)*SS)+DR
 10   CONTINUE
C  Stoerungen durch Venus
      DO 20 K=5,44
C  Argument J*Venus+ I*Erde
      DBLARG=X(1,K)*GV + X(2,K)*EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      CS=COS(ARG)
      SS=SIN(ARG)
      DL=(X(3,K)*CS + X(4,K)*SS)+DL
      DR=(X(5,K)*CS + X(6,K)*SS)+DR
      DB=(X(7,K)*CS + X(8,K)*SS)+DB
 20   CONTINUE
C  Stoerungen durch Mars
      DO 30 K=45,89
C  Argument J*Mars + I*Erde
      DBLARG=X(1,K)*GMA+X(2,K)*EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      CS=COS(ARG)
      SS=SIN(ARG)
      DL=(X(3,K)*CS + X(4,K)*SS)+DL
      DR=(X(5,K)*CS + X(6,K)*SS)+DR
      DB=(X(7,K)*CS + X(8,K)*SS)+DB
 30   CONTINUE
C  Stoerungen durch Jupiter
      DO 40 K=90,110
C  Argument J*Jupiter + I*Erde
      DBLARG=X(1,K)*GJ+X(2,K)*EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      CS=COS(ARG)
      SS=SIN(ARG)
      DL=(X(3,K)*CS + X(4,K)*SS)+DL
      DR=(X(5,K)*CS + X(6,K)*SS)+DR
      DB=(X(7,K)*CS + X(8,K)*SS)+DB
 40   CONTINUE
C  Stoerungen durch Saturn
      DO 50 K=111,121
C  Argument J*Saturn + I*Erde
      DBLARG=X(1,K)*GS+X(2,K)*EL(10)
      DBLARG=DMOD(DBLARG,ZWEIPI)
      ARG=DBLARG
      CS=COS(ARG)
      SS=SIN(ARG)
      DL=(X(3,K)*CS + X(4,K)*SS)+DL
      DR=(X(5,K)*CS + X(6,K)*SS)+DR
      DB=(X(7,K)*CS + X(8,K)*SS)+DB
 50   CONTINUE

C  Teil 6:  Berechnung der ekliptikalen Koordinaten

      DBR=EL(11)*(10.D0**(DBR*1.D-9)-1.D0)
      SR=EL(11)*10.D0**(DR*1.D-9)
      SL=DMOD((DL+DG+EL(12))/STR + EL(9),ZWEIPI)
      SB=DB/STR

      RETURN
      END

      SUBROUTINE STERNZ(T0,S)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  UP zur Berechnung der Greenwicher Ortssternzeit S zum Zeitpunkt T0 (UT)
      T=(T0-2415020.D0)/36525.D0
      S=DMOD((8640184.542D0*T+0.0929D0*T*T)/86400.D0+DMOD(T0,1.D0)+0.776
     &9194D0,1.D0)
      S=S*8.D0*DATAN(1.D0)
      RETURN
      END

      SUBROUTINE POLREC(A,B,R,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  UP zum berechnen von rechtwinkligen aus polaren Koordinaten
C  A,B,R = Azimut, Hoehe und Distanz,  X,Y,Z = rechtwinklige Koordinaten
      X=R*DCOS(B)*DCOS(A)
      Y=R*DCOS(B)*DSIN(A)
      Z=R*DSIN(B)
      RETURN
      END

      SUBROUTINE EKLAEQ(Y,Z,EPS,Y0,Z0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  UP zum Umrechnen von ekliptikalen Koordinaten in aequatoriale
C  EPS ist die Neigung der Ekliptik
      Y0=Y*DCOS(EPS)-Z*DSIN(EPS)
      Z0=Y*DSIN(EPS)+Z*DCOS(EPS)
      RETURN
      END

      SUBROUTINE NUTAT(TA,DE,DO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  UP zur Berechnung der Nutation in der Laenge des Aequinoktiums
C  und der Neigung der Ekliptik
      DIMENSION X(9,69),X1(180),X2(180),X3(81)
      DIMENSION XA(9,10),XB(9,10),X1A(90),X1B(90),X2A(90),X2B(90)
      EQUIVALENCE (X(1,21),X1(1)),(X(1,41),X2(1)),(X(1,61),X3(1))
      DATA XA/   0., 0., 0., 0., 1.,-172327.,-173.7,  92100., 9.1,
     &           0., 0., 0., 0., 2.,   2088.,    .2,   -904., 0.4,
     &          -2., 0., 2., 0., 1.,     45.,    .0,    -24., 0.0,
     &           2., 0.,-2., 0., 0.,     10.,    .0,      0., 0.0,
     &           0.,-2., 2.,-2., 1.,     -4.,    .0,      2., 0.0,
     &          -2., 0., 2., 0., 2.,     -3.,    .0,      2., 0.0,
     &           1.,-1., 0.,-1., 0.,     -2.,    .0,      0., 0.0,
     &           0., 0., 2.,-2., 2., -12729.,  -1.3,   5522.,-2.9,
     &           0., 1., 0., 0., 0.,   1261.,  -3.1,      0., 0.0,
     &           0., 1., 2.,-2., 2.,   -497.,   1.2,    216.,-0.6/
      DATA XB/   0.,-1., 2.,-2., 2.,    214.,  -0.5,    -93., 0.3,
     &           0., 0., 2.,-2., 1.,    124.,    .1,    -66., 0.0,
     &           2., 0., 0.,-2., 0.,     45.,    .0,      0., 0.0,
     &           0., 0., 2.,-2., 0.,    -21.,    .0,      0., 0.0,
     &           0., 2., 0., 0., 0.,     16.,  -0.1,      0., 0.0,
     &           0., 1., 0., 0., 1.,    -15.,    .0,      8., 0.0,
     &           0., 2., 2.,-2., 2.,    -15.,    .1,      7., 0.0,
     &           0.,-1., 0., 0., 1.,    -10.,    .0,      5., 0.0,
     &          -2., 0., 0., 2., 1.,     -5.,    .0,      3., 0.0,
     &           0.,-1., 2.,-2., 1.,     -5.,    .0,      3., 0.0/
      DATA X1A/  2., 0., 0.,-2., 1.,      4.,    .0,     -2., 0.0,
     &           0., 1., 2.,-2., 1.,      3.,    .0,     -2., 0.0,
     &           1., 0., 0.,-1., 0.,     -3.,    .0,      0., 0.0,
     &           0., 0., 2., 0., 2.,  -2037.,  -0.2,    884.,-0.5,
     &           1., 0., 0., 0., 0.,    675.,    .1,      0., 0.0,
     &           0., 0., 2., 0., 1.,   -342.,  -0.4,    183., 0.0,
     &           1., 0., 2., 0., 2.,   -261.,    .0,    113.,-0.1,
     &           1., 0., 0.,-2., 0.,   -149.,    .0,      0., 0.0,
     &          -1., 0., 2., 0., 2.,    114.,    .0,    -50., 0.0,
     &           0., 0., 0., 2., 0.,     60.,    .0,      0., 0.0/
      DATA X1B/  1., 0., 0., 0., 1.,     58.,    .0,    -31., 0.0,
     &          -1., 0., 0., 0., 1.,    -57.,    .0,     30., 0.0,
     &          -1., 0., 2., 2., 2.,    -52.,    .0,     22., 0.0,
     &           1., 0., 2., 0., 1.,    -44.,    .0,     23., 0.0,
     &           0., 0., 2., 2., 2.,    -32.,    .0,     14., 0.0,
     &           2., 0., 0., 0., 0.,     28.,    .0,      0., 0.0,
     &           1., 0., 2.,-2., 2.,     26.,    .0,    -11., 0.0,
     &           2., 0., 2., 0., 2.,    -26.,    .0,     11., 0.0,
     &           0., 0., 2., 0., 0.,     25.,    .0,      0., 0.0,
     &          -1., 0., 2., 0., 1.,     19.,    .0,    -10., 0.0/
      DATA X2A/ -1., 0., 0., 2., 1.,     14.,    .0,     -7., 0.0,
     &           1., 0., 0.,-2., 1.,    -13.,    .0,      7., 0.0,
     &          -1., 0., 2., 2., 1.,     -9.,    .0,      5., 0.0,
     &           1., 1., 0.,-2., 0.,     -7.,    .0,      0., 0.0,
     &           0., 1., 2., 0., 2.,      7.,    .0,     -3., 0.0,
     &           1., 0., 0., 2., 0.,      6.,    .0,      0., 0.0,
     &           0., 0., 0., 2., 1.,     -6.,    .0,      3., 0.0,
     &           0.,-1., 2., 0., 2.,     -6.,    .0,      3., 0.0,
     &           1., 0., 2., 2., 2.,     -6.,    .0,      3., 0.0,
     &           2., 0., 2.,-2., 2.,      6.,    .0,     -2., 0.0/
      DATA X2B/  0., 0., 0.,-2., 1.,     -5.,    .0,      3., 0.0,
     &           0., 0., 2., 2., 1.,     -5.,    .0,      3., 0.0,
     &           1., 0., 2.,-2., 1.,      5.,    .0,     -3., 0.0,
     &           0., 0., 0., 1., 0.,     -4.,    .0,      0., 0.0,
     &           0., 1., 0.,-2., 0.,     -4.,    .0,      0., 0.0,
     &           1.,-1., 0., 0., 0.,      4.,    .0,      0., 0.0,
     &           1., 0.,-2., 0., 0.,      4.,    .0,      0., 0.0,
     &           2., 0., 2., 0., 1.,     -4.,    .0,      2., 0.0,
     &           1., 0., 2., 0., 0.,      3.,    .0,      0., 0.0,
     &           1., 1., 0., 0., 0.,     -3.,    .0,      0., 0.0/
      DATA X3/   1.,-1., 2., 0., 2.,     -3.,    .0,      1., 0.0,
     &          -2., 0., 0., 0., 1.,     -2.,    .0,      1., 0.0,
     &          -1., 0., 2.,-2., 1.,     -2.,    .0,      1., 0.0,
     &           2., 0., 0., 0., 1.,      2.,    .0,     -1., 0.0,
     &          -1.,-1., 2., 2., 2.,     -2.,    .0,      1., 0.0,
     &           0.,-1., 2., 2., 2.,     -2.,    .0,      1., 0.0,
     &           1., 0., 0., 0., 2.,     -2.,    .0,      1., 0.0,
     &           1., 1., 2., 0., 2.,      2.,    .0,     -1., 0.0,
     &           3., 0., 2., 0., 2.,     -2.,    .0,      1., 0.0/

C  Feldzuweisungen
      DO 1988 I=1,90
      J=90+I
      X1(I)=X1A(I)
      X1(J)=X1B(I)
      X2(I)=X2A(I)
      X2(J)=X2B(I)
 1988 CONTINUE
      DO 1987 I=1,10
      K=10+I
      DO 1986 J=1,9
      X(J,I)=XA(J,I)
      X(J,K)=XB(J,I)
 1986 CONTINUE
 1987 CONTINUE

      T=TA-2415020.D0
      Z=T/10000.D0
      P=296.104608D0+13.0649924465D0*T+6.89D-4*Z**2+29.5D-8*Z**3
      PI=358.475833D0+0.9856002669D0*T-0.112D-4*Z**2-6.8D-8*Z**3
      O=11.250889D0+13.2293504490D0*T-2.407D-4*Z**2-0.7D-8*Z**3
      W=350.737486D0+12.1907491914D0*T-1.076D-4*Z**2+3.9D-8*Z**3
      Q=259.183275D0-0.0529539222D0*T+1.557D-4*Z**2+4.6D-8*Z**3
      R=8.D0*DATAN(1.D0)/360.D0
      P=DMOD(P,360.D0)*R
      PI=DMOD(PI,360.D0)*R
      O=DMOD(O,360.D0)*R
      W=DMOD(W,360.D0)*R
      Q=DMOD(Q,360.D0)*R
      DE=0.D0
      DO=0.D0
      T=T/36525.D0
      DO 2 I=1,69
      ARG=X(1,I)*P+X(2,I)*PI+X(3,I)*O+X(4,I)*W+X(5,I)*Q
      DE=DE+(X(6,I)+X(7,I)*T)*DSIN(ARG)
      DO=DO+(X(8,I)+X(9,I)*T)*DCOS(ARG)
 2    CONTINUE
      DE=DE/2062648062.470964D0
      DO=DO/2062648062.470964D0
      RETURN
      END


      SUBROUTINE PRAEZ(I,T1,R1,D1,T2,R2,D2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  UP zum Berechnen der Praezession nach den Formeln von ANDOYER (CR 174,506)
C  T1,R1,D1 ; T2,R2,D2 = Zeit (trop. Jahr), Rektaszension und Deklination des
C  alten und neuen Aequinoktiums.  Bei  i=1  Praezessionskonstante nach
C  NEWCOMB, bei i=2  nach LIESKE (FK5).
      DIMENSION EP(2)
      DATA EP,R /1900.D0,2000.D0,206264.8D0/
      S=(T1-EP(I))/100.D0
      T=(T2-T1)/100.D0
      IF(DABS(T).GT.1.D-7)GOTO 1
      R2=R1
      D2=D1
      RETURN
 1    CONTINUE
      IF(I.EQ.2)GOTO 2
C  Praezession nach NEWCOMB (bei FK4 - Koordinaten)
      A=((2304.253+1.3973*S+0.00006*S*S)*T
     &  +(0.3023-0.00027*S)*T*T+0.01800*T**3)/R
      B=((2004.685-0.8533*S-0.00037*S*S)*T
     &  -(0.4267+0.00037*S)*T*T-0.04180*T**3)/R
      C=((4608.506+2.7945*S+0.00012*S*S)*T
     &  +(1.3973+0.00012*S)*T*T+0.03632*T**3)/R
      GOTO 3
 2    CONTINUE
C  Praezession nach LIESKE (bei FK5 - Koordinaten)
      A=((2306.2181+1.39656*S-0.000139*S*S)*T
     &  +(0.30188-0.000345*S)*T*T+0.017998*T**3)/R
      B=((2004.3109-0.85330*S-0.000217*S*S)*T
     &  -(0.42665+0.000217*S)*T*T-0.041833*T**3)/R
      C=((4612.4362+2.79312*S-0.000278*S*S)*T
     &  +(1.39656-0.000279*S)*T*T+0.036201*T**3)/R
 3    CONTINUE
      W=DASIN(DSIN(B)*DSIN(R1+A))/2.D0
      P=DATAN(DTAN(B)*DCOS(R1+A))/2.D0
      F=DATAN(DSIN(2.D0*W)*DTAN(D1+2.D0*P))/2.D0
      R2=R1+(F-DATAN(DTAN(W)*DTAN(P)))*2.D0+C
      D2=D1+(P-DATAN(DTAN(W)*DTAN(F)))*2.D0
      RETURN
      END

      SUBROUTINE RECPOL(U,V,W,R,A,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  UP zum Umformen von rechtwinkligen Koordinaten in Polarkoordinaten
C  U,V,W = rechtwinklige Koordinaten; A=Hoehe,B=Azimut,R=Distanz
      ZWEIPI=8.D0*DATAN(1.D0)
      R=DSQRT(U*U+V*V+W*W)
      A=DATAN(W/DSQRT(U*U+V*V))
      B=DMOD(DATAN2(V,U)+ZWEIPI,ZWEIPI)
      RETURN
      END


      SUBROUTINE YDAT(T,YD)
C Berechnung des Julianischen Datums aus dem Buergerlichen.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      AJ=T/1.D4
      J=INT(AJ)
      AM=(AJ-AINT(AJ))*1.D2
      M=INT(AM)
      D=(AM-AINT(AM))*1.D2
      IF(M.LE.2)THEN
      J=J-1
      M=M+12
      END IF
      B=0.D0
      IF(T.GE.15821015.D0)THEN
      IA=INT(J/100)
      B=2-IA+INT(IA/4)
      END IF
      YD=AINT(365.25D0*J)+AINT(30.6001D0*(M+1))+1720994.5D0+B+D
      RETURN
      END

      DOUBLE PRECISION FUNCTION DATUM(YD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      AN=YD+0.5D0
      IZ=AN
      F=AN-IZ
      IF(IZ.LT.2299161)THEN
      IA=IZ
      GOTO 1
      END IF
      IA0=(IZ-1867216.25D0)/36524.25D0
      IA=IZ+1+IA0-IA0/4
    1 IB=IA+1524
      IC=(IB-122.1D0)/365.25D0
      K=365.25D0*IC
      IE=(IB-K)/30.6001D0
      TAG=IB-K-INT(30.6001D0*IE)+F
      IF(IE.GT.13)THEN
      MONAT=IE-13
      JAHR=IC-4715
      ELSE
      MONAT=IE-1
      JAHR=IC-4716
      END IF
      DATUM=JAHR*1.D4+MONAT*1.D2+TAG
      RETURN
      END

