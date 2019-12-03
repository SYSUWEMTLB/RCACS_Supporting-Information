      MODULE ALLOCRCALFS
       REAL, ALLOCATABLE :: FLXMB(:,:,:,:)
       REAL, ALLOCATABLE :: FLYMB(:,:,:,:)
      END MODULE ALLOCRCALFS
      SUBROUTINE RCALFS
      USE ALLOCRCALFS

C     RCALF IS THE MASTER SUBROUTINE USED TO INTEGRATE THE DIFFERENTIAL
C           EQUATIONS OF THE WATER QUALITY MODEL.
C           IT CONTROLS THE NUMERICAL INTEGRATION USING A CENTRAL
C           DIFFERENCING SCHEME FOR TIME AND A SMOLARKIEWICZ CORRECTED
C           UPWIND SCHEME FOR SPACE.
  
      SAVE
      INCLUDE 'RCACM' 
      INTEGER   NEGSEG(99,3)
      REAL   C(NX,NY,NZ,NOSYS)
     .    ,CB(NX,NY,NZ,NOSYS),CF(NX,NY,NZ,NOSYS)
     .    ,DM(NX,NY,NZ),MF(NX,NY,NZ)
     .    ,BVOLB(NX,NY,NZ),BVOLF(NX,NY,NZ),TOTMAS
     .    ,AL(NX,NY,NZ),AD(NX,NY,NZ),AU(NX,NY,NZ)
     .    ,BETA(NX,NY,NZ),GAMMA(NX,NY,NZ)
 
      REAL   CONCS(NX,NY,NZ),CSLAVE(NX,NY,NZ,NOSYS)
     .      ,DXBAR1(NX,NY),DYBAR1(NX,NY)
     .      ,DXBAR2(NX,NY),DYBAR2(NX,NY)
     .      ,D(NX,NY)
     .      ,ANTIQX(NX,NY,NZ),QXFLX(NX,NY,NZ)
     .      ,ANTIQY(NX,NY,NZ),QYFLX(NX,NY,NZ)
     .      ,ANTIQZ(NX,NY,NZ),QZFLX(NX,NY,NZ)
     .      ,RXFLX(NX,NY,NZ),RYFLX(NX,NY,NZ)
     .      ,U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ)
     .      ,FLUXIN(NX,NY,NZ),FLUXOUT(NX,NY,NZ)
     .      ,SYSMASS(NOSYS+NOKINSYS),SYSLOADS(4,NOSYS)
   
 
C  INITIALIZATION
      INITB=0
      IFLAG_TBRK=1
      IFLAG_TBC=1
      ITRAK=1 
      ITIMESECS=0
      ITIMEWQSECS=0
      TIME=0.0
      NXPRTD=IPRNTDSECS
      NXPRTG=IPRNTGSECS
      IWRTADDM=0
      IADDMDEBUG=0
      INITMB=1
      IF(MASSBAL.EQ.1) THEN
        IF(ISMBSECS.EQ.0) THEN
          NXPRTMB=0
          INITMB=0
        ELSE
          NXPRTMB=ISMBSECS
          INITMB=1
        ENDIF
      ENDIF
      R=1.E-10
C
      IF(TZERO.NE.0.) THEN
         TIME = TZERO
         CALL RCA08
         ITIMESECS = 86400.*TZERO
         ITIMEWQSECS = 86400.*TZERO
         NXPRTD = ITIMESECS+IPRNTDSECS
         NXPRTG = ITIMESECS+IPRNTGSECS
         IF(MASSBAL.EQ.1) THEN
           IF(ITIMESECS.LT.ISMBSECS) THEN
             NXPRTMB=ISMBSECS
             INITMB=1
           ELSE
             NXPRTMB = ITIMESECS
             INITMB=0
           ENDIF
         ENDIF
      ENDIF
C        INITIALIZE ARRAYS
      DO 10 IZ=1,NZ
        DO 10 IY=1,NY
          DO 10 IX=1,NX
            CONCS(IX,IY,IZ) = 0.0
            ANTIQX(IX,IY,IZ) = 0.0
            ANTIQY(IX,IY,IZ) = 0.0
            ANTIQZ(IX,IY,IZ) = 0.0
            QXFLX(IX,IY,IZ) = 0.0
            QYFLX(IX,IY,IZ) = 0.0
            QZFLX(IX,IY,IZ) = 0.0
            RXFLX(IX,IY,IZ) = 0.0
            RYFLX(IX,IY,IZ) = 0.0
   10 CONTINUE
      IF(MASSBAL.EQ.1) THEN
       ALLOCATE(FLXMB(NX,NY,NZ,NOSYS+NOKINSYS),STAT=ISTAT)
       ALLOCATE(FLYMB(NX,NY,NZ,NOSYS+NOKINSYS),STAT=ISTAT)
       DO 12 ISYS=1,NOSYS+NOKINSYS
        DO 12 IZ=1,NZ
         DO 12 IY=1,NY
          DO 12 IX=1,NX
           FLXMB(IX,IY,IZ,ISYS)=0.
           FLYMB(IX,IY,IZ,ISYS)=0.
   12  CONTINUE
       DO 14 ISYS=1,NOSYS
         DO I=1,4
          SYSLOADS(I,ISYS)=0.0
         ENDDO
   14  CONTINUE
       IF(INITMB.EQ.0) THEN
        DO 16 ISYS=1,NOSYS
         SYSMASS(ISYS)=0.0
         DO 16 IZ=1,NZ
          DO 16 IY=1,NY
           DO 16 IX=1,NX
            IF(FSM(IX,IY).EQ.1.)  SYSMASS(ISYS) = SYSMASS(ISYS) +
     .                            BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
  16    CONTINUE
        DO 17 ISYS=1,NOKINSYS
         SYSMASS(NOSYS+ISYS)=0.0
         DO 17 IZ=1,NZ
          DO 17 IY=1,NY
           DO 17 IX=1,NX
            IF(FSM(IX,IY).EQ.1.)
     .                 SYSMASS(NOSYS+ISYS) = SYSMASS(NOSYS+ISYS) +
     .                 BVOL(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS)
  17    CONTINUE
        WRITE(17)  TIME,SYSMASS,SYSLOADS,FLXMB,FLYMB
        NXPRTMB=NXPRTMB+IPRNTMBSECS
        INITMB=1
       ENDIF
      ENDIF

      DO 25 ISYS=1,NOSYS
       DO 25 IZ=1,NZ
        DO 25 IY=1,NY
         DO 25 IX=1,NX
           CB(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)
           C(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)
           CSLAVE(IX,IY,IZ,ISYS) = 0.0
   25 CONTINUE

C        CALCULATE CONSTANTS USED FOR SMOLARKIEWICZ CALCULATIONS
      DO 30 IY=2,NY
       DO 30 IX=2,NX
           DXBAR1(IX,IY) = 0.5*(DX(IX,IY)+DX(IX-1,IY))
           DXBAR2(IX,IY) = 0.5*(DX(IX,IY)+DX(IX,IY-1))
           DYBAR1(IX,IY) = 0.5*(DY(IX,IY)+DY(IX-1,IY))
           DYBAR2(IX,IY) = 0.5*(DY(IX,IY)+DY(IX,IY-1))
   30 CONTINUE

C        SET INITIAL VOLUMES AT TIME 'N-1'
       DO 35 IZ=1,NZ
        DO 35 IY=1,NY
         DO 35 IX=1,NX
          BVOLB(IX,IY,IZ) = BVOL(IX,IY,IZ)
   35 CONTINUE

C        SET UP -AVECT- ARRAYS FOR TIME = 0.0 (TZERO)
      CALL RCAEXP1
C        PRINT INITIAL CONDITIONS
      CALL RCA09
      IDISK = 3 

C  BEGIN LEAP-FROG STEP
  40  CONTINUE

C        GET TOTAL MASS PER SYSTEM IF MASS/FLUX BALANCES REQUESTED
      IF(MASSBAL.EQ.1 .AND. ITIMESECS.EQ.NXPRTMB) THEN
       DO 41 ISYS=1,NOSYS
        SYSMASS(ISYS)=0.0
        DO 41 IZ=1,NZ
         DO 41 IY=1,NY
          DO 41 IX=1,NX
           IF(FSM(IX,IY).EQ.1.)  SYSMASS(ISYS) = SYSMASS(ISYS) +
     .                           BVOLB(IX,IY,IZ)*CB(IX,IY,IZ,ISYS)
  41   CONTINUE
       DO 42 ISYS=1,NOKINSYS
        SYSMASS(NOSYS+ISYS)=0.0
        DO 42 IZ=1,NZ
         DO 42 IY=1,NY
          DO 42 IX=1,NX
           IF(FSM(IX,IY).EQ.1.)
     .                SYSMASS(NOSYS+ISYS) = SYSMASS(NOSYS+ISYS) +
     .                BVOLB(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS)
  42   CONTINUE
      ENDIF


C        IF PIECEWISE LINEAR BOUNDARY CONDITIONS THEN
C           GET CONCENTRATIONS AT CURRENT TIME LEVEL
      IF(IBCPWLOPT.EQ.1) THEN
        DO 50 ISYS=1,NOSYS
          IF(NOBC(ISYS).EQ.0) GO TO 50
          DELTBC = TIME - NXBCT
          DO 45 I=1,NOBC(ISYS)
            IX=IBC(1,I,ISYS)
            IY=IBC(2,I,ISYS)
            IZ=IBC(3,I,ISYS)
            CARAY(IX,IY,IZ,ISYS) = DELTBC*SBC(I,ISYS) + BBC(I,ISYS)
            CB(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)
   45     CONTINUE
   50  CONTINUE
      ENDIF

C        COMPUTE DEPTHS AT TIME LEVEL N
      DELT = TIME - (NXHYDTSECS-IHYDDTSECS)/86400.
      DO 55 IY=1,NY
        DO 55 IX=1,NX
            D(IX,IY) = H(IX,IY)+XAZ(IX,IY)+DELT*DETA(IX,IY)*86400.
  55  CONTINUE

C        CALCULATE HORIZONTAL AND VERTICAL VELOCITIES AT TIME 'N'
      DO 60 IZ=1,NZ
        DO 60 IY=2,NY
          DO 60 IX=2,NX
            U(IX,IY,IZ)=QX(IX,IY,IZ)/(DZ(IZ)*0.5*(D(IX,IY)+D(IX-1,IY))*
     .                  DYBAR1(IX,IY))
            V(IX,IY,IZ)=QY(IX,IY,IZ)/(DZ(IZ)*0.5*(D(IX,IY)+D(IX,IY-1))*
     .                  DXBAR2(IX,IY))
            W(IX,IY,IZ)=QZ(IX,IY,IZ)/XAZ(IX,IY)
 60   CONTINUE

C        EVALUATE KINETIC DERIVATIVES (IF NECESSARY)
      IF(DTWQ.EQ.0.0 .OR. ITIMESECS.GE.ITIMEWQSECS) THEN
        CALL TUNER
        IF(INITB.EQ.0) THEN
          DO ISYS=1,NOSYS
           DO IZ=1,NZ
            DO IY=1,NY
             DO IX=1,NX
              C(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)
              CB(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ENDIF
        IDISK = 0
        ITIMEWQSECS=ITIMEWQSECS+IDTWQSECS
        TIMEWQ=ITIMEWQSECS/86400.
      ENDIF
      INITB = 1
 
C        COMPUTE VOLUMES AT TIME LEVEL 'N+1'
      DO 80 IZ=1,NZ
       DO 80 IY=1,NY
        IS = IXS(IY)
        IE = IXE(IY)
        IF(IS.EQ.0) GO TO 80
        DO 79 IX=IS,IE
          BVOLF(IX,IY,IZ) = BVOLB(IX,IY,IZ) + 2.*DT*VDER(IX,IY,IZ)
  79    CONTINUE
  80  CONTINUE

      DO 320 ISYS=1,NOSYS 
C        CHECK FOR SYSTEM BYPASS
      IF(SYSBY(ISYS).EQ.1) GO TO 320
        INEGS=0
        TOTMAS=0.0
        DO 100 IZ=1,NZ
         DO 100 IY=1,NY
          DO 100 IX=1,NX
           DM(IX,IY,IZ) = 0.0
           MF(IX,IY,IZ) = 0.0
  100   CONTINUE

C        THE 120 THROUGH 131 DO LOOPS TAKE AN EULER INTEGRATION STEP
C        STEP FOR HORIZONTAL AND VERTICAL ADVECTIVE TRANSPORT
C           (-DIAG- AND -AVECT- CONTAIN JUST  Q  TERMS)
 
C        COMPUTE C AT TIME LEVEL N+1
       IF (NZ.GT.1) THEN
C      SURFACE LAYER
        DO 120 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0) GO TO 120
          DO 119 IX=IS,IE
            IF(FSM(IX,IY).LE.0.) GO TO 119
            DM(IX,IY,1) = - DIAG(IX,IY,1)*CB(IX,IY,1,ISYS)
     .        + AVECT(IX,IY,1,1)*CB(IX-1,IY,1,ISYS)
     .         + AVECT(IX,IY,1,2)*CB(IX+1,IY,1,ISYS)
     .          + AVECT(IX,IY,1,3)*CB(IX,IY-1,1,ISYS)
     .           + AVECT(IX,IY,1,4)*CB(IX,IY+1,1,ISYS)
     .            + AVECT(IX,IY,1,6)*CB(IX,IY,2,ISYS)
  119     CONTINUE
  120   CONTINUE
C      LAYERS 2 TO NZ-1
        DO 125 IZ=2,NZ-1
          DO 125 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 125
            DO 124 IX=IS,IE
              IF(FSM(IX,IY).LE.0.) GO TO 124
              DM(IX,IY,IZ) = - DIAG(IX,IY,IZ)*CB(IX,IY,IZ,ISYS)
     .         + AVECT(IX,IY,IZ,1)*CB(IX-1,IY,IZ,ISYS)
     .          + AVECT(IX,IY,IZ,2)*CB(IX+1,IY,IZ,ISYS)
     .           + AVECT(IX,IY,IZ,3)*CB(IX,IY-1,IZ,ISYS)
     .            + AVECT(IX,IY,IZ,4)*CB(IX,IY+1,IZ,ISYS)
     .             + AVECT(IX,IY,IZ,5)*CB(IX,IY,IZ-1,ISYS)
     .              + AVECT(IX,IY,IZ,6)*CB(IX,IY,IZ+1,ISYS)
  124       CONTINUE
  125   CONTINUE
C      BOTTOM LAYER
         DO 130 IY=1,NY
           IS = IXS(IY)
           IE = IXE(IY)
           IF(IS.EQ.0) GO TO 130
           DO 129 IX=IS,IE
             IF(FSM(IX,IY).LE.0.) GO TO 129
             DM(IX,IY,NZ) =
     .        - DIAG(IX,IY,NZ)*CB(IX,IY,NZ,ISYS)
     .         + AVECT(IX,IY,NZ,1)*CB(IX-1,IY,NZ,ISYS)
     .          + AVECT(IX,IY,NZ,2)*CB(IX+1,IY,NZ,ISYS)
     .           + AVECT(IX,IY,NZ,3)*CB(IX,IY-1,NZ,ISYS)
     .            + AVECT(IX,IY,NZ,4)*CB(IX,IY+1,NZ,ISYS)
     .             + AVECT(IX,IY,NZ,5)*CB(IX,IY,NZ-1,ISYS)
  129      CONTINUE
  130   CONTINUE

       ELSE

        IZ=1
        DO 132 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0) GO TO 132
          DO 131 IX=IS,IE
            IF(FSM(IX,IY).LE.0.) GO TO 131
            DM(IX,IY,IZ) = - DIAG(IX,IY,IZ)*CB(IX,IY,IZ,ISYS)
     .       + AVECT(IX,IY,IZ,1)*CB(IX-1,IY,IZ,ISYS)
     .        + AVECT(IX,IY,IZ,2)*CB(IX+1,IY,IZ,ISYS)
     .         + AVECT(IX,IY,IZ,3)*CB(IX,IY-1,IZ,ISYS)
     .          + AVECT(IX,IY,IZ,4)*CB(IX,IY+1,IZ,ISYS)
  131     CONTINUE
  132   CONTINUE

       ENDIF
 
C        ADD SOURCES/SINKS DUE TO KINETICS
        DO 135 IZ=1,NZ
          DO 135 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 135
            DO 134 IX=IS,IE
              DM(IX,IY,IZ) = DM(IX,IY,IZ) + CDARAY(IX,IY,IZ,ISYS)
 134        CONTINUE
 135    CONTINUE

C        SET MASS BALANCE FLAG FOR LOADS
       IMBWK=0
       IF(MASSBAL.EQ.1 .AND. ITIMESECS.GE.ISMBSECS) THEN
        IF(IMBDOPT.EQ.0) THEN
         IF(ITIMESECS.EQ.NXPRTMB) IMBWK=1
         DTMB=1.
        ELSE
         IMBWK=1
         DTMB=DT
        ENDIF
       ENDIF

C        ADD LOADS
        IF(NOPS(ISYS).GT.0) THEN
          DELTPS = TIME - NXPST
          DO 140 I=1,NOPS(ISYS) 
           WK = DELTPS*SPS(I,ISYS) + BPS(I,ISYS)
           IF(IMBWK.EQ.1)  SYSLOADS(1,ISYS) = SYSLOADS(1,ISYS)+DTMB*WK
           IX = IPS(1,I,ISYS)
           IY = IPS(2,I,ISYS)
           IWK = IPS(3,I,ISYS)
           DO IZ=1,NZ
            DM(IX,IY,IZ) = DM(IX,IY,IZ) + ZFRACPS(IZ,IWK)*WK
           ENDDO
  140     CONTINUE
        ENDIF
        IF(NONPS(ISYS).GT.0) THEN
          DELTNPS = TIME - NXNPST
          DO 145 I=1,NONPS(ISYS)
           WK = DELTNPS*SNPS(I,ISYS) + BNPS(I,ISYS)
           IF(IMBWK.EQ.1)  SYSLOADS(2,ISYS) = SYSLOADS(2,ISYS)+DTMB*WK
           IX = INPS(1,I,ISYS)
           IY = INPS(2,I,ISYS)
           IWK = INPS(3,I,ISYS)
           DO IZ=1,NZ
            DM(IX,IY,IZ) = DM(IX,IY,IZ) + ZFRACNPS(IZ,IWK)*WK
           ENDDO
  145     CONTINUE
        ENDIF
        IF(NOFL(ISYS).GT.0) THEN
          DELTFL = TIME - NXFLT
          DO 150 I=1,NOFL(ISYS)
           WK = DELTFL*SFL(I,ISYS) + BFL(I,ISYS)
           IF(IMBWK.EQ.1)  SYSLOADS(3,ISYS) = SYSLOADS(3,ISYS)+DTMB*WK
           IX = IFL(1,I,ISYS)
           IY = IFL(2,I,ISYS)
           IWK = IFL(3,I,ISYS)
           DO IZ=1,NZ
            DM(IX,IY,IZ) = DM(IX,IY,IZ) + ZFRACFL(IZ,IWK)*WK
           ENDDO
  150     CONTINUE
        ENDIF
        IF(NOATM(ISYS).GT.0) THEN
          DELTATM = TIME - NXATMT
          DO 155 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 155
            DO 154 IX=IS,IE
              IF(FSM(IX,IY).LE.0.) GO TO 154
              DM(IX,IY,1) = DM(IX,IY,1) + 
     .          (DELTATM*SATM(IX,IY,ISYS) + BATM(IX,IY,ISYS))*XAZ(IX,IY)
              IF(IMBWK.EQ.1)  SYSLOADS(4,ISYS) = SYSLOADS(4,ISYS) +
     .       DTMB*(DELTATM*SATM(IX,IY,ISYS)+BATM(IX,IY,ISYS))*XAZ(IX,IY)
  154       CONTINUE
  155     CONTINUE
        ENDIF

C  CALCULATE CONCENTRATIONS: UPWIND ADVECTION, KINETICS, LOADS 
        DO 160 IZ=1,NZ
          DO 160 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 160
            DO 159 IX=IS,IE
              IF(FSM(IX,IY).LE.0.) GO TO 159
              CONCS(IX,IY,IZ) = (BVOLB(IX,IY,IZ)*CB(IX,IY,IZ,ISYS)
     .                           + 2.*DT*DM(IX,IY,IZ))/BVOLF(IX,IY,IZ)
 159        CONTINUE
 160    CONTINUE

 
C  BOUNDARY ELEMENTS ARE UNCHANGED
      DO 165 I=1,NOBCALL
        IX=IBCALL(1,I)
        IY=IBCALL(2,I)
        IZ=IBCALL(3,I)
        CONCS(IX,IY,IZ)=CB(IX,IY,IZ,ISYS)
 165  CONTINUE

C        CHECK FOR MASS/FLUX BALANCE COMPUTATIONS
      IF(MASSBAL.EQ.0) GO TO 175
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0 .AND. ITIMESECS.NE.NXPRTMB) GO TO 175
      IF(ITIMESECS.LT.ISMBSECS) GO TO 175
       DO 170 IZ=1,NZ
        DO 170 IY=1,NY
         DO 170 IX=1,NX
          IF(FSM(IX,IY).NE.1.) GO TO 170
          FLXMB(IX,IY,IZ,ISYS) = FLXMB(IX,IY,IZ,ISYS) + DTMB*
     .         (AMAX1(0.,QX(IX,IY,IZ)*CB(IX-1,IY,IZ,ISYS))
     .        + AMIN1(0.,QX(IX,IY,IZ)*CB(IX,IY,IZ,ISYS)))
          FLYMB(IX,IY,IZ,ISYS) = FLYMB(IX,IY,IZ,ISYS) + DTMB*
     .         (AMAX1(0.,QY(IX,IY,IZ)*CB(IX,IY-1,IZ,ISYS))
     .        + AMIN1(0.,QY(IX,IY,IZ)*CB(IX,IY,IZ,ISYS)))
  170  CONTINUE
 
C
C******************************************************************
C
C           SMOLARKIEWICZ METHOD FOR SECOND-ORDER ACCURACY
 
C        CALCULATE ANTIDIFFUSIION FLOWS ( = VELOCITIES * AREAS)
  175 IF(ISMOLAR .EQ. 0) THEN
      DO 180 IZ=1,NZ
        DO 180 IY=2,NY
          DO 180 IX=2,NX
          IF(FSM(IX,IY).EQ.0. .OR. FSM(IX,IY).EQ.-1.) GO TO 180
          IF(FSM(IX,IY).EQ.-2. .AND. ISMOLBCOPT.EQ.0) GO TO 180
            IF(CONCS(IX,IY,IZ).LE.1.0E-09 .OR.
     .            CONCS(IX-1,IY,IZ).LE.1.0E-09) THEN
               ANTIQX(IX,IY,IZ) = 0.0
            ELSE IF(ABS(U(IX,IY,IZ)).LT.(U(IX,IY,IZ)*U(IX,IY,IZ)*
     .                 2.*DT/DXBAR1(IX,IY))) THEN
               ANTIQX(IX,IY,IZ) = 0.0
            ELSE
               ANTIQX(IX,IY,IZ) = ABS(U(IX,IY,IZ))*
     .          (1.-ABS(U(IX,IY,IZ))*2.*DT/DXBAR1(IX,IY))*
     .               (CONCS(IX,IY,IZ)-CONCS(IX-1,IY,IZ))/
     .                    (CONCS(IX,IY,IZ)+CONCS(IX-1,IY,IZ)+R)
            ENDIF


            IF(CONCS(IX,IY,IZ).LE.1.0E-09 .OR.
     .            CONCS(IX,IY-1,IZ).LE.1.0E-09) THEN
               ANTIQY(IX,IY,IZ) = 0.0
            ELSE IF(ABS(V(IX,IY,IZ)).LT.(V(IX,IY,IZ)*V(IX,IY,IZ)*
     .                 2.*DT/DYBAR2(IX,IY))) THEN
               ANTIQY(IX,IY,IZ) = 0.0
            ELSE
               ANTIQY(IX,IY,IZ) = ABS(V(IX,IY,IZ))*
     .          (1.-ABS(V(IX,IY,IZ))*2.*DT/DYBAR2(IX,IY))*
     .               (CONCS(IX,IY,IZ)-CONCS(IX,IY-1,IZ))/
     .                    (CONCS(IX,IY,IZ)+CONCS(IX,IY-1,IZ)+R)
            ENDIF

            IF(FSM(IX,IY).GT.0.) THEN
            IF (IZ.GT.1) THEN
               IF(CONCS(IX,IY,IZ).LE.1.0E-09 .OR.
     .               CONCS(IX,IY,IZ-1).LE.1.0E-09) THEN
                  ANTIQZ(IX,IY,IZ) = 0.0
               ELSE IF(ABS(W(IX,IY,IZ)).LT.(W(IX,IY,IZ)*W(IX,IY,IZ)*
     .                    2.*DT/(D(IX,IY)*DZZ(IZ-1)))) THEN
                  ANTIQZ(IX,IY,IZ) = 0.0
               ELSE
                 ANTIQZ(IX,IY,IZ) = ABS(W(IX,IY,IZ))*
     .                 (1.-ABS(W(IX,IY,IZ))*2.*DT/(D(IX,IY)*DZZ(IZ-1)))*
     .                  (CONCS(IX,IY,IZ-1)-CONCS(IX,IY,IZ))/
     .                       (CONCS(IX,IY,IZ-1)+CONCS(IX,IY,IZ)+R)
               ENDIF
            ENDIF
          ENDIF

 180  CONTINUE
      ENDIF

C-----------------------------------------------------------------------
C        RECURSIVE SMOLARKIEWICZ (ISMOLAR.EQ.1)
      IF(ISMOLAR .EQ. 1) THEN
c$doacross local(iz,iy,ix,aa)
        DO 181 IZ=1,NZ
        DO 181 IY=2,NY-1
        DO 181 IX=2,NX
          IF(CONCS(IX,IY,IZ).LT.1.E-9.OR.CONCS(IX-1,IY,IZ).LT.1.E-9)THEN
            ANTIQX(IX,IY,IZ)=0.0 
          ELSE
            AA=(CONCS(IX,IY,IZ)-CONCS(IX-1,IY,IZ))/
     .             (CONCS(IX-1,IY,IZ)+CONCS(IX,IY,IZ)+1.0E-15)
            ANTIQX(IX,IY,IZ)=ABS(U(IX,IY,IZ))*AA/(1.-ABS(AA)+1.0E-15) 
          END IF
  181   CONTINUE

c$doacross local(iz,iy,ix,aa)
        DO 182 IZ=1,NZ
        DO 182 IY=2,NY
        DO 182 IX=2,NX-1
          IF(CONCS(IX,IY,IZ).LT.1.E-9.OR.CONCS(IX,IY-1,IZ).LT.1.E-9)THEN
            ANTIQY(IX,IY,IZ)=0.0  
          ELSE
            AA=(CONCS(IX,IY,IZ)-CONCS(IX,IY-1,IZ))/
     .             (CONCS(IX,IY-1,IZ)+CONCS(IX,IY,IZ)+1.0E-15)
            ANTIQY(IX,IY,IZ)=ABS(V(IX,IY,IZ))*AA/(1.-ABS(AA)+1.0E-15)  
          END IF
  182   CONTINUE

c$doacross local(iz,iy,ix,aa)
        DO 183 IZ=2,NZ
        DO 183 IY=2,NY-1
        DO 183 IX=2,NX-1
          IF(CONCS(IX,IY,IZ).LT.1.E-9.OR.CONCS(IX,IY,IZ-1).LT.1.E-9)THEN
            ANTIQZ(IX,IY,IZ)=0.0 
          ELSE
            AA=(CONCS(IX,IY,IZ-1)-CONCS(IX,IY,IZ))/
     .             (CONCS(IX,IY,IZ)+CONCS(IX,IY,IZ-1)+1.0E-15)
            ANTIQZ(IX,IY,IZ)=ABS(W(IX,IY,IZ))*AA/(1.-ABS(AA)+1.0E-15)
          END IF
  183   CONTINUE

C
C------------- ADJUST FOR RIVER/WALL FLUXES -------------
        DO 184 IZ=1,NZ
        DO 184 IY=2,NY
        DO 184 IX=2,NX
        IF(FSM(IX,IY).EQ.1..AND.FSM(IX-1,IY).EQ.-1.)ANTIQX(IX,IY,IZ)=0.0
        IF(FSM(IX,IY).EQ.-1..AND.FSM(IX-1,IY).EQ.1.)ANTIQX(IX,IY,IZ)=0.0
        IF(FSM(IX,IY).EQ.1..AND.FSM(IX,IY-1).EQ.-1.)ANTIQY(IX,IY,IZ)=0.0
        IF(FSM(IX,IY).EQ.-1..AND.FSM(IX,IY-1).EQ.1.)ANTIQY(IX,IY,IZ)=0.0
  184   CONTINUE

        DO 185 IZ=1,NZ
        DO 185 IY=1,NY
        DO 185 IX=1,NX
          FLUXIN(IX,IY,IZ)=0.0
          FLUXOUT(IX,IY,IZ)=0.0
  185   CONTINUE

c$doacross local(iz,iy,ix), share(fluxin,fluxout)
        DO 186 IZ=1,NZ
        DO 186 IY=2,NY-1
         DO IX=2,NX
          IF(ANTIQX(IX,IY,IZ).GE.0.)THEN
            FLUXIN( IX,IY,IZ)=FLUXIN(IX,IY,IZ)+
     .             ANTIQX(IX,IY,IZ)/DXBAR1(IX,IY)*CONCS(IX-1,IY,IZ)
          ELSE
            FLUXOUT(IX,IY,IZ)=FLUXOUT(IX,IY,IZ)-
     .             ANTIQX(IX,IY,IZ)/DXBAR1(IX,IY)*CONCS(IX,IY,IZ)
          ENDIF
         ENDDO
         DO IX=1,NX-1
          IF(ANTIQX(IX+1,IY,IZ).GE.0.)THEN
            FLUXOUT(IX,IY,IZ)=FLUXOUT(IX,IY,IZ)+
     .             ANTIQX(IX+1,IY,IZ)/DXBAR1(IX+1,IY)*CONCS(IX,IY,IZ)
          ELSE
            FLUXIN(IX,IY,IZ)=FLUXIN(IX,IY,IZ)-
     .             ANTIQX(IX+1,IY,IZ)/DXBAR1(IX+1,IY)*CONCS(IX+1,IY,IZ)
          ENDIF
         ENDDO
  186   CONTINUE
c$doacross local(iz,iy,ix), share(fluxin,fluxout)
        DO 187 IZ=1,NZ
        DO 187 IX=2,NX-1
         DO  IY=2,NY
          IF(ANTIQY(IX,IY,IZ).GE.0.)THEN
            FLUXIN(IX,IY,IZ)=FLUXIN(IX,IY,IZ)+
     .             ANTIQY(IX,IY,IZ)/DYBAR2(IX,IY)*CONCS(IX,IY-1,IZ)
          ELSE
            FLUXOUT(IX,IY,IZ)=FLUXOUT(IX,IY,IZ)-
     .             ANTIQY(IX,IY,IZ)/DYBAR2(IX,IY)*CONCS(IX,IY,IZ)
          ENDIF
         ENDDO
         DO IY=1,NY-1
          IF(ANTIQY(IX,IY+1,IZ).GE.0.)THEN
            FLUXOUT(IX,IY,IZ)=FLUXOUT(IX,IY,IZ)+
     .             ANTIQY(IX,IY+1,IZ)/DYBAR2(IX,IY+1)*CONCS(IX,IY,IZ)
          ELSE
            FLUXIN(IX,IY,IZ)=FLUXIN(IX,IY,IZ)-
     .             ANTIQY(IX,IY+1,IZ)/DYBAR2(IX,IY+1)*CONCS(IX,IY+1,IZ)
          ENDIF
         ENDDO
  187   CONTINUE
c for layer 1 or layer NLAYER
c$doacross local(iz,iy,ix), share(fluxin,fluxout)
        DO 188 IY=2,NY-1
        DO 188 IX=2,NX-1
          IF(ANTIQZ(IX,IY,2).GE.0.)THEN
            FLUXIN(IX,IY,1)=FLUXIN(IX,IY,1)+
     .             ANTIQZ(IX,IY,2)/(DZZ(1)*D(IX,IY))*CONCS(IX,IY,2)
          ELSE
            FLUXOUT(IX,IY,1)=FLUXOUT(IX,IY,1)-
     .             ANTIQZ(IX,IY,2)/(DZZ(1)*D(IX,IY))*CONCS(IX,IY,1)
          ENDIF
          IF(ANTIQZ(IX,IY,NZ).GE.0.)THEN
           FLUXOUT(IX,IY,NZ)=FLUXOUT(IX,IY,NZ)+ANTIQZ(IX,IY,NZ)
     .                     /(DZZ(NZ-1)*D(IX,IY))*CONCS(IX,IY,NZ)
          ELSE
            FLUXIN(IX,IY,NZ)=FLUXIN(IX,IY,NZ)-ANTIQZ(IX,IY,NZ)
     .                   /(DZZ(NZ-1)*D(IX,IY))*CONCS(IX,IY,NZ-1)
          ENDIF
  188   CONTINUE
c for layer 2 through layer NLAYER-1
c$doacross local(iz,iy,ix), share(fluxin,fluxout)
        DO 189 IZ=2,NZ-1
        DO 189 IY=2,NY-1
        DO 189 IX=2,NX-1
          IF(ANTIQZ(IX,IY,IZ).GE.0.)THEN
            FLUXOUT(IX,IY,IZ)=FLUXOUT(IX,IY,IZ)+ANTIQZ(IX,IY,IZ)
     .                         /(DZZ(IZ-1)*D(IX,IY))*CONCS(IX,IY,IZ)
          ELSE
            FLUXIN(IX,IY,IZ)=FLUXIN(IX,IY,IZ)-ANTIQZ(IX,IY,IZ)
     .                         /(DZZ(IZ-1)*D(IX,IY))*CONCS(IX,IY,IZ-1)
          ENDIF
          IF(ANTIQZ(IX,IY,IZ+1).GE.0.)THEN
            FLUXIN(IX,IY,IZ)=FLUXIN(IX,IY,IZ)+ANTIQZ(IX,IY,IZ+1)
     .                         /(DZZ(IZ)*D(IX,IY))*CONCS(IX,IY,IZ+1)
          ELSE
            FLUXOUT(IX,IY,IZ)=FLUXOUT(IX,IY,IZ)-ANTIQZ(IX,IY,IZ+1)
     .                         /(DZZ(IZ)*D(IX,IY))*CONCS(IX,IY,IZ)
          ENDIF
  189   CONTINUE

c assign caray and concs to be 0.0 temporally 
        DO IZ=1,NZ 
        DO IY=1,NY
        DO IX=1,NX
          IF(FSM(IX,IY).EQ.0.)CB(IX,IY,IZ,ISYS)=0.0
          IF(FSM(IX,IY).EQ.0.)CONCS(IX,IY,IZ)=0.0
        ENDDO
        ENDDO
        ENDDO
        DO 191 IZ=1,NZ 
        DO 191 IY=1,NY
        DO 191 IX=1,NX
          FMAX=AMAX1(CB(IX,IY,IZ,ISYS),
     .     CB(MAX(1,IX-1),IY,IZ,ISYS),CB(MIN(NX,IX+1),IY,IZ,ISYS),
     .     CB(IX,MAX(1,IY-1),IZ,ISYS),CB(IX,MIN(NY,IY+1),IZ,ISYS),
     .     CB(IX,IY,MAX(1,IZ-1),ISYS),CB(IX,IY,MIN(NZ,IZ+1),ISYS),
     .     CONCS(IX,IY,IZ),
     .     CONCS(MAX(1,IX-1),IY,IZ),CONCS(MIN(NX,IX+1),IY,IZ),
     .     CONCS(IX,MAX(1,IY-1),IZ),CONCS(IX,MIN(NY,IY+1),IZ),
     .     CONCS(IX,IY,MAX(1,IZ-1)),CONCS(IX,IY,MIN(NZ,IZ+1)))
          FLUXIN(IX,IY,IZ)=(FMAX-CONCS(IX,IY,IZ))/
     .                 (2.*DT*FLUXIN(IX,IY,IZ)+1.0E-15)
  191   CONTINUE

        DO IZ=1,NZ 
        DO IY=1,NY
        DO IX=1,NX
          IF(FSM(IX,IY).EQ.0.)CB(IX,IY,IZ,ISYS)=1.E+20
          IF(FSM(IX,IY).EQ.0.)CONCS(IX,IY,IZ)=1.E+20
        ENDDO
        ENDDO
        ENDDO
        DO 192 IZ=1,NZ 
        DO 192 IY=1,NY
        DO 192 IX=1,NX
          FMIN=AMIN1(CB(IX,IY,IZ,ISYS),
     .     CB(MAX(1,IX-1),IY,IZ,ISYS),CB(MIN(NX,IX+1),IY,IZ,ISYS),
     .     CB(IX,MAX(1,IY-1),IZ,ISYS),CB(IX,MIN(NY,IY+1),IZ,ISYS),
     .     CB(IX,IY,MAX(1,IZ-1),ISYS),CB(IX,IY,MIN(NZ,IZ+1),ISYS),
     .     CONCS(IX,IY,IZ),
     .     CONCS(MAX(1,IX-1),IY,IZ),CONCS(MIN(NX,IX+1),IY,IZ),
     .     CONCS(IX,MAX(1,IY-1),IZ),CONCS(IX,MIN(NY,IY+1),IZ),
     .     CONCS(IX,IY,MAX(1,IZ-1)),CONCS(IX,IY,MIN(NZ,IZ+1)))
          FLUXOUT(IX,IY,IZ)=(CONCS(IX,IY,IZ)-FMIN)/
     .                  (2.*DT*FLUXOUT(IX,IY,IZ)+1.0E-15)
  192   CONTINUE
 
        DO 193 IZ=1,NZ 
        DO 193 IY=1,NY
        DO 193 IX=1,NX
          IF(FSM(IX,IY).EQ.0.)CB(IX,IY,IZ,ISYS)=0.0
          IF(FSM(IX,IY).EQ.0.)CONCS(IX,IY,IZ)=0.0
          IF(FSM(IX,IY).EQ.-1.)FLUXIN(IX,IY,IZ)=0.
          IF(FSM(IX,IY).EQ.-1.)FLUXOUT(IX,IY,IZ)=0.
  193   CONTINUE
        DO 194 IZ=1,NZ
        DO 194 IY=2,NY
        DO 194 IX=2,NX
          IF(ANTIQX(IX,IY,IZ).GE.0.)THEN
            ANTIQX(IX,IY,IZ)=ANTIQX(IX,IY,IZ)*
     .            (AMIN1(1.,FLUXIN(IX,IY,IZ),FLUXOUT(IX-1,IY,IZ)))
          ELSE
            ANTIQX(IX,IY,IZ)=ANTIQX(IX,IY,IZ)*
     .            (AMIN1(1.,FLUXIN(IX-1,IY,IZ),FLUXOUT(IX,IY,IZ)))
          ENDIF
          IF(ANTIQY(IX,IY,IZ).GE.0.)THEN
            ANTIQY(IX,IY,IZ)=ANTIQY(IX,IY,IZ)*
     .            (AMIN1(1.,FLUXIN(IX,IY,IZ),FLUXOUT(IX,IY-1,IZ)))
          ELSE
            ANTIQY(IX,IY,IZ)=ANTIQY(IX,IY,IZ)*
     .            (AMIN1(1.,FLUXIN(IX,IY-1,IZ),FLUXOUT(IX,IY,IZ)))
          ENDIF
          IF(IZ.GE.2)THEN
            IF(ANTIQZ(IX,IY,IZ).GE.0.)THEN
              ANTIQZ(IX,IY,IZ)=ANTIQZ(IX,IY,IZ)*
     .            (AMIN1(1.,FLUXIN(IX,IY,IZ-1),FLUXOUT(IX,IY,IZ)))
            ELSE
              ANTIQZ(IX,IY,IZ)=ANTIQZ(IX,IY,IZ)*
     .            (AMIN1(1.,FLUXIN(IX,IY,IZ),FLUXOUT(IX,IY,IZ-1)))
            ENDIF
          ENDIF
  194   CONTINUE

      ENDIF

C  END OF EVALUATION OF TERMS REQUIRED FOR SMOLARKIEWICZ CORRECTION

C-----------------------------------------------------------------------
C  NOW USE TERMS TO CALCULATE ANTIDIFFUSION FLUXES

      DO 196 IZ=1,NZ
        DO 196 IY=2,NY
          DO 196 IX=2,NX
              QXFLX(IX,IY,IZ) = 0.5*(ANTIQX(IX,IY,IZ)*(CONCS(IX,IY,IZ)+
     .           CONCS(IX-1,IY,IZ))-ABS(ANTIQX(IX,IY,IZ))*
     .           (CONCS(IX,IY,IZ)-CONCS(IX-1,IY,IZ)))
     .           *DYBAR1(IX,IY)*0.5*(D(IX,IY)+D(IX-1,IY))*DZ(IZ)
              QYFLX(IX,IY,IZ) = 0.5*(ANTIQY(IX,IY,IZ)*(CONCS(IX,IY,IZ)+
     .          CONCS(IX,IY-1,IZ))-ABS(ANTIQY(IX,IY,IZ))*
     .          (CONCS(IX,IY,IZ)-CONCS(IX,IY-1,IZ)))      
     .           *DXBAR2(IX,IY)*0.5*(D(IX,IY)+D(IX,IY-1))*DZ(IZ)
            IF(NZ.GT.1.AND.IZ.GT.1.AND.FSM(IX,IY).GT.0.) 
     .        QZFLX(IX,IY,IZ) = 0.5*
     .        (ANTIQZ(IX,IY,IZ)*(CONCS(IX,IY,IZ-1)+CONCS(IX,IY,IZ))-
     .          ABS(ANTIQZ(IX,IY,IZ))*(CONCS(IX,IY,IZ-1)-
     .          CONCS(IX,IY,IZ)))*XAZ(IX,IY)

            IF(FSM(IX,IY).EQ.0. .OR. FSM(IX,IY).EQ.-1.) THEN
             QXFLX(IX,IY,IZ) = 0.0
             QYFLX(IX,IY,IZ) = 0.0
             QZFLX(IX,IY,IZ) = 0.0
            ENDIF
 196  CONTINUE
 
C  SET HORIZONTAL ANTIDIFFUSION FLUXES IN BOUNDARY ELEMENTS TO ZERO
        IF(ISMOLBCOPT.EQ.0) THEN
          DO 198 I=1,NOBCALL
            IX=IBCALL(1,I)
            IY=IBCALL(2,I)
            IZ=IBCALL(3,I)
            QXFLX(IX,IY,IZ)=0.0
            QXFLX(IX+1,IY,IZ)=0.0
            QYFLX(IX,IY,IZ)=0.0
            QYFLX(IX,IY+1,IZ)=0.0
 198      CONTINUE
        ENDIF


C        TAKE CORRECTION STEP, PLACING RESULTS
C          BACK INTO CONCS NOW AT '(N+1)'
        IF (NZ.GT.1) THEN
          DO 200 IZ=1,NZ-1
            DO 200 IY=1,NY
              IS = IXS(IY)
              IE = IXE(IY)
              IF(IS.EQ.0) GO TO 200
              DO 199 IX=IS,IE
                IF(FSM(IX,IY).LE.0.) GO TO 199
                CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)
     .             -2.*DT*(QXFLX(IX+1,IY,IZ)-QXFLX(IX,IY,IZ)
     .                    +QYFLX(IX,IY+1,IZ)-QYFLX(IX,IY,IZ)
     .                    +QZFLX(IX,IY,IZ)-QZFLX(IX,IY,IZ+1))/
     .                                             BVOLF(IX,IY,IZ)
 199          CONTINUE
 200      CONTINUE
          IZ=NZ
          DO 205 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 205
            DO 204 IX=IS,IE
              IF(FSM(IX,IY).LE.0.) GO TO 204
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)
     .           -2.*DT*(QXFLX(IX+1,IY,IZ)-QXFLX(IX,IY,IZ)
     .                  +QYFLX(IX,IY+1,IZ)-QYFLX(IX,IY,IZ)
     .                  +QZFLX(IX,IY,IZ))/BVOLF(IX,IY,IZ)
 204        CONTINUE
 205      CONTINUE
        ELSE
          IZ=1
          DO 210 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 210
            DO 209 IX=IS,IE
              IF(FSM(IX,IY).LE.0.) GO TO 209
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)
     .           -2.*DT*(QXFLX(IX+1,IY,IZ)-QXFLX(IX,IY,IZ)
     .                  +QYFLX(IX,IY+1,IZ)-QYFLX(IX,IY,IZ))/
     .                                             BVOLF(IX,IY,IZ)
 209        CONTINUE
 210      CONTINUE
        ENDIF

C  END OF SMOLARKIEWICZ CORRECTION CODE

C        CHECK FOR MASS/FLUX BALANCE COMPUTATIONS
      IF(MASSBAL.EQ.0) GO TO 213
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0 .AND. ITIMESECS.NE.NXPRTMB) GO TO 213
      IF(ITIMESECS.LT.ISMBSECS) GO TO 213
      DO IZ=1,NZ
       DO IY=1,NY
        DO IX=1,NX
          FLXMB(IX,IY,IZ,ISYS)=FLXMB(IX,IY,IZ,ISYS)+DTMB*QXFLX(IX,IY,IZ)
          FLYMB(IX,IY,IZ,ISYS)=FLYMB(IX,IY,IZ,ISYS)+DTMB*QYFLX(IX,IY,IZ)
        ENDDO
       ENDDO
      ENDDO 

C**********************************************************************

C  CALCULATE HORIZONTAL DIFFUSION FLUXES
 
C        CHECK FOR SIGMA-LEVEL CORRECTION OPTION
  213   IF(SLCOPT.EQ.1 .AND.
     .      MOD((ITIMESECS-IFIX(86400.*TZERO)),IDTSLCSECS).LE.0.) THEN
         CALL RCAMPROF(CB(1,1,1,ISYS),CSLAVE(1,1,1,ISYS))
        ENDIF
 
C  CALCULATE HORIZONTAL DIFFUSION FLUXES
C  - USING CONCENTRATIONS AT TIME LEVEL "N-1" -
c$doacross local(iz,iy,ix,is,ie)
        DO 220 IZ=1,NZ
        DO 220 IY=2,NY
        DO 220 IX=2,NX
              IF(FSM(IX,IY).GE.1. .OR. FSM(IX,IY).LE.-2.)THEN
                IF(FSM(IX-1,IY).GE.1. .OR. FSM(IX-1,IY).LE.-2.)THEN
                  RXFLX(IX,IY,IZ)=RX(IX,IY,IZ)*
     .             ((CB(IX-1,IY,IZ,ISYS)-CSLAVE(IX-1,IY,IZ,ISYS)) -
     .             (CB(IX,IY,IZ,ISYS)-CSLAVE(IX,IY,IZ,ISYS)))
                ENDIF
                IF(FSM(IX,IY-1).GE.1. OR. FSM(IX,IY-1).LE.-2.)THEN
                  RYFLX(IX,IY,IZ)=RY(IX,IY,IZ)*
     .             ((CB(IX,IY-1,IZ,ISYS)-CSLAVE(IX,IY-1,IZ,ISYS)) -
     .             (CB(IX,IY,IZ,ISYS)-CSLAVE(IX,IY,IZ,ISYS)))
                ENDIF
              ENDIF
 220    CONTINUE
 
C  CALCULATE CONCENTRATIONS: HORIZONTAL DIFFUSION
        DO 230 IZ=1,NZ
          DO 230 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0) GO TO 230
            DO 229 IX=IS,IE
              IF(FSM(IX,IY).LE.0.) GO TO 229
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)+2.*DT*
     .           (RXFLX(IX,IY,IZ)-RXFLX(IX+1,IY,IZ)+RYFLX(IX,IY,IZ)
     .           -RYFLX(IX,IY+1,IZ))/BVOLF(IX,IY,IZ)
 229        CONTINUE
 230    CONTINUE
 
C  BOUNDARY ELEMENTS ARE UNCHANGED
        DO 235 I=1,NOBCALL
          IX=IBCALL(1,I)
          IY=IBCALL(2,I)
          IZ=IBCALL(3,I)
          CONCS(IX,IY,IZ)=CB(IX,IY,IZ,ISYS)  !!!!!!!!!!!!!!!!!!!!!!!!!!
  235   CONTINUE

C        CHECK FOR MASS/FLUX BALANCE COMPUTATIONS
      IF(MASSBAL.EQ.0) GO TO 237
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0 .AND. NXPRTMB.NE.ITIMESECS) GO TO 237
      IF(ITIMESECS.LT.ISMBSECS) GO TO 237
      DO IZ=1,NZ
       DO IY=1,NY
        DO IX=1,NX
         FLXMB(IX,IY,IZ,ISYS)=FLXMB(IX,IY,IZ,ISYS)+DTMB*RXFLX(IX,IY,IZ)
         FLYMB(IX,IY,IZ,ISYS)=FLYMB(IX,IY,IZ,ISYS)+DTMB*RYFLX(IX,IY,IZ)
        ENDDO
       ENDDO
      ENDDO
 
C**********************************************************************
C  COMPLETE REST OF INTEGRATION STEP
 
C        IMPLICIT VERTICAL DIFFUSION STEP (IF NZ > 1)
 237    IF(NZ.EQ.1) GO TO 280
        DO 240 IZ=1,NZ
         DO 240 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0) GO TO 240
          DO 239 IX=IS,IE
            MF(IX,IY,IZ) = CONCS(IX,IY,IZ)*BVOLF(IX,IY,IZ)
 239      CONTINUE
 240    CONTINUE

C        SET UP DIAGONAL AND OFF-DIAGONAL ELEMENTS
        DO 250 IY=1,NY
         IS = IXS(IY)
         IE = IXE(IY)
         IF(IS.EQ.0) GO TO 250
         DO 249 IX=IS,IE
          IF(FSM(IX,IY).LE.0) GO TO 249
          AL(IX,IY,1) = 0.
          AD(IX,IY,1) = BVOLF(IX,IY,1) + 2.*DT*RZ(IX,IY,2)
          AU(IX,IY,1) = -2.*DT*RZ(IX,IY,2)
           DO 245 IZ=2,NZ-1
            AL(IX,IY,IZ) = -2.*DT*RZ(IX,IY,IZ)
            AD(IX,IY,IZ) = BVOLF(IX,IY,IZ) +
     .                        2.*DT*(RZ(IX,IY,IZ) + RZ(IX,IY,IZ+1))
            AU(IX,IY,IZ) = -2.*DT*RZ(IX,IY,IZ+1)
  245      CONTINUE
          AL(IX,IY,NZ) = -2.*DT*RZ(IX,IY,NZ)
          AD(IX,IY,NZ) = BVOLF(IX,IY,NZ) +
     .                        2.*DT*RZ(IX,IY,NZ)
          AU(IX,IY,NZ) = 0.0
  249    CONTINUE
  250   CONTINUE

C        FORWARD SWEEP OF TRIDIAGONAL SCHEME
        DO 260 IY=1,NY
         IS = IXS(IY)
         IE = IXE(IY)
         IF(IS.EQ.0) GO TO 260
         DO 259 IX=IS,IE
          IF(FSM(IX,IY).LE.0) GO TO 259
          BETA(IX,IY,1) = AD(IX,IY,1)
          GAMMA(IX,IY,1) = MF(IX,IY,1)/BETA(IX,IY,1)
          DO 255 IZ=2,NZ
           BETA(IX,IY,IZ) = AD(IX,IY,IZ)
     .                 - AL(IX,IY,IZ)*AU(IX,IY,IZ-1)/BETA(IX,IY,IZ-1)
           GAMMA(IX,IY,IZ) = (MF(IX,IY,IZ)
     .                 - AL(IX,IY,IZ)*GAMMA(IX,IY,IZ-1))/BETA(IX,IY,IZ)
  255     CONTINUE
  259    CONTINUE
  260   CONTINUE

C        BACKWARD SWEEP FOR SOLUTION
        DO 270 IY=1,NY
         IS = IXS(IY)
         IE = IXE(IY)
         IF(IS.EQ.0) GO TO 270
         DO 269 IX=IS,IE
          IF(FSM(IX,IY).LE.0) GO TO 269
          CF(IX,IY,NZ,ISYS) = GAMMA(IX,IY,NZ)
          DO 265 ILYR=2,NZ
           IZ = NZ+1-ILYR
           CF(IX,IY,IZ,ISYS) = GAMMA(IX,IY,IZ)
     .               - AU(IX,IY,IZ)*CF(IX,IY,IZ+1,ISYS)/BETA(IX,IY,IZ)
  265     CONTINUE
  269    CONTINUE
  270   CONTINUE

C  BOUNDARY ELEMENTS ARE UNCHANGED
      DO 275 I=1,NOBCALL
        IX=IBCALL(1,I)
        IY=IBCALL(2,I)
        IZ=IBCALL(3,I)
        CF(IX,IY,IZ,ISYS)=CB(IX,IY,IZ,ISYS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  275 CONTINUE

  280   CONTINUE


C  FILTER SOLUTIONS
        DO 300 IZ=1,NZ
         DO 300 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0) GO TO 300
          DO 299 IX=IS,IE
           IF(FSM(IX,IY).LE.0.) GO TO 299
           C(IX,IY,IZ,ISYS) = C(IX,IY,IZ,ISYS) + 0.025 *
     .        (CF(IX,IY,IZ,ISYS)+CB(IX,IY,IZ,ISYS)-2.*C(IX,IY,IZ,ISYS))
  299     CONTINUE
  300   CONTINUE

  320 CONTINUE

C        FLUX BALANCES FOR MISC KINETIC VARIABLES
      IF(MASSBAL.EQ.0 .OR. NOKINSYS.EQ.0) GO TO 330
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0) THEN 
       IF(ITIMESECS.NE.NXPRTMB) GO TO 330
       DO 325 ISYS=1,NOKINSYS
        DO 325 IZ=1,NZ
         DO 325 IY=1,NY
          DO 325 IX=1,NX
           IF(FSM(IX,IY).NE.1.) GO TO 325
            FLXMB(IX,IY,IZ,NOSYS+ISYS) = 
     .         (AMAX1(0.,QX(IX,IY,IZ)*CKINARRAY(IX-1,IY,IZ,ISYS))
     .        + AMIN1(0.,QX(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS))
     .        + RX(IX,IY,IZ)
     .           *(CKINARRAY(IX-1,IY,IZ,ISYS)-CKINARRAY(IX,IY,IZ,ISYS)))
            FLYMB(IX,IY,IZ,NOSYS+ISYS) = 
     .         (AMAX1(0.,QY(IX,IY,IZ)*CKINARRAY(IX,IY-1,IZ,ISYS))
     .        + AMIN1(0.,QY(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS))
     .        + RY(IX,IY,IZ)
     .           *(CKINARRAY(IX,IY-1,IZ,ISYS)-CKINARRAY(IX,IY,IZ,ISYS)))
  325  CONTINUE
      ELSE
       IF(ITIMESECS.LT.ISMBSECS) GO TO 330
       DO 327 ISYS=1,NOKINSYS
        DO 327 IZ=1,NZ
         DO 327 IY=1,NY
          DO 327 IX=1,NX
           IF(FSM(IX,IY).NE.1.) GO TO 327
            FLXMB(IX,IY,IZ,NOSYS+ISYS) = FLXMB(IX,IY,IZ,NOSYS+ISYS) 
     .        +DT*(AMAX1(0.,QX(IX,IY,IZ)*CKINARRAY(IX-1,IY,IZ,ISYS))
     .             + AMIN1(0.,QX(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS))
     .             + RX(IX,IY,IZ)
     .           *(CKINARRAY(IX-1,IY,IZ,ISYS)-CKINARRAY(IX,IY,IZ,ISYS)))
            FLYMB(IX,IY,IZ,NOSYS+ISYS) = FLYMB(IX,IY,IZ,NOSYS+ISYS) 
     .        +DT*(AMAX1(0.,QY(IX,IY,IZ)*CKINARRAY(IX,IY-1,IZ,ISYS))
     .             + AMIN1(0.,QY(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS))
     .             + RY(IX,IY,IZ)
     .           *(CKINARRAY(IX,IY-1,IZ,ISYS)-CKINARRAY(IX,IY,IZ,ISYS)))
  327  CONTINUE
      ENDIF 


C  PREPARE FOR NEXT TIME STEP
  330 CONTINUE
      ITIMESECS = ITIMESECS + IDTSECS
      TIME = ITIMESECS/86400.
      DO 340 ISYS=1,NOSYS
       DO 340 IZ=1,NZ
        DO 340 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0) GO TO 340
          DO 339 IX=IS,IE
            IF(FSM(IX,IY).LE.0.) GO TO 339
            CB(IX,IY,IZ,ISYS) = C(IX,IY,IZ,ISYS)
            C(IX,IY,IZ,ISYS) = CF(IX,IY,IZ,ISYS)
  339     CONTINUE
  340 CONTINUE
       DO 345 IZ=1,NZ
         DO 345 IY=1,NY
           DO 345 IX=1,NX
          BVOLB(IX,IY,IZ) = BVOL(IX,IY,IZ)
          BVOL(IX,IY,IZ) = BVOLF(IX,IY,IZ)
  345 CONTINUE

      DO 370 ISYS=1,NOSYS
        DO 370 IZ=1,NZ
          DO 370 IY=1,NY
            DO 370 IX=1,NX
            IF(FSM(IX,IY).LE.0.) GO TO 370
            CARAY(IX,IY,IZ,ISYS) = MAX(CB(IX,IY,IZ,ISYS),CMIN(ISYS))
            CB(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)
  370 CONTINUE

C        CHECK TO SEE IF IT IS TIME TO PRINT RESULTS
C        FIRST DETAILED DUMPS
      IF(ITIMESECS.LT.NXPRTD) GO TO 385
C        DETAILED PRINT TIME
C        SET IDISK FOR DETAILED DUMPS
      IDISK=2
      NXPRTD = NXPRTD + IPRNTDSECS
C        SECOND GLOBAL DUMPS
 385  IF(ITIMESECS.LT.NXPRTG) GO TO 401
C        GLOBAL PRINT TIME
C        CALL PRINT ROUTINE
      CALL RCA09
      IDISK=IDISK+1
      NXPRTG = NXPRTG + IPRNTGSECS

C  NOW PERFORM STABILITY CHECK
      DO 400 ISYS=1,NOSYS 
        IF(SYSBY(ISYS).EQ.1) GO TO 400
        CKMAX = CMAX(ISYS)
        DO 398 IZ=1,NZ
         DO 398 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0) GO TO 398
          DO 397 IX=IS,IE
           IF(FSM(IX,IY).LE.0) GO TO 397
           IF(CF(IX,IY,IZ,ISYS).GE.CKMAX) GO TO 430
  397     CONTINUE
  398   CONTINUE
  400 CONTINUE

C  MASS/FLUX BALANCE DUMPS
  401 IF(MASSBAL.EQ.0) GO TO 420
      IF(ITIMESECS.LT.NXPRTMB) GO TO 420
       IF(IMBDOPT.EQ.0) THEN
        WRITE(17)  TIME,SYSMASS,SYSLOADS,FLXMB,FLYMB
       ELSE
        TIMEAVE=FLOAT(IPRNTMBSECS)/86400.
        DO 405 ISYS=1,NOSYS+NOKINSYS
         DO 405 IZ=1,NZ
          DO 405 IY=1,NY
           DO 405 IX=1,NX
            FLXMB(IX,IY,IZ,ISYS)=FLXMB(IX,IY,IZ,ISYS)/TIMEAVE
            FLYMB(IX,IY,IZ,ISYS)=FLYMB(IX,IY,IZ,ISYS)/TIMEAVE
  405   CONTINUE
        DO ISYS=1,NOSYS
         DO I=1,4
          SYSLOADS(I,ISYS)=SYSLOADS(I,ISYS)/TIMEAVE
         ENDDO
        ENDDO
        WRITE(17)  TIME,SYSMASS,SYSLOADS,FLXMB,FLYMB
       ENDIF
C        RESET ARRAYS AND COUNTER
      DO 410 ISYS=1,NOSYS
       DO 410 I=1,4
        SYSLOADS(I,ISYS)=0.
  410 CONTINUE
      DO 415 ISYS=1,NOSYS+NOKINSYS
       DO 415 IZ=1,NZ
        DO 415 IY=1,NY
         DO 415 IX=1,NX
          FLXMB(IX,IY,IZ,ISYS)=0.0
          FLYMB(IX,IY,IZ,ISYS)=0.0
  415 CONTINUE
      NXPRTMB = NXPRTMB + IPRNTMBSECS
      IF(NXPRTMB.GT.IEMBSECS) MASSBAL=0

C        CHECK TO SEE IF NECESSARY TO UPDATE TIME-VARIABLE FUNCTIONS
C        POINT SOURCE LOADS 
  420 IF(IPSOPT.GT.1.AND.TIME.GE.NXPST)   
     .   CALL RCA10(SPS,BPS,MXWK,NOPS,NXPST,33,IPSPWLOPT,SCALPS)
C        NONPOINT SOURCE LOADS 
      IF(INPSOPT.GT.1.AND.TIME.GE.NXNPST)   
     .   CALL RCA10(SNPS,BNPS,MXWK,NONPS,NXNPST,34,INPSPWLOPT,SCALNPS)
C        FALL-LINE LOADS 
      IF(IFLOPT.GT.1.AND.TIME.GE.NXFLT)   
     .   CALL RCA10(SFL,BFL,MXWK,NOFL,NXFLT,35,IFLPWLOPT,SCALFL)
C        ATMOSPHERIC LOADS 
      IF(IATMOPT.GT.1.AND.TIME.GE.NXATMT)   
     .   CALL RCA10(SATM,BATM,NX*NY,NOATM,NXATMT,36,IATMPWLOPT,
     .              SCALATM)
C        NEXT BOUNDARY CONDITIONS 
      IFLAG_TBC=1
      IF((IBCOPT.EQ.2.OR.IBCOPT.EQ.4).AND.TIME.GE.NXBCT) THEN
         IFLAG_TBC=0
         CALL RCA11 
      ENDIF

      IF(IBCPWLOPT.NE.1) THEN
        DO 425 ISYS=1,NOSYS
        DO 425 I=1,NOBCALL
          IX=IBCALL(1,I)
          IY=IBCALL(2,I)
          IZ=IBCALL(3,I)
          CB(IX,IY,IZ,ISYS)=CARAY(IX,IY,IZ,ISYS) !! stepwise
 425    CONTINUE
      ENDIF

C        NEXT TRANSPORT FIELDS
         IFLAG_TBRK=1
C     IF(ITIMESECS.GE.HYDBRK(ITIMHYD)) THEN
      IF(ITIMESECS.EQ.NXHYDTSECS) THEN
         IFLAG_TBRK=0
         CALL RCA03A
         CALL RCAEXP1
C        ITIMHYD = ITIMHYD+1

      ENDIF
   
       IF(IFLAG_TBRK.EQ.0)THEN

       DO 428 IZ=1,NZ
         DO 428 IY=1,NY
           IS = IXS(IY)
           IE = IXE(IY)
           IF(IS.EQ.0) GO TO 428
           DO 427 IX=IS,IE
             BVOLB(IX,IY,IZ) = BVOL(IX,IY,IZ)
  427      CONTINUE
  428  CONTINUE

       ENDIF

C        CHECK FOR END OF SIMULATION...
C           IF NOT GO BACK AND TAKE NEXT INTEGRATION STEP
C
      IF(TIME.LE.TEND) GO TO 40 
      ITRAK = ITRAK+1 
      IF(ITRAK.LE.NSTEP) THEN
         IDTSECS=ISTEP(ITRAK)
         DT=ISTEP(ITRAK)
         TEND = TBRK(ITRAK) 
         IF(DT.EQ.0.0) THEN
            CALL RCAMESS(3,DT)
            IDISK=1
            INITB=2
            CALL TUNER
            CALL EXIT
         ENDIF
      ELSE
         TEND = 0.0
      ENDIF 

      IF(TEND.GT.0.) GO TO 40 
 
C        FINISHED...
      WRITE(OUT,8000)
 8000 FORMAT(///30X,'RCALFS FINISHED INTEGRATION'/
     .          30X,'USER DUMPS TO FOLLOW'//)
      IF(MASSBAL.EQ.1) THEN
        DEALLOCATE(FLXMB)
        DEALLOCATE(FLYMB)
      ENDIF
      RETURN

C
C         STABILITY CRITERIA VIOLATED...ABEND 
C
  430 CALL RCAMESS(1,CKMAX) 
      WRITE(OUT,9000)  BVOLB(IX,IY,IZ),DIAG(IX,IY,IZ)
 9000 FORMAT(10X,'VOLUME =',E15.6,'  ,DIAGONAL =',E13.4)
      IREC = IREC+1 
      IDISK = 1 
      INITB = 2
      CALL TUNER
      RETURN
      END 
