      MODULE ALLOCRCAEXPS
       REAL, ALLOCATABLE :: FLXMB(:,:,:,:)
       REAL, ALLOCATABLE :: FLYMB(:,:,:,:)
      END MODULE ALLOCRCAEXPS
      SUBROUTINE RCAEXPS
      USE ALLOCRCAEXPS
C
C        RCAEXPS IS THE MASTER ROUTINE USED TO INTEGRATE THE DIFFERENTIAL
C                EQUATIONS OF THE WATER QUALITY MODEL.
C                TIME DERIVATIVE: EXPLICIT - 1ST ORDER EULER
C                SPACE DERIVATIVE: UPWIND WITH SMOLARKIEWICZ CORRECTOR
C
C************************************************************************
C
      SAVE
      INCLUDE 'RCACM'            
      INTEGER  NEGSEG(99,3)
      REAL   CONCS(NX,NY,0:NZ),CSLAVE(NX,NY,NZ,NOSYS)
     .      ,DM(NX,NY,NZ),BVOLN(NX,NY,NZ)
     .      ,DXBAR1(NX,NY),DYBAR1(NX,NY)
     .      ,DXBAR2(NX,NY),DYBAR2(NX,NY)
     .      ,D(NX,NY),U(NX,NY,NZ),V(NX,NY,NZ)
     .      ,W(NX,NY,NZ)
     .      ,FLUXIN(NX,NY,NZ),FLUXOUT(NX,NY,NZ)
     .      ,ANTIQX(NX,NY,NZ),QXFLX(NX,NY,NZ)
     .      ,ANTIQY(NX,NY,NZ),QYFLX(NX,NY,NZ)
     .      ,ANTIQZ(NX,NY,NZ),QZFLX(NX,NY,NZ)
     .      ,RXFLX(NX,NY,NZ),RYFLX(NX,NY,NZ)
     .      ,AL(NX,NY,NZ),AD(NX,NY,NZ),AU(NX,NY,NZ)
     .      ,BETA(NX,NY,NZ),GAMMA(NX,NY,NZ)
     .      ,MF(NX,NY,NZ)
     .      ,SYSMASS(NOSYS+NOKINSYS),SYSLOADS(4,NOSYS)
     .      ,ADDMASS(NX,NY,NZ,NOSYS)
       COMMON SKIPHOURS
C        INITIALIZATION
      INITB=0
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
      IF(TZERO.NE.0.)   THEN
         CALL RCA08
         ITIMESECS = 86400.*TZERO
         ITIMEWQSECS = 86400.*TZERO
         TIME = TZERO
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
       ALLOCATE(FLXMB(NX,NY,NZ,31),STAT=ISTAT)   !ldh
       ALLOCATE(FLYMB(NX,NY,NZ,31),STAT=ISTAT)   !ldh
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

c$doacross local(isys,iz,iy,ix) , share(cslave)
      DO 25 ISYS=1,NOSYS
       DO 25 IZ=1,NZ
        DO 25 IY=1,NY
         DO 25 IX=1,NX
            CSLAVE(IX,IY,IZ,ISYS) = 0.0
   25 CONTINUE

C        SET UP -AVECT- ARRAYS FOR TIME = 0.0 (TZERO)
      CALL RCAEXP1

C        PRINT INITIAL CONDITIONS
      CALL RCA09
      IDISK = 3 

C        BEGINNING OF FULL INTEGRATION TIMESTEP LOOP
   30 CONTINUE

C        GET TOTAL MASS PER SYSTEM IF MASS/FLUX BALANCES REQUESTED
      IF(MASSBAL.EQ.1 .AND. ITIMESECS.EQ.NXPRTMB) THEN
       DO 35 ISYS=1,NOSYS
        SYSMASS(ISYS)=0.0
        DO 35 IZ=1,NZ
         DO 35 IY=1,NY
          DO 35 IX=1,NX
           IF(FSM(IX,IY).EQ.1.)  SYSMASS(ISYS) = SYSMASS(ISYS) +
     .                           BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
  35   CONTINUE
       DO 40 ISYS=1,NOKINSYS
        SYSMASS(NOSYS+ISYS)=0.0
        DO 40 IZ=1,NZ
         DO 40 IY=1,NY
          DO 40 IX=1,NX
           IF(FSM(IX,IY).EQ.1.)
     .                SYSMASS(NOSYS+ISYS) = SYSMASS(NOSYS+ISYS) +
     .                BVOL(IX,IY,IZ)*CKINARRAY(IX,IY,IZ,ISYS)
  40   CONTINUE
      ENDIF

C        IF PIECEWISE LINEAR BOUNDARY CONDITIONS THEN
C           GET CONCENTRATIONS AT CURRENT TIME LEVEL
      IF(IBCPWLOPT.EQ.1)   THEN
c$doacross local(isys,deltbc,iz,iy,ix,i)
        DO 50 ISYS=1,NOSYS
          IF(NOBC(ISYS).EQ.0)  GO TO 50
          DELTBC = TIME - NXBCT
          DO 45 I=1,NOBC(ISYS)
            IX=IBC(1,I,ISYS)
            IY=IBC(2,I,ISYS)
            IZ=IBC(3,I,ISYS)
            CARAY(IX,IY,IZ,ISYS) = DELTBC*SBC(I,ISYS) + BBC(I,ISYS)
   45     CONTINUE
   50  CONTINUE
      ENDIF
 
C        EVALUATE KINETIC DERIVATIVES
      IF (DTWQ.EQ.0.0.OR.ITIMESECS.GE.ITIMEWQSECS) THEN
        CALL TUNER
        IDISK = 0
        ITIMEWQSECS=ITIMEWQSECS+IDTWQSECS
        TIMEWQ=ITIMEWQSECS/86400.
      ENDIF
      INITB = 1
 
C        COMPUTE DEPTHS AT TIME LEVEL N
      DELT = TIME - (NXHYDTSECS-IHYDDTSECS)/86400.
c$doacross local(iy,ix) , share(d)
      DO 60 IY=1,NY
       DO 60 IX=1,NX
        D(IX,IY) = H(IX,IY)+ETA(IX,IY)+DELT*DETA(IX,IY)*86400.
  60  CONTINUE

C        CALCULATE CONSTANTS USED FOR SMOLARKIEWICZ METHOD
c$doacross local(iy,ix) , share(dxbar1,dxbar2,dybar1,dybar2)
      DO 65 IY=2,NY
       DO 65 IX=2,NX
        DXBAR1(IX,IY) = 0.5*(DX(IX,IY)+DX(IX-1,IY))
        DXBAR2(IX,IY) = 0.5*(DX(IX,IY)+DX(IX,IY-1))
        DYBAR1(IX,IY) = 0.5*(DY(IX,IY)+DY(IX-1,IY))
        DYBAR2(IX,IY) = 0.5*(DY(IX,IY)+DY(IX,IY-1))
   65 CONTINUE
c$doacross local(iz,iy,ix) , share(u,v,w)
      DO 70 IZ=1,NZ
       DO 70 IY=2,NY
        DO 70 IX=2,NX
          U(IX,IY,IZ)=QX(IX,IY,IZ)/
     .       (DZ(IZ)*0.5*(D(IX,IY)+D(IX-1,IY))*DYBAR1(IX,IY))
          V(IX,IY,IZ)=QY(IX,IY,IZ)/
     .       (DZ(IZ)*0.5*(D(IX,IY)+D(IX,IY-1))*DXBAR2(IX,IY))
          W(IX,IY,IZ)=QZ(IX,IY,IZ)/XAZ(IX,IY)
 70   CONTINUE

C             BEGIN INTEGRATION STEP

C        COMPUTE VOLUMES AT TIME LEVEL N+1
c$doacross local(iz,iy,ix,is,ie)
      DO 80 IZ=1,NZ
       DO 80 IY=1,NY
        IS = IXS(IY)
        IE = IXE(IY)
        IF(IS.EQ.0) GO TO 80
        DO 79 IX=IS,IE
          BVOLN(IX,IY,IZ) = BVOL(IX,IY,IZ) + DT*VDER(IX,IY,IZ)
  79    CONTINUE
  80  CONTINUE

C        MAIN SYSTEM LOOP

      DO 400 ISYS=1,NOSYS 
C        CHECK FOR SYSTEM BYPASS
        IF(SYSBY(ISYS).EQ.1)   GO TO 400
        INEGS=0
        TOTMAS=0.0
 
C        THE 120 THROUGH 132 DO LOOPS TAKE AN EULER INTEGRATION STEP
C        STEP FOR HORIZONTAL AND VERTICAL ADVECTIVE TRANSPORT
C           (-DIAG- AND -AVECT- CONTAIN JUST  Q  TERMS)
        IF (NZ.GT.1) THEN
C      SURFACE LAYER
c$doacross local(iy,ix,is,ie)
        DO 120 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0)   GO TO 120
          DO 119 IX=IS,IE
            IF(FSM(IX,IY).LE.0.)  GO TO 119
            DM(IX,IY,1) = - DIAG(IX,IY,1)*CARAY(IX,IY,1,ISYS)
     .        + AVECT(IX,IY,1,1)*CARAY(IX-1,IY,1,ISYS)
     .         + AVECT(IX,IY,1,2)*CARAY(IX+1,IY,1,ISYS)
     .          + AVECT(IX,IY,1,3)*CARAY(IX,IY-1,1,ISYS)
     .           + AVECT(IX,IY,1,4)*CARAY(IX,IY+1,1,ISYS)
     .            + AVECT(IX,IY,1,6)*CARAY(IX,IY,2,ISYS)
  119     CONTINUE
  120   CONTINUE
C      LAYERS 2 TO NZ-1
c$doacross local(iy,ix,iz,is,ie)
        DO 125 IZ=2,NZ-1
          DO 125 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 125
            DO 123 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 123
              DM(IX,IY,IZ) = - DIAG(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
     .         + AVECT(IX,IY,IZ,1)*CARAY(IX-1,IY,IZ,ISYS)
     .          + AVECT(IX,IY,IZ,2)*CARAY(IX+1,IY,IZ,ISYS)
     .           + AVECT(IX,IY,IZ,3)*CARAY(IX,IY-1,IZ,ISYS)
     .            + AVECT(IX,IY,IZ,4)*CARAY(IX,IY+1,IZ,ISYS)
     .             + AVECT(IX,IY,IZ,5)*CARAY(IX,IY,IZ-1,ISYS)
     .              + AVECT(IX,IY,IZ,6)*CARAY(IX,IY,IZ+1,ISYS)
  123     CONTINUE
  125   CONTINUE
C      BOTTOM LAYER
c$doacross local(iy,ix,is,ie)
        DO 130 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0)   GO TO 130
          DO 129 IX=IS,IE
            IF(FSM(IX,IY).LE.0.)  GO TO 129
            DM(IX,IY,NZ) =
     .        - DIAG(IX,IY,NZ)*CARAY(IX,IY,NZ,ISYS)
     .         + AVECT(IX,IY,NZ,1)*CARAY(IX-1,IY,NZ,ISYS)
     .          + AVECT(IX,IY,NZ,2)*CARAY(IX+1,IY,NZ,ISYS)
     .           + AVECT(IX,IY,NZ,3)*CARAY(IX,IY-1,NZ,ISYS)
     .            + AVECT(IX,IY,NZ,4)*CARAY(IX,IY+1,NZ,ISYS)
     .             + AVECT(IX,IY,NZ,5)*CARAY(IX,IY,NZ-1,ISYS)
  129     CONTINUE
  130   CONTINUE

      ELSE

        IZ=1
c$doacross local(iy,is,ie,ix)
        DO 132 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0)   GO TO 132
          DO 131 IX=IS,IE
            IF(FSM(IX,IY).LE.0.)  GO TO 131
            DM(IX,IY,IZ) = - DIAG(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
     .       + AVECT(IX,IY,IZ,1)*CARAY(IX-1,IY,IZ,ISYS)
     .        + AVECT(IX,IY,IZ,2)*CARAY(IX+1,IY,IZ,ISYS)
     .         + AVECT(IX,IY,IZ,3)*CARAY(IX,IY-1,IZ,ISYS)
     .          + AVECT(IX,IY,IZ,4)*CARAY(IX,IY+1,IZ,ISYS)
  131     CONTINUE
  132   CONTINUE

      ENDIF

C        ADD SOURCES/SINKS DUE TO KINETICS
c$doacross local(iy,ix,iz,is,ie)
c lb 20170326
c      OPEN(UNIT=428,FILE='FSM',FORM='UNFORMATTED')
c      WRITE(428) FSM
      DO 135 IZ=1,NZ
       DO 135 IY=1,NY
        IS = IXS(IY)
        IE = IXE(IY)
        IF(IS.EQ.0)   GO TO 135
        DO 134 IX=IS,IE
         DM(IX,IY,IZ) = DM(IX,IY,IZ) + CDARAY(IX,IY,IZ,ISYS)
  134   CONTINUE
  135 CONTINUE

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
      IF(NOPS(ISYS).GT.0)  THEN
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
  140  CONTINUE
      ENDIF
      IF(NONPS(ISYS).GT.0)   THEN
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
  145  CONTINUE
      ENDIF
      IF(NOFL(ISYS).GT.0)   THEN
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
  150  CONTINUE
      ENDIF
      IF(NOATM(ISYS).GT.0)   THEN
       DELTATM = TIME - NXATMT
       DO 155 IY=1,NY
        IS = IXS(IY)
        IE = IXE(IY)
        IF(IS.EQ.0)   GO TO 155
        DO 154 IX=IS,IE
         IF(FSM(IX,IY).LE.0.)  GO TO 154
         DM(IX,IY,1) = DM(IX,IY,1) + (DELTATM*SATM(IX,IY,ISYS)
     .               + BATM(IX,IY,ISYS))*XAZ(IX,IY)
         IF(IMBWK.EQ.1)  SYSLOADS(4,ISYS) = SYSLOADS(4,ISYS) +
     .      DTMB*(DELTATM*SATM(IX,IY,ISYS)+BATM(IX,IY,ISYS))*XAZ(IX,IY)
  154   CONTINUE
  155  CONTINUE
      ENDIF

C  CALCULATE C AT TIME LEVEL N+1
C    PREDICTOR STEP (UPWIND ADVECTION, KINETICS, LOADS)
 
c$doacross local(iz,iy,ix,is,ie)
        DO 160 IZ=1,NZ
          DO 160 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 160
            DO 159 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 159
              CONCS(IX,IY,IZ) = (BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
     .                           + DT*DM(IX,IY,IZ))/BVOLN(IX,IY,IZ)
  159       CONTINUE
  160   CONTINUE

C  BOUNDARY ELEMENTS ARE UNCHANGED
      DO 165 I=1,NOBCALL
        IX=IBCALL(1,I)
        IY=IBCALL(2,I)
        IZ=IBCALL(3,I)
        CONCS(IX,IY,IZ)=CARAY(IX,IY,IZ,ISYS)
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
     .        (AMAX1(0.,QX(IX,IY,IZ)*CARAY(IX-1,IY,IZ,ISYS))
     .       + AMIN1(0.,QX(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)))
         FLYMB(IX,IY,IZ,ISYS) = FLYMB(IX,IY,IZ,ISYS) + DTMB*
     .        (AMAX1(0.,QY(IX,IY,IZ)*CARAY(IX,IY-1,IZ,ISYS))
     .       + AMIN1(0.,QY(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)))
  170 CONTINUE
 
C***********************************************************************

C           SMOLARKIEWICZ CORRECTIONS
 
C-----------------------------------------------------------------------
C        SECOND-ORDER ACCURATE SCHEME (ISMOLAR.EQ.0)
  175 IF(ISMOLAR.EQ.0) THEN
C        CALCULATE ANTIDIFFUSIION FLOWS ( = VELOCITIES * AREAS)
c$doacross local(iz,iy,ix,is,ie)
       DO 180 IZ=1,NZ
        DO 180 IY=1,NY
         DO 180 IX=1,NX
          IF(FSM(IX,IY).EQ.0. .OR. FSM(IX,IY).EQ.-1.) GO TO 180
          IF(FSM(IX,IY).EQ.-2. .AND. ISMOLBCOPT.EQ.0) GO TO 180
            IF(CONCS(IX,IY,IZ).LE.1.0E-09 .OR.
     .            CONCS(IX-1,IY,IZ).LE.1.0E-09)  THEN
               ANTIQX(IX,IY,IZ) = 0.0
            ELSE IF(ABS(U(IX,IY,IZ)).LT.(U(IX,IY,IZ)*U(IX,IY,IZ)*
     .                 DT/DXBAR1(IX,IY)))   THEN
               ANTIQX(IX,IY,IZ) = 0.0
            ELSE
               ANTIQX(IX,IY,IZ) = ABS(U(IX,IY,IZ))*
     .          (1.-ABS(U(IX,IY,IZ))*DT/DXBAR1(IX,IY))*
     .               (CONCS(IX,IY,IZ)-CONCS(IX-1,IY,IZ))/
     .                    (CONCS(IX,IY,IZ)+CONCS(IX-1,IY,IZ)+R)
            ENDIF

            IF(CONCS(IX,IY,IZ).LE.1.0E-09 .OR.
     .            CONCS(IX,IY-1,IZ).LE.1.0E-09)  THEN
               ANTIQY(IX,IY,IZ) = 0.0
            ELSE IF(ABS(V(IX,IY,IZ)).LT.(V(IX,IY,IZ)*V(IX,IY,IZ)*
     .                 DT/DYBAR2(IX,IY)))   THEN
               ANTIQY(IX,IY,IZ) = 0.0
            ELSE
               ANTIQY(IX,IY,IZ) = ABS(V(IX,IY,IZ))*
     .          (1.-ABS(V(IX,IY,IZ))*DT/DYBAR2(IX,IY))*
     .               (CONCS(IX,IY,IZ)-CONCS(IX,IY-1,IZ))/
     .                    (CONCS(IX,IY,IZ)+CONCS(IX,IY-1,IZ)+R)
            ENDIF

          IF(FSM(IX,IY).GT.0.) THEN
            IF (IZ.GT.1) THEN
               IF(CONCS(IX,IY,IZ).LE.1.0E-09 .OR.
     .               CONCS(IX,IY,IZ-1).LE.1.0E-09)  THEN
                  ANTIQZ(IX,IY,IZ) = 0.0
               ELSE IF(ABS(W(IX,IY,IZ)).LT.(W(IX,IY,IZ)*W(IX,IY,IZ)*
     .                    DT/(D(IX,IY)*DZZ(IZ-1))))   THEN
                  ANTIQZ(IX,IY,IZ) = 0.0
               ELSE
                 ANTIQZ(IX,IY,IZ) = ABS(W(IX,IY,IZ))*
     .                (1.-ABS(W(IX,IY,IZ))*DT/(D(IX,IY)*DZZ(IZ-1)))*
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

C        ASSIGN CARAY AND CONCS TO BE 0.0 TEMPORALLY
        DO IZ=1,NZ 
        DO IY=1,NY
        DO IX=1,NX
          IF(FSM(IX,IY).EQ.0.)CARAY(IX,IY,IZ,ISYS)=0.0
          IF(FSM(IX,IY).EQ.0.)CONCS(IX,IY,IZ)=0.0
        ENDDO
        ENDDO
        ENDDO
        DO 191 IZ=1,NZ 
        DO 191 IY=1,NY
        DO 191 IX=1,NX
          FMAX=AMAX1(CARAY(IX,IY,IZ,ISYS),
     .     CARAY(MAX(1,IX-1),IY,IZ,ISYS),CARAY(MIN(NX,IX+1),IY,IZ,ISYS),
     .     CARAY(IX,MAX(1,IY-1),IZ,ISYS),CARAY(IX,MIN(NY,IY+1),IZ,ISYS),
     .     CARAY(IX,IY,MAX(1,IZ-1),ISYS),CARAY(IX,IY,MIN(NZ,IZ+1),ISYS),
     .     CONCS(IX,IY,IZ),
     .     CONCS(MAX(1,IX-1),IY,IZ),CONCS(MIN(NX,IX+1),IY,IZ),
     .     CONCS(IX,MAX(1,IY-1),IZ),CONCS(IX,MIN(NY,IY+1),IZ),
     .     CONCS(IX,IY,MAX(1,IZ-1)),CONCS(IX,IY,MIN(NZ,IZ+1)))
          FLUXIN(IX,IY,IZ)=(FMAX-CONCS(IX,IY,IZ))/
     .                 (DT*FLUXIN(IX,IY,IZ)+1.0E-15)
  191   CONTINUE

        DO IZ=1,NZ 
        DO IY=1,NY
        DO IX=1,NX
          IF(FSM(IX,IY).EQ.0.)CARAY(IX,IY,IZ,ISYS)=1.E+20
          IF(FSM(IX,IY).EQ.0.)CONCS(IX,IY,IZ)=1.E+20
        ENDDO
        ENDDO
        ENDDO
        DO 192 IZ=1,NZ 
        DO 192 IY=1,NY
        DO 192 IX=1,NX
          FMIN=AMIN1(CARAY(IX,IY,IZ,ISYS),
     .     CARAY(MAX(1,IX-1),IY,IZ,ISYS),CARAY(MIN(NX,IX+1),IY,IZ,ISYS),
     .     CARAY(IX,MAX(1,IY-1),IZ,ISYS),CARAY(IX,MIN(NY,IY+1),IZ,ISYS),
     .     CARAY(IX,IY,MAX(1,IZ-1),ISYS),CARAY(IX,IY,MIN(NZ,IZ+1),ISYS),
     .     CONCS(IX,IY,IZ),
     .     CONCS(MAX(1,IX-1),IY,IZ),CONCS(MIN(NX,IX+1),IY,IZ),
     .     CONCS(IX,MAX(1,IY-1),IZ),CONCS(IX,MIN(NY,IY+1),IZ),
     .     CONCS(IX,IY,MAX(1,IZ-1)),CONCS(IX,IY,MIN(NZ,IZ+1)))
          FLUXOUT(IX,IY,IZ)=(CONCS(IX,IY,IZ)-FMIN)/
     .                  (DT*FLUXOUT(IX,IY,IZ)+1.0E-15)
  192   CONTINUE
 
        DO 193 IZ=1,NZ 
        DO 193 IY=1,NY
        DO 193 IX=1,NX
          IF(FSM(IX,IY).EQ.0.)CARAY(IX,IY,IZ,ISYS)=0.0
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

c$doacross local(iz,iy,ix,is,ie)
      DO 200 IZ=1,NZ
        DO 200 IY=2,NY
          DO 200 IX=2,NX
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
  200 CONTINUE
 
C  SET HORIZONTAL ANTIDIFFUSION FLUXES IN BOUNDARY ELEMENTS TO ZERO
        IF(ISMOLBCOPT.EQ.0) THEN
          DO 210 I=1,NOBCALL
            IX=IBCALL(1,I)
            IY=IBCALL(2,I)
            IZ=IBCALL(3,I)
            QXFLX(IX,IY,IZ)=0.0
            QXFLX(IX+1,IY,IZ)=0.0
            QYFLX(IX,IY,IZ)=0.0
            QYFLX(IX,IY+1,IZ)=0.0
 210      CONTINUE
        ENDIF

C        TAKE CORRECTION STEP, PLACING RESULTS
C          BACK INTO CONCS NOW AT '(N+1)'
        IF (NZ.GT.1) THEN
c$doacross local(iz,iy,ix,is,ie)
          DO 230 IZ=1,NZ-1
           DO 230 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 230
            DO 229 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 229
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)
     .           +DT*(QXFLX(IX,IY,IZ)-QXFLX(IX+1,IY,IZ)
     .               +QYFLX(IX,IY,IZ)-QYFLX(IX,IY+1,IZ)
     .               +QZFLX(IX,IY,IZ+1)-QZFLX(IX,IY,IZ))/
     .                                           BVOLN(IX,IY,IZ)
 229        CONTINUE
 230      CONTINUE

          IZ=NZ
c$doacross local(iy,ix,is,ie)
          DO 235 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 235
            DO 234 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 234
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)
     .           +DT*(QXFLX(IX,IY,IZ)-QXFLX(IX+1,IY,IZ)
     .               +QYFLX(IX,IY,IZ)-QYFLX(IX,IY+1,IZ)
     .               -QZFLX(IX,IY,IZ))/BVOLN(IX,IY,IZ)
 234        CONTINUE
 235      CONTINUE
        ELSE
          IZ=1
          DO 240 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 240
            DO 239 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 239
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)
     .           +DT*(QXFLX(IX,IY,IZ)-QXFLX(IX+1,IY,IZ)
     .               +QYFLX(IX,IY,IZ)-QYFLX(IX,IY+1,IZ))/
     .                                             BVOLN(IX,IY,IZ)
 239      CONTINUE
 240      CONTINUE
        ENDIF

C  END OF SMOLARKIEWICZ CORRECTION CODE

C        CHECK FOR MASS/FLUX BALANCE COMPUTATIONS
      IF(MASSBAL.EQ.0) GO TO 333
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0 .AND. ITIMESECS.NE.NXPRTMB) GO TO 333
      IF(ITIMESECS.LT.ISMBSECS) GO TO 333
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
  333   IF(SLCOPT.EQ.1 .AND.
     .      MOD((ITIMESECS-IFIX(86400.*TZERO)),IDTSLCSECS).LE.0.) THEN
          CALL RCAMPROF(CARAY(1,1,1,ISYS),CSLAVE(1,1,1,ISYS))
        ENDIF

c$doacross local(iz,iy,ix,is,ie)
        DO 250 IZ=1,NZ
        DO 250 IY=2,NY
        DO 250 IX=2,NX
              IF(FSM(IX,IY).GE.1. .OR. FSM(IX,IY).LE.-2.)THEN
                IF(FSM(IX-1,IY).GE.1. .OR. FSM(IX-1,IY).LE.-2.)THEN
                  RXFLX(IX,IY,IZ)=RX(IX,IY,IZ)*
     .             ((CARAY(IX-1,IY,IZ,ISYS)-CSLAVE(IX-1,IY,IZ,ISYS)) -
     .             (CARAY(IX,IY,IZ,ISYS)-CSLAVE(IX,IY,IZ,ISYS)))
                ENDIF
                IF(FSM(IX,IY-1).GE.1. OR. FSM(IX,IY-1).LE.-2.)THEN
                  RYFLX(IX,IY,IZ)=RY(IX,IY,IZ)*
     .             ((CARAY(IX,IY-1,IZ,ISYS)-CSLAVE(IX,IY-1,IZ,ISYS)) -
     .             (CARAY(IX,IY,IZ,ISYS)-CSLAVE(IX,IY,IZ,ISYS)))
                ENDIF
              ENDIF
 250    CONTINUE
 
C  CALCULATE CONCENTRATIONS: HORIZONTAL DIFFUSION
c$doacross local(iz,iy,ix,is,ie)
        DO 260 IZ=1,NZ
          DO 260 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 260
            DO 259 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 259
              CONCS(IX,IY,IZ)=CONCS(IX,IY,IZ)+DT*
     .           (RXFLX(IX,IY,IZ)-RXFLX(IX+1,IY,IZ)+RYFLX(IX,IY,IZ)
     .           -RYFLX(IX,IY+1,IZ))/BVOLN(IX,IY,IZ)
 259      CONTINUE
 260    CONTINUE

C  BOUNDARY ELEMENTS ARE UNCHANGED
        DO 265 I=1,NOBCALL
          IX=IBCALL(1,I)
          IY=IBCALL(2,I)
          IZ=IBCALL(3,I)
          CONCS(IX,IY,IZ)=CARAY(IX,IY,IZ,ISYS)
  265   CONTINUE

C        CHECK FOR MASS/FLUX BALANCE COMPUTATIONS
      IF(MASSBAL.EQ.0) GO TO 280
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0 .AND. ITIMESECS.NE.NXPRTMB) GO TO 280
      IF(ITIMESECS.LT.ISMBSECS) GO TO 280
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
 
C        IMPLICIT VERTICAL DIFFUSION STEP  (IF NZ > 1)
  280 IF (NZ.EQ.1)  GO TO 380
c$doacross local(iy,ix,iz,is,ie)
      DO 340 IZ=1,NZ
       DO 340 IY=1,NY
        IS = IXS(IY)
        IE = IXE(IY)
        IF(IS.EQ.0)   GO TO 340
        DO 339 IX=IS,IE
          MF(IX,IY,IZ) = CONCS(IX,IY,IZ)*BVOLN(IX,IY,IZ)
 339    CONTINUE
 340  CONTINUE

C        SET UP DIAGONAL AND OFF-DIAGONAL ELEMENTS
c$doacross local(iz,iy,ix,is,ie)
      DO 350 IY=1,NY
       IS = IXS(IY)
       IE = IXE(IY)
       IF(IS.EQ.0)   GO TO 350
       DO 349 IX=IS,IE
        IF(FSM(IX,IY).LE.0.)  GO TO 349
        AL(IX,IY,1) = 0.
        AD(IX,IY,1) = BVOLN(IX,IY,1) + DT*RZ(IX,IY,2)
        AU(IX,IY,1) = -DT*RZ(IX,IY,2)
        DO 345 IZ=2,NZ-1
          AL(IX,IY,IZ) = -DT*RZ(IX,IY,IZ)
          AD(IX,IY,IZ) = BVOLN(IX,IY,IZ) +
     .                      DT*(RZ(IX,IY,IZ) + RZ(IX,IY,IZ+1))
          AU(IX,IY,IZ) = -DT*RZ(IX,IY,IZ+1)
  345   CONTINUE
        AL(IX,IY,NZ) = -DT*RZ(IX,IY,NZ)
        AD(IX,IY,NZ) = BVOLN(IX,IY,NZ) + DT*RZ(IX,IY,NZ)
        AU(IX,IY,NZ) = 0.0
  349  CONTINUE
  350 CONTINUE

C        FORWARD SWEEP OF TRIDIAGONAL SCHEME
c$doacross local(iz,iy,ix,is,ie)
      DO 360 IY=1,NY
       IS = IXS(IY)
       IE = IXE(IY)
       IF(IS.EQ.0)   GO TO 360
       DO 359 IX=IS,IE
        IF(FSM(IX,IY).LE.0.)  GO TO 359
        BETA(IX,IY,1) = AD(IX,IY,1)
        GAMMA(IX,IY,1) = MF(IX,IY,1)/BETA(IX,IY,1)
        DO 355 IZ=2,NZ
         BETA(IX,IY,IZ) = AD(IX,IY,IZ)
     .               - AL(IX,IY,IZ)*AU(IX,IY,IZ-1)/BETA(IX,IY,IZ-1)
         GAMMA(IX,IY,IZ) = (MF(IX,IY,IZ)
     .               - AL(IX,IY,IZ)*GAMMA(IX,IY,IZ-1))/BETA(IX,IY,IZ)
  355   CONTINUE
  359  CONTINUE
  360 CONTINUE

C        BACKWARD SWEEP FOR SOLUTION
c$doacross local(iz,iy,ix,is,ie,ilyr)
      DO 370 IY=1,NY
       IS = IXS(IY)
       IE = IXE(IY)
       IF(IS.EQ.0)   GO TO 370
       DO 369 IX=IS,IE
        IF(FSM(IX,IY).LE.0.)  GO TO 369
        CONCS(IX,IY,NZ) = GAMMA(IX,IY,NZ)
        DO 365 ILYR=2,NZ
         IZ = NZ+1-ILYR
         CONCS(IX,IY,IZ) = GAMMA(IX,IY,IZ)
     .           - AU(IX,IY,IZ)*CONCS(IX,IY,IZ+1)/BETA(IX,IY,IZ)
  365   CONTINUE
  369  CONTINUE
  370 CONTINUE

C  BOUNDARY ELEMENTS ARE UNCHANGED
      DO 375 I=1,NOBCALL
        IX=IBCALL(1,I)
        IY=IBCALL(2,I)
        IZ=IBCALL(3,I)
        CONCS(IX,IY,IZ)=CARAY(IX,IY,IZ,ISYS)
  375 CONTINUE
 
  380 CONTINUE

C        CHECK TO SEE IF NEGATIVE SOLUTIONS PERMITTED 
c$doacross local(iz,iy,ix,is,ie,inegs,totmas) , share(negseg)
      DO 390 IZ=1,NZ
        DO 390 IY=1,NY
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0)   GO TO 390
          DO 389 IX=IS,IE
            IF(FSM(IX,IY).LE.0.)  GO TO 389
            IF (NEGSLN.EQ.1)   THEN
              CARAY(IX,IY,IZ,ISYS) = CONCS(IX,IY,IZ)
            ELSE
C        IF NOT ... CHECK TO MAKE SURE SOLUTION IS NOT
C                   FORCED THRU ZERO (1.0E-38)
              IF(CONCS(IX,IY,IZ).GE.CMIN(ISYS)) THEN
                CARAY(IX,IY,IZ,ISYS) = CONCS(IX,IY,IZ)
              ELSE
C          PREVENT NEG SOLN BY QUARTERING PRESENT CONC AND PRINT WARNING
                CARAY(IX,IY,IZ,ISYS) = 
     .                   MAX(CARAY(IX,IY,IZ,ISYS)/4.,CMIN(ISYS))
C$              call m_lock
                INEGS = INEGS+1
                IF(INEGS.LE.99)  THEN
                  NEGSEG(INEGS,1) = IX
                  NEGSEG(INEGS,2) = IY
                  NEGSEG(INEGS,3) = IZ
                ENDIF
                TOTMAS = TOTMAS + BVOLN(IX,IY,IZ)*
     .                   (CARAY(IX,IY,IZ,ISYS)-CONCS(IX,IY,IZ))
                ADDMASS(IX,IY,IZ,ISYS) = ADDMASS(IX,IY,IZ,ISYS) +
     .           BVOLN(IX,IY,IZ)*(CARAY(IX,IY,IX,ISYS)-CONCS(IX,IY,IZ))
                IWRTADDM=1
C$              call m_unlock
              ENDIF
            ENDIF
 389    CONTINUE
 390  CONTINUE

       IF(IADDMDEBUG.EQ.1) THEN
        IF(INEGS.GT.0)  THEN
C$            call mp_setlock
          CALL RCAMESS(2,0.001*TOTMAS)
          WRITE(OUT,1000)  ((NEGSEG(I,J),J=1,3),I=1,MIN(INEGS,99))
 1000     FORMAT(30X,12I5)
C$            call mp_unsetlock
        ENDIF
       ENDIF

  400 CONTINUE

C        SEE IF MASS HAS BEEN ADDED TO THE SYSTEM
      IF(IWRTADDM.EQ.1) THEN
        OPEN(19,FILE='RCAFADDM',FORM='UNFORMATTED',STATUS='REPLACE')
        WRITE(19)  ADDMASS
        CLOSE(19)
        IWRTADDM=0
      ENDIF

C        FLUX BALANCES FOR MISC KINETIC VARIABLES
      IF(MASSBAL.EQ.0 .OR. NOKINSYS.EQ.0) GO TO 410
C        INSTANTANEOUS OR AVERAGED
      IF(IMBDOPT.EQ.0) THEN 
       IF(ITIMESECS.NE.NXPRTMB) GO TO 410
       DO 405 ISYS=1,NOKINSYS
        DO 405 IZ=1,NZ
         DO 405 IY=1,NY
          DO 405 IX=1,NX
           IF(FSM(IX,IY).NE.1.) GO TO 405
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
  405  CONTINUE
      ELSE
       IF(ITIMESECS.LT.ISMBSECS) GO TO 410
       DO 407 ISYS=1,NOKINSYS
        DO 407 IZ=1,NZ
         DO 407 IY=1,NY
          DO 407 IX=1,NX
           IF(FSM(IX,IY).NE.1.) GO TO 407
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
  407  CONTINUE
      ENDIF 

c$doacross local(iz,iy,ix)
  410 DO 415 IZ=1,NZ
       DO 415 IY=1,NY
        DO 415 IX=1,NX
         BVOL(IX,IY,IZ)=BVOLN(IX,IY,IZ)
  415 CONTINUE
C    LDH FOR OUTPUT SETTING
       IF(ITIMESECS.EQ.OUTPUTNUM*86400)THEN       
       IRECORD=IRECORD+1
       WRITE(FILESTR,'(I2)')IRECORD
       WRITE(SIC,'(I2)')IRECORD 
       CLOSE(11)
       CLOSE(524)
       CLOSE(525)
       CLOSE(14)            
       OPEN(UNIT=11,FILE=TRIM('RCAF11')//'_'//TRIM(ADJUSTL(FILESTR))
     .      ,FORM='UNFORMATTED')
       OPEN(UNIT=524,FILE=TRIM('RCAFLB')//'_'//TRIM(ADJUSTL(FILESTR))
     .      ,FORM='UNFORMATTED')
       OPEN(UNIT=525,FILE=TRIM('SEDFLB')//'_'//TRIM(ADJUSTL(FILESTR))
     .      ,FORM='UNFORMATTED')
       OPEN(UNIT=14,FILE=TRIM('RCAF14')//'_'//TRIM(ADJUSTL(FILESTR))
     .      ,FORM='UNFORMATTED')
      IF(NDMPS.GT.0)THEN
      CLOSE(13)
      OPEN(UNIT=13,FILE=TRIM('RCAF13')//'_'//TRIM(ADJUSTL(FILESTR))
     .      ,FORM='UNFORMATTED')       
      ENDIF
      OPEN(UNIT=1515,FILE=TRIM('RESTART_RCAFIC')//'_'//TRIM(ADJUSTL
     .      (SIC)),FORM='UNFORMATTED')
      WRITE(1515)CARAY      
      CLOSE(1515)      
      OUTPUTNUM=OUTPUTNUM+OUTPUTNUM
      WRITE(OUT,1333)'CREATING FILE OF RCAF11,RCAF13 AND RCAF14 AT'
     .       ,OUTPUTNUM ,'DAYS FOR FILESTR',IRECORD
 1333 FORMAT(/,A45,I3,A16,I2,/)
      ENDIF

C  PREPARE FOR NEXT TIME STEP
      ITIMESECS = ITIMESECS + IDTSECS
      TIME = ITIMESECS/86400.
      WRITE(OUT,'(A,10X,I15)')'ITIMESECS = ',ITIMESECS
C        CHECK TO SEE IF IT IS TIME TO PRINT RESULTS
C  FIRST DETAILED DUMPS
      IF(ITIMESECS.LT.NXPRTD)  GO TO 435
C  DETAILED PRINT TIME
C  SET IDISK FOR DETAILED DUMPS
      IDISK=2
      NXPRTD = NXPRTD + IPRNTDSECS
C  SECOND GLOBAL DUMPS
 435  IF(ITIMESECS.LT.NXPRTG)  GO TO 450
C  GLOBAL PRINT TIME
C  CALL PRINT ROUTINE
      CALL RCA09
      IDISK=IDISK+1
      NXPRTG = NXPRTG + IPRNTGSECS

C  NOW PERFORM STABILITY CHECK
      DO 445 ISYS=1,NOSYS 
        IF(SYSBY(ISYS).EQ.1)   GO TO 445
        CKMAX = CMAX(ISYS)
        DO 440 IZ=1,NZ
          DO 440 IY=1,NY
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 440
            DO 439 IX=IS,IE
              IF(FSM(IX,IY).LE.0.)  GO TO 439
              IF(CARAY(IX,IY,IZ,ISYS).GE.CKMAX)   GO TO 500
  439       CONTINUE
  440   CONTINUE
  445 CONTINUE

C  MASS/FLUX BALANCE DUMPS
  450 IF(MASSBAL.EQ.0) GO TO 470
      IF(ITIMESECS.LT.NXPRTMB) GO TO 470
       IF(IMBDOPT.EQ.0) THEN
        WRITE(17)  TIME,SYSMASS,SYSLOADS,FLXMB,FLYMB
       ELSE
        TIMEAVE=FLOAT(IPRNTMBSECS)/86400.
        DO 455 ISYS=1,NOSYS+NOKINSYS
         DO 455 IZ=1,NZ
          DO 455 IY=1,NY
           DO 455 IX=1,NX
            FLXMB(IX,IY,IZ,ISYS)=FLXMB(IX,IY,IZ,ISYS)/TIMEAVE
            FLYMB(IX,IY,IZ,ISYS)=FLYMB(IX,IY,IZ,ISYS)/TIMEAVE
  455   CONTINUE
        DO ISYS=1,NOSYS
         DO I=1,4
          SYSLOADS(I,ISYS)=SYSLOADS(I,ISYS)/TIMEAVE
         ENDDO
        ENDDO
        WRITE(17)  TIME,SYSMASS,SYSLOADS,FLXMB,FLYMB
       ENDIF
C        RESET ARRAYS AND COUNTER
      DO 460 ISYS=1,NOSYS
       DO 460 I=1,4
        SYSLOADS(I,ISYS)=0.
  460 CONTINUE
      DO 465 ISYS=1,NOSYS+NOKINSYS
       DO 465 IZ=1,NZ
        DO 465 IY=1,NY
         DO 465 IX=1,NX
          FLXMB(IX,IY,IZ,ISYS)=0.0
          FLYMB(IX,IY,IZ,ISYS)=0.0
  465 CONTINUE
      NXPRTMB = NXPRTMB + IPRNTMBSECS
      IF(NXPRTMB.GT.IEMBSECS) MASSBAL=0

C        CHECK TO SEE IF NECESSARY TO UPDATE TIME-VARIABLE FUNCTIONS
C        FIRST POINT SOURCE LOADS 
  470 IF(IPSOPT.GT.1.AND.TIME.GE.NXPST)   
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
      !IF(IBCOPT.GT.1.AND.TIME.GE.NXBCT) CALL RCA11    ! ---original code
      IF(IBCOPT.EQ.2.OR.IBCOPT.EQ.4.AND.TIME.GE.NXBCT) CALL RCA11      ! by ldh
C        NEXT TRANSPORT FIELDS
C     IF(ITIMESECS.GE.HYDBRK(ITIMHYD))  THEN
      IF(ITIMESECS.EQ.NXHYDTSECS)  THEN    
        IF(86400*TBRK(ITRAK).GT.ITIMESECS) CALL RCA03A
        CALL RCAEXP1
C        ITIMHYD = ITIMHYD+1
      ENDIF      
C        CHECK FOR END OF SIMULATION...
C           IF NOT GO BACK AND TAKE NEXT INTEGRATION STEP             
      IF(TIME.LE.TEND)   GO TO 30 
      ITRAK = ITRAK+1 
      IF(ITRAK.LE.NSTEP)  THEN
         IDTSECS=ISTEP(ITRAK)
         DT=ISTEP(ITRAK)/86400.
         TEND = TBRK(ITRAK) 
         IF(DT.EQ.0.0)   THEN
            CALL RCAMESS(3,DT)
            IDISK=1
            INITB=2
            CALL TUNER
            CALL EXIT
         ENDIF
      ELSE
         TEND = 0.0
      ENDIF 

      IF(TEND.GT.0.)   GO TO 30 
C        FINISHED...
      WRITE(OUT,8000)
 8000 FORMAT(///30X,'RCAEXPS FINISHED INTEGRATION'/
     .          30X,'USER DUMPS TO FOLLOW'//)
      IF(MASSBAL.EQ.1) THEN
        DEALLOCATE(FLXMB)
        DEALLOCATE(FLYMB)
      ENDIF
      RETURN

C         STABILITY CRITERIA VIOLATED...ABEND 
  500 CALL RCAMESS(1,CKMAX) 
      WRITE(OUT,9000)  BVOLN(IX,IY,IZ),DIAG(IX,IY,IZ)
 9000 FORMAT(10X,'VOLUME =',E13.4,'  ,DIAGONAL =',E13.4)
      IREC = IREC+1 
      IDISK = 1 
      INITB = 2
      CALL TUNER
      RETURN
      END 
