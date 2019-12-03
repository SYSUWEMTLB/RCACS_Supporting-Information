      SUBROUTINE  RCAEXP1
C 
C        RCAEXP1 SETS UP THE -DIAG- AND -AVECT- VECTORS
C                 IF SMOLARKIEWICZ CORRECTION THEN -DIAG- AND -AVECT-
C                 CONTAIN JUST ADVECTIVE TERMS
C                 IF NO SMOLARKIEWICZ CORRECTION THEN -DIAG- AND -AVECT-
C                 CONTAIN ADVECTIVE AND HORIZONTAL DISPERSIVE TERMS
C 
      SAVE
      INCLUDE  'RCACM'

      DO 100 IT=1,6
       DO 100 IZ=0,NZ+1
        DO 100 IY=1,NY
         DO 100 IX=1,NX
          AVECT(IX,IY,IZ,IT) = 0.
 100  CONTINUE

      DO 110 IZ=1,NZ
       DO 110 IY=1,NY
        DO 110 IX=1,NX
         DIAG(IX,IY,IZ) = 0.
 110  CONTINUE

C        CHECK FOR SMOLARKIEWICZ OPTION
      IF(INTGRTYP.LT.4) THEN
       IF(NZ.GT.1) THEN
        DO 150 IZ=1,NZ
          DO 150 IY=2,NY-1
            IS = IXS(IY)
            IE = IXE(IY)
            IF(IS.EQ.0)   GO TO 150
            DO 149 IX=IS,IE
              IF (FSM(IX,IY).LE.0.0)  GO TO 149
              IF (FSM(IX-1,IY).NE.0.0) AVECT(IX,IY,IZ,1)=RX(IX,IY,IZ)
              IF (FSM(IX+1,IY).NE.0.0) AVECT(IX,IY,IZ,2)=RX(IX+1,IY,IZ)
              IF (FSM(IX,IY-1).NE.0.0) AVECT(IX,IY,IZ,3)=RY(IX,IY,IZ)
              IF (FSM(IX,IY+1).NE.0.0) AVECT(IX,IY,IZ,4)=RY(IX,IY+1,IZ)
              DIAG(IX,IY,IZ)=0.0
              IF (FSM(IX-1,IY).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                                RX(IX,IY,IZ)
              IF (FSM(IX+1,IY).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                                RX(IX+1,IY,IZ)
              IF (FSM(IX,IY-1).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                                RY(IX,IY,IZ)
              IF (FSM(IX,IY+1).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                                RY(IX,IY+1,IZ)
 149        CONTINUE
 150    CONTINUE
       ELSE
        IZ=1
        DO 160 IY=2,NY-1
          IS = IXS(IY)
          IE = IXE(IY)
          IF(IS.EQ.0)   GO TO 160
          DO 159 IX=IS,IE
            IF (FSM(IX,IY).LE.0.0)  GO TO 159
            IF (FSM(IX-1,IY).NE.0.0) AVECT(IX,IY,IZ,1)=RX(IX,IY,IZ)
            IF (FSM(IX+1,IY).NE.0.0) AVECT(IX,IY,IZ,2)=RX(IX+1,IY,IZ)
            IF (FSM(IX,IY-1).NE.0.0) AVECT(IX,IY,IZ,3)=RY(IX,IY,IZ)
            IF (FSM(IX,IY+1).NE.0.0) AVECT(IX,IY,IZ,4)=RY(IX,IY+1,IZ)
            DIAG(IX,IY,IZ)=0.0
            IF (FSM(IX-1,IY).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                              RX(IX,IY,IZ)
            IF (FSM(IX+1,IY).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                              RX(IX+1,IY,IZ)
            IF (FSM(IX,IY-1).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                              RY(IX,IY,IZ)
            IF (FSM(IX,IY+1).NE.0.0) DIAG(IX,IY,IZ)=DIAG(IX,IY,IZ)+
     +                                              RY(IX,IY+1,IZ)
 159      CONTINUE
 160    CONTINUE
       ENDIF
      ENDIF
C
      DO 200 IZ=1,NZ
       DO 200 IY=2,NY-1
        IS = IXS(IY)
        IE = IXE(IY)
        IF(IS.EQ.0)   GO TO 200
        DO 199 IX=IS,IE

         IF (FSM(IX,IY).LE.0.0)  GO TO 199

         IF(QX(IX,IY,IZ).GT.0.)  THEN
            AVECT(IX,IY,IZ,1) = AVECT(IX,IY,IZ,1) + QX(IX,IY,IZ)
         ELSE
            DIAG(IX,IY,IZ) = DIAG(IX,IY,IZ) - QX(IX,IY,IZ)
         ENDIF

         IF(QX(IX+1,IY,IZ).GT.0.)  THEN
            DIAG(IX,IY,IZ) = DIAG(IX,IY,IZ) + QX(IX+1,IY,IZ)
         ELSE
            AVECT(IX,IY,IZ,2) = AVECT(IX,IY,IZ,2) - QX(IX+1,IY,IZ)
         ENDIF

         IF(QY(IX,IY,IZ).GT.0.)  THEN
            AVECT(IX,IY,IZ,3) = AVECT(IX,IY,IZ,3) + QY(IX,IY,IZ)
         ELSE
            DIAG(IX,IY,IZ) = DIAG(IX,IY,IZ) - QY(IX,IY,IZ)
         ENDIF

         IF(QY(IX,IY+1,IZ).GT.0.)  THEN
            DIAG(IX,IY,IZ) = DIAG(IX,IY,IZ) + QY(IX,IY+1,IZ)
         ELSE
            AVECT(IX,IY,IZ,4) = AVECT(IX,IY,IZ,4) - QY(IX,IY+1,IZ)
         ENDIF

         IF (NZ.EQ.1)  GO TO 199

         IF(QZ(IX,IY,IZ).GT.0.)  THEN
            DIAG(IX,IY,IZ) = DIAG(IX,IY,IZ) + QZ(IX,IY,IZ)
         ELSE
            AVECT(IX,IY,IZ,5) = AVECT(IX,IY,IZ,5) - QZ(IX,IY,IZ)
         ENDIF

         IF(QZ(IX,IY,IZ+1).GT.0.)  THEN
            AVECT(IX,IY,IZ,6) = AVECT(IX,IY,IZ,6) + QZ(IX,IY,IZ+1)
         ELSE
            DIAG(IX,IY,IZ) = DIAG(IX,IY,IZ) - QZ(IX,IY,IZ+1)
         ENDIF

  199   CONTINUE
  200 CONTINUE
C
C         ZERO TRANSPORT TERMS FOR BOUNDARY SEGMENTS
C
      DO 210 I=1,NOBCALL
        IX = IBCALL(1,I)
        IY = IBCALL(2,I)
        IZ = IBCALL(3,I)
        DIAG(IX,IY,IZ) = 0.0
        DO 220 IT=1,6
          AVECT(IX,IY,IZ,IT) = 0.0
  220   CONTINUE
  210 CONTINUE
 
      RETURN
      END
