      SUBROUTINE RCAMPROF(CSIGMA,CSLAVE)
CC
C        RCAMPROF COMPUTES DOMAIN-AVERAGED WATER QUALITY CONCENTRATIONS
C                 FOR USE IN SIGMA-COORDINATE CORRECTION
C

      SAVE
      INCLUDE 'RCACM'
      REAL  CSIGMA(NX,NY,NZ),CSLAVE(NX,NY,NZ)
      REAL  ZM(NZ),CI(NZ+1),CA(NSLC),CSL(NX,NY,NSLC)

C        CSIGMA - INPUT ARRAY OF SIGMA-LEVEL CONCENTRATIONS
C        CSLAVE - OUTPUT ARRAY OF SIGMA-LEVEL CONCENTRATIONS
C                   MINUS THE DOMAIN-AVERAGED CONCNETRATIONS,
C                   I.E. THE SIGMA-LEVEL RESIDULALS

C        INTERPOLATE CSIGMA ONTO STANDARD LEVELS

ccc c$doacross local(iz,iy,ix) , share(ci,csl,zm)
      DO 40 IY=1,NY
        DO 30 IX=1,NX
          IF (FSM(IX,IY).NE.0.0) THEN
            DO 10 IZ=1,NZ
              ZM(IZ) = -ZZ(IZ) * HBAR(IX,IY)
              CI(IZ) = CSIGMA(IX,IY,IZ)
   10       CONTINUE
            CALL SINTER(ZM,CI,SLCDPTH,CA,NZ,NOSLC)
            DO 20 IZ=1,NOSLC
              CSL(IX,IY,IZ) = CA(IZ)
   20       CONTINUE
          END IF
   30   CONTINUE
   40 CONTINUE

C        FIND MEAN CONC AT EACH STANDARD LEVEL, BY AREALLY
C        INTEGRATING AND DIVIDING BY THE AREA AT THAT DEPTH

ccc c$doacross local(iz,iy,ix,srfarea) , share(ca)
      DO 70 IZ=1,NOSLC
        CA(IZ) = 0.0
        SRFAREA = 0.0
        DO 60 IY=1,NY
          DO 50 IX=1,NX
            IF (FSM(IX,IY).NE.0.0) THEN
              IF(-H(IX,IY).LE.SLCDPTH(IZ)) THEN
                SRFAREA = SRFAREA + XAZ(IX,IY)
                CA(IZ) = CA(IZ) + CSL(IX,IY,IZ)*XAZ(IX,IY)
              END IF
            END IF
   50     CONTINUE
   60   CONTINUE
        IF (SRFAREA.GT.0.0)  CA(IZ) = CA(IZ)/SRFAREA
   70 CONTINUE

C        INTERPOLATE DOMAIN-AVERAGED PROFILE BACK ONTO SIGMA GRID
  
ccc c$doacross local(iz,iy,ix) , share(cslave)
      DO 140 IY=1,NY
        DO 130 IX=1,NX
          IF (FSM(IX,IY).NE.0.0) THEN
            DO 110 IZ=1,NZ
              ZM(IZ) = -ZZ(IZ) * HBAR(IX,IY)
  110       CONTINUE
            CALL SINTER(SLCDPTH,CA,ZM,CI,NOSLC,NZ)
            DO 120 IZ=1,NZ
              CSLAVE(IX,IY,IZ) = CI(IZ)
  120       CONTINUE
          END IF
  130   CONTINUE
  140 CONTINUE

c     call exit
      RETURN
      END
