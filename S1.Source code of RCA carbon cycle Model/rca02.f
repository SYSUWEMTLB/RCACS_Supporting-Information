      SUBROUTINE RCA02
C 
C        RCA02 READS PRINT CONTROL PARAMETERS AND INTEGR HISTORY
C 
      SAVE
      INCLUDE  'RCACM'
      INTEGER*2  DDGDPBY(40)
      CHARACTER*4  TWARPSLC
      DATA  DDGDPBY/40*1/

C        Integration Option and Integration Timestep Information

C        INTGRTYP -  INTEGRATION METHOD
C                     = 1  -  EXPLICIT - UPWIND
C                     = 3  -  SPLIT TIME STEP - UPWIND
C                     = 4  -  EXPLICIT - UPWIND/SMOLARKIEWICZ
C                     = 5  -  LEAPFROG - UPWIND/SMOLARKIEWICZ
C                     = 6  -  SPLIT TIME STEP - UPWIND/SMOLARKIEWICZ
C           NEGSLN =  0  POSITIVE SOLUTIONS ONLY
C                     1  NEGATIVE SOLUTIONS PERMITTED
C           SLCOPT =  0  DO NOT PERFORM SIGMA-LEVEL CORRECTION FOR
C                        HORIZONTAL DIFFUSION
C                     1  PERFORM SIGMA-LEVEL CORRECTION FOR
C                        HORIZONTAL DIFFUSION
C          ISMOLAR =  0  "ECOM3D" SMOLAR_2 - SECOND ORDER ACCURATE
C                     1  "ECOM3D" SMOLAR_R - RECURSIVE VERSION
C       ISMOLBCOPT =  0  DO NOT APPLY SMOLARKIEWICZ CORRECTIONS AT BOUNDARIES
C                     1  APPLY SMOLARKIEWICZ CORRECTIONS AT BOUNDARIES

      READ(IN,1000) COMMENT
 1000 FORMAT(A)
      READ(IN,2000,ERR=950) INTGRTYP,NEGSLN,SLCOPT,ISMOLAR,ISMOLBCOPT
 2000 FORMAT(5I10) 
      IF(INTGRTYP.EQ.1)    WRITE(OUT,2010) 
 2010 FORMAT(////10X,'EXPLICIT INTEGRATION SCHEME SELECTED')
      IF(INTGRTYP.EQ.2) THEN
         WRITE(OUT,2020) 
 2020    FORMAT(////10X,'LEAFFROG SCHEME NOT A VALID OPTION'/
     .              10X,'RCA TERMINATED')
        CALL EXIT
      ENDIF
      IF(INTGRTYP.EQ.3)    WRITE(OUT,2030) 
 2030 FORMAT(////10X,'SPLIT TIMESTEP INTEGRATION SCHEME SELECTED')
      IF(INTGRTYP.EQ.4)    WRITE(OUT,2040) 
 2040 FORMAT(////10X,'EXPLICIT-SMOLARKIEWICZ SCHEME SELECTED')
      IF(INTGRTYP.EQ.5)    WRITE(OUT,2050) 
 2050 FORMAT(////10X,'LEAPFROG-SMOLARKIEWICZ SELECTED')
      IF(INTGRTYP.EQ.6)    WRITE(OUT,2060) 
 2060 FORMAT(////10X,'SPLIT TIMESTEP-SMOLARKIEWICZ SCHEME SELECTED')
      IF(ISMOLAR.EQ.0) WRITE(OUT,2061)
 2061 FORMAT(10X,'SECOND ORDER ACCURATE SCHEME SELECTED')
      IF(ISMOLAR.EQ.1) WRITE(OUT,2062)
 2062 FORMAT(10X,'RECURSIVE SCHEME SELECTED')
      IF(ISMOLBCOPT.EQ.0) WRITE(OUT,2063)
 2063 FORMAT(10X,
     .   'SMOLARKIEWICZ CORRECTIONS WILL NOT BE APPLIED AT BOUNDARIES')
      IF(ISMOLBCOPT.EQ.1) WRITE(OUT,2064)
 2064 FORMAT(10X,
     .   'SMOLARKIEWICZ CORRECTIONS WILL BE APPLIED AT BOUNDARIES')
      IF(NEGSLN.EQ.1)   WRITE(OUT,2070) 
 2070 FORMAT(/10X,'USER HAS CHOSEN TO PERMIT NEGATIVE SOLUTIONS'/)
      IF(SLCOPT.EQ.0)    WRITE(OUT,2080) 
 2080 FORMAT(/10X,'NO SIGMA-LEVEL CORRECTION FOR HORIZONTAL DIFFUSION')
      IF(SLCOPT.EQ.1)  THEN
        WRITE(OUT,2090) 
 2090   FORMAT(///10X,'SIGMA-LEVEL CORRECTION FOR HORIZONTAL',
     .    ' DIFFUSION ENABLED')
        READ(IN,1000)   COMMENT
        READ(IN,2100)  IDTSLCSECS,TWARPSLC,NOSLC,(SLCDPTH(I),I=1,NOSLC)
 2100   FORMAT(I10,6X,A4,I10/(8F10.0))
        ISCALT=IUNITCHECK(TWARPSLC,'TWARPSLC')
        IDTSLCSECS=IDTSLCSECS*ISCALT
        DTSLC=IDTSLCSECS/86400.
        WRITE(OUT,2110)  IDTSLCSECS,DTSLC,(SLCDPTH(I),I=1,NOSLC)
 2110   FORMAT(10X,'TIME INTERVAL BETWEEN CALCULATING SIGMA-LEVEL',
     .    ' DOMAIN AVERAGES = ',I5,' SECONDS (= ',F7.5,' DAYS)'/10X,
     .    'STANDARD LEVEL DEPTHS TO BE USED FOR DOMAIN AVERAGING'/
     .    (10X,10F7.2))
        DO 20 I=1,NOSLC
          SLCDPTH(I) = -SLCDPTH(I)
   20   CONTINUE
      ENDIF

C     Time Wrap Scale Fact, Simulation Start Time, Integration Interval for WQ
C       Note: Units for -TZERO- are DAYS ... Irrespective of -TWARP-
      READ(IN,1000)   COMMENT
      READ(IN,2200,ERR=950)  TWARP,TZERO,IDTWQ
 2200 FORMAT(6X,A4,E10.4,I10)
      ISCALT=IUNITCHECK(TWARP,'TWARP   ')
      SCALT = ISCALT
      WRITE(OUT,2150)    TWARP 
 2150 FORMAT(10X,'TIME WARP FACTOR ',A4) 
      IF(TZERO.NE.0.)   WRITE(OUT,2175)   TZERO 
 2175 FORMAT(10X,'TIME =',F7.2,' DAYS IS THE START OF THE SIMULATION'/)
      IF(IDTWQ.NE.0.)   WRITE(OUT,2180)   IDTWQ,TWARP
 2180 FORMAT(
     .  10X,'INTEGRATION STEP-SIZE TO BE USED FOR THE KINETIC SUBROUTINE
     . -TUNER- IS',I7,1X,A4)
      IDTWQSECS=ISCALT*IDTWQ
      DTWQ=FLOAT(IDTWQSECS)/86400.
 
C        Note: Units for -TEND- and -TBRK- are DAYS ... Irrespective of -TWARP-
      IF(INTGRTYP.EQ.3 .OR. INTGRTYP.EQ.6)  THEN
         READ(IN,1000)   COMMENT
         READ(IN,2201,ERR=950)  IDTSPLIT,IDTFULL,TEND
 2201    FORMAT(2I10,F10.0)
         WRITE(OUT,2505)  IDTSPLIT,TWARP,IDTFULL,TWARP,TEND
 2505    FORMAT(10X,'SPLIT TIMESTEP INTEGRATION'/
     .     10X,'DTSPLIT =',I5,1X,A4,' DTFULL =',I5,1X,A4,
     .     ' TEND =',F6.1,' DAYS')
         IDTSPLITSECS=ISCALT*IDTSPLIT
         DTSPLIT = FLOAT(IDTSPLITSECS)/86400.
         IDTFULLSECS=ISCALT*IDTFULL
         DTFULL = FLOAT(IDTFULLSECS)/86400.
         IDTSECS = IDTFULLSECS
         DT = DTFULL
         IF(MOD(IDTFULL,IDTSPLIT).NE.0) THEN
           WRITE(OUT,2510)  IDTSPLIT,IDTFULL
 2510      FORMAT(///10X,'ERROR...IDTSPLIT =',I5,
     .       ' IS NOT AN EXACT MULTIPLE OF IDTFULL =',I5
     .       /10X,'RCA TERMINATED')
           CALL EXIT
         ENDIF
         IF(MOD(IDTWQ,IDTFULL).NE.0) THEN
           WRITE(OUT,2511)  IDTWQ,IDTFULL
 2511      FORMAT(///10X,'ERROR...IDTWQ =',I5,
     .       ' IS NOT AN EXACT MULTIPLE OF IDTFULL =',I5
     .       /10X,'RCA TERMINATED')
           CALL EXIT
         ENDIF
         NSTEP = 1
      ELSE
         READ(IN,1000)   COMMENT
         READ(IN,2210,ERR=950)  NSTEP
 2210    FORMAT(I10)
         READ(IN,1000)   COMMENT
         READ(IN,2202,ERR=950) (ISTEP(I),TBRK(I),I=1,NSTEP)         
 2202    FORMAT(4(I10,F10.0))
         WRITE(OUT,2350)  TWARP,TWARP,TWARP,TWARP
 2350    FORMAT(/5X,4('   DT(',A4,')    T(DAYS)')/)
         WRITE(OUT,2500) (ISTEP(I),TBRK(I),I=1,NSTEP)
 2500    FORMAT(5X,4(I11,F11.2))
         DO I=1,NSTEP
          ISTEP(I)=ISCALT*ISTEP(I)
         ENDDO
         TEND = TBRK(1)
         IDTSECS = ISTEP(1)
         DT = FLOAT(ISTEP(1))/86400.
        DO I=1,NSTEP
         IF(MOD(IDTWQ,ISTEP(I)).NE.0) THEN
           WRITE(OUT,2512)  IDTWQ,I,ISTEP(I)
 2512      FORMAT(///10X,'ERROR...IDTWQ =',I5,
     .       ' IS NOT AN EXACT MULTIPLE OF ISTEP(',I3,') =',I5
     .       /10X,'RCA TERMINATED')
           CALL EXIT
         ENDIF
      ENDDO
      IF(DT.EQ.0.0)   THEN
         WRITE(OUT,2310)
 2310    FORMAT(///10X,'ERROR...INTEGRATION STEP SIZE = 0.0')
         CALL EXIT
      ENDIF
      ENDIF

C     Global and Detailed Dump Intervals
      READ(IN,1000)   COMMENT
      READ(IN,1105,ERR=950) IPRNTG,IPRNTD,NDMPS,TWARPP,IGDOPT,IDDOPT
 1105 FORMAT(2I10,I10,6X,A4,2I10)
      ISCALP=IUNITCHECK(TWARPP,'TWARPP  ')
      IF(NDMPS.GT.0)   THEN
         OPEN(UNIT=12,FILE='RCAF12',FORM='UNFORMATTED')
         OPEN(UNIT=13,FILE=TRIM('RCAF13')//'_'//TRIM(ADJUSTL('1'))
     .      ,FORM='UNFORMATTED')
C        Segments for Detailed Dumps
         READ(IN,1000)   COMMENT
         READ(IN,1100)  ((IFDMPS(I,J),J=1,3),I=1,NDMPS)
 1100    FORMAT(18I4)
C        Bypass options for Detailed Dumps
         READ(IN,1000)   COMMENT
         READ(IN,1150)  (DDMPSBY(ISYS),ISYS=1,NOSYS)
 1150    FORMAT(40I2)
         DO 10 ISYS=1,NOSYS
           DDGDPBY(ISYS) = DDMPSBY(ISYS)
   10    CONTINUE
         WRITE(12)  NDMPS
         WRITE(12)  ((IFDMPS(I,J),J=1,3),I=1,NDMPS)
         WRITE(12)  DDGDPBY
      ENDIF
 
      IPRNTGSECS=ISCALP*IPRNTG
      IPRNTDSECS=ISCALP*IPRNTD
      PRNTG=FLOAT(IPRNTGSECS)/86400.
      PRNTD=FLOAT(IPRNTDSECS)/86400.
      WRITE(OUT,1200)   IPRNTGSECS,PRNTG
 1200 FORMAT(//10X,'GLOBAL PRINT INTERVAL',I7,' SECONDS (',
     .    E12.3,' DAYS)')
      IF(IGDOPT.EQ.0)  WRITE(OUT,1205)
 1205 FORMAT(10X,'USER SELECTED NO GLOBAL DUMP AVERAGING OPTION')
      IF(IGDOPT.EQ.1)  WRITE(OUT,1210)
 1210 FORMAT(10X,'USER SELECTED GLOBAL DUMP AVERAGING OPTION')
      IF(NDMPS.GT.0)  THEN
         WRITE(OUT,1220)   IPRNTDSECS,PRNTD
         WRITE(OUT,1235)   ((IFDMPS(I,J),J=1,3),I=1,NDMPS)
      ENDIF
 1220 FORMAT(/10X,'THE FOLLOWING SEGMENTS WILL BE SAVED USING A PRINT 
     .INTERVAL OF',I7,' SECONDS (',E12.3,' DAYS)')
 1235 FORMAT(10X,'(',2I4,I3,')',2X,'(',2I4,I3,')',2X,'(',2I4,I3,')',2X
     .          ,'(',2I4,I3,')',2X,'(',2I4,I3,')',2X,'(',2I4,I3,')',2X
     .          ,'(',2I4,I3,')',2X)
      IF(IDDOPT.EQ.0)  WRITE(OUT,1240)
 1240 FORMAT(10X,'USER SELECTED NO DETAILED DUMP AVERAGING OPTION')
      IF(IDDOPT.EQ.1)  WRITE(OUT,1245)
 1245 FORMAT(10X,'USER SELECTED DETAILED DUMP AVERAGING OPTION')

C     Segments to be displayed for intermediate dumps
      READ(IN,1000) COMMENT
      READ(IN,1100,ERR=950) ((IDUMP(I,J),J=1,3),I=1,6) 

C        Mass Balance Checks
      READ(IN,1000) COMMENT
      READ(IN,1300,ERR=950) MASSBAL,IPRNTMB,TWARPMB,IMBDOPT,
     .                      ISTARTMB,IENDMB     
 1300 FORMAT(2I10,6X,A4,3I10)
      IF(MASSBAL.EQ.0)  THEN
       WRITE(OUT,3000)
 3000  FORMAT(//10X,'MASS BALANCE CHECKS WILL NOT BE PERFORMED')
      ELSE
       ISCALPMB=IUNITCHECK(TWARPMB,'TWARPMB ')
       IPRNTMBSECS=ISCALPMB*IPRNTMB
       ISMBSECS=ISCALPMB*ISTARTMB
       IEMBSECS=ISCALPMB*IENDMB
       IF(IEMBSECS.EQ.0) THEN
         WRITE(OUT,3050) 
 3050    FORMAT(10X,'INPUT ERROR ... ZERO VALUE SPECIFIED FOR IENDMB'/
     .          10X,'RCA TERMINATED')
         CALL EXIT
       ENDIF
       IF(IPRNTMBSECS.EQ.0) THEN
         WRITE(OUT,3060) 
 3060    FORMAT(10X,'INPUT ERROR ... ZERO VALUE SPECIFIED FOR IPRNTMB'/
     .          10X,'RCA TERMINATED')
         CALL EXIT
       ENDIF
       PRNTMB=FLOAT(IPRNTMBSECS)/86400.
       OPEN(UNIT=17,FILE='RCAFMB',FORM='UNFORMATTED')
       WRITE(OUT,3100)  IPRNTMBSECS,PRNTMB
 3100  FORMAT(//
     .  10X,'USER HAS REQUESTED MASS BALANCE/FLUX BALANCE COMPUTATIONS'/
     .  10X,'BALANCES PRINT INTERVAL',I7,' SECONDS (',
     .       E12.3,' DAYS)')
       IF(IMBDOPT.EQ.0)  WRITE(OUT,3200)
 3200  FORMAT(10X,'BALANCES WILL BE INSTANTANEOUS VALUES')
       IF(IMBDOPT.EQ.1)  WRITE(OUT,3210)
 3210  FORMAT(10X,'BALANCES WILL BE AVERAGED VALUES')
      ENDIF

      RETURN

  950 CALL FMTER
      CALL EXIT 
      END 
