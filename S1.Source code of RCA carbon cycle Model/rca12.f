      SUBROUTINE RCA12
C
C        RCA12 PRINTS USER REQUESTED DUMPS AT END OF SIMULATION
C
      SAVE
      INCLUDE  'RCACM'
      CHARACTER   GDNAMES(50)*8,LYRDMP*3
      REAL   T(9000)
      INTEGER   IDMP(3,8)

      EQUIVALENCE
     .      (IDMP(1,1),IX1) , (IDMP(2,1),IY1) , (IDMP(3,1),IZ1)
     .    , (IDMP(1,2),IX2) , (IDMP(2,2),IY2) , (IDMP(3,2),IZ2)
     .    , (IDMP(1,3),IX3) , (IDMP(2,3),IY3) , (IDMP(3,3),IZ3)
     .    , (IDMP(1,4),IX4) , (IDMP(2,4),IY4) , (IDMP(3,4),IZ4)
     .    , (IDMP(1,5),IX5) , (IDMP(2,5),IY5) , (IDMP(3,5),IZ5)
     .    , (IDMP(1,6),IX6) , (IDMP(2,6),IY6) , (IDMP(3,6),IZ6)
     .    , (IDMP(1,7),IX7) , (IDMP(2,7),IY7) , (IDMP(3,7),IZ7)
     .    , (IDMP(1,8),IX8) , (IDMP(2,8),IY8) , (IDMP(3,8),IZ8)
 
C        READ TIME RECORD (SKIPPING -SYSBY- HEADER RECORD)
      REWIND 10
      READ(10)  INX,INY,INZ,INOSYS,NGDMP
      READ(10)  (GDNAMES(I),I=1,NGDMP)
      READ(10)
      READ(10)
      DO 30 I=1,9000
   30 READ(10,END=40)  T(I)
   40 IF(I.GT.40)  WRITE(OUT,1100)
 1100 FORMAT(/////10X,'DUE TO DIMENSION LIMITS, THE FOLLOWING DUMPS WILL
     . BE LIMITED TO THE FIRST 40 TIME RECORDS'//)
      IREC=MIN(40,I-1)
      IF(IREC.LE.16)   WRITE(OUT,1110)
 1110 FORMAT('1') 
      LASTVAR=0

      READ(IN,1234) COMMENT
 1234 FORMAT(A)

  100   READ(IN,1000,ERR=950,END=975)    IVAR,((IDMP(J,I),J=1,3),
     .     I=1,8),LYRDMP,LYR 
 1000   FORMAT(25I3,A3,I3)
        IF(IVAR.LE.0) GO TO 975
        IF(IVAR.GT.NGDMP) GO TO 970
        IF(IVAR.NE.LASTVAR) THEN
C         PRODUCE A SYSTEM ONLY DUMP FILE FROM MASTER DUMP FILE 
          CALL RCABYSS (NGDMP,IVAR)
          LASTVAR=IVAR
        ENDIF

        IF(LYRDMP.EQ.'LYR')   THEN
         WRITE(OUT,1200)  IVAR,GDNAMES(IVAR),LYR
 1200    FORMAT('1'//20X,'CONCS FOR SYSTEM',I3,' == ',A8,'  LAYER ',I3/)
         DO 110 IT=1,IREC
         IF( (T(2)-T(1)).LT.1.0 )  WRITE(OUT,3001)   T(IT)
 3001    FORMAT(//5X,'TIME =',F6.3)
         IF( (T(2)-T(1)).GE.1.0 )  WRITE(OUT,3002)   T(IT)
 3002    FORMAT(//5X,'TIME =',F6.1)
         CALL RCAPRNT(SCRATCH_ARRAY(1,1,1,IT),LYR,LYR)
  110    CONTINUE

        ELSE
  
         IF(IX1.EQ.0)   GO TO 400
         IF(IREC.GT.16)  WRITE(OUT,2200)  IVAR,GDNAMES(IVAR)
 2200    FORMAT(/////20X,'CONCS FOR SYSTEM',I3,' == ',A8/)
         IF(IREC.LE.16)  WRITE(OUT,2201)  IVAR,GDNAMES(IVAR)
 2201    FORMAT(/////20X,'CONCS FOR SYSTEM',I3,' == ',A8/)
         WRITE(OUT,2500)  ((IDMP(J,I),J=1,3),I=1,8)
 2500    FORMAT(1X,//4X,'TIME',1X,8(1X,3I4)) 

         DO 200 I=1,IREC 
          IF( (T(2)-T(1)).LT.1.0 )  WRITE(OUT,3000) T(I),
     .     SCRATCH_ARRAY(IX1,IY1,IZ1,I),SCRATCH_ARRAY(IX2,IY2,IZ2,I),
     .     SCRATCH_ARRAY(IX3,IY3,IZ3,I),SCRATCH_ARRAY(IX4,IY4,IZ4,I),
     .     SCRATCH_ARRAY(IX5,IY5,IZ5,I),SCRATCH_ARRAY(IX6,IY6,IZ6,I),
     .     SCRATCH_ARRAY(IX7,IY7,IZ7,I),SCRATCH_ARRAY(IX8,IY8,IZ8,I) 
 3000    FORMAT(1X,F8.5,8E13.4)
          IF( (T(2)-T(1)).GE.1.0 )  WRITE(OUT,3010) T(I),
     .      SCRATCH_ARRAY(IX1,IY1,IZ1,I),SCRATCH_ARRAY(IX2,IY2,IZ2,I),
     .      SCRATCH_ARRAY(IX3,IY3,IZ3,I),SCRATCH_ARRAY(IX4,IY4,IZ4,I),
     .      SCRATCH_ARRAY(IX5,IY5,IZ5,I),SCRATCH_ARRAY(IX6,IY6,IZ6,I),
     .      SCRATCH_ARRAY(IX7,IY7,IZ7,I),SCRATCH_ARRAY(IX8,IY8,IZ8,I) 
 3010    FORMAT(1X,F8.1,8E13.4)
  200    CONTINUE

        ENDIF

      GO TO 100

  400 RETURN

  950 CALL FMTER
      CALL EXIT
  970 WRITE(OUT,4500)  IVAR
 4500 FORMAT(//10X,'ILLEGAL VARIABLE NUMBER SPECIFIED FOR GLOBAL DUMPS'/
     .         10X,'IVAR = ',I3,5X,'RCA TERMINATED')

  975 RETURN
      END