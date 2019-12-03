cOPTION CROSSREF ON
cOPTION SYMTABLE ON
cPAGEWIDTH 132
      SUBROUTINE SED_READ
      SAVE
C
C*******************************************************************************
C
C             VARIABLE STOICHIOMETRY SEDIMENT SUBMODEL - VSSEDMODL
C             ----------------------------------------------------
C
C                      ORIGINALLY DEVELOPED FOR THE CHESAPEAKE BAY
C
C   ORGANIC FLUX TO SEDIMENT ==> DIAGENESIS ==> INORGANIC FLUX TO WATER COLUMN
C             3-G MODEL - G1=LABILE, G2=REFRACTORY, G3=CONSERVATIVE
C
C*******************************************************************************
C
C  INPUTS
C  ------
C             REQUIRED INPUTS FOR SEDIMENT SUB-MODEL
C
C  A.  PASSED TO SEDIMENT SUBROUTINE FROM WATER QUALITY SUBROUTINE
C
C    1.  OVERLYING WATER COLUMN CONCENTRATIONS OF AMMONIA, NITRATE,
C        SALINITY, DISSOLVED OXYGEN, PHOSPHORUS AND SILICA
C    2.  OVERLYING WATER COLUMN SEGMENT VOLUME (V1)
C    3.  OVERLYING WATER COLUMN SEGMENT DEPTH (BL(ISEG,3))
C    4.  OVERLYING WATER COLUMN SEGMENT TEMPERATURE
C    8.  ALGAL PHOSPHORUS TO CARBON, NITROGEN TO CARBON, AND SILICA
C        TO CARBON RATIOS
C    9.  DEPTH OF OVERLYING WATER COLUMN
C
C  B.  INPUTS SUPPLIED VIA DIRECT INPUT TO THE SEDIMENT SUBROUTINE
C
C     VARIABLE NAMES        DESCRIPTION                            UNITS
C     --------------   ----------------------------------------    -----
C           HSED       DEPTH OF SEDIMENT LAYER                        CM
C          DIFFT       WATER COLUMN-SEDIMENT LAYER DIFFUSION     CM2/SEC
C                      COEFFICIENT FOR TEMPERATURE
C         SALTSW       SALINITY SWITCH: (ARCHAIC: FORMERLY USED TO   PPT
C                      DETERMINE WHETHER METHANE OR SULFIDE SOD
C                      FORMULATION IS TO BE USED); NOW AFFECTS
C                      BENTHIC NITRIF AND DE-NITRIF RATES
C         XXSNET       NET SETTLING VELOCITY FROM WATER COLUMN
C
C  POM FLUX RELATED
C                      FRACTIONS OF G1, G2, AND G3 CONTAINED IN ...
C        FRPPH1(3)     ALGAL GROUP NO.1 PHOSPHORUS
C        FRPPH2(3)     ALGAL GROUP NO.2 PHOSPHORUS
C        FRPPH3(3)     ALGAL GROUP NO.3 PHOSPHORUS
C        FRPOP(3)      NON-ALGAL PARTICULATE ORGANIC PHOSPHORUS
C        FRNPH1(3)     ALGAL GROUP NO.1 NITROGEN
C        FRNPH2(3)     ALGAL GROUP NO.2 NITROGEN
C        FRNPH3(3)     ALGAL GROUP NO.3 NITROGEN
C        FRPON(3)      NON-ALGAL PARTICULATE ORGANIC NITROGEN
C        FRCPH1(3)     ALGAL GROUP NO.1 CARBON
C        FRCPH2(3)     ALGAL GROUP NO.2 CARBON
C        FRCPH3(3)     ALGAL GROUP NO.3 CARBON
C        FRPOC(3)      NON-ALGAL PARTICULATE ORGANIC CARBON
C
C  DIAGENESIS RELATED
C        KPDIAG(3)     REACTION RATES FOR POP G1, G2, AND G3        /DAY
C        DPTHTA(3)     TEMPERATURE THETAS FOR POP G1, G2, AND G3
C        KNDIAG(3)     REACTION RATES FOR PON G1, G2, AND G3        /DAY
C        DNTHTA(3)     TEMPERATURE THETAS FOR PON G1, G2, AND G3
C        KCDIAG(3)     REACTION RATES FOR POC G1, G2, AND G3        /DAY
C        DCTHTA(3)     TEMPERATURE THETAS FOR POC G1, G2, AND G3
C                      
C
C  FLUX-FLUX RELATED
C        VPMIX(NSED)   SOLID-PHASE MIXING RATE                  M**2/DAY
C        VDMIX(NSED)   DISSOLVED-PHASE MIXING RATE              M**2/DAY
C           M1         CONCENTRATION OF SS IN LAYER 1            KG SS/L
C           M2         CONCENTRATION OF SS IN LAYER 2            KG SS/L
C          PIE1S       PARTITION COEFFICIENT ASSOCIATED             L/KG
C                      SULFIDE SORPTION IN LAYER 1
C          PIE2S       PARTITION COEFFICIENT ASSOCIATED             L/KG
C                      SULFIDE SORPTION IN LAYER 2
C          KLN2        NITRATE TRANSFER COEFFICIENT BETWEEN        M/DAY
C                        SEDIMENT LAYERS 1 AND 2
C        REACTION VELOCITIES FOR...
C          KAPPD1      DISSOLVED SULFIDE IN LAYER 1                M/DAY
C          KAPPP1      PARTICULATE SULFIDE IN LAYER 1              M/DAY
C          KAPPNH4     AMMONIA NITROGEN                            M/DAY
C          KP1NO3      AEROBIC DENITRIFICATION                     M/DAY
C          KP2NO3      ANAEROBIC DENITRIFICATION                   M/DAY
C          KAPPC       METHANE                                     M/DAY
C          KAPPD       METHANE DIFFUSION MASS TRANSFER             M/DAY
C
C
C  OUTPUT
C  ------
C     THE SUBROUTINE RETURNS FLUXES FOR....
C      JSOD, JAQSOD, JCH4AQ AND JCH4G [GM O2*/M2-DAY]
C      JNH4, JPO4, JNO3 AND JSI [MG/M2-DAY]
C
C***********************************************************************
C
      INCLUDE  'RCACM'
      INCLUDE  'SEDCM'
      INCLUDE  'SEDLB'

      CHARACTER   SSNAME(25)*25,FRNAME(9)*24,SEDFILNA*40
     .           ,COMMENT*80,TWARPSED*4,DUMPNAMES(76)*8
      INTEGER  ISEDPRNT,PIESIOPT,PIEPO4OPT
      REAL  PGFRAC(3),NGFRAC(3),CGFRAC(3),H1S(NX,NY)

      EXTERNAL  SEDF

C        WATER COLUMN STATE-VARIABLES
      REAL
     .      SAL(NX,NY,NZ)   , PHYT1(NX,NY,NZ)
     .    , PHYT2(NX,NY,NZ) , PHYT3(NX,NY,NZ)
     .    , RPOP(NX,NY,NZ)  , LPOP(NX,NY,NZ)
     .    , RDOP(NX,NY,NZ)  , LDOP(NX,NY,NZ)
     .    , PO4T(NX,NY,NZ)
     .    , RPON(NX,NY,NZ)  , LPON(NX,NY,NZ)
     .    , RDON(NX,NY,NZ)  , LDON(NX,NY,NZ)
     .    , NH4T(NX,NY,NZ)  , NO23(NX,NY,NZ)
     .    , BSI(NX,NY,NZ)   , SIT(NX,NY,NZ)
     .    , RPOC(NX,NY,NZ)  , LPOC(NX,NY,NZ)
     .    , RDOC(NX,NY,NZ)  , LDOC(NX,NY,NZ)
     .    , EXDOC(NX,NY,NZ) , REPOC(NX,NY,NZ)
     .    , REDOC(NX,NY,NZ)
     .    , O2EQ(NX,NY,NZ)  , DO(NX,NY,NZ)
     .    , DIC(NX,NY,NZ)  , TALK(NX,NY,NZ) 
 
C        GLOBAL DUMP AVERAGE VARIABLES
      REAL
     .      CTEMP_GDA(NX,NY)
     .    , POP1R_GDA(NX,NY)     , POP2R_GDA(NX,NY)
     .    , POP3R_GDA(NX,NY)     , POPR_GDA(NX,NY)
     .    , PON1R_GDA(NX,NY)     , PON2R_GDA(NX,NY)
     .    , PON3R_GDA(NX,NY)     , PONR_GDA(NX,NY)
     .    , POC1R_GDA(NX,NY)     , POC2R_GDA(NX,NY)
     .    , POC3R_GDA(NX,NY)     , POCR_GDA(NX,NY)
     .    , PO4T2R_GDA(NX,NY)    , HST2R_GDA(NX,NY)
     .    , SIT2R_GDA(NX,NY)     , PSIAVR_GDA(NX,NY)
     .    , FLXPOP_GDA(NX,NY)    , FLXPON_GDA(NX,NY)
     .    , FLXPOC_GDA(NX,NY)    , FLXPOS_GDA(NX,NY)
     .    , FLXPO4T2S_GDA(NX,NY) , FLXSIT2S_GDA(NX,NY)
      REAL
     .      FLXPOP1_GDA(NX,NY) , FLXPOP2_GDA(NX,NY) , FLXPOP3_GDA(NX,NY)
     .    , FLXPON1_GDA(NX,NY) , FLXPON2_GDA(NX,NY) , FLXPON3_GDA(NX,NY)
     .    , FLXPOC1_GDA(NX,NY) , FLXPOC2_GDA(NX,NY) , FLXPOC3_GDA(NX,NY)
      REAL
     .      O20_GDA(NX,NY)     , CSOD_GDA(NX,NY)
     .    , SOD_GDA(NX,NY)     , S_GDA(NX,NY)
     .    , JP_GDA(NX,NY)      , JN_GDA(NX,NY)
     .    , JC_GDA(NX,NY)      , JO2NH4_GDA(NX,NY)
     .    , XJCO2AV_GDA(NX,NY) , XJC1AV_GDA(NX,NY)
      REAL
     .      JPO4_GDA(NX,NY)    , JNH4_GDA(NX,NY)
     .    , JNO3_GDA(NX,NY)    , JH2S_GDA(NX,NY)
     .    , JSI_GDA(NX,NY)     , JCH4AQ_GDA(NX,NY)
     .    , JCH4G_GDA(NX,NY)   , JN2_GDA(NX,NY)
     .    , H1_GDA(NX,NY)
      REAL
     .      PO40_GDA(NX,NY)    , PO41_GDA(NX,NY)
     .    , PO42_GDA(NX,NY)    , PO4T2_GDA(NX,NY)
     .    , NH40_GDA(NX,NY)    , NH41_GDA(NX,NY)
     .    , NH42_GDA(NX,NY)    , NH4T2_GDA(NX,NY)
     .    , NO30_GDA(NX,NY)    , NO31_GDA(NX,NY)
     .    , NO32_GDA(NX,NY)    , NO3T2_GDA(NX,NY)
     .    , HS1_GDA(NX,NY)     , HS2_GDA(NX,NY)
     .    , HST2_GDA(NX,NY)
     .    , SI0_GDA(NX,NY)     , SI1_GDA(NX,NY)
     .    , SI2_GDA(NX,NY)     , SIT2_GDA(NX,NY)
      REAL
     .      JN2(NX,NY)
     .    , CH41_GDA(NX,NY)    , CH42_GDA(NX,NY)
     .    , CH4T2_GDA(NX,NY)   , SO41_GDA(NX,NY)
     .    , SO42_GDA(NX,NY)    , SO4T2_GDA(NX,NY)

 
C        REAL AND LABELED COMMON FOR LINKING WATER COLUMN AND SEDIMENT SUBMODELS
C        LABELED COMMON FOR EUTROPHICATION ALGAL GROWTH ROUTINES
      REAL
     .    PO4(NX,NY,NZ),NH4(NX,NY,NZ),SI(NX,NY,NZ)
      EQUIVALENCE
     .    (CKINARRAY(1,1,1,1),PO4(1,1,1))
     .   ,(CKINARRAY(1,1,1,2),NH4(1,1,1))
     .   ,(CKINARRAY(1,1,1,3),SI(1,1,1))
      REAL
     .   PCRB1(NX,NY,NZ),PCRB2(NX,NY,NZ)
     .  ,PCRB3(NX,NY,NZ)
     .  ,NCRB1(NX,NY,NZ),NCRB2(NX,NY,NZ)
     .  ,NCRB3(NX,NY,NZ)
     .  ,SCRB1(NX,NY,NZ),SCRB2(NX,NY,NZ)
     .  ,SCRB3(NX,NY,NZ)
      COMMON /EUTRO/ 
     .   PCRB1,PCRB2,PCRB3,NCRB1,NCRB2,NCRB3,SCRB1,SCRB2,SCRB3
      REAL
     .   SODS(NX,NY),FLXHS(NX,NY),FLXNH4(NX,NY)
     .  ,FLXPO4(NX,NY),FLXNO3(NX,NY),FLXSI(NX,NY)
     .  ,FLXCH4AQ(NX,NY),FLXCH4G(NX,NY)
     .  ,FLXDIC(NX,NY),FLXTALK(NX,NY)
     .  ,FLXPOCN(NX,NY) 
      COMMON /FLUXES/
     .   DEPFLUX(NX,NY,15),CFLUXS(NX,NY,11)
      EQUIVALENCE
     .   (CFLUXS(1,1,1),SODS(1,1))     , (CFLUXS(1,1,2),FLXHS(1,1))  
     .  ,(CFLUXS(1,1,3),FLXCH4AQ(1,1)) , (CFLUXS(1,1,4),FLXCH4G(1,1))
     .  ,(CFLUXS(1,1,5),FLXNH4(1,1))   , (CFLUXS(1,1,6),FLXNO3(1,1)) 
     .  ,(CFLUXS(1,1,7),FLXPO4(1,1))   , (CFLUXS(1,1,8),FLXSI(1,1))
     .  ,(CFLUXS(1,1,9),FLXDIC(1,1))   , (CFLUXS(1,1,10),FLXTALK(1,1))
     .  ,(CFLUXS(1,1,11),FLXPOCN(1,1))

      EQUIVALENCE
     .   (CARAY(1,1,1,1),SAL(1,1,1))    , (CARAY(1,1,1,2),PHYT1(1,1,1))
     . , (CARAY(1,1,1,3),PHYT2(1,1,1))  , (CARAY(1,1,1,4),PHYT3(1,1,1)) 
     . , (CARAY(1,1,1,5),RPOP(1,1,1))   , (CARAY(1,1,1,6),LPOP(1,1,1)) 
     . , (CARAY(1,1,1,7),RDOP(1,1,1))   , (CARAY(1,1,1,8),LDOP(1,1,1))
     . , (CARAY(1,1,1,9),PO4T(1,1,1))
     . , (CARAY(1,1,1,10),RPON(1,1,1))  , (CARAY(1,1,1,11),LPON(1,1,1))
     . , (CARAY(1,1,1,12),RDON(1,1,1))  , (CARAY(1,1,1,13),LDON(1,1,1))
     . , (CARAY(1,1,1,14),NH4T(1,1,1))  , (CARAY(1,1,1,15),NO23(1,1,1))
     . , (CARAY(1,1,1,16),BSI(1,1,1))   , (CARAY(1,1,1,17),SIT(1,1,1))
     . , (CARAY(1,1,1,18),RPOC(1,1,1))  , (CARAY(1,1,1,19),LPOC(1,1,1))
     . , (CARAY(1,1,1,20),RDOC(1,1,1))  , (CARAY(1,1,1,21),LDOC(1,1,1))
     . , (CARAY(1,1,1,22),EXDOC(1,1,1)) , (CARAY(1,1,1,23),REPOC(1,1,1))
     . , (CARAY(1,1,1,24),REDOC(1,1,1)) , (CARAY(1,1,1,25),O2EQ(1,1,1))
     . , (CARAY(1,1,1,26),DO(1,1,1))
     . , (CARAY(1,1,1,27),DIC(1,1,1))  
     . , (CARAY(1,1,1,28),TALK(1,1,1)) 

      DATA  SSNAME/
     .  '    Temperature [Deg C]    ' , '     G1 POP [mg P/m**3]    ' ,
     .  '     G1 PON [mg N/m**3]    ' , '     G1 POC [mg C/m**3]    ' ,
     .  '     G2 POP [mg P/m**3]    ' , '     G2 PON [mg N/m**3]    ' ,
     .  '     G2 POC [mg C/m**3]    ' , '     G3 POP [mg P/m**3]    ' ,
     .  '     G3 PON [mg N/m**3]    ' , '     G3 POC [mg C/m**3]    ' ,
     .  '  Biogenic Si [mg Si/m**3] ' , ' Dissolved PO4 [mg P/m**3] ' ,
     .  ' Dissolved NH4 [mg N/m**3] ' , ' Dissolved NO3 [mg N/m**3] ' ,
     .  'Dissolved H2S [mg O2*/m**3]' , ' Dissolved Si [mg Si/m**3] ' ,
     .  '      Benthic Stress       ' , ' Dissolved CH4 [mg O2/m**3]' ,
     .  '    Sulfate [mg O2/m**3]   ' , '      DIC [ umol C/m**3]   ',
     .  '      TALK [umol C/m**3]    ' , '                           ' ,
     .  '                           ' , '                           ' ,
     .  '                           '/

      DATA   FRNAME/
     .    'WINTER PHYTO PHOSPHORUS ','SUMMER PHYTO PHOSPHORUS '
     .   ,'DETRITAL ORG PHOSPHORUS '
     .   ,'WINTER PHYTO NITROGEN   ','SUMMER PHYTO NITROGEN   '
     .   ,'DETRITAL ORG NITROGEN   '
     .   ,'WINTER PHYTO CARBON     ','SUMMER PHYTO CARBON     '
     .   ,'DETRITAL ORG CARBON     '/
      DATA  DUMPNAMES  /
     .   ' sedtemp' , '  G1 POP' , '  G2 POP' , '  G3 POP' , 'Totl POP' 
     . , '  G1 PON' , '  G2 PON' , '  G3 PON' , 'Totl PON' , '  G1 POC' 
     . , '  G2 POC' , '  G3 POC' , 'Totl POC' , 'Ly2 TPO4' , 'Ly2 TH2S' 
     . , 'Lyr2 TSi' , 'BiogenSi' , '    JPOP' , '    JPON' , '    JPOC' 
     . , '     O20' , '    CSOD' , '     SOD' , 's-MasTrf' , 'JPdiagen' 
     . , 'JNdiagen' , 'JCdiagen' , '  JO2NH4' , ' XJCO2AV' , '  XJC1AV' 
     . , '    JPO4' , '    JNH4' , '    JNO3' , '     JHS' , '     JSI' 
     . , '  JCH4AQ' , '   JCH4G' , '      H1' , '    PO40' , '    PO41' 
     . , '    PO42' , '   PO4T2' , '    NH40' , '    NH41' , '    NH42' 
     . , '   NH4T2' , '    NO30' , '    NO31' , '    NO32' , '   NO3T2' 
     . , '     HS1' , '     HS2' , '    HST2' , '     SI0' , '     SI1' 
     . , '     SI2' , '    SIT2' , '     JN2' , '    CH41' , '    CH42' 
     . , '   CH4T2' , '    SO41' , '    SO42' , '   SO4T2' , '  JPOPG1' 
     . , '  JPOPG2' , '  JPOPG3' , '  JPONG1' , '  JPONG2' , '  JPONG3' 
     . , '  JPOCG1' , '  JPOCG2' , '  JPOCG3' , '   JPBSi' , ' JdepPO4' 
     . , '  JdepSi' /
C             INITIALIZATION

C        OPEN FILES FOR SEDIMENT DUMPS

      !OPEN(14,FILE='RCAF14',FORM='UNFORMATTED')
      OPEN(16,FILE='RCAFICSED',FORM='UNFORMATTED')
      IF(MASSBAL.EQ.1)  OPEN(18,FILE='RCAFSEDBAL',FORM='UNFORMATTED')

C        WRITE -FSM- (LAND-MASK) TO RCAF14
      WRITE(14)  FSM,76
c     do 1069 iy = 1,NY
c     WRITE(14,1068)  (FSM(ix,iy),ix=1,NX)
c1068 format(53f4.0)
c1069 continue     
C        SET NUMBER OF LAYERS FOR SEDIMENT MODEL => 2
C        (AEROBIC + 1 ANAEROBIC)
      
      nsedlyr = 2
      xnsedlyr = nsedlyr

      OPEN(40,FILE=KINFILNA(1),FORM='FORMATTED')

C        PRINT CONTROLS
      READ(40,1080)   COMMENT
 1080 FORMAT(A80)
      READ(40,1100)   ISEDPRNT,IPRNTSED,TWARPSED,IGDSEDOPT
      ISEDPRNT = 0.0
 1100 FORMAT(2I10,6X,A4,I10)
      IF(TWARPSED.NE.'SECS' .AND. TWARPSED.NE.'secs' .AND.
     .   TWARPSED.NE.'MINS' .AND. TWARPSED.NE.'mins' .AND.
     .   TWARPSED.NE.'HRS ' .AND. TWARPSED.NE.'hrs ' .AND.
     .   TWARPSED.NE.' HRS' .AND. TWARPSED.NE.' hrs' .AND.
     .   TWARPSED.NE.'DAYS' .AND. TWARPSED.NE.'days')  THEN
       WRITE(OUT,9001)  TWARPSED
 9001  FORMAT(/5X,'THE UNITS CHOSEN FOR TWARPSED ',A4,' ARE NOT VALID'/
     .   5X,'RCA TERMINATED')
       CALL EXIT
      ENDIF
      IF(TWARPSED.EQ.'SECS'.OR.TWARPSED.EQ.'secs')  ISCALPSED=1
      IF(TWARPSED.EQ.'MINS'.OR.TWARPSED.EQ.'mins')  ISCALPSED=60
      IF(TWARPSED.EQ.'HRS '.OR.TWARPSED.EQ.'hrs ')  ISCALPSED=3600
      IF(TWARPSED.EQ.' HRS'.OR.TWARPSED.EQ.' hrs')  ISCALPSED=3600
      IF(TWARPSED.EQ.'DAYS'.OR.TWARPSED.EQ.'days')  ISCALPSED=86400
      IPRNTSEDSECS=ISCALPSED*IPRNTSED
      PRNTSED=IPRNTSEDSECS/86400.
      WRITE(OUT,1110)   IPRNTSEDSECS,PRNTSED
 1110 FORMAT(//5X,'SEDIMENT DUMP PRINT INTERVAL',I7,' SECONDS (',
     .    E12.3,' DAYS)')
      IF(IGDSEDOPT.EQ.0)  WRITE(OUT,1120)
 1120 FORMAT(5X,'USER SELECTED NO GLOBAL DUMP AVERAGING OPTION')
      IF(IGDSEDOPT.EQ.1)  WRITE(OUT,1130)
 1130 FORMAT(5X,'USER SELECTED GLOBAL DUMP AVERAGING OPTION')
      IGDCNT=0
      NXPRTSED = IPRNTSEDSECS
      IF(TZERO.NE.0.)  NXPRTSED = ITIMESECS+IPRNTSEDSECS
C        Initialize variables for Global Dump Averaging, if required
      IF(IGDSEDOPT.EQ.1) THEN
        DO IY = 1,NY
          DO IX = 1,NX
            CTEMP_GDA(IX,IY) = 0.
            POP1R_GDA(IX,IY) = 0.
            POP2R_GDA(IX,IY) = 0.
            POP3R_GDA(IX,IY) = 0.
            POPR_GDA(IX,IY) = 0.
            PON1R_GDA(IX,IY) = 0.
            PON2R_GDA(IX,IY) = 0.
            PON3R_GDA(IX,IY) = 0.
            PONR_GDA(IX,IY) = 0.
            POC1R_GDA(IX,IY) = 0.
            POC2R_GDA(IX,IY) = 0.
            POC3R_GDA(IX,IY) = 0.
            POCR_GDA(IX,IY) = 0.
            PO4T2R_GDA(IX,IY) = 0.
            HST2R_GDA(IX,IY) = 0.
            SIT2R_GDA(IX,IY) = 0.
            PSIAVR_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            FLXPOP_GDA(IX,IY) = 0.
            FLXPOP1_GDA(IX,IY) = 0.
            FLXPOP2_GDA(IX,IY) = 0.
            FLXPOP3_GDA(IX,IY) = 0.
            FLXPON_GDA(IX,IY) = 0.
            FLXPON1_GDA(IX,IY) = 0.
            FLXPON2_GDA(IX,IY) = 0.
            FLXPON3_GDA(IX,IY) = 0.
            FLXPOC_GDA(IX,IY) = 0.
            FLXPOC1_GDA(IX,IY) = 0.
            FLXPOC2_GDA(IX,IY) = 0.
            FLXPOC3_GDA(IX,IY) = 0.
            FLXPOS_GDA(IX,IY) = 0.
            FLXPO4T2S(IX,IY) = 0.
            FLXSIT2S(IX,IY) = 0.
            O20_GDA(IX,IY) = 0.
            CSOD_GDA(IX,IY) = 0.
            SOD_GDA(IX,IY) = 0.
            S_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            JP_GDA(IX,IY) = 0.
            JN_GDA(IX,IY) = 0.
            JC_GDA(IX,IY) = 0.
            JO2NH4_GDA(IX,IY) = 0.
            XJCO2AV_GDA(IX,IY) = 0.
            XJC1AV_GDA(IX,IY) = 0.
            JPO4_GDA(IX,IY) = 0.
            JNH4_GDA(IX,IY) = 0.
            JNO3_GDA(IX,IY) = 0.
            JH2S_GDA(IX,IY) = 0.
            JSI_GDA(IX,IY) = 0.
            JCH4AQ_GDA(IX,IY) = 0.
            JCH4G_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            H1_GDA(IX,IY) = 0.
            PO40_GDA(IX,IY) = 0.
            PO41_GDA(IX,IY) = 0.
            PO42_GDA(IX,IY) = 0.
            PO4T2_GDA(IX,IY) = 0.
            NH40_GDA(IX,IY) = 0.
            NH41_GDA(IX,IY) = 0.
            NH42_GDA(IX,IY) = 0.
            NH4T2_GDA(IX,IY) = 0.
            NO30_GDA(IX,IY) = 0.
            NO31_GDA(IX,IY) = 0.
            NO32_GDA(IX,IY) = 0.
            NO3T2_GDA(IX,IY) = 0.
            HS1_GDA(IX,IY) = 0.
            HS2_GDA(IX,IY) = 0.
            HST2_GDA(IX,IY) = 0.
            SI0_GDA(IX,IY) = 0.
            SI1_GDA(IX,IY) = 0.
            SI2_GDA(IX,IY) = 0.
            SIT2_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            JN2_GDA(IX,IY) = 0.
            CH41_GDA(IX,IY) = 0.
            CH42_GDA(IX,IY) = 0.
            CH4T2_GDA(IX,IY) = 0.
            SO41_GDA(IX,IY) = 0.
            SO42_GDA(IX,IY) = 0.
            SO4T2_GDA(IX,IY) = 0.
          ENDDO
        ENDDO
      ENDIF

C             GET MODEL COEFFICIENTS REQUIRED FOR SEDIMENT MODEL

C        SEDIMENT DEPTH
      READ(40,1080)   COMMENT
      DO 10 IY=1,NY
       READ(40,1200,ERR=900) (HSED(IX,IY),IX=1,NX)
 1200 FORMAT(8F10.0)
   10 CONTINUE

C        COMPUTE SEDIMENT VOLUMES ==> M**3

      DO 20 IY=1,NY
       DO 20 IX=1,NX

C        SEDIMENT VOLUME (M**3)
        BSVOL(IX,IY) = 0.01*HSED(IX,IY)*XAZ(IX,IY)

C        SET BFORMAX
        BFORMAXS(IX,IY) = 0.0
        ISWBNTHS(IX,IY) = 0

   20 CONTINUE

      WRITE(OUT,1010)
 1010 FORMAT(///34X,
     .'SEDIMENT-WATER COLUMN LINKAGES AND SEDIMENT DEPTHS AND VOLUMES'/)
       WRITE(OUT,1020)
 1020    FORMAT(10X,'SEDIMENT SEGMENT DEPTHS (CM)')
C       IF(ISEDPRNT.EQ.1) CALL RCAPRNT(HSED,1,1)
C       WRITE(OUT,1021)
 1021    FORMAT(10X,'SEDIMENT VOLUMES (M**3)')
C       IF(ISEDPRNT.EQ.1) CALL RCAPRNT(BSVOL,1,1)

C        CONVERT DEPTHS TO METERS

      DO 30 IY=1,NY
       DO 30 IX=1,NX
        HSED(IX,IY) = 0.01*HSED(IX,IY)
   30 CONTINUE
 
C        READ SEDIMENT INITIAL CONDITIONS
      WRITE(OUT,1040)
 1040 FORMAT(////38X,
     .  'S E D I M E N T   I N I T I A L   C O N D I T I O N S'/)

C        INITIAL CONDITIONS

      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CTEMP,SSNAME(1))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOP(1,1,1),SSNAME(2))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPON(1,1,1),SSNAME(3))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOC(1,1,1),SSNAME(4))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOP(1,1,2),SSNAME(5))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPON(1,1,2),SSNAME(6))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOC(1,1,2),SSNAME(7))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOP(1,1,3),SSNAME(8))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPON(1,1,3),SSNAME(9))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOC(1,1,3),SSNAME(10))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CPOS,SSNAME(11))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,PO4T2TM1S,SSNAME(12))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,NH4T2TM1S,SSNAME(13))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,NO3T2TM1S,SSNAME(14))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,HST2TM1S,SSNAME(15))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,SIT2TM1S,SSNAME(16))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,BNTHSTR1S,SSNAME(17))

      CALL ICREAD(NX,NY,OUT,ISEDPRNT,PO4T1TM1S,SSNAME(12))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,NH4T1TM1S,SSNAME(13))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,NO3T1TM1S,SSNAME(14))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,SIT1TM1S,SSNAME(16))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,HST1TM1S,SSNAME(15))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CH4T1TM1S,SSNAME(18))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,CH4T2TM1S,SSNAME(18))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,SO4T2TM1S,SSNAME(19))

      CALL ICREAD(NX,NY,OUT,ISEDPRNT,DICT1TM1S,SSNAME(20))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,DICT2TM1S,SSNAME(20))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,TAT1TM1S,SSNAME(21))
      CALL ICREAD(NX,NY,OUT,ISEDPRNT,TAT2TM1S,SSNAME(21))
c by lb 20171019 DICT1TM1S,DICT2TM1S,TAT1TM1S,TAT2TM1S
      do iy=1,ny
       do ix=1,nx
        ch41tm1s(ix,iy)=ch4t1tm1s(ix,iy)
       enddo
      enddo

      if(cycle.eq.1)    then
         rewind(16)
         READ(16)   CTEMP,CPOP,CPON,CPOC,CPOS
     .           ,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S
     .           ,HST2TM1S,SIT2TM1S,BNTHSTR1S
     .           ,PO4T1TM1S,NH4T1TM1S,NO3T1TM1S
     .           ,SIT1TM1S,HST1TM1S
     .           ,CH4T1TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S
     .           ,PO41TM1S,NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S
     .           ,O20TM1S,SODTM1S,DD0TM1S
     .           ,BFORMAXS,ISWBNTHS
     .           ,DICT1TM1S,DICT2TM1S,TAT1TM1S,TAT2TM1S
     .           ,DIC1TM1S,TA1TM1S
c by lb 20171019 DICT1TM1S,DICT2TM1S,TAT1TM1S,TAT2TM1S
c         dic0 = DICT1TM1S(87,159)
c         dic1 = DIC1TM1S(87,159)
c         ta0 = TAT1TM1S(87,159)
c         ta1 = TA1TM1S(87,159)
c       do iy=1,ny
c       do ix=1,nx
c        DICT1TM1S(ix,iy) = dic0/2.
c        TAT1TM1S(ix,iy) = ta0/2.
c        DIC1TM1S(ix,iy) = dic1/2.
c        TA1TM1S(ix,iy) = ta1/2.
c       enddo
c      enddo
      endif
c      dic0 = 0.0
c      dic1 = 0.0
c      ta0 = 0.0
c      ta1 = 0.0

C        TRANSFER I.C. TO APPROPRIATE VECTORS
      CALL LOADSED

C        TEMPERATURE DIFFUSION COEFFICIENT

      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  DIFFT
      DIFFT = 0.0001*86400.*DIFFT
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  SALTSW
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  FRPPH1
      IF(CONST(2).GE.2) READ(40,1080)   COMMENT
      IF(CONST(2).GE.2) READ(40,1200,ERR=900)  FRPPH2
      IF(CONST(2).GE.3) READ(40,1080)   COMMENT
      IF(CONST(2).GE.3) READ(40,1200,ERR=900)  FRPPH3
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  FRPOP
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  FRNPH1
      IF(CONST(2).GE.2) READ(40,1080)   COMMENT
      IF(CONST(2).GE.2) READ(40,1200,ERR=900)  FRNPH2
      IF(CONST(2).GE.3) READ(40,1080)   COMMENT
      IF(CONST(2).GE.3) READ(40,1200,ERR=900)  FRNPH3
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  FRPON
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  FRCPH1
      IF(CONST(2).GE.2) READ(40,1080)   COMMENT
      IF(CONST(2).GE.2) READ(40,1200,ERR=900)  FRCPH2
      IF(CONST(2).GE.3) READ(40,1080)   COMMENT
      IF(CONST(2).GE.3) READ(40,1200,ERR=900)  FRCPH3
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  FRPOC
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  (KPDIAG(I),DPTHTA(I),I=1,3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  (KNDIAG(I),DNTHTA(I),I=1,3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  (KCDIAG(I),DCTHTA(I),I=1,3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  KSI,THTASI

      WRITE(OUT,1210)  DIFFT
      WRITE(OUT,1215)  SALTSW            ! New CH4 formation mechanism
      WRITE(OUT,1220)
      WRITE(OUT,1230)  FRNAME(1),FRPPH1
      WRITE(OUT,1230)  FRNAME(2),FRPPH2
      WRITE(OUT,1230)  FRNAME(3),FRPOP
      WRITE(OUT,1230)  FRNAME(4),FRNPH1
      WRITE(OUT,1230)  FRNAME(5),FRNPH2
      WRITE(OUT,1230)  FRNAME(6),FRPON
      WRITE(OUT,1230)  FRNAME(7),FRCPH1
      WRITE(OUT,1230)  FRNAME(8),FRCPH2
      WRITE(OUT,1230)  FRNAME(9),FRPOC
      WRITE(OUT,1240)  (KPDIAG(I),DPTHTA(I),I=1,3),(KNDIAG(I),
     .   DNTHTA(I),I=1,3),(KCDIAG(I),DCTHTA(I),I=1,3),KSI,THTASI

 1210 FORMAT(
     .   //30X,'TEMPERATURE DIFFUSION COEFFICIENT ',E10.3,' CM**2/SEC')
 1215 FORMAT(
     .   //30X,'SALINITY SWITCH:      ',F10.3,' PPT')
 1220 FORMAT(//30X,'PARTICULATE ORGANIC MATTER G-MODEL SPLITS'/ 10X,
     .   'FRACTION OF....',5X,'RECYCLED TO',5X,'G1',5X,'G2',5X,'G3')
 1230 FORMAT(6X,A24,11X,3F7.2)
 1240 FORMAT(//30X,'DIAGENESIS RATES (/DAY)  TEMP CORR FACTOR'/
     .   30X,'PHOSPHORUS'/
     .   39X,'G1',E11.3,5X,F7.3/39X,'G2',E11.3,5X,F7.3/
     .   39X,'G3',E11.3,5X,F7.3/
     .   30X,'NITROGEN'/
     .   39X,'G1',E11.3,5X,F7.3/39X,'G2',E11.3,5X,F7.3/
     .   39X,'G3',E11.3,5X,F7.3/
     .   30X,'CARBON'/
     .   39X,'G1',E11.3,5X,F7.3/39X,'G2',E11.3,5X,F7.3/
     .   39X,'G3',E11.3,5X,F7.3/
     .   30X,'SILICA'/
     .   41X,E11.3,5X,F7.3)

C        READ SEDIMENTATION AND PARTICLE MIXING RATES
      WRITE(OUT,1050)
 1050 FORMAT(////38X,
     .  'S E D I M E N T A T I O N   A N D   P A R T I C L E   M I X I N
     . G   R A T E S'//)

      CALL VREAD(NX,NY,OUT,ISEDPRNT,VSED,1)
      CALL VREAD(NX,NY,OUT,ISEDPRNT,VPMIX,2)
      CALL VREAD(NX,NY,OUT,ISEDPRNT,VDMIX,3)

C        COMPUTE G FRACTION SPLITS FOR REFRACTORY POM
      DO 35 I=2,3
       PGFRAC(I) = FRPOP(I)/(FRPOP(2)+FRPOP(3))
       NGFRAC(I) = FRPON(I)/(FRPON(2)+FRPON(3))
       CGFRAC(I) = FRPOC(I)/(FRPOC(2)+FRPOC(3))
   35 CONTINUE

C        CONVERT CM/YEAR TO M/DAY

      DO 37 IY=1,NY
       DO 37 IX=1,NX
        VSED(IX,IY) = 2.73791E-5*VSED(IX,IY)
   37 CONTINUE

      READ(40,1080)   COMMENT
      READ(40,1103,ERR=900)  m1m2opt,piesiopt,piepo4opt
 1103 FORMAT(3I10)
      IF(M1M2OPT.EQ.0) THEN
        READ(40,1080)   COMMENT
        READ(40,1200,ERR=900)  m1g,m2g,thtaDp,thtaDd
        WRITE(OUT,1269)   M1g,M2g,thtaDp,thtaDd
 1269   FORMAT(//35X,'ADDITIONAL CONSTANTS'/
     .    /30X,'m1........',F8.2,' KG/L'
     .    /30X,'m2........',F8.2,' KG/L'
     .    /30X,'thtaDp....',F8.3,
     .    /30X,'thtaDd....',F8.3)
        DO 38 IY=1,NY
         DO 38 IX=1,NX
          IF(FSM(IX,IY).EQ.1) THEN
            M1(IX,IY) = M1G
            M2(IX,IY) = M2G
          ENDIF
 38     CONTINUE
      ELSE
        CALL VREAD(NX,NY,OUT,ISEDPRNT,M1,4)
        CALL VREAD(NX,NY,OUT,ISEDPRNT,M2,5)
        READ(40,1080)   COMMENT
        READ(40,1200,ERR=900)  thtaDp,thtaDd
        WRITE(OUT,1270)   thtaDp,thtaDd
 1270   FORMAT(
     .    /30X,'thtaDp....',F8.3,
     .    /30X,'thtaDd....',F8.3)
      ENDIF
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  Dd0,thtaDd0
      WRITE(OUT,1271)   Dd0,thtaDd0
 1271 FORMAT(
     .    /30X,'Dd0.......',F8.3,' M/DAY'
     .    /30X,'thtaDd0...',F8.3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  kappnh4s,pienh4,thtanh4s,kmnh4,thtakmnh4,
     .     kmnh4o2
      WRITE(OUT,1272)  kappnh4s,pienh4,thtanh4s,kmnh4,thtakmnh4,kmnh4o2
 1272 FORMAT(
     .    /30X,'Saltwater values'
     .    /30X,'kappnh4s..',F8.3,' M/DAY'
     .    /30X,'pienh4....',F8.3,' L/KG'
     .    /30X,'thtanh4s..',F8.3,
     .    /30X,'kmnh4.....',F8.3,' MG N/M**3'
     .    /30X,'thtakmnh4.',F8.3,
     .    /30X,'kmnh4o2...',F8.3,' MG O2/L')
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  kappnh4f,thtanh4f
      WRITE(OUT,1273)  kappnh4f,thtanh4f
 1273 FORMAT(
     .    /30X,'Freshwater values'
     .    /30X,'kappnh4f..',F8.3,' M/DAY'
     .    /30X,'thtanh4f..',F8.3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  kapp1no3s,k2no3s,thtano3s
      WRITE(OUT,1274)  kapp1no3s,k2no3s,thtano3s
 1274 FORMAT(
     .    /30X,'Saltwater values'
     .    /30X,'kapp1no3s.',F8.3,' M/DAY'
     .    /30X,'k2no3s....',F8.3,' /DAY'
     .    /30X,'thtano3s..',F8.3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  kapp1no3f,k2no3f,thtano3f
      WRITE(OUT,1275)  kapp1no3f,k2no3f,thtano3f
 1275 FORMAT(
     .    /30X,'Freshwater values'
     .    /30X,'kapp1no3f.',F8.3,' M/DAY'
     .    /30X,'k2no3f....',F8.3,' /DAY'
     .    /30X,'thtano3f..',F8.3)
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  kappd1,kappp1,pie1s,pie2s,thtapd1,kmhso2
      WRITE(OUT,1276)  kappd1,kappp1,pie1s,pie2s,thtapd1,kmhso2
 1276 FORMAT(
     .    /30X,'kappd1....',F8.3,' M/DAY'
     .    /30X,'kappp1....',F8.3,' M/DAY'
     .    /30X,'pie1s.....',F8.3,' L/KG'
     .    /30X,'pie2s.....',F8.3,' L/KG'
     .    /30X,'thtapd1...',F8.3,
     .    /30X,'kmhso2....',F8.3,' MG O2/L')
      IF(PIESIOPT.EQ.0) THEN
        READ(40,1080)   COMMENT
        READ(40,1200,ERR=900)  csisat,pie1sig,pie2sig,ksi,thtasi,kmpsi,
     .       o2critsi,jsidetr
        WRITE(OUT,1278)  csisat,pie1sig,pie2sig,ksi,thtasi,kmpsi,
     .       o2critsi,jsidetr
 1278 FORMAT(
     .    /30X,'csisat....',F8.3,' MG SI/M**3'
     .    /30X,'pie1si....',F8.3,' L/KG'
     .    /30X,'pie2si....',F8.3,' L/KG'
     .    /30X,'ksi.......',F8.3,' /DAY'
     .    /30X,'thtasi....',F8.3,
     .    /30X,'kmpsi.....',F8.3,' MG SI/M**3'
     .    /30X,'o2critsi..',F8.3,' MG O2/L'
     .    /30X,'jsidetr...',F8.1,' MG SI/M**2')
        DO 39 IY=1,NY
         DO 39 IX=1,NX
          IF(FSM(IX,IY).EQ.1) THEN
           pie1si(IX,IY) = pie1sig
           pie2si(IX,IY) = pie2sig
          ENDIF
 39     CONTINUE
      ELSE
        CALL VREAD(NX,NY,OUT,ISEDPRNT,pie1si,6)
        CALL VREAD(NX,NY,OUT,ISEDPRNT,pie2si,7)
        READ(40,1080)   COMMENT
        READ(40,1200,ERR=900)  csisat,ksi,thtasi,kmpsi,
     .     o2critsi,jsidetr
        WRITE(OUT,1378)  csisat,ksi,thtasi,kmpsi,
     .     o2critsi,jsidetr
 1378 FORMAT(
     .    /30X,'csisat....',F8.3,' MG SI/M**3'
     .    /30X,'ksi.......',F8.3,' /DAY'
     .    /30X,'thtasi....',F8.3,
     .    /30X,'kmpsi.....',F8.3,' MG SI/M**3'
     .    /30X,'o2critsi..',F8.3,' MG O2/L'
     .    /30X,'jsidetr...',F8.1,' MG SI/M**2')
      ENDIF
      IF(PIEPO4OPT.EQ.0) THEN
        READ(40,1080)   COMMENT
        READ(40,1200,ERR=900)  pie1po4mg,pie1po4ng,o2crit,kmo2Dp
        WRITE(OUT,1280)  pie1po4mg,pie1po4ng,o2crit,kmo2Dp
 1280   FORMAT(
     .    /30X,'pie1po4m..',F8.3,' L/KG'
     .    /30X,'pie1po4n..',F8.3,' L/KG'
     .    /30X,'o2crit....',F8.3,' MG O2/L'
     .    /30X,'kmo2Dp....',F8.3,' MG O2/L')
        DO 41 IY=1,NY
         DO 41 IX=1,NX
          IF(FSM(IX,IY).EQ.1) THEN
           pie1po4m(IX,IY) = pie1po4mg
           pie1po4n(IX,IY) = pie1po4ng
          ENDIF
 41     CONTINUE
      ELSE
        CALL VREAD(NX,NY,OUT,ISEDPRNT,pie1po4m,8)
        CALL VREAD(NX,NY,OUT,ISEDPRNT,pie1po4n,9)
        READ(40,1080)   COMMENT
        READ(40,1200,ERR=900)  o2crit,kmo2Dp
        WRITE(OUT,1380)  o2crit,kmo2Dp
 1380   FORMAT(
     .    /30X,'o2crit....',F8.3,' MG O2/L'
     .    /30X,'kmo2Dp....',F8.3,' MG O2/L')
      ENDIF
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  tempbnth,kbnthstr,klbnth,Dpmin
      WRITE(OUT,1282)  tempbnth,kbnthstr,klbnth,Dpmin
 1282 FORMAT(
     .  /30X,'tempbnth..',F8.3,' DEG C'
     .  /30X,'kbnthstr..',F8.3,' /DAY'
     .  /30X,'klbnth....',F8.3,
     .  /30X,'Dpmin.....',F8.3,' M**2/DAY')
      READ(40,1080)   COMMENT
      READ(40,1200,ERR=900)  kappch4,thtach4,KMCH4O2,KMSO4   !CH4-SO4 Code
      WRITE(OUT,1284)  kappch4,thtach4,KMCH4O2,KMSO4         !CH4-SO4 Code
 1284 FORMAT(
     .  /30X,'kappch4...',F8.3,' M/DAY'
     .  /30X,'thtach4...',F8.3,
     .  /30X,'KmCH4O2...',F8.3,
     .  /30X,'KmSO4  ...',F8.3)

C        SET MIN AND MAX RANGES FOR ROOT FINDER

      SODMIN = 0.0001
      SODMAX = 100.
      ITER = 50
      EPS = 0.00005

      DO 46 IY=1,NY
       DO 46 IX=1,NX
        IF(FSM(IX,IY).EQ.1.)  THEN
C        SET UP REACTION RATES IN TABLE LOOK-UP FORM
C        (ASSUMES ALL SEDIMENT DEPTHS AND MIXING RATES ARE THE SAME)
         H20 = HSED(IX,IY)
         DELTAZ = HSED(IX,IY)
         w2 = VSED(IX,IY)
         Dp = VPMIX(IX,IY)
         Dd = VDMIX(IX,IY)
         CALL SEDTS2TH
C        RETURN TO MAIN PROGRAM AFTER INITIALIZATION
         RETURN
        ENDIF
   46 CONTINUE

   47 RETURN


      ENTRY  SED_FLUX


C     Increment counter
      IGDCT = IGDCT + 1


C        PASS WQ TIME-STEP (IN DAYS) TO SEDIMENT SUBR
      IF(DTWQ.NE.0.0)  THEN
        DELTAT = DTWQ
      ELSE
        DELTAT = DT
      ENDIF

C      COMPUTE TEMPERATURE AND BENTHIC-STRESS DEPENDENT DEPOSITION RATES
C     DO 55 IY=1,NY
C     DO 55 IX=1,NX
C     IF(FSM(IX,IY).NE.1.) GO TO 55
C     TEMPD = CTEMP(IX,IY)
C     TEMPD = MIN(CTEMP(IX,IY),25.)
C     ITEMP = 10.*TEMPD+1
C     TERM = (1.-KBNTHSTR*BNTHSTR1S(IX,IY))*ZW12NOM(ITEMP)
C  55 CONTINUE

C             PARTICULATE ORGANIC FLUX TO SEDIMENT
c  c$doacross local(iy,ix,term,ppflx,pnflx,pcflx,sarea)

      DO 100 IY=1,NY
      DO 100 IX=1,NX
      IF(FSM(IX,IY).NE.1.) GO TO 100

C        FLUX RATE FROM WATER COLUMN = DEPOSITION RATE * CONC
C                    M/DAY * MG/L * 1000 L/M**3 ==> MG/M2-DAY 

      TERM = 1000./XAZ(IX,IY)
      SSFWS(IX,IY) = 0.
      PCFWS(IX,IY) = 0.
      PNFWS(IX,IY) = 0.
      PPFWS(IX,IY) = 0.
      PSFWS(IX,IY) = 0.

      DO 60 I=1,3

C USE LABILE PARTICULATES AS G1, SPLIT REFRACTORY BETWEEN G2, G3

      IF (I .EQ. 1) THEN
        PPFLX = TERM*DEPFLUX(IX,IY,5)
        PNFLX = TERM*DEPFLUX(IX,IY,7)
        PCFLX = TERM*(DEPFLUX(IX,IY,10)+DEPFLUX(IX,IY,11))
      ELSE
        PPFLX = TERM*DEPFLUX(IX,IY,4)*PGFRAC(I)
        PNFLX = TERM*DEPFLUX(IX,IY,6)*NGFRAC(I)
        PCFLX = TERM*DEPFLUX(IX,IY,9)*CGFRAC(I)
      END IF

C          PHOSPHORUS

      FLXPOP(IX,IY,I)
     .        = TERM*PCRB1(IX,IY,NZ)*FRPPH1(I)*DEPFLUX(IX,IY,1)
     .        + TERM*PCRB2(IX,IY,NZ)*FRPPH2(I)*DEPFLUX(IX,IY,2)
     .        + TERM*PCRB3(IX,IY,NZ)*FRPPH3(I)*DEPFLUX(IX,IY,3)
     .        + PPFLX

C        NITROGEN

      FLXPON(IX,IY,I)
     .        = TERM*NCRB1(IX,IY,NZ)*FRNPH1(I)*DEPFLUX(IX,IY,1)
     .        + TERM*NCRB2(IX,IY,NZ)*FRNPH2(I)*DEPFLUX(IX,IY,2)
     .        + TERM*NCRB3(IX,IY,NZ)*FRNPH3(I)*DEPFLUX(IX,IY,3)
     .        + PNFLX

C        CARBON

      FLXPOC(IX,IY,I) = TERM*FRCPH1(I)*DEPFLUX(IX,IY,1)
     .                + TERM*FRCPH2(I)*DEPFLUX(IX,IY,2)
     .                + TERM*FRCPH3(I)*DEPFLUX(IX,IY,3)
     .                + PCFLX

   60 CONTINUE

C        SILICA

      FLXPOS(IX,IY) = TERM*SCRB1(IX,IY,NZ)*DEPFLUX(IX,IY,1)
     .              + TERM*SCRB2(IX,IY,NZ)*DEPFLUX(IX,IY,2)
     .              + TERM*SCRB3(IX,IY,NZ)*DEPFLUX(IX,IY,3)
     .              + TERM*DEPFLUX(IX,IY,8)

C WILL ADD ADSORBED PHOSPHATE AND SILICA TO LAYER 2 INORGANIC POOL

      FLXPO4T2S(IX,IY) = TERM*DEPFLUX(IX,IY,12)
      FLXSIT2S(IX,IY) = TERM*DEPFLUX(IX,IY,13)
      
      FLXDICT2S(IX,IY) = TERM*DEPFLUX(IX,IY,14)
      FLXTAT2S(IX,IY) = TERM*DEPFLUX(IX,IY,15)

C SUM PARTICULATE FLUXES TO SEDIMENTS, NEGATIVE INTO SEDIMENTS

      SAREA=1./XAZ(IX,IY)
c     SAREA=XAZ(IX,IY)/1000.
      DO 62 I = 1,3
        PPFWS(IX,IY) = PPFWS(IX,IY) - FLXPOP(IX,IY,I)*SAREA
        PNFWS(IX,IY) = PNFWS(IX,IY) - FLXPON(IX,IY,I)*SAREA
        PCFWS(IX,IY) = PCFWS(IX,IY) - FLXPOC(IX,IY,I)*SAREA
 62   CONTINUE
        PSFWS(IX,IY) = PSFWS(IX,IY) - FLXPOS(IX,IY)*SAREA

100   CONTINUE

      DO 200 IY=1,NY
      DO 200 IX=1,NX
      IF(FSM(IX,IY).NE.1.)  GO TO 200

C        LOAD CONCS AT TIME LEVEL   T - DELTAT

      CALL LOADTM1
C       -DIAGX- UNITS ARE [GM/M**3-DAY]
C       -JDIAGX- UNITS ARE [MG/M2/DAY]

C        SEDIMENT DEPTH
      H2 = HSED(IX,IY)
      DELTAZ = H2

C        SEDIMENTATION AND MIXING RATES
      w2 = VSED(IX,IY)
      Dp = VPMIX(IX,IY)
      Dd = VDMIX(IX,IY)

C        SEDIMENT TEMPERATURE
      TEMPD = MAX(CTEMP(IX,IY),0.0)
      TEMPD = MIN(TEMPD,34.9)
      STP20 = TEMPD - 20.

C        CONVERT OVERLYING WATER COLUMN CONCENTRATIONS INTO MG/M**3
      PO40 = PO4(IX,IY,NZ)*1000.
      NH40 = NH4(IX,IY,NZ)*1000.
      NO30 = NO23(IX,IY,NZ)*1000.
      SI0  = SI(IX,IY,NZ)*1000.
      DIC0 = DIC(IX,IY,NZ)*1000. 
c     mmol/l to  mol/l
      TA0 = TALK(IX,IY,NZ)*1000.

C        OVERLYING WATER COLUMN OXYGEN (WITH 0.001 MG/L FLOOR)
      O20 = AMAX1(DO(IX,IY,NZ),0.001)
      HS0 = 0.

C        OVERLYING WATER COLUMN SALINITY
      SAL0 = SAL(IX,IY,NZ)

C       Regression function to get SO4 concentration from SAL0
C       [SO4] = 20 mg/L            for        [Cl] < 6 mg/L
C             = 20 + (27/190)[Cl]  for        [Cl] > 6 mg/L
C       1 ppt = 607.445 mg/L Cl
        IF (SAL0 .GT. 0.0099) THEN
           SO40MG = 20.0 + (27./190.)*607.445*SAL0
        ELSE
           SO40MG = 20.0
        ENDIF

C        DEPTH OF OVERLYING WATER COLUMN
      OWCDPH = HBAR(IX,IY)
c     OWCDPH = ZD(IX,IY,NZ)+HBAR(IX,IY)

C        METHANE SATURATION
      CH4SAT = 99.0*(1.+(OWCDPH+HSED(IX,IY))/10.)*0.9759**STP20

C        SORBED PO4 AND SI FLUXES
      FLXPO4T2 = FLXPO4T2S(IX,IY)
      FLXSIT2 = FLXSIT2S(IX,IY)
      
      FLXDICT2 = FLXDICT2S(IX,IY)
      FLXTAT2  = FLXTAT2S(IX,IY)

c-- evaluate the temperature dependent coefficients
      error=sedf(sodmin,0)
c-- diagenesis computation
      call diagen

d        write(out,9911)  time,ctemp(ix,iy),pop1r,pop2r,pop3r,popr
d9911    format(/1x,'entry Zbrent time,ctemp,pop1r,pop2r,pop3r,popr'/
d    .          8e10.3)
d        write(out,9901)  pon1r,pon2r,pon3r,ponr,poc1r,poc2r,poc3r,pocr
d        write(out,9902)  po4t2r,hst2r,sit2r,psiavr
d        write(out,9903)
d    .        (flxpop(ix,iy,1)+flxpop(ix,iy,2)+flxpop(ix,iy,3))
d    .       ,(flxpon(ix,iy,1)+flxpon(ix,iy,2)+flxpon(ix,iy,3))
d    .       ,(flxpoc(ix,iy,1)+flxpoc(ix,iy,2)+flxpoc(ix,iy,3))
d        write(out,9904)  o20,csod,sod,s
d        write(out,9905)  jp,jn,jc,jo2nh4,xjco2av,xjc1av
d        write(out,9906)  jpo4,jnh4,jno3,jhs,jsi,jch4aq,jch4g,bnthstr
d        write(out,9907)  po40,po41,po42,po4t2,nh40,nh41,nh42,nh4t2
d        write(out,9908)  no30,no31,no32,no3t2,hs1,hs2,hst2
d        write(out,9909)  si0,si1,si2,sit2
c-- solve the nh4,no3,sod equations
      sod=zbrent(sedf,sodmin,sodmax,eps,ierr)
      if(ierr.ne.0)  then
         write(out,9000)   time,ierr,ix,iy,sodmin,sodmax
 9000    format(/
     .      5x,'Zbrent failure at time =',f8.3,' with ierr=',i2/
     .      5x,'in sediment segment ix=',i5,' iy=',i5/
     .      5x,'(sodmin,sodmax=',f6.4,f6.1,')'/
     .      5x,'Program termination follows diagnostic dumps')
         if(ierr.eq.2)  then
         write(out,9900)  time,ctemp(ix,iy),pop1r,pop2r,pop3r,popr
 9900    format(/1x,' time,ctemp(ix,iy),pop1r,pop2r,pop3r,popr'/8e10.3)
         write(out,9901)  pon1r,pon2r,pon3r,ponr,poc1r,poc2r,poc3r,pocr
 9901    format(/1x,' pon1r,pon2r,pon3r,ponr,poc1r,poc2r,poc3r,pocr'/
     .          8e10.3)
         write(out,9902)  po4t2r,hst2r,sit2r,psiavr
 9902    format(/1x,' po4t2r,hst2r,sit2r,psiavr'/8e10.3)
         write(out,9903)
     .        (flxpop(ix,iy,1)+flxpop(ix,iy,2)+flxpop(ix,iy,3))
     .       ,(flxpon(ix,iy,1)+flxpon(ix,iy,2)+flxpon(ix,iy,3))
     .       ,(flxpoc(ix,iy,1)+flxpoc(ix,iy,2)+flxpoc(ix,iy,3))
 9903    format(/1x,' flxpop,flxpon,flxpoc'/8e10.3)
         write(out,9904)  o20,csod,sod,s,h1,h2,h1dot,h2dot
 9904    format(/1x,' o20,csod,sod,s,h1,h2,h1dot,h2dot'/8e10.3)
         write(out,9905)  jp,jn,jc,jo2nh4,xjco2av,xjc1av
 9905    format(/1x,' jp,jn,jc,jo2nh4,xjco2av,xjc1av'/8e10.3)
         write(out,9906)  jpo4,jnh4,jno3,jhs,jsi,jch4aq,jch4g,bnthstr
 9906    format(/1x,' jpo4,jnh4,jno3,jhs,jsi,jch4aq,jch4g,bnthstr'/
     .          8e10.3)
         write(out,9907)  po40,po41,po42,po4t2,nh40,nh41,nh42,nh4t2
 9907    format(/1x,' po40,po41,po42,po4t2,nh40,nh41,nh42,nh4t2'/8e10.3)
         write(out,9908)  no30,no31,no32,no3t2,hs1,hs2,hst2
 9908    format(/1x,' no30,no31,no32,no3t2,hs1,hs2,hst2'/8e10.3)
         write(out,9909)  si0,si1,si2,sit2
 9909    format(/1x,' si0,si1,si2,sit2'/8e10.3)
         call exit
         else
         error=sedf(sodmin,1)
         write(out,9889)  time,sodmin,error
 9889    format(/5x,'Zbrent diagnostics at time =',f8.3,
     .     ' for sod =',f8.4,' error =',e12.3/)
         write(out,9900)  time,ctemp(ix,iy),pop1r,pop2r,pop3r,popr
         write(out,9901)  pon1r,pon2r,pon3r,ponr,poc1r,poc2r,poc3r,pocr
         write(out,9902)  po4t2r,hst2r,sit2r,psiavr
         write(out,9903)
     .        (flxpop(ix,iy,1)+flxpop(ix,iy,2)+flxpop(ix,iy,3))
     .       ,(flxpon(ix,iy,1)+flxpon(ix,iy,2)+flxpon(ix,iy,3))
     .       ,(flxpoc(ix,iy,1)+flxpoc(ix,iy,2)+flxpoc(ix,iy,3))
         write(out,9904)  o20,csod,sod,s,h1,h2,h1dot,h2dot
         write(out,9905)  jp,jn,jc,jo2nh4,xjco2av,xjc1av
         write(out,9906)  jpo4,jnh4,jno3,jhs,jsi,jch4aq,jch4g,bnthstr
         write(out,9907)  po40,po41,po42,po4t2,nh40,nh41,nh42,nh4t2
         write(out,9908)  no30,no31,no32,no3t2,hs1,hs2,hst2
         write(out,9909)  si0,si1,si2,sit2
         error=sedf(sodmax,1)
         write(out,9889)  time,sodmax,error
         write(out,9900)  time,ctemp(ix,iy),pop1r,pop2r,pop3r,popr
         write(out,9901)  pon1r,pon2r,pon3r,ponr,poc1r,poc2r,poc3r,pocr
         write(out,9902)  po4t2r,hst2r,sit2r,psiavr
         write(out,9903)
     .        (flxpop(ix,iy,1)+flxpop(ix,iy,2)+flxpop(ix,iy,3))
     .       ,(flxpon(ix,iy,1)+flxpon(ix,iy,2)+flxpon(ix,iy,3))
     .       ,(flxpoc(ix,iy,1)+flxpoc(ix,iy,2)+flxpoc(ix,iy,3))
         write(out,9904)  o20,csod,sod,s,h1,h2,h1dot,h2dot
         write(out,9905)  jp,jn,jc,jo2nh4,xjco2av,xjc1av
         write(out,9906)  jpo4,jnh4,jno3,jhs,jsi,jch4aq,jch4g,bnthstr
         write(out,9907)  po40,po41,po42,po4t2,nh40,nh41,nh42,nh4t2
         write(out,9908)  no30,no31,no32,no3t2,hs1,hs2,hst2
         write(out,9909)  si0,si1,si2,sit2
         call exit
         endif
      endif
      h1s(ix,iy)=h1
c-- and evaluate the po4,si equations
      error=sedf(sod,2)

c-- replace the t minus 1 concentrations
      call storetm1
C Burial Flux Calculation added by LB
      DO IG = 1,3        
           IF(FSM(IX,IY).EQ.1.)  THEN
            CHDIC(IX,IY,IG+2)=VSED(IX,IY)*CPOC(IX,IY,IG)
           ENDIF
      ENDDO
C ASSIGN FLUX-FLUX RESULTS TO WQM ARRAYS

      SODS(IX,IY) = SOD
      FLXHS(IX,IY) = JHS
      FLXCH4AQ(IX,IY) = JCH4AQ
      FLXCH4G(IX,IY) = JCH4G
      FLXNH4(IX,IY) = JNH4/1000.
      FLXPO4(IX,IY) = JPO4/1000.
      FLXNO3(IX,IY) = JNO3/1000.
      FLXSI(IX,IY) = JSI/1000.
      FLXDIC(IX,IY) = JDIC/1000.
      FLXPOCN(IX,IY) = (FLXPOC(IX,IY,1)+FLXPOC(IX,IY,2)+FLXPOC(IX,IY,3))
     .                 /XAZ(IX,IY)
C      if(jdic.gt.100.0) then
C          pause
C      endif
      FLXTALK(IX,IY) = JTALK/1000.
C    dic0单位mmol/m3或umol/L,CDDIC中SA单位m3,需要JDIC单位为mmol/m2/day,所以不需要换算 

c- denitrification flux
      itemp = 10*tempd+1
      if(sal0.gt.saltsw)  then
         xapp1no3=zhtano3s(itemp)
         xk2no3=zhta2no3s(itemp)*h2
      else
         xapp1no3=zhtano3f(itemp)
         xk2no3=zhta2no3f(itemp)*h2
      endif
      bnthden = (xapp1no3*xapp1no3/s*no31 + xk2no3*no32)

C        Perform Global Dump Sums, if required
      IF(IGDSEDOPT.EQ.1) THEN
            M2G=M2(IX,IY)
            CTEMP_GDA(IX,IY)  = CTEMP_GDA(IX,IY) + CTEMP(IX,IY)
            POP1R_GDA(IX,IY)  = POP1R_GDA(IX,IY) + POP1AV/1.0E6/M2G
            POP2R_GDA(IX,IY)  = POP2R_GDA(IX,IY) + POP2AV/1.0E6/M2G
            POP3R_GDA(IX,IY)  = POP3R_GDA(IX,IY) + POP3AV/1.0E6/M2G
            POPR_GDA(IX,IY)   = POPR_GDA(IX,IY) + (POP1AV+POP2AV+POP3AV)
     .                          /1.0E6/M2G
            PON1R_GDA(IX,IY)  = PON1R_GDA(IX,IY) + PON1AV/1.0E6/M2G
            PON2R_GDA(IX,IY)  = PON2R_GDA(IX,IY) + PON2AV/1.0E6/M2G
            PON3R_GDA(IX,IY)  = PON3R_GDA(IX,IY) + PON3AV/1.0E6/M2G
            PONR_GDA(IX,IY)   = PONR_GDA(IX,IY) + (PON1AV+PON2AV+PON3AV)
     .                          /1.0E6/M2G
            POC1R_GDA(IX,IY)  = POC1R_GDA(IX,IY) + POC1AV/1.0E6/M2G
            POC2R_GDA(IX,IY)  = POC2R_GDA(IX,IY) + POC2AV/1.0E6/M2G
            POC3R_GDA(IX,IY)  = POC3R_GDA(IX,IY) + POC3AV/1.0E6/M2G
            POCR_GDA(IX,IY)   = POCR_GDA(IX,IY) + (POC1AV+POC2AV+POC3AV)
     .                          /1.0E6/M2G
            PO4T2R_GDA(IX,IY) = PO4T2R_GDA(IX,IY) + 
     .                          PO4T2AV/1000./M2G
            HST2R_GDA(IX,IY)  = HST2R_GDA(IX,IY) + HST2AV/1000./M2G
            SIT2R_GDA(IX,IY)  = SIT2R_GDA(IX,IY) + (PSIAV+SIT2AV)
     .                          /1.0E6/M2G
            PSIAVR_GDA(IX,IY) = PSIAVR_GDA(IX,IY) + PSIAV/1.0E6/M2G
            FLXPOP_GDA(IX,IY) = FLXPOP_GDA(IX,IY) + FLXPOP(IX,IY,1)
     .                        + FLXPOP(IX,IY,2) + FLXPOP(IX,IY,3)
            FLXPOP1_GDA(IX,IY) = FLXPOP1_GDA(IX,IY) + FLXPOP(IX,IY,1)
            FLXPOP2_GDA(IX,IY) = FLXPOP2_GDA(IX,IY) + FLXPOP(IX,IY,2)
            FLXPOP3_GDA(IX,IY) = FLXPOP3_GDA(IX,IY) + FLXPOP(IX,IY,3)
            FLXPON_GDA(IX,IY) = FLXPON_GDA(IX,IY) + FLXPON(IX,IY,1)
     .                        + FLXPON(IX,IY,2) + FLXPON(IX,IY,3)
            FLXPON1_GDA(IX,IY) = FLXPON1_GDA(IX,IY) + FLXPON(IX,IY,1)
            FLXPON2_GDA(IX,IY) = FLXPON2_GDA(IX,IY) + FLXPON(IX,IY,2)
            FLXPON3_GDA(IX,IY) = FLXPON3_GDA(IX,IY) + FLXPON(IX,IY,3)
            FLXPOC_GDA(IX,IY) = FLXPOC_GDA(IX,IY) + FLXPOC(IX,IY,1)
     .                        + FLXPOC(IX,IY,2) + FLXPOC(IX,IY,3)
            FLXPOC1_GDA(IX,IY) = FLXPOC1_GDA(IX,IY) + FLXPOC(IX,IY,1)
            FLXPOC2_GDA(IX,IY) = FLXPOC2_GDA(IX,IY) + FLXPOC(IX,IY,2)
            FLXPOC3_GDA(IX,IY) = FLXPOC3_GDA(IX,IY) + FLXPOC(IX,IY,3)
            FLXPOS_GDA(IX,IY) = FLXPOS_GDA(IX,IY) + FLXPOS(IX,IY)
            FLXPO4T2S_GDA(IX,IY) = FLXPO4T2S_GDA(IX,IY)+FLXPO4T2S(IX,IY)
            FLXSIT2S_GDA(IX,IY) = FLXSIT2S_GDA(IX,IY) + FLXSIT2S(IX,IY)
            O20_GDA(IX,IY)    = O20_GDA(IX,IY) + O20
            CSOD_GDA(IX,IY)   = CSOD_GDA(IX,IY) + CSOD
            SOD_GDA(IX,IY)    = SOD_GDA(IX,IY) + SOD
            S_GDA(IX,IY)      = S_GDA(IX,IY) + S
            JP_GDA(IX,IY)     = JP_GDA(IX,IY) + JP
            JN_GDA(IX,IY)     = JN_GDA(IX,IY) + JN
            JC_GDA(IX,IY)     = JC_GDA(IX,IY) + JC
            JO2NH4_GDA(IX,IY) = JO2NH4_GDA(IX,IY) + JO2NH4
            XJCO2AV_GDA(IX,IY) = XJCO2AV_GDA(IX,IY) + XJCO2AV
            XJC1AV_GDA(IX,IY) = XJC1AV_GDA(IX,IY) + XJC1AV
            JPO4_GDA(IX,IY)   = JPO4_GDA(IX,IY) + JPO4
            JNH4_GDA(IX,IY)   = JNH4_GDA(IX,IY) + JNH4
            JNO3_GDA(IX,IY)   = JNO3_GDA(IX,IY) + JNO3
            JH2S_GDA(IX,IY)    = JH2S_GDA(IX,IY) + JHS
            JSI_GDA(IX,IY)    = JSI_GDA(IX,IY) + JSI
            JCH4AQ_GDA(IX,IY) = JCH4AQ_GDA(IX,IY) + JCH4AQ
            JCH4G_GDA(IX,IY)  = JCH4G_GDA(IX,IY) + JCH4G
            H1_GDA(IX,IY)     = H1_GDA(IX,IY) + H1
            PO40_GDA(IX,IY)   = PO40_GDA(IX,IY) + PO40
            PO41_GDA(IX,IY)   = PO41_GDA(IX,IY) + PO41
            PO42_GDA(IX,IY)   = PO42_GDA(IX,IY) + PO42
            PO4T2_GDA(IX,IY)  = PO4T2_GDA(IX,IY) + PO4T2
            NH40_GDA(IX,IY)   = NH40_GDA(IX,IY) + NH40
            NH41_GDA(IX,IY)   = NH41_GDA(IX,IY) +NH41
            NH42_GDA(IX,IY)   = NH42_GDA(IX,IY) + NH42
            NH4T2_GDA(IX,IY)  = NH4T2_GDA(IX,IY) + NH4T2
            NO30_GDA(IX,IY)   = NO30_GDA(IX,IY) + NO30
            NO31_GDA(IX,IY)   = NO31_GDA(IX,IY) + NO31
            NO32_GDA(IX,IY)   = NO32_GDA(IX,IY) + NO32
            NO3T2_GDA(IX,IY)  = NO3T2_GDA(IX,IY) + NO3T2
            HS1_GDA(IX,IY)    = HS1_GDA(IX,IY) + HS1
            HS2_GDA(IX,IY)    = HS2_GDA(IX,IY) + HS2
            HST2_GDA(IX,IY)   = HST2_GDA(IX,IY) + HST2
            SI0_GDA(IX,IY)    = SI0_GDA(IX,IY) + SI0
            SI1_GDA(IX,IY)    = SI1_GDA(IX,IY) + SI1
            SI2_GDA(IX,IY)    = SI2_GDA(IX,IY) + SI2
            SIT2_GDA(IX,IY)   = SIT2_GDA(IX,IY) + SIT2
            JN2_GDA(IX,IY) = JN2_GDA(IX,IY) + BNTHDEN
      ENDIF
      if(itimesecs.ge.nxprtsed)    then
        iseddisk=1

c- output solid phase concentrations in mg/g (O2) mg/kg (po4) mg/g (Si)
        m2g=m2(ix,iy)
        hst2r=hst2av/1000./m2g
        po4t2r=po4t2av/1000./m2g
        sit2r=(psiav+sit2av)/1.0e6/m2g
        psiavr=(psiav)/1.0e6/m2g
c- output POC,PON,etc in mg/g
        pon1r=pon1av/1.0e6/m2g
        pon2r=pon2av/1.0e6/m2g
        pon3r=pon3av/1.0e6/m2g
        ponr=pon1r+pon2r+pon3r
        poc1r=poc1av/1.0e6/m2g
        poc2r=poc2av/1.0e6/m2g
        poc3r=poc3av/1.0e6/m2g
        pocr=poc1r+poc2r+poc3r
        pop1r=pop1av/1.0e6/m2g
        pop2r=pop2av/1.0e6/m2g
        pop3r=pop3av/1.0e6/m2g
        popr=pop1r+pop2r+pop3r
c- output dissolved concentrations in mg/L
        nh40z=nh40/1000.
        nh41z=nh41/1000.
        nh42z=nh42av/1000.
        nh4t2z=nh4t2av/1000.
        no30z=no30/1000.
        no31z=no31/1000.
        no32z=no32av/1000.
        no3t2z=no3t2av/1000.
        si0z=si0/1000.
        si1z=si1/1000.
        si2z=si2av/1000.
        sit2z=sit2av/1000.
        po40z=po40/1000.
        po41z=po41/1000.
        po42z=po42av/1000.
        po4t2z=po4t2av/1000.
        hs1z=hs1
        hs2z=hs2av
        hst2z=hst2av

        IF(IGDSEDOPT.EQ.0)  THEN
c        Instantaneous Dumps
         write(14)  time,ctemp(ix,iy),pop1r,pop2r,pop3r,popr
     .    ,pon1r,pon2r,pon3r,ponr,poc1r,poc2r,poc3r,pocr
     .    ,po4t2r,hst2r,sit2r,psiavr
     .    ,(flxpop(ix,iy,1)+flxpop(ix,iy,2)+flxpop(ix,iy,3))
     .    ,(flxpon(ix,iy,1)+flxpon(ix,iy,2)+flxpon(ix,iy,3))
     .    ,(flxpoc(ix,iy,1)+flxpoc(ix,iy,2)+flxpoc(ix,iy,3))
     .    ,o20,csod,sod,s,jp,jn,jc,jo2nh4,xjco2av,xjc1av
     .    ,jpo4,jnh4,jno3,jhs,jsi,jch4aq,jch4g,h1
     .    ,po40,po41,po42,po4t2,nh40,nh41,nh42,nh4t2
     .    ,no30,no31,no32,no3t2,hs1,hs2,hst2,si0,si1,si2,sit2
     .    ,bnthden,ch41,ch42,ch4t2,so41,so42,so4t2
     .    ,flxpop(ix,iy,1),flxpop(ix,iy,2),flxpop(ix,iy,3)
     .    ,flxpon(ix,iy,1),flxpon(ix,iy,2),flxpon(ix,iy,3)
     .    ,flxpoc(ix,iy,1),flxpoc(ix,iy,2),flxpoc(ix,iy,3)
     .    ,flxpos(ix,iy),flxpo4t2s(ix,iy),flxsit2s(ix,iy) 
         if(iseddisk.eq.1) then
            jpo4_gda(ix,iy)=jpo4
            jnh4_gda(ix,iy)=jnh4
            jno3_gda(ix,iy)=jno3
            jn2_gda(ix,iy)=bnthden
            jsi_gda(ix,iy)=jsi
            jh2s_gda(ix,iy)=jhs
            jch4aq_gda(ix,iy)=jch4aq
            jch4g_gda(ix,iy)=jch4g
            sod_gda(ix,iy)=sod
         endif
        ELSE
C        Global Dumps
           CTEMP_GDA(IX,IY) = CTEMP_GDA(IX,IY)/FLOAT(IGDCNT)
           POP1R_GDA(IX,IY) = POP1R_GDA(IX,IY)/FLOAT(IGDCNT)
           POP2R_GDA(IX,IY) = POP2R_GDA(IX,IY)/FLOAT(IGDCNT)
           POP3R_GDA(IX,IY) = POP3R_GDA(IX,IY)/FLOAT(IGDCNT)
           POPR_GDA(IX,IY)  = POPR_GDA(IX,IY)/FLOAT(IGDCNT)
           PON1R_GDA(IX,IY) = PON1R_GDA(IX,IY)/FLOAT(IGDCNT)
           PON2R_GDA(IX,IY) = PON2R_GDA(IX,IY)/FLOAT(IGDCNT)
           PON3R_GDA(IX,IY) = PON3R_GDA(IX,IY)/FLOAT(IGDCNT)
           PONR_GDA(IX,IY)  = PONR_GDA(IX,IY)/FLOAT(IGDCNT)
           POC1R_GDA(IX,IY) = POC1R_GDA(IX,IY)/FLOAT(IGDCNT)
           POC2R_GDA(IX,IY) = POC2R_GDA(IX,IY)/FLOAT(IGDCNT)
           POC3R_GDA(IX,IY) = POC3R_GDA(IX,IY)/FLOAT(IGDCNT)
           POCR_GDA(IX,IY)  = POCR_GDA(IX,IY)/FLOAT(IGDCNT)
           PO4T2R_GDA(IX,IY) = PO4T2R_GDA(IX,IY)/FLOAT(IGDCNT)
           HST2R_GDA(IX,IY) = HST2R_GDA(IX,IY)/FLOAT(IGDCNT)
           SIT2R_GDA(IX,IY) = SIT2R_GDA(IX,IY)/FLOAT(IGDCNT)
           PSIAVR_GDA(IX,IY) = PSIAVR_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOP_GDA(IX,IY) = FLXPOP_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOP1_GDA(IX,IY) = FLXPOP1_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOP2_GDA(IX,IY) = FLXPOP2_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOP3_GDA(IX,IY) = FLXPOP3_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPON_GDA(IX,IY) = FLXPON_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPON1_GDA(IX,IY) = FLXPON1_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPON2_GDA(IX,IY) = FLXPON2_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPON3_GDA(IX,IY) = FLXPON3_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOC_GDA(IX,IY) = FLXPOC_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOC1_GDA(IX,IY) = FLXPOC1_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOC2_GDA(IX,IY) = FLXPOC2_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOC3_GDA(IX,IY) = FLXPOC3_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPOS_GDA(IX,IY) = FLXPOS_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXPO4T2S_GDA(IX,IY) = FLXPO4T2S_GDA(IX,IY)/FLOAT(IGDCNT)
           FLXSIT2S_GDA(IX,IY) = FLXSIT2S_GDA(IX,IY)/FLOAT(IGDCNT)
           O20_GDA(IX,IY)   = O20_GDA(IX,IY)/FLOAT(IGDCNT)
           CSOD_GDA(IX,IY)  = CSOD_GDA(IX,IY)/FLOAT(IGDCNT)
           SOD_GDA(IX,IY)   = SOD_GDA(IX,IY)/FLOAT(IGDCNT)
           S_GDA(IX,IY)     = S_GDA(IX,IY)/FLOAT(IGDCNT)
           JP_GDA(IX,IY)    = JP_GDA(IX,IY)/FLOAT(IGDCNT)
           JN_GDA(IX,IY)    = JN_GDA(IX,IY)/FLOAT(IGDCNT)
           JC_GDA(IX,IY)    = JC_GDA(IX,IY)/FLOAT(IGDCNT)
           JO2NH4_GDA(IX,IY) = JO2NH4_GDA(IX,IY)/FLOAT(IGDCNT)
           XJCO2AV_GDA(IX,IY) = XJCO2AV_GDA(IX,IY)/FLOAT(IGDCNT)
           XJC1AV_GDA(IX,IY) = XJC1AV_GDA(IX,IY)/FLOAT(IGDCNT)
           JPO4_GDA(IX,IY)  = JPO4_GDA(IX,IY)/FLOAT(IGDCNT)
           JNH4_GDA(IX,IY)  = JNH4_GDA(IX,IY)/FLOAT(IGDCNT)
           JNO3_GDA(IX,IY)  = JNO3_GDA(IX,IY)/FLOAT(IGDCNT)
           JH2S_GDA(IX,IY)   = JH2S_GDA(IX,IY)/FLOAT(IGDCNT)
           JSI_GDA(IX,IY)   = JSI_GDA(IX,IY)/FLOAT(IGDCNT)
           JCH4AQ_GDA(IX,IY) = JCH4AQ_GDA(IX,IY)/FLOAT(IGDCNT)
           JCH4G_GDA(IX,IY) = JCH4G_GDA(IX,IY)/FLOAT(IGDCNT)
           H1_GDA(IX,IY)    = H1_GDA(IX,IY) /FLOAT(IGDCNT)
           PO40_GDA(IX,IY)  = PO40_GDA(IX,IY)/FLOAT(IGDCNT)
           PO41_GDA(IX,IY)  = PO41_GDA(IX,IY)/FLOAT(IGDCNT)
           PO42_GDA(IX,IY)  = PO42_GDA(IX,IY)/FLOAT(IGDCNT)
           PO4T2_GDA(IX,IY) = PO4T2_GDA(IX,IY)/FLOAT(IGDCNT)
           NH40_GDA(IX,IY)  = NH40_GDA(IX,IY)/FLOAT(IGDCNT)
           NH41_GDA(IX,IY)  = NH41_GDA(IX,IY)/FLOAT(IGDCNT)
           NH42_GDA(IX,IY)  = NH42_GDA(IX,IY)/FLOAT(IGDCNT)
           NH4T2_GDA(IX,IY) = NH4T2_GDA(IX,IY)/FLOAT(IGDCNT)
           NO30_GDA(IX,IY)  = NO30_GDA(IX,IY)/FLOAT(IGDCNT)
           NO31_GDA(IX,IY)  = NO31_GDA(IX,IY)/FLOAT(IGDCNT)
           NO32_GDA(IX,IY)  = NO32_GDA(IX,IY)/FLOAT(IGDCNT)
           NO3T2_GDA(IX,IY) = NO3T2_GDA(IX,IY)/FLOAT(IGDCNT)
           HS1_GDA(IX,IY)   = HS1_GDA(IX,IY)/FLOAT(IGDCNT)
           HS2_GDA(IX,IY)   = HS2_GDA(IX,IY)/FLOAT(IGDCNT)
           HST2_GDA(IX,IY)  = HST2_GDA(IX,IY)/FLOAT(IGDCNT)
           SI0_GDA(IX,IY)   = SI0_GDA(IX,IY)/FLOAT(IGDCNT)
           SI1_GDA(IX,IY)   = SI1_GDA(IX,IY)/FLOAT(IGDCNT)
           SI2_GDA(IX,IY)   = SI2_GDA(IX,IY)/FLOAT(IGDCNT)
           SIT2_GDA(IX,IY)  = SIT2_GDA(IX,IY)/FLOAT(IGDCNT)
           JN2_GDA(IX,IY) = JN2_GDA(IX,IY)/FLOAT(IGDCNT)
           CH41_GDA(IX,IY) = CH41_GDA(IX,IY)/FLOAT(IGDCNT)
           CH42_GDA(IX,IY) = CH42_GDA(IX,IY)/FLOAT(IGDCNT)
           CH4T2_GDA(IX,IY) = CH4T2_GDA(IX,IY)/FLOAT(IGDCNT)
           SO41_GDA(IX,IY) = SO41_GDA(IX,IY)/FLOAT(IGDCNT)
           SO42_GDA(IX,IY) = SO42_GDA(IX,IY)/FLOAT(IGDCNT)
           SO4T2_GDA(IX,IY) = SO4T2_GDA(IX,IY)/FLOAT(IGDCNT)
          WRITE(14)  TIME,CTEMP_GDA(IX,IY),POP1R_GDA(IX,IY),
     .     POP2R_GDA(IX,IY),POP3R_GDA(IX,IY),POPR_GDA(IX,IY),
     .     PON1R_GDA(IX,IY),PON2R_GDA(IX,IY),PON3R_GDA(IX,IY),
     .     PONR_GDA(IX,IY),POC1R_GDA(IX,IY),POC2R_GDA(IX,IY),
     .     POC3R_GDA(IX,IY),POCR_GDA(IX,IY),PO4T2R_GDA(IX,IY),
     .     HST2R_GDA(IX,IY),SIT2R_GDA(IX,IY),PSIAVR_GDA(IX,IY),
     .     FLXPOP_GDA(IX,IY),FLXPON_GDA(IX,IY),FLXPOC_GDA(IX,IY),
     .     O20_GDA(IX,IY),CSOD_GDA(IX,IY),SOD_GDA(IX,IY),S_GDA(IX,IY),
     .     JP_GDA(IX,IY),JN_GDA(IX,IY),JC_GDA(IX,IY),JO2NH4_GDA(IX,IY),
     .     XJCO2AV_GDA(IX,IY),XJC1AV_GDA(IX,IY),JPO4_GDA(IX,IY),
     .     JNH4_GDA(IX,IY),JNO3_GDA(IX,IY),JH2S_GDA(IX,IY),
     .     JSI_GDA(IX,IY),JCH4AQ_GDA(IX,IY),JCH4G_GDA(IX,IY),
     .     H1_GDA(IX,IY),PO40_GDA(IX,IY),PO41_GDA(IX,IY),
     .     PO42_GDA(IX,IY),PO4T2_GDA(IX,IY),NH40_GDA(IX,IY),
     .     NH41_GDA(IX,IY),NH42_GDA(IX,IY),NH4T2_GDA(IX,IY),
     .     NO30_GDA(IX,IY),NO31_GDA(IX,IY),NO32_GDA(IX,IY),
     .     NO3T2_GDA(IX,IY),HS1_GDA(IX,IY),HS2_GDA(IX,IY),
     .     HST2_GDA(IX,IY),SI0_GDA(IX,IY),SI1_GDA(IX,IY),
     .     SI2_GDA(IX,IY),SIT2_GDA(IX,IY),JN2_GDA(IX,IY),
     .     CH41_GDA(IX,IY),CH42_GDA(IX,IY),CH4T2_GDA(IX,IY),
     .     SO41_GDA(IX,IY),SO42_GDA(IX,IY),SO4T2_GDA(IX,IY),
     .     FLXPOP1_GDA(IX,IY),FLXPOP2_GDA(IX,IY),FLXPOP3_GDA(IX,IY),
     .     FLXPON1_GDA(IX,IY),FLXPON2_GDA(IX,IY),FLXPON3_GDA(IX,IY),
     .     FLXPOC1_GDA(IX,IY),FLXPOC2_GDA(IX,IY),FLXPOC3_GDA(IX,IY),
     .     FLXPOS_GDA(IX,IY),FLXPO4T2S_GDA(IX,IY),FLXSIT2S_GDA(IX,IY)
          ENDIF

      endif    ! endif for itimesecs .ge. nxprtsed
C
      call storesed

  200 CONTINUE
      IF(IGDSEDOPT.EQ.1) IGDCNT=IGDCNT+1

C             TAKE TEMPERATURE INTEGRATION STEP
      DO 300 IY=1,NY
      DO 300 IX=1,NX
      IF(FSM(IX,IY).NE.1.)  GO TO 300
         CTEMP(IX,IY) = CTEMP(IX,IY)+DELTAT*DIFFT/HSED(IX,IY)/
     .      HSED(IX,IY)*(HYDTEMP(IX,IY,NZ)-CTEMP(IX,IY))
C        CTEMP(IX,IY)=MIN(CTEMP(IX,IY),25.)
  300 CONTINUE

      IF(MASSBAL.EQ.1)  THEN
        TOTORGP=0.
        TOTORGN=0.
        TOTORGC=0.
        TOTBSI=0.
        TOTPO4=0.
        TOTNH4=0.
        TOTNO3=0.
        TOTSI=0.
        TOTH2S=0.
        TOTCH4=0.
        TOTSO4=0.
        DO 320 IG=1,3
         DO 320 IY=1,NY
          DO 320 IX=1,NX
           IF(FSM(IX,IY).EQ.1.)  THEN
            TOTORGP=TOTORGP+CPOP(IX,IY,IG)*BSVOL(IX,IY)
            TOTORGN=TOTORGN+CPON(IX,IY,IG)*BSVOL(IX,IY)
            TOTORGC=TOTORGC+CPOC(IX,IY,IG)*BSVOL(IX,IY)
           ENDIF
  320   CONTINUE
        DO 325 IY=1,NY
         DO 325 IX=1,NX
          IF(FSM(IX,IY).EQ.1.)  THEN
           H1=H1S(IX,IY)
           TOTBSI=TOTBSI+CPOS(IX,IY)*BSVOL(IX,IY)
           TOTPO4=TOTPO4+XAZ(IX,IY)*
     .                   (H1*PO4T1TM1S(IX,IY)+(H2-H1)*PO4T2TM1S(IX,IY))
           TOTNH4=TOTNH4+XAZ(IX,IY)*
     .                   (H1*NH4T1TM1S(IX,IY)+(H2-H1)*NH4T2TM1S(IX,IY))
           TOTNO3=TOTNO3+XAZ(IX,IY)*
     .                   (H1*NO3T1TM1S(IX,IY)+(H2-H1)*NO3T2TM1S(IX,IY))
           TOTSI=TOTSI+XAZ(IX,IY)*
     .                   (H1*SIT1TM1S(IX,IY)+(H2-H1)*SIT2TM1S(IX,IY))
           TOTH2S=TOTH2S+XAZ(IX,IY)*
     .                   (H1*HST1TM1S(IX,IY)+(H2-H1)*HST2TM1S(IX,IY))
           TOTCH4=TOTCH4+XAZ(IX,IY)*
     .                   (H1*PO41TM1S(IX,IY)+(H2-H1)*PO4T2TM1S(IX,IY))
           TOTSO4=TOTSO4+XAZ(IX,IY)*
     .                   (H1*SO41TM1S(IX,IY)+(H2-H1)*SO4T2TM1S(IX,IY))
          ENDIF
  325   CONTINUE
        JPOP=0.
        JPON=0.
        JPOC=0.
        JBSI=0.
        JPO4SS=0.
        JSISS=0.
        BURIALORGP=0.
        BURIALORGN=0.
        BURIALORGC=0.
        BURIALBSI=0.
        BURIALPO4=0.
        BURIALNH4=0.
        BURIALNO3=0.
        BURIALSI=0.
        BURIALH2S=0.
        BURIALCH4=0.
        BURIALSO4=0.
        DO 330 IG=1,3
         DO 330 IY=1,NY                !ldh-1
          DO 330 IX=1,NX               !ldh-1
           IF(FSM(IX,IY).EQ.1.)  THEN
            BURIALORGP=BURIALORGP+VSED(IX,IY)*CPOP(IX,IY,IG)
            BURIALORGN=BURIALORGN+VSED(IX,IY)*CPON(IX,IY,IG)
            BURIALORGC=BURIALORGC+VSED(IX,IY)*CPOC(IX,IY,IG)
           ENDIF
  330   CONTINUE
        BURIALBSI=BURIALBSI+VSED(IX,IY)*CPOS(IX,IY)
        BURIALPO4=BURIALPO4+VSED(IX,IY)*PO4T2TM1S(IX,IY)
        BURIALNH4=BURIALNH4+VSED(IX,IY)*NH4T2TM1S(IX,IY)
        BURIALNO3=BURIALNO3+VSED(IX,IY)*NO3T2TM1S(IX,IY)
        BURIALSI=BURIALSI+VSED(IX,IY)*SIT2TM1S(IX,IY)
        BURIALH2S=BURIALH2S+VSED(IX,IY)*HST2TM1S(IX,IY)
        BURIALCH4=BURIALCH4+VSED(IX,IY)*CH4T2TM1S(IX,IY)
        BURIALSO4=BURIALSO4+VSED(IX,IY)*SO4T2TM1S(IX,IY)
        DO 340 IY=1,NY
         DO 340 IX=1,NX
          IF(FSM(IX,IY).EQ.1.)  THEN
           JPOP=JPOP+FLXPOP(IX,IY,1)+FLXPOP(IX,IY,2)+FLXPOP(IX,IY,3)
           JPON=JPON+FLXPON(IX,IY,1)+FLXPON(IX,IY,2)+FLXPON(IX,IY,3)
           JPOC=JPOC+FLXPOC(IX,IY,1)+FLXPOC(IX,IY,2)+FLXPOC(IX,IY,3)
           JBSI=JBSI+FLXPOS(IX,IY)
           JPO4SS=JPO4SS+FLXPO4T2S(IX,IY)
           JSISS=JSISS+FLXSIT2S(IX,IY)
           IF(IGDSEDOPT.EQ.0)  THEN
            TOTJPO4=TOTJPO4+JPO4_GDA(IX,IY)
            TOTJNH4=TOTJNH4+JNH4_GDA(IX,IY)
            TOTJNO3=TOTJNO3+JNO3_GDA(IX,IY)
            TOTJN2=TOTJN2+JN2_GDA(IX,IY)
            TOTJSI=TOTJSI+JSI_GDA(IX,IY)
            TOTJH2S=TOTJH2S+JH2S_GDA(IX,IY)
            TOTJCH4AQ=TOTJCH4AQ+JCH4AQ_GDA(IX,IY)
            TOTJCH4G=TOTJCH4G+JCH4G_GDA(IX,IY)
            TOTJSOD=TOTJSOD+SOD_GDA(IX,IY)
           ELSE
            TOTJPO4=TOTJPO4+JPO4_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJNH4=TOTJNH4+JNH4_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJNO3=TOTJNO3+JNO3_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJN2=TOTJN2+JN2_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJSI=TOTJSI+JSI_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJH2S=TOTJH2S+JH2S_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJCH4AQ=TOTJCH4AQ+JCH4AQ_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJCH4G=TOTJCH4G+JCH4G_GDA(IX,IY)/FLOAT(IGDCNT)
            TOTJSOD=TOTJSOD+SOD_GDA(IX,IY)/FLOAT(IGDCNT)
           ENDIF
          ENDIF
  340   CONTINUE
        WRITE(18)  TIME,TOTORGP,TOTORGN,TOTORGC,TOTBSI,TOTPO4,TOTNH4
     .            ,TOTNO3,TOTSI,TOTH2S,TOTCH4,TOTSO4
     .            ,JPOP,JPON,JPOC,JBSI,JPO4SS,JSISS
     .            ,BURIALORGP,BURIALORGN,BURIALORGC
     .            ,BURIALBSI,BURIALPO4,BURIALNH4,BURIALNO3,BURIALSI
     .            ,BURIALH2S,BURIALCH4,BURIALSO4,TOTJPO4,TOTJNH4
     .            ,TOTJNO3,TOTJN2,TOTJSI,TOTJH2S,TOTJCH4AQ,TOTJCH4G
     .            ,TOTSOD
      ENDIF

      IF(ISEDDISK.EQ.1) THEN
       IF(IGDSEDOPT.EQ.1) THEN
        IGDCNT=0
        DO IY = 1,NY
          DO IX = 1,NX
            CTEMP_GDA(IX,IY) = 0.
            POP1R_GDA(IX,IY) = 0.
            POP2R_GDA(IX,IY) = 0.
            POP3R_GDA(IX,IY) = 0.
            POPR_GDA(IX,IY) = 0.
            PON1R_GDA(IX,IY) = 0.
            PON2R_GDA(IX,IY) = 0.
            PON3R_GDA(IX,IY) = 0.
            PONR_GDA(IX,IY) = 0.
            POC1R_GDA(IX,IY) = 0.
            POC2R_GDA(IX,IY) = 0.
            POC3R_GDA(IX,IY) = 0.
            POCR_GDA(IX,IY) = 0.
            PO4T2R_GDA(IX,IY) = 0.
            HST2R_GDA(IX,IY) = 0.
            SIT2R_GDA(IX,IY) = 0.
            PSIAVR_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            FLXPOP_GDA(IX,IY) = 0.
            FLXPOP1_GDA(IX,IY) = 0.
            FLXPOP2_GDA(IX,IY) = 0.
            FLXPOP3_GDA(IX,IY) = 0.
            FLXPON_GDA(IX,IY) = 0.
            FLXPON1_GDA(IX,IY) = 0.
            FLXPON2_GDA(IX,IY) = 0.
            FLXPON3_GDA(IX,IY) = 0.
            FLXPOC_GDA(IX,IY) = 0.
            FLXPOC1_GDA(IX,IY) = 0.
            FLXPOC2_GDA(IX,IY) = 0.
            FLXPOC3_GDA(IX,IY) = 0.
            FLXPOS_GDA(IX,IY) = 0.
            FLXPO4T2S_GDA(IX,IY) = 0.
            FLXSIT2S_GDA(IX,IY) = 0.
            O20_GDA(IX,IY) = 0.
            CSOD_GDA(IX,IY) = 0.
            SOD_GDA(IX,IY) = 0.
            S_GDA(IX,IY) = 0.
            JP_GDA(IX,IY) = 0.
            JN_GDA(IX,IY) = 0.
            JC_GDA(IX,IY) = 0.
            JO2NH4_GDA(IX,IY) = 0.
            XJCO2AV_GDA(IX,IY) = 0.
            XJC1AV_GDA(IX,IY) = 0.
            JPO4_GDA(IX,IY) = 0.
            JNH4_GDA(IX,IY) = 0.
            JNO3_GDA(IX,IY) = 0.
            JH2S_GDA(IX,IY) = 0.
            JSI_GDA(IX,IY) = 0.
            JCH4AQ_GDA(IX,IY) = 0.
            JCH4G_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            H1_GDA(IX,IY) = 0.
            PO40_GDA(IX,IY) = 0.
            PO41_GDA(IX,IY) = 0.
            PO42_GDA(IX,IY) = 0.
            PO4T2_GDA(IX,IY) = 0.
            NH40_GDA(IX,IY) = 0.
            NH41_GDA(IX,IY) = 0.
            NH42_GDA(IX,IY) = 0.
            NH4T2_GDA(IX,IY) = 0.
            NO30_GDA(IX,IY) = 0.
            NO31_GDA(IX,IY) = 0.
            NO32_GDA(IX,IY) = 0.
            NO3T2_GDA(IX,IY) = 0.
            HS1_GDA(IX,IY) = 0.
            HS2_GDA(IX,IY) = 0.
            HST2_GDA(IX,IY) = 0.
            SI0_GDA(IX,IY) = 0.
            SI1_GDA(IX,IY) = 0.
            SI2_GDA(IX,IY) = 0.
            SIT2_GDA(IX,IY) = 0.
          ENDDO
          DO IX = 1,NX
            JN2_GDA(IX,IY) = 0.
            CH41_GDA(IX,IY) = 0.
            CH42_GDA(IX,IY) = 0.
            CH4T2_GDA(IX,IY) = 0.
            SO41_GDA(IX,IY) = 0.
            SO42_GDA(IX,IY) = 0.
            SO4T2_GDA(IX,IY) = 0.
          ENDDO
        ENDDO
       ENDIF
       iseddisk=0
      ENDIF

      if(idisk.ge.1)    then
       rewind(16)
       write(16)   CTEMP,CPOP,CPON,CPOC,CPOS
     .        ,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S
     .        ,HST2TM1S,SIT2TM1S,BNTHSTR1S
     .        ,PO4T1TM1S,NH4T1TM1S,NO3T1TM1S
     .        ,SIT1TM1S,HST1TM1S
     .        ,CH4T1TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S
     .        ,PO41TM1S,NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S
     .        ,O20TM1S,SODTM1S,DD0TM1S
     .        ,BFORMAXS,ISWBNTHS
     .        ,DICT1TM1S,DICT2TM1S,TAT1TM1S,TAT2TM1S
     .        ,DIC1TM1S,TA1TM1S  
      endif

      if(itimesecs.ge.nxprtsed)  nxprtsed=nxprtsed+iprntsedsecs
      IF(ITIMESECS.EQ.OUTPUTNUM*86400)THEN
      WRITE(SIC,'(I2)')IRECORD 
      OPEN(1616,FILE=TRIM('RESTART_RCAFICSED')//'_'//TRIM(ADJUSTL(SIC))
     . ,FORM='UNFORMATTED')      
        WRITE(1616)   CTEMP,CPOP,CPON,CPOC,CPOS
     .        ,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S
     .        ,HST2TM1S,SIT2TM1S,BNTHSTR1S
     .        ,PO4T1TM1S,NH4T1TM1S,NO3T1TM1S
     .        ,SIT1TM1S,HST1TM1S
     .        ,CH4T1TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S
     .        ,PO41TM1S,NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S
     .        ,O20TM1S,SODTM1S,DD0TM1S
     .        ,BFORMAXS,ISWBNTHS
     .        ,DICT1TM1S,DICT2TM1S,TAT1TM1S,TAT2TM1S
     .        ,DIC1TM1S,TA1TM1S
C by lb 20171019
      CLOSE(1616)
      END IF
      RETURN

  900 WRITE(OUT,9005)
 9005 FORMAT(/' READ ERROR IN SEDIMENT INPUT DECK')
      STOP
      END
