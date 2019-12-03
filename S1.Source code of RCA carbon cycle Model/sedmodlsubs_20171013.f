cOPTION CROSSREF ON
cOPTION SYMTABLE ON
cPAGEWIDTH 132
      SUBROUTINE DIAGEN
c-- Diagenesis conputation
c-    frpon1 = G1 fraction of jnin
c-    frpon2 = G2 fraction of jnin
c-    frpon3 = G3 fraction of jnin
c-    pon1(nl) = G1 pon concentration
c-    pon2(nl) = G2 pon concentration
c-    pon3(nl) = G3 pon concentration
c-    kpon1 = G1 reaction rate
c-    kpon2 = G2 reaction rate
c-    kpon3 = G3 reaction rate
c-    thtapon1 = theta for kpon1
c-    thtapon2 = theta for kpon2
c-    thtapon3 = theta for kpon3
c-- Similarly for poc and pop
c-    psi(nl) = biogenic si
c-    ksi = specific reaction rate for dissolution
c-    thtaksi = theta for ksi

      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON

      real
     .     kpop1  , kpop2  , kpop3
     .   , kpon1  , kpon2  , kpon3
     .   , kpoc1  , kpoc2  , kpoc3
      equivalence
     .     (kpdiag(1),kpop1) , (kpdiag(2),kpop2) , (kpdiag(3),kpop3)
     .   , (kndiag(1),kpon1) , (kndiag(2),kpon2) , (kndiag(3),kpon3)
     .   , (kcdiag(1),kpoc1) , (kcdiag(2),kpoc2) , (kcdiag(3),kpoc3)
     .  , (dpthta(1),thtapop1),(dpthta(2),thtapop2),(dpthta(3),thtapop3)
     .  , (dnthta(1),thtapon1),(dnthta(2),thtapon2),(dnthta(3),thtapon3)
     .  , (dcthta(1),thtapoc1),(dcthta(2),thtapoc2),(dcthta(3),thtapoc3)
      itemp = 10.*tempd + 1
      xkpop1=zhtapop1(itemp)*deltaz
      xkpop2=zhtapop2(itemp)*deltaz
      xkpop3=zhtapop3(itemp)*deltaz
      xkpon1=zhtapon1(itemp)*deltaz
      xkpon2=zhtapon2(itemp)*deltaz
      xkpon3=zhtapon3(itemp)*deltaz
      xkpoc1=zhtapoc1(itemp)*deltaz
      xkpoc2=zhtapoc2(itemp)*deltaz
      xkpoc3=zhtapoc3(itemp)*deltaz
      xksi=zhtasi(itemp)*deltaz
      xkjsi=zhtajsi(itemp)

c- pon1
      pon1=(flxpon(IX,IY,1)*deltat/h20+pon1tm1)
     .       /(1.0+(xkpon1+w2)*deltat/h20)
      pon1av=pon1

c- pon2
      pon2=(flxpon(IX,IY,2)*deltat/h20+pon2tm1)
     .       /(1.0+(xkpon2+w2)*deltat/h20)
      pon2av=pon2

c- pon3
      pon3=(flxpon(IX,IY,3)*deltat/h20+pon3tm1)
     .       /(1.0+(xkpon3+w2)*deltat/h20)
      pon3av=pon3

      iprt=0

c- jn
      xjn=xkpon1*pon1+xkpon2*pon2+xkpon3*pon3
      jn=xjn

c- poc1
      poc1=(flxpoc(IX,IY,1)*deltat/h20+poc1tm1)
     .       /(1.0+(xkpoc1+w2)*deltat/h20)
      poc1av=poc1

c- poc2
      poc2=(flxpoc(IX,IY,2)*deltat/h20+poc2tm1)
     .       /(1.0+(xkpoc2+w2)*deltat/h20)
      poc2av=poc2

c- poc3
      poc3=(flxpoc(IX,IY,3)*deltat/h20+poc3tm1)
     .       /(1.0+(xkpoc3+w2)*deltat/h20)
      poc3av=poc3
      iprt=0

c- jc
      xjc=xkpoc1*poc1+xkpoc2*poc2+xkpoc3*poc3
      jc=xjc

c- pop1
      pop1=(flxpop(IX,IY,1)*deltat/h20+pop1tm1)
     .       /(1.0+(xkpop1+w2)*deltat/h20)
      pop1av=pop1

c- pop2
      pop2=(flxpop(IX,IY,2)*deltat/h20+pop2tm1)
     .       /(1.0+(xkpop2+w2)*deltat/h20)
      pop2av=pop2

c- pop3
      pop3=(flxpop(IX,IY,3)*deltat/h20+pop3tm1)
     .       /(1.0+(xkpop3+w2)*deltat/h20)
      pop3av=pop3

      iprt=0

c- jp
      xjp=xkpop1*pop1+xkpop2*pop2+xkpop3*pop3
      jp=xjp

c- psi
c     Reaction kinetics
c     d psi/ dt  = -xksi*psi/(psi+kmpsi)*(csisat-si2) + jsiin + jsidetr
c
c     Two sources: jsiin computed from C/Si ratio, and jsidetr, detritial
c     source of Si.

      fd2=1.0/(1.0+m2(ix,iy)*pie2si(ix,iy))
      k3=xksi*(csisat-fd2*sit2tm1)
     .       /(psitm1+kmpsi)
      psi=((flxpos(IX,IY)+jsidetr)*deltat/h20+psitm1)
     .       /(1.0+(k3+w2)*deltat/h20)

      iprt=0
d     if(idisk.ge.1)  iprt=0
d     ixp=4
d     iyp=6

d       if(iprt.eq.1.and.IX.eq.ixp.and.IY.eq.iyp)then
d         print *,'DIAGEN ix,iy',ix,iy
d         print *,'    '
d         print *,'psi,flxpos,deltat,h2,psitm1,k3,w2'
d         print *, psi,flxpos(ix,iy),deltat,h2,psitm1,k3,w2 
d         print *,'    '
d       endif

      psiav=psi
      return
      end

      FUNCTION SEDF(SOD1,IUSE)
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
      REAL*8 AD(4,4), BX(4), G(2), HX(2,2)                       !CH4-SO4 Code
      REAL*8 DBLSO41, DBLSO42, RA0, RA1, RA2, R1, R2, DISC, SN1  !CH4-SO4 Code
cLIST ON

      iprt=0
d     if(idisk.ge.1)  iprt=0
d     ixp=4
d     iyp=6
d     if(iprt.eq.1.and.ix.eq.ixp.and.iy.eq.iyp .and. iuse.eq.1)  
d    .   print *,'  SOD iteration'
d     if(iprt.eq.1.and.ix.eq.ixp.and.iy.eq.iyp .and. iuse.ne.1) 
d    .   print *,'  FLUX iteration'

c-- if iuse=1 then this is an SOD iteration
c   if iuse<>1 then this is a full evaluation of all the fluxes

      if(iuse.eq.0) then
c-- Temperature corrections
      itemp = 10.*tempd + 1
      xkmnh4=zhtakmnh4(itemp)
      xappd1=zhtad1(itemp)
      xappp1=zhtap1(itemp)

      if(sal0.ge.saltsw)  then
         xappnh4=zhtanh4s(itemp)
         xapp1no3=zhtano3s(itemp)
         xk2no3=zhta2no3s(itemp)*DELTAZ
      else
         xappnh4=zhtanh4f(itemp)
         xapp1no3=zhtano3f(itemp)
         xk2no3=zhta2no3f(itemp)*DELTAZ
      endif

      xksi=zhtasi(itemp)*deltaz
      xDd0=zhtaDd0(itemp)
      xappch4=zhtach4(itemp)
c- layer 1-2 transport
      kl12nom=Dd/deltaz*zl12nom(itemp)
c- include organic carbon term to simulate benthic biomass
      w12nom=Dp/deltaz*zw12nom(itemp)*poc1av/1.0e5
c-- Benthic stress equation
c--           iswbnth=0 -> low temperature period
c--           iswbnth=1 -> high temperature period
      if(iswbnth.eq.0) then
c-- low temp
              if(tempd.ge.tempbnth) then
c-- change to high temp
              iswbnth=1
              bformax=0.0
              endif
         bfor=kmo2Dp/(kmo2Dp+o20)
         else
c-- high temp
              if(tempd.lt.tempbnth) then
c-- change to low temp
              iswbnth=0
              endif
         bformax=amax1(kmo2Dp/(kmo2Dp+o20),bformax)
         bfor=bformax
         endif

         bnthstr=(bnthstr1+deltat*bfor)
     .           /(1.0+kbnthstr*deltat)
         bnth=(1.0-kbnthstr*bnthstr)
c-- w12min= Dpmin/deltaz is physical mixing
         w12min= Dpmin/deltaz
         w12=w12nom*bnth + w12min

c*-add in the enhancement due to benthic activity
         kl12=kl12nom + klbnth*w12

c-*****
      return
      endif


c-- Entry for computing the NH4,NO3, and SOD fluxes
      if(iuse.eq.1) then

c-- surface mass transfer
      s = sod1/o20

c-- h1dot evaluation
      h1=xDd0*o20/sod1
c   limit to the full depth
      h1=amin1(h1,h20)
      h2=h20-h1
c        prevent h2 from going to 0
      h2=amax1(h2,(0.001*h1))
      h1=h20-h2
      Dd0o20dot=(xDd0*o20-Dd0tm1*o20tm1)/deltat
      h1dot=(-h1*(sod1-sodtm1)/deltat + Dd0o20dot)/sod1
      h2dot=-h1dot
      DEW(IX,IY) = H20
      DE1(IX,IY) = H1

c- ammonia flux

      k0h1p=0.
      k1h1p=0.
      k2h2d=0.
      k2h2p=0.
      if(kmnh4.ne.0.0) then
         k0h1d=xappnh4**2/s*xkmnh4*(o20/(kmnh4o2+o20))
         k1h1d=s
      else
         k1h1d=xappnh4**2/s*(o20/(kmnh4o2+o20))+s
         k0h1d=0.
      endif
      j1=s*nh40
      k3=0.0
      j2=jn
      pie1=pienh4
      pie2=pienh4
      kmc1=xkmnh4
d       if(iprt.eq.1.and.IX.eq.ixp.and.IY.eq.iyp)then
d          print *,' k0h1d,k1h1d,j1,j2,pie1,pie2,kmc1'
d          print *,  k0h1d,k1h1d,j1,j2,pie1,pie2,kmc1 
d       endif
      call sedtsfnl
     .(nh41,nh42,nh42av,nh4t1,nh4t2,nh4t2av,nh41tm1,nh4t1tm1,nh4t2tm1,1)
      jnh4=s*(nh41-nh40)
d       if(iprt.eq.1.and.IX.eq.ixp.and.IY.eq.iyp)then
d         print *,' nh41,nh42,nh42av'
d         print *,  nh41,nh42,nh42av
d         print *,' nh4t1,nh4t2,nh4t2av'
d         print *,  nh4t1,nh4t2,nh4t2av 
d         print *,' nh41tm1,nh4t2tm1'
d         print *,  nh41tm1,nh4t2tm1 
d         print *,' jnh4=s*(nh41-nh40)'
d         print *,  jnh4,s,nh41,nh40
d       endif

c- Oxygen consumed by nitrification
c        a1 = 64/14 * 1/1000 - mole ratio and mg/m2-day to gm/m2-day
         a1 = 0.004571429
      if(kmnh4.ne.0.0) then
         jo2nh4=a1*k0h1d*nh41/(xkmnh4+nh41tm1)
      else
         jo2nh4=a1*(k1h1d-s)*nh41
      endif
c- DIC consumed by nitrification
c        a1dic  mmol C/m2-day
         a1dic = 1./14.*(-1.98+1.88)
      if(kmnh4.ne.0.0) then
         jdicnh4=a1dic*k0h1d*nh41/(xkmnh4+nh41tm1)
      else
         jdicnh4=a1dic*(k1h1d-s)*nh41
      endif
      NIDIC(IX,IY,1) = JDICNH4
      NIDIC(IX,IY,2) = 0.0
c- Talk consumed by nitrification
c        a1ta  mmol C/m2-day
         a1ta = 1./14.*(-1.98)
      if(kmnh4.ne.0.0) then
         jtanh4=a1ta*k0h1d*nh41/(xkmnh4+nh41tm1)
      else
         jtanh4=a1ta*(k1h1d-s)*nh41
      endif
      NITA(IX,IY,1) = JTANH4
      NITA(IX,IY,2) = 0.0

c- denitrification

      k0h1d=0.
      k0h1p=0.
      kmc1=0.0
      k1h1d=xapp1no3**2/s+s
      k1h1p=0.
      k2h2d=xk2no3
      k2h2p=0.
      if(xkmnh4.ne.0.0) then
        j1=s*no30+
     .   xappnh4**2/s*xkmnh4*(o20/(kmnh4o2+o20))*nh41/(xkmnh4+nh41tm1)
      else
         j1=s*no30+
     .      xappnh4**2/s*(o20/(kmnh4o2+o20))*nh41
      endif
      k3=0.0
      j2=0.0
      pie1=0.
      pie2=0.
      call sedtsfnl
     .(no31,no32,no32av,no3t1,no3t2,no3t2av,no31tm1,no3t1tm1,no3t2tm1,1)
      jno3=s*(no31-no30)

c- Sulfide/methane oxidation
c        diagenesis consumed by denitrification
c        a2 = 10/8 * 32/14 * 1/1000
         a2 = 0.00285714
         xjcno31 = a2*xapp1no3**2/s*no31
         xjcno3 = a2*xk2no3*no32
c        DIC consumed by denitrification
         a2dic = 10./8./14.
         jdicno31 = a2dic*xapp1no3**2/s*no31
         jdicno3 = a2dic*xk2no3*no32
         DNDIC(IX,IY,1) = -1.0*JDICNO31
         DNDIC(IX,IY,2) = -1.0*JDICNO3
c        Talk consumed by denitrification
         a2ta = 10./8.*1./14.*0.8
         jtano31 = a2ta*xapp1no3**2/s*no31
         jtano3 = a2ta*xk2no3*no32
         DNTA(IX,IY,1) = -1.0*JTANO31
         DNTA(IX,IY,2) = -1.0*JTANO3
c- add the aerobic and first anaerobic layer to keep mass balance
         xjcno3=xjcno31+xjcno3

d       if(iprt.eq.1.and.IX.eq.ixp.and.IY.eq.iyp)then
d         print *,'ixp,iyp,xjcno31,a2,xapp1no3,s,no31'
d         print *, ixp,iyp,xjcno31,a2,xapp1no3,s,no31 
d         print *,'ixp,iyp,xjcno3 ,a2,xk2no3,no32'
d         print *, ixp,iyp,xjcno3 ,a2,xk2no3,no32 
d       endif

c- Convert carbon diagenesis flux to o2*  units and decrement jco2

c        a0= 32./12./1000.
         a0=2.666666666e-3
         xjco2av=0.0
         xjc1av=0.0
         xjco2 = a0*jc
         xjc1 = amax1(xjco2-xjcno3,0.0)
         xjco2av=xjco2av+xjco2
         xjc1av=xjc1av+xjc1

d       if(iprt.eq.1.and.IX.eq.ixp.and.IY.eq.iyp)then
d         print *,'ixp,iyp',ixp,iyp
d         print *,'xjco2,a0,xjc,xjc1,xjcno3,xjco2av,xjco2av,xjc1av'
d         print *, xjco2,a0,xjc,xjc1,xjcno3,xjco2av,xjco2av,xjc1av 
d       endif

C**** **********************************************************
C**** New code for methane formation.  CH4 starts forming
C**** once all sulfate is used up.
C**** **********************************************************

C**** Sulfide and sulfate in O2 equivalents
C**** units: so4 in o2 equivalents
C     SO4 (mg so4/L)* 1 mmol SO4 /98 mg SO4 * 2 mmol O2/ 1 mmol SO4
C     . 32 mg O2 / mmol O2= 0.65306122

      SO40=SO40MG*0.65306122
      K0H1D=0.
      K0H1P=0.
      KMC1=0.0
      K1H1D=XAPPD1**2/S*(O20/KMHSO2) + S
      K1H1P=XAPPP1**2/S*(O20/KMHSO2)
      K2H2D=0.
      K2H2P=0.
      J1=0.
      K3=0.0
      J2=XJC1
      PIE1=PIE1S
      PIE2=PIE2S

C**** Set KL12 using H for SO4
      ITEMP = 10.*TEMPD+1
      DDSO4 = ZL12NOM(ITEMP)*H20
      IF (XJC1 .LE. 0.0000001) THEN
        HSO4 = H20
      ELSE
        HSO4 = SQRT(2.*DDSO4*SO40*H20/XJC1)
      ENDIF

C**** No deeper than H20
      IF(HSO4.GT.H20) HSO4=H20
      KL12SO4=KL12*H20/HSO4

C**** Fractions and overall decay reaction velocity
      FD1=1./(1.+M1(IX,IY)*PIE1)
      FP1=M1(IX,IY)*PIE1/(1.+M1(IX,IY)*PIE1)
      FD2=1./(1.+M2(IX,IY)*PIE2)
      FP2=M2(IX,IY)*PIE2/(1.+M2(IX,IY)*PIE2)
      FP1SO4=FP1
      FP2SO4=FP2
      KHS_1=FP1*XAPPP1**2/S*(O20/KMHSO2)+FD1*XAPPD1**2/S*(O20/KMHSO2)
c     HS and FeS oxidation ratio      
      
      BX(1) = DBLE(S)*DBLE(SO40)
      BX(2) = DBLE(H20)*DBLE(SO4T2TM1)/DBLE(DELTAT)
      BX(3) = DBLE(HS0)*DBLE(S)
      BX(4) = DBLE(H20)*DBLE(HST2TM1)/DBLE(DELTAT)

      AD(1,1) = -DBLE(S)-DBLE(KL12SO4)
      AD(1,2) = DBLE(KL12SO4)
      AD(1,3) = DBLE(KHS_1)
      AD(2,1) = DBLE(KL12SO4)
      AD(2,2) = -(DBLE(DELTAT)*DBLE(KL12SO4)+DBLE(H20))/DBLE(DELTAT)
      AD(3,3) = -DBLE(W2)-DBLE(FP1)*DBLE(W12)-DBLE(FD1)*DBLE(S)
     .          -DBLE(FD1)*DBLE(KL12SO4)-DBLE(KHS_1)
      AD(3,4) = DBLE(FP2)*DBLE(W12)+DBLE(FD2)*DBLE(KL12SO4)
      AD(4,3) = DBLE(W2)+DBLE(FP1)*DBLE(W12)+DBLE(FD1)*DBLE(KL12SO4)
      AD(4,4) = -(DBLE(DELTAT)*DBLE(FP2)*DBLE(W12)
     .          +DBLE(DELTAT)*DBLE(FD2)*DBLE(KL12SO4)
     .          +DBLE(DELTAT)*DBLE(W2)+DBLE(H20)) /DBLE(DELTAT)


      G(1) = ((BX(1)*AD(3,3) - AD(1,3)*BX(3))*AD(4,4) -
     .       BX(1)*AD(3,4)*AD(4,3) + AD(1,3)*AD(3,4)*BX(4) +
     .       AD(1,3)*BX(2)*AD(3,4))/(AD(1,3)*AD(3,4))

      G(2) = ((BX(1)*AD(3,3) - AD(1,3)*BX(3))*AD(4,4) -
     .       BX(1)*AD(3,4)*AD(4,3) + AD(1,3)*AD(3,4)*BX(4))
     .                             /(AD(1,3)*AD(3,4))

      HX(1,1) = (AD(1,1)*AD(3,3)*AD(4,4)-AD(1,1)*AD(3,4)*AD(4,3)+AD(1,3)
     .     *AD(2,1)*AD(3,4))/(AD(1,3)*AD(3,4))
      HX(1,2) = (AD(1,2)*AD(3,3)*AD(4,4)-AD(1,2)*AD(3,4)*AD(4,3)+AD(1,3)
     .     *AD(2,2)*AD(3,4))/(AD(1,3)*AD(3,4))
      HX(2,1) = (AD(1,1)*AD(3,3)*AD(4,4)-AD(1,1)*AD(3,4)*
     .     AD(4,3))/(AD(1,3)*AD(3,4))
      HX(2,2) = (AD(1,2)*AD(3,3)*AD(4,4)-AD(1,2)*AD(3,4)*
     .     AD(4,3))/(AD(1,3)*AD(3,4))

      RA0 = (HX(1,1)*G(2)-G(1)*HX(2,1))*DBLE(KMSO4)
      RA1 = - G(1)*HX(2,1) + HX(1,1)*G(2) +
     .      (HX(1,1)*HX(2,2)-HX(1,2)*HX(2,1))*DBLE(KMSO4) + HX(1,1)*J2
      RA2 = HX(1,1)*HX(2,2)-HX(1,2)*HX(2,1)

      SN1 = 1.                          !solution of a2*q^2+a1*x+a0
      IF (RA1.LE.0.0) SN1 = -1.         !see Num Rec p178
      DISC = -(RA1+SN1*DSQRT(DBLE(RA1)**2-DBLE(RA2)*DBLE(RA0)*4.) )/2.
      IF (DISC .EQ. 0.) THEN
        PRINT*,'Singularity in CH4-SO4 Code'
C        STOP
C LB 20170721
      ELSEIF (DISC .NE. DISC) THEN
        PRINT*,'DISC is NaN'
C        STOP
      ENDIF
      R1 = DISC / RA2
      R2 = RA0 / DISC

      DBLSO42 = R1
      IF (DBLSO42 .LT. 0.) DBLSO42 = R2

      DBLSO41 = -(HX(1,2)*DBLSO42+G(1))/HX(1,1)
      HST1 = -(AD(1,2)*DBLSO42+AD(1,1)*DBLSO41+BX(1))/AD(1,3)
      HST2 = (AD(1,2)*AD(3,3)*DBLSO42+AD(1,1)*AD(3,3)*DBLSO41+BX(1)
     .       *AD(3,3)-AD(1,3)*BX(3))/(AD(1,3)*AD(3,4))
      HS1=FD1*HST1
      HS2=FD2*HST2
c      cdic = (HS1+HS2)*70.
c      ctalk = (HS1+HS2)*106.
      HS2AV=FD2*HST2
      SO42=DBLSO42
      SO42AV=SO42
      SO4T2 = SO42
      SO41=DBLSO41
      XJ2=J2*KMSO4/(SO42+KMSO4)
      SADIC(IX,IY,1) = 0.0
      SATA(IX,IY,1) = 0.0
      XJ2CH4=XJ2
      X1J2=J2*DBLSO42/(SO42+KMSO4)
      JHS=S*(HS1-HS0)
      CSODHS=(XAPPD1**2/S*FD1 + XAPPP1**2/S*FP1)*(O20/KMHSO2)*HST1
      
      SADIC(IX,IY,2) = -1.0*XJ2CH4/64.0
      SATA(IX,IY,2) = 2*XJ2CH4/64.0
      
      SIDIC(IX,IY,1) = 0.0
      SIDIC(IX,IY,2) = 0.0
      SITA(IX,IY,1) = -2*CSODHS/64
      SITA(IX,IY,2) = 0.0

C**** Methane
      CH40 = 0.
      K0H1P=0.
      K1H1P=0.
      K2H2D=0.
      K2H2P=0.
      K1H1D=XAPPCH4**2/S*(O20/(KMCH4O2+O20))+S
      K0H1D=0.
      J1=S*CH40
      K3=0.0
      J2=XJ2
      PIE1=0.0
      PIE2=0.0
      KMC1=0.0

      CALL SEDTSFNL
     .(CH41,CH42,CH42AV,CH4T1,CH4T2,CH4T2AV,CH41TM1,CH4T1TM1,CH4T2TM1,1)

      IF(CH42.GT.CH4SAT) THEN
         CH42=CH4SAT
         CH41 = (CH40*S**2+CH42*KL12*S)/(S**2+KL12*S+
     .           XAPPCH4**2*(O20/(KMCH4O2+O20)))
      ENDIF

C**** Calculate changes in CH4 and HS stored in the sediment
      DCH4T2 = (CH4T2 - CH4T2TM1)*H20/DELTAT
      DHST2  = (HST2 - HST2TM1)*H20/DELTAT

C**** Calculate CSOD
      CSODCH4 = XAPPCH4**2/S*(O20/(KMCH4O2+O20))*CH41
      MEDIC(IX,IY,1) = -1*CSODCH4
      MEDIC(IX,IY,2) = 0.0
      META(IX,IY,1) = 0.0*CSODCH4
      META(IX,IY,2) = 0.0
      
      CSOD    = CSODCH4+CSODHS

C**** Calculate Fluxes
      IF(HST1.NE.0.0) THEN
         FD1=HS1/HST1
      ELSE
         FD1=0.0
      ENDIF
      JCH4      = S*(CH41-CH40)
      JCH4AQ    = S*CH41
      FLUXHS    = S*FD1*HS1
      FLUXHSCH4 = JCH4AQ + FLUXHS

C**** If not flux or SOD or stored then it must escape as gas flux
      JCH4G = 0.
      IF (CH42.EQ.CH4SAT) THEN
         JCH4G = XJ2CH4 - DCH4T2 - CSODCH4 - JCH4AQ
      ENDIF

C**** Volumetric methane and total gas flux (L/m2-d)
      VJCH4G=22.4/64.*JCH4G
c     JGAS=JN2GAS+VJCH4G                   ! jn2gas not computed
C**** **********************************************************

c- SOD function - Evaluate the error
      sod = csod + jo2nh4
      soderr=sod-sod1
      sedf = sod - sod1
      
d       if(iprt.eq.1.and.IX.eq.ixp.and.IY.eq.iyp)then
d         print *,'ixp,iyp',ixp,iyp
d         print *,'sod,csod,jo2nh4,soderr,sod1,sedf'
d         print *, sod,csod,jo2nh4,soderr,sod1,sedf 
d       endif

c-- if this is sod iteration then return
      return
      endif

c-- Entry for Silica and PO4 flux

c-  Silica
      if(iuse.eq.2) then
      k0h1d=0.
      k0h1p=0.
      kmc1=0.0
      k1h1d=s
      k1h1p=0.
      k2h2d=0.
      k2h2p=0.
      j1=s*si0

c-- Oxygen dependency of pie1
      if(o20.lt.o2critsi) then
         pie1=pie2si(ix,iy)*pie1si(ix,iy)**(o20/o2critsi)
      else
         pie1=pie2si(ix,iy)*pie1si(ix,iy)
      endif
      pie2=pie2si(ix,iy)
c- Silica dissolution kinetics
c     d psi/ dt  = -xksi*psi/(psi+kmpsi)*(csisat-si2)
c     d sit2/ dt = +xksi*psi/(psi+kmpsi)*csisat
c                  -xksi*psi/(psi+kmpsi)*si2
      fd2=1./(1.+m2(ix,iy)*pie2)
      k3=xksi*psi/(psitm1+kmpsi)*fd2
      j2=xksi*psi/(psitm1+kmpsi)*csisat + flxsit2
      call sedtsfnl
     .(si1,si2,si2av,sit1,sit2,sit2av,si1tm1,sit1tm1,sit2tm1,1)
      jsi=s*(si1-si0)

c- DIC TALK flux
      cdic1 = nidic(ix,iy,1)+dndic(ix,iy,1)+sadic(ix,iy,1)
     .      + sidic(ix,iy,1)+medic(ix,iy,1)
      cdic2 = nidic(ix,iy,2)+dndic(ix,iy,2)+sadic(ix,iy,2)
     .      + sidic(ix,iy,2)+medic(ix,iy,2)
      cta1 = nita(ix,iy,1)+dnta(ix,iy,1)+sata(ix,iy,1)
     .      + sidic(ix,iy,1)+medic(ix,iy,1)
      cta2 = nita(ix,iy,2)+dnta(ix,iy,2)+sata(ix,iy,2)
     .      + sita(ix,iy,2)+meta(ix,iy,2)
c      if(ix.eq.75)then
c          if(iy.eq.65)then
c              pause
c              endif
c      endif
      k0h1d=0.
      k0h1p=0.
      kmc1=0.0
      k1h1d=0.00010195/DE1(IX,IY)
      k1h1p=0.
      k2h2d=0.
      k2h2p=0.
      j1=0.00010195/DE1(IX,IY)*dic0+cdic1
      if(o20.lt.o2critsi) then
         pie1=pie2si(ix,iy)*pie1si(ix,iy)**(o20/o2critsi)/4.*0.0+250.0
      else
         pie1=pie2si(ix,iy)*pie1si(ix,iy)/4.*0.0+250.0
      endif
      pie2=pie2si(ix,iy)/4.*0.0+80.0
      fd2=1./(1.+m2(ix,iy)*pie2)
      k3=xksi*(poc1+poc2+poc3)/(poc1tm1+poc2tm1+poc3tm1+kmpsi)*fd2
      j2= xksi*(poc1+poc2+poc3)/(poc1tm1+poc2tm1+poc3tm1+kmpsi)*csisat
     .  + flxdict2 + cdic2
      KMC1=0.0
      CALL SEDTSFNL
     .(DIC1,DIC2,DIC2AV,DICT1,DICT2,DICT2AV,
     . DIC1TM1,DICT1TM1,DICT2TM1,1)
c      dic1 = dic0+100
      jdic=0.00010195/DE1(IX,IY)*(dic1-dic0)
C      IF(IX.EQ.107) THEN
C          IF(IY.EQ.174) THEN
C              J1DIC(IX,IY,1) = J1
C          ENDIF
C      ENDIF
C      if(jdic.gt.100.0) then
C          pause
C      endif
      J1DIC(IX,IY,1) = J1
      J1DIC(IX,IY,2) = J2
      JCDIC(IX,IY,1) = 0.0
      JCDIC(IX,IY,2) = JC
      
      k0h1d=0.
      k0h1p=0.
      kmc1=0.0
      k1h1d=0.00010195/DE1(IX,IY)
      k1h1p=0.
      k2h2d=0.
      k2h2p=0.
      j1=0.00010195/DE1(IX,IY)*ta0+cta1
      if(o20.lt.o2critsi) then
         pie1=pie2si(ix,iy)*pie1si(ix,iy)**(o20/o2critsi)/4.*0.0+250.0
      else
         pie1=pie2si(ix,iy)*pie1si(ix,iy)/4.*0.0+250.0
      endif
      pie2=pie2si(ix,iy)/4.*0.0+80.0
      fd2=1./(1.+m2(ix,iy)*pie2)
      k3=xksi*(poc1+poc2+poc3)/(poc1tm1+poc2tm1+poc3tm1+kmpsi)*fd2
      j2= xksi*(poc1+poc2+poc3)/(poc1tm1+poc2tm1+poc3tm1+kmpsi)*csisat
     .  + flxtalkt2 + cta2
      KMC1=0.0
      CALL SEDTSFNL
     .(TA1,TA2,TA2AV,TAT1,TAT2,TAT2AV,
     . TA1TM1,TAT1TM1,TAT2TM1,1)
c      ta1 = ta0+100
      jtalk=0.00010195/DE1(IX,IY)*(ta1-ta0)
      JCTA(IX,IY,1) = 0.0
      JCTA(IX,IY,2) = JC
      J1TA(IX,IY,1) = J1
      J1TA(IX,IY,2) = J2
      
      CHDIC(IX,IY,1) = dic1
      CHDIC(IX,IY,2) = dic2
      CHTA(IX,IY,1) = ta1
      CHTA(IX,IY,2)= ta2
c-- Phosphate

      k0h1d=0.
      k0h1p=0.
      kmc1=0.0
      k1h1d=s
      k1h1p=0.
      k2h2d=0.
      k2h2p=0.
      j1=s*po40
      k3=0.0
      j2=jp + flxpo4t2
c-- Oxygen dependency of pie1
      if(o20.lt.o2crit) then
         pie1=pie1po4n(ix,iy)*pie1po4m(ix,iy)**(o20/o2crit)
      else
         pie1=pie1po4n(ix,iy)*pie1po4m(ix,iy)
      endif
      pie2=pie1po4n(ix,iy)
      call sedtsfnl
     .(po41,po42,po42av,po4t1,po4t2,po4t2av,po41tm1,po4t1tm1,po4t2tm1,1)
      jpo4=s*(po41-po40)

c-- printing for mixing output

      if(ibnth.eq.1)  then
        k1h1d=xappd1**2/s
        k1h1p=xappp1**2/s*(o20/kmhso2)
        psod=w12*(hst2av-hst1)
        dsod=kl12*(hs2av-hs1)
        sodjhs=jhs
        sodpld=sodjhs+dsod
        sodplp=sodpld+psod
      endif
      return

      endif

      return
      end

      SUBROUTINE SEDTSFNL(C1,C2,C2AV,CT1,CT2,CT2AV,C1TM1,
     .                    CT1TM1,CT2TM1,ITYPE)
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON
      real    ans(2),c1,c2,c2av,ct1,ct2,ct2av,c1tm1,ct2tm1

c-  Notation
c     c0 = olw concentration
c     c1 = layer 1 dissolved conc.
c     c2(i) = layer 2-nsedlyr dissolved conc.
c     c2av = layer 2-nsedlyr average conc. for output
c     ct1 = layer 1 total conc.
c     ct2(i) = layer 2-nsedlyr total conc.
c     ct2av = later 2-nsedlyr average conc. for output
c     c1tm1 = c1 at time level t - deltat
c     ct1tm1 = ct1 at time level t - deltat
c     ct2tm1 = ct2 at time level t - deltat
c     h1 =  aerobic layer depth
c     h2 =  anaerobic layer depth
c     h1dot =  d(h1)dt
c     h2dot =  d(h2)dt
c- For c = nh4,no3,hs,po4,si
c-    if itype=1 then return all variables,
c        else return only ct1 and ct2 - for diagenesis computation

      fd1=1./(1.+m1(ix,iy)*pie1)
      fp1=m1(ix,iy)*pie1/(1.+m1(ix,iy)*pie1)
      fd2=1./(1.+m2(ix,iy)*pie2)
      fp2=m2(ix,iy)*pie2/(1.+m2(ix,iy)*pie2)

c- Transport and Decay terms
      f12 = w12*fp1 + kl12*fd1
      f21 = w12*fp2 + kl12*fd2
c     evaluate the MM term at time level t-1
      if(kmc1.ne.0.0) then
        xk0 = (k0h1d*fd1 + k0h1p*fp1)/(kmc1+c1tm1)
      else
        xk0=0.0
      endif
      xk1 = xk0 + k1h1d*fd1 + k1h1p*fp1
      xk2 = k2h2d*fd2 + k2h2p*fp2

c-- aerobic layer displacement flux
      h1dotp=0.5*(h1dot+abs(h1dot))
      h1dotm=-0.5*(h1dot-abs(h1dot))

c-- linear equation coefficients
c     a11 = -f12 -xk1 - w2 + h1dotm - h1dot - h1/deltat
      a11 = -h1dotm -h1dot -h1/deltat -f12 -xk1 -w2
      a12 = f21 + h1dotp
      a21 = f12 + w2 + h1dotm
      b1 = -j1 - h1/deltat*ct1tm1
c     a22= -f21 -xk2 - w2 - k3 -h1dotp -h2dot - h2/deltat
      a22= -h1dotp -h2dot -h2/deltat -f21 -xk2 -w2 -k3
      b2= -j2 - h2/deltat*ct2tm1

c-- Solve the 2x2 set of linear equations
      delta=a11*a22-a12*a21
      if(delta.eq.0.0) then
        print *,'twod is singular: a11,a12,a21,a22'
        print *,a11,a12,a21,a22
        stop
      endif
      ans(1)=(b1*a22-b2*a12)/delta
      ans(2)=(b2*a11-b1*a21)/delta
c        print *,' ans(1),ans(2)...',ans(1),ans(2)

        ct2=ans(2)
        ct2av=ct2

c-- Evaluate other terms if not diagenesis solution
        if(itype.eq.1) then
        ct1=ans(1)
        c1=fd1*ct1
        c2=fd2*ans(2)
        c2av=c2
        endif

        return
        end

      SUBROUTINE SEDTS2TH
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON

      real
     .     kpop1  , kpop2  , kpop3
     .   , kpon1  , kpon2  , kpon3
     .   , kpoc1  , kpoc2  , kpoc3
      equivalence
     .     (kpdiag(1),kpop1) , (kpdiag(2),kpop2) , (kpdiag(3),kpop3)
     .   , (kndiag(1),kpon1) , (kndiag(2),kpon2) , (kndiag(3),kpon3)
     .   , (kcdiag(1),kpoc1) , (kcdiag(2),kpoc2) , (kcdiag(3),kpoc3)
     .  , (dpthta(1),thtapop1),(dpthta(2),thtapop2),(dpthta(3),thtapop3)
     .  , (dnthta(1),thtapon1),(dnthta(2),thtapon2),(dnthta(3),thtapon3)
     .  , (dcthta(1),thtapoc1),(dcthta(2),thtapoc2),(dcthta(3),thtapoc3)

      do 10 itemp=1,350
      temp=float(itemp-1)/10.0 + 0.05
      temp20=temp-20.0
      temp202=temp20/2.0
      zhtanh4s(itemp)=kappnh4s*thtanh4s**temp202
      zhtanh4f(itemp)=kappnh4f*thtanh4f**temp202
      zhtakmnh4(itemp)=kmnh4*thtakmnh4**temp202
      zhtad1(itemp)=kappd1*thtapd1**temp202
      zhtap1(itemp)=kappp1*thtapd1**temp202
      zhtano3s(itemp)=kapp1no3s*thtano3s**temp202
      zhta2no3s(itemp)=k2no3s*thtano3s**temp20
      zhtano3f(itemp)=kapp1no3f*thtano3f**temp202
      zhta2no3f(itemp)=k2no3f*thtano3f**temp20
      zhtaDd0(itemp)=Dd0*thtaDd0**temp20
      zl12nom(itemp)=thtaDd**temp20
      zw12nom(itemp)=thtaDp**temp20
      zhtapon1(itemp)=kpon1*thtapon1**temp20
      zhtapon2(itemp)=kpon2*thtapon2**temp20
      zhtapon3(itemp)=kpon3*thtapon3**temp20
      zhtapoc1(itemp)=kpoc1*thtapoc1**temp20
      zhtapoc2(itemp)=kpoc2*thtapoc2**temp20
      zhtapoc3(itemp)=kpoc3*thtapoc3**temp20
      zhtapop1(itemp)=kpop1*thtapop1**temp20
      zhtapop2(itemp)=kpop2*thtapop2**temp20
      zhtapop3(itemp)=kpop3*thtapop3**temp20
      zhtasi(itemp)=ksi*thtasi**temp20
      zhtajsi(itemp)=thtasi**temp20
      zhtach4(itemp) = kappch4*thtach4**temp202
10    continue
        return
        end

      FUNCTION ZBRENT(FUNC,X1,X2,TOL,IERR)
      external func
c-- modified for two argument FUNC to accomodate sedf
      PARAMETER (ITMAX=100,EPS=3.E-8)
      IERR=0
      A=X1
      B=X2
      FA=FUNC(A,1)
      FB=FUNC(B,1)
C        Check --- Root must bracket ZBRENT
      IF(FB*FA.GT.0.) THEN
         IERR=1
         RETURN
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          FB=FUNC(B,1)
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B,1)
11    CONTINUE
C        ZBRENT exceeded maximum iterations
      IERR=2
      ZBRENT=B
      RETURN
      END

      SUBROUTINE STORETM1
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON

c        dissolved concentrations in layer 1
      nh41tm1s(IX,IY)=nh41
      no31tm1s(IX,IY)=no31
      hs1tm1s(IX,IY)=hs1
      si1tm1s(IX,IY)=si1
      po41tm1s(IX,IY)=po41
c        benthic stress
      bnthstr1s(IX,IY) = bnthstr
c        total concentrations in layer 1
      nh4t1tm1s(IX,IY)=nh4t1
      no3t1tm1s(IX,IY)=no3t1
      hst1tm1s(IX,IY)=hst1
      sit1tm1s(IX,IY)=sit1
      po4t1tm1s(IX,IY)=po4t1
c        total concentrations in layer 2
      nh4t2tm1s(IX,IY)=nh4t2
      no3t2tm1s(IX,IY)=no3t2
      hst2tm1s(IX,IY)=hst2
      sit2tm1s(IX,IY)=sit2
      po4t2tm1s(IX,IY)=po4t2
c        POM concentrations in layer 2
      pon1tm1s(IX,IY)=pon1
      pon2tm1s(IX,IY)=pon2
      pon3tm1s(IX,IY)=pon3
      poc1tm1s(IX,IY)=poc1
      poc2tm1s(IX,IY)=poc2
      poc3tm1s(IX,IY)=poc3
      pop1tm1s(IX,IY)=pop1
      pop2tm1s(IX,IY)=pop2
      pop3tm1s(IX,IY)=pop3
      psitm1s(IX,IY)=psi
c        h1 variables
      sodtm1s(IX,IY)=sod
      o20tm1s(IX,IY)=o20
      Dd0tm1s(IX,IY)=xDd0
      CH4T1TM1S(IX,IY) = CH4T1
      CH4T2TM1S(IX,IY) = CH4T2
      CH41TM1S(IX,IY)  = CH41
      SO4T2TM1S(IX,IY) = SO4T2
C     DIC与TA均为溶解态，因此不区分dic1和dict1之前的错误就是有误dict1没有定义造成的
      dic1tm1s(ix,iy) = dic1
      ta1tm1s(ix,iy) = ta1
      dict1tm1s(ix,iy) = dict1
      dict2tm1s(ix,iy) = dict2
      tat1tm1s(ix,iy) = tat1
      tat2tm1s(ix,iy) = tat2

d     iprt=0
d     ixp=4
d     iyp=6
d     if(idisk.ge.0)  iprt=0
d     if(iprt.eq.1 .and. ix.eq.ixp .and. iy.eq.iyp)  then
d        print *,' storetm1...poc1tm1s=poc1',poc1tm1s(ixp,iyp),poc1
d        print *,' storetm1...pon1tm1s=pon1',pon1tm1s(ixp,iyp),pon1
d        print *,' storetm1...pop1tm1s=pop1',pop1tm1s(ixp,iyp),pop1
d        print *,' storetm1...psitm1s=psi',psitm1s(ixp,iyp),psi
d     endif

C NEW
      BFORMAXS(IX,IY)=BFORMAX
      ISWBNTHS(IX,IY)=ISWBNTH
C END NEW

      return
      end

      SUBROUTINE LOADTM1
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON

      nh41tm1=nh41tm1s(IX,IY)
      no31tm1=no31tm1s(IX,IY)
      hs1tm1=hs1tm1s(IX,IY)
      si1tm1=si1tm1s(IX,IY)
      po41tm1=po41tm1s(IX,IY)
      bnthstr1 = bnthstr1s(IX,IY)
      nh4t1tm1=nh4t1tm1s(IX,IY)
      no3t1tm1=no3t1tm1s(IX,IY)
      hst1tm1=hst1tm1s(IX,IY)
      sit1tm1=sit1tm1s(IX,IY)
      po4t1tm1=po4t1tm1s(IX,IY)
      nh4t2tm1=nh4t2tm1s(IX,IY)
      no3t2tm1=no3t2tm1s(IX,IY)
      hst2tm1=hst2tm1s(IX,IY)
      sit2tm1=sit2tm1s(IX,IY)
      po4t2tm1=po4t2tm1s(IX,IY)
      pon1tm1=pon1tm1s(IX,IY)
      pon2tm1=pon2tm1s(IX,IY)
      pon3tm1=pon3tm1s(IX,IY)
      poc1tm1=poc1tm1s(IX,IY)
      poc2tm1=poc2tm1s(IX,IY)
      poc3tm1=poc3tm1s(IX,IY)
      pop1tm1=pop1tm1s(IX,IY)
      pop2tm1=pop2tm1s(IX,IY)
      pop3tm1=pop3tm1s(IX,IY)
      psitm1=psitm1s(IX,IY)
      sodtm1=sodtm1s(IX,IY)
      o20tm1=o20tm1s(IX,IY)
      Dd0tm1=Dd0tm1s(IX,IY)
      CH4T1TM1 = CH4T1TM1S(IX,IY)
      CH4T2TM1 = CH4T2TM1S(IX,IY)
      CH41TM1  = CH41TM1S(IX,IY)
      SO4T2TM1 = SO4T2TM1S(IX,IY)
      dic1tm1 = dic1tm1s(ix,iy)
      ta1tm1 = ta1tm1s(ix,iy)
      dict1tm1 = dict1tm1s(ix,iy)
      dict2tm1 = dict2tm1s(ix,iy)
      tat1tm1 = tat1tm1s(ix,iy)
      tat2tm1 = tat2tm1s(ix,iy)

d     iprt=0
d     ixp=4
d     iyp=6
d     if(idisk.ge.0)  iprt=0
d     if(iprt.eq.1 .and. ix.eq.ixp .and. iy.eq.iyp)  then
d        print *,' loadtm1...poc1tm1=poc1tm1s',poc1tm1,poc1tm1s(ixp,iyp)
d        print *,' loadtm1...pon1tm1=pon1tm1s',pon1tm1,pon1tm1s(ixp,iyp)
d        print *,' loadtm1...pop1tm1=pop1tm1s',pop1tm1,pop1tm1s(ixp,iyp)
d        print *,' loadtm1...psitm1=psitm1s',psitm1,psitm1s(ixp,iyp)
d     endif

C NEW
      BFORMAX=BFORMAXS(IX,IY)
      ISWBNTH=ISWBNTHS(IX,IY)
C END NEW

      return
      end

      SUBROUTINE STORESED
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON

C     print *,' entry storesed...idisk',idisk
      do 30 IYY=1,NY
      do 30 IXX=1,NX
C CPOC ETC NOW IN MG/M**3
      cpop(IXX,IYY,1)=pop1tm1s(IXX,IYY)
      cpop(IXX,IYY,2)=pop2tm1s(IXX,IYY)
      cpop(IXX,IYY,3)=pop3tm1s(IXX,IYY)
      cpon(IXX,IYY,1)=pon1tm1s(IXX,IYY)
      cpon(IXX,IYY,2)=pon2tm1s(IXX,IYY)
      cpon(IXX,IYY,3)=pon3tm1s(IXX,IYY)
      cpoc(IXX,IYY,1)=poc1tm1s(IXX,IYY)
      cpoc(IXX,IYY,2)=poc2tm1s(IXX,IYY)
      cpoc(IXX,IYY,3)=poc3tm1s(IXX,IYY)
      cpos(IXX,IYY)=psitm1s(IXX,IYY)

d     iprt=0
d     ixp=4
d     iyp=6
d     if(idisk.ge.0)  iprt=0
d     if(iprt.eq.1 .and.ixx.eq.ixp .and. iyy.eq.iyp)  then
d        print *,' storesed...poc(ix,iy)=poc1tm1s',cpoc(ixx,iyy,1)
d    .          ,poc1tm1s(ixx,iyy)
d        print *,' storesed...pon(ix,iy)=pon1tm1s',cpon(ixx,iyy,1)
d    .          ,pon1tm1s(ixx,iyy)
d        print *,' storesed...pop(ix,iy)=pop1tm1s',cpop(ixx,iyy,1)
d    .          ,pop1tm1s(ixx,iyy)
d        print *,' storesed...pos(ix,iy)=psitm1s',cpos(ixx,iyy)
d    .          ,psitm1s(ixx,iyy)
d     endif

   30 continue
      return
      end

      SUBROUTINE LOADSED
      SAVE                                                              
cLIST OFF
      include 'RCACM'
      include 'SEDCM'
      include 'SEDLB'
cLIST ON

d     print *,' entry loadsed...idisk',idisk
      do 30 IYY=1,NY
      do 30 IXX=1,NX
      pop1tm1s(IXX,IYY)=cpop(IXX,IYY,1)
      pop2tm1s(IXX,IYY)=cpop(IXX,IYY,2)
      pop3tm1s(IXX,IYY)=cpop(IXX,IYY,3)
      pon1tm1s(IXX,IYY)=cpon(IXX,IYY,1)
      pon2tm1s(IXX,IYY)=cpon(IXX,IYY,2)
      pon3tm1s(IXX,IYY)=cpon(IXX,IYY,3)
      poc1tm1s(IXX,IYY)=cpoc(IXX,IYY,1)
      poc2tm1s(IXX,IYY)=cpoc(IXX,IYY,2)
      poc3tm1s(IXX,IYY)=cpoc(IXX,IYY,3)
      psitm1s(IXX,IYY)=cpos(IXX,IYY)

d     iprt=0
d     ixp=4
d     iyp=6
d     if(idisk.ge.0)  iprt=0
d     if(iprt.eq.1 .and.ixx.eq.ixp .and. iyy.eq.iyp)  then
d        print *,' loadsed...poc1tm1s=poc(ix,iy)',poc1tm1s(ixx,iyy)
d    .          ,cpoc(ixx,iyy,1)
d        print *,' loadsed...pon1tm1s=pon(ix,iy)',pon1tm1s(ixx,iyy)
d    .          ,cpon(ixx,iyy,1)
d        print *,' loadsed...pop1tm1s=pop(ix,iy)',pop1tm1s(ixx,iyy)
d    .          ,cpop(ixx,iyy,1)
d        print *,' loadsed...psitm1s=pos(ix,iy)',psitm1s(ixx,iyy)
d    .          ,cpos(ixx,iyy)
d     endif

   30 continue
      return
      end

      SUBROUTINE ICREAD(NX,NY,IOUT,ISEDPRNT,ARRAY,SSNAME)
      CHARACTER   SSNAME*25,COMMENT*80
      REAL  ARRAY(NX,NY)
   
      READ(40,1080)  COMMENT
 1080 FORMAT(A80)
      READ(40,1000)  SCALE
 1000 FORMAT(8F10.0)
      DO 50 IY=1,NY
       READ(40,1000)   (ARRAY(IX,IY),IX=1,NX)
   50 CONTINUE

      WRITE(IOUT,1200)  SSNAME,SCALE
 1200 FORMAT(//50X,'SEDIMENT INITIAL CONDITIONS FOR',/53X,A25/
     .   53X,'SCALE FACTOR =',E10.3//3X,'COLUMN      ROW -->'/
     .   10X,7X,'1',11X,'2',11X,'3',11X,'4',11X,'5',11X,'6',11X,'7'
     .   11X,'8',11X,'9',10X,'10'/)
      DO 150 IY=1,NY
       IF(ISEDPRNT.GT.0) WRITE(IOUT,1300)  IY,(ARRAY(IX,IY),IX=1,NX)
 1300  FORMAT(4X,I3,3X,10E12.3,/(10X,10E12.3))
       DO 130 IX=1,NX
         ARRAY(IX,IY)=SCALE*ARRAY(IX,IY)
  130  CONTINUE
  150 CONTINUE
      RETURN
      END

      SUBROUTINE VREAD(NX,NY,IOUT,ISEDPRNT,ARRAY,INUM)

      CHARACTER   COMMENT*80
      REAL  ARRAY(NX,NY)
   
      READ(40,1080)  COMMENT
 1080 FORMAT(A80)
      READ(40,1000)  SCALE
 1000 FORMAT(8F10.0)
      DO 50 IY=1,NY
       READ(40,1000)   (ARRAY(IX,IY),IX=1,NX)
   50 CONTINUE

      IF(INUM.EQ.1)   WRITE(IOUT,1250)
 1250 FORMAT(//52X,'SEDIMENTATION RATES (CM/YR)'/)
      IF(INUM.EQ.2)   WRITE(IOUT,1260)
 1260 FORMAT(//43X,'SEDIMENT SOLID-PHASE MIXING RATES (M**2/DAY)'/)
      IF(INUM.EQ.3)   WRITE(IOUT,1270)
 1270 FORMAT(//41X,'SEDIMENT DISSOLVED-PHASE MIXING RATES (M**2/DAY)'/)
      IF(INUM.EQ.4)   WRITE(IOUT,1271)
 1271 FORMAT(//41X,'SEDIMENT SOLIDS CONCENTRATION - AEROBIC (KG/L)'/)
      IF(INUM.EQ.5)   WRITE(IOUT,1272)
 1272 FORMAT(//41X,'SEDIMENT SOLIDS CONCENTRATION - ANAEROBIC (KG/L)'/)
      IF(INUM.EQ.6)   WRITE(IOUT,1273)
 1273 FORMAT(//41X,'AEROBIC LAYER SI PARTITION COEFFICIENT (L/KG)'/)
      IF(INUM.EQ.7)   WRITE(IOUT,1274)
 1274 FORMAT(//41X,'ANAEROBIC LAYER SI PARTITION COEFFICIENT (L/KG)'/)
      IF(INUM.EQ.8)   WRITE(IOUT,1275)
 1275 FORMAT(//41X,'AEROBIC LAYER PO4 PARTITION COEFFICIENT (L/KG)'/)
      IF(INUM.EQ.9)   WRITE(IOUT,1276)
 1276 FORMAT(//41X,'ANAEROBIC LAYER PO4 PARTITION COEFFICIENT (L/KG)'/)
      WRITE(IOUT,1300)  SCALE
 1300 FORMAT(53X,'SCALE FACTOR =',E10.3//3X,'COLUMN      ROW -->'/
     .   10X,7X,'1',11X,'2',11X,'3',11X,'4',11X,'5',11X,'6',11X,'7'
     .   11X,'8',11X,'9',10X,'10'/)
      DO 150 IY=1,NY
       IF(ISEDPRNT.GT.0) WRITE(IOUT,1350)  IY,(ARRAY(IX,IY),IX=1,NX)
 1350  FORMAT(4X,I3,3X,10E12.3,/(10X,10E12.3))
       DO 130 IX=1,NX
         ARRAY(IX,IY)=SCALE*ARRAY(IX,IY)
  130  CONTINUE
  150 CONTINUE
      RETURN
      END
