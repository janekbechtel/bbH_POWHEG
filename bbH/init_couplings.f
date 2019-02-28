      subroutine init_couplings
      implicit none
      include "coupl.inc"
      include 'PhysPars.h'
      include "pwhg_physpar.h"

      real * 8 powheginput
      external powheginput
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      data called/.false./
      save called
      if(called) then
         return
      else
         called=.true.
      endif


*********************************************************
***********         MADGRAPH                 ************
*********************************************************
c Parameters are read from the MadGraph param_card.dat,
c except the strong coupling constant, which is defined
c somewhere else
      call setpara
      call madtophys
 
c     Set here lepton and quark masses for momentum reshuffle in the LHE event file
      physpar_ml(1) = 0.51099891d-3
      physpar_ml(2) = 0.1056583668d0
      physpar_ml(3) = 1.77684d0
      physpar_mq(1) = 0.33d0     ! down
      physpar_mq(2) = 0.33d0     ! up
      physpar_mq(3) = 0.50d0     ! strange
      physpar_mq(4) = 1.50d0     ! charm
      physpar_mq(5) = 4.5d0      ! bottom

      end


      subroutine lh_readin
c overrides the lh_readin subroutine in MODEL/couplings.f;
      implicit none
      include 'coupl.inc'
      include 'PhysPars.h'
      double precision  Two, Four, Rt2, Pi
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Rt2   = 1.414213562d0 )
      parameter( Pi = 3.14159265358979323846d0 )
 
c     Common to lh_readin and printout
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb

      real*8 sthw,cthw,g2

      real*8 mbr
      integer imsbar
      common/mb_running/mbr,imsbar

      real * 8 powheginput
      external powheginput

      integer kmin,kmax
      common/yuksum/kmin,kmax

c consider top and bottom Yukawa in loops? 
c (k=1 for top, k=2 for bottom)  
      kmax = 2 ! do not change
      if(powheginput('#toploop').eq.0) then
         kmin = 2               ! no top loop contribution
         print*,'top loop effects disregarded'
      else
         kmin = 1               ! include top loop contribution
         print*,'top loop effects included'
      endif   

c decide about renormalization scheme for b-quarks:
c      imsbar = 1 ! (1=MSbar;0=OnShell) 
c      imsbar = 0 ! (1=MSbar;0=OnShell)       
      if(powheginput('#msbar').eq.1) then
         print*,'MSbar renormalization used for b-quark'
         imsbar = 1
      elseif(powheginput('#msbar').eq.0) then
         print*,'OnShell renormalization used for b-quark'
         imsbar = 0 
      else
         stop 'check setting for renormalization scheme'
      endif   
      
      alfas = 0.119d0 !initial dummy value

      lmass = 0d0
      mtaMS = 1.777d0
      cmass = 0d0
      mcMS = 0d0

      tmass = 173d0
      mtMS = 173d0
      twidth = 0d0 ! this value is not used in MEs, but needs to be initialized

      bmass = 4.75d0
      mbMS = 4.75d0 ! initial value only (overruled for imsbar=1)

      hmass = powheginput('hmass')
      hwidth = powheginput('hwidth')

cccccccccccccccccccccccccccccccccccccccccccccccc
c check if chosen value of hmass can be handled:
      if ((hmass.lt.100d0.or.hmass.gt.200).or.
     &    (hwidth.gt.2d0)) then
         print*,' ################################################'
         print*,' this code is applicable only for bbH production '
         print*,' with a Higgs mass in the range 100 < mH < 200 GeV, '
         print*,' and decay width Gamma_H < 2 GeV, '
         print*,' where Higgs decays can be treated '
         print*,' using a narrow-width approximation '
         print*,' and off-shell effects are small'
         print*,''
         print*,' you entered a Higgs mass (in GeV) of ',hmass
         print*,' and Higgs width (in GeV) of ', hwidth
         print*,' please change your settings'
         print*,' ################################################'
         stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
          
      zwidth=2.441d0
      wwidth=2.0476d0

      zmass = 91.1876d0 
      wmass = 80.385d0
c
ccc
c
      ph_Wmass = wmass
      ph_Zmass = zmass
      
      gfermi = 0.1166390d-4
      ph_sthw2 = 1d0 - (ph_Wmass/ph_Zmass)**2
      sthw = SQRT(ph_sthw2)
      cthw = SQRT(1.d0 -ph_sthw2 )
      G2 = SQRT(8.d0*GFERMI/SQRT(2.d0))*ph_Zmass*cthw
      alpha = g2**2*ph_sthw2/(4.d0*PI)

c
c CKM matrix set to unit in MEs:

      Vud=1d0
      Vus=1d-10
      Vub=1d-10
      Vcd=1d-10
      Vcs=1d0
      Vcb=1d-10
      Vtd=1d-10
      Vts=1d-10
      Vtb=1d0

      end
 
      subroutine set_ebe_couplings
      implicit none
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include "coupl.inc"

      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud

c running Yukawa coupling      
      double precision  Zero
      parameter( Zero = 0.0d0 )

      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 pwhg_alphas0
      external pwhg_alphas0
      real*8 b0,b1,g1,a1
      real * 8 mb,mbr,mbsq,alphasq,alphas
      integer imsbar,nlf
      common/number_lflav/nlf
      common/mb_running/mbr,imsbar

      real *8 powheginput
      external powheginput
     
      logical init
      data init/.true./
      save init

      logical yuk_comp
      data yuk_comp/.false./
      save yuk_comp

      integer nloop
      data nloop/2/
      save nloop
cc
c
c QCD coupling constant:
c
c st_alphas is computed by powheg and can be used here
c
      G=sqrt(st_alpha*4d0*pi)
      GG(1)=-G
      GG(2)=-G
      
      if (imsbar.eq.0) then 
         mbr = bmass
         return                 !OnShell renormalization
      endif   

c for MSbar scheme compute runnning bmass (used in Yukawa coupling only):  
      if (init) then 
         if(powheginput('#runningscales').eq.1) then
            yuk_comp=.true.
         endif   
         if(powheginput('#bornonly').eq.1) then
            print*,'mb(MSbar): one-loop running'
            nloop = 1
         else   
            print*,'mb(MSbar): two-loop running'
         endif   
      endif   

      if (init.or.yuk_comp) then

         init=.false.

cccccccccc
c   
c for MSbar scheme, compute running b-mass and Yukawa:
         nlf=4 
         mbsq = bmass**2
         mb = bmass

         alphasq = pwhg_alphas(mbsq,st_lambda5MSB,st_nlight)
         alphas = st_alpha

         b0=(11d0-2d0/3d0*nlf)/4d0
         b1=(102d0-38d0/3d0*nlf)/16d0
         g1=(202d0/3d0-20d0/9d0*nlf)/16d0
         a1=-b1/b0/b0+g1/b0

c running of mb(mu): solution of RGE  
         if (nloop.eq.1) then 
c one-loop running for LO calculation
            mbr=mb*(alphas/alphasq)**(1d0/b0)
         else
c two-loop running for NL0 calculation (expanded)
            mbr=mb*(alphas/alphasq)**(1d0/b0)*
     $           ((1d0+a1*(alphas-alphasq)/pi))*
     $           (1d0-4d0/3d0*alphasq/pi)
         endif !loop-running

        mbMS = mbr

         if(mbMS.gt.1d0) then
            ghbot(1) = dcmplx( -mbMS/v, Zero )
         else
            ghbot(1) = dcmplx( Zero, Zero )
         endif
         ghbot(2) = ghbot(1)
      
      endif

      return
      end


      subroutine madtophys
      implicit none
      include 'coupl.inc'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      real * 8 e_em,g_weak
      real*8 masswindow
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb

      e_em=gal(1)
      ph_unit_e=e_em
      ph_alphaem=e_em**2/(4*pi)
      ph_sthw2=1-(wmass/zmass)**2
      ph_sthw=sqrt(ph_sthw2)
      g_weak=e_em/ph_sthw
      ph_gfermi=sqrt(2d0)*g_weak**2/(8*wmass**2)

      ph_Hmass = hmass
      ph_Zwidth = zwidth
      ph_Wwidth = wwidth
      ph_Hwidth = hwidth
 
      ph_tmass = tmass
      ph_twidth = twidth !1.5083d0 ! this value is used for top decays

      ph_bmass = bmass

      ph_HmHw = ph_Hmass * ph_Hwidth
      ph_WmWw = ph_Wmass * ph_Wwidth
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_Wmass2 = ph_Wmass**2
      ph_Zmass2 = ph_Zmass**2
      ph_Hmass2 = ph_Hmass**2

c     set mass windows around the resonance peak 
c     It is used in the generation of the Born phase space
      masswindow = 1000d0
      ph_Hmass2low=max(0d0,ph_Hmass-masswindow*ph_Hwidth)
      ph_Hmass2low=ph_Hmass2low**2
      ph_Hmass2high=(ph_Hmass+masswindow*ph_Hwidth)**2

c     CKM from PDG 2010 (eq. 11.27)
      ph_CKM(1,1)=Vud
      ph_CKM(1,2)=Vus
      ph_CKM(1,3)=Vub
      ph_CKM(2,1)=Vcd
      ph_CKM(2,2)=Vcs
      ph_CKM(2,3)=Vcb
      ph_CKM(3,1)=Vtd
      ph_CKM(3,2)=Vts
      ph_CKM(3,3)=Vtb

      end


