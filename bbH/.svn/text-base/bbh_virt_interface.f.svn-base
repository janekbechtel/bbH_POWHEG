      subroutine bbh_virtual_interface(nlo,p,flav,born,virt)
      implicit none

      include 'nlegborn.h'
      include 'coupl.inc'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'

      real*8 alphas2
      real*8 QCDL4,QCDL5,TMAS
      common/W50512/QCDL4,QCDL5
      integer NFL,LO
      common/W50511/NFL,LO,TMAS

      real*8 mtop,alphas,mb
      common/topas/mtop,alphas,mb
      real*8 mh
      common/smhiggs/mh
      real*8 gmu
      common/gf/gmu
      real*8 delta,muedr
      common/dimreg/delta,muedr
      integer nlf
      common/number_lflav/nlf

      integer nlo
      real * 8 p(0:3,nlegborn)
      integer flav(nlegborn)
      real * 8 born,virt
      integer flav_id
      real*8 sig0qq(2),sig0gg(2),dum
      real*8 sigvqq(2),sigvgg(2)
c coefficients of IR single and double poles (for testing only)
      real*8 sigv_irpole(2,2),sigvqq_ir(2,2),sigvgg_ir(2,2)
*
c
c four-momenta and invariants
c
      real*8 pt(4),ptb(4),ph(4),b1(4),b2(4)
      real*8 sp,t1(2),t2(2),u1(2),u2(2)

      logical ps_test
      parameter (ps_test=.false.) 

c
c pwhg version:      
      mtop=tmass          
      mb=bmass
      mh=hmass

      alphas=g**2/4d0/pi !should be identical to st_alpha
c
      gmu = ph_gfermi

ccccccccccccccc
* 
      delta=1d0
c renormalization scale
      muedr = dsqrt(st_muren2) !value set in POWHEG
c
c initialize:
      born = 0d0
      virt = 0d0
c
c determine which channel we need:
c
      if (flav(1).eq.0) then !gg
         flav_id = 2
      elseif (flav(1).gt.0) then 
         flav_id = 1 !gives qqbar
      else 
         flav_id = -1 !gives qbar q
      endif  

c four momenta:
c bottom
      pt(1:3) = p(1:3,4)
      pt(4)   = p(0,4)
c
c anti-bottom
      ptb(1:3) = p(1:3,5)
      ptb(4)   = p(0,5)
c
c Higgs 
      ph(1:3) = p(1:3,3)
      ph(4)   = p(0,3)
c
c beam
      if (flav_id.ge.0) then !gg or qqbar
c quark:
      b1(1:3) = p(1:3,1)
      b1(4)   = p(0,1)
c
c anti-quark
      b2(1:3) = p(1:3,2)
      b2(4)   = p(0,2)

      else !qbar q
c quark:
      b2(1:3) = p(1:3,1)
      b2(4)   = p(0,1)
c
c anti-quark
      b1(1:3) = p(1:3,2)
      b1(4)   = p(0,2)
c
      endif
c

ccccccc
c
c test only:
      if (ps_test) then 

         flav_id = 3
c
         mtop=174d0
         mb=4.6d0
         mh=120d0

c for alphas-test only:
c         mh = 125d0
*     
         delta=1d0
c     renormalization scale
         muedr=mb+mh/2d0
c     alphas is calculated at 2-loop order. For 1-loop running use LO=1 and
c     the corresponding LO choices for QCDL4 and QCDL5
         NFL=4
         nlf=NFL
         LO=2
         tmas=mtop
c     QCDL4 and QCDL5 have to be chosen the same as used in the PDFs (here: CTEQ6)
c     NLO
         if(LO.eq.2)then
            QCDL4=0.326d0
            QCDL5=0.226d0
         else
c     LO
            QCDL4=0.215d0
            QCDL5=0.165d0
         endif
         alphas=alphas2(muedr)
         write(6,*)'bbh_interface: alphas=',alphas
c
c four-momenta:
c bottom
      pt(1) =396.0826028513954d0   
      pt(2) =72.54310572811833d0     
      pt(3) =-229.9239044928934d0     
      pt(4) =463.7133730306625d0     
c anti-bottom
      ptb(1) = -53.02844924829163d0     
      ptb(2) = -106.3429914857804d0     
      ptb(3) = 67.86970358133982d0     
      ptb(4) = 136.9244497233506d0     
c Higgs
      ph(1) =-343.0541536031038d0     
      ph(2) = 33.79988575766206d0       
      ph(3) =  162.0542009115537d0     
      ph(4) =  399.3621772459870d0    
c beam momenta
      b1(1) =0.00000000000000E+00
      b1(2) =0.00000000000000E+00
      b1(3) =0.50000000000000E+03
      b1(4) =0.50000000000000E+03
      b2(1) =0.00000000000000E+00
      b2(2) =0.00000000000000E+00
      b2(3) =-0.50000000000000E+03
      b2(4) =0.50000000000000E+03
c
c     
      endif !ps_test
c
c calculate the invariants :
c
      sp=(b1(4)+b2(4))**2-((b1(1)+b2(1))**2+
     $     (b1(2)+b2(2))**2+(b1(3)+b2(3))**2)
      t1(1)=(b1(4)-pt(4))**2-((b1(1)-pt(1))**2+
     $     (b1(2)-pt(2))**2+(b1(3)-pt(3))**2)
      t2(1)=(b2(4)-ptb(4))**2-((b2(1)-ptb(1))**2+
     $     (b2(2)-ptb(2))**2+(b2(3)-ptb(3))**2)
      u1(1)=(b2(4)-pt(4))**2-((b2(1)-pt(1))**2+
     $     (b2(2)-pt(2))**2+(b2(3)-pt(3))**2)
      u2(1)=(b1(4)-ptb(4))**2-((b1(1)-ptb(1))**2+
     $     (b1(2)-ptb(2))**2+(b1(3)-ptb(3))**2)
*
      t1(2)=(b2(4)-pt(4))**2-((b2(1)-pt(1))**2+
     $     (b2(2)-pt(2))**2+(b2(3)-pt(3))**2)
      t2(2)=(b1(4)-ptb(4))**2-((b1(1)-ptb(1))**2+
     $     (b1(2)-ptb(2))**2+(b1(3)-ptb(3))**2)
      u1(2)=(b1(4)-pt(4))**2-((b1(1)-pt(1))**2+
     $     (b1(2)-pt(2))**2+(b1(3)-pt(3))**2)
      u2(2)=(b2(4)-ptb(4))**2-((b2(1)-ptb(1))**2+
     $     (b2(2)-ptb(2))**2+(b2(3)-ptb(3))**2)

c     
c calculate the matrix elements squared at the parton level:	
c qqg=1 -> qqbar , qqg=2 -> gg, qqg=3 -> qqbar and gg
c
      if (ps_test) then 
         call bbhvirt(nlo,3,sp,t1,t2,u1,u2,sig0qq,sig0gg,sigvqq,sigvgg,
     $        sigvqq_ir,sigvgg_ir)
c         if (nlo.eq.0) then 
            print*,'nlo=',nlo
            write(6,*)'LO, qq:',sig0qq
            write(6,*)'LO, gg:',sig0gg
c         elseif (nlo.eq.1) then 
            write(6,*)'virt, qq:',sigvqq
            write(6,*)'virt, gg:',sigvgg
            write(6,*)'single pole, qq:',sigvqq_ir(1,1),sigvqq_ir(1,2)
            write(6,*)'double pole, qq:',sigvqq_ir(2,1),sigvqq_ir(2,2)
            write(6,*)'single pole, gg:',sigvgg_ir(1,1),sigvgg_ir(1,2)
            write(6,*)'double pole, gg:',sigvgg_ir(2,1),sigvgg_ir(2,2)
c         endif
         
         stop
      endif !pstest

      call bbhvirt(nlo,abs(flav_id),
     $     sp,t1,t2,u1,u2,sig0qq,sig0gg,sigvqq,sigvgg,
     $     sigvqq_ir,sigvgg_ir)

      if (abs(flav_id).eq.1) then !qqbar
         born = sig0qq(1)
         virt = sigvqq(1)
      elseif (abs(flav_id).eq.2) then    !gg
         born = sig0gg(1)
         virt = sigvgg(1)
      else 
         print*,'illegal flav_id in bbh_virt:',flav_id
         stop
      endif   

c      if (ps_test) then
c         write(6,*)'LO:',sig0qq,sig0gg
c         write(6,*)'virt:',sigvqq,sigvgg
c         print*
c         stop
c      endif

c output for these momenta:a with nlf=5 in modified MSbar, ie the
c top quark is decoupled as heavy quark (OS scheme for bottom):
c alphas=  0.1244681550355927     
c top-Yukawa: -0.7066895038715881     
c bottom-Yukawa: -1.8682596079363822E-002
c LO, qq:  1.9814310862278459E-009  1.9814310862278459E-009
c LO, gg:  3.2568570893942059E-009  3.2568570893941953E-009
c without finite parts from the expansion of the pole parts:
c virt, qq: -9.1627999026167970E-009 -7.8095368281486550E-009
c virt, gg: -2.0887184184681491E-008 -2.0887184184681473E-008
c with the finite parts from the expansion from the pole parts:
c virt, qq:  4.9438362700313760E-010  8.1792865078980022E-010
c virt, gg:  2.2540286485345509E-009  2.2540286485344780E-009
c single pole, qq:  9.6729380746375746E-010  8.2503071037576845E-010
c double pole, qq: -1.0467093130873819E-010 -1.0467093130873819E-010
c single pole, gg:  1.8401591126621714E-009  1.8401591126621685E-009
c double pole, gg: -3.8710460377631705E-010 -3.8710460377631485E-010
*
c output for these momenta with nlf=4 in modified MSbar, ie both the
c top and bottom quarks are decoupled as heavy quarks (OS scheme for bottom):
c alphas=  0.1196560344822110     
c top-Yukawa: -0.7066895038715881     
c bottom-Yukawa: -1.8682596079363822E-002
c LO, qq:  1.8311827120681732E-009  1.8311827120681732E-009
c LO, gg:  3.0098954433631775E-009  3.0098954433631659E-009
c virt, qq:  5.6208069358288171E-010  8.4953166828598601E-010
c virt, gg:  1.9279199072803790E-009  1.9279199072803175E-009
c single pole, qq:  8.5938440506836944E-010  7.3299190042212636E-010
c double pole, qq: -9.2994047244619047E-011 -9.2994047244619047E-011
c single pole, gg:  1.5966612363422225E-009  1.5966612363422206E-009
c double pole, gg: -3.4391997245159835E-010 -3.4391997245159622E-010
c
      return 
      end
      
