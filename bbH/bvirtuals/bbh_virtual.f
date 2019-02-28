
      subroutine bbhvirt(nlo,qqg,sp,t1,t2,u1,u2,
     $     sig0qq,sig0gg,sigvqq,sigvgg,sigvqq_ir,sigvgg_ir)
***********************************************************
c virtual O(alpha_s) corrections to bbbarH production
c extracted from bbh code in June 2014, D.W.
***********************************************************
      implicit none
      integer nlo
      integer i,j
      real*8 sig0(2),sig0qq(2),sig0gg(2)
      real*8 sigv(2),sigvqq(2),sigvgg(2)
      real*8 sp,t1(2),t2(2),u1(2),u2(2)
      real*8 mtop,alphas,mb
      common/topas/mtop,alphas,mb
      real*8 mh
      common/smhiggs/mh
      real*8 pi,pi2,w2
      common/para/pi,pi2,w2
      real*8 delta,muedr
      common/dimreg/delta,muedr
      real*8 ceps1(2,2),ceps2(2,2)
      common/irpoles/ceps1,ceps2
      real*8 irfinite(2,2)
      real*8 gmu,yukt,yukb,yukf
      common/gf/gmu
      common/yukawa/yukt,yukb,yukf
      integer kmin,kmax
      common/yuksum/kmin,kmax
c running Yukawa coupling
      real*8 b0,b1,g1,a1
      real*8 mbr
      integer imsbar
      common/mb_running/mbr,imsbar
      integer nlf
      common/number_lflav/nlf
      real*8 alphasq,alphas2
c
c o_fac: overall factor (coupling constant)
c
      real*8 o_fac
c matrix elements squared
      real*8 matb_gg,matb_qq
c coefficients of IR single and double poles
      real*8 sigv_irpole(2,2),sigvqq_ir(2,2),sigvgg_ir(2,2)
c switch for qqbar and/or gluon fusion
      integer qqg
c
c initialization
      do i=1,2
         sig0qq(i)=0d0
         sig0gg(i)=0d0
         sigvqq(i)=0d0
         sigvgg(i)=0d0
         do j=1,2
            sigvqq_ir(i,j)=0d0
            sigvgg_ir(i,j)=0d0
         enddo
      enddo
c constants
      pi=4d0*datan(1d0)
      pi2=pi*pi
      w2=dsqrt(2d0)

c set in init_couplings:
c consider top and bottom Yukawa in loops? 
c (k=1 for top, k=2 for bottom)
c      kmax = 2 ! do not change
c      kmin = 1 ! can be set to 2 -> no top loop contribution
c
c overall factor and Yukawa couplings
c
c MSM-Higgs-fermion coupling:
c          eta
*          -----
c     yukt=-mtop/2d0/sw/mw*dsqrt(4d0*pi*alpha)
c      gmu=1.16639d-5
      yukt=-mtop*dsqrt(gmu*w2)
      yukb=-mb*dsqrt(gmu*w2)
c MSSM
c bbh0
c      yukt=-mtop*dsqrt(gmu*w2)*dcos(alpha)/dsin(beta)
c      yukb=mb*dsqrt(gmu*w2)*dsin(alpha)/dcos(beta)
c bbH0
c       yukt=-mtop*dsqrt(gmu*w2)*dsin(alpha)/dsin(beta)
c       yukb=-mb*dsqrt(gmu*w2)*dcos(alpha)/dcos(beta)


      yukb=yukb*mbr/mb
c      write(6,*)'top-Yukawa:',yukt
c      write(6,*)'bottom-Yukawa:',yukb
c      write(6,*)'alpha_s=',alphas
      yukf=yukb

*
      o_fac=yukf**2*(4d0*pi*alphas)**2

c     
c Born-matrix element squared:
c     
c qqbar annihilation
      do i=1,2
         sig0(i)=0d0
      enddo
      if (qqg.eq.1.or.qqg.eq.3) then
         sig0(1)=matb_qq(sp,t1(1),t2(1),u1(1),u2(1),mb,mh)
         sig0(2)=matb_qq(sp,t1(2),t2(2),u1(2),u2(2),mb,mh)
         do i=1,2
            sig0qq(i)=o_fac*sig0(i)
         enddo
      end if
c gluon fusion
      do i=1,2
         sig0(i)=0d0
      enddo
      if (qqg.eq.2.or.qqg.eq.3) then
         sig0(1)=matb_gg(sp,t1(1),t2(1),u1(1),u2(1),mb,mh)
         sig0(2)=matb_gg(sp,t1(2),t2(2),u1(2),u2(2),mb,mh)
         do i=1,2
            sig0gg(i)=o_fac*sig0(i)
         enddo
      endif

      if (nlo.eq.0) return
c
c virtual O(alpha_s) corrections 
c
c renormalization constants      
c      do j=1,3
c         muedr=10d0**(j**2)
      call renorm(mb,mtop,mh)
c qqbar annihilation:
      do i=1,2
         sigv(i)=0d0
         do j=1,2
            sigv_irpole(i,j)=0d0
         enddo
      enddo
      if (qqg.eq.1.or.qqg.eq.3) then
         call mvirt_qq(sp,t1(1),t2(1),u1(1),u2(1),mb,mh,sigv(1),
     $        sigv_irpole(1,1),sigv_irpole(2,1))
         call mvirt_qq(sp,t1(2),t2(2),u1(2),u2(2),mb,mh,sigv(2),
     $        sigv_irpole(1,2),sigv_irpole(2,2))
c add extra term when using the MSbar scheme for the b Yukawa coupling:
         if(imsbar.eq.1) then
            sigv(1)=sigv(1)+sig0qq(1)/o_fac*
     $           8d0*4d0/3d0*(1d0-3d0/2d0*dlog(mb/muedr))
            sigv(2)=sigv(2)+sig0qq(2)/o_fac*
     $           8d0*4d0/3d0*(1d0-3d0/2d0*dlog(mb/muedr))
         endif
         do i=1,2
            sigvqq(i)=o_fac*sigv(i)*alphas/4d0/pi
c IR poles
c single
            sigvqq_ir(1,i)=o_fac*sigv_irpole(1,i)*alphas/4d0/pi
c double
            sigvqq_ir(2,i)=o_fac*sigv_irpole(2,i)*alphas/4d0/pi
         enddo
      endif
c gluon fusion:
      do i=1,2
         sigv(i)=0d0
         do j=1,2
            sigv_irpole(i,j)=0d0
         enddo
      enddo
      if (qqg.eq.2.or.qqg.eq.3) then
         call mvirt_gg(sp,t1(1),t2(1),u1(1),u2(1),mb,mh,sigv(1),
     $        sigv_irpole(1,1),sigv_irpole(2,1))
         call mvirt_gg(sp,t1(2),t2(2),u1(2),u2(2),mb,mh,sigv(2),
     $        sigv_irpole(1,2),sigv_irpole(2,2))
         if(imsbar.eq.1) then
            sigv(1)=sigv(1)+sig0gg(1)/o_fac*
     $           8d0*4d0/3d0*(1d0-3d0/2d0*dlog(mb/muedr))
            sigv(2)=sigv(2)+sig0gg(2)/o_fac*
     $           8d0*4d0/3d0*(1d0-3d0/2d0*dlog(mb/muedr))
         endif
         do i=1,2
            sigvgg(i)=o_fac*sigv(i)*alphas/4d0/pi
c           sigvgg(i)=sigvgg(i)-2d0*sig0gg(i)*alphas/4d0/pi*
c     $        (-2d0/3d0*nlf+11d0/3d0*3d0)*dlog(muedr**2/sp)
c IR poles
c single
            sigvgg_ir(1,i)=o_fac*sigv_irpole(1,i)*alphas/4d0/pi
c double
            sigvgg_ir(2,i)=o_fac*sigv_irpole(2,i)*alphas/4d0/pi
         enddo
      endif
c      write(6,*)muedr,sigvqq(1),sigvgg(1)
c      enddo
      return
      end
