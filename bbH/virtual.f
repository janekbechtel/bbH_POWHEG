      subroutine setvirtual(p,vflav,virtual)
c Virtual needs to be provided by the user and put here
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real*8 virt_amp2,born,dum
      real * 8 virtual

      real *8 powheginput
      external powheginput 

      logical, save :: firsttime = .true. 

      integer fakevirt
      save fakevirt 

      virtual = 0d0

      if (firsttime) then
         fakevirt=powheginput("#fakevirt")
         if (fakevirt == 1) write(*,*) 'WARNING: Using fakevirt !'
         firsttime = .false.
      endif

      if(fakevirt.eq.1) then  

c calculate the RW matrix elements squared at the parton level:	
         call bbh_virtual_interface(0,p,vflav,born,dum)
         virtual = 0.2d0*born
 
      else !true virtuals   

c calculate the RW matrix elements squared at the parton level:	
         call bbh_virtual_interface(1,p,vflav,born,virt_amp2)
         virtual = virt_amp2

c convert to POWHEG normalization:
         virtual = virtual*(2d0*pi/st_alpha)

      endif   

      end
