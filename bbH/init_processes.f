      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include "coupl.inc"
      include 'LesHouches.h'
      integer i,nf_loc

      common /nf_proc/nf_loc

      integer hdecaymode
      real * 8 powheginput
      external powheginput

c******************************************************
c select Higgs decay mode:
c******************************************************
c     decay products of the Higgs boson:
      hdecaymode=powheginput('#hdecaymode')
      if (hdecaymode.lt.0) then
c     default: no Higgs boson decay
         hdecaymode=-1
      endif
c     change the LHUPI id of the process according to vector boson id
c     and decay
      lprup(1)=10000+hdecaymode ! 10000+idup of Higgs decay product 

c******************************************************
c     index of the first light coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=6
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Set here the number of light flavours
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      st_nlight=5 ! this number is used for computation of pwhg_alphas
      print*,'number of quark flavors for alpha_s =',st_nlight

      nf_loc = 4 ! number of incoming light flavors 
      if (nf_loc.gt.4) then
         print*,'Nf > 4 is not supported by this code'
         stop
      endif

 
      call init_processes_born
      call init_processes_real

      call init_couplings
      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.nf_loc) then
            flst_lightpart=i
            exit
         endif
      enddo
 
      end
 
 
 
      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      integer nf_loc
      common /nf_proc/nf_loc
 
      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=           1
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=           5
      flst_born(   5,   1)=          -5
 
      flst_born(   1,   2)=           1
      flst_born(   2,   2)=          -1
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=           5
      flst_born(   5,   2)=          -5
 
      flst_born(   1,   3)=          -2
      flst_born(   2,   3)=           2
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=           5
      flst_born(   5,   3)=          -5
 
      flst_born(   1,   4)=           2
      flst_born(   2,   4)=          -2
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=           5
      flst_born(   5,   4)=          -5
 
      flst_born(   1,   5)=          -4
      flst_born(   2,   5)=           4
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=           5
      flst_born(   5,   5)=          -5
 
      flst_born(   1,   6)=           4
      flst_born(   2,   6)=          -4
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=           5
      flst_born(   5,   6)=          -5
 
      flst_born(   1,   7)=          -3
      flst_born(   2,   7)=           3
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=           5
      flst_born(   5,   7)=          -5
 
      flst_born(   1,   8)=           3
      flst_born(   2,   8)=          -3
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=           5
      flst_born(   5,   8)=          -5
 
      flst_born(   1,   9)=           0
      flst_born(   2,   9)=           0
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=           5
      flst_born(   5,   9)=          -5
 
      flst_nborn=           9

      return
      if (nf_loc.eq.5) then
 
      flst_born(   1,   10)=          -5
      flst_born(   2,   10)=           5
      flst_born(   3,   10)=          25
      flst_born(   4,   10)=           5
      flst_born(   5,   10)=          -5
 
      flst_born(   1,   11)=           5
      flst_born(   2,   11)=          -5
      flst_born(   3,   11)=          25
      flst_born(   4,   11)=           5
      flst_born(   5,   11)=          -5


      flst_nborn=           11

      endif
 
      end
 
 
 
      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      integer nf_loc
      common /nf_proc/nf_loc
 
      flst_real(   1,   1)=          -1
      flst_real(   2,   1)=           1
      flst_real(   3,   1)=          25
      flst_real(   4,   1)=           5
      flst_real(   5,   1)=          -5
      flst_real(   6,   1)=           0
 
      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=           0
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=           5
      flst_real(   5,   2)=          -5
      flst_real(   6,   2)=          -1
 
      flst_real(   1,   3)=           1
      flst_real(   2,   3)=          -1
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=           5
      flst_real(   5,   3)=          -5
      flst_real(   6,   3)=           0
 
      flst_real(   1,   4)=           1
      flst_real(   2,   4)=           0
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=           5
      flst_real(   5,   4)=          -5
      flst_real(   6,   4)=           1
 
      flst_real(   1,   5)=          -2
      flst_real(   2,   5)=           2
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=           5
      flst_real(   5,   5)=          -5
      flst_real(   6,   5)=           0
 
      flst_real(   1,   6)=          -2
      flst_real(   2,   6)=           0
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=           5
      flst_real(   5,   6)=          -5
      flst_real(   6,   6)=          -2
 
      flst_real(   1,   7)=           2
      flst_real(   2,   7)=          -2
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=           5
      flst_real(   5,   7)=          -5
      flst_real(   6,   7)=           0
 
      flst_real(   1,   8)=           2
      flst_real(   2,   8)=           0
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=           5
      flst_real(   5,   8)=          -5
      flst_real(   6,   8)=           2
 
      flst_real(   1,   9)=          -4
      flst_real(   2,   9)=           4
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=           5
      flst_real(   5,   9)=          -5
      flst_real(   6,   9)=           0
 
      flst_real(   1,  10)=          -4
      flst_real(   2,  10)=           0
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=           5
      flst_real(   5,  10)=          -5
      flst_real(   6,  10)=          -4
 
      flst_real(   1,  11)=           4
      flst_real(   2,  11)=          -4
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=           5
      flst_real(   5,  11)=          -5
      flst_real(   6,  11)=           0
 
      flst_real(   1,  12)=           4
      flst_real(   2,  12)=           0
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=           5
      flst_real(   5,  12)=          -5
      flst_real(   6,  12)=           4
 
      flst_real(   1,  13)=          -3
      flst_real(   2,  13)=           3
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=           5
      flst_real(   5,  13)=          -5
      flst_real(   6,  13)=           0
 
      flst_real(   1,  14)=          -3
      flst_real(   2,  14)=           0
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=           5
      flst_real(   5,  14)=          -5
      flst_real(   6,  14)=          -3
 
      flst_real(   1,  15)=           3
      flst_real(   2,  15)=          -3
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=           5
      flst_real(   5,  15)=          -5
      flst_real(   6,  15)=           0
 
      flst_real(   1,  16)=           3
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=           5
      flst_real(   5,  16)=          -5
      flst_real(   6,  16)=           3
 
      flst_real(   1,  17)=           0
      flst_real(   2,  17)=          -1
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=           5
      flst_real(   5,  17)=          -5
      flst_real(   6,  17)=          -1
 
      flst_real(   1,  18)=           0
      flst_real(   2,  18)=           1
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=           5
      flst_real(   5,  18)=          -5
      flst_real(   6,  18)=           1
 
      flst_real(   1,  19)=           0
      flst_real(   2,  19)=          -2
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=           5
      flst_real(   5,  19)=          -5
      flst_real(   6,  19)=          -2
 
      flst_real(   1,  20)=           0
      flst_real(   2,  20)=           2
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=           5
      flst_real(   5,  20)=          -5
      flst_real(   6,  20)=           2
 
      flst_real(   1,  21)=           0
      flst_real(   2,  21)=          -4
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=           5
      flst_real(   5,  21)=          -5
      flst_real(   6,  21)=          -4
 
      flst_real(   1,  22)=           0
      flst_real(   2,  22)=           4
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=           5
      flst_real(   5,  22)=          -5
      flst_real(   6,  22)=           4
 
      flst_real(   1,  23)=           0
      flst_real(   2,  23)=          -3
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=           5
      flst_real(   5,  23)=          -5
      flst_real(   6,  23)=          -3
 
      flst_real(   1,  24)=           0
      flst_real(   2,  24)=           3
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=           5
      flst_real(   5,  24)=          -5
      flst_real(   6,  24)=           3
 
      flst_real(   1,  25)=           0
      flst_real(   2,  25)=           0
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=           5
      flst_real(   5,  25)=          -5
      flst_real(   6,  25)=           0
 
      flst_nreal=          25

      return
      if (nf_loc.eq.5) then
 
      flst_real(   1,   26)=          -5
      flst_real(   2,   26)=           5
      flst_real(   3,   26)=          25
      flst_real(   4,   26)=           5
      flst_real(   5,   26)=          -5
      flst_real(   6,   26)=           0

      flst_real(   1,   27)=          -5
      flst_real(   2,   27)=           0
      flst_real(   3,   27)=          25
      flst_real(   4,   27)=           5
      flst_real(   5,   27)=          -5
      flst_real(   6,   27)=          -5

      flst_real(   1,  28)=           0
      flst_real(   2,  28)=           5
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=           5
      flst_real(   5,  28)=          -5
      flst_real(   6,  28)=           5

      flst_real(   1,   29)=           5
      flst_real(   2,   29)=          -5
      flst_real(   3,   29)=          25
      flst_real(   4,   29)=           5
      flst_real(   5,   29)=          -5
      flst_real(   6,   29)=           0
 
      flst_real(   1,   30)=           5
      flst_real(   2,   30)=           0
      flst_real(   3,   30)=          25
      flst_real(   4,   30)=           5
      flst_real(   5,   30)=          -5
      flst_real(   6,   30)=           5
 
      flst_real(   1,  31)=           0
      flst_real(   2,  31)=          -5
      flst_real(   3,  31)=          25
      flst_real(   4,  31)=           5
      flst_real(   5,  31)=          -5
      flst_real(   6,  31)=          -5

      flst_nreal=          31

      endif !nf=5   
 
      return
      end
 
