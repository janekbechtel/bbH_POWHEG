      subroutine sreal_proc(p,legs,wgt)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision p(0:3,nexternal),wgt
      integer legs(nexternal),lstr
      character*140 str
      double precision P1(0:3,nexternal)
      integer i,ic(nexternal),legs1(nexternal)
      logical mtc,even
      
      do i=1,nexternal
         ic(i)=i
      enddo
      mtc=.false.
 10   call nexper(nexternal- 3,ic( 3+1),mtc,even)
      do i= 3+1,nexternal
         ic(i)=ic(i)+ 3
      enddo
      CALL SWITCHMOM (P,P1,IC,NEXTERNAL)
      CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL)
      
      call convert_to_string(nexternal,legs1,str,lstr)

      
      ! qbq->Hbbg:
      if ((legs(1).lt.0).and.(legs(2).gt.0).and.(legs(6).eq.0)) then 
         call srealmtrx_001(p1,wgt)         
         goto 20
      ! qqb->Hbbg:
      elseif ((legs(1).gt.0).and.(legs(2).lt.0).and.(legs(6).eq.0)) then 
         call srealmtrx_003(p1,wgt)         
         goto 20
      ! qbg->Hbbqb:
      elseif ((legs(1).lt.0).and.(legs(2).eq.0).and.(legs(6).lt.0)) then 
         call srealmtrx_002(p1,wgt)         
         goto 20
      ! qg->Hbbq:
      elseif ((legs(1).gt.0).and.(legs(2).eq.0).and.(legs(6).gt.0)) then 
         call srealmtrx_004(p1,wgt)         
         goto 20
      ! gqb->Hbbqb:
      elseif ((legs(1).eq.0).and.(legs(2).lt.0).and.(legs(6).lt.0)) then 
         call srealmtrx_017(p1,wgt)         
         goto 20
      ! gq->Hbbq:
      elseif ((legs(1).eq.0).and.(legs(2).gt.0).and.(legs(6).gt.0)) then 
         call srealmtrx_018(p1,wgt)         
         goto 20
      ! gg->Hbbg:
      elseif ((legs(1).eq.0).and.(legs(2).eq.0).and.(legs(6).eq.0)) then 
         call srealmtrx_025(p1,wgt)        
         goto 20
      endif   
      
      
      do while(mtc)
         do i= 3+1,nexternal
            ic(i)=ic(i)- 3
         enddo
         goto 10
      enddo
      if (.not.mtc) then
         write (*,*) "Error #1, in sreal_proc.f"
         stop
      endif
      
 20   continue
      return
      end
      
      
      subroutine real_color(legs,color)
      implicit none
      include "nexternal.inc"
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2001(maxamps), jamp2001(0:maxamps)
      common/to_Ramps_001/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxamps)
      common/to_Ramps_002/amp2002,jamp2002
      Double Precision amp2003(maxamps), jamp2003(0:maxamps)
      common/to_Ramps_003/amp2003,jamp2003
      Double Precision amp2004(maxamps), jamp2004(0:maxamps)
      common/to_Ramps_004/amp2004,jamp2004

      Double Precision amp2017(maxamps), jamp2017(0:maxamps)
      common/to_Ramps_017/amp2017,jamp2017
      Double Precision amp2018(maxamps), jamp2018(0:maxamps)
      common/to_Ramps_018/amp2018,jamp2018

      Double Precision amp2025(maxamps), jamp2025(0:maxamps)
      common/to_Ramps_025/amp2025,jamp2025
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal,maxamps)
      integer color(2,nexternal),color1(2,nexternal)
      double precision random,xtarget
      external random
      integer legs(nexternal),lstr,i,j
      character*140 str
      integer ic(nexternal),legs1(nexternal)
      integer iflow,ifl
      logical mtc,even
      
      do i=1,nexternal
         ic(i)=i
      enddo
      mtc=.false.
 10   call nexper(nexternal- 3,ic( 3+1),mtc,even)
      do i= 3+1,nexternal
         ic(i)=ic(i)+ 3
      enddo
      CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL)
      
      call convert_to_string(nexternal,legs1,str,lstr)
      
      ! qbq->Hbbg:
      if ((legs(1).lt.0).and.(legs(2).gt.0).and.(legs(6).eq.0)) then 
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      ! qqb->Hbbg:
      elseif ((legs(1).gt.0).and.(legs(2).lt.0).and.(legs(6).eq.0)) then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo 
         goto 20
      ! qbg->Hbbqb:
      elseif ((legs(1).lt.0).and.(legs(2).eq.0).and.(legs(6).lt.0)) then 
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      ! qg->Hbbq:
      elseif ((legs(1).gt.0).and.(legs(2).eq.0).and.(legs(6).gt.0)) then 
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      ! gqb->Hbbqb:
      elseif ((legs(1).eq.0).and.(legs(2).lt.0).and.(legs(6).lt.0)) then 
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      ! gq->Hbbq:
      elseif ((legs(1).eq.0).and.(legs(2).gt.0).and.(legs(6).gt.0)) then 
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      ! gg->Hbbg:
      elseif ((legs(1).eq.0).and.(legs(2).eq.0).and.(legs(6).eq.0)) then 
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      endif   

      
      do while(mtc)
         do i= 3+1,nexternal
            ic(i)=ic(i)- 3
         enddo
         goto 10
      enddo
      if (.not.mtc) then
         write (*,*) "Error #1, in sborn_proc.f"
         stop
      endif
      
 20   continue
      xtarget=jamp2cum(iflow)*random()
      ifl=1
      do while (jamp2cum(ifl).lt.xtarget)
         ifl=ifl+1
      enddo
      do i=1,2
         do j=1,nexternal
            color1(i,j)=ICOLUP(i,j,ifl)
         enddo
      enddo
      call switchcolor(color1,color,
     &     ic,nexternal)
      
      return
      end
      
      
      
      
