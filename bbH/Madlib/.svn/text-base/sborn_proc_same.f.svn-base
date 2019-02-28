c
c this routine is adpated from sborn_proc.f : 
c     only genuinely different contributions are retained
c     in sborn_proc and born_color routines
c
      subroutine sborn_proc(p_born,legs,wgt,wgtjk,wgtmunu)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision wgt
      double precision p_born(0:3,nexternal-1),wgt2(nexternal-1),
     &   wgtmunu(0:3,0:3,nexternal-1),wgtjk(nexternal-1,nexternal-1)
      double precision p_born1(0:3,nexternal-1),wgt1(0:nexternal-1),
     &   wgtmunu1(0:3,0:3,nexternal-1),wgtjk1(nexternal-1,nexternal-1)
      integer legs(nexternal-1),lstr,i
      character*140 str
      integer ic(nexternal-1),legs1(nexternal-1)
      logical mtc,even
      logical calculatedBorn
      integer skip
      common/cBorn/calculatedBorn,skip
 
      calculatedBorn=.false.
      
      do i=1,nexternal-1
         ic(i)=i
      enddo
      mtc=.false.
 10   call nexper(nexternal-1- 3,ic( 3+1),mtc,even)
      do i= 3+1,nexternal-1
         ic(i)=ic(i)+ 3
      enddo
      CALL SWITCHMOM (P_born,P_born1,IC,NEXTERNAL-1)
      CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL-1)
      
      call convert_to_string(nexternal-1,legs1,str,lstr)

      if (legs(1).lt.0) then !qbq->Hbb
         call sborn_cl_001(p_born1,wgtmunu1,wgt1)
         call sborn_sf_001(p_born1,wgtjk1)
         goto 20
      elseif (legs(1).gt.0) then !qqb->Hbb
         call sborn_cl_002(p_born1,wgtmunu1,wgt1)
         call sborn_sf_002(p_born1,wgtjk1)
         goto 20
      else ! gg->Hbb
         call sborn_cl_009(p_born1,wgtmunu1,wgt1)
         call sborn_sf_009(p_born1,wgtjk1)
         goto 20
      endif    
      
      do while(mtc)
         do i= 3+1,nexternal-1
            ic(i)=ic(i)- 3
         enddo
         goto 10
      enddo
      if (.not.mtc) then
         write (*,*) "Error #1, in sborn_proc.f"
         stop
      endif
      
 20   wgt=0d0
      call switchborns(wgt1(1),wgt2,wgtjk1,wgtjk,wgtmunu1,wgtmunu,
     &     ic,nexternal-1)
      do i=1,nexternal-1
         if(wgt.eq.0d0 .and. wgt2(i).ne.0d0) then
            wgt=wgt2(i)
         elseif (wgt.ne.0d0 .and. wgt2(i).ne.0d0 .and.
     &           abs((wgt-wgt2(i))/wgt).gt.1d-7 ) then
            write (*,*) "Error #2 in sborn_proc ",i,wgt2
            stop
         endif
      enddo
      
      end
      
      
      
      
      
      subroutine born_color(legs,color)
      implicit none
      include "nexternal.inc"
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2001(maxamps), jamp2001(0:maxamps)
      common/to_amps_001/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxamps)
      common/to_amps_002/amp2002,jamp2002
      Double Precision amp2009(maxamps), jamp2009(0:maxamps)
      common/to_amps_009/amp2009,jamp2009
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal-1,maxamps)
      integer color(2,nexternal-1),color1(2,nexternal-1)
      double precision random,xtarget
      external random
      integer legs(nexternal-1),lstr,i,j
      character*140 str
      integer ic(nexternal-1),legs1(nexternal-1)
      integer iflow,ifl
      logical mtc,even
      
      do i=1,nexternal-1
         ic(i)=i
      enddo
      mtc=.false.
 10   call nexper(nexternal-1- 3,ic( 3+1),mtc,even)
      do i= 3+1,nexternal-1
         ic(i)=ic(i)+ 3
      enddo
      CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL-1)
      
      call convert_to_string(nexternal-1,legs1,str,lstr)

      if (legs(1).lt.0) then !qbq->Hbb
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif (legs(1).gt.0) then !qqb->Hbb
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      else ! gg->Hbb
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      endif    
      
      do while(mtc)
         do i= 3+1,nexternal-1
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
         do j=1,nexternal-1
            color1(i,j)=ICOLUP(i,j,ifl)
         enddo
      enddo
      call switchcolor(color1,color,
     &     ic,nexternal-1)
      
      end
      
      
      
      
      subroutine convert_to_string(npart,id,string,lstring)
      implicit none
      integer npart,lstring,i
      integer id(npart)
      character*140 string
      character*3 s
      
      do i=1,140
         string(i:i)=' '
      enddo
      lstring=0
      do i=1,npart
         if (id(i).eq.21) id(i)=0
         if (abs(id(i)).le.9) then
            s=char(abs(id(i))+48)
         elseif(abs(id(i)).le.99)then
            s=char(abs(id(i))/10+48)
     &           //char(mod(abs(id(i)),10)+48)
               elseif(abs(id(i)).le.999) then
                  s=char(abs(id(i))/100+48)
     &           //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &           //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         else
            write (*,*) 'error, particle ID is too large',abs(id(i))
         endif
         if (id(i).ge.0) then
            if (id(i).le.9) then
               string(lstring+1:lstring+1)=s
               lstring=lstring+1
            elseif (id(i).le.99) then
               string(lstring+1:lstring+2)=s
               lstring=lstring+2
            elseif (id(i).le.999) then
               string(lstring+1:lstring+3)=s
               lstring=lstring+3
            endif
         else
            if (abs(id(i)).le.9) then
               string(lstring+1:lstring+2)='-'//s
               lstring=lstring+2
            elseif (abs(id(i)).le.99) then
               string(lstring+1:lstring+3)='-'//s
               lstring=lstring+3
            elseif (abs(id(i)).le.999) then
               string(lstring+1:lstring+4)='-'//s
               lstring=lstring+4
            endif
         endif
      enddo
      end
