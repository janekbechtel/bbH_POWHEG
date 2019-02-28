c adapted from analysis of Wbbj directory

c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer i
      real * 8 dy,dpt,dr,dm
      real*8 ptmin,ptmax,ymin,ymax,rmin,rmax,mmin,mmax
      real*8 ptzoom_max,dpt_zoom
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      integer maxjet
      parameter (maxjet=3)
      integer nptmin
      parameter (nptmin=1)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)
      data cptmin/  '-025'/
      data ptminarr/   25d0 /
      common/infohist/ptminarr,cnum,cptmin
      save /infohist/
      real * 8 powheginput
      external powheginput

      call inihists

      ymin = -8d0
      ymax =  8d0
      dy = 0.25d0

      ptmin = 0d0
      ptmax = 800d0
      dpt = 5d0

      ptzoom_max = 100d0
      dpt_zoom = 1d0

      rmin = 0d0
      rmax = 10d0
      dr = 0.25d0

      mmin = 0d0
      mmax = 1000d0
      dm = 10d0

      call bookupeqbins('sigtot(no-cuts)',1d0,0.5d0,1.5d0)


      do i=1,nptmin

cccc H+1bjet analysis:
      call bookupeqbins('sig(1bj-cuts)'//cptmin(i),1d0,0.5d0,1.5d0)

      call bookupeqbins('H-y(1bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('H-eta(1bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('H-pt(1bj-cuts)'//cptmin(i),dpt,ptmin,ptmax)
      call bookupeqbins('H-m(1bj-cuts)'//cptmin(i),1d0,120d0,130d0)

      call bookupeqbins('b1-y(1bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('b1-eta(1bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('b1-pt(1bj-cuts)'//cptmin(i),dpt,ptmin,ptmax)
      call bookupeqbins('b1-ptzoom(1bj-cuts)'//cptmin(i),
     &     dpt_zoom,ptmin,ptzoom_max)

      call bookupeqbins('j1-y(1bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('j1-eta(1bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('j1-pt(1bj-cuts)'//cptmin(i),dpt,ptmin,ptmax)
      call bookupeqbins('j1-ptzoom(1bj-cuts)'//cptmin(i),
     &    dpt_zoom,ptmin,ptzoom_max)

      call bookupeqbins('Rb1j(1bj-cuts)'//cptmin(i),dr,rmin,rmax)


cccc H+2bjet analysis:
 
      call bookupeqbins('sig(2bj-cuts)'//cptmin(i),1d0,0.5d0,1.5d0)

      call bookupeqbins('H-y(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('H-eta(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('H-pt(2bj-cuts)'//cptmin(i),dpt,ptmin,ptmax)
      call bookupeqbins('H-m(2bj-cuts)'//cptmin(i),1d0,120d0,130d0)

      call bookupeqbins('b1-y(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('b1-eta(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('b1-pt(2bj-cuts)'//cptmin(i),dpt,ptmin,ptmax)
      call bookupeqbins('b1-ptzoom(2bj-cuts)'//cptmin(i),
     &     dpt_zoom,ptmin,ptzoom_max)

      call bookupeqbins('b2-y(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('b2-eta(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('b2-pt(2bj-cuts)'//cptmin(i),dpt,ptmin,400d0)
      call bookupeqbins('b2-ptzoom(2bj-cuts)'//cptmin(i),2d0,1d0,151d0)

      call bookupeqbins('Rb1b2(2bj-cuts)'//cptmin(i),dr,rmin,rmax) 
      call bookupeqbins('Mb1b2(2bj-cuts)'//cptmin(i),dm,mmin,mmax) 


      call bookupeqbins('j1-y(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('j1-eta(2bj-cuts)'//cptmin(i),dy,ymin,ymax)
      call bookupeqbins('j1-pt(2bj-cuts)'//cptmin(i),dpt,ptmin,ptmax)
      call bookupeqbins('j1-ptzoom(2bj-cuts)'//cptmin(i),
     &    dpt_zoom,ptmin,ptzoom_max)

      call bookupeqbins('Rb1j(2bj-cuts)'//cptmin(i),dr,rmin,rmax)
      call bookupeqbins('Rb2j(2bj-cuts)'//cptmin(i),dr,rmin,rmax) 

      enddo
      end
     
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_math.h' 
      include 'pwhg_weights.h'
      integer isthep_loc(NMXHEP)  ! local copy of isthep
      logical ini
      data ini/.true./
      save ini
      integer   maxjet,mjets,njets,numjets,ntracks
      parameter (maxjet=2048)
      real * 8 pj(4,maxjet)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),ihep_of_track(maxtrack)
      character * 1 cnum(9)
      integer nptmin
      parameter (nptmin=1)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)      
      real * 8 ptb1min,ptb2min,ybmax,etabmax,yjmax
      common/infohist/ptminarr,cnum,cptmin
      save /infohist/
      integer j,i,k
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      real * 8 pw(4)
      real * 8 y,eta,pt,m
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp
      integer maxnumlep
      parameter (maxnumlep=10)
      integer mu
      real * 8 ptminfastjet,R
      real * 8 palg
      real * 8 dsig(7)
      integer nweights
      logical inimulti
      data inimulti/.true./
      save inimulti
      logical found_hardjet,found_nexthardjet,found_bjet1,found_bjet2
      real * 8 phardjet(4),pnexthardjet(4),pbjet1(4),pbjet2(4),px(4)
      logical is_B_hadron,is_BBAR_hadron
      external is_B_hadron,is_BBAR_hadron
      real * 8 p_b(4,maxnumlep),p_bbar(4,maxnumlep)
      integer nbjet_array(maxjet),
     $     nbbarjet_array(maxjet),jetinfo(maxjet),id,nb,
     $     nbbar,nbjet,nbbarjet,typeb1,typeb2,wrong_bb_sequence
      real * 8 pthardjet, ptbjet1, ptbjet2
      real*8 dy,deta,dphi,dr
      real*8 y_b1,eta_b1,pt_b1,m_b1
      real*8 y_b2,eta_b2,pt_b2,m_b2
      real*8 y_j1,eta_j1,pt_j1,m_j1
      real*8 pbjet1s(4)
      real*8 pbsum(4),mass_bb
      integer dummy
      integer ih,nhiggs
      real*8 ph(4)
      logical write_cuts
      parameter (write_cuts=.true.)
      logical pass_cuts1,pass_cuts2
      logical pass_cuts

      if(inimulti) then
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         inimulti=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
         nweights=1
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
          nweights=weights_num
      endif

      if(sum(abs(dsig)).eq.0) return


C     ---------------------------------------------------------------
C     Cut parameters
C     ---------------------------------------------------------------
c
c general jet parameters:
      r = 0.5d0
      ptminfastjet = 25d0

c b-jet parameters:
      ptb1min=25d0
      ptb2min=25d0
      etabmax=2.5d0

      yjmax=4.5d0

      pass_cuts1 = .true.
      pass_cuts2 = .true.

      pass_cuts = .true.

      if (ini) then
      if(write_cuts) then
         write(*,*) 'start analysis for WHCPRG =',WHCPRG
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) '                ANALYSIS CUTS               '
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) ''
         write(*,*) 'parameters for jet algorithm:'
         write(*,*) 'pt_min = ',ptminfastjet
         write(*,*) 'R = ',r       
         write(*,*) ''  
         write(*,*) 'jet parameters for analysis'   
         write(*,*) 'pt_min(hard jet) = ',ptminarr
         write(*,*) 'abs(yj) < ',yjmax
         write(*,*) ''  
         write(*,*) 'b-JETS:'
         write(*,*) 'pt(b1) > ',ptb1min
         write(*,*) '2-bjet-analysis: pt(b2) > ',ptb2min
         write(*,*) 'abs(etab) < ',etabmax
         write(*,*) ''
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
      endif
      
      ini=.false.
      endif

c init:
      ph = 0d0
      p_b = 0d0
      p_bbar = 0d0
      pj = 0d0

      ptrack = 0d0
      phardjet = 0d0
      pnexthardjet = 0d0

      pbjet1 = 0d0
      pbjet2 = 0d0


      do ihep=1,nhep  
         isthep_loc(ihep) = isthep(ihep)
      enddo
      
      if ((WHCPRG.eq.'NLO   ').or.(WHCPRG.eq.'LHE   ')) then
         nb=0
         nbbar=0
         nhiggs = 0
         do ihep=1,nhep            
            if(idhep(ihep).eq.25) then !H
               ph=phep(1:4,ihep)
               ih=ihep
               nhiggs=nhiggs+1
            elseif(is_B_hadron(idhep(ihep))) then
               nb=nb+1
               do mu=1,4
                  p_b(mu,nb)=phep(mu,ihep)
               enddo     
            elseif(is_BBAR_hadron(idhep(ihep))) then
               nbbar=nbbar+1
               do mu=1,4
                  p_bbar(mu,nbbar)=phep(mu,ihep)
               enddo              
            endif
         enddo !ihep

      elseif((WHCPRG.eq.'HERWIG').or.(WHCPRG.eq.'PYTHIA')) then
         nb=0
         nbbar=0
         nhiggs = 0
         do ihep=1,nhep
            if(idhep(ihep).eq.25.and.isthep_loc(ihep).eq.1) then
               ph=phep(1:4,ihep)
               ih=ihep
               nhiggs = nhiggs + 1
            endif
            if (isthep_loc(ihep).eq.1.and.(is_B_hadron(idhep(ihep))))
     1           then  
               nb=nb+1
               do mu=1,4
                  p_b(mu,nb)=phep(mu,ihep)
               enddo  
            endif
            if (isthep_loc(ihep).eq.1.and.(is_BBAR_hadron(idhep(ihep))))
     1           then  
               nbbar=nbbar+1
               do mu=1,4
                  p_bbar(mu,nbbar)=phep(mu,ihep)
               enddo  
            endif
         enddo !ihep
      endif ! NLO / LHEF/ PYTHIA
      
c check that exactly 1 Higgs is selected
      if (nhiggs.ne.1) then
         write(*,*) "Error in pwhg_analysis: ",nhiggs,"Higgses"
         call printleshouches
         call exit(1)
      endif
 

      if (nb*nbbar.eq.0) then
         if (((nb.eq.0).and.(nbbar.ne.0)).or.
     $       ((nb.ne.0).and.(nbbar.eq.0))) then
            write(*,*) 'SEVERE WARNING: ***** One b is missing ******' 
            return
         endif
      endif
        
      ntracks=0
      mjets=0
      ihep_of_track = 10000
c     Loop over final state particles to find jets, exclude Higgs! 
      do ihep=1,nhep
         if ((isthep_loc(ihep).eq.1).and.(idhep(ihep).ne.25)
     &        ) then
           if (ntracks.eq.maxtrack) then
              write(*,*) 'Too many particles. Increase maxtrack.'//
     #             ' PROGRAM ABORTS'
              call exit(1)
           endif
c     copy momenta to construct jets 
           ntracks=ntracks+1
           ihep_of_track(ntracks)=ihep
           do mu=1,4
              ptrack(mu,ntracks)=phep(mu,ihep)
           enddo
        endif
      enddo !ihep


      jetvec = 0
      numjets = 0
      if (ntracks.eq.0) then
         numjets=0
      else
c     palg=1 is standard kt, -1 is antikt
         palg = -1d0
         call fastjetppgenkt(ptrack,ntracks,R,palg,ptminfastjet,
     $        pj,numjets,jetvec)
c         call fastjetktwhich(ptrack,ntracks,ptminfastjet,R,
c     $        pj,mjets,jetvec) 

      endif
      
c     find in which ptrack the B hadrons ended up
      nbjet_array = 0
      nbbarjet_array = 0
      nbjet=0
      nbbarjet=0
c     loop over tracks
      do i=1,ntracks
         id=idhep(ihep_of_track(i))
         if (is_B_hadron(id)) then
            nbjet=nbjet+1
            nbjet_array(nbjet)=jetvec(i)            
         elseif (is_BBAR_hadron(id)) then   
            nbbarjet=nbbarjet+1
            nbbarjet_array(nbbarjet)=jetvec(i)                        
         endif
      enddo

c     we want at least two(one) b-type jets for the 2bjet (1bjet) analysis:
      if (nbjet+nbbarjet.lt.2) pass_cuts2 = .false. 
      if (nbjet+nbbarjet.lt.1) pass_cuts1 = .false. 

c     jets are ordered in decreasing pt. Set up array of info on jets
c     if jetinfo=0 then non-b jet
c     if jetinfo=5 then b jet
c     if jetinfo=-5 then bbar jet
c      do i=1,numjets
c         jetinfo(i)=0
c      enddo
      
      do i=1,numjets
         jetinfo(i)=0
         do j=1,nbjet
            if (i.eq.nbjet_array(j)) then
               jetinfo(i)=5
            endif
         enddo
         do j=1,nbbarjet
            if (i.eq.nbbarjet_array(j)) then
               jetinfo(i)=-5
             endif
         enddo
      enddo
       

      found_hardjet=.false.
      found_nexthardjet=.false.
      found_bjet1=.false.
      found_bjet2=.false.
      typeb1=0
      typeb2=0
      do i=1,numjets
         if (jetinfo(i).eq.0) then
            if (.not.found_hardjet) then
               found_hardjet=.true.
               do mu=1,4
                  phardjet(mu)=pj(mu,i)
               enddo            
               goto 111
            elseif (found_hardjet.and..not.found_nexthardjet) then
               found_nexthardjet=.true.
               do mu=1,4
                  pnexthardjet(mu)=pj(mu,i)
               enddo      
               goto 111
            endif
         elseif (abs(jetinfo(i)).eq.5) then
            if (.not.found_bjet1) then
               found_bjet1=.true.
               do mu=1,4
                  pbjet1(mu)=pj(mu,i)
               enddo 
c     keep track of b flavor of the 1st jet (quark or antiquark)
               typeb1=jetinfo(i)
               goto 111
            endif
            if (.not.found_bjet2.and.found_bjet1) then
               found_bjet2=.true.
               do mu=1,4
                  pbjet2(mu)=pj(mu,i)
               enddo 
c     keep track of b flavor of the 2nd jet (quark or antiquark)
               typeb2=jetinfo(i)
            endif
            if (found_bjet1.and.found_bjet2) then
c     they must come from a b-bbar couple. Otherwise, return warning:
               if (typeb1*typeb2.gt.0) then
                  wrong_bb_sequence = wrong_bb_sequence + 1                  
                  if (mod(wrong_bb_sequence,1).eq.0) then
                     write(*,*) 'WARNING: 2 b or 2 bbar in sequence ', 
     $                    wrong_bb_sequence
                  endif
c                  return ! return, if you require b-bbar couple
               endif
            endif
         endif
 111     continue
      enddo

      call filld('sigtot(no-cuts)',1d0,dsig)         

c if there are not one (two) bjets, set cuts1 (cuts2) to false:
      if (.not.(found_bjet1.and.found_bjet2)) then
         pass_cuts2 = .false.
      endif
      if (.not.(found_bjet1.or.found_bjet2)) then
         pass_cuts1 = .false.
         return
      endif
     
! init:
      pt = 0d0
      eta = 0d0
      y = 0d0
      m = 0d0

      pt_j1 = 0d0
      eta_j1 = 0d0
      y_j1 = 0d0
      m_j1 = 0d0

      pt_b1 = 0d0
      eta_b1 = 0d0
      y_b1 = 0d0
      m_b1 = 0d0

      pt_b2 = 0d0
      eta_b2 = 0d0
      y_b2 = 0d0
      m_b2 = 0d0

      dy = 0d0
      deta = 0d0
      dphi= 0d0
      dr = 0d0

      pbsum = 0d0
      mass_bb = 0d0

      if (found_bjet1) call getyetaptmass(pbjet1,y_b1,eta_b1,pt_b1,m_b1)      
      if ((pt_b1.lt.ptb1min).or.(abs(eta_b1).gt.etabmax)) then 
         pass_cuts2 = .false.      
         if (found_bjet2) 
     &   call getyetaptmass(pbjet2,y_b2,eta_b2,pt_b2,m_b2)
         if ((pt_b2.lt.ptb2min).or.(abs(eta_b2).gt.etabmax)) then 
            pass_cuts1 = .false.
            return
         else
            pbjet1s = pbjet2
         endif
      else ! bjet1 passes cuts
         pbjet1s = pbjet1     
         if (found_bjet2) 
     &   call getyetaptmass(pbjet2,y_b2,eta_b2,pt_b2,m_b2)
         if ((pt_b2.lt.ptb2min).or.(abs(eta_b2).gt.etabmax)) then 
            pass_cuts2 = .false.
         endif
      endif    

      if (.not.pass_cuts1) return

c*****************************************************
      
      do i=1,nptmin        

         call filld('sig(1bj-cuts)'//cptmin(i),1d0,dsig)         

cccc H+1jet analysis:
c
c     Higgs
         call getyetaptmass(ph,y,eta,pt,m)
         call filld('H-y(1bj-cuts)'//cptmin(i),    y, dsig)
         call filld('H-eta(1bj-cuts)'//cptmin(i),eta, dsig)
         call filld('H-pt(1bj-cuts)'//cptmin(i),  pt, dsig)
         call filld('H-m(1bj-cuts)'//cptmin(i), m, dsig)

c     hardest b jet
         call getyetaptmass(pbjet1s,y_b1,eta_b1,pt_b1,m_b1)
         call filld('b1-y(1bj-cuts)'//cptmin(i),     y_b1, dsig)
         call filld('b1-eta(1bj-cuts)'//cptmin(i), eta_b1, dsig)
         call filld('b1-pt(1bj-cuts)'//cptmin(i),   pt_b1, dsig)
         call filld('b1-ptzoom(1bj-cuts)'//cptmin(i),pt_b1, dsig)

c     hardest jet         
         if (found_hardjet) then
            call getyetaptmass(phardjet,y_j1,eta_j1,pt_j1,m_j1)

            if ((pt_j1.gt.ptminarr(i)).and.(abs(y_j1).lt.yjmax)) then
               call filld('j1-y(1bj-cuts)'//cptmin(i),     y_j1, dsig)
               call filld('j1-eta(1bj-cuts)'//cptmin(i), eta_j1, dsig)
               call filld('j1-pt(1bj-cuts)'//cptmin(i),   pt_j1, dsig)
               call filld('j1-ptzoom(1bj-cuts)'//cptmin(i), pt_j1, dsig)

c bj separation:
               call getdydetadphidr(pbjet1s,phardjet,dy,deta,dphi,dr)
               call filld('Rb1j(1bj-cuts)'//cptmin(i),     dr, dsig)

            elseif (found_nexthardjet) then  
               found_hardjet = .false.
               call getyetaptmass(pnexthardjet,y_j1,eta_j1,pt_j1,m_j1)
               if ((pt_j1.gt.ptminarr(i)).and.(abs(y_j1).lt.yjmax)) then
                  found_nexthardjet = .false.
                 call filld('j1-y(1bj-cuts)'//cptmin(i),     y_j1, dsig)
                 call filld('j1-eta(1bj-cuts)'//cptmin(i), eta_j1, dsig)
                 call filld('j1-pt(1bj-cuts)'//cptmin(i),   pt_j1, dsig)
                 call filld('j1-ptzoom(1bj-cuts)'//cptmin(i),pt_j1,dsig)
c bj separation:
                  call getdydetadphidr(
     &                 pbjet1s,pnexthardjet,dy,deta,dphi,dr)
                  call filld('Rb1j(1bj-cuts)'//cptmin(i),     dr, dsig)
               endif
            endif !pt,y
         endif !hard jet

      

         if (.not.pass_cuts2) goto 123

cccc H+2jet analysis:
c   
         call filld('sig(2bj-cuts)'//cptmin(i),1d0,dsig)         

c     Higgs
         call getyetaptmass(ph,y,eta,pt,m)
         call filld('H-y(2bj-cuts)'//cptmin(i),    y, dsig)
         call filld('H-eta(2bj-cuts)'//cptmin(i),eta, dsig)
         call filld('H-pt(2bj-cuts)'//cptmin(i),  pt, dsig)
         call filld('H-m(2bj-cuts)'//cptmin(i), m, dsig)

c     hardest b jet
         call getyetaptmass(pbjet1,y_b1,eta_b1,pt_b1,m_b1)
         call filld('b1-y(2bj-cuts)'//cptmin(i),     y_b1, dsig)
         call filld('b1-eta(2bj-cuts)'//cptmin(i), eta_b1, dsig)
         call filld('b1-pt(2bj-cuts)'//cptmin(i),   pt_b1, dsig)
         call filld('b1-ptzoom(2bj-cuts)'//cptmin(i),pt_b1, dsig)
        
c     next-to-hardest b jet
         call getyetaptmass(pbjet2,y_b2,eta_b2,pt_b2,m_b2)
         call filld('b2-y(2bj-cuts)'//cptmin(i),     y_b2, dsig)
         call filld('b2-eta(2bj-cuts)'//cptmin(i), eta_b2, dsig)
         call filld('b2-pt(2bj-cuts)'//cptmin(i),   pt_b2, dsig)
         call filld('b2-ptzoom(2bj-cuts)'//cptmin(i),   pt_b2, dsig)


c bb separation:
         call getdydetadphidr(pbjet1,pbjet2,dy,deta,dphi,dr)
         call filld('Rb1b2(2bj-cuts)'//cptmin(i),     dr, dsig)

         pbsum = pbjet1+pbjet2
         call pwhg_getinvmass(pbsum,mass_bb)
         call filld('Mb1b2(2bj-cuts)'//cptmin(i),mass_bb, dsig)

         dy = 0d0
         deta = 0d0
         dphi= 0d0
         dr = 0d0

c     hardest jet         
         if (found_hardjet) then
            call getyetaptmass(phardjet,y_j1,eta_j1,pt_j1,m_j1)

            if ((pt_j1.gt.ptminarr(i)).and.(abs(y_j1).lt.yjmax)) then
               call filld('j1-y(2bj-cuts)'//cptmin(i),     y_j1, dsig)
               call filld('j1-eta(2bj-cuts)'//cptmin(i), eta_j1, dsig)
               call filld('j1-pt(2bj-cuts)'//cptmin(i),   pt_j1, dsig)
               call filld('j1-ptzoom(2bj-cuts)'//cptmin(i),pt_j1, dsig)

c bj separation:
               call getdydetadphidr(pbjet1,phardjet,dy,deta,dphi,dr)
               call filld('Rb1j(2bj-cuts)'//cptmin(i),     dr, dsig)
               call getdydetadphidr(pbjet2,phardjet,dy,deta,dphi,dr)
               call filld('Rb2j(2bj-cuts)'//cptmin(i),     dr, dsig)

            elseif (found_nexthardjet) then  
               found_hardjet = .false.
               call getyetaptmass(pnexthardjet,y_j1,eta_j1,pt_j1,m_j1)
               if ((pt_j1.gt.ptminarr(i)).and.(abs(y_j1).lt.yjmax)) then
                  found_nexthardjet = .false.
                 call filld('j1-y(2bj-cuts)'//cptmin(i),     y_j1, dsig)
                 call filld('j1-eta(2bj-cuts)'//cptmin(i), eta_j1, dsig)
                 call filld('j1-pt(2bj-cuts)'//cptmin(i),   pt_j1, dsig)
                 call filld('j1-ptzoom(2bj-cuts)'//cptmin(i),pt_j1,dsig)
c bj separation:
                  call getdydetadphidr(
     &                 pbjet1,pnexthardjet,dy,deta,dphi,dr)
                  call filld('Rb1j(2bj-cuts)'//cptmin(i),     dr, dsig)
                  call getdydetadphidr(
     &                 pbjet2,pnexthardjet,dy,deta,dphi,dr)
                  call filld('Rb2j(2bj-cuts)'//cptmin(i),     dr, dsig)

               endif
            endif !pt,y         
         endif !hard jet

 123     continue

      enddo
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass
      call pwhg_getrapidity(p,y)      
      pt=sqrt(p(1)**2+p(2)**2)
      call pwhg_getpseudorapidity(p,eta)
      call pwhg_getinvmass(p,mass)
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end


      function is_B_hadron(id)
      implicit none
      logical is_B_hadron
      integer id
      is_B_hadron=((id.gt.-600).and.(id.lt.-500)).or.
     $     ((id.gt.5000).and.(id.lt.6000)).or.(id.eq.5)
      end

      function is_BBAR_hadron(id)
      implicit none
      logical is_BBAR_hadron
      integer id
      is_BBAR_hadron=((id.gt.500).and.(id.lt.600)).or.
     $     ((id.gt.-6000).and.(id.lt.-5000)).or.(id.eq.-5)
      end

