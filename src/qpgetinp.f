      subroutine qpgetinp(unit)
      implicit none
      integer unit
c
      include 'qpglobal.h'
c
c     work space
c
      integer i,j,l,ir,ig,isg,is,is1,flen,iswap,nhypo
      double precision twindow,twinout,suppress,munit,sdfsel
      double precision strike,dip,rake,depdif,dswap(11)
      character*800 grndir,outfile,fswap
      character*800 line
c
c     uniform receiver depth
c     ======================
c
      call skip_comments(unit)
      read(unit,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call skip_comments(unit)
      read(unit,*)twindow,dt
      ntcut=1+idnint(twindow/dt)
      nt=2
100   nt=2*nt
      if(nt.lt.ntcut)goto 100
      nf=nt/2
      if(nf.gt.nfmax)then
        stop 'Error: nfmax defined in qpglobal.h too small!'
      endif
      df=1.d0/(dble(nt)*dt)
c
      call skip_comments(unit)
      read(unit,*)fcut
      nfcut=min0(nf,1+idnint(fcut/df))
      fcut=dble(nfcut-1)*df
      call skip_comments(unit)
      read(unit,*)slwmax
      if(slwmax.le.0.d0)then
        stop ' Error: bad selection of max. slowness!'
      else
        slwmax=slwmax/KM2M
      endif
c
      call skip_comments(unit)
      read(unit,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        fi=0.d0
      else
        fi=dlog(suppress)*df/PI2
      endif
c
c     cutoffs of spectra
c     ==================
c
      call skip_comments(unit)
      read(unit,*)fgr,ldeggr
      if(fgr.lt.0.d0)fgr=0.d0
      if(ldeggr.lt.0)ldeggr=0
      if(fgr.gt.0.d0.and.ldeggr.le.0.or.
     &   fgr.le.0.d0.and.ldeggr.gt.0)then
        stop ' Bad fgr and ldeggr combination!'
      endif
      nogravity=fgr*dble(ldeggr).le.0.d0
c
      call skip_comments(unit)
      read(unit,'(a)') line
      read(line,*)i,j
      selpsv=i.eq.1
      selsh=j.eq.1
      if(.not.(selpsv.or.selsh))then
        stop ' Error: none of PSV and SH is selected!'
      endif
      read(line,*,end=110)i,j,ldegmin,ldegcut
110   ldegcut=min0(ldegmax-ndmax-1,ldegcut)
      if(ldegcut.lt.ldegmin)ldegcut=ldegmin
      ldegmin=min0(max0(1+ndmax,ldegmin),ldegcut,ldegmax)
c
c     Green's function files
c     ======================
c
      call skip_comments(unit)
      read(unit,*)ngrn,rr0,grndir
      if(ngrn.le.0)then
        stop ' bad number of source depths!'
      else if(ngrn.gt.ngrnmax)then
        stop ' number of source depths exceeds the maximum!'
      endif
      rr0=rr0*KM2M
c
      do ig=1,ngrn
        call skip_comments(unit)
        read(unit,*)grndep(ig),grnfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop ' bad Green function selection!'
        endif
        grndep(ig)=grndep(ig)*KM2M
      enddo
c
c     sort green function files by source depth
c
      do i=1,ngrn
        do j=i+1,ngrn
          if(grndep(j).lt.grndep(i))then
            dswap(1)=grndep(i)
            fswap=grnfile(i)
            iswap=grnsel(i)
c
            grndep(i)=grndep(j)
            grnfile(i)=grnfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=dswap(1)
            grnfile(j)=fswap
            grnsel(j)=iswap
          endif
        enddo
      enddo
c
      do flen=800,1,-1
        if(grndir(flen:flen).ne.' ')goto 200
      enddo
200   continue
      do ig=1,ngrn
        rgrnfile(ig)=grndir(1:flen)//'R_'//grnfile(ig)
        tgrnfile(ig)=grndir(1:flen)//'T_'//grnfile(ig)
        pgrnfile(ig)=grndir(1:flen)//'P_'//grnfile(ig)
        sgrnfile(ig)=grndir(1:flen)//'S_'//grnfile(ig)
        ggrnfile(ig)=grndir(1:flen)//'G_'//grnfile(ig)
        vgrnfile(ig)=grndir(1:flen)//'V_'//grnfile(ig)
        wgrnfile(ig)=grndir(1:flen)//'W_'//grnfile(ig)
      enddo
c
c     multi-event source parameters
c     =============================
c
      call skip_comments(unit)
      read(unit,*)ns,sdfsel
      if(ns.gt.nsmax)then
        stop ' Error: too many subevents'
      endif
      if(sdfsel.eq.1)then
        do is=1,ns
          call skip_comments(unit)
c
c         the six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
c
          read(unit,*)munit,mrr(is),mtt(is),mpp(is),
     &                          mrt(is),mpr(is),mtp(is),
     &                    lats(is),lons(is),deps(is),togs(is),trss(is)

          mtt(is)=mtt(is)*munit
          mpp(is)=mpp(is)*munit
          mrr(is)=mrr(is)*munit
          mtp(is)=mtp(is)*munit
          mpr(is)=mpr(is)*munit
          mrt(is)=mrt(is)*munit
          deps(is)=deps(is)*KM2M
        enddo
      else if(sdfsel.eq.2)then
        do is=1,ns
          call skip_comments(unit)
          read(unit,*)munit,strike,dip,rake,
     &                    lats(is),lons(is),deps(is),togs(is),trss(is)
          call moments(munit,strike,dip,rake,
     &                 mtt(is),mpp(is),mrr(is),
     &                 mtp(is),mpr(is),mrt(is))
          deps(is)=deps(is)*KM2M
        enddo
      else
        stop ' bad selection for the source data format!'
      endif
c
      togsmin=togs(1)
      do is=1,ns
        togsmin=dmin1(togsmin,togs(is))
      enddo
      do is=1,ns
        togs(is)=togs(is)-togsmin
      enddo
c
c     sort sub-events by depth
c
      do i=1,ns
        do j=i+1,ns
          if(deps(j).lt.deps(i))then
            dswap(1)=lats(i)
            dswap(2)=lons(i)
            dswap(3)=deps(i)
            dswap(4)=mtt(i)
            dswap(5)=mpp(i)
            dswap(6)=mrr(i)
            dswap(7)=mtp(i)
            dswap(8)=mpr(i)
            dswap(9)=mrt(i)
            dswap(10)=togs(i)
            dswap(11)=trss(i)
c
            lats(i)=lats(j)
            lons(i)=lons(j)
            deps(i)=deps(j)
            mtt(i)=mtt(j)
            mpp(i)=mpp(j)
            mrr(i)=mrr(j)
            mtp(i)=mtp(j)
            mpr(i)=mpr(j)
            mrt(i)=mrt(j)
            togs(i)=togs(j)
            trss(i)=trss(j)
c
            lats(j)=dswap(1)
            lons(j)=dswap(2)
            deps(j)=dswap(3)
            mtt(j)=dswap(4)
            mpp(j)=dswap(5)
            mrr(j)=dswap(6)
            mtp(j)=dswap(7)
            mpr(j)=dswap(8)
            mrt(j)=dswap(9)
            togs(j)=dswap(10)
            trss(j)=dswap(11)
          endif
        enddo
      enddo
c
      isg1(1)=1
      is1=1
      do ig=1,ngrn-1
        depdif=0.5d0*(grndep(ig)+grndep(ig+1))
        isg2(ig)=isg1(ig)-1
        do is=is1,ns
          if(deps(is).lt.depdif)then
            isg2(ig)=is
          endif
        enddo
        isg1(ig+1)=isg2(ig)+1
        is1=isg1(ig+1)
      enddo
      isg2(ngrn)=ns
c
      do ig=1,ngrn
        nsg(ig)=max0(0,1+isg2(ig)-isg1(ig))
      enddo
c
c     receiver parameters
c     ===================
c
      call skip_comments(unit)
      read(unit,*)outfile,ioutform
      if(ioutform.lt.1.or.ioutform.gt.2)then
        stop ' Error: bad selection of output format!'
      endif
      call skip_comments(unit)
      read(unit,*)twinout
      ntcutout=min0(nt,1+idnint(twinout/dt))
      call skip_comments(unit)
      read(unit,*)nlpf,f1corner,f2corner
      call skip_comments(unit)
      read(unit,*)slwlwcut,slwupcut
      if(slwlwcut.gt.slwupcut)then
        slwlwcut=0.d0
        slwupcut=slwmax
      else
        slwlwcut=slwlwcut/KM2M
        slwupcut=slwupcut/KM2M
      endif
      call skip_comments(unit)
      read(unit,*)nr
      if(nr.gt.nrmax)then
        stop ' Error: too many receivers'
      endif
c
      if(ioutform.eq.2)then
        ishypo=1
        do is=2,ns
          if(togs(is).lt.togs(ishypo))then
            ishypo=is
          endif
        enddo
        nhypo=0
        do is=1,ns
          if(togs(is).le.togs(ishypo))then
            nhypo=nhypo+1
          endif
        enddo
        if(nhypo.gt.1)then
          print *,' Warning: hypocenter cannot be uniquely determined!'
        endif
      endif
c
      do ir=1,nr
        call skip_comments(unit)
        read(unit,*)latr(ir),lonr(ir),rname(ir),tred(ir)
      enddo
c
      do flen=800,1,-1
        if(outfile(flen:flen).ne.' ')goto 300
      enddo
300   continue
      if(ioutform.eq.1)then
        uxout=outfile(1:flen)//'.un'
        uyout=outfile(1:flen)//'.ue'
        uzout=outfile(1:flen)//'.uz'
c
        vxout=outfile(1:flen)//'.vn'
        vyout=outfile(1:flen)//'.ve'
        vzout=outfile(1:flen)//'.vz'
c
        axout=outfile(1:flen)//'.an'
        ayout=outfile(1:flen)//'.ae'
        azout=outfile(1:flen)//'.az'
c
        tvout=outfile(1:flen)//'.tn'
        twout=outfile(1:flen)//'.te'
      else
        uxout=outfile(1:flen)//'.ut'
        uyout=outfile(1:flen)//'.up'
        uzout=outfile(1:flen)//'.ur'
c
        vxout=outfile(1:flen)//'.vt'
        vyout=outfile(1:flen)//'.vp'
        vzout=outfile(1:flen)//'.vr'
c
        axout=outfile(1:flen)//'.at'
        ayout=outfile(1:flen)//'.ap'
        azout=outfile(1:flen)//'.ar'
c
        tvout=outfile(1:flen)//'.tt'
        twout=outfile(1:flen)//'.tp'
      endif
      sdout=outfile(1:flen)//'.sd'
      grout=outfile(1:flen)//'.gr'
c
c     multilayered model parameters
c     =============================
c
      call skip_comments(unit)
      read(unit,*)l,i
      if(l.ge.lymax-2)then
        stop ' Error: lymax defined too small!'
      endif
      if(i.eq.1)then
        dispersion=.true.
      else
        dispersion=.false.
      endif
c
      qsmin=10000.d0
      depatmos=0.d0
      rratmos=REARTH
      do i=1,l
        call skip_comments(unit)
        read(unit,*)j,dp0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        if(i.eq.1)then
          depatmos=-KM2M*dp0(1)
          rratmos=REARTH+depatmos
        endif
        if(i.eq.l.and.KM2M*dp0(i).ne.REARTH)then
          if(KM2M*dp0(i).gt.1.001d0*REARTH)then
            stop ' Error: earth radius larger than pre-defined!'
          else if(KM2M*dp0(i).lt.0.999d0*REARTH)then
            stop ' Error: earth radius smaller than pre-defined!'
          else
            print *,' Warning: earth radius changed to pre-defined!'
            dp0(i)=REARTH/KM2M
          endif
        endif
c
        dp0(i)=KM2M*dp0(i)+depatmos
        vp0(i)=KM2M*vp0(i)
        vs0(i)=KM2M*vs0(i)
        ro0(i)=KM2M*ro0(i)
        if(vs0(i).gt.0.d0)qsmin=dmin1(qsmin,qs0(i))
        if(i.gt.1)then
          if(dp0(i).lt.dp0(i-1))then
            stop ' Error: bad layering of earth model!'
          endif
        endif
      enddo
c
      dpr=dpr+depatmos
      do i=1,ngrn
        grndep(i)=grndep(i)+depatmos
      enddo
      do is=1,ns
        deps(is)=deps(is)+depatmos
      enddo
c
      if(dp0(1).ne.0.d0)then
        stop ' Error: bad start depth!'
      else if(dp0(l).gt.rratmos)then
        stop ' Error: bad definition of earth radius!'
      else if(dp0(l).lt.rratmos)then
        l=l+1
        if(l.ge.lymax-2)stop ' Error: lymax defined too small!'
        dp0(l)=rratmos
        vp0(l)=vp0(l-1)
        vs0(l)=vs0(l-1)
        ro0(l)=ro0(l-1)
        qp0(l)=qp0(l-1)
        qs0(l)=qs0(l-1)
      endif
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then
          l0=l0+1
          dp0up(l0)=dp0(i-1)
          vp0up(l0)=vp0(i-1)
          vs0up(l0)=vs0(i-1)
          ro0up(l0)=ro0(i-1)
          qp0up(l0)=qp0(i-1)
          qs0up(l0)=qs0(i-1)
c
          dp0lw(l0)=dp0(i)
          vp0lw(l0)=vp0(i)
          vs0lw(l0)=vs0(i)
          ro0lw(l0)=ro0(i)
          qp0lw(l0)=qp0(i)
          qs0lw(l0)=qs0(i)
        endif
      enddo
c
c     end of inputs
c     =============
c
      return
      end
