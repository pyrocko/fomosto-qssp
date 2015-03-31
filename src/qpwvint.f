      subroutine qpwvint(ierr)
      implicit none
      integer ierr
c
      include 'qpglobal.h'
c
      integer i,id,is,ir,ig,nd,nt0,nf0,ntcut0,nfcut0,ly,ishift
      integer lf,lf1,istp,ldeg,ldeg0,ldegf,ldegf0
      integer istat,ldegup,ldeglw,ldegneed
      integer ldegtap(4,nsmax,nrmax)
      double precision depsarc,rrs,slws,slwr,twin,anorm
      double precision f,dt0,df0,rn,re,azi,bazi,bazi0,f20
      double precision tap(0:ldegmax),ldf(nsmax,nrmax)
      double precision f0(nsmax,nrmax)
      double complex cp0,cp1,cp2,wavelet
      double complex cfac,ca,cb,dur,dut,dup,dus,dug,duv,duw
      double complex expl(nsmax),clvd(nsmax),ss12(nsmax)
      double complex ss11(nsmax),ds31(nsmax),ds23(nsmax)
c
      do ldeg=0,ldegmax
        tap(ldeg)=1.d0
      enddo
c
      do is=1,ns
        depsarc=deps(is)/REARTH
        rrs=REARTH-deps(is)
        lys=lyob
        do ly=lyob+1,min0(lycm-1,ly0)
          if(rrs.le.rrup(ly))then
            lys=ly
          endif
        enddo
        slws=dmin1(slwupcut,slwmax)*REARTH
        do ir=1,nr
          call disazi(1.d0,lats(is),lons(is),
     &                     latr(ir),lonr(ir),rn,re)
c
          dis(is,ir)=dsqrt(rn**2+re**2)
          twin=dble(ntcutout)*dt+tred(ir)-togs(is)
          if(dis(is,ir).gt.0.d0)then
            slwr=twin/dis(is,ir)
          else
            slwr=slws
          endif
c
          if(slwr.ge.slws)then
            ldf(is,ir)=0.d0
            idr(is,ir)=0
          else
            ldf(is,ir)=PI2*slwr
            idr(is,ir)=min0(ndmax,
     &             idnint(dlog(dis(is,ir)/depsarc)/dlog(10.d0)))
          endif
c
c         f0 is the frequency, at which the distance is equal to 10 times
c         of cut-off wavelength
c
          f0(is,ir)=10.d0/(slwr*dsqrt(dis(is,ir)**2+depsarc**2))
c
          ssd(is,ir)=dcmplx(dsin(dis(is,ir)),0.d0)
          ssf(is,ir)=dcmplx(2.d0*dsin(0.5d0*dis(is,ir))**2,0.d0)
c
c         azi = station azimuth (from south to east)
c
          azi=datan2(re,-rn)
          ssa(is,ir)=dcmplx(dsin(azi),0.d0)
          csa(is,ir)=dcmplx(dcos(azi),0.d0)
          ss2a(is,ir)=dcmplx(dsin(2.d0*azi),0.d0)
          cs2a(is,ir)=dcmplx(dcos(2.d0*azi),0.d0)
        enddo
      enddo
c
      if(ioutform.eq.1)then
        do ir=1,nr
          do is=1,ns
            call disazi(1.d0,latr(ir),lonr(ir),
     &                       lats(is),lons(is),rn,re)
            bazi=datan2(-re,-rn)
            ssb(is,ir)=dcmplx(dsin(bazi),0.d0)
            csb(is,ir)=dcmplx(dcos(bazi),0.d0)
          enddo
        enddo
      else
        do ir=1,nr
          call disazi(1.d0,latr(ir),lonr(ir),
     &                     lats(ishypo),lons(ishypo),rn,re)
          bazi0=datan2(-re,-rn)
          do is=1,ns
            call disazi(1.d0,latr(ir),lonr(ir),
     &                       lats(is),lons(is),rn,re)
            bazi=datan2(-re,-rn)
            ssb(is,ir)=dcmplx(dsin(bazi-bazi0),0.d0)
            csb(is,ir)=dcmplx(dcos(bazi-bazi0),0.d0)
          enddo
        enddo
      endif
c
      do is=1,ns
        expl(is)=dcmplx((mtt(is)+mpp(is)+mrr(is))/3.d0,0.d0)
        clvd(is)=dcmplx(mrr(is),0.d0)-expl(is)
        ss12(is)=dcmplx(mtp(is),0.d0)
        ss11(is)=dcmplx((mtt(is)-mpp(is))/2.d0,0.d0)
        ds31(is)=dcmplx(mrt(is),0.d0)
        ds23(is)=dcmplx(mpr(is),0.d0)
      enddo
c
c     initiation
c
      do lf=1,nf
        do ir=1,nr
          uz(lf,ir)=(0.d0,0.d0)
          ux(lf,ir)=(0.d0,0.d0)
          uy(lf,ir)=(0.d0,0.d0)
          us(lf,ir)=(0.d0,0.d0)
          ug(lf,ir)=(0.d0,0.d0)
          uv(lf,ir)=(0.d0,0.d0)
          uw(lf,ir)=(0.d0,0.d0)
        enddo
      enddo
c
      do ig=1,ngrn
        if(nsg(ig).le.0)goto 500
        write(*,'(a)')' '
        write(*,'(a,i4,a,f5.1,a)')' processing ',1+isg2(ig)-isg1(ig),
     &    ' point source(s) at depth ',(grndep(ig)-depatmos)/KM2M,' km'
        write(*,'(a)')' open Green function data base: '
     &              //grnfile(ig)(1:40)
        write(*,'(a)')' ... please wait ...'
c
        open(21,file=rgrnfile(ig),
     &       form='unformatted',status='old')
        open(22,file=tgrnfile(ig),
     &       form='unformatted',status='old')
        open(23,file=pgrnfile(ig),
     &       form='unformatted',status='old')
        open(24,file=sgrnfile(ig),
     &       form='unformatted',status='old')
        open(25,file=ggrnfile(ig),
     &       form='unformatted',status='old')
        open(26,file=vgrnfile(ig),
     &       form='unformatted',status='old')
        open(27,file=wgrnfile(ig),
     &       form='unformatted',status='old')
c
        read(21)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(22)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(23)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(24)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(25)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(26)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(27)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
c
        if(ntcut0.ne.ntcut.or.dt0.ne.dt.or.
     &     nfcut0.ne.nfcut.or.df0.ne.df)then
          print *,' Error in qpwvint: t/f sampling'
     &          //' inconsistent with Green functions!'
          write(*,'(a)')'               '//'  ntcut             dt'
     &                                   //'  nfcut             df'
          write(*,'(a,2(i7,f16.8))')' Current input:',
     &                               ntcut,dt,nfcut,df
          write(*,'(a,2(i7,f16.8))')'     Data base:',
     &                               ntcut0,dt0,nfcut0,df0
          stop
        endif
c
        nd=0
        do is=isg1(ig),isg2(ig)
          do ir=1,nr
            nd=max0(nd,idr(is,ir))
          enddo
        enddo
c
        if(ldegup+nd.gt.ldegmax)then
          write(*,'(a,i6,a,i6)')' Error in qpwvint: '
     &     //'max. harmonic degree required = ',ldegup+nd,
     &     ' > ldegmax defined: ',ldegmax
          stop
        endif
c
        ldegf0=0
        do lf=1,nfcut
          f=dble(lf-1)*df
          read(21)ldegf
          read(22)ldegf
          read(23)ldegf
          read(24)ldegf
          read(25)ldegf
          read(26)ldegf
          read(27)ldegf
c
          read(21)((yr(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(22)((yt(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(23)((yp(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
          read(24)((ys(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(25)((yg(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(26)((yv(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(27)((yw(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
c
          if(lf.eq.1)ldegf0=ldegf
          ldegneed=0
          ldeglw=ldegf
          ldeg0=ldegf
          do is=isg1(ig),isg2(ig)
            do ir=1,nr
              id=idr(is,ir)
              if(ldf(is,ir).le.0.d0.or.dis(is,ir).le.0.d0)then
                ldegtap(4,is,ir)=max0(0,ldegf-id)
              else
                ldegtap(4,is,ir)=max0(0,min0(ldegf-id,
     &             ldegf0+idnint(dsqrt(f**2+f0(is,ir)**2)*ldf(is,ir))))
              endif
              ldegtap(3,is,ir)=max0(0,min0(ldegtap(4,is,ir)*4/5,
     &                            max0(0,ldegtap(4,is,ir)-40)))
              ldegneed=max0(ldegneed,ldegtap(4,is,ir))
              ldeglw=min0(ldeglw,ldegtap(4,is,ir))
              ldegtap(2,is,ir)=min0(ldegtap(3,is,ir),
     &                              idnint(REARTH*PI2*f*slwlwcut))
              ldegtap(1,is,ir)=min0(ldegtap(3,is,ir),
     &                              idnint(0.8d0*REARTH*PI2*f*slwlwcut))
              ldeg0=min0(ldeg0,ldegtap(1,is,ir))
            enddo
          enddo
          ldegneed=min0(ldegneed+nd,ldegf)          
c
c         use differential filter to suppress spatial aliasing
c
          do id=0,nd
            yt(0,1,id)=(0.d0,0.d0)
            yt(0,4,id)=(0.d0,0.d0)
            yv(0,1,id)=(0.d0,0.d0)
            yv(0,4,id)=(0.d0,0.d0)
c
            yr(0,3,id)=(0.d0,0.d0)
            yt(0,3,id)=(0.d0,0.d0)
            yp(0,3,id)=(0.d0,0.d0)
            ys(0,3,id)=(0.d0,0.d0)
            yg(0,3,id)=(0.d0,0.d0)
            yv(0,3,id)=(0.d0,0.d0)
            yw(0,3,id)=(0.d0,0.d0)
c
            yr(0,2,id)=(0.d0,0.d0)
            yt(0,2,id)=(0.d0,0.d0)
            yp(0,2,id)=(0.d0,0.d0)
            ys(0,2,id)=(0.d0,0.d0)
            yg(0,2,id)=(0.d0,0.d0)
            yv(0,2,id)=(0.d0,0.d0)
            yw(0,2,id)=(0.d0,0.d0)
c
            yr(1,2,id)=(0.d0,0.d0)
            yt(1,2,id)=(0.d0,0.d0)
            yp(1,2,id)=(0.d0,0.d0)
            ys(1,2,id)=(0.d0,0.d0)
            yg(1,2,id)=(0.d0,0.d0)
            yv(1,2,id)=(0.d0,0.d0)
            yw(1,2,id)=(0.d0,0.d0)
          enddo
c
          do id=1,nd
c
c           m = 0
c
            ca=dcmplx(1.d0/3.d0,0.d0)
            yr(0,1,id)=yr(0,1,id-1)-ca*yr(1,1,id-1)
            ys(0,1,id)=ys(0,1,id-1)-ca*ys(1,1,id-1)
            yg(0,1,id)=yg(0,1,id-1)-ca*yg(1,1,id-1)
            yr(0,4,id)=yr(0,4,id-1)-ca*yr(1,4,id-1)
            ys(0,4,id)=ys(0,4,id-1)-ca*ys(1,4,id-1)
            yg(0,4,id)=yg(0,4,id-1)-ca*yg(1,4,id-1)
            do ldeg=1,ldegneed-id
              ca=dcmplx(dble(ldeg+1)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg)/dble(2*ldeg-1),0.d0)
              yr(ldeg,1,id)=yr(ldeg,1,id-1)-ca*yr(ldeg+1,1,id-1)
     &                     -cb*yr(ldeg-1,1,id-1)
              ys(ldeg,1,id)=ys(ldeg,1,id-1)-ca*ys(ldeg+1,1,id-1)
     &                     -cb*ys(ldeg-1,1,id-1)
              yg(ldeg,1,id)=yg(ldeg,1,id-1)-ca*yg(ldeg+1,1,id-1)
     &                     -cb*yg(ldeg-1,1,id-1)
              yr(ldeg,4,id)=yr(ldeg,4,id-1)-ca*yr(ldeg+1,4,id-1)
     &                     -cb*yr(ldeg-1,4,id-1)
              ys(ldeg,4,id)=ys(ldeg,4,id-1)-ca*ys(ldeg+1,4,id-1)
     &                     -cb*ys(ldeg-1,4,id-1)
              yg(ldeg,4,id)=yg(ldeg,4,id-1)-ca*yg(ldeg+1,4,id-1)
     &                     -cb*yg(ldeg-1,4,id-1)
            enddo
c
c           m = 1
c
            do ldeg=1,ldegneed-id
              ca=dcmplx(dble(ldeg+2)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-1)/dble(2*ldeg-1),0.d0)
              yt(ldeg,1,id)=yt(ldeg,1,id-1)-ca*yt(ldeg+1,1,id-1)
     &                     -cb*yt(ldeg-1,1,id-1)
              yt(ldeg,4,id)=yt(ldeg,4,id-1)-ca*yt(ldeg+1,4,id-1)
     &                     -cb*yt(ldeg-1,4,id-1)
              yv(ldeg,1,id)=yv(ldeg,1,id-1)-ca*yv(ldeg+1,1,id-1)
     &                     -cb*yv(ldeg-1,1,id-1)
              yv(ldeg,4,id)=yv(ldeg,4,id-1)-ca*yv(ldeg+1,4,id-1)
     &                     -cb*yv(ldeg-1,4,id-1)
              yr(ldeg,3,id)=yr(ldeg,3,id-1)-ca*yr(ldeg+1,3,id-1)
     &                     -cb*yr(ldeg-1,3,id-1)
              yt(ldeg,3,id)=yt(ldeg,3,id-1)-ca*yt(ldeg+1,3,id-1)
     &                     -cb*yt(ldeg-1,3,id-1)
              yp(ldeg,3,id)=yp(ldeg,3,id-1)-ca*yp(ldeg+1,3,id-1)
     &                     -cb*yp(ldeg-1,3,id-1)
              ys(ldeg,3,id)=ys(ldeg,3,id-1)-ca*ys(ldeg+1,3,id-1)
     &                     -cb*ys(ldeg-1,3,id-1)
              yg(ldeg,3,id)=yg(ldeg,3,id-1)-ca*yg(ldeg+1,3,id-1)
     &                     -cb*yg(ldeg-1,3,id-1)
              yv(ldeg,3,id)=yv(ldeg,3,id-1)-ca*yv(ldeg+1,3,id-1)
     &                     -cb*yv(ldeg-1,3,id-1)
              yw(ldeg,3,id)=yw(ldeg,3,id-1)-ca*yw(ldeg+1,3,id-1)
     &                     -cb*yw(ldeg-1,3,id-1)
            enddo
c
c           m = 2
c
            do ldeg=2,ldegneed-id
              ca=dcmplx(dble(ldeg+3)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-2)/dble(2*ldeg-1),0.d0)
              yr(ldeg,2,id)=yr(ldeg,2,id-1)-ca*yr(ldeg+1,2,id-1)
     &                     -cb*yr(ldeg-1,2,id-1)
              yt(ldeg,2,id)=yt(ldeg,2,id-1)-ca*yt(ldeg+1,2,id-1)
     &                     -cb*yt(ldeg-1,2,id-1)
              yp(ldeg,2,id)=yp(ldeg,2,id-1)-ca*yp(ldeg+1,2,id-1)
     &                     -cb*yp(ldeg-1,2,id-1)
              ys(ldeg,2,id)=ys(ldeg,2,id-1)-ca*ys(ldeg+1,2,id-1)
     &                     -cb*ys(ldeg-1,2,id-1)
              yg(ldeg,2,id)=yg(ldeg,2,id-1)-ca*yg(ldeg+1,2,id-1)
     &                     -cb*yg(ldeg-1,2,id-1)
              yv(ldeg,2,id)=yv(ldeg,2,id-1)-ca*yv(ldeg+1,2,id-1)
     &                     -cb*yv(ldeg-1,2,id-1)
              yw(ldeg,2,id)=yw(ldeg,2,id-1)-ca*yw(ldeg+1,2,id-1)
     &                     -cb*yw(ldeg-1,2,id-1)
            enddo
          enddo
c
          do is=isg1(ig),isg2(ig)

            wavelet=(1.d0,0.d0)/dcmplx(1.d0,pi2*f*trss(is))**2
     &             *cdexp(-dcmplx(-fi,f)*dcmplx(PI2*togs(is),0.d0))
            do ir=1,nr
              id=idr(is,ir)
c
              call taper(ldegtap(1,is,ir),ldegtap(4,is,ir),tap(0))
c
              call qplegendre(ldegtap(4,is,ir),dis(is,ir))
c
              do ldeg=ldegtap(1,is,ir),ldegtap(4,is,ir)
                cp0=dcmplx(plm(ldeg,0),0.d0)
                cp1=dcmplx(plm(ldeg,1),0.d0)
                cp2=dcmplx(plm(ldeg,2),0.d0)
                dur=(expl(is)*yr(ldeg,1,id)+clvd(is)*yr(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *yr(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *yr(ldeg,3,id)*cp1*ssd(is,ir)
                dus=(expl(is)*ys(ldeg,1,id)+clvd(is)*ys(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *ys(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *ys(ldeg,3,id)*cp1*ssd(is,ir)
                dug=(expl(is)*yg(ldeg,1,id)+clvd(is)*yg(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *yg(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *yg(ldeg,3,id)*cp1*ssd(is,ir)
                dut=(expl(is)*yt(ldeg,1,id)+clvd(is)*yt(ldeg,4,id))
     &             *cp1*ssd(is,ir)
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *yt(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *yt(ldeg,3,id)*cp1
                dup=(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &             *yp(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &             *yp(ldeg,3,id)*cp1
                duv=(expl(is)*yv(ldeg,1,id)+clvd(is)*yv(ldeg,4,id))
     &             *cp1*ssd(is,ir)
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *yv(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *yv(ldeg,3,id)*cp1
                duw=(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &             *yw(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &             *yw(ldeg,3,id)*cp1
                cfac=dcmplx(tap(ldeg),0.d0)*wavelet
                if(id.gt.0)cfac=cfac/ssf(is,ir)**id
                uz(lf,ir)=uz(lf,ir)+dur*cfac
                us(lf,ir)=us(lf,ir)+dus*cfac
                ug(lf,ir)=ug(lf,ir)+dug*cfac
                if(ioutform.eq.1)then
                  ux(lf,ir)=ux(lf,ir)
     &                     +(dut*csb(is,ir)+dup*ssb(is,ir))*cfac
                  uy(lf,ir)=uy(lf,ir)
     &                     +(dut*ssb(is,ir)-dup*csb(is,ir))*cfac
                  uv(lf,ir)=uv(lf,ir)
     &                     +(duv*csb(is,ir)+duw*ssb(is,ir))*cfac
                  uw(lf,ir)=uw(lf,ir)
     &                     +(duv*ssb(is,ir)-duw*csb(is,ir))*cfac
                else
                  ux(lf,ir)=ux(lf,ir)
     &                     +(dut*csb(is,ir)+dup*ssb(is,ir))*cfac
                  uy(lf,ir)=uy(lf,ir)
     &                     -(dut*ssb(is,ir)-dup*csb(is,ir))*cfac
                  uv(lf,ir)=uv(lf,ir)
     &                     +(duv*csb(is,ir)+duw*ssb(is,ir))*cfac
                  uw(lf,ir)=uw(lf,ir)
     &                     -(duv*ssb(is,ir)-duw*csb(is,ir))*cfac
                endif
              enddo
            enddo
          enddo
          write(*,'(i6,a,f10.4,3(a,i5))')lf,'.',1.0d+03*f,
     &                        ' mHz: spectra read: ',ldegf,
     &                        ', used: ',ldeglw,' - ',ldegneed
        enddo
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(27)
        write(*,'(i6,a)')lf-1,' spectra read from '
     &                      //grnfile(ig)(1:40)
500     continue
      enddo
      return
      end
