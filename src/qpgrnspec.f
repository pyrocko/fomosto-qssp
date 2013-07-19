      subroutine qpgrnspec(ig)
      implicit none
      integer ig
c
      include 'qpglobal.h'
c
      integer i,il,istp,ly,lf,ldeg,ldeg0,ldegf,ldegup,ldegpre
      double precision f,ksp2,dll1,omi,expo
      double precision fl,depst1,depst2,dys2,dxs
      double precision kcut1(4),kcut2(4)
      double precision disk(0:ldegmax),fdisk(0:ldegmax)
      double complex ca,cb,cag,cs1,cs2,cs3,cs4,ct1,ct2,cll1
      double complex ys1(4,0:2),ys2(4,0:2),ys3(4,0:2)
      double complex ys4(4,0:2),ys5(4,0:2)
      double complex yt1(4,0:2),yt2(4,0:2)
      double complex ypsv(6,4),ypsvg(6,4),ysh(2,2)
      double complex ysgr(4),yspr(4),ystt(4),yttt(4)
c
      double precision expos,expoa,expof
      double complex c2,c3,c4
      data expos,expoa,expof/12.d0,12.d0,6.d0/
      data c2,c3,c4/(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
c     Initiation
c
      if(rr0.gt.0.d0)then
c
c       tappered disk source:
c       2*[cos(theta)-cos(alpha)]/[1-cos(alpha)]^2
c
        dxs=dcos(rr0/rrup(lys))
        dys2=(dsin(0.5d0*rr0/rrup(lys)))**2
        call qpdlegendre(ldegmax,rr0/rrup(lys))
        fdisk(0)=1.d0
        fdisk(1)=(2.d0+dxs)/3.d0
        do ldeg=2,ldegmax
          fdisk(ldeg)=0.5d0/dys2**2
     &               /dble(ldeg)/dble(ldeg+1)/dble(ldeg+2)
     &               *(2.d0*dys2*(2.d0*dble(ldeg)*dys2+1.d0)
     &                +dble(ldeg)*dxs**2*plm(ldeg,0)
     &                -dble(2*ldeg+1)*dxs*plm(ldeg-1,0)
     &                +dble(ldeg+1)*plm(ldeg-2,0))
        enddo
        open(19,file='disk.dat',status='unknown')
        write(19,'(a)')'     n        fdisk'
        do ldeg=0,ldegmax
          write(19,'(i6,E16.8)')ldeg,fdisk(ldeg)
        enddo
        close(19)
        do ldeg=0,ldegmax
          disk(ldeg)=fdisk(ldeg)*dble(2*ldeg+1)/(4.d0*PI*rrup(lys)**2)
        enddo
      else
        do ldeg=0,ldegmax
          fdisk(ldeg)=1.d0
          disk(ldeg)=dble(2*ldeg+1)/(4.d0*PI*rrup(lys)**2)
        enddo
      endif
c
      do ldeg=0,ldegmax
        fdisk(ldeg)=dabs(fdisk(ldeg))
        do il=ldeg+1,ldegmax
          fdisk(ldeg)=dmax1(fdisk(ldeg),dabs(fdisk(il)))
        enddo
      enddo
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      open(21,file=rgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(22,file=tgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(23,file=pgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(24,file=sgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(25,file=ggrnfile(ig),
     &     form='unformatted',status='unknown')
      open(26,file=vgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(27,file=wgrnfile(ig),
     &     form='unformatted',status='unknown')
c
      ldeg0=1+ndmax
      do ldeg0=1+ndmax,ldegmin
        dll1=dble(ldeg0)*dble(ldeg0+1)
        expo=0.d0
        do ly=min0(lys,lyr),max0(lys,lyr)-1
          expo=expo+dsqrt(dll1)*dlog(rrup(ly)/rrlw(ly))
        enddo
        if(expo.gt.expos+dlog(dsqrt(dll1)*fdisk(ldeg0)))goto 10
      enddo
10    continue
c
      omi=PI2*fcut
      ldegup=ldeg0
      do ldegup=ldeg0,min0(ldegmax,ldeg0+idint(REARTH*omi*slwmax))
        dll1=dble(ldegup)*dble(ldegup+1)
        expo=0.d0
        do ly=min0(lys,lyr),max0(lys,lyr)-1
          if(vsup(ly).gt.0.d0)then
            ksp2=(omi*rrup(ly)/vsup(ly))**2
          else
            ksp2=(omi*rrup(ly)/vpup(ly))**2
          endif
          if(dll1.gt.ksp2)then
            expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
          endif
        enddo
        if(expo.gt.expof+dlog(dsqrt(dll1)*fdisk(ldegup)))goto 20
      enddo
20    continue
c
      if(ldegup.ge.ldegmax-ndmax-1)then
        print *,' Warning in qpgrnspec: cutoff degree = ',ldegup,
     &          ' exceeds the limit!'
        ldegup=ldegmax-ndmax-1
        ldeg0=min0(ldeg0,ldegup)
      endif
c
      write(*,*)' '
      write(*,'(a,i3,a,f7.2,a)')' ... calculate Green functions for ',
     &        ig,'. source at depth ',(grndep(ig)-depatmos)/KM2M,' km'
	write(*,'(a,i5)')'   max. harmonic degree: L_max = ',ldegup
c
      write(21)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(22)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(23)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(24)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(25)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(26)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(27)nt,ntcut,dt,nf,nfcut,df,ldegup
c
      do istp=1,4
        do ldeg=0,1
          yr(ldeg,istp,0)=(0.d0,0.d0)
          yt(ldeg,istp,0)=(0.d0,0.d0)
          yp(ldeg,istp,0)=(0.d0,0.d0)
          ys(ldeg,istp,0)=(0.d0,0.d0)
          yg(ldeg,istp,0)=(0.d0,0.d0)
          yv(ldeg,istp,0)=(0.d0,0.d0)
          yw(ldeg,istp,0)=(0.d0,0.d0)
        enddo
      enddo
c
      ldegpre=ldeg0
      do lf=1,nfcut
        f=dble(lf-1)*df
        omi=PI2*f
        call qpqmodel(f)
c
        do ldegf=ldegpre,min0(ldegup,ldeg0+idint(REARTH*omi*slwmax))
          dll1=dble(ldegf)*dble(ldegf+1)
          expo=0.d0
          do ly=min0(lys,lyr),max0(lys,lyr)-1
            if(vsup(ly).gt.0.d0)then
              ksp2=(omi*rrup(ly)/vsup(ly))**2
            else
              ksp2=(omi*rrup(ly)/vpup(ly))**2
            endif
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
          enddo
          if(expo.gt.expos+dlog(dsqrt(dll1)*fdisk(ldegf)))goto 50
        enddo
50      ldegf=min0(ldegf,ldegup)
        ldegpre=ldegf
c
        do ldeg=0,ldegf+1
          dll1=dble(ldeg)*dble(ldeg+1)
c
c         determine degree dependent starting layer number
c         of sh solution
c
          lylwsh(ldeg)=min0(lycm,ly0)
          expo=0.d0
          do ly=max0(lys,lyr,lyob),min0(lycm,ly0)-1
            ksp2=(omi*rrup(ly)/vsup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwsh(ldeg)=ly
              goto 100
            endif
          enddo
100       continue
c
c         determine degree dependent starting layer number
c         of psv solution
c
          lylwpsv(ldeg)=ly0
          expo=0.d0
          do ly=max0(lys,lyr,lyob),min0(lycm,ly0)-1
            ksp2=(omi*rrup(ly)/vsup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwpsv(ldeg)=ly
              goto 200
            endif
          enddo
          do ly=lycm,min0(lycc,ly0)-1
            ksp2=(omi*rrup(ly)/vpup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwpsv(ldeg)=ly
              goto 200
            endif
          enddo
          do ly=lycc,ly0-1
            ksp2=(omi*rrup(ly)/vsup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwpsv(ldeg)=ly
              goto 200
            endif
          enddo
200       continue
c
          lyupatm(ldeg)=1
          expo=0.d0
          do ly=lyob-1,1,-1
            ksp2=(omi*rrup(ly)/vpup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expoa)then
              lyupatm(ldeg)=ly
              goto 300
            endif
          enddo
300       continue
        enddo
c
        if(selpsv.or.lys.lt.lyob)then
          depst1=(rratmos-depatmos-rrup(lylwpsv(0)))/KM2M
          depst2=(rratmos-depatmos-rrup(lylwpsv(ldegf+1)))/KM2M
        else
          depst1=(rratmos-depatmos-rrup(lylwsh(1)))/KM2M
          depst2=(rratmos-depatmos-rrup(lylwsh(ldegf+1)))/KM2M
        endif
c
c       determine layer dependent max. harmonic degree
c       of sh solution
c
        do ly=lyob,min0(lycm-1,ly0)
          ldegsh(ly)=1
          do ldeg=1,ldegf+1
            if(lylwsh(ldeg).ge.ly)then
              ldegsh(ly)=ldeg
            endif
          enddo
        enddo
c
c       determine layer dependent max. harmonic degree
c       of psv solution
c
        do ly=1,lyob-1
          ldegpsv(ly)=0
          do ldeg=0,ldegf+1
            if(lyupatm(ldeg).le.ly)then
              ldegpsv(ly)=ldeg
            endif
          enddo
        enddo
c
        do ly=lyob,ly0
          ldegpsv(ly)=0
          do ldeg=0,ldegf+1
            if(lylwpsv(ldeg).ge.ly)then
              ldegpsv(ly)=ldeg
            endif
          enddo
        enddo
c
        do ly=1,ly0
          ksmallpsv(ly)=.true.
          ksmallsh(ly)=.true.
        enddo
c
        do il=0,2
          do istp=1,4
            ys1(istp,il)=(0.d0,0.d0)
            ys2(istp,il)=(0.d0,0.d0)
            ys3(istp,il)=(0.d0,0.d0)
            ys4(istp,il)=(0.d0,0.d0)
            ys5(istp,il)=(0.d0,0.d0)
            yt1(istp,il)=(0.d0,0.d0)
            yt2(istp,il)=(0.d0,0.d0)
          enddo
        enddo
c
        do ldeg=0,ldegf+1
          do istp=1,4
            ys1(istp,0)=ys1(istp,1)
            ys2(istp,0)=ys2(istp,1)
            ys3(istp,0)=ys3(istp,1)
            ys4(istp,0)=ys4(istp,1)
            ys5(istp,0)=ys5(istp,1)
c
            ys1(istp,1)=ys1(istp,2)
            ys2(istp,1)=ys2(istp,2)
            ys3(istp,1)=ys3(istp,2)
            ys4(istp,1)=ys4(istp,2)
            ys5(istp,1)=ys5(istp,2)
c
            yt1(istp,0)=yt1(istp,1)
            yt2(istp,0)=yt2(istp,1)
c
            yt1(istp,1)=yt1(istp,2)
            yt2(istp,1)=yt2(istp,2)
          enddo
c
          if(.not.selpsv.or.lys.lt.lyob.and.f.le.0.d0)then
            do istp=1,4
	        do i=1,6
                ypsv(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else if(nogravity)then
            call qppsvkern(f,ldeg,ypsv)
          else
            fl=fgr*dble(ldeg)/dble(ldeggr)
            if(f.le.fgr-fl)then
              call qppsvkerng(f,ldeg,ypsv)
            else if(f.ge.fgr*(1.d0+FLTAPER)-fl)then
              call qppsvkern(f,ldeg,ypsv)
            else
              call qppsvkerng(f,ldeg,ypsvg)
              call qppsvkern(f,ldeg,ypsv)
              ca=dcmplx(dsin(((f+fl)/fgr-1.d0)/FLTAPER)**2,0.d0)
              cb=(1.d0,0.d0)-ca
              do istp=1,4
                do i=1,6
                  ypsv(i,istp)=ca*ypsv(i,istp)+cb*ypsvg(i,istp)
                enddo
              enddo
            endif
          endif
c
          if(.not.selsh.or.lys.lt.lyob)then
            do istp=1,2
	        do i=1,2
                ysh(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else
            call qpshkern(f,ldeg,ysh)
          endif
c          
          if(lyr.le.lyob)then
            do istp=1,4
              yspr(istp)=-ypsv(2,istp)
            enddo
          else
            ca=(claup(lyr)+cmuup(lyr)*c2/c3)/(claup(lyr)+cmuup(lyr)*c2)
            cb=cmuup(lyr)/dcmplx(rrup(lyr),0.d0)
            cll1=dcmplx(dble(ldeg)*dble(ldeg+1),0.d0)
            do istp=1,4
              yspr(istp)=-ca*(ypsv(2,istp)+c4*cb*ypsv(1,istp)
     &                   -c2*cll1*cb*ypsv(3,istp))
            enddo
          endif
c
          comi2=dcmplx(PI2*f,PI2*fi)**2
          ca=dcmplx(dble(ldeg+1)/rrup(lyr),0.d0)
          cb=dcmplx(freeairgrd,0.d0)-comi2
          do istp=1,4
            ysgr(istp)=-ypsv(6,istp)+ca*ypsv(5,istp)+cb*ypsv(1,istp)
          enddo
c
          do istp=1,4
            ystt(istp)=(ypsv(5,istp)/cgrup(lyr)-ypsv(1,istp))/crrup(lyr)
     &                +comi2*ypsv(3,istp)/cgrup(lyr)
          enddo
c
          do istp=1,2
            yttt(istp)=comi2*ysh(1,istp)/cgrup(lyr)
          enddo
c
c         1. Explosion (M11=M22=M33=1)
c
          cs1=dcmplx( disk(ldeg)/(roup(lys)*vpup(lys)**2),0.d0)
          cs2=dcmplx(-disk(ldeg)*4.d0*(vsup(lys)/vpup(lys))**2
     &                /rrup(lys),0.d0)
          cs4=dcmplx( disk(ldeg)*2.d0*(vsup(lys)/vpup(lys))**2
     &                /rrup(lys),0.d0)
          ys1(1,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
          ys2(1,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
          ys3(1,2)=cs1*yspr(1)+cs2*yspr(2)+cs4*yspr(4)
          ys4(1,2)=cs1*ysgr(1)+cs2*ysgr(2)+cs4*ysgr(4)
          ys5(1,2)=cs1*ystt(1)+cs2*ystt(2)+cs4*ystt(4)
          yt1(1,2)=(0.d0,0.d0)
          yt2(1,2)=(0.d0,0.d0)
c
c         2. Strike-slip (M12=M21=1)
c
          if(ldeg.lt.2.or.lys.lt.lyob)then
            ys1(2,2)=(0.d0,0.d0)
            ys2(2,2)=(0.d0,0.d0)
            ys3(2,2)=(0.d0,0.d0)
            ys4(2,2)=(0.d0,0.d0)
            ys5(2,2)=(0.d0,0.d0)
            yt1(2,2)=(0.d0,0.d0)
            yt2(2,2)=(0.d0,0.d0)
          else
            ct2=dcmplx(disk(ldeg)
     &                /(dble(ldeg)*dble(ldeg+1)*rrup(lys)),0.d0)
            cs4=-ct2
            ys1(2,2)=cs4*ypsv(1,4)
            ys2(2,2)=cs4*ypsv(3,4)
            ys3(2,2)=cs4*yspr(4)
            ys4(2,2)=cs4*ysgr(4)
            ys5(2,2)=cs4*ystt(4)
            yt1(2,2)=ct2*ysh(1,2)
            yt2(2,2)=ct2*yttt(2)
          endif
c
c         3. Dip-slip (M13=M31=1)
c
          if(ldeg.lt.1.or.lys.lt.lyob)then
            ys1(3,2)=(0.d0,0.d0)
            ys2(3,2)=(0.d0,0.d0)
            ys3(3,2)=(0.d0,0.d0)
            ys4(3,2)=(0.d0,0.d0)
            ys5(3,2)=(0.d0,0.d0)
            yt1(3,2)=(0.d0,0.d0)
            yt2(3,2)=(0.d0,0.d0)
          else
            ct1=dcmplx(disk(ldeg)/(dble(ldeg)*dble(ldeg+1)
     &                *roup(lys)*vsup(lys)**2),0.d0)
            cs3=ct1
            ys1(3,2)=cs3*ypsv(1,3)
            ys2(3,2)=cs3*ypsv(3,3)
            ys3(3,2)=cs3*yspr(3)
            ys4(3,2)=cs3*ysgr(3)
            ys5(3,2)=cs3*ystt(3)
            yt1(3,2)=ct1*ysh(1,1)
            yt2(3,2)=ct1*yttt(1)
          endif
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c
          if(lys.lt.lyob)then
            ys1(4,2)=(0.d0,0.d0)
            ys2(4,2)=(0.d0,0.d0)
            ys3(4,2)=(0.d0,0.d0)
            ys4(4,2)=(0.d0,0.d0)
            ys5(4,2)=(0.d0,0.d0)
          else
            cs1=dcmplx( disk(ldeg)/(roup(lys)*vpup(lys)**2),0.d0)
            cs2=dcmplx( disk(ldeg)*(3.d0-4.d0*(vsup(lys)/vpup(lys))**2)
     &                  /rrup(lys),0.d0)
            cs4=-(0.5d0,0.d0)*cs2
            ys1(4,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
            ys2(4,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
            ys3(4,2)=cs1*yspr(1)+cs2*yspr(2)+cs4*yspr(4)
            ys4(4,2)=cs1*ysgr(1)+cs2*ysgr(2)+cs4*ysgr(4)
            ys5(4,2)=cs1*ystt(1)+cs2*ystt(2)+cs4*ystt(4)
          endif
          yt1(4,2)=(0.d0,0.d0)
          yt2(4,2)=(0.d0,0.d0)
c
c         convert to (r,t,p,g)-system
c         ===========================
c
c         1. Explosion (M11=M22=M33=1)
c            yr,yg normalized by Plm(l,0,cos(t))
c            yt normalized by Plm(l,1,cos(t)) = -dPlm/dt
c
          yr(ldeg,1,0)=ys1(1,2)
          yt(ldeg,1,0)=-ys2(1,2)
          ys(ldeg,1,0)=ys3(1,2)
          yg(ldeg,1,0)=ys4(1,2)
          yv(ldeg,1,0)=-ys5(1,2)
c
c         2. Strike-slip (M12=M21=1)
c            yr,yg normalized by Plm(l,2,cos(t))*sin(2p)
c            yt normalized by Plm(l,2,cos(t))*sin(2p)/sin(t)
c            yp normalized by Plm(l,2,cos(t))*cos(2p)/sin(t)
c
          yr(ldeg,2,0)=ys1(2,2)
          if(ldeg.ge.3)then
            ca=dcmplx(dble(ldeg-2)*dble(ldeg-3)/dble(2*ldeg-3),0.d0)
            cb=dcmplx(dble(ldeg+1)*dble(ldeg+2)/dble(2*ldeg+1),0.d0)
            yt(ldeg-1,2,0)= ca*ys2(2,0)-cb*ys2(2,2)-(2.d0,0.d0)*yt1(2,1)
            yp(ldeg-1,2,0)=-ca*yt1(2,0)+cb*yt1(2,2)+(2.d0,0.d0)*ys2(2,1)
            yv(ldeg-1,2,0)= ca*ys5(2,0)-cb*ys5(2,2)-(2.d0,0.d0)*yt2(2,1)
            yw(ldeg-1,2,0)=-ca*yt2(2,0)+cb*yt2(2,2)+(2.d0,0.d0)*ys5(2,1)
          endif
          ys(ldeg,2,0)=ys3(2,2)
          yg(ldeg,2,0)=ys4(2,2)
c
c         3. Dip-slip (M13=M31=1)
c            yr,yg normalized by Plm(l,2,cos(t))*cos(p)
c            yt normalized by Plm(l,1,cos(t))*cos(p)/sin(t)
c            yp normalized by Plm(l,1,cos(t))*sin(p)/sin(t)
c
          yr(ldeg,3,0)=ys1(3,2)
          if(ldeg.ge.2)then
            ca=dcmplx(dble(ldeg-2)*dble(ldeg-2)/dble(2*ldeg-3),0.d0)
            cb=dcmplx(dble(ldeg+1)*dble(ldeg+1)/dble(2*ldeg+1),0.d0)
            yt(ldeg-1,3,0)= ca*ys2(3,0)-cb*ys2(3,2)+yt1(3,1)
            yp(ldeg-1,3,0)=-ca*yt1(3,0)+cb*yt1(3,2)-ys2(3,1)
            yv(ldeg-1,3,0)= ca*ys5(3,0)-cb*ys5(3,2)+yt2(3,1)
            yw(ldeg-1,3,0)=-ca*yt2(3,0)+cb*yt2(3,2)-ys5(3,1)
          endif
          ys(ldeg,3,0)=ys3(3,2)
          yg(ldeg,3,0)=ys4(3,2)
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c            yr,yg normalized by Plm(l,0,cos(t))
c            yt normalized by Plm(l,0,cos(t)) = -dPlm/dt
c
          yr(ldeg,4,0)=ys1(4,2)
          yt(ldeg,4,0)=-ys2(4,2)
          ys(ldeg,4,0)=ys3(4,2)
          yg(ldeg,4,0)=ys4(4,2)
          yv(ldeg,4,0)=-ys5(4,2)
        enddo
c
        write(*,'(i6,a,f12.4,a,i5,a,2(f7.2,a))')lf,'.',1.0d+03*f,
     &       ' mHz: cut-off degree = ',ldegf,
     &       ', start depth = ',depst1,' - ',depst2,' km'
c
        write(21)ldegf
        write(22)ldegf
        write(23)ldegf
        write(24)ldegf
        write(25)ldegf
        write(26)ldegf
        write(27)ldegf
        write(21)((yr(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(22)((yt(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(23)((yp(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
        write(24)((ys(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(25)((yg(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(26)((yv(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(27)((yw(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
      enddo
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      return
      end
