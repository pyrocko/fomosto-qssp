      subroutine qpfftinv(ierr)
      implicit none
      integer ierr
c
      include 'qpglobal.h'
c
      integer lf,mf,ir,it
      double precision t,sr,si
      double precision f(nfmax),y0(nrmax),y(nrmax)
      double complex s
      double complex swap(2*nfmax),lpf(nfmax)
c
c      open(39,file='static.dat',status='unknown')
c      write(39,'(a,$)')'        Lat         Lon'
c      if(ioutform.eq.1)then
c        write(39,'(a)')'         Un          Ue          Uz          Gr'
c      else
c        write(39,'(a)')'         Ut          Up          Ur          Gr'
c      endif
c      do ir=1,nr
c        write(39,'(2f12.6,4E16.8)')latr(ir),lonr(ir),dreal(ux(1,ir)),
c     &      dreal(uy(1,ir)),dreal(uz(1,ir)),dreal(ug(1,ir))
c      enddo
c      close(39)
c
      do lf=1,nf
        f(lf)=dble(lf-1)*df
      enddo
c
c     low-pass filter
c
      if(fcorner.gt.0.d0.and.nlpf.gt.0)then
        call butterworth(nlpf,fcorner,df,nf,lpf)
        do lf=1,nf
          do ir=1,nr
            ux(lf,ir)=ux(lf,ir)*lpf(lf)
            uy(lf,ir)=uy(lf,ir)*lpf(lf)
            uz(lf,ir)=uz(lf,ir)*lpf(lf)
            us(lf,ir)=us(lf,ir)*lpf(lf)
            ug(lf,ir)=ug(lf,ir)*lpf(lf)
            uv(lf,ir)=uv(lf,ir)*lpf(lf)
            uw(lf,ir)=uw(lf,ir)*lpf(lf)
          enddo
        enddo
      endif
c
c     calculate north or radial component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=ux(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf  ))*dexp(-PI2*fi*t)*df
          ux(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=uxout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(ux(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(ux(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
c
      open(21,file=vxout,status='unknown')
      write(21,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(21,'(a16,$)')rname(ir)
      enddo
      write(21,'(a16)')rname(nr)
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(21,'(f12.3,$)')t
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=dreal(ux(lf,ir))
          enddo
        else
          do ir=1,nr
            y(ir)=dimag(ux(lf,ir))
          enddo
        endif
        do ir=1,nr-1
          write(21,'(E16.8,$)')y(ir)
        enddo
        write(21,'(E16.8)')y(nr)
      enddo
      close(21)
c
      open(22,file=axout,status='unknown')
      write(22,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(22,'(a16,$)')rname(ir)
      enddo
      write(22,'(a16)')rname(nr)
      do ir=1,nr
        y0(ir)=dreal(ux(1,ir))
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=dreal(ux(lf,ir))
          enddo
        else
          do ir=1,nr
            y(ir)=dimag(ux(lf,ir))
          enddo
        endif
        write(22,'(f12.3,$)')t
        do ir=1,nr-1
          write(22,'(E16.8,$)')(y(ir)-y0(ir))/dt
        enddo
        write(22,'(E16.8)')(y(nr)-y0(nr))/dt
        do ir=1,nr
          y0(ir)=y(ir)
        enddo
      enddo
      close(22)
c
c     calculate east or tangential component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=uy(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf)  )*dexp(-PI2*fi*t)*df
          uy(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=uyout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(uy(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(uy(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
c
      open(21,file=vyout,status='unknown')
      write(21,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(21,'(a16,$)')rname(ir)
      enddo
      write(21,'(a16)')rname(nr)
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(21,'(f12.3,$)')t
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=dreal(uy(lf,ir))
          enddo
        else
          do ir=1,nr
            y(ir)=dimag(uy(lf,ir))
          enddo
        endif
        do ir=1,nr-1
          write(21,'(E16.8,$)')y(ir)
        enddo
        write(21,'(E16.8)')y(nr)
      enddo
      close(21)
c
      open(22,file=ayout,status='unknown')
      write(22,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(22,'(a16,$)')rname(ir)
      enddo
      write(22,'(a16)')rname(nr)
      do ir=1,nr
        y0(ir)=dreal(uy(1,ir))
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=dreal(uy(lf,ir))
          enddo
        else
          do ir=1,nr
            y(ir)=dimag(uy(lf,ir))
          enddo
        endif
        write(22,'(f12.3,$)')t
        do ir=1,nr-1
          write(22,'(E16.8,$)')(y(ir)-y0(ir))/dt
        enddo
        write(22,'(E16.8)')(y(nr)-y0(nr))/dt
        do ir=1,nr
          y0(ir)=y(ir)
        enddo
      enddo
      close(22)
c
c     calculate vertical component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=uz(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf  ))*dexp(-PI2*fi*t)*df
          uz(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=uzout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(uz(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(uz(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
c
      open(21,file=vzout,status='unknown')
      write(21,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(21,'(a16,$)')rname(ir)
      enddo
      write(21,'(a16)')rname(nr)
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(21,'(f12.3,$)')t
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=dreal(uz(lf,ir))
          enddo
        else
          do ir=1,nr
            y(ir)=dimag(uz(lf,ir))
          enddo
        endif
        do ir=1,nr-1
          write(21,'(E16.8,$)')y(ir)
        enddo
        write(21,'(E16.8)')y(nr)
      enddo
      close(21)
c
      open(22,file=azout,status='unknown')
      write(22,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(22,'(a16,$)')rname(ir)
      enddo
      write(22,'(a16)')rname(nr)
      do ir=1,nr
        y0(ir)=dreal(uz(1,ir))
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=dreal(uz(lf,ir))
          enddo
        else
          do ir=1,nr
            y(ir)=dimag(uz(lf,ir))
          enddo
        endif
        write(22,'(f12.3,$)')t
        do ir=1,nr-1
          write(22,'(E16.8,$)')(y(ir)-y0(ir))/dt
        enddo
        write(22,'(E16.8)')(y(nr)-y0(nr))/dt
        do ir=1,nr
          y0(ir)=y(ir)
        enddo
      enddo
      close(22)
c
c     calculate infrasound component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=us(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf  ))*dexp(-PI2*fi*t)*df
          us(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=sdout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(us(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(us(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
c
c     calculate gravity component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=ug(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf  ))*dexp(-PI2*fi*t)*df
          ug(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=grout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(ug(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(ug(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
c
c     calculate north or radial tilt component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=uv(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf  ))*dexp(-PI2*fi*t)*df
          uv(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=tvout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(uv(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(uv(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
c
c     calculate east or tangential tilt component
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,f(lf))
     &     *dcmplx(PI2*tred(ir),0.d0)
          swap(lf)=uw(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          swap(lf)=dconjg(swap(mf))
        enddo
        swap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1(swap,2*nf,+1)
        do lf=1,nf
          t=dble(2*lf-2)*dt
          sr=dreal(swap(2*lf-1))*dexp(-PI2*fi*t)*df
          t=dble(2*lf-1)*dt
          si=dreal(swap(2*lf  ))*dexp(-PI2*fi*t)*df
          uw(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=twout,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      do ir=1,nr
        y(ir)=0.d0
      enddo
      do it=1,ntcutout
        t=togsmin+dble(it-1)*dt
        write(20,'(f12.3,$)')t
        do ir=1,nr-1
          write(20,'(E16.8,$)')y(ir)
        enddo
        write(20,'(E16.8)')y(nr)
        lf=(it+1)/2
        if(it.gt.2*(it/2))then
          do ir=1,nr
            y(ir)=y(ir)+dreal(uw(lf,ir))*dt
          enddo
        else
          do ir=1,nr
            y(ir)=y(ir)+dimag(uw(lf,ir))*dt
          enddo
        endif
      enddo
      close(20)
      return
      end
