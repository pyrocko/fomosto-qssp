      subroutine brunewvlet(tau,df,mm2,wvf)
      implicit none
      integer mm2
      double precision tau,df,fi
      double complex wvf(mm2)
c
      integer l
      double precision f
c
      double precision pi2
      double complex c1
      data pi2/6.28318530717959d0/
      data c1/(1.d0,0.d0)/
c
c     for wavelet: Brune's source time function
c
      do l=1,mm2
        f=df*dble(l-1)
        wvf(l)=c1/dcmplx(1.d0,pi2*f*tau)**2
      enddo
      return
      end