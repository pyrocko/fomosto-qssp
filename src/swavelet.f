      subroutine swavelet(tau,df,nf,wvf)
      implicit none
      integer nf
      double precision tau,df
      double complex wvf(nf)
c
      integer l,n
      double precision f,omi,x,dt0
      double complex alfa,beta,gamma,eta
c
      double precision pi,pi2,eps
      data pi,pi2,eps/3.14159265358979d0,6.28318530717959d0,1.0d-04/
c
c       for wavelet: normalized square half-sinus
c
      do l=1,nf
        f=df*dble(l-1)
        x=f*tau
        if(x.eq.0.d0)then
          wvf(l)=(1.d0,0.d0)
        else if(x.ge.1.d0-eps.and.x.le.1.d0+eps)then
          wvf(l)=dcmplx(-1.d0/x/(1+x),0.d0)
        else
          wvf(l)=dcmplx(0.d0,1.d0/(pi2*x*(1.d0+x)*(1.d0-x)))
     &             *(cdexp(dcmplx(0.d0,-pi2*x))-(1.d0,0.d0))
        endif
      enddo
      return
      end
