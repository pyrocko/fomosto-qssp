      subroutine qplegendre(ldeg,raddis)
      implicit none
      integer ldeg
      double precision raddis
c
      include 'qpglobal.h'
c
c     calculate associated Plm(l,m,x)/(1-x^2)^(m/2)
c     where Plm are Legendre polynomials
c
      integer l,m
      double precision x
c
      x=dcos(raddis)
      plm(0,0)=1.d0
      plm(0,1)=0.d0
      plm(1,1)=1.d0
      plm(0,2)=0.d0
      plm(1,2)=0.d0
      plm(2,2)=3.d0
      do m=0,2
        plm(m+1,m)=dble(2*m+1)*x*plm(m,m)
        do l=m+2,ldeg
          plm(l,m)=(dble(2*l-1)*x*plm(l-1,m)
     &             -dble(l+m-1)*plm(l-2,m))/dble(l-m)
        enddo
      enddo
      return
      end
