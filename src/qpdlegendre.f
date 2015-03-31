      subroutine qpdlegendre(ldeg,raddis)
      implicit none
      integer ldeg
      double precision raddis
c
      include 'qpglobal.h'
c
c     x = cos(raddis)
c     define: Flm(l,m,x) = Plm(l,m,x)/(1-x^2)^(m/2), where Plm are Legendre polynomials
c     This subroutine calculate associated dFlm(l,m,x) = Flm(l,m,x) - Flm(l,m,1)
c
      integer l,m
      double precision x,y
c
      x=dcos(raddis)
      y=-2.d0*(dsin(0.5d0*raddis))**2
c
      plm(0,0)=0.d0
      plm(1,0)=y
      do l=2,ldeg
        plm(l,0)=(dble(2*l-1)*y+dble(2*l-1)*x*plm(l-1,0)
     &             -dble(l-1)*plm(l-2,0))/dble(l)
      enddo
c
      plm(0,1)=0.d0
      plm(1,1)=0.d0
      plm(2,1)=3.d0*y
      do l=3,ldeg
        plm(l,1)=0.5d0*dble(l)*dble(2*l-1)*y
     &          +(dble(2*l-1)*x*plm(l-1,1)
     &           -dble(l)*plm(l-2,1))/dble(l-1)
      enddo
c
      plm(0,2)=0.d0
      plm(1,2)=0.d0
      plm(2,2)=0.d0
      plm(3,2)=15.d0*y
      do l=4,ldeg
        plm(l,2)=0.125d0*dble(l-1)*dble(l)*dble(l+1)*dble(2*l-1)*y
     &          +(dble(2*l-1)*x*plm(l-1,2)
     &           -dble(l+1)*plm(l-2,2))/dble(l-2)
      enddo
c
      return
      end
