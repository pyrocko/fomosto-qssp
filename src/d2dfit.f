      subroutine d2dfit(n,dis,p,b,ndeg,disfit)
      implicit none
      integer n,ndeg
      double precision dis(n),disfit(n),p(n,0:ndeg),b(0:ndeg)
c
      integer i,ideg
      double precision x,dx
c
      dx=2.d0/dble(n-1)
      do i=1,n
        x=-1.d0+dble(i-1)*dx
        p(i,0)=1.d0
        p(i,1)=x
        do ideg=2,ndeg
          p(i,ideg)=(dble(2*ideg-1)*x*p(i,ideg-1)
     &              -dble(ideg-1)*p(i,ideg-2))/dble(ideg)
        enddo
      enddo
c
      do ideg=0,ndeg
        b(ideg)=0.5d0*(dis(1)*p(1,ideg)+dis(n)*p(n,ideg))
        do i=2,n-1
          b(ideg)=b(ideg)+dis(i)*p(i,ideg)
        enddo
        b(ideg)=b(ideg)*dx*0.5d0*dble(2*ideg+1)
      enddo
c
      do i=1,n
        x=-1.+dble(i-1)*dx
        disfit(i)=0.d0
        do ideg=0,ndeg
          disfit(i)=disfit(i)+b(ideg)*p(i,ideg)
        enddo
      enddo
      return
      end