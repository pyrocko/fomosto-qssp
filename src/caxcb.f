      subroutine caxcb(a,b,n,l,m,c)
      implicit none
c
c     calculate c=a*b
c
      integer n,l,m
      double complex a(n,l),b(l,m),c(n,m)
c
      integer i,j,k
c
      do j=1,m
        do i=1,n
          c(i,j)=dcmplx(0.d0,0.d0)
          do k=1,l
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
          enddo
        enddo
      enddo
      return
      end
