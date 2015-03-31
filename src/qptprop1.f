      subroutine qptprop1(ysh,lylw)
      implicit none
c
c     calculation of response to sh source
c     ysh(2,2): solution vector (complex)
c
      integer lylw
      double complex ysh(2,2)
c
      include 'qpglobal.h'
c
c     work space
c
      integer i,istp,ly,key
      double complex momsum,momyup,momylw
      double complex y0(2),yup(2),ylw(2),c(2)
      double complex mat(2,2),mai(2,2),coef(2,2),b(2,2)
c
c===============================================================================
c
c     propagation from surface to source
c
      momsum=(0.d0,0.d0)
      do ly=lyob,lycm-1
        momsum=momsum
     &        +cro(ly)*(crrup(ly)**5-crrlw(ly)**5)/(5.d0,0.d0)
      enddo
c
      yup(1)=(1.d0,0.d0)
      yup(2)=(0.d0,0.d0)
c
      if(lyr.eq.lyob)call cmemcpy(yup,y0,2)
c
      do ly=lyob,lys-1
        call caxcb(mat2x2inv(1,1,ly),yup,2,2,1,c)
        c(1)=c(1)*cdexp(cpt(1,ly))
        c(2)=c(2)*cdexp(cpt(2,ly))
        call caxcb(mat2x2lw(1,1,ly),c,2,2,1,yup)
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,2)
      enddo
c     yup(1)=yup(1)
      yup(2)=yup(2)/crrup(lys)
c
c===============================================================================
c
c     propagation from bottom to source
c
      if(lylw.eq.lycm)then
        ylw(1)=(1.d0,0.d0)
        ylw(2)=(0.d0,0.d0)
      else
        ylw(1)=mat2x2up(1,1,lylw)
        ylw(2)=mat2x2up(2,1,lylw)
      endif
      if(lylw.eq.lyr.and.lylw.gt.lys)then
        call cmemcpy(ylw,y0,2)
      endif
c
      do ly=lylw-1,lys,-1
        call caxcb(mat2x2inv(1,1,ly),ylw,2,2,1,c)
        c(1)=c(1)*cdexp(cpt(1,ly))
        c(2)=c(2)*cdexp(cpt(2,ly))
        call caxcb(mat2x2up(1,1,ly),c,2,2,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,2)
      enddo
c     ylw(1)=ylw(1)
      ylw(2)=ylw(2)/crrup(lys)
c
c     y0(1)=y0(1)
      y0(2)=y0(2)/crrup(lyr)
c
c     momyup = [a^3*Y2(a)-r_s^3*Y2(r_s)]/(-omega^2)
c     momylw = [r_s^3*Y2(r_s)-r_cm^3*Y2(r_cm)]/(-omega^2)
c
      momyup= yup(2)*crrup(lys)**3/comi2
      momylw=-ylw(2)*crrup(lys)**3/comi2
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qptprop1: anormal exit from cdgemp!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ysh(i,istp)=b(1,istp)*y0(i)
          enddo
          ysh(1,istp)=ysh(1,istp)
     &      -(b(1,istp)*momyup+b(2,istp)*momylw)*crrup(lyr)/momsum
        enddo
      else
        do istp=1,2
          do i=1,2
            ysh(i,istp)=b(2,istp)*y0(i)
          enddo
          ysh(1,istp)=ysh(1,istp)
     &      -(b(1,istp)*momyup+b(2,istp)*momylw)*crrup(lyr)/momsum
        enddo
      endif
      return
      end