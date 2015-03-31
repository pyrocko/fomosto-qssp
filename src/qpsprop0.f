      subroutine qpsprop0(ypsv,lyup,lylw)
      implicit none
c
c     calculation of response to psv source (l=0)
c     ypsv(6,4): solution vector (complex)
c
      integer lyup,lylw
      double complex ypsv(6,4)
c
      include 'qpglobal.h'
c
c     work space
c
      integer i,istp,ly,key
      double complex y0(3),c(3),yup(3),ylw(3)
      double complex wave(3),coef(2,2),b(2,2)
c
c===============================================================================
c
c     propagation from surface to source
c
      yup(1)=(1.d0,0.d0)
      yup(2)=(0.d0,0.d0)
      yup(3)=(0.d0,0.d0)
c
      if(lyr.eq.lyup)call cmemcpy(yup,y0,3)
c
      do ly=lyup,lys-1
        call caxcb(mas3x3inv(1,1,ly),yup,3,3,1,c)
        wave(1)=cdexp( cps(1,ly))
        wave(2)=cdexp(-cps(2,ly))
        wave(3)=(1.d0,0.d0)
        if(ly.ge.lyr)then
          do i=1,3
            y0(i)=y0(i)*wave(2)
          enddo
        endif
c
	  c(1)=c(1)*wave(1)*wave(2)
c       c(2)=c(2)
        c(3)=c(3)*wave(3)*wave(2)
c
        call caxcb(mas3x3lw(1,1,ly),c,3,3,1,yup)
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,3)
      enddo
      yup(1)=yup(1)/crrup(lys)
      yup(2)=yup(2)/crrup(lys)**2
c     yup(3)=yup(3)
c
c===============================================================================
c
c     propagation from bottom to source
c
      do i=1,3
        ylw(i)=mas3x3up(i,1,lylw)
      enddo
      if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,3)
c
      do ly=lylw-1,lys,-1
        call caxcb(mas3x3inv(1,1,ly),ylw,3,3,1,c)
        wave(1)=cdexp(-cps(1,ly))
        wave(2)=cdexp( cps(2,ly))
        wave(3)=(1.d0,0.d0)
        if(ly.lt.lyr)then
          do i=1,3
            y0(i)=y0(i)*wave(1)
          enddo
        endif
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
        c(3)=c(3)*wave(1)*wave(3)
c
        call caxcb(mas3x3up(1,1,ly),c,3,3,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,3)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      y0(1)=y0(1)/crrup(lyr)
      y0(2)=y0(2)/crrup(lyr)**2
c     y0(3)=y0(3)
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
        print *,' Warning in qpsprop0: anormal exit from cdgemp!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(1,istp)*y0(i)
          enddo
          ypsv(5,istp)=b(1,istp)*y0(3)
          ypsv(6,istp)=ypsv(5,istp)/crrup(lyr)
        enddo
      else
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(2,istp)*y0(i)
          enddo
c
c         add a constant to the potential below the source so that
c         it becomes continuous through the source level.
c
          ypsv(5,istp)=b(2,istp)*(y0(3)-ylw(3))+b(1,istp)*yup(3)
          ypsv(6,istp)=ypsv(5,istp)/crrup(lyr)
        enddo
      endif
      return
      end
