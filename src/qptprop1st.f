      subroutine qptprop1st(ysh)
      implicit none
c
c     calculation of static response to sh source of degree 1
c     ysh(2,2): solution vector (complex)
c
      double complex ysh(2,2)
c
      include 'qpglobal.h'
c
c     work space
c
      integer ly
      double complex momup,momlw
c
c===============================================================================
c
c     rigid rotation: y ~ r
c
      momup=(0.d0,0.d0)
      do ly=lyob,lys-1
        momup=momup+cro(ly)*(crrup(ly)**5-crrlw(ly)**5)/(5.d0,0.d0)
      enddo
c
      momlw=(0.d0,0.d0)
      do ly=lys,lycm-1
        momlw=momlw+cro(ly)*(crrup(ly)**5-crrlw(ly)**5)/(5.d0,0.d0)
      enddo
c
      if(lyr.le.lys)then
        ysh(1,1)= momlw*(crrup(lyr)/crrup(lys))/(momup+momlw)
      else
        ysh(1,1)=-momup*(crrup(lyr)/crrup(lys))/(momup+momlw)
      endif
      ysh(2,1)=(0.d0,0.d0)
      ysh(1,2)=(0.d0,0.d0)
      ysh(2,2)=(0.d0,0.d0)
c
      return
      end