      subroutine qpsmat0(ly,lwup)
      implicit none
c
c     calculate 3x3 spheroidal layer matrix for a solid shell in case of degree l = 0
c
      integer ly,lwup
c
      include 'qpglobal.h'
c
      integer i,j
      double complex cmudr,c2mu,cxi,cx2,cr2,cdet
      double complex ph0jup,ph1jup,ph0yup,ph1yup,ps0jup,ps0yup
      double complex ph0jlw,ph1jlw,ph0ylw,ph1ylw,ps0jlw,ps0ylw
      double complex spb(2,2)
      double complex spbphj,spbphy,spbpsj,spbpsy
c
      double complex ci,c1,c2,c3,c6
      data ci,c1,c2,c3,c6/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),
     &                    (3.d0,0.d0),(6.d0,0.d0)/
c
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
c
      cx2=(kp(ly)*crrup(ly))**2
      if(ksmallpsv(ly))then
        ph0jup=spbphj(0,kp(ly)*crrup(ly))
        ph1jup=spbphj(1,kp(ly)*crrup(ly))
        ps0jup=spbpsj(0,kp(ly)*crrup(ly))
c
        ph0yup=spbphy(0,kp(ly)*crrup(ly))
        ph1yup=spbphy(1,kp(ly)*crrup(ly))
        ps0yup=spbpsy(0,kp(ly)*crrup(ly))
c
        cr2=crrup(ly)**2
        mas3x3up(1,1,ly)=cr2*ph1jup/c3
        mas3x3up(2,1,ly)=cxi*cr2*ph0jup-c2*c2mu*cr2*ph1jup/c3
        mas3x3up(3,1,ly)=cga(ly)*cr2*ps0jup/c6
        mas3x3up(1,2,ly)=ph1yup
        mas3x3up(2,2,ly)=cxi*cx2*ph0yup-c2*c2mu*ph1yup
        mas3x3up(3,2,ly)=cga(ly)*cr2*ps0yup/c2
      else
        spb(1,1)=c1-zjup(0,ly,1)
        spb(2,1)=-ci*zjup(0,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhup(0,ly,1)
c
        do j=1,2
          mas3x3up(1,j,ly)=spb(2,j)
          mas3x3up(2,j,ly)=cxi*cx2*spb(1,j)-c2*c2mu*spb(2,j)
          mas3x3up(3,j,ly)=-cga(ly)*spb(1,j)
        enddo
      endif
      mas3x3up(1,3,ly)=(0.d0,0.d0)
      mas3x3up(2,3,ly)=(0.d0,0.d0)
      mas3x3up(3,3,ly)=c1
c
      if(ly.eq.lylwpsv(0))return
c
      cx2=(kp(ly)*crrlw(ly))**2
      if(ksmallpsv(ly))then
        ph0jlw=spbphj(0,kp(ly)*crrlw(ly))
        ph1jlw=spbphj(1,kp(ly)*crrlw(ly))
        ps0jlw=spbpsj(0,kp(ly)*crrlw(ly))
c
        ph0ylw=spbphy(0,kp(ly)*crrlw(ly))
        ph1ylw=spbphy(1,kp(ly)*crrlw(ly))
        ps0ylw=spbpsy(0,kp(ly)*crrlw(ly))
c
        cr2=crrlw(ly)**2
        mas3x3lw(1,1,ly)=cr2*ph1jlw/c3
        mas3x3lw(2,1,ly)=cxi*cr2*ph0jlw-c2*c2mu*cr2*ph1jlw/c3
        mas3x3lw(3,1,ly)=cga(ly)*cr2*ps0jlw/c6
        mas3x3lw(1,2,ly)=ph1ylw
        mas3x3lw(2,2,ly)=cxi*cx2*ph0ylw-c2*c2mu*ph1ylw
        mas3x3lw(3,2,ly)=cga(ly)*cr2*ps0ylw/c2
      else
        spb(1,1)=c1-zjlw(0,ly,1)
        spb(2,1)=-ci*zjlw(0,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhlw(0,ly,1)
c
        do j=1,2
          mas3x3lw(1,j,ly)=spb(2,j)
          mas3x3lw(2,j,ly)=cxi*cx2*spb(1,j)-c2*c2mu*spb(2,j)
          mas3x3lw(3,j,ly)=-cga(ly)*spb(1,j)
        enddo
      endif
      mas3x3lw(1,3,ly)=(0.d0,0.d0)
      mas3x3lw(2,3,ly)=(0.d0,0.d0)
      mas3x3lw(3,3,ly)=c1
c
      if(lwup.eq.0)then
c
c       calculate inverse masrix at upper radius
c
        cdet=mas3x3up(1,1,ly)*mas3x3up(2,2,ly)
     &      -mas3x3up(1,2,ly)*mas3x3up(2,1,ly)
        mas3x3inv(1,1,ly)=mas3x3up(2,2,ly)/cdet
        mas3x3inv(1,2,ly)=-mas3x3up(1,2,ly)/cdet
        mas3x3inv(1,3,ly)=(0.d0,0.d0)
        mas3x3inv(2,1,ly)=-mas3x3up(2,1,ly)/cdet
        mas3x3inv(2,2,ly)=mas3x3up(1,1,ly)/cdet
        mas3x3inv(2,3,ly)=(0.d0,0.d0)
        mas3x3inv(3,1,ly)=(mas3x3up(2,1,ly)*mas3x3up(3,2,ly)
     &                    -mas3x3up(3,1,ly)*mas3x3up(2,2,ly))/cdet
        mas3x3inv(3,2,ly)=(mas3x3up(1,2,ly)*mas3x3up(3,1,ly)
     &                    -mas3x3up(1,1,ly)*mas3x3up(3,2,ly))/cdet
        mas3x3inv(3,3,ly)=c1
      else
c
c       calculate inverse masrix at lower radius
c
        cdet=mas3x3lw(1,1,ly)*mas3x3lw(2,2,ly)
     &      -mas3x3lw(1,2,ly)*mas3x3lw(2,1,ly)
        mas3x3inv(1,1,ly)=mas3x3lw(2,2,ly)/cdet
        mas3x3inv(1,2,ly)=-mas3x3lw(1,2,ly)/cdet
        mas3x3inv(1,3,ly)=(0.d0,0.d0)
        mas3x3inv(2,1,ly)=-mas3x3lw(2,1,ly)/cdet
        mas3x3inv(2,2,ly)=mas3x3lw(1,1,ly)/cdet
        mas3x3inv(2,3,ly)=(0.d0,0.d0)
        mas3x3inv(3,1,ly)=(mas3x3lw(2,1,ly)*mas3x3lw(3,2,ly)
     &                    -mas3x3lw(3,1,ly)*mas3x3lw(2,2,ly))/cdet
        mas3x3inv(3,2,ly)=(mas3x3lw(1,2,ly)*mas3x3lw(3,1,ly)
     &                    -mas3x3lw(1,1,ly)*mas3x3lw(3,2,ly))/cdet
        mas3x3inv(3,3,ly)=c1
      endif
c
      return
      end