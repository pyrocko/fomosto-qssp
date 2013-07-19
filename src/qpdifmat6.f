      subroutine qpdifmat6(ly,ldeg,rr,mat)
      implicit none
      integer ly,ldeg
      double precision rr
      double complex mat(6,6)
c
      include 'qpglobal.h'
c
      double complex cldeg,clp1,cll1,cup,clw
      double complex crr,drr,crorr,clarr,cmurr,cxirr,cgarr,cgrrr
c
      double complex c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      cll1=cldeg*clp1
      crr=dcmplx(rr,0.d0)
      cup=(crr-rrlw(ly))/(crrup(ly)-crrlw(ly))
      clw=c1-cup
      crorr=cup*croup(ly)+clw*crolw(ly)
      clarr=cup*claup(ly)+clw*clalw(ly)
      cmurr=cup*cmuup(ly)+clw*cmulw(ly)
      cgarr=cup*cgaup(ly)+clw*cgalw(ly)
      cgrrr=cup*cgrup(ly)+clw*cgrlw(ly)
      cxirr=clarr+c2*cmurr
c
      mat(1,1)=(c1-c2*clarr/cxirr)/crr
      mat(1,2)=c1/cxirr/crr
      mat(1,3)=cll1*clarr/cxirr/crr
      mat(1,4)=(0.d0,0.d0)
      mat(1,5)=(0.d0,0.d0)
      mat(1,6)=(0.d0,0.d0)
c
      mat(2,1)=c4*(cmurr*(c1+c2*clarr/cxirr)/crr-crorr*cgrrr)
     &        -crorr*comi2*crr
      mat(2,2)=c2*clarr/cxirr/crr
      mat(2,3)=cll1*(crorr*cgrrr
     &        -c2*cmurr*(c1+c2*clarr/cxirr)/crr)
      mat(2,4)=cll1/crr
      mat(2,5)=crorr*clp1*crr
      mat(2,6)=-crorr*crr
c
      mat(3,1)=-c1/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=c2/crr
      mat(3,4)=c1/cmurr/crr
      mat(3,5)=(0.d0,0.d0)
      mat(3,6)=(0.d0,0.d0)
c
      mat(4,1)=crorr*cgrrr-c2*cmurr*(c1+c2*clarr/cxirr)/crr
      mat(4,2)=-clarr/cxirr/crr
      mat(4,3)=c2*cmurr*(c2*cll1*(c1-cmurr/cxirr)-c1)/crr
     &        -crorr*comi2*crr
      mat(4,4)=-c1/crr
      mat(4,5)=-crorr*crr
      mat(4,6)=(0.d0,0.d0)
c
      mat(5,1)=cgarr/crr
      mat(5,2)=(0.d0,0.d0)
      mat(5,3)=(0.d0,0.d0)
      mat(5,4)=(0.d0,0.d0)
      mat(5,5)=-clp1/crr
      mat(5,6)=c1/crr
c
      mat(6,1)=cgarr*clp1/crr
      mat(6,2)=(0.d0,0.d0)
      mat(6,3)=-cgarr*cll1/crr
      mat(6,4)=(0.d0,0.d0)
      mat(6,5)=(0.d0,0.d0)
      mat(6,6)=cldeg/crr
      return
      end