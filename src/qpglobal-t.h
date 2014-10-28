c     CONSTANTS
c     =========
      double precision PI,PI2
      parameter(PI=3.14159265358979d0,PI2=6.28318530717959d0)
      double precision DEG2RAD,KM2M
      parameter(DEG2RAD=1.745329251994328d-02,KM2M=1.0d+03)
      double precision REARTH,BIGG
      parameter(REARTH=6.371d+06,BIGG=6.6732d-11)
      double precision RESOLUT
      parameter(RESOLUT=0.01d0)
      double precision FLTAPER
      parameter(FLTAPER=0.2d0)
      double precision FSBUP,FSBLW,FSBREF
      parameter(FSBUP=0.1d+02,FSBLW=0.5d-03,FSBREF=1.0d+00)

c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     lymax: max. number of layers
c     nrmax: max. number of receivers
c     nsmax: max. number of sub-events
c     ngrnmax: max. number of source depths (Green's function data bases)
c     nfmax: max. number of frequency samples
c     ldegmax: max. degree of Legendre polynomials
c
      integer lymax,nrmax,nsmax,ngrnmax,nfmax,ldegmax,ndmax
      parameter(nrmax=21,nsmax=1500,ngrnmax=50)
      parameter(nfmax=4096,ldegmax=10000)
      parameter(lymax=220+ngrnmax)
      parameter(ndmax=4)
c
c     GREEN FUNCTION PARAMETERS
c     =========================
c
c     ntcut = number of time samples
c     ntcutout = number of time samples of output
c     nt = power of 2 integer nearest to ntcut
c     nf = nt/2
c     nfcut = int(fcut/dt)
c     lygrn = layer number of Green's function source
c     grnsel = selection of Green's function to be calculated
c     dt = time sampling interval
c     df = frequency sampling interval
c     fi = imaginary frequency derived from anti-aliasing factor
c     fcut = nfcut*df
c     fgr = critical frequency from "with" to "without gravitation"
c     ldeggr = critical harmonic degree
c     comi = complex angular frequency
c     comi2 = comi**2
c     yr,yt,yp,ys,yg,yv,yw = Green's function spectra
c     (r,t,p = displacement, s = sound, g = gravity, v,w = tilt)
c     grnfile = file name of Green's function spectra
c
      logical nogravity,selpsv,selsh
      integer ngrn,nt,ntcut,ntcutout,nf,nfcut,nlpf,ldeggr,ldegmin
      integer lygrn(ngrnmax),grnsel(ngrnmax)
      integer ldegpsv(lymax),ldegsh(lymax)
      double precision dt,df,fi,fcut,fgr,rratmos,depatmos,rr0
      double precision grndep(ngrnmax)
      double complex comi,comi2
      complex yr(0:ldegmax,4,0:ndmax),yt(0:ldegmax,4,0:ndmax)
      complex yp(0:ldegmax,4,0:ndmax),ys(0:ldegmax,4,0:ndmax)
      complex yg(0:ldegmax,4,0:ndmax),yv(0:ldegmax,4,0:ndmax)
      complex yw(0:ldegmax,4,0:ndmax)
      character*800 grnfile(ngrnmax),
     &             rgrnfile(ngrnmax),tgrnfile(ngrnmax),
     &             pgrnfile(ngrnmax),sgrnfile(ngrnmax),
     &             ggrnfile(ngrnmax),vgrnfile(ngrnmax),
     &             wgrnfile(ngrnmax)
      common /lgreen/ nogravity,selpsv,selsh
      common /igreen/ ngrn,nt,ntcut,ntcutout,nf,nfcut,nlpf,
     &                ldeggr,ldegmin,lygrn,grnsel,ldegpsv,ldegsh
      common /dgreen/ dt,df,fi,fcut,fgr,rratmos,depatmos,rr0,
     &                grndep,comi,comi2
      common /rgreen/ yr,yt,yp,ys,yg,yv,yw
      common /cgreen/ rgrnfile,tgrnfile,pgrnfile,sgrnfile,ggrnfile,
     &                vgrnfile,wgrnfile,grnfile
c
c     slowness cut-offs
c
      double precision slwmax,slwlwcut,slwupcut,fcorner
      common /dslwcutoffs/ slwmax,slwlwcut,slwupcut,fcorner
c
c     RECEIVER PARAMETERS
c     ===================
c
c     nr = number of receivers
c     ioutform = output format:
c                1 = vertical(z)/north(n)/east(e)
c                2 = vertical(z)/radial(r)/tangential(t)
c     latr,lonr = geographic receiver coordinates
c     tred = time reduction
c     rname = receiver name
c     dpr = uniform receiver depth
c     uxout, ... = synthetic seismogram outputs
c          (u = displacement,v = velocity,a = acceleration,gr = gravimetric)
c
      integer nr,ioutform,ishypo
      double precision dpr,freeairgrd
      double precision latr(nrmax),lonr(nrmax),tred(nrmax)
      double complex uz(nfmax,nrmax),ux(nfmax,nrmax)
      double complex uy(nfmax,nrmax),us(nfmax,nrmax)
      double complex ug(nfmax,nrmax),uv(nfmax,nrmax)
      double complex uw(nfmax,nrmax)
      character*10 rname(nrmax)
      character*800 uxout,uyout,uzout,
     &             vxout,vyout,vzout,
     &             axout,ayout,azout,
     &             sdout,grout,tvout,
     &             twout
      common /ireceiver/ nr,ioutform,ishypo
      common /dreceiver/ dpr,freeairgrd,latr,lonr,tred,
     &                   uz,ux,uy,us,ug,uv,uw
      common /creceiver/ rname,sdout,grout,
     &             uxout,uyout,uzout,
     &             vxout,vyout,vzout,
     &             axout,ayout,azout,
     &             tvout,twout
c
c     SUBEVENT PARAMETERS
c     ===================
c
c     ns = number of sub-events
c     isg1,isg2 = the respective Green's function number of subevents
c     mtt,... = moment tensor
c     lats,lons = geographic coordinates
c     togs = origin time
c     trss = rise time
c
      integer ns
      integer isg1(ngrnmax),isg2(ngrnmax),nsg(ngrnmax)
      double precision togsmin
      double precision mtt(nsmax),mpp(nsmax),mrr(nsmax),
     &                 mtp(nsmax),mpr(nsmax),mrt(nsmax),
     &                 lats(nsmax),lons(nsmax),deps(nsmax),
     &                 togs(nsmax),trss(nsmax)
      common /isubevents/ns,isg1,isg2,nsg
      common /dsubevents/mtt,mpp,mrr,mtp,mpr,mrt,
     &                   lats,lons,deps,togs,trss,togsmin
c
c     SOURCE-RECEIVER CONFIGURATION PARAMETERS
c     ========================================
c
      integer idr(nsmax,nrmax)
      double precision dis(nsmax,nrmax)
      double complex ssa(nsmax,nrmax),ss2a(nsmax,nrmax)
      double complex csa(nsmax,nrmax),cs2a(nsmax,nrmax)
      double complex ssb(nsmax,nrmax),csb(nsmax,nrmax)
      double complex ssd(nsmax,nrmax),ssf(nsmax,nrmax)
      common /isourcerec/ idr
      common /dsourcerec/ dis,ssa,ss2a,csa,cs2a,ssb,csb,ssd,ssf
c
c     LEGENDRE POLYNOMIALS TABLES
c     ===========================
c
c     plm = associated Legendre polynomials divided by sin(x)**m
c
      double precision plm(0:ldegmax,0:2)
      common /dlegendre/ plm
c
c     ORIGINAL MODEL PARAMETERS
c     =========================
c
c     disperion = yes/no
c     lys = layer number of source
c     lyr = layer number of receiver
c     lyob = layer number of ocean/atmosphere bottom
c     lycm = layer number of core-mantle boundary
c     lycc = layer number of inner and outer core boundary
c     ly0 = max. layer number for integration
c     lylwpsv = starting layer number for psv solution
c     lylwsh = starting layer number for sh solution
c
      logical dispersion
      integer lys,lyr,lyos,lyob,lycm,lycc,ly0
      integer lylwsh(0:ldegmax),lylwpsv(0:ldegmax)
      integer lyupatm(0:ldegmax)
      common /ldispers/ dispersion
      common /ilnumber/ lys,lyr,lyos,lyob,lycm,lycc,ly0,
     &                  lylwsh,lylwpsv,lyupatm
c
c     qsmin = min. quality factor for shear wave
c     dp = depth (up = top of layer, lw = bottom of layer)
c     vp, vs = p and s wave velocity
c     ro = density
c     qp, qs = quality factors of p and s waves
c
      integer l0
      double precision qsmin
      double precision dp0(2*lymax),dp0up(2*lymax),dp0lw(2*lymax)
      double precision vp0(2*lymax),vp0up(2*lymax),vp0lw(2*lymax)
      double precision vs0(2*lymax),vs0up(2*lymax),vs0lw(2*lymax)
      double precision ro0(2*lymax),ro0up(2*lymax),ro0lw(2*lymax)
      double precision qp0(2*lymax),qp0up(2*lymax),qp0lw(2*lymax)
      double precision qs0(2*lymax),qs0up(2*lymax),qs0lw(2*lymax)
      common /imodel0/ l0
      common /dmodel0/ qsmin,
     &                 dp0,dp0up,dp0lw,vp0,vp0up,vp0lw,
     &                 vs0,vs0up,vs0lw,ro0,ro0up,ro0lw,
     &                 qp0,qp0up,qp0lw,qs0,qs0up,qs0lw
c
c     rr = radius (up = top of layer, lw = bottom of layer)
c     else see above.
c
      double precision rrup(lymax),rrlw(lymax)
      double precision vpup(lymax),vplw(lymax)
      double precision vsup(lymax),vslw(lymax)
      double precision roup(lymax),rolw(lymax)
      double precision qpup(lymax),qplw(lymax)
      double precision qsup(lymax),qslw(lymax)
      common /dmodel/ rrup,rrlw,vpup,vplw,
     &                vsup,vslw,roup,rolw,
     &                qpup,qplw,qsup,qslw
c
c     LAYER MATRICES
c     ==============
c
c     la, mu = Lame constants
c     gr = gravity
c     ga = 4*pi*G*rho
c     (up = top of layer, lw = bottom of layer, without = average)
c
      double complex crrup(lymax),crrlw(lymax)
      double complex cro(lymax),croup(lymax),crolw(lymax)
      double complex cla(lymax),claup(lymax),clalw(lymax)
      double complex cmu(lymax),cmuup(lymax),cmulw(lymax)
      double complex cvp(lymax),cvpup(lymax),cvplw(lymax)
      double complex cvs(lymax),cvsup(lymax),cvslw(lymax)
      double complex cgr(lymax),cgrup(lymax),cgrlw(lymax)
	  double complex cga(lymax),cgaup(lymax),cgalw(lymax)
      common /cmodel/ crrup,crrlw,cro,croup,crolw,
     &                cla,claup,clalw,cmu,cmuup,cmulw,
     &                cvp,cvpup,cvplw,cvs,cvsup,cvslw,
     &                cgr,cgrup,cgrlw,cga,cgaup,cgalw
c
c     ksmallpsv = case for very small enough psv wave number
c     ksmallsh = case for very small enough sh wave number
c
      logical ksmallpsv(lymax),ksmallsh(lymax)
      common /spblist/ ksmallpsv,ksmallsh
c
c     kp = omi/vp, ks = omi/vs (sv), kt = omi/vs (sh)
c     cps = phase term of psv propagator
c     cpt = phase term of sh propagator
c     mat2x2 = 2x2 toroidal solution matrix
c     mas2x2 = 2x2 spheroidal solution matrix (l = 0)
c     mas4x4 = 4x4 spheroidal solution matrix (l > 0, in liquid)
c     mas6x6 = 6x6 spheroidal solution matrix (l > 0, in solid)
c     mas(t)inv = inverse solution matrix
c
      double complex kp(lymax),ks(lymax),kt(lymax)
      double complex cps(6,lymax),cpt(2,lymax)
      double complex mat2x2up(2,2,lymax),mat2x2lw(2,2,lymax)
      double complex mat2x2inv(2,2,lymax)
	  double complex mas2x2up(2,2,lymax),mas2x2lw(2,2,lymax)
      double complex mas2x2inv(2,2,lymax)
	  double complex mas4x4up(4,4,lymax),mas4x4lw(4,4,lymax)
      double complex mas4x4inv(4,4,lymax)
	  double complex mas6x6up(6,6,lymax),mas6x6lw(6,6,lymax)
      double complex mas6x6inv(6,6,lymax)
      common /matrix/ kp,ks,kt,cps,cpt,
     &                mat2x2up,mat2x2lw,mat2x2inv,
     &                mas2x2up,mas2x2lw,mas2x2inv,
     &                mas4x4up,mas4x4lw,mas4x4inv,
     &                mas6x6up,mas6x6lw,mas6x6inv
c
c     SPHERICAL BESSEL FUNCTIONS
c     ==========================
c
c     zj(x)=x*j_(n+1)(x)/[j_n(x)+ix*j_(n+1)(x)]
c     wj(xa,xb)=ln{[j_n(xa)+i*xa*j_(n+1)(xa)]/[j_n(xb)+i*xb*j_(n+1)(xb)]}
c     zh(x)=x*h_(n+1)(x)/h_n(x), where h_n(x) = j_n(x) + sign[im(x)]iy_n(x)
c     wh(xa,xb)=ln[h_n(xa)/h_n(xb)]
c
      double complex zjup(0:ldegmax,lymax,3),zjlw(0:ldegmax,lymax,3)
      double complex zhup(0:ldegmax,lymax,3),zhlw(0:ldegmax,lymax,3)
      double complex wj(0:ldegmax,lymax,3),wh(0:ldegmax,lymax,3)
      common /spbessels0/ zjup,zjlw,zhup,zhlw,wj,wh
c
c     same as above but with gravity effect
c
      double complex zjupg(0:ldegmax),zjlwg(0:ldegmax)
      double complex zhupg(0:ldegmax),zhlwg(0:ldegmax)
      double complex wjg(0:ldegmax),whg(0:ldegmax)
      common /spbesselsg/ zjupg,zjlwg,zhupg,zhlwg,wjg,whg
c
c     cua = solution describing rigid motion of the earth
c
      double complex cua(8,6),cypnorm(6,lymax)
      common /masscenteracc/ cua,cypnorm
c
c     mmantle = inertia moment of mantle
c     mshell = inertia moment of shells in mantle
c
      double precision mmantle,mshell(lymax)
      common /mantlemoment/ mmantle,mshell
