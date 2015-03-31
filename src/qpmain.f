      program qssp
      implicit none
c
      include 'qpglobal.h'
c
c     work space
c
      integer ig,ierr,runtime
      integer time
      character*80 inputfile
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#    AA    TTTTT  M   M   QQQ    SSSS   SSSS  PPPP   #'
      print *,'#   A  A     T    MM MM  Q   Q  S      S      P   P  #'
      print *,'#   AAAA     T    M M M  Q   Q   SSS    SSS   PPPP   #'
      print *,'#  A    A    T    M   M  Q  QQ      S      S  P      #'
      print *,'#  A    A    T    M   M   QQQQ  SSSS   SSSS   P      #'
      print *,'#                                                    #'
      print *,'#                  (Version 2010)                    #'
      print *,'#                                                    #'
      print *,'#              synthetic seismograms                 #'
      print *,'#                      based on                      #'
      print *,'#          a spherically symmetric earth model       #'
      print *,'#                                                    #'
      print *,'#                  (Version 2010)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#            last modified: April 2010               #'
      print *,'######################################################'
      print *,'                          '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qpgetinp(10)
      close(10)
c
      call qpsublayer(ierr)
      do ig=1,ngrn
        if(grnsel(ig).eq.1)then
          lys=lygrn(ig)
          call qpgrnspec(ig)
        endif
      enddo
      call qpwvint(ierr)
      call qpfftinv(ierr)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with atmqssp     #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
