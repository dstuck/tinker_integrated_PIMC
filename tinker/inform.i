c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  inform.i  --  control values for I/O and program flow  ##
c     ##                                                         ##
c     #############################################################
c
c
c     digits    decimal places output for energy and coordinates
c     iprint    steps between status printing (0=no printing)
c     iwrite    steps between coordinate dumps (0=no dumps)
c     isend     steps between socket communication (0=no sockets)
c     verbose   logical flag to turn on extra information
c     debug     logical flag to turn on full debug printing
c     holdup    logical flag to wait for carriage return on exit
c     abort     logical flag to stop execution at next chance
c
c
      integer digits,iprint,iwrite,isend
      logical verbose,debug,holdup,abort
      common /inform/ digits,iprint,iwrite,isend,verbose,debug,holdup,
     &                abort
