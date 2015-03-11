c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine dstuckOptimize                              ##
c     ##                                                         ##
c     #############################################################
c
c
c
c
      subroutine dstuckOptimize(typeVec, xyzCoord, connMat, modeNum, 
     $partNum)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
      include 'rigid.i'
      include 'vdwpot.i'
c DES: from initial.f
c      include 'sizes.i'
      include 'align.i'
      include 'atoms.i'
      include 'bath.i'
c      include 'bound.i'
      include 'cell.i'
      include 'files.i'
      include 'group.i'
      include 'inform.i'
c      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'minima.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'neigh.i'
      include 'openmp.i'
      include 'output.i'
      include 'params.i'
      include 'pdb.i'
      include 'precis.i'
c      include 'rigid.i'
      include 'scales.i'
      include 'sequen.i'
      include 'socket.i'
      include 'warp.i'
      include 'zclose.i'
c DES: from optimize
      include 'usage.i'
c DES out
c DES: from mechanic.f misc calls
      include 'boxes.i'
c DES out
c DES: from readxyz.f
      include 'atmtyp.i'
      include 'couple.i'
      include 'titles.i'
c DES out
      integer modeNum
      integer partNum
      integer, dimension(partNum) :: typeVec
      real*8, dimension(partNum*3) :: xyzCoord
      real*8, dimension(partNum*8) :: connMat
      real*8, dimension(modeNum) :: redMass
c DES: from initial.f
      real*8 precise
c DES out
c DES: from getxyz.f
      integer ixyz
      integer freeunit
      logical exist
      character*120 xyzfile
c DES out
c DES: from optimize
      integer i,j,imin,nvar
      integer next
      real*8 minimum,optimiz1
      real*8 grdmin,gnorm,grms
      real*8 energy,eps
      real*8, allocatable :: xx(:)
      real*8, allocatable :: derivs(:,:)
      logical analytic
      character*20 keyword
      character*120 minfile
      character*120 record
      character*120 string
      external energy
      external optimiz1
      external optsave
c DES out
c DES: for projtrm
      integer hcounter
      real*8 xc(3,N)
      real*8 hfull(3*N,3*N)
      real*8 trvec(3*N,3*N)
      real*8 temp(3*N,3*N)
c DES out
c DES: from readxyz.f
      integer k,m
      logical clash
c DES out

c DES: from initial.f
c     default unit numbers for input and output
c
      input = 5
      iout = 6
c     cores, thread count and options for OpenMP
c
      nproc = 1
      nthread = 1
!$    nproc = omp_get_num_procs ()
!$    nthread = nproc
!$    call omp_set_num_threads (nthread)
!$    call omp_set_nested (.true.)
c
c     Intel compiler extensions to OpenMP standard
c
!$    call kmp_set_stacksize (2**28)
!$    call kmp_set_blocktime (0)
c
c     values of machine precision constants
c
      tiny = precise (1)
      small = precise (2)
      huge = precise (3)
c
c     number of lines in the keyfile
c
      nkey = 0
c
c     number of lines in the parameter file
c
      nprm = 0
c
c     number of atoms in the system
c
      n = 0
c
c     number of molecules in the system
c
      nmol = 0
c
c     number of unit cell replicates
c
      ncell = 0
c
c     number of atoms used in superposition
c
      nfit = 0
c
c     number of mutated atoms in the system
c
      nmut = 0
c
c     number of bonds added or deleted from Z-matrix
c
      nadd = 0
      ndel = 0
c
c     number of atoms in Protein Data Bank format
c
      npdb = 0
c
c     number of residues and chains in biopolymer sequence
c
      nseq = 0
      nchain = 0
c
c     highest numbered previous cycle file
c
      nprior = 0
c
c     flags for information levels within the program
c
      verbose = .false.
      debug = .false.
      abort = .false.
c
c     flag for use of atom groups
c
      use_group = .false.
c
c     flags for periodic boundaries
c
      use_bounds = .false.
      use_replica = .false.
      use_polymer = .false.
c
c     flags for temperature and pressure baths
c
      isothermal = .false.
      isobaric = .false.
c
c     flags for rebuilding of neighbor lists
c
      dovlst = .true.
      doclst = .true.
      domlst = .true.
c
c     flag for use of rigid bodies
c
      use_rigid = .false.
c
c     flag to show setting of optimization scale factors
c
      set_scale = .false.
c
c     flags for external Java socket communication
c
      skt_init = .false.
      use_socket = .false.
c
c     flags for potential energy smoothing
c
      use_smooth = .false.
      use_dem = .false.
      use_gda = .false.
      use_tophat = .false.
      use_stophat = .false.
c
c     type of coordinates file
c
      coordtype = 'NONE'
c
c     atomic symbols for elements
c
      call initatom
c
c     default parameters used by optimizations
c
      fctmin = 0.0d0
      maxiter = 0
      nextiter = 0
      iprint = -1
      iwrite = -1
      stpmax = 0.0d0
c DES out

c DES: from getxyz.f
      coordtype = 'CARTESIAN'
c DES: from getxyz.f
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'CARTESIAN'
c DES: readxyz -> replaced by directly passing info
c
c     initialize the total number of atoms in the system
c
      n = partNum
      title = ' '
      ltitle = 0
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         tag(i) = 0
         name(i) = '   '
         x(i) = xyzCoord(3*(i-1)+1)
         y(i) = xyzCoord(3*(i-1)+2)
         z(i) = xyzCoord(3*(i-1)+3)
         type(i) = typeVec(i)
         do j = 1, 8
            i12(j,i) = connMat(j+(i-1)*8)
         end do
      end do
      
c
c     for each atom, count and sort its attached atoms
c
      do i = 1, n
         n12(i) = 0
         do j = 8, 1, -1
            if (i12(j,i) .ne. 0) then
               n12(i) = j
               goto 90
            end if
         end do
   90    continue
         call sort (n12(i),i12(1,i))
      end do
c
c     check for atom pairs with identical coordinates
c
      clash = .false.
      if (n .le. 10000)  call chkxyz (clash)
c
c     make sure that all connectivities are bidirectional
c
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            do m = 1, n12(k)
               if (i12(m,k) .eq. i)  goto 120
            end do
            write (iout,110)  k,i
  110       format (/,' READXYZ  --  Check Connection of Atom',
     &                 i6,' to Atom',i6)
            call fatal
  120       continue
         end do
      end do
c DES end of readxyz
c DES out

c DES from optimize.f
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
c
c     search the keywords for output frequency parameters
c
c      do i = 1, nkey
c         next = 1
c         record = keyline(i)
c         call gettext (record,keyword,next)
c         call upcase (keyword)
c         string = record(next:120)
c         if (keyword(1:9) .eq. 'PRINTOUT ') then
c            read (string,*,err=10,end=10)  iprint
c         else if (keyword(1:9) .eq. 'WRITEOUT ') then
c            read (string,*,err=10,end=10)  iwrite
c         end if
c   10    continue
c      end do
c
c     get termination criterion as RMS gradient per atom
c
      grdmin = 0.0001d0
c
c     write out a copy of coordinates for later update
c
c      imin = freeunit ()
c      minfile = filename(1:leng)//'.xyz'
c      call version (minfile,'new')
c      open (unit=imin,file=minfile,status='new')
c      call prtxyz (imin)
c      close (unit=imin)
c      outfile = minfile
c
c     set scaling parameter for function and derivative values;
c     use square root of median eigenvalue of typical Hessian
c
      set_scale = .true.
      nvar = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               nvar = nvar + 1
               scale(nvar) = 12.0d0
            end do
         end if
      end do
c
c     check for too many parameters to be optimized
c
      if (nvar .gt. maxopt) then
         write (iout,50)
   50    format (/,' OPTIMIZE  --  Too many Parameters,',
     &              ' Increase the Value of MAXOPT')
         call fatal
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvar))
      allocate (derivs(3,n))
c
c     scale the coordinates of each active atom
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i) * scale(nvar)
            nvar = nvar + 1
            xx(nvar) = y(i) * scale(nvar)
            nvar = nvar + 1
            xx(nvar) = z(i) * scale(nvar)
         end if
      end do
c
c     make the call to the optimization routine
c
      iprint=0
      call ocvm (nvar,xx,minimum,grdmin,optimiz1,optsave)
c
c     unscale the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     compute the final function and RMS gradient values
c
      if (analytic) then
         call gradient (minimum,derivs)
      else
         minimum = energy ()
         call numgrad (energy,derivs,eps)
      end if
      gnorm = 0.0d0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               gnorm = gnorm + derivs(j,i)**2
            end do
         end if
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(nvar/3))
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (derivs)
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         xyzCoord(3*(i-1)+1) = x(i)
         xyzCoord(3*(i-1)+2) = y(i)
         xyzCoord(3*(i-1)+3) = z(i)
      end do
c
c     write out the final function and gradient values
c
c      if (digits .ge. 8) then
c         if (grms .gt. 1.0d-8) then
c            write (iout,60)  minimum,grms,gnorm
c   60       format (/,' Final Function Value :',2x,f20.8,
c     &              /,' Final RMS Gradient :',4x,f20.8,
c     &              /,' Final Gradient Norm :',3x,f20.8)
c         else
c            write (iout,70)  minimum,grms,gnorm
c   70       format (/,' Final Function Value :',2x,f20.8,
c     &              /,' Final RMS Gradient :',4x,d20.8,
c     &              /,' Final Gradient Norm :',3x,d20.8)
c         end if
c      else if (digits .ge. 6) then
c         if (grms .gt. 1.0d-6) then
c            write (iout,80)  minimum,grms,gnorm
c   80       format (/,' Final Function Value :',2x,f18.6,
c     &              /,' Final RMS Gradient :',4x,f18.6,
c     &              /,' Final Gradient Norm :',3x,f18.6)
c         else
c            write (iout,90)  minimum,grms,gnorm
c   90       format (/,' Final Function Value :',2x,f18.6,
c     &              /,' Final RMS Gradient :',4x,d18.6,
c     &              /,' Final Gradient Norm :',3x,d18.6)
c         end if
c      else
c         if (grms .gt. 1.0d-4) then
c            write (iout,100)  minimum,grms,gnorm
c  100       format (/,' Final Function Value :',2x,f16.4,
c     &              /,' Final RMS Gradient :',4x,f16.4,
c     &              /,' Final Gradient Norm :',3x,f16.4)
c         else
c            write (iout,110)  minimum,grms,gnorm
c  110       format (/,' Final Function Value :',2x,f16.4,
c     &              /,' Final RMS Gradient :',4x,d16.4,
c     &              /,' Final Gradient Norm :',3x,d16.4)
c         end if
c      end if
c
c     write the final coordinates into a file
c
c      imin = freeunit ()
c      open (unit=imin,file=minfile,status='old')
c      rewind (unit=imin)
c      call prtxyz (imin)
c      close (unit=imin)
c
c     perform any final tasks before program exit
c
c      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function optimiz1  --  energy and gradient for optimize  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "optimiz1" is a service routine that computes the energy and
c     gradient for optimally conditioned variable metric optimization
c     in Cartesian coordinate space
c
c
      function optimiz1 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'scales.i'
      include 'usage.i'
      integer i,nvar
      real*8 optimiz1,e
      real*8 energy,eps
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
      logical analytic
      external energy
c
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      if (analytic) then
         call gradient (e,derivs)
      else
         e = energy ()
         call numgrad (energy,derivs,eps)
      end if
      optimiz1 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(2,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(3,i) / scale(nvar)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
cc
cc  DES: from liblas/prntmat.F
cc
c      SUBROUTINE PRNTMAT(N,NROW,NCOL,A)
c      REAL*8 A(NROW,NCOL)
c      PARAMETER (MAXCOL=6)
c      NP = N
c      IF(NP.GT.NCOL) NP = NCOL  ! CAN'T PRINT MORE THAN NCOL
c      NT = NP/MAXCOL
c      IF(NT.EQ.0) GO TO 30
c      DO 20 I=1,NT
c         IMIN = (I-1)*MAXCOL + 1
c         IMAX = I*MAXCOL
c         WRITE(6,1000)
c         DO 10 J=1,NROW
c            WRITE(6,1100) (A(J,K),K=IMIN,IMAX)
c 10      CONTINUE
c 20   CONTINUE
c 30   CONTINUE
c      NS = NT*MAXCOL
c      NLEFT = NP - NS
c      IF(NLEFT.EQ.0) RETURN
c      WRITE(6,1000)
c      DO 40 J=1,NROW
c         WRITE(6,1100) (A(J,K),K=NS+1,NP)
c 40   CONTINUE
c      RETURN
c 1000 FORMAT(/)
c 1100 FORMAT(1X,6F12.6)
c      END

