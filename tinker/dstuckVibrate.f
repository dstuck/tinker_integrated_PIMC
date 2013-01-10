c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  function energy  --  evaluates energy terms and total  ##
c     ##  subroutine dstuckenergy  --  evaluates energy total    ##
c     ##                                                         ##
c     #############################################################
c
c
c     "energy" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy
c
c
      subroutine dstuckvibrate(typeVec, xyzCoord, connMat, modeNum, 
     $partNum, redMass, evectors, evalues)
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
c DES: from vibrate
      include 'atmtyp.i'
      include 'hescut.i'
      include 'math.i'
      include 'units.i'
      include 'usage.i'
c DES out
c DES: from mechanic.f misc calls
c      include 'boxes.i'
c DES out
c DES: from readxyz.f
      include 'couple.i'
      include 'titles.i'
c DES out
      integer modeNum
      integer partNum
      integer, dimension(partNum) :: typeVec
      real*8, dimension(partNum*3) :: xyzCoord
      real*8, dimension(partNum*8) :: connMat
      real*8, dimension(modeNum) :: redMass
      real*8  tempInvMass
      real*8, dimension(modeNum*partnum*3) :: evectors
      real*8, dimension(modeNum) :: evalues
c DES: from initial.f
      real*8 precise
c DES out
c DES: from getxyz.f
      integer ixyz
      integer freeunit
      logical exist
      character*120 xyzfile
c DES out
c DES: from vibrate
      integer i,j,k,m
      integer ihess
      integer lext
      integer nfreq,ndummy
      integer nvib,ivib
      integer nview,next
      integer nlist,ilist
      integer, allocatable :: list(:)
      integer, allocatable :: iv(:)
      integer, allocatable :: hindex(:)
      integer, allocatable :: hinit(:,:)
      integer, allocatable :: hstop(:,:)
      real*8 factor,vnorm
      real*8 sum,ratio
c      real*8 scale
      real*8, allocatable :: xref(:)
      real*8, allocatable :: yref(:)
      real*8, allocatable :: zref(:)
      real*8, allocatable :: mass2(:)
      real*8, allocatable :: h(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: matrix(:)
      real*8, allocatable :: hdiag(:,:)
      real*8, allocatable :: vects(:,:)
      logical query
      character*1 letter
      character*7 ext
      character*120 record
      character*120 string
c DES out
c DES: for projtrm
      integer hcounter
      real*8 xc(3,N)
      real*8 hfull(3*N,3*N)
      real*8 trvec(3*N,3*N)
      real*8 temp(3*N,3*N)
c DES out
c DES: from readxyz.f
c      integer i,j,k,m
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

c DES from vibrate.f
c
c     perform dynamic allocation of some local arrays
c
      nfreq = 3 * nuse
      allocate (mass2(n))
      allocate (hinit(3,n))
      allocate (hstop(3,n))
      allocate (hdiag(3,n))
      allocate (hindex(nfreq*(nfreq-1)/2))
      allocate (h(nfreq*(nfreq-1)/2))
      allocate (matrix(nfreq*(nfreq+1)/2))
c
c     initialize various things needed for vibrations
c
      ndummy = 0
      do i = 1, n
         if (use(i) .and. atomic(i).eq.0) then
            ndummy = ndummy + 1
            mass(i) = 0.001d0
         end if
         mass2(i) = sqrt(mass(i))
      end do
      nvib = nfreq - 3*ndummy
c
c     calculate the Hessian matrix of second derivatives
c
      hesscut = 0.0d0
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     store upper triangle of the Hessian in "matrix"
c
      ihess = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               ihess = ihess + 1
               matrix(ihess) = hdiag(j,i)
               do k = hinit(j,i), hstop(j,i)
                  m = (hindex(k)+2) / 3
                  if (use(m)) then
                     ihess = ihess + 1
                     matrix(ihess) = h(k)
                  end if
               end do
            end do
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (eigen(nfreq))
      allocate (vects(nfreq,nfreq))
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nfreq,nfreq,matrix,eigen,vects)
c      write (iout,10)
c   10 format (/,' Eigenvalues of the Hessian Matrix :',/)
c      write (iout,20)  (i,eigen(i),i=1,nvib)
c   20 format (5(i5,f10.3))
c
c     store upper triangle of the mass-weighted Hessian matrix
c
      ihess = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               ihess = ihess + 1
               matrix(ihess) = hdiag(j,i) / mass(i)
               do k = hinit(j,i), hstop(j,i)
                  m = (hindex(k)+2) / 3
                  if (use(m)) then
                     ihess = ihess + 1
                     matrix(ihess) = h(k) / (mass2(i)*mass2(m))
                  end if
               end do
            end do
         end if
      end do
c
c  DES: Vectorize the way projtrm wants its inputs
c
      do i=1, N
         xc(1,i) = x(i)
         xc(2,i) = y(i)
         xc(3,i) = z(i)
      end do
c
c  DES: Fill out hessian from upper triangle
c
c      hcounter = N*(N+1)/2;
c      do i=3*N, 1, -1
c         hfull(i,i) = matrix(hcounter)
c         hcounter = hcounter-1
c         do j=i-1, 1, -1
cc TODO: Is this backwards?
c            hfull(j,i) = matrix(hcounter)
c            hfull(i,j) = hfull(j,i)
c            hcounter = hcounter-1
c         end do
c      end do
c
      hcounter = 1;
      do i=1, 3*N 
         hfull(i,i) = matrix(hcounter)
         hcounter = hcounter+1
         do j=i+1, 3*N
            hfull(j,i) = matrix(hcounter)
            hfull(i,j) = hfull(j,i)
            hcounter = hcounter+1
         end do
      end do
c  DES: Project rotation and translation from the Hessian
c
c      write(*,*) "hfull before projection"
c      call prntmat(9,9,9,hfull)
      call projtrm(n,1,xc,mass,temp,trvec,hfull)
c      write(*,*) "hfull"
c      call prntmat(9,9,9,hfull)
c
c  DES: Fill out upper triangle from hfull
c
c      hcounter = N*(N+1)/2;
c      do i=3*N, 1, -1
c         matrix(hcounter) = hfull(i,i)
c         hcounter = hcounter-1
c         do j=i-1, 1, -1
cc TODO: Is this backwards?
c            matrix(hcounter) = hfull(j,i)
c            hcounter = hcounter-1
c         end do
c      end do
      hcounter = 1;
      do i=1, 3*N 
         matrix(hcounter) = hfull(i,i)
         hcounter = hcounter+1
         do j=i+1, 3*N
            matrix(hcounter) = hfull(j,i)
            hcounter = hcounter+1
         end do
      end do
c
c     diagonalize to get vibrational frequencies and normal modes
c
      call diagq (nfreq,nfreq,matrix,eigen,vects)
c      write(*,*) "vects"
c      call prntmat(9,9,9,vects)
c      write(*,*) "eigen"
c      call prntmat(6,9,1,eigen)
      factor = sqrt(convert) / (2.0d0*pi*lightspd)
      do i = 1, nvib
         eigen(i) = factor * sign(1.0d0,eigen(i)) * sqrt(abs(eigen(i)))
      end do
c      write (iout,30)
c   30 format (/,' Vibrational Frequencies (cm-1) :',/)
c      write (iout,40)  (i,eigen(i),i=1,nvib)
c   40 format (5(i5,f10.3))
c  DES: Pass frequencies in cm-1 back to V_Tinker
      do i = 1, modenum
         evalues(modenum+1-i) = eigen(nfreq+1-i)
      end do
c  DES: Form reduced mass and pass back to V_Tinker
      do i = 1, modenum
         tempInvMass = 0.0
         do j = 1, partnum
            do k = 1, 3
c              tempInvMass = tempInvMass + evectors(k+3*(j-1)+
c     $3*partnum*(i-1))**2/mass(j)
              tempInvMass = tempInvMass + vects(k+3*(j-1),i+
     $(nfreq-modenum))**2/mass(j)
            end do
         end do
         redMass(i) = tempInvMass**(-1)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (hinit)
      deallocate (hstop)
      deallocate (hdiag)
      deallocate (h)
      deallocate (matrix)
c
c     form Cartesian coordinate displacements from normal modes
c
      do i = 1, nvib
         vnorm = 0.0d0
         do j = 1, nfreq
            k = iuse((j+2)/3)
            vects(j,i) = vects(j,i) / mass2(k)
            vnorm = vnorm + vects(j,i)**2
         end do
         vnorm = sqrt(vnorm)
         do j = 1, nfreq
            vects(j,i) = vects(j,i) / vnorm
         end do
      end do
c  DES: Pass normalModes in cm-1 back to V_Tinker
      do i = 1, modenum
         do j = 1, partnum
            do k = 1, 3
               evectors(k + 3*(j-1) + 3*partnum*(i-1)) = vects(k+3*(j-1)
     $,i+(nfreq-modenum))
            end do
         end do
      end do

c
c     perform dynamic allocation of some local arrays
c
      deallocate (mass2)
      deallocate (eigen)
      deallocate (vects)
c  DES note: may need to call final.f??
c  DES out

      return
      end
c
c  DES: Snagged from qchem! Modified to fix tinker
c
c       $Id: projtrm.F,v 1.1 1997/06/20 20:41:37 johnson Rel $
c
      SUBROUTINE PROJTRM(NATOMS, IPRNT,  XC,     ATMASS,  P,
     $                   TRVEC,  HESS)                      
      IMPLICIT REAL*8(A-H,O-Z)                              
      REAL*8 XC(3,NATOMS),ATMASS(NATOMS),P(3*NATOMS,*),     
     $       TRVEC(3*NATOMS,*),HESS(3*NATOMS,3*NATOMS),
     $       eVecMat(3,3)
      REAL*8 T(9),PMOM(3),TUpper(6)
      PARAMETER (ONE=1.0D0)                                 
      CALL COM(NATOMS,ATMASS,XC,CX,CY,CZ)                   
      CALL DZERO(9,T)
      DO 20 I=1,NATOMS                                      
         X = XC(1,I)
         Y = XC(2,I)
         Z = XC(3,I)
         AMI = ATMASS(I)                                       
         T(1) = T(1) + AMI*(Y*Y + Z*Z)                         
         T(5) = T(5) + AMI*(X*X + Z*Z)                         
         T(9) = T(9) + AMI*(X*X + Y*Y)                         
         T(2) = T(2) - AMI*X*Y                                 
         T(3) = T(3) - AMI*X*Z                                 
         T(6) = T(6) - AMI*Y*Z                                 
 20   CONTINUE 
      T(4) = T(2)
      T(7) = T(3)
      T(8) = T(6)
c      write(*,*) "Printing T"
c      call prntmat(3,3,3,T)
c DES: changing DIAGMAT call
c      CALL DIAGMAT(T,3,TRVEC,P,PMOM,IERR)    ! TRVEC & P USED AS SCRATCH
      do i=1, 6
         iOff = 0
         if(i.gt.3) then
            iOff = iOff+1
         end if
         if(i.gt.5) then
            iOff = iOff+2
         end if
         TUpper(i) = T(i+iOff)
      end do
      call diagq (3,3,TUpper,PMOM,eVecMat)
c      write(*,*) "Printing PMOM"
c      call prntmat(3,3,1,PMOM)
c      write(*,*) "Printing eVecMat"
c      call prntmat(3,3,3,eVecMat)
      do i=1, 3
         do j=1, 3
c  TODO: Is this reversed?
            T(j+(i-1)*3) = eVecMat(j,i)
         end do
      end do

      IF(IERR.NE.0) THEN                         
       WRITE(6,1000)                             
       CALL OPTEXIT(9)                           
      ENDIF                                      
      NAT3 = 3*NATOMS                            
      CALL DZERO(NAT3*6,TRVEC)                   
      CALL FORMTRM(NATOMS,ATMASS,XC,T,TRVEC)     
      IF(IPRNT.GT.5) THEN                        
       WRITE(6,1100)                             
       CALL PRNTMAT(6,NAT3,6,TRVEC)              
      ENDIF                                      
      CALL DZERO(NAT3*NAT3,P)                    
      DO 30 K=1,6                                
      DO 30 J=1,NAT3                             
      DO 30 I=1,NAT3                             
      P(I,J) = P(I,J) - TRVEC(I,K)*TRVEC(J,K)    
 30   CONTINUE                                   
      DO 40 I=1,NAT3                             
      P(I,I) = ONE + P(I,I)                      
 40   CONTINUE                                   
      IF(IPRNT.GT.5) THEN                        
       WRITE(6,1200)                             
       CALL PRNTMAT(NAT3,NAT3,NAT3,P)            
      ENDIF
      CALL MATAB(NAT3,NAT3,NAT3,HESS,P,TRVEC,0)  
      CALL MATAB(NAT3,NAT3,NAT3,P,TRVEC,HESS,0)  
      IF(IPRNT.GT.1) WRITE(6,1300)               
      DO 50 I=1,NATOMS
      XC(1,I) = XC(1,I) + CX
      XC(2,I) = XC(2,I) + CY
      XC(3,I) = XC(3,I) + CZ
 50   CONTINUE                                   
c      write(*,*) "Hess in Proj"
c      call prntmat(9,9,9,hess)
      RETURN                                     
 1000 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize inertia tensor') 
 1100 FORMAT(/,' Vectors for Translations and Rotations')             
 1200 FORMAT(/,' The Projection Matrix')                              
 1300 FORMAT(/,' Translations and Rotations Projected Out of Hessian')
      END                            
*Deck formtrm                                                                              
      SUBROUTINE FormTRM(NATOMS,ATMASS,XC,T,V)  
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ATMASS(NATOMS),XC(3,NATOMS),T(9),V(3,NATOMS,6)
      PARAMETER (ONE=1.0d0,TollZERO=1.0d-8)
      NAT3 = 3*NATOMS
      DO 10 I=1,NATOMS
         X = XC(1,I)
         Y = XC(2,I)
         Z = XC(3,I)
         ami = DSQRT(ATMASS(I))
         CX = ami*(X*T(1) + Y*T(2) + Z*T(3))
         CY = ami*(X*T(4) + Y*T(5) + Z*T(6))
         CZ = ami*(X*T(7) + Y*T(8) + Z*T(9))
         V(1,I,1) = ami*ONE
         V(2,I,2) = ami*ONE
         V(3,I,3) = ami*ONE
         V(1,I,4) = CY*T(7) - CZ*T(4)
         V(2,I,4) = CY*T(8) - CZ*T(5)
         V(3,I,4) = CY*T(9) - CZ*T(6)
         V(1,I,5) = CZ*T(1) - CX*T(7)
         V(2,I,5) = CZ*T(2) - CX*T(8)
         V(3,I,5) = CZ*T(3) - CX*T(9)
         V(1,I,6) = CX*T(4) - CY*T(1)
         V(2,I,6) = CX*T(5) - CY*T(2)
         V(3,I,6) = CX*T(6) - CY*T(3)
 10   CONTINUE
      DO 20 I=1,6
         skal = SProd(NAT3,V(1,1,I),V(1,1,I))
         IF(skal.GT.TollZERO) THEN
            skal = ONE/SQRT(skal)
            CALL VScal(NAT3,skal,V(1,1,I))
         ENDIF
 20   CONTINUE
      RETURN
      END 
c       $Id: com.F,v 1.1 1997/06/20 20:40:53 johnson Rel $                                                                                                                             
c
c  DES: vibman.com
c
      SUBROUTINE COM(NATOMS,ATMASS,XC,X,Y,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ATMASS(NATOMS),XC(3,NATOMS)
      PARAMETER (ZERO=0.0D0)
      TOTMAS = ZERO
      X = ZERO
      Y = ZERO
      Z = ZERO
      DO 10 IATM=1,NATOMS
      AMI = ATMASS(IATM)
      TOTMAS = TOTMAS + AMI
      X = X + AMI*XC(1,IATM)
      Y = Y + AMI*XC(2,IATM)
      Z = Z + AMI*XC(3,IATM)
 10   CONTINUE
      X = X/TOTMAS
      Y = Y/TOTMAS
      Z = Z/TOTMAS
      DO 20 IATM=1,NATOMS
      XC(1,IATM) = XC(1,IATM) - X
      XC(2,IATM) = XC(2,IATM) - Y
      XC(3,IATM) = XC(3,IATM) - Z
 20   CONTINUE
      RETURN
      END
c
c  DES: from libopt/dzero.F
c
      SUBROUTINE DZERO(N,A)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(N)
      PARAMETER (ZERO=0.0D0)
      DO 10 I=1,N
      A(I) = ZERO
 10   CONTINUE
      RETURN
      END
c
c  DES: from libopt/vscal.F
c
      SUBROUTINE VSCAL(N,SKAL,V)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 V(N)
      DO 10 I=1,N
      V(I) = V(I)*SKAL
 10   CONTINUE
      RETURN
      END 
c
c  DES: from libopt/matab.F
c
      SUBROUTINE MATAB(NR,NC,NI,A,B,C,IACC)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NR,NI),B(NI,NC),C(NR,NC)
      IF(IACC.EQ.0) CALL DZERO(NR*NC,C)
      DO 20 J=1,NC
      DO 20 K=1,NI
      BVAL = B(K,J)
      DO 10 I=1,NR
      C(I,J) = C(I,J) + A(I,K)*BVAL
 10   CONTINUE
 20   CONTINUE
      RETURN
      END
c
c  DES: from libopt/sprod.F
c
      DOUBLE PRECISION FUNCTION SPROD(N,V1,V2)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 V1(N),V2(N)
      PARAMETER (ZERO=0.0D0)
      SPROD = ZERO
      DO 10 I=1,N
      SPROD = SPROD + V1(I)*V2(I)
 10   CONTINUE
      RETURN
      END
c
c  DES: from liblas/prntmat.F
c
      SUBROUTINE PRNTMAT(N,NROW,NCOL,A)
      REAL*8 A(NROW,NCOL)
      PARAMETER (MAXCOL=6)
      NP = N
      IF(NP.GT.NCOL) NP = NCOL  ! CAN'T PRINT MORE THAN NCOL
      NT = NP/MAXCOL
      IF(NT.EQ.0) GO TO 30
      DO 20 I=1,NT
         IMIN = (I-1)*MAXCOL + 1
         IMAX = I*MAXCOL
         WRITE(6,1000)
         DO 10 J=1,NROW
            WRITE(6,1100) (A(J,K),K=IMIN,IMAX)
 10      CONTINUE
 20   CONTINUE
 30   CONTINUE
      NS = NT*MAXCOL
      NLEFT = NP - NS
      IF(NLEFT.EQ.0) RETURN
      WRITE(6,1000)
      DO 40 J=1,NROW
         WRITE(6,1100) (A(J,K),K=NS+1,NP)
 40   CONTINUE
      RETURN
 1000 FORMAT(/)
 1100 FORMAT(1X,6F12.6)
      END

