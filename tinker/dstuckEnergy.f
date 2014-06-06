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
      subroutine dstuckenergy(numAt, typeVec, xyzCoord, connMat,
     $prmName, prmLen, eOut, initialize)
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
      include 'align.i'
      include 'atoms.i'
      include 'bath.i'
      include 'cell.i'
      include 'files.i'
      include 'group.i'
      include 'inform.i'
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
      include 'scales.i'
      include 'sequen.i'
      include 'socket.i'
      include 'warp.i'
      include 'zclose.i'
c DES out
c DES: from mechanic.f misc calls
      include 'boxes.i'
c DES out
c DES: from readxyz.f
      include 'atmtyp.i'
      include 'couple.i'
      include 'titles.i'
c DES out
      integer numAt
      integer, dimension(numAt) :: typeVec
      real*8, dimension(numAt*3) :: xyzCoord
      real*8, dimension(numAt*8) :: connMat
      character*120 prmName
      integer prmLen
      real*8 eOut
c      logical initialize
      integer initialize
      real*8 energy
      real*8 cutoff
c DES: from initial.f
      real*8 precise
c DES out
c DES: from getxyz.f
      integer freeunit
      character*120 xyzfile
c DES out
c DES: from getprm.f
      logical useprm
      integer iprm
c      character*120 prmfile
      character(len=prmLen) prmfile
      character*120 record
c DES out 
c DES: from readxyz.f
      integer i,j,k,m
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
c DES: readxyz -> replaced by directly passing info
c
c     initialize the total number of atoms in the system
c
      n = numAt
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

      if(initialize.NE.0) then
c DES: from mechanic.f
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     find unit cell type, lattice parameters and cutoff values
c
c      call unitcell
      use_bounds = .false.
      use_replica = .false.
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
      orthogonal = .false.
      monoclinic = .false.
      triclinic = .false.
      octahedron = .false.
      spacegrp = '          '
c      nosymm = .false.
c      call lattice
c      call polymer
      use_polymer = .false.
      polycut = 5.5d0
      call cutoffs
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     get the force field parameters and assign atom types
c
c      call field
      use_bond = .true.
      use_angle = .true.
      use_strbnd = .true.
      use_urey = .true.
      use_angang = .true.
      use_opbend = .true.
      use_opdist = .true.
      use_improp = .true.
      use_imptor = .true.
      use_tors = .true.
      use_pitors = .true.
      use_strtor = .true.
      use_tortor = .true.
      use_vdw = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
      use_mpole = .true.
      use_polar = .true.
      use_rxnfld = .false.
      use_solv = .true.
      use_metal = .false.
      use_geom = .true.
      use_extra = .true.
c  DES: from getprm.f
      useprm = .true.
c      prmfile = filename(1:leng)//'.prm'
c      prmName = '../sulfate_water.prm'
c      prmLen = 20
      prmfile(1:prmLen) = prmName
      call initprm
      nprm = 0
      if (useprm) then
         iprm = freeunit ()
c         write (*,*) "prmfile = ", prmfile
         open (unit=iprm,file=prmfile,status='old')
         rewind (unit=iprm)
         do while (.true.)
            read (iprm,30,err=50,end=50)  record
   30       format (a120)
            nprm = nprm + 1
            prmline(nprm) = record
            if (nprm .ge. maxprm) then
               write (iout,40)
   40          format (/,' GETPRM  --  Parameter File Too Large;',
     &                    ' Increase MAXPRM')
               call fatal
            end if
         end do
   50    continue
         close (unit=iprm)
      end if
c
c     get control and parameter values from the parameter file
c
      if (useprm)  call readprm
c DES out
c
      call katom
c
c     assign atoms to molcules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
c
c     assign bond, angle and cross term potential parameters
c
      if (use_bond .or. use_strbnd .or. use_strtor
     &    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
      if (use_angle .or. use_strbnd .or. use_angang)  call kangle
      if (use_strbnd)  call kstrbnd
      if (use_urey)  call kurey
      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
      if (use_angle .or. use_opbend)  call kopbend
      if (use_angle .or. use_opdist)  call kopdist
      if (use_improp)  call kimprop
      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
      if (use_pitors)  call kpitors
      if (use_strtor)  call kstrtor
      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole
      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld)  call kmpole
      if (use_polar .or. use_solv)  call kpolar
      if (use_ewald)  call kewald
c
c     assign solvation, metal, pisystem and restraint parameters
c
      if (use_solv)  call ksolv
      if (use_metal)  call kmetal
      if (use_orbit)  call korbit
      if (use_geom)  call kgeom
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,11)
   11    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
c DES out
c DES: End of initialize if block
      else

c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond
      if (use_angle)  call eangle
      if (use_strbnd)  call estrbnd
      if (use_urey)  call eurey
      if (use_angang)  call eangang
      if (use_opbend)  call eopbend
      if (use_opdist)  call eopdist
      if (use_improp)  call eimprop
      if (use_imptor)  call eimptor
      if (use_tors)  call etors
      if (use_pitors)  call epitors
      if (use_strtor)  call estrtor
      if (use_tortor)  call etortor
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge)  call echarge
      if (use_chgdpl)  call echgdpl
      if (use_dipole)  call edipole
      if (use_mpole .or. use_polar)  call empole
      if (use_rxnfld)  call erxnfld
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv
      if (use_geom)  call egeom
      if (use_metal)  call emetal
      if (use_extra)  call extra
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
      energy = esum
      eOut = energy
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' ENERGY  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      end if
      return
      end
