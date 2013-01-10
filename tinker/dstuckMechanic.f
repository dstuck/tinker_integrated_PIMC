
      subroutine dstuckmechanic(prmName, prmLen)
      implicit none
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'sizes.i'
      include 'files.i'
      include 'keys.i'
      include 'params.i'
      include 'bound.i'
      include 'boxes.i'

      character*120 prmName
      integer prmLen
      logical useprm
      integer iprm
      character*120 prmfile
      character*120 record
      integer freeunit
      
      write(*,*) "in mechanic"

c DES: from mechanic.f
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      write(*,*) "after attach"
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      write(*,*) "before bonds"
      call bonds
      call angles
      call torsions
      call bitors
      call rings
      write(*,*) "after rings"
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
      write(*,*) "before getprm.f"
c  DES: from getprm.f
      useprm = .true.
c      prmfile = filename(1:leng)//'.prm'
      prmName = '../sulfate_water.prm'
      prmLen = 20
      prmfile(1:prmLen) = prmName
      write(*,*) "before initprm"
      call initprm
      write(*,*) "after initprm"
      nprm = 0
      if (useprm) then
         iprm = freeunit ()
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
      write(*,*) "before readprm"
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
      write(*,*) "leaving mechanic"
c DES out
      return
      end

