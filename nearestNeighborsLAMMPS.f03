! Title: Nearest Neigbhor List
! Purpose: Implement nearest neigbhor searches and creates a list of nearest atoms
!           in a structure. Look for the explanation at the bottom, as well as TODO
! Author: Axel-Jose Persinger

program nearestNeighborsLAMMPS

  !START VARIABLE DECLARATIONS
  ! -filename is the LAMMPS dumpfile with a ridiculous memory allocation
  ! -interactionR is the interaction radius of the atoms, given by user.
  ! -linesToSkipFirst is how many lines to skip to the toal atoms ASSUMING you
  !        have a normal dump file
  ! -ID is the integer that keeps track of the atom ID to store in the atomArray
  ! -totalAtoms  is how many total atoms in the dump file
  ! -xBoundCiel is the upper ammount any number can be on x-axis (important for
  !        memory allocation)
  ! -yBoundCiel is the upper ammount any number can be on y-axis (important for
  !        memory allocation)
  ! -zBoundCiel is the upper ammount any number can be on z-axis (important for
  !        memory allocation)
  ! -xBoundFloor is the lower ammount any number can be on x-axis (important for
  !        memory allocation)
  ! -yBoundFloor is the lower ammount any number can be on y-axis (important for
  !        memory allocation)
  ! -zBoundFloor is the lower ammount any number can be on z-axis (important for
  !        memory allocation)
  ! -i is used for looping
  !-uselessChar is used in case I have to read in character data I don't need
  !        from a file
  !-uselessInt is used in case I have to read in integer data I don't need
  !        from a file
  !-atomArray stores the atom data.First is the atom ID, then it's the 3
  !        dimensional (x,y,z coords) allocatable so we can store the size of
  !        it after we read the totalAtoms from dump
  !END  VARIABLE  DECLARATIONS


  character(LEN=1000) :: filename
  integer :: linesToSkipFirst = 3
  real :: xBoundCiel, yBoundCiel, zBoundCiel, xBoundFloor, yBoundFloor, zBoundFloor
  integer:: i, totalAtoms, ID
  real, dimension(:,:), allocatable :: atomArray
  ! This gets the LAMMPS dumpfile
  print *, "Please enter the filename of the LAMMPS dumpfile (should be in directory for ease): "
  !read  *, filename
  filename = "out.Pt_3nm.20000"
  print *, "Using dumpfile: "// trim(filename)

  print *, "Please enter the lattice parameter: "
  read *, interactionR
  interactionR = interactionR * .85

  ! Opening the dumpfile to parse through
  ! Using '1' as the unit for the file because FORTRAN ignores any normal coding conventions
  open (unit = 1, file = trim(filename))

  ! Read through unnecessary header data until we get to the total atoms
  do i = 1, linesToSkipFirst
    read(1,*)
  end do

  ! Read the total number of atoms and convert the atomArray to appropriate size
  read(1,*) totalAtoms
  print *, "Total number of atoms: "
  totalAtoms = int(totalAtoms)
  print *, totalAtoms
  allocate(atomArray(totalAtoms, 3))

  ! Gets the boundaries of the simulation box. Okay what now?
  read(1,*) ! Skips through unnecessary header data
  read(1,*) xBoundFloor, xBoundCiel
  read(1,*) yBoundFloor, yBoundCiel
  read(1,*) zBoundFloor, zBoundCiel

  ! Debugging purposes, just printing the information... Ya boy still got it
  print *, xBoundFloor
  print *, xBoundCiel
  print *, yBoundFloor
  print *, yBoundCiel
  print *, zBoundFloor
  print *, zBoundCiel

  read(1,*) ! Skips through unnecessary header data

  ! LET'S GET THE ATOMS, HOMEBOY
  do i = 1, totalAtoms, 1
    read(1,*) ID, uselessInt, atomArray (ID, 1), atomArray (ID, 2), atomArray (ID, 3)
  end do



end program nearestNeighborsLAMMPS




!START PROGRAM EXPLANATION
! The program is broken up into X parts:
!  1. Get LAMMPS dump file and interaction radius
!  2. Read in header data (Boundaries, Total Atoms)
!  3. Read in atom data into an array
!  4. Calculate sphereical coordinates of interactions
!END  PROGRAM  EXPLANATION

!TODO: Accomadate for timesteps. Good luck.
