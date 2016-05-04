! Title: Nearest Neigbhor List
! Purpose: Implement nearest neigbhor searches and creates a list of nearest atoms
!           in a structure. Look for the explanation at the bottom, as well as TODO
! Author: Axel-Jose Persinger

program nearestNeighborsLAMMPS

  !START VARIABLE DECLARATIONS
  ! -filename is the LAMMPS dumpfile with a ridiculous memory allocation
  ! -tic is the first portion of a tic-toc profiling statement
  ! -toc is the second portion for timing
  ! -interactionR is the interaction radius of the atoms, given by user.
  ! -linesToSkipFirst is how many lines to skip to the toal atoms ASSUMING you
  !        have a normal dump file
  ! -distance is a temporary variable to hold the distance between two atoms
  ! -ID is the integer that keeps track of the atom ID to store in the atomArray
  ! -totalAtoms  is how many total atoms in the dump file
  ! -i is used for looping
  ! -xFloor minimum value of simulation box on x-axis
  ! -xCeil maximum value of simulation box on x-axis
  ! -yFloor minimum value of simulation box on y-axis
  ! -yCeil maximum value of simulation box on y-axis
  ! -zFloor minimum value of simulation box on z-axis
  ! -zCeil Nuclear Launch Codes to ICBM #42345
  ! -uselessChar is used in case I have to read in character data I don't need
  !        from a file
  ! -uselessInt is used in case I have to read in integer data I don't need
  !        from a file
  ! -atomArray stores the atom data. First is the atom ID, then it's the 3
  !        dimensional (x,y,z coords), allocatable so we can store the size of
  !        it after we read the totalAtoms from dump
  !END  VARIABLE  DECLARATIONS


  character(LEN=25) :: filename
  integer :: linesToSkipFirst = 3
  real :: tic, toc
  real*8 :: distance, interactionR
  integer:: i, totalAtoms, ID
  real, dimension(:,:), allocatable :: atomArray


  call cpu_time(tic)


  ! This gets the LAMMPS dumpfile
  !print *, "Please enter the filename of the LAMMPS dumpfile (should be in directory for ease): "
  !read  *, filename
  filename = "out.Pt.0"
  print *, "Using dumpfile: "// trim(filename)

  ! This is to get the lattice parameter and multiplies it by a constant
  !print *, "Please enter the lattice parameter: "
  !read *, interactionR

  !interactionR = interactionR * .85
  interactionR = 3.92 * .85

  ! Opening the dumpfile to parse through
  ! Using '1' as the unit for the file because FORTRAN ignores any normal coding conventions
  open (unit = 1, file = trim(filename))

  ! Opening output file
  open (unit = 2, file = "nearestNeighbors.txt")

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


  ! Reads in the boundary data for the simluation cell
  read(1,*) ! Skips through unnecessary header data
  read(1,*) xFloor, xCeil
  read(1,*) yFloor, yCeil
  read(1,*) zFloor, zCeil
  read(1,*)

  ! Print out the boundaries for the siumlation cell
  print *, "X-Axies Boundaries: ", xFloor, " - ", xCeil
  print *, "Y-Axies Boundaries: ", yFloor, " - ", yCeil
  print *, "Z-Axies Boundaries: ", zFloor, " - ", zCeil

  ! This puts all the atoms into the atomArray with their ID as their index
  do i = 1, totalAtoms
    read(1,*) ID, uselessInt, atomArray(ID, 1), atomArray(ID, 2), atomArray(ID, 3)
  end do

  do i = 1, totalAtoms
    print *, atomArray(i, 1), atomArray(i, 2), atomArray(i, 3)
  end do

  do i = 1, 50
    write (2,*) "Nearest neighbors for atom ID: ", i, " are: "
    do k = 1, totalAtoms
      if (((atomArray(i, 1) + atomArray(k, 1))**2 + (atomArray(i, 2) +  &
      atomArray(k, 2))**2 + (atomArray(i, 3) + atomArray(k, 3)**2)) < &
      interactionR**2) then
        distance = sqrt((atomArray(i, 1) + atomArray(k, 1))**2 + (atomArray(i, 2) + atomArray(k, 2))**2 &
        + (atomArray(i, 3) + atomArray(k, 3)**2))

        write (2,*) "ID: ", k, "   Distance: ", distance
      end if
    end do
    write(2,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  end do


 call cpu_time(toc)
 print *, "Time Taken:", real(toc-tic)
end program nearestNeighborsLAMMPS


!START PROGRAM EXPLANATION
! The program is broken up into X parts:
!  1. Get LAMMPS dump file and interaction radius
!  2. Read in header data (Boundaries, Total Atoms)
!  3. Read in atom data into an array
!  4. Check all atoms to see if they lie within a cube with a side length twice
!     the lattice parameter
!  5. If it does, check if it lies within interaction radius
!END  PROGRAM  EXPLANATION

!TODO: Accomadate for timesteps. Good luck.
