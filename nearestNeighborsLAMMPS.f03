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
  ! -nearestIterable is an interator to count up the array for the atoms
  ! -linesToSkipFirst is how many lines to skip to the toal atoms ASSUMING you
  !        have a normal dump file
  ! -ID is the integer that keeps track of the atom ID to store in the atomArray
  ! -totalAtoms  is how many total atoms in the dump file
  ! -i is used for looping
  ! -uselessChar is used in case I have to read in character data I don't need
  !        from a file
  ! -uselessInt is used in case I have to read in integer data I don't need
  !        from a file
  ! -atomArray stores the atom data. First is the atom ID, then it's the 3
  !        dimensional (x,y,z coords) allocatable so we can store the size of
  !        it after we read the totalAtoms from dump
  ! -tempArray stores the nearest neighbors
  ! -nearestNeighbors is the final array that stores all the nearest neighbors
  !END  VARIABLE  DECLARATIONS


  character(LEN=1000) :: filename
  integer :: linesToSkipFirst = 3
  integer :: nearestIterable = 1
  real :: tic, toc
  integer:: i, totalAtoms, ID
  real, dimension(:,:), allocatable :: atomArray
  integer, dimension(:), allocatable :: tempArray
  integer, dimension(:), allocatable :: nearestNeighbors


  call cpu_time(tic)


  ! This gets the LAMMPS dumpfile
  print *, "Please enter the filename of the LAMMPS dumpfile (should be in directory for ease): "
  !read  *, filename
  filename = "out.Pt_3nm.20000"
  print *, "Using dumpfile: "// trim(filename)

  ! This is to get the lattice parameter and multiplies it by a constant
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
  allocate(tempArray(totalAtoms))



  ! Skips through unnecessary header data
  read(1,*) ! Skips through unnecessary header data
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)

  ! This puts all the atoms into the atomArray with their ID as their index
  do i = 1, totalAtoms
    read(1,*) ID, uselessInt, atomArray(ID, 1), atomArray(ID, 2), atomArray(ID, 3)
  end do


  ! Test if the atoms are within a cube with a side length twice the lattice
  !     parameter.
  ! For x axis test: (atomArray(i, 1) <= (atomArray(1, 1) + (2 * interactionR)))
  ! For y axis test: (atomArray(i, 2) <= (atomArray(1, 2) + (2 * interactionR)))
  ! For z axis test: (atomArray(i, 2) <= (atomArray(1, 2) + (2 * interactionR)))
  do i = 2, totalAtoms
    if (atomArray(i, 1) <= (atomArray(1, 1) + (2 * interactionR))) then
      if (atomArray(i, 2) <= (atomArray(1, 2) + (2 * interactionR))) then
        if (atomArray(i, 2) <= (atomArray(1, 2) + (2 * interactionR))) then
          tempArray(nearestIterable) = i
          nearestIterable = nearestIterable + 1
        end if
      end if
    end if
  end do

  allocate(nearestNeighbors(nearestIterable))

  do i = 0, nearestIterable
    if ((sqrt((atomArray(1, 1) - atomArray(tempArray(i), 1)) ** 2 + &
              (atomArray(1, 2) - atomArray(tempArray(i), 2)) ** 2 + &
              (atomArray(1, 3) - atomArray(tempArray(i), 3)) ** 2   &
              )) > interactionR) then
              nearestNeighbors(i) = tempArray(i)
    end if
  end do

  print *, nearestNeighbors




  ! To read a specific entry just use atomArray (idnumber, entry OR :)


 call cpu_time(toc)
 print *, "Time Taken -->", real(toc-tic)
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
