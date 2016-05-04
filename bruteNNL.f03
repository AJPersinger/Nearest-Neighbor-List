! Title: Nearest Neigbhor List
! Purpose: Implement nearest neigbhor searches and creates a list of nearest atoms
!           in a structure. Look for the explanation in the README
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
  real :: tic, toc, x1, y1, z1, x2, y2, z2, distance
  real*8 :: interactionR
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


  ! Nearest Neighbor calculations
  do i = 1, totalAtoms
    write (2,*) "Nearest neighbors for atom ID: ", i, " are: "
    ! Assign coordinates to variables (Just to make the code less cluttered)
    x1 = atomArray(i, 1)
    y1 = atomArray(i, 2)
    z1 = atomArray(i, 3)

    do k = 1, totalAtoms
      ! Assign coordinates to variables (Just to make the code less cluttered)
      x2 = atomArray(k, 1)
      y2 = atomArray(k, 2)
      z2 = atomArray(k, 3)

      ! Check to see if the distances squared is within the interactionR
      distance = pnearestNeighborCalc(xFloor, xCeil, yFloor, yCeil, zFloor, zCeil, x1, y1, z1, x2, y2, z2, interactionR)
      if (distance /= -1) then
        ! Write the data to the file
        write (2,*) "ID: ", k, "   Distance: ", sqrt(distance)
      end if
    end do
    write(2,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  end do



 ! Close the file to give permission to modify to whatever other program
 close (2, status='KEEP')

 call cpu_time(toc)
 print *, "Time Taken:", real(toc-tic)
end program nearestNeighborsLAMMPS

!START FUNCTION DECLARATIONS
! -squareDistance:
!     return type: real
!     purpose: calculate the distance of two points before taking the sqrt
!
! -isNN:
!     return type: logical
!     purpose: determine if two atoms are neighbors
!
!END FUNCTION DECLARATIONS



real function squareDistance(x1, y1, z1, x2, y2, z2)
  real x1, y1, z1, x2, y2, z2
  squareDistance = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2
  return
end function

! Why is there a p? Because fortran though it was a good idea for intrinsic variable names to be a thing
real function pnearestNeighborCalc(xFloor, xCeil, yFloor, yCeil, zFloor, zCeil, x1, y1, z1, x2, y2, z2, interactionR)

  if ((x1 + interactionR) > xCeil .or. (x1  - interactionR) < xFloor) then
    if ((x1 + interactionR) > xCeil) then
      x1 = x1 + interactionR
    else if ((x1 - interactionR) < xFloor) then
      x1 = x1 - interactionR
    end if
    pnearestNeighborCalc = squareDistance(x1, y1, z1, x2, y2, z2)

  else if ((y1 + interactionR) > yCeil .or. (y1  - interactionR) < yFloor) then
    if ((y1 + interactionR) > yCeil) then
      x1 = x1 + interactionR
    else if ((y1 - interactionR) < yFloor) then
      y1 = y1 - interactionR
    end if
    pnearestNeighborCalc = squareDistance(x1, y1, z1, x2, y2, z2)

  else if ((z1 + interactionR) > xCeil .or. (z1  - interactionR) < zFloor) then
    if ((z1 + interactionR) > zCeil) then
      z1 = z1 + interactionR
    else if ((z1 - interactionR) < zFloor) then
      z1 = z1 - interactionR
    end if
    pnearestNeighborCalc = squareDistance(x1, y1, z1, x2, y2, z2)

  else
    if (squareDistance(x1, y1, z1, x2, y2, z2) <= interactionR**2) then
      pnearestNeighborCalc = squareDistance(x1, y1, z1, x2, y2, z2)
    else
      pnearestNeighborCalc = -1
    end if

  end if
  return
end function
