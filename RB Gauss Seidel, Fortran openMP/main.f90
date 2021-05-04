module precision

  implicit none
  
  integer, parameter :: dp = kind(1.d0)
  
end module precision

module routines

  use precision
  use omp_lib
  implicit none
  
contains
  
! ------------------------------------------------------------------------ !

  subroutine printmat(mat,string)
  real(dp), intent(in), dimension(:,:,:) :: mat
  character(*) :: string

  integer :: i,j,k

  print *, string
  DO i = 1, size(mat,1)
    DO j = 1, size(mat,2)
      DO k = 1, size(mat,3)
        write (*, "(F8.2)", ADVANCE="NO") mat(i,j,k)     
      END DO
    print *
    END DO
    print *
    print *
  END DO
  end subroutine printmat


! ------------------------------------------------------------------------ !

  subroutine write_vtk(u)
    ! This routines is taken from 3-weeks course in HPC!
    ! This is C-interoperability using bindings
    use, intrinsic :: iso_c_binding

    interface
      subroutine write_vtk_c(n, u) BIND(C, name="write_vtk")
        import c_ptr, c_int
        implicit none
        integer(c_int), value :: n
        type(c_ptr), value :: u
      end subroutine
    end interface

    ! Array to write-out
    real(dp), intent(in), target :: u(:,:,:)

    ! Local variables
    integer(c_int) :: N

    ! Get size of the array
    N = size(u, 1)

    call write_vtk_c(n, c_loc(u(1,1,1)))

  end subroutine write_vtk

! -------------------------------------------------------------------------!

  subroutine Gauss_Seidel_redblack_OMP(u,delta,N,it)
  real(dp), dimension(:,:,:), intent(inout) :: u
  integer,  intent(in) :: N,it
  real(dp), intent(in) :: Delta
  
  integer :: i,j,k
  integer :: iter
  integer :: N_loc

  !$OMP parallel default(none) private(i,j,k,iter) shared(u,delta,N,it,N_loc) 

  ! used for debugging
  print*, 'threads:',omp_get_num_threads()

  ! For the evalRadiator function:
  N_loc = size(u,1)

  ! iterating over the volume
  do iter = 1,it
    
    !$OMP do schedule(static)
    ! RED DOTS (Even coordinte sums)
    do i = 2,size(u,1)-1
      do j = 2,size(u,2)-1
        do k = 2+mod(i+j,2),size(u,3)-1,2
          u(i,j,k) =(u(i-1,j,k)+u(i+1,j,k)&!X-direction
                    +u(i,j-1,k)+u(i,j+1,k)&!Y-direction
                    +u(i,j,k-1)+u(i,j,k+1)&!Z-direction
                    +delta**2._dp*evalRadiator(i,j,k,N_loc,N)&
                     )/6._dp
        end do
      end do
    end do
    !$OMP end do

    !$OMP do schedule(static)
    ! BLACK DOTS (Uneven coordinate sums)
    do i = 2,size(u,1)-1
      do j = 2,size(u,2)-1
        do k = 2+mod(i+j+1,2),size(u,3)-1,2
          u(i,j,k) =(u(i-1,j,k)+u(i+1,j,k)&!X-direction
                    +u(i,j-1,k)+u(i,j+1,k)&!Y-direction
                    +u(i,j,k-1)+u(i,j,k+1)&!Z-direction
                    +delta**2._dp*evalRadiator(i,j,k,N_loc,N)&
                     )/6._dp
        end do
      end do
    end do
    !$OMP end do

  end do

  !$OMP end parallel
  end subroutine Gauss_Seidel_redblack_OMP
  
! ------------------------------------------------------------------------ !

  real(dp) function evalRadiator(i,j,k,N_loc,N)
    integer,  intent(in) :: i,j,k,N_loc,N
    real(dp) :: x,y,z

    x = (i+N_loc*0+1._dp)/N
    y = (j+N_loc*0+1._dp)/N
    z = (k+N_loc*0+1._dp)/N

    if (0._dp <= x .and. x <= 5._dp/16._dp .and. &
  3._dp/4._dp <= y .and. y <= 1._dp        .and. &
  1._dp/6._dp <= z .and. z <= 1._dp/2._dp) then
      evalRadiator = 200._dp !C/m^2
    else
      evalRadiator = 0._dp   !C/m^2
    end if
  end function evalRadiator

end module routines

! ------------------------------------------------------------------------ !
! --------------------------- MAIN PROGRAM START ------------------------- !
! ------------------------------------------------------------------------ !

program main

  use precision
  use routines
  use mpi
  use omp_lib
  
  implicit none

  ! Inputs:
  integer :: N! = 64 ! side length of matrix
  integer :: P! = 8 ! total nr. of workers
  integer :: c! = 2 ! collumn
  integer :: r! = 2 ! row
  integer :: l! = 2 ! layer
  integer :: algo ! selects the algorithm used.

  ! Arguments from cmd
  integer :: nargs, iargs(6)
  character(len=64) :: line

  ! File variables
  character(125) :: filename

  ! etc.
  integer :: i, flag = 0

  ! field variables  
  real(dp), allocatable, dimension(:,:,:) :: uloc
  real(dp) :: Delta

  ! timing
  real(dp) t1, t2
  real(dp), dimension(:), allocatable :: times

  integer :: ierr

  ! Defining inputs:
  nargs = command_argument_count()
  if ( nargs < 6 ) then
    write(*,'(a,i0)') 'Expected 6 arguments, only got ', nargs
    write(*,'(a)') "     N  np_r np_c np_l <algo>"
    stop
  else
    do nargs = 1, 6
      line = ' '
      call get_command_argument(nargs, line, ierr)
      read(line, *) iargs(nargs)
    end do
  end if

  N = iargs(1); P = iargs(2); r = iargs(3); c = iargs(4); l = iargs(5);
  algo = iargs(6);  

  ! using openMP to define the number of threads
  !$OMP parallel
  ! Printout the setup:
  !$OMP single
  write(*,'(8(a,i0),a)') "project: algo[", algo, &
          "] ranks=", 1, " N=", N, " Prc(",r,',',c,',',l,")"

   write(filename,'(7(a,i0),a)') 'HPCdata/N',N,'_P',P,'_px',r,'_py',c,&
         '_pz',l,'_algo',algo,'_numthreads',omp_get_num_threads(),'.dat'

  open(420,file=filename)!,status="new")
  write(420,'(8(a,i0),a)') "project: algo[", algo, &
            "] ranks=", 1, " matrix=", N," Threads=",omp_get_num_threads()," Prc(",r,',',c,',',l,")"
  !$OMP end single
  !$OMP end parallel  


  ! Allocating u
  allocate(uloc(N+2,N+2,N+2))
  
  ! Applying boundary conditions
  uloc(:,:,:)   = 0
  uloc(1,:,:)   = 20
  uloc(N+2,:,:) = 20
  uloc(:,1,:)   = 20
  uloc(:,N+2,:) = 0
  uloc(:,:,1)   = 20
  uloc(:,:,N+2) = 20

  !Defining Delta (2 is to define the 2 meters of box length):
  Delta = 2._dp/(N-1._dp)
  
  t1 = MPI_Wtime()
  select case (algo)
  case(0)
    ! Calling the Gauss Seidel routine:
    call Gauss_Seidel_redblack_OMP(uloc,Delta,N+2,1000)
  case default
    write(*,'(A,I0)') 'Unknown algorithm?? = ',algo
  end select
  t2 = MPI_Wtime()-t1  

  ! print to terminal
  write(*,*) 'Wallclock: ', t2

  ! write this to file:
  allocate(times(1))
  times(1) = t2
  write(420,*) times(i)
  close(420)
  
  ! Prints the global matrix to a VTK File.
  call write_vtk(uloc)
  
  deallocate(uloc)
end program main
