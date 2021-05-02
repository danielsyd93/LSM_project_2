module precision

  implicit none
  
  integer, parameter :: dp = kind(1.d0)
  
end module precision

module routines

  use precision
  use mpi
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

  subroutine ApplyBC(dim,pos,BCval,coords,irank,u,comm)
  integer, intent(in) :: dim,pos,BCval,coords(3),irank,comm
  real(dp), intent(inout), dimension(:,:,:) :: u
  
  integer :: ierr
  integer :: max_coords,size_u

  ! Finding the max dimension:
  call MPI_allreduce(maxval(coords), max_coords, 1, MPI_DOUBLE_PRECISION, &
                     MPI_MAX, comm, ierr)


  if (pos == 0 .AND. coords(dim+1) == 0) then
    if     (dim == 0) then
      u(1,:,:) = BCval
    elseif (dim == 1) then
      u(:,1,:) = BCval
    elseif (dim == 2) then
      u(:,:,1) = BCval
    else
      print*, 'Error, dim is 0-indexed i.e. dim is [0,2].'
    end if
  end if

  if (pos .ne. 0 .and. coords(dim+1) == max_coords) then
    size_u = size(u,dim+1)
    
    if     (dim == 0) then
      u(size_u,:,:) = BCval
    elseif (dim == 1) then
      u(:,size_u,:) = BCval
    elseif (dim == 2) then
      u(:,:,size_u) = BCval
    else
      print*, 'Error, dim is 0-indexed i.e. dim is [0,2].'
    end if

  end if

  end subroutine ApplyBC
! ------------------------------------------------------------------------ !

  subroutine vtk_u_global(uloc,coords,irank,isize,comm)
    real(dp),intent(in),dimension(:,:,:) :: uloc
    integer :: coords(3),irank,isize,comm

    real(dp), dimension(:,:,:), allocatable :: u_glob,uloc_r
    integer :: max_coords(3),cr(3)
    integer :: i,j,k
    integer :: si,sj,sk
    integer :: ierr
    integer :: stat(MPI_STATUS_SIZE) 

    ! Finding the max dimension:
    call MPI_allreduce(coords, max_coords, 3, MPI_INTEGER, &
                       MPI_MAX, comm, ierr)

    call mpi_barrier(comm,ierr)

    if (irank == 0) then
      si = size(uloc,1)-2; sj = size(uloc,2)-2; sk = size(uloc,3)-2

      ! Defining the global solution
      allocate(u_glob(si*(max_coords(1)+1)+2,&
                      sj*(max_coords(2)+1)+2,&
                      sk*(max_coords(3)+1)+2))

      allocate(uloc_r(si+2,sj+2,sk+2))
      
      ! Applying BC's according with the local elements
      u_glob = 20
      u_glob(:,size(u_glob,2),:) = 0
      u_glob(2:si+1,2:sj+1,2:sk+1) = &
        uloc(2:si+1,2:sj+1,2:sk+1)
      
      do i = 1,isize-1
        ! recieve coordinates:
        call MPI_RECV(cr, 3, MPI_INTEGER, i, 1, comm,&
                      stat, ierr)

        ! Recieve local temp field
        call MPI_RECV(uloc_r, size(uloc), MPI_DOUBLE_PRECISION, i, 2, &
                      comm, stat, ierr)

        ! Inserting the the recieved mat into the global temp field
        u_glob(cr(1)*(si)+2:(cr(1)+1)*(si)+1,&
               cr(2)*(sj)+2:(cr(2)+1)*(sj)+1,&
               cr(3)*(sk)+2:(cr(3)+1)*(sk)+1) = &
               uloc_r(2:si+1,2:sj+1,2:sk+1)
        !end if
      end do   

      ! Prints the global heat field to terminal  
      !call printmat(u_glob,'Global heat field')

      ! finally creating a VTK file:
      call write_vtk(u_glob)


    else
      do i = 1,isize-1
        if (irank == i) then
          ! send coordinates:
          call MPI_SEND(coords, 3, MPI_INTEGER, 0, 1, comm, ierr)
          
          ! send the local matrix:
          call MPI_SEND(uloc,size(uloc),MPI_DOUBLE_PRECISION,0,2,comm,ierr)
        end if
      end do 
    end if    
  end subroutine vtk_u_global
  
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

  subroutine Gauss_Seidel_redblack(u,comm,delta,irank,coords,N)
  real(dp), dimension(:,:,:), intent(inout) :: u
  integer,  intent(in) :: comm,irank,coords(3),N
  real(dp), intent(in) :: Delta
  
  integer :: i,j,k
  integer :: N_loc
  integer :: ierr

  ! For the evalRadiator function:
  N_loc = size(u,1)

  ! RED DOTS (Even coordinate sums)
  do i = 2,size(u,1)-1
    do j = 2,size(u,2)-1
      do k = 2+mod(i+j,2),size(u,3)-1,2
        u(i,j,k) =(u(i-1,j,k)+u(i+1,j,k)&!X-direction
                  +u(i,j-1,k)+u(i,j+1,k)&!Y-direction
                  +u(i,j,k-1)+u(i,j,k+1)&!Z-direction
                  +delta**2._dp*evalRadiator(i,j,k,N_loc,coords,N)&
                   )/6._dp
      end do
    end do
  end do

  ! Update walls:
  call updateComm(u,irank,comm)

  ! BLACK DOTS (Uneven coordinate sums)
  do i = 2,size(u,1)-1
    do j = 2,size(u,2)-1
      do k = 2+mod(i+j+1,2),size(u,3)-1,2
        u(i,j,k) =(u(i-1,j,k)+u(i+1,j,k)&!X-direction
                  +u(i,j-1,k)+u(i,j+1,k)&!Y-direction
                  +u(i,j,k-1)+u(i,j,k+1)&!Z-direction
                  +delta**2._dp*evalRadiator(i,j,k,N_loc,coords,N)&
                   )/6._dp
      end do
    end do
  end do
  
  ! Update walls:
  call updateComm(u,irank,comm)
  end subroutine Gauss_Seidel_redblack
  
! ------------------------------------------------------------------------ !

subroutine Jacobi(u,u_old,coords,irank,N)
real(dp), dimension(:,:,:), intent(inout) :: u
real(dp), dimension(:,:,:), intent(inout) :: u_old
integer,  intent(in) :: comm,irank,coords(3),N
real(dp), intent(in) :: Delta

integer :: i,j,k
integer :: N_loc
integer :: ierr
real(dp) :: temp

do i = 2,size(u,1)-1
    do j = 2,size(u,2)-1
      do k = 2+mod(i+j,2),size(u,3)-1,2
        u_old(i,j,k) =(u(i-1,j,k)+u(i+1,j,k)&!X-direction
                  +u(i,j-1,k)+u(i,j+1,k)&!Y-direction
                  +u(i,j,k-1)+u(i,j,k+1)&!Z-direction
                  +delta**2._dp*evalRadiator(i,j,k,N_loc,coords,N)&
                   )/6._dp
      end do
    end do
  end do


temp=u
u=u_old
u_old=temp

end subroutine Jacobi

! ------------------------------------------------------------------------ !
  real(dp) function evalRadiator(i,j,k,N_loc,coords,N)
    integer,  intent(in) :: i,j,k,N_loc,coords(3),N
    real(dp) :: x,y,z

    x = (i+N_loc*coords(1)+1._dp)/N
    y = (j+N_loc*coords(2)+1._dp)/N
    z = (k+N_loc*coords(3)+1._dp)/N

    if (0._dp <= x .and. x <= 5._dp/16._dp .and. &
  3._dp/4._dp <= y .and. y <= 1._dp        .and. &
  1._dp/6._dp <= z .and. z <= 1._dp/2._dp) then
      evalRadiator = 200._dp !C/m^2
    else
      evalRadiator = 0._dp   !C/m^2
    end if
  end function evalRadiator

! ------------------------------------------------------------------------ !
! ---------------------------- COMMUNICATION ----------------------------- !
! ------------------------------------------------------------------------ !
  subroutine updateComm(u,irank,comm)
  real(dp), dimension(:,:,:), intent(inout) :: u
  integer, intent(in) :: irank,comm

  real(dp), dimension(:,:,:), allocatable :: sBuff,rBuff
  integer :: rp, rn ! rank prev. and rank next
  integer :: NumNe ! number of neigbors
  integer :: l1,lend
  integer :: i,j,k,c

  integer :: ierr,errc
  integer, dimension(:),   allocatable :: sreq,rreq
  integer, dimension(:,:), allocatable :: sstats,rstats

  ! Number of Neighbors:
  NumNe = 0
  
  ! Getting number of neighbors:
  do i = 1,3
    call MPI_Cart_shift(comm, i-1, 1, rp, rn, ierr)

    if (rp .ne. MPI_PROC_NULL) NumNe = NumNe + 1

    if (rn .ne. MPI_PROC_NULL) NumNe = NumNe + 1
  end do

  ! Allocating the buffer for the communication:
  allocate(sBuff(size(u,1),size(u,2),NumNe),&
           rBuff(size(u,1),size(u,2),NumNe))

  ! Allocating the requests
  allocate(rreq(NumNe),sreq(NumNe))  

  ! Allocating the statuses
  allocate(sstats(MPI_STATUS_SIZE,NumNe),rstats(MPI_STATUS_SIZE,NumNe))

  ! Now the buffers are filled
  c = 0 ! Counter for the neigbors 
  do i = 1,3
    call MPI_Cart_shift(comm, i-1, 1, rp, rn, ierr)
    ! Saving Rank_prev into buffer
    if (rp .NE. MPI_PROC_NULL) then
      c = c + 1

      if     (i == 1) then      
        sBuff(:,:,c) = u(2,:,:)
      elseif (i == 2) then
        sBuff(:,:,c) = u(:,2,:)
      elseif (i == 3) then
        sBuff(:,:,c) = u(:,:,2)
      end if

      ! Sending the data:
      call Mpi_Isend(sBuff(:,:,c), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rp, 0, comm, sreq(c), ierr)

      ! Recieving the data:
      call Mpi_Irecv(rBuff(:,:,c), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rp, 1, comm, rreq(c), ierr)
    end if
    
    ! Saving Rank_next into buffer
    if (rn .NE. MPI_PROC_NULL) then
      c = c + 1

      if     (i == 1) then      
        sBuff(:,:,c) = u(size(u,1)-1,:,:)
      elseif (i == 2) then
        sBuff(:,:,c) = u(:,size(u,2)-1,:)
      elseif (i == 3) then
        sBuff(:,:,c) = u(:,:,size(u,3)-1)
      end if
      ! Sending the data:
      call Mpi_Isend(sBuff(:,:,c), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rn, 1, comm, sreq(c), ierr)

      ! Recieving the data:
      call Mpi_Irecv(rBuff(:,:,c), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rn, 0, comm, rreq(c), ierr)
    end if
  end do

  ! Waiting for all messages to be sent
  call Mpi_Waitall(c, sreq, sstats, ierr)

  ! Waiting for all messages to be recieved
  call Mpi_Waitall(c, rreq, rstats, ierr)

  ! Back inserting the step:
  c = 0
  do i = 1,3
    call MPI_Cart_shift(comm, i-1, 1, rp, rn, ierr)
    
    if (rp .NE. MPI_PROC_NULL) then
      c = c + 1
      if     (i == 1) then      
        u(1,:,:) = rBuff(:,:,c)
      elseif (i == 2) then
        u(:,1,:) = rBuff(:,:,c)
      elseif (i == 3) then
        u(:,:,1) = rBuff(:,:,c)
      end if
    end if
    
    if (rn .NE. MPI_PROC_NULL) then
      c = c + 1
      if     (i == 1) then      
        u(size(u,1),:,:) = rBuff(:,:,c)
      elseif (i == 2) then
        u(:,size(u,2),:) = rBuff(:,:,c)
      elseif (i == 3) then
        u(:,:,size(u,3)) = rBuff(:,:,c)
      end if
    end if

  end do

  end subroutine updateComm
  
  subroutine assert(condition,string,flag)
  logical, intent(in)      :: condition
  character(*), intent(in) :: string
  integer, intent(out)     :: flag

  if (condition) then
    print*, string
    flag = 1
  end if

  end subroutine assert
end module routines

! ------------------------------------------------------------------------ !
! --------------------------- MAIN PROGRAM START ------------------------- !
! ------------------------------------------------------------------------ !

program main

  use precision
  use routines
  use mpi
  
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

  ! MPI Variables
  integer :: isize,irank,ierr,errc

  ! Cartesian
  integer cart_comm, ndims
  integer dim_size(3)
  logical periods(3), reorder
  integer coords(3)

  ! timing
  real(dp) t1, t2
  real(dp), dimension(:), allocatable :: times

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, isize,ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, irank,ierr)

  ! Defining inputs:
  nargs = command_argument_count()
  if ( nargs < 6 ) then
    call MPI_Finalize(ierr)
    write(*,'(a,i0)') 'Expected 6 arguments, only got ', nargs
    write(*,'(a)') "     N  np_r np_c np_l <algo>"
  else
    do nargs = 1, 6
      line = ' '
      call get_command_argument(nargs, line, ierr)
      read(line, *) iargs(nargs)
    end do
  end if

  ! Distribute arguments
  call MPI_Bcast(iargs, 6, MPI_Integer, 0, MPI_COMM_WORLD, ierr)

  N = iargs(1); P = iargs(2); r = iargs(3); c = iargs(4); l = iargs(5);
  algo = iargs(6);  

  call MPI_barrier(MPI_COMM_WORLD,ierr)

  ! Printout the setup:
  if (irank == 0) then
    write(*,'(7(a,i0),a)') "project: algo[", algo, &
            "] ranks=", isize, " N=", N, " Prc(",r,',',c,',',l,")"

     write(filename,'(6(a,i0),a)') 'HPCdata/N',N,'_P',P,'_px',r,'_py',c,&
                                   '_pz',l,'_algo',algo,'.dat'

    open(420,file=filename,status="new")
    write(420,'(7(a,i0),a)') "project: algo[", algo, &
              "] ranks=", isize, " matrix=", N, " Prc(",r,',',c,',',l,")"
  end if

  ! check that the correct conditions are given.
  if (irank == 0) then
    call assert((.not. p==(r*c*l)), 'ERROR this should hold: p==(r*c*l)' ,flag)
    call assert((.not. p==isize),   'ERROR this should hold: p==isize'   ,flag)
    call assert((.not. mod(N,p)==0),'ERROR this should hold: mod(N,p)==0',flag)
    call assert((.not. mod(N,c)==0),'ERROR this should hold: mod(N,c)==0',flag)
    call assert((.not. mod(N,r)==0),'ERROR this should hold: mod(N,r)==0',flag)
    call assert((.not. mod(N,l)==0),'ERROR this should hold: mod(N,r)==0',flag)
    if (flag == 1) then
      write(*,*) 'An error occoured in the assesment of the settings!'
      write(420,*) 'An error occoured in the assesment of the settings!'
      call MPI_abort(MPI_COMM_WORLD,errc,ierr)
      stop
    end if
  end if

  

  ! First a Cartesian grid is generated:
  ndims = 3
  dim_size = 2
  periods = .false.
  reorder = .true.
  CALL MPI_Cart_create(MPI_COMM_WORLD,ndims,dim_size,periods,reorder&
                       ,cart_comm,ierr)

  ! Now getting the new reordered iranks:
  CALL MPI_Comm_rank(cart_comm, irank,ierr)
  ! and the corresponding coords:
  CALL MPI_Cart_Coords(cart_comm, irank, 3,coords,ierr)

  ! Allocating u
  allocate(uloc(N/l+2,N/c+2,N/r+2))
  
  ! Applying boundary conditions
    ! u(:,:,:) = 0
  uloc(:,:,:) = 0
    ! u(1,:,:) = 0
  call ApplyBC(0,0,20,coords,irank,uloc,cart_comm)
    ! u(N,:,:) = 0
  call ApplyBC(0,N,20,coords,irank,uloc,cart_comm)
    ! u(:,1,:) = 0
  call ApplyBC(1,0,20,coords,irank,uloc,cart_comm)
    ! u(:,:,1) = 0
  call ApplyBC(2,0,20,coords,irank,uloc,cart_comm)
    ! u(:,:,N) = 0
  call ApplyBC(2,N,20,coords,irank,uloc,cart_comm)
    ! u(:,N,:) = 20
  call ApplyBC(1,N,0,coords,irank,uloc,cart_comm)

  !Defining Delta (2 is to define the 2 meters of box length):
  Delta = 2._dp/(N-1._dp)

  t1 = MPI_Wtime()
  do i = 1,1000
    select case (algo)
    case(0)
      ! Calling the Gauss Seidel routine:
      call Gauss_Seidel_redblack(uloc,cart_comm,Delta,irank,coords,N+2)
    case default
      write(*,'(A,I0)') 'Unknown algorithm?? = ',algo
    end select
  end do
  t2 = MPI_Wtime()-t1  

  ! print to terminal
  do i = 0,isize-1  
    if (irank==i) then
      write(*,*) 'Wallclock: ', t2
    end if
    call mpi_barrier(cart_comm,ierr)
  end do

  ! write this to file:
  allocate(times(isize))
  call MPI_GATHER(t2, 1, MPI_DOUBLE_PRECISION, times, 1,&
                  MPI_DOUBLE_PRECISION, 0, cart_comm, ierr)
  if (irank == 0) then
    do i = 1,isize
      write(420,*) times(i)
    end do    
    close(420)
  end if
  

  ! Prints the global matrix to a VTK File.
  !call vtk_u_global(uloc,coords,irank,isize,cart_comm)
  
  deallocate(uloc)
  CALL MPI_Finalize(ierr)
end program main
