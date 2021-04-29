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
        !if (i == 1) then
        !print*,cr
        !print*,si,sj,sk
        !print*,cr(1)*(si)+2,(cr(1)+1)*(si)+1
        !print*,cr(2)*(sj)+2,(cr(2)+1)*(sj)+1
        !print*,cr(3)*(sk)+2,(cr(3)+1)*(sk)+1
        ! Inserting the the recieved mat into the global temp field
        u_glob(cr(1)*(si)+2:(cr(1)+1)*(si)+1,&
               cr(2)*(sj)+2:(cr(2)+1)*(sj)+1,&
               cr(3)*(sk)+2:(cr(3)+1)*(sk)+1) = &
               uloc_r(2:si+1,2:sj+1,2:sk+1)
        !end if
      end do   

      ! Prints the global heat field to terminal  
      call printmat(u_glob,'Global heat field')

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

  subroutine Gauss_Seidel_redblack(u,comm,delta,irank,coords)
  real(dp), dimension(:,:,:), intent(inout) :: u
  integer,  intent(in) :: comm,irank,coords(3)
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
                  +delta**2._dp*evalRadiator(i,j,k,N_loc,delta,coords)&
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
                  +delta**2._dp*evalRadiator(i,j,k,N_loc,delta,coords)&
                   )/6._dp
      end do
    end do
  end do
  
  ! Update walls:
  call updateComm(u,irank,comm)
  end subroutine Gauss_Seidel_redblack
  
! ------------------------------------------------------------------------ !

  real(dp) function evalRadiator(i,j,k,N_loc,delta,coords)
    integer,  intent(in) :: i,j,k,N_loc,coords(3)
    real(dp), intent(in) :: delta
    real(dp) :: x,y,z

    x = i*delta+(N_loc-2)*delta*coords(1)
    y = j*delta+(N_loc-2)*delta*coords(2)
    z = k*delta+(N_loc-2)*delta*coords(3)

    if (-1._dp <= x .and. x <= -3._dp/8._dp .and. &
        -1._dp <= y .and. y <= -1._dp/2._dp .and. &
      -2._dp/3._dp <= z .and. z <= 0._dp) then
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
  integer :: rank_prev, rank_next
  integer :: i,j,k
  integer :: req(12),stats(MPI_STATUS_SIZE,12),ierr,errc

  ! Allocating the buffer for the communication:
  allocate(sBuff(size(u,1),size(u,2),6),rBuff(size(u,1),size(u,2),6))

  call mpi_barrier(comm,ierr)
  if (irank == 1) then
    call printmat(u,'uloc for irank == 1')
  end if
  

  ! Send the calculated points:
  do i = 1,3
    call MPI_Cart_shift(comm, i-1, 1, rank_prev, rank_next, ierr)
    
    
    ! ---------------------------- SENDING ------------------------------- !  

    ! Saving Rank_prev into buffer
    if (rank_prev .NE. MPI_PROC_NULL) then
      if     (i == 1) then      
        sBuff(:,:,i) = u(2,:,:)
      elseif (i == 2) then
        sBuff(:,:,i) = u(:,2,:)
      elseif (i == 3) then
        sBuff(:,:,i) = u(:,:,2)
      end if

      ! Sending the data:
      call Mpi_Isend(sBuff(:,:,i), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rank_prev, 0, comm, req(i), ierr)
    else 
      req(i) = MPI_REQUEST_NULL
    end if
    
    ! Saving Rank_next into buffer
    if (rank_next .NE. MPI_PROC_NULL) then
      if     (i == 1) then      
        sBuff(:,:,i+3) = u(size(u,1)-1,:,:)
      elseif (i == 2) then
        sBuff(:,:,i+3) = u(:,size(u,2)-1,:)
      elseif (i == 3) then
        sBuff(:,:,i+3) = u(:,:,size(u,3)-1)
      end if
      ! Sending the data:
      call Mpi_Isend(sBuff(:,:,i+3), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rank_next, 1, comm, req(i+3), ierr)
    else 
      req(i+3) = MPI_REQUEST_NULL
    end if
    
    ! --------------------------- RECIEVING ------------------------------ !  

    ! Recieving Rank_prev into u
    if (rank_prev .NE. MPI_PROC_NULL) then
      call Mpi_Irecv(rBuff(:,:,i), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rank_prev, 1, comm, req(i+6), ierr)
    else 
      req(i+6) = MPI_REQUEST_NULL
    end if
    
    ! Saving Rank_next into buffer
    if (rank_next .NE. MPI_PROC_NULL) then
      call Mpi_Irecv(rBuff(:,:,i+3), size(u,1)**2, MPI_DOUBLE_PRECISION, &
                     rank_prev, 1, comm, req(i+9), ierr)
    else 
      req(i+9) = MPI_REQUEST_NULL
    end if
  end do


  ! Waiting for all messages to be recieved
  call Mpi_Waitall(12, req, stats, ierr)

  ! Once data is recieved are the information inserted into the array.
  do i = 1,3
    call MPI_Cart_shift(comm, i-1, 1, rank_prev, rank_next, ierr)
    
    if (rank_prev .NE. MPI_PROC_NULL) then
      if     (i == 1) then      
        u(1,:,:) = rBuff(:,:,i)
      elseif (i == 2) then
        u(:,1,:) = rBuff(:,:,i)
      elseif (i == 3) then
        u(:,:,1) = rBuff(:,:,i)
      end if
    end if
    
    if (rank_next .NE. MPI_PROC_NULL) then
      if     (i == 1) then      
        u(size(u,1),:,:) = rBuff(:,:,i+3)
      elseif (i == 2) then
        u(:,size(u,2),:) = rBuff(:,:,i+3)
      elseif (i == 3) then
        u(:,:,size(u,3)) = rBuff(:,:,i+3)
      end if
    end if

  end do

  call mpi_barrier(comm,ierr)
  if (irank == 1) then
    call printmat(u,'uloc for irank == 1, after!')
    call MPI_abort(comm,errc,ierr)
  end if

  end subroutine updateComm
end module routines

! ------------------------------------------------------------------------ !
! --------------------------- MAIN PROGRAM START ------------------------- !
! ------------------------------------------------------------------------ !

program main

  use precision
  use routines
  use mpi
  
  implicit none

  integer :: N = 8 ! side length of matrix
  integer :: P = 8 ! total nr. of workers
  integer :: c = 2 ! collumn
  integer :: r = 2 ! row
  integer :: l = 2 ! layer

  integer :: i

  ! field variables  
  real(dp), allocatable, dimension(:,:,:) :: uloc
  real(dp) :: Delta

  ! MPI Variables
  integer :: isize,irank,ierr

  ! Cartesian
  integer cart_comm, ndims
  integer dim_size(3)
  logical periods(3), reorder
  integer coords(3)

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, isize,ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, irank,ierr)

  if (irank == 0) then
    if (.not. p==(r*c*l))  print*, 'ERROR this should hold: p==(r*c*l)'
    if (.not. p==isize)    print*, 'ERROR this should hold: p==isize'
    if (.not. mod(N,p)==0) print*, 'ERROR this should hold: mod(N,p)==0'
    if (.not. mod(N,c)==0) print*, 'ERROR this should hold: mod(N,c)==0'
    if (.not. mod(N,r)==0) print*, 'ERROR this should hold: mod(N,r)==0'    
    if (.not. mod(N,l)==0) print*, 'ERROR this should hold: mod(N,r)==0'
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

  !do i = 1,1000
    ! Calling the Gauss Seidel routine:
    call Gauss_Seidel_redblack(uloc,cart_comm,Delta,irank,coords)
  !end do

  ! Prints the global matrix to a VTK File.
  call vtk_u_global(uloc,coords,irank,isize,cart_comm)
  
  deallocate(uloc)
  CALL MPI_Finalize(ierr)

end program main
