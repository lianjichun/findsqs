program pos2atat
  implicit none
  real(8) :: a(3,3)
  real(8),allocatable :: x(:,:)
  integer(2),allocatable :: natom(:)
  character(len=2),allocatable :: symbol(:)
  character(len=20) :: poscar
  integer(2) :: i,j

  read(*,*) poscar
  call readpos(a,x,symbol,natom,poscar)
  write(*,'(3F10.6)') (a(:,i),i=1,3)
  write(*,'(A)') '  1.000000  0.000000  0.000000'
  write(*,'(A)') '  0.000000  1.000000  0.000000'
  write(*,'(A)') '  0.000000  0.000000  1.000000'
  j=1
  do i=1,sum(natom)
    write(*,'(3F10.6,A4)') x(:,i),symbol(j)
    if ( i == sum(natom(1:j)) ) j=j+1
  end do
  write(*,'(A)') 'end'
  write(*,*)
contains

subroutine readpos(a,x,symbol,natom,poscar)
  implicit none
  real(8) :: a(3,3)
  real(8),allocatable :: x(:,:)
  integer(1) :: ntype,io_error
  integer(2) :: satom,i
  integer(2),allocatable :: natom(:)
  real(8) :: scal
  character(len=2),allocatable :: symbol(:)
  character(len=1) :: coordinate
  character(len=*) :: poscar

  open(10,file=poscar)
  read(10,*)
  read(10,*) scal
  read(10,*) a
  a=scal*a
  read(10,*)
  ntype=1
  io_error=0
  do while (io_error == 0)
    allocate(natom(ntype))
    read(10,*,iostat=io_error) natom
    ntype=ntype+1
    deallocate(natom)
    backspace(10)
  end do
  ntype=ntype-2
  allocate(symbol(ntype))
  allocate(natom(ntype))
  backspace(10)
  backspace(10)
  read(10,*) symbol
  read(10,*) natom
  satom=sum(natom)
  allocate(x(3,satom))
  read(10,*) coordinate
  if ( coordinate == 'S' ) then
    read(10,*) coordinate
  end if
  do i=1,satom
    read(10,*) x(:,i)
  end do
  close(10)
  if ( coordinate == 'C' ) then
    x=matmul(inverse(a),x)
  end if
  return
end subroutine readpos

function inverse(A)
  implicit none
  real(8) :: inverse(3,3)
  real(8) :: A(3,3)
  integer(1) :: ind(2,3),i,j

  ind(:,1)=(/ 2, 3 /)
  ind(:,2)=(/ 1, 3 /)
  ind(:,3)=(/ 1, 2 /)
  do i=1,3
    do j=1,3
      inverse(j,i)=(-1)**(i+j)*det2(A(ind(:,i),ind(:,j)))
    end do
  end do
  inverse=inverse/det3(A)
  return
end function inverse

function det2(A)
  implicit none
  real(8) :: det2,A(2,2)

  det2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
  return
end function det2

function det3(A)
  implicit none
  real(8) :: det3,A(3,3)

  det3=     A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  det3=det3-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  det3=det3+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  return
end function det3

end program pos2atat
