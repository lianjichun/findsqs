module outfiles
  implicit none
contains

subroutine outposcar(a,x,atom,natom,symbols,k,site,iconf,nl,l)
  use structure
  use functions
  implicit none
  real(8) :: a(3,3),x(:,:)
  integer(2) :: natom(:),k(:),iconf(:,:)
  integer(1) :: site
  character(len=2) :: atom(:),symbols(:)
  real(8),allocatable :: x0(:,:),x1(:,:),x2(:,:),x_new(:,:)
  character(len=2),allocatable :: atom_new(:)
  integer(2),allocatable :: natom_new(:),nx(:),nx0(:),nx1(:),nx2(:),aconf(:)
  integer(1) :: ntype,ntype_new,nd,nk,nKw,j
  integer(2) :: na,na0,na1,na2,na_new,n,nl,l
  integer(8) :: nc,i
  character(len=30) :: fm,tail,poscar

  na0=natom(site)
  allocate(nx0(na0))
  allocate(x0(3,na0))
  forall(i=1:na0) nx0(i)=i
  nx0=nx0+sum(natom(1:site-1))
  x0=x(:,nx0)

  ntype=size(natom)
  na=sum(natom)
  na1=na-na0
  allocate(nx1(na1))
  allocate(nx(na))
  allocate(x1(3,na1))
  nx=complement(nx0,na)
  nx1=nx(1:na1)
  x1=x(:,nx1)

  nk=size(k)
  ntype_new=ntype+nk-1
  if ( any(symbols == 'Kw') ) then
    ntype_new=ntype_new-1
  end if
  allocate(natom_new(ntype_new))
  allocate(atom_new(ntype_new))
  j=0
  do i=1,ntype
    if ( i /= site ) then
      j=j+1
      natom_new(j)=natom(i)
      atom_new(j)=atom(i)
    end if
  end do
  nKw=0
  do i=1,nk
    if ( symbols(i) /= 'Kw' ) then
      j=j+1
      natom_new(j)=k(i)
      atom_new(j)=symbols(i)
    else
      nKw=i
    end if
  end do

  na_new=sum(natom_new)
  allocate(x_new(3,na_new))
  na2=na_new-na1
  allocate(x2(3,na2))
  allocate(nx2(na2))

  call system('mkdir poscar 2> /dev/null')
  nc=size(iconf,2)
  nd=width(i_8=nc)
  write(fm,*) nd
  na=sum(k)
  allocate(aconf(na))
  do i=1,nc
    if ( nKw == nk ) then
      nx2=iconf(:,i)
    else
      aconf=complement(iconf(:,i),na)
      n=na-k(nk)
      aconf(n+1:)=aconf(1:k(nk))
      aconf(1:n)=iconf(:,i)
      if ( nkW /= 0 ) then
        n=sum(k(1:nKw))
        aconf(n-k(nKw)+1:na2)=aconf(n+1:)
      end if
      nx2=aconf(1:na2)
    end if
    x2=x0(:,nx2)
    x_new(:,1:na1)=x1
    x_new(:,na1+1:)=x2
    write(tail,'(I'//fm//'.'//fm//')') i
    poscar='poscar/POSCAR-'//tail//''
    call writepos(a,x_new,atom_new,natom_new,poscar)
  end do
  nd=width(i_2=nl)
  write(fm,*) nd
  write(tail,'(I'//fm//'.'//fm//')') l
  poscar='poscar-'//tail//''
  call system('mv poscar '//poscar//'')
  return
end subroutine outposcar

end module outfiles
