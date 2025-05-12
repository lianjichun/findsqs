module configurations
  use functions
  implicit none
contains

subroutine irrconfig(iconf,eqamat,group,k,clus2,deg2,clus3,deg3,tol,npmc,nsqs,tcorr)
  use cluster
  implicit none
  type :: datalink
    integer(2),allocatable :: iconf(:)
    type(datalink),pointer :: next
  end type datalink
  real(8),parameter :: eps=1D-5
  real(8) :: tcorr(:,:),tol
  real(8),allocatable :: corr(:,:),dcorr(:)
  integer(2) :: eqamat(:,:),group(:),k(:),deg3(:)
  integer(2) :: clus2(:,:,:),deg2(:),clus3(:,:,:)
  integer(2),allocatable :: aconf(:),oconf(:),lconf(:)
  integer(2),allocatable :: ib(:),ie(:),m(:),iconf(:,:)
  integer(8),allocatable :: E(:,:),ncm(:),nct(:)
  integer(1) :: nk,nsubc,npmc
  integer(2) :: na,ns,no,ng,g,o,ncorr,nsqs
  integer(8) :: i2,i3,i4,i5,i6,ic,i,nirr
  character(len=20) :: fm
  logical(1) :: break
  logical(1),allocatable :: occ2(:),occ3(:),occ4(:)
  logical(1),allocatable :: occ5(:),occ6(:),valid(:,:)
  type(datalink),pointer :: p,head,next

  nk=size(k)-1
  na=sum(k)
  ns=na-k(nk+1)
  no=size(eqamat,2)
  ng=size(group)-1
  allocate(ncm(nk))
  allocate(nct(ng))
  allocate(m(nk))
  allocate(ib(nk))
  allocate(ie(nk))
  allocate(aconf(ns))
  allocate(oconf(ns))
  allocate(lconf(ns))
  allocate(valid(no,nk))
  call confnum(ncm,nct,ib,ie,m,k,group,ng)
  call matrixE(E,m,k)
  call lastconf(lconf,m,k)
  if ( nk > 1 ) allocate(occ3(ncm(2)))
  if ( nk > 2 ) allocate(occ4(ncm(3)))
  if ( nk > 3 ) allocate(occ5(ncm(4)))
  if ( nk > 4 ) allocate(occ6(ncm(5)))
  allocate(head)
  nullify(head%next)
  p=>head
  ncorr=size(tcorr,1)
  nsubc=nk*(nk+1)/2 !2**nk-1
  allocate(corr(ncorr,nsubc))
  allocate(dcorr(ncorr))
  write(fm,*) ncorr
  nirr=0

  select case(nk)
  case(1)
    call binary
  case(2)
    call ternary
  case(3)
    call quaternary
  case(4)
    call quinary
  case(5)
    call senary
  end select

  allocate(iconf(ns,nirr))
  nullify(p%next)
  i=0
  p=>head
  do while ( associated(p%next) )
    p=>p%next
    i=i+1
    iconf(:,i)=p%iconf
  end do
  p=>head
  do while ( associated(p) )
    next=>p%next
    deallocate(p)
    p=>next
  end do
  return
contains

subroutine binary
  do g=1,ng
    allocate(occ2(nct(g)))
    occ2=.true.
    do i2=1,nct(g)
      if ( occ2(i2) .eqv. .false. ) cycle
      call component1
      call calcorr(1,aconf(ib(1):ie(1)))
      if (  break ) cycle
      call trialclose
      if ( break ) return
    end do
    deallocate(occ2)
  end do
  return
end subroutine binary

subroutine ternary
  do g=1,ng
    allocate(occ2(nct(g)))
    occ2=.true.
    do i2=1,nct(g)
      if ( occ2(i2) .eqv. .false. ) cycle
      call component1
      call calcorr(1,aconf(ib(1):ie(1)))
      if (  break ) cycle
      occ3=.true.
      do i3=1,ncm(2)
        if ( occ3(i3) .eqv. .false. ) cycle
        call components(2,i3,occ3,ib(2),ie(2))
        call calcorr(2,aconf(ib(2):ie(2)))
        if ( break ) cycle
        call calcorr(3,aconf(ib(1):ie(2)))
        if ( break ) cycle
        call trialclose
        if ( break ) return
      end do
    end do
    deallocate(occ2)
  end do
  return
end subroutine ternary

subroutine quaternary
  do g=1,ng
    allocate(occ2(nct(g)))
    occ2=.true.
    do i2=1,nct(g)
      if ( occ2(i2) .eqv. .false. ) cycle
      call component1
      call calcorr(1,aconf(ib(1):ie(1)))
      if ( break ) cycle
      occ3=.true.
      do i3=1,ncm(2)
        if ( occ3(i3) .eqv. .false. ) cycle
        call components(2,i3,occ3,ib(2),ie(2))
        call calcorr(2,aconf(ib(2):ie(2)))
        if ( break ) cycle
        call calcorr(3,aconf(ib(1):ie(2)))
        if ( break ) cycle
        occ4=.true.
        do i4=1,ncm(3)
          if ( occ4(i4) .eqv. .false. ) cycle
          call components(3,i4,occ4,ib(3),ie(3))
          call calcorr(4,aconf(ib(3):ie(3)))
          if ( break ) cycle
          call calcorr(5,aconf(ib(2):ie(3)))
          if ( break ) cycle
          call calcorr(6,aconf(ib(1):ie(3)))
          if ( break ) cycle
          call trialclose
          if ( break ) return
        end do
      end do
    end do
    deallocate(occ2)
  end do
  return
end subroutine quaternary

subroutine quinary
  do g=1,ng
    allocate(occ2(nct(g)))
    occ2=.true.
    do i2=1,nct(g)
      if ( occ2(i2) .eqv. .false. ) cycle
      call component1
      call calcorr(1,aconf(ib(1):ie(1)))
      if ( break ) cycle
      occ3=.true.
      do i3=1,ncm(2)
        if ( occ3(i3) .eqv. .false. ) cycle
        call components(2,i3,occ3,ib(2),ie(2))
        call calcorr(2,aconf(ib(2):ie(2)))
        if ( break ) cycle
        call calcorr(3,aconf(ib(1):ie(2)))
        if ( break ) cycle
        occ4=.true.
        do i4=1,ncm(3)
          if ( occ4(i4) .eqv. .false. ) cycle
          call components(3,i4,occ4,ib(3),ie(3))
          call calcorr(4,aconf(ib(3):ie(3)))
          if ( break ) cycle
          call calcorr(5,aconf(ib(2):ie(3)))
          if ( break ) cycle
          call calcorr(6,aconf(ib(1):ie(3)))
          if ( break ) cycle
          occ5=.true.
          do i5=1,ncm(4)
            if ( occ5(i5) .eqv. .false. ) cycle
            call components(4,i5,occ5,ib(4),ie(4))
            call calcorr(7,aconf(ib(4):ie(4)))
            if ( break ) cycle
            call calcorr(8,aconf(ib(3):ie(4)))
            if ( break ) cycle
            call calcorr(9,aconf(ib(2):ie(4)))
            if ( break ) cycle
            call calcorr(10,aconf(ib(1):ie(4)))
            if ( break ) cycle
            call trialclose
            if ( break ) return
          end do
        end do
      end do
    end do
    deallocate(occ2)
  end do
  return
end subroutine quinary

subroutine senary
  do g=1,ng
    allocate(occ2(nct(g)))
    occ2=.true.
    do i2=1,nct(g)
      if ( occ2(i2) .eqv. .false. ) cycle
      call component1
      call calcorr(1,aconf(ib(1):ie(1)))
      if ( break ) cycle
      occ3=.true.
      do i3=1,ncm(2)
        if ( occ3(i3) .eqv. .false. ) cycle
        call components(2,i3,occ3,ib(2),ie(2))
        call calcorr(2,aconf(ib(2):ie(2)))
        if ( break ) cycle
        call calcorr(3,aconf(ib(1):ie(2)))
        if ( break ) cycle
        occ4=.true.
        do i4=1,ncm(3)
          if ( occ4(i4) .eqv. .false. ) cycle
          call components(3,i4,occ4,ib(3),ie(3))
          call calcorr(4,aconf(ib(3):ie(3)))
          if ( break ) cycle
          call calcorr(5,aconf(ib(2):ie(3)))
          if ( break ) cycle
          call calcorr(6,aconf(ib(1):ie(3)))
          if ( break ) cycle
          occ5=.true.
          do i5=1,ncm(4)
            if ( occ5(i5) .eqv. .false. ) cycle
            call components(4,i5,occ5,ib(4),ie(4))
            call calcorr(7,aconf(ib(4):ie(4)))
            if ( break ) cycle
            call calcorr(8,aconf(ib(3):ie(4)))
            if ( break ) cycle
            call calcorr(9,aconf(ib(2):ie(4)))
            if ( break ) cycle
            call calcorr(10,aconf(ib(1):ie(4)))
            if ( break ) cycle
            occ6=.true.
            do i6=1,ncm(5)
              if ( occ6(i6) .eqv. .false. ) cycle
              call components(5,i6,occ6,ib(5),ie(5))
              call calcorr(11,aconf(ib(5):ie(5)))
              if ( break ) cycle
              call calcorr(12,aconf(ib(4):ie(5)))
              if ( break ) cycle
              call calcorr(13,aconf(ib(3):ie(5)))
              if ( break ) cycle
              call calcorr(14,aconf(ib(2):ie(5)))
              if ( break ) cycle
              call calcorr(15,aconf(ib(1):ie(5)))
              if ( break ) cycle
              call trialclose
              if ( break ) return
            end do
          end do
        end do
      end do
    end do
    deallocate(occ2)
  end do
  return
end subroutine senary

subroutine component1
  valid(:,1)=.false.
  aconf(1)=0
  call int2cconf(i2,aconf(2:k(1)),nct(g),E,m(1)-group(g),k(1)-1_2,lconf(2:k(1))-group(g))
  aconf(1:k(1))=aconf(1:k(1))+group(g)
  do o=1,no
    oconf(1:k(1))=eqamat(aconf(1:k(1)),o)
    if ( all(oconf(1:k(1)) /= group(g)) ) cycle
    call sort(oconf(1:k(1)),1_2,k(1))
    call cconf2int(ic,oconf(2:k(1))-group(g),lconf(2:k(1))-group(g),nct(g),E,k(1)-1_2)
    if ( ic == i2 ) valid(o,1)=.true.
    if ( occ2(ic) .eqv. .true. ) then
      occ2(ic)=.false.
    end if
  end do
  return
end subroutine component1

subroutine components(n,im,occ,ib,ie)
  integer(1) :: n
  integer(8) :: im
  logical(1) :: occ(:)
  integer(2) :: ib,ie

  valid(:,n)=.false.
  call int2cconf(im,aconf(ib:ie),ncm(n),E,m(n),k(n),lconf(ib:ie))
  call cconf2aconf(aconf(ib:ie),aconf(1:ib-1),na)
  do o=1,no
    if ( valid(o,n-1) .eqv. .false. ) cycle
    oconf(ib:ie)=eqamat(aconf(ib:ie),o)
    call sort(oconf(ib:ie),1_2,k(n))
    call aconf2cconf(oconf(ib:ie),aconf(1:ib-1),na)
    call cconf2int(ic,oconf(ib:ie),lconf(ib:ie),ncm(n),E,k(n))
    if ( ic == im ) valid(o,n)=.true.
    if ( occ(ic) .eqv. .true. ) then
      occ(ic)=.false.
    end if
  end do
  return
end subroutine components

subroutine calcorr(n,conf)
  integer(1) :: n
  integer(2) :: conf(:)

  break=.false.
  call correction(corr(:,n),clus2,deg2,clus3,deg3,na,conf)
  dcorr=dabs(corr(:,n)-tcorr(:,n))
  if ( any(dcorr(1:npmc) > eps) .or. any(dcorr(npmc+1:) > tol) ) break=.true.
  return
end subroutine calcorr

subroutine trialclose
  do i=1,nsubc
    write(*,'('//fm//'F8.4)') corr(:,i)
  end do
  write(*,*)
  nirr=nirr+1
  if ( nirr == nsqs ) break=.true. 
  allocate(p%next)
  p=>p%next
  allocate(p%iconf(ns))
  p%iconf=aconf
  return
end subroutine trialclose

end subroutine irrconfig

subroutine confnum(ncm,nct,ib,ie,m,k,group,ng)
  implicit none
  integer(2) :: i,j,nk,ng,group(:)
  integer(2) :: k(:),m(:),ib(:),ie(:)
  integer(8) :: ncm(:),nct(:)

  nk=size(k)-1
  m(1)=sum(k)
  do i=1,nk-1
    m(i+1)=m(i)-k(i)
  end do
  do i=1,nk
    ncm(i)=nchoosek(m(i),k(i))
  end do
  j=m(1)-1
  do i=1,ng
    nct(i)=nchoosek(j,k(1)-1_2)
    j=j+group(i)-group(i+1)
  end do
  ib(1)=1
  ie(1)=k(1)
  do i=1,nk-1
    ib(i+1)=ib(i)+k(i)
    ie(i+1)=ie(i)+k(i+1)
  end do
  return
end subroutine confnum

subroutine matrixE(E,m,k)
  implicit none
  integer(8),allocatable :: E(:,:)
  integer(8),allocatable :: D(:,:)
  integer(2) :: k(:),i,j,m(:),ni,nj
  integer(1) :: nk

  nk=size(k)-1
  ni=maxval(m-k(1:nk))+1
  nj=maxval(k(1:nk))
  allocate(D(ni,nj))
  allocate(E(ni,nj))
  D(1,:)=0
  do i=2,ni
    do j=1,nj
      D(i,j)=nchoosek(nj-j+i-2_2,nj-j)
    end do
  end do
  do i=1,ni
    E(i,1)=sum(D(1:i,1))
  end do
  E(:,2:)=D(:,1:nj-1)
  return
end subroutine matrixE

subroutine lastconf(lconf,m,k)
  implicit none
  integer(2) :: lconf(:),m(:),k(:)
  integer(2) :: ib,j
  integer(1) :: i,nk

  nk=size(k)-1
  ib=0
  do i=1,nk
    forall(j=1:k(i)) lconf(ib+j)=m(i)-k(i)+j
    ib=ib+k(i)
  end do
end subroutine lastconf

subroutine int2cconf(ic,cconf,nc,E,m,k,lconf)
  implicit none
  integer(8) :: ic,nc,n
  integer(8),target :: E(:,:)
  integer(2) :: m,k,ni,nj,i
  integer(2) :: cconf(k),lconf(:)

  n=nc-ic
  ni=m-k+1
  nj=size(E,2)-k
  do i=1,k-1
    cconf(i)=binsearch_8(n,E(1:ni,nj+i))
    ni=cconf(i)
    n=n-E(cconf(i),nj+i)
  end do
  cconf(i)=n+1
  cconf=lconf-cconf+1
  return
end subroutine int2cconf

subroutine cconf2aconf(conf,fconf,na)
  implicit none
  integer(2) :: conf(:),fconf(:)
  integer(2) :: na
  integer(2) :: label(na)

  label=complement(fconf,na)
  conf=label(conf)
  return
end subroutine cconf2aconf

subroutine aconf2cconf(conf,fconf,na)
  implicit none
  integer(2) :: i,j,l,na
  integer(2) :: conf(:),fconf(:),c(na)

  c=1
  c(fconf)=0
  c(conf)=2
  j=0
  l=0
  do i=1,maxval(conf)
    j=j+c(i)
    if ( c(i) == 2 ) then
      l=l+1
      j=j-1
      conf(l)=j
    end if
  end do
  return
end subroutine aconf2cconf

subroutine cconf2int(ic,cconf,lconf,nc,E,k)
  implicit none
  integer(2) :: cconf(:),lconf(:),k
  integer(8) :: ic,nc,E(:,:)
  integer(2) :: n,j,conf(k)

  ic=0
  n=size(E,2)-k
  conf=lconf-cconf+1
  do j=1,k
    ic=ic+E(conf(j),n+j)
  end do
  ic=nc-ic
  return
end subroutine cconf2int

end module configurations
