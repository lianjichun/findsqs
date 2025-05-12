module supercell
  use functions
  implicit none
contains

subroutine buildsupercell(a,x,as,xs,M)
  implicit none
  real(8) :: a(3,3),as(3,3),x(:,:)
  real(8),allocatable :: xt(:,:),xs(:,:)
  real(8),allocatable :: T(:,:),Pd(:,:)
  integer(1) :: M(3,3),bd(2,3),vertex(3,8)
  integer(1),allocatable :: Pc(:,:)
  integer(2) :: i,j,k,n,nt,np,na

  nt=abs(det(M))
  allocate(T(3,nt))
  vertex(:,1)=0
  vertex(:,2:4)=M
  vertex(:,5)=M(:,1)+M(:,2)
  vertex(:,6)=M(:,1)+M(:,3)
  vertex(:,7)=M(:,2)+M(:,3)
  vertex(:,8)=M(:,1)+M(:,2)+M(:,3)
  do i=1,3
    bd(1,i)=maxval(vertex(i,:))
    if ( bd(1,i) == 0 ) bd(1,i)=bd(1,i)+1
    bd(2,i)=minval(vertex(i,:))
  end do
  np=1_2*(bd(1,1)-bd(2,1))*(bd(1,2)-bd(2,2))*(bd(1,3)-bd(2,3))
  allocate(Pc(3,np))
  allocate(Pd(3,np))
  n=0
  do i=bd(2,1),bd(1,1)-1
    do j=bd(2,2),bd(1,2)-1
      do k=bd(2,3),bd(1,3)-1
        n=n+1
        Pc(:,n)=(/ i, j, k /)
      end do
    end do
  end do
  Pd=matmul(inverse(dble(M)),dble(Pc))
  n=0
  do i=1,np
    if ( all(Pd(:,i) > -0.01) .and. all(Pd(:,i) < 0.99) ) then
      n=n+1
      T(:,n)=Pd(:,i)
    end if
  end do
  na=size(x,2)
  allocate(xt(3,na))
  xt=matmul(a,x)
  as=matmul(a,dble(M))
  xt=matmul(inverse(as),xt)
  n=0
  allocate(xs(3,na*nt))
  do i=1,na
    do j=1,nt
      n=n+1
      xs(:,n)=xt(:,i)+T(:,j)
    end do
  end do
  xs=xs-floor(xs+1D-5)
  return
end subroutine buildsupercell

subroutine sqscell(a,nlat,rot,n)
  implicit none
  real(8) :: a(3,3)
  integer(2) :: nlat,i,nm
  integer(1) :: rot(:,:,:),n
  real(8),allocatable :: lens(:)  
  integer(1),allocatable :: M(:,:,:)
  integer(2),allocatable :: id(:)

  call genlattices(rot,n,M)
  nm=size(M,3)
  allocate(lens(nm))
  do i=1,nm
    call delaunary(a,M(:,:,i),lens(i))
  end do
  allocate(id(nm))
  call sortindices(lens,id)
  M=M(:,:,id)
  if ( nlat < 1 .or. nlat > nm ) nlat=nm
  open(11,file='SQSCELL')
  write(11,'(I6)') nlat
  do i=1,nlat
    write(11,*)
    write(11,'(3I4)') M(:,1,i)
    write(11,'(3I4)') M(:,2,i)
    write(11,'(3I4)') M(:,3,i)
  end do
  close(11)
  stop
end subroutine sqscell

subroutine genlattices(rot,n,M)
  implicit none
  integer(1) :: rot(:,:,:),n
  integer(1),allocatable :: fac(:),pow(:)
  integer(1),allocatable :: Mi(:,:,:),M(:,:,:)
  integer(1) :: Mt(3,3),Bi(3,3)
  integer(1) :: b,c,d,e,f,g,o,nf,nb,no,inpnum
  integer(2) :: ns,i,j,k
  real(8),allocatable :: R(:,:,:),Minv(:,:,:)
  real(8) :: RB(3,3),H(3,3),dH(3,3)
  logical(1),allocatable :: occ(:)

  inpnum=n
  call prime(fac,pow,inpnum)
  nf=size(fac)
  ns=1
  do i=1,nf
    ns=ns*(((fac(i)**(pow(i)+2)-1)*(fac(i)**(pow(i)+1)-1))/((fac(i)-1)**2*(fac(i)+1)))
  end do
  allocate(M(3,3,ns))
  i=0
  do b=1,n
    if ( mod(n,b) /= 0 ) cycle
    nb=n/b
    do d=1,nb
      if ( mod(nb,d) /= 0 ) cycle
      g=nb/d
      do c=0,d-1
        do e=0,g-1
          do f=0,g-1
            Mt(1,:)=(/ b, 0_1, 0_1 /)
            Mt(2,:)=(/ c,  d,  0_1 /)
            Mt(3,:)=(/ e,  f,   g  /)
            if ( det(Mt) == n ) then
              i=i+1
              M(:,:,i)=Mt
            end if
          end do
        end do
      end do
    end do
  end do
  no=size(rot,3)
  allocate(R(3,3,no))
  do i=1,no
    R(:,:,i)=inverse(dble(rot(:,:,i)))
  end do
  allocate(Minv(3,3,ns))
  do i=1,ns
    Minv(:,:,i)=inverse(dble(M(:,:,i)))
  end do
  allocate(occ(ns))
  allocate(Mi(3,3,ns))
  occ=.false.
  k=0
  do i=1,ns
    if ( occ(i) ) cycle
    Bi=M(:,:,i)
    k=k+1
    Mi(:,:,k)=Bi
    do o=1,no
      RB=matmul(R(:,:,o),dble(Bi))
      do j=i+1,ns
        if ( occ(j) ) cycle
        H=matmul(Minv(:,:,j),RB)
        dH=dabs(H-nint(H))
        if ( all(dH < 1D-5) ) occ(j)=.true.
      end do
    end do
  end do
  deallocate(M)
  allocate(M(3,3,k))
  M=Mi(:,:,1:k)
  return
end subroutine genlattices

subroutine prime(factors,powers,inpnum)
  implicit none
  integer(1),allocatable :: factors(:),powers(:)
  integer(1) :: factors_t(3),powers_t(3)
  integer(1) :: inpnum,uplimit,fac,pow,n

  uplimit=sqrt(real(inpnum))
  n=0
  do fac=2,uplimit
    if ( isprime(fac) ) then
      pow=0
      do while ( mod(inpnum,fac) == 0 )
        pow=pow+1
        inpnum=inpnum/fac
      end do
      if ( pow /= 0 ) then
        n=n+1
        factors_t(n)=fac
        powers_t(n)=pow
      end if
    end if
  end do
  if ( inpnum /= 1 ) then
    n=n+1
    factors_t(n)=inpnum
    powers_t(n)=1
  end if
  allocate(factors(n))
  allocate(powers(n))
  factors=factors_t(1:n)
  powers=powers_t(1:n)
  return
end subroutine prime

function isprime(num)
  implicit none
  logical(1) :: isprime
  integer(1) :: num,sqr,i

  if ( num == 1 ) then
    isprime=.false.
    return
  end if
  sqr=sqrt(real(num))
  do i=2,sqr
    if ( mod(num,i) == 0 ) then
      isprime=.false.
      return
    end if
  end do
  isprime=.true.
  return
end function isprime

subroutine delaunary(a,M,len2)
  implicit none
  real(8) :: a(3,3),as(3,3),volume,vol
  real(8) :: v(3,7),b(3,4),g(6),s(7),d,len2
  real(8),allocatable :: lens(:)
  integer(1) :: M(3,3),vec(3,35),i,j,k,k1,k2,n,l(1)

  as=matmul(a,dble(M))
  volume=dot_product(cross(as(:,1),as(:,2)),as(:,3))
  b(:,1:3)=as
  b(:,4)=-(as(:,1)+as(:,2)+as(:,3))
  do while ( .true. )
    g(1)=dot_product(b(:,1),b(:,2))
    g(2)=dot_product(b(:,1),b(:,3))
    g(3)=dot_product(b(:,1),b(:,4))
    g(4)=dot_product(b(:,2),b(:,3))
    g(5)=dot_product(b(:,2),b(:,4))
    g(6)=dot_product(b(:,3),b(:,4))
    d=maxval(g)
    l=maxloc(g)
    if ( d < 1D-6 ) then
      exit
    end if
    select case (l(1))
      case (1)
        k1=1;k2=2
      case (2)
        k1=1;k2=3
      case (3)
        k1=1;k2=4
      case (4)
        k1=2;k2=3
      case (5)
        k1=2;k2=4
      case (6)
        k1=3;k2=4
    end select
    b(:,k1)=-b(:,k1)
    do i=1,4
      if ( i /= k1 .and. i /= k2 ) then
        b(:,i)=b(:,i)-b(:,k1)
      end if
    end do
  end do
  v(:,1:4)=b
  v(:,5)=b(:,1)+b(:,2)
  v(:,6)=b(:,2)+b(:,3)
  v(:,7)=b(:,3)+b(:,1)
  do i=1,7
    s(i)=norm(v(:,i))
  end do
  n=0
  do i=1,5
    do j=i+1,6
      do k=j+1,7
        vol=dot_product(cross(v(:,i),v(:,j)),v(:,k))
        if ( dabs(dabs(vol)-volume) < 1D-5 ) then
          n=n+1
          if ( vol > 0 ) then
            vec(:,n)=(/ i, j, k /)
          else
            vec(:,n)=(/ i, k, j /)
          end if
        end if
      end do
    end do
  end do
  allocate(lens(n))
  do i=1,n
    lens(i)=s(vec(1,i))+s(vec(2,i))+s(vec(3,i))
  end do
  l=minloc(lens)
  n=l(1)
  as(:,1)=v(:,vec(1,n))
  as(:,2)=v(:,vec(2,n))
  as(:,3)=v(:,vec(3,n))
  len2=norm(as(:,1))**2+norm(as(:,2))**2+norm(as(:,3))**2+norm(as(:,1)+as(:,2)+as(:,3))**2
  M=nint(matmul(inverse(a),as))
  return
end subroutine delaunary

end module supercell
