module cluster
  use functions
  implicit none
contains

subroutine clusterprim(tcorr,pairs,trips,deg2,deg3,xc,a,x,rot,tra,cut,conc,symprec)
  implicit none
  real(8) :: a(3,3),x(:,:),tra(:,:),cut(2),symprec
  integer(1) :: rot(:,:,:)
  real(8) :: la,lb,lc,sinab,sinbc,sinca,r,ab,bc,ca,aa(2),bb(2),cc(2),xt(3)
  real(8),allocatable :: d(:),xc(:,:),xc1(:,:),dis2_t(:),dis2(:)
  real(8),allocatable :: dis3_t(:),dis3(:),tcorr(:,:)
  integer(1),allocatable :: xf1(:,:)
  integer(2),allocatable :: ip(:),eqamat(:,:,:),deg2_t(:),deg2(:)
  integer(2),allocatable :: deg3_t(:),deg3(:),trips_t(:,:),con(:)
  integer(2),allocatable :: pairs_t(:,:),pairs(:,:),trips(:,:)
  integer(1) :: ma,mb,mc,maxsz(3),minsz(3),lensz(3),o,no,xf(3),nk,nsubc
  integer(2) :: satom,na,i,j,k,m,n,ns,pair(2),trip(3),conc(:)
  integer(4) :: nc
  logical(1),allocatable :: occ2(:,:),occ3(:,:,:)

  la=norm(a(:,1))
  lb=norm(a(:,2))
  lc=norm(a(:,3))
  sinab=norm(cross(a(:,1),a(:,2)))/(la*lb)
  sinbc=norm(cross(a(:,2),a(:,3)))/(lb*lc)
  sinca=norm(cross(a(:,3),a(:,1)))/(lc*la)
  r=maxval(cut)
  ab=r/sinab
  bc=r/sinbc
  ca=r/sinca
  aa=(/ ab, ca /)
  bb=(/ ab, bc /)
  cc=(/ bc, ca /)
  ma=ceiling(maxval(aa)/la)
  mb=ceiling(maxval(bb)/lb)
  mc=ceiling(maxval(cc)/lc)
  satom=size(x,2)
  allocate(d(satom))
  allocate(xc1(3,satom*(2*ma+1)*(2*mb+1)*(2*mc+1)))
  xc1(:,1:satom)=x
  na=satom
  do i=-ma,ma+1
    do j=-mb,mb+1
      do k=-mc,mc+1
        xf=(/ i,j,k /)
        if ( all(xf == 0) ) cycle
        do m=1,satom
          xt=xf+x(:,m)
          do n=1,satom
            d(n)=distance(a,x(:,n),xt)
          end do
          if ( any(d < r+symprec) ) then
            na=na+1
            xc1(:,na)=xt
          end if
        end do
      end do
    end do
  end do
  allocate(xc(3,na))
  xc=xc1(:,1:na)

  allocate(xf1(3,na))
  xf1=floor(xc+symprec)
  maxsz=maxval(xf1,2)
  minsz=minval(xf1,2)
  lensz=maxsz-minsz+1
  ns=1_2*lensz(1)*lensz(2)*lensz(3)
  allocate(ip(ns*satom))
  do i=1,na
    n=pos2int(x,xc(:,i),xf1(:,i),minsz,lensz,ns,symprec)
    ip(n)=i
  end do
  no=size(rot,3)
  allocate(eqamat(na,na,no))
  deallocate(xc1)
  allocate(xc1(3,na))
  do i=1,no
    do j=1,na
      xc1(:,j)=matmul(rot(:,:,i),xc(:,j))+tra(:,i)
    end do
    do j=1,na
      xf=floor(xc1(:,j)+symprec)
      xc1=xc1-spread(xf,2,na)
      do m=1,na
        xf=floor(xc1(:,m)+symprec)
        if ( any(xf-minsz < 0) .or. any(xf-minsz > lensz-1) ) cycle
        eqamat(m,j,i)=ip(pos2int(x,xc1(:,m),xf,minsz,lensz,ns,symprec))
      end do
    end do
  end do

  allocate(occ2(na,satom))
  allocate(dis2_t(na*satom))
  allocate(deg2_t(na*satom))
  allocate(pairs_t(2,na*satom))
  deallocate(d)
  allocate(d(1))
  m=0
  n=0
  occ2=.false.
  do i=1,satom
    do j=1,na
      if ( occ2(j,i) ) cycle
      d=distance(a,xc(:,i),xc(:,j))
      if ( d(1) > cut(1)+symprec .or. d(1) < symprec ) cycle
      n=n+1
      dis2_t(n)=d(1)
      deg2_t(n)=0
      do o=1,no
        pair=(/ i, j /)
        pair=eqamat(pair,i,o)
        if ( occ2(pair(2),pair(1)) == .false. ) then
          m=m+1
          deg2_t(n)=deg2_t(n)+1
          pairs_t(:,m)=pair
          occ2(pair(2),pair(1))=.true.
        end if
      end do
      do o=1,no
        pair=(/ j, i /)
        pair=eqamat(pair,j,o)
        if ( occ2(pair(2),pair(1)) == .false. ) then
          m=m+1
          deg2_t(n)=deg2_t(n)+1
          pairs_t(:,m)=pair
          occ2(pair(2),pair(1))=.true.
        end if
      end do
    end do
  end do
  allocate(dis2(n))
  allocate(deg2(n))
  allocate(pairs(2,m))
  dis2=dis2_t(1:n)
  deg2=deg2_t(1:n)
  pairs=pairs_t(:,1:m)

  nc=nchoosek(na,2_2)
  allocate(occ3(na,na,satom))
  allocate(dis3_t(nc*satom))
  allocate(deg3_t(nc*satom))
  allocate(trips_t(3,nc*satom))
  deallocate(d)
  allocate(d(3))
  m=0
  n=0
  occ3=.false.
  do i=1,satom
    do j=1,na-1
      do k=j+1,na
        if ( occ3(k,j,i) ) cycle
        d(1)=distance(a,xc(:,i),xc(:,j))
        d(2)=distance(a,xc(:,j),xc(:,k))
        d(3)=distance(a,xc(:,k),xc(:,i))
        if ( any(d > cut(2)+symprec) .or. any(d < symprec) ) cycle
        n=n+1
        dis3_t(n)=maxval(d)
        deg3_t(n)=0
        do o=1,no
          trip=(/ i, j, k /)
          trip=eqamat(trip,i,o)
          call sort(trip(2:3),1_2,2_2)
          if ( occ3(trip(3),trip(2),trip(1)) == .false. ) then
            m=m+1
            deg3_t(n)=deg3_t(n)+1
            trips_t(:,m)=trip
            occ3(trip(3),trip(2),trip(1))=.true.
          end if
        end do
        do o=1,no
          trip=(/ j, i, k /)
          trip=eqamat(trip,j,o)
          call sort(trip(2:3),1_2,2_2)
          if ( occ3(trip(3),trip(2),trip(1)) == .false. ) then
            m=m+1
            deg3_t(n)=deg3_t(n)+1
            trips_t(:,m)=trip
            occ3(trip(3),trip(2),trip(1))=.true.
          end if
        end do
        do o=1,no
          trip=(/ k, i, j /)
          trip=eqamat(trip,k,o)
          call sort(trip(2:3),1_2,2_2)
          if ( occ3(trip(3),trip(2),trip(1)) == .false. ) then
            m=m+1
            deg3_t(n)=deg3_t(n)+1
            trips_t(:,m)=trip
            occ3(trip(3),trip(2),trip(1))=.true.
          end if
        end do
      end do
    end do
  end do
  allocate(dis3(n))
  allocate(deg3(n))
  allocate(trips(3,m))
  dis3=dis3_t(1:n)
  deg3=deg3_t(1:n)
  trips=trips_t(:,1:m)

  call sortcluster(dis2,deg2,pairs)
  call sortcluster(dis3,deg3,trips)

  m=size(deg2)
  n=size(deg3)
  na=sum(conc)
  nk=size(conc)
  nsubc=nk*(nk-1)/2 !2**(nk-1)-1
  allocate(con(nsubc))
  allocate(tcorr(m+n,nsubc))

  select case(nk)
  case(2)
    con=conc(1)
  case(3)
    con(1:2)=conc(1:2)
    con(3)=sum(conc(1:2))
  case(4)
    con(1:2)=conc(1:2)
    con(3)=sum(conc(1:2))
    con(4)=conc(3)
    con(5)=sum(conc(2:3))
    con(6)=sum(conc(1:3))
  case(5)
    con(1:2)=conc(1:2)
    con(3)=sum(conc(1:2))
    con(4)=conc(3)
    con(5)=sum(conc(2:3))
    con(6)=sum(conc(1:3))
    con(7)=conc(4)
    con(8)=sum(conc(3:4))
    con(9)=sum(conc(2:4))
    con(10)=sum(conc(1:4))
  case(6)
    con(1:2)=conc(1:2)
    con(3)=sum(conc(1:2))
    con(4)=conc(3)
    con(5)=sum(conc(2:3))
    con(6)=sum(conc(1:3))
    con(7)=conc(4)
    con(8)=sum(conc(3:4))
    con(9)=sum(conc(2:4))
    con(10)=sum(conc(1:4))
    con(11)=conc(5)
    con(12)=sum(conc(4:5))
    con(13)=sum(conc(3:5))
    con(14)=sum(conc(2:5))
    con(15)=sum(conc(1:5))
  end select

  do i=1,nsubc
    do j=1,m
      tcorr(j,i)=(1-2.D0*con(i)/na)**2
    end do
    do j=m+1,m+n
      tcorr(j,i)=(1-2.D0*con(i)/na)**3
    end do
  end do
  return
contains

subroutine sortcluster(dis,deg,tuple)
  implicit none
  real(8) :: dis(:)
  integer(2) :: n,i,m1,n1,m2,n2
  integer(2) :: deg(:),tuple(:,:)
  integer(2),allocatable :: tuple_t(:,:),id(:)

  n=size(dis)
  allocate(id(n))
  call sortindices(dis,id)
  allocate(tuple_t(size(tuple,1),size(tuple,2)))
  n1=0
  do i=1,n
    m1=n1+1
    n1=m1+deg(id(i))-1
    m2=sum(deg(1:id(i)-1))+1
    n2=m2+deg(id(i))-1
    tuple_t(:,m1:n1)=tuple(:,m2:n2)
  end do
  tuple=tuple_t
  deg=deg(id)  
  return
end subroutine sortcluster

function distance(a,x0,x1)
  real(8) :: distance,a(3,3),x0(3),x1(3)

  distance=norm(matmul(a,x1-x0))
  return
end function distance

function pos2int(x,xc,xf,minsz,lensz,ns,symprec)
  real(8) :: x(:,:),xc(3),xc1(3),symprec
  integer(1) :: xf(3),xt(3),minsz(3),lensz(3),m
  integer(2) :: pos2int,ns,n

  xt=xf-minsz
  n=1_2*xt(1)*lensz(2)*lensz(3)+1_2*xt(2)*lensz(3)+xt(3)+1
  xc1=xc-xf
  do m=1,size(x,2)
    if ( all(dabs(xc1-x(:,m)) < symprec) ) exit
  end do
  pos2int=(m-1)*ns+n
  return
end function pos2int

end subroutine clusterprim

subroutine clustersuper(clus2,clus3,a,x,xc0,M,tra,pairs,trips,symprec)
  implicit none
  real(8),allocatable :: xc(:,:)
  integer(2),allocatable :: clus2(:,:,:),clus3(:,:,:)
  integer(2) :: pairs(:,:),trips(:,:)
  real(8) :: a(3,3),x(:,:),xc0(:,:),tra(:,:),xt(3),symprec
  integer(2),allocatable :: eqamat(:)
  integer(1) :: nt,i,M(3,3)
  integer(2) :: na,nx,np,nq,n,j

  nt=size(tra,2)
  na=size(xc0,2)
  nx=size(x,2)
  np=size(pairs,2)
  nq=size(trips,2)
  allocate(xc(3,na))
  allocate(eqamat(na))
  allocate(clus2(2,np,nt))
  allocate(clus3(3,nq,nt))
  xc=matmul(inverse(dble(M)),xc0)
  do i=1,nt
    do j=1,na
      xt=xc(:,j)+tra(:,i)
      xt=xt-floor(xt+symprec)
      do n=1,nx
        if ( all(dabs(matmul(a,xt-x(:,n))) < symprec) ) then
          eqamat(j)=n
          exit
        end if
      end do
    end do
    do j=1,np
      clus2(:,j,i)=eqamat(pairs(:,j))
    end do
    do j=1,nq
      clus3(:,j,i)=eqamat(trips(:,j))
    end do
  end do
  return
end subroutine clustersuper

subroutine correction(corr,clus2,deg2,clus3,deg3,na,conf)
  implicit none
  real(8) :: corr(:),corrt
  integer(2) :: clus2(:,:,:),clus3(:,:,:),deg2(:),deg3(:)
  integer(2) :: na,conf(:),np,nq,m,n,j,k
  integer(1) :: nt,i,sigma(na)

  sigma=1
  sigma(conf)=-1
  np=size(deg2)
  nq=size(deg3)
  nt=size(clus2,3)
  do i=1,np
    corrt=0D0
    m=sum(deg2(1:i-1))+1
    n=sum(deg2(1:i))
    do j=1,nt
      do k=m,n
        corrt=corrt+sigma(clus2(1,k,j))*sigma(clus2(2,k,j))
      end do
    end do
    corr(i)=corrt/deg2(i)
  end do

  do i=1,nq
    corrt=0D0
    m=sum(deg3(1:i-1))+1
    n=sum(deg3(1:i))
    do j=1,nt
      do k=m,n
        corrt=corrt+sigma(clus3(1,k,j))*sigma(clus3(2,k,j))*sigma(clus3(3,k,j))
      end do
    end do
    corr(np+i)=corrt/deg3(i)
  end do
  corr=corr/nt
  return
end subroutine correction

end module cluster
