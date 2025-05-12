program findsqs
  use functions
  use structure
  use symmetry
  use supercell
  use cluster
  use configurations
  use outfiles
  implicit none
  real(8) :: a(3,3),as(3,3),cuts(2),prec,tole
  real(8),allocatable :: x(:,:),xs(:,:),xc(:,:)
  real(8),allocatable :: tcorr(:,:),tra(:,:)
  integer(1) :: M(3,3),cell,comp,site,npmc
  integer(2) :: nlat,nsqs,i,conc(6)
  integer(1),allocatable :: rot(:,:,:),lat(:,:,:)
  integer(2),allocatable :: deg3(:),deg2(:),natom(:)
  integer(2),allocatable :: pairs(:,:),trips(:,:)
  integer(2),allocatable :: clus2(:,:,:),clus3(:,:,:)
  integer(2),allocatable :: k(:),group(:),eqamat(:,:)
  integer(2),allocatable :: iconf(:,:)
  integer(2),allocatable :: natoms(:)
  real(4) :: time0,time1
  integer(4) :: time
  integer(1) :: hour,mins,secs
  character(len=2) :: symb(6)
  character(len=2),allocatable :: symbol(:),symbs(:)
  character(len=30) :: fm,tail,poscar
  logical(1) :: alive
  namelist /insqs/ comp,cell,conc,nlat,nsqs,prec,symb,tole,cuts,npmc

  call cpu_time(time0)

  open(11,file='INSQS')
  read(11,nml=insqs)
  close(11)
  allocate(k(comp))
  allocate(symbs(comp))
  k=conc(1:comp)
  symbs=symb(1:comp)
  call readpos(a,x,symbol,natom,'POSCAR')
  call rotation(eqamat,rot,tra,a,x,natom,prec)
  inquire(file='SQSCELL',exist=alive)
  if ( .not. alive ) call sqscell(a,nlat,rot,cell)
  call clusterprim(tcorr,pairs,trips,deg2,deg3,xc,a,x,rot,tra,cuts,k,prec)
  natom=natom*cell
  allocate(group(2))
  group(1)=1
  group(2)=sum(natom)+1
  site=1

  call system('rm -r poscar-* 2> /dev/null')
  open(13,file='SQSCELL')
  read(13,*) nlat
  do i=1,nlat
    read(13,*);read(13,*) M
    deallocate(eqamat,tra)
    call buildsupercell(a,x,as,xs,M)
    call eqamatrix(eqamat,tra,as,xs,natom,prec,site)
    call clustersuper(clus2,clus3,as,xs,xc,M,tra,pairs,trips,prec)
    call irrconfig(iconf,eqamat,group,k,clus2,deg2,clus3,deg3,tole,npmc,nsqs,tcorr)
    if (size(iconf,2) /= 0 ) then
      call outposcar(as,xs,symbol,natom,symbs,k,site,iconf,nlat,i)
    end if 
    deallocate(clus2,clus3,xs,iconf)
    write(*,*) i
  end do
  close(13)

  call cpu_time(time1)
  time=nint(time1-time0)
  hour=time/3600
  mins=(time-hour*3600)/60
  secs=time-hour*3600-mins*60
  write(*,'(I2.2,A,I2.2,A,I2.2,A,/)') hour,' hour ',mins,' min ',secs,' sec'
end program findsqs
