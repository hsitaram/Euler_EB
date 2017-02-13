program driver
  implicit none
  integer :: nnodes,ntri
  real*8, allocatable :: x(:,:)
  integer, allocatable :: conn(:,:),iblank(:,:,:)
  real*8 :: xlo(3),ds
  integer :: j,k,l
  integer :: jmax,kmax,lmax,holeCount
  character*1 :: hash
  character*40 :: arg
  !
  if (iargc() < 2) then
   write(6,*) 'Usage : holeCut <filename>'
   write(6,*) 'example : ./holeCut wedge.plt'
   stop
  endif
  call getarg(1,arg)
  open(unit=1,file=arg,form='formatted')
  read(1,*) hash,nnodes,ntri
  allocate(x(3,nnodes))
  allocate(conn(3,ntri))
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) x
  read(1,*) conn
  close(1)
  !
  xlo=(/-3.,-3.,-3./)
  ds=0.06
  jmax=101
  kmax=101
  lmax=101
  !
  allocate(iblank(jmax,kmax,lmax))
  iblank=1
  call holeCutSurf(jmax,kmax,lmax,iblank,xlo,ds,x,conn,0,nnodes,ntri,0,holeCount)
  write(6,*) 'holeCount=',holeCount
  !
  open(unit=1,file='check.p3d',form='unformatted')
  write(1) jmax,kmax,lmax
  write(1) (((xlo(1)+(j-1)*ds,j=1,jmax),k=1,kmax),l=1,lmax),&
          (((xlo(2)+(k-1)*ds,j=1,jmax),k=1,kmax),l=1,lmax),&
          (((xlo(3)+(l-1)*ds,j=1,jmax),k=1,kmax),l=1,lmax),&
          (((iblank(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
  close(1)
  deallocate(x,conn,iblank)
  !
end program driver
