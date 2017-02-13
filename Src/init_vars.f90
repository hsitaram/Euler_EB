!===============================================================================================
subroutine initialize_sod(dens,velx,vely,velz,pres,lo,hi,ng,dx,prob_lo,prob_hi) bind(C, name="initialize_sod")

  implicit none
  include 'f_constants.h'

  integer          :: lo(THREEDIM), hi(THREEDIM), ng
 
  double precision :: dens(lo(DIRX)-ng:hi(DIRX)+ng,lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)
  double precision :: velx(lo(DIRX)-ng:hi(DIRX)+ng,lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)
  double precision :: vely(lo(DIRX)-ng:hi(DIRX)+ng,lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)
  double precision :: velz(lo(DIRX)-ng:hi(DIRX)+ng,lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)
  double precision :: pres(lo(DIRX)-ng:hi(DIRX)+ng,lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)

  double precision :: dx(THREEDIM) 
  double precision :: prob_lo(THREEDIM) 
  double precision :: prob_hi(THREEDIM) 

  integer          :: i,j,k
  double precision :: x,y,z,xlen
  
  xlen = prob_hi(DIRX)-prob_lo(DIRX)

  velx = 0.d0
  vely = 0.d0
  velz = 0.d0

  do k = lo(DIRZ), hi(DIRZ)
     z = prob_lo(DIRZ) + (dble(k)+0.5d0) * dx(DIRZ)

     do j = lo(DIRY), hi(DIRY)
     	y = prob_lo(DIRY) + (dble(j)+0.5d0) * dx(DIRY)

        do i = lo(DIRX), hi(DIRX)
           x = prob_lo(DIRX) + (dble(i)+0.5d0) * dx(DIRX)

	   if(x < 0.5*xlen) then
           	dens(i,j,k) = 1.d0
		pres(i,j,k) = 1.d0
	   else
		dens(i,j,k) = 0.125
		pres(i,j,k) = 0.1
	   endif

        end do
     end do
  end do

end subroutine initialize_sod
!===============================================================================================
subroutine init_blanking_field(phi,lo,hi,ng,dx,prob_lo,prob_hi,findices,sindices,nfindices,nsindices) bind(C,name="init_blanking_field")

implicit none

integer :: lo(3),hi(3),ng
double precision :: dx(3),prob_lo(3),prob_hi(3)
double precision :: phi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
integer :: nfindices,nsindices
integer :: findices(nfindices),sindices(nsindices)

integer :: i,j,k
integer :: f_ind,c_ind  !fortran and c index :)

phi(:,:,:)=1.d0

do f_ind=1,nsindices

	c_ind=f_ind-1
	i = sindices(3*c_ind+1)
	j = sindices(3*c_ind+2)
	k = sindices(3*c_ind+3)

	phi(i,j,k) = 0.d0

enddo
	
end subroutine init_blanking_field
!===============================================================================================
subroutine print_var(phi, lo, hi, ng, dx, prob_lo, prob_hi) bind(C, name="print_var")

  implicit none

  integer          :: lo(3), hi(3), ng
  double precision :: phi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: dx(3) 
  double precision :: prob_lo(3) 
  double precision :: prob_hi(3) 

  integer          :: i,j,k
  double precision :: x,y,z,r2

  do k = lo(3)-ng, hi(3)+ng
     do j = lo(2)-ng, hi(2)+ng
        do i = lo(1)-ng, hi(1)+ng
		print *,"i,j,k,phi:",i,j,k,phi(i,j,k)
        end do
     end do
  end do

end subroutine print_var
!==============================================================================================
