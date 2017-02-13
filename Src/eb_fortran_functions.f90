!==================================================================
subroutine spherefindindices(lo, hi,dx, prob_lo, prob_hi,rad,pos,fldindices,sldindices, &
nindmax,nfldindices,nsldindices) bind(C, name="spherefindindices")

implicit none

integer :: lo(3),hi(3)
integer :: nindmax
double precision :: dx(3),prob_lo(3),prob_hi(3)
double precision :: rad,pos(3)

integer   :: fldindices(nindmax)
integer   :: sldindices(nindmax)
integer   :: nfldindices,nsldindices

integer :: i,j,k
double precision :: x,y,z,r2

nfldindices=0
nsldindices=0


 do k = lo(3), hi(3)
     z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)

     do j = lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)

        do i = lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

           r2 = ((x-pos(1))**2 + (y-pos(2))**2 + (z-pos(3))**2)

	   if(r2 < rad**2) then
		sldindices(3*nsldindices+1) = i
		sldindices(3*nsldindices+2) = j
		sldindices(3*nsldindices+3) = k
	        nsldindices=nsldindices+1
	   else
		fldindices(3*nfldindices+1) = i
		fldindices(3*nfldindices+2) = j
		fldindices(3*nfldindices+3) = k
	        nfldindices=nfldindices+1
	   endif
		

        end do
     end do
  end do

end subroutine spherefindindices
!==================================================================
subroutine find_volfrac_in_box_forsphere(vfrac,lo,hi,ng,dx,problo,probhi,rad,pos) bind(C,name="find_volfrac_in_box_forsphere")

	implicit none
	include 'f_constants.h'

	integer          :: lo(THREEDIM),hi(THREEDIM),ng
	double precision :: dx(THREEDIM),pos(THREEDIM)
	double precision :: problo(THREEDIM),probhi(THREEDIM)	
	double precision :: vfrac(lo(DIRX)-ng:hi(DIRX)+ng,&
			    lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)
	double precision :: rad

	
	integer :: i,j,k
	double precision :: x,y,z,r2


	vfrac = 1.d0

 	do k = lo(DIRZ), hi(DIRZ)
     		z = problo(DIRZ) + (dble(k)+0.5d0) * dx(DIRZ)

     		do j = lo(DIRY), hi(DIRY)
        		y = problo(DIRY) + (dble(j)+0.5d0) * dx(DIRY)

        		do i = lo(DIRX), hi(DIRX)
           			x = problo(DIRX) + (dble(i)+0.5d0) * dx(DIRX)

           			r2 = ((x-pos(DIRX))**2 + (y-pos(DIRY))**2 + (z-pos(DIRZ))**2)

	   			if(r2 < rad**2) then
					vfrac(i,j,k)=0.d0
				endif

			enddo
		enddo
	enddo


end subroutine find_volfrac_in_box_forsphere
!==================================================================
subroutine find_volfrac_in_box(vfrac,lo,hi,ng,dx,problo,probhi,surfcoord,conn,nnodes,ntri) bind(C,name="find_volfrac_in_box")

	implicit none
	include 'f_constants.h'

	integer          :: lo(THREEDIM),hi(THREEDIM),ng
	double precision :: dx(THREEDIM),pos(THREEDIM)
	double precision :: problo(THREEDIM),probhi(THREEDIM)
	
	double precision :: vfrac(lo(DIRX)-ng:hi(DIRX)+ng,&
			    lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)
  	integer :: nnodes,ntri
  	double precision :: surfcoord(THREEDIM,nnodes)
  	integer :: conn(THREEDIM,ntri)
	integer :: i,j,k,solidpoints

	vfrac = 1.d0
  	call holeCutSurf(lo,hi,problo,dx,ng,surfcoord,conn,nnodes,ntri,solidpoints,vfrac)

end subroutine find_volfrac_in_box
!==================================================================
subroutine read_triangulation_file(surfcoord,conn,nnodes,ntri) bind(C,name="read_triangulation_file")

	use paraviewfilewrite
	implicit none
	include 'f_constants.h'

	integer :: nnodes,ntri
	double precision :: surfcoord(THREEDIM,nnodes)
	integer :: conn(THREEDIM,ntri)
  	character*1 :: hash
	integer :: fptr

	fptr=31
  	open(unit=fptr,file="surface_triangulation.plt",form='formatted')

  	read(fptr,*) 
  	read(fptr,*)
  	read(fptr,*)
  	read(fptr,*)
  	read(fptr,*) surfcoord
  	read(fptr,*) conn

  	close(fptr)

	!write vtu file
  	call Writeumeshfile("surface.vtu",surfcoord,conn,nnodes,ntri)

end subroutine read_triangulation_file
!==================================================================
