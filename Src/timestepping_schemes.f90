!===============================================================================================
subroutine advance_fwdeuler(consvars_nm1,consvars,residual,vfraction,lo,hi,ng,dx,tstep,tfactor,problo,probhi) bind(C,name="advance_fwdeuler")

	implicit none
	include 'f_constants.h'
	
	integer :: lo(THREEDIM),hi(THREEDIM),ng
	double precision :: dx(THREEDIM),tstep,tfactor
	double precision :: problo(THREEDIM),probhi(THREEDIM)

	double precision :: consvars_nm1(lo(DIRX)-ng:hi(DIRX)+ng,   &
					 lo(DIRY)-ng:hi(DIRY)+ng,   &
					 lo(DIRZ)-ng:hi(DIRZ)+ng,NCVRS)

	double precision :: consvars(lo(DIRX)-ng:hi(DIRX)+ng,  &
			             lo(DIRY)-ng:hi(DIRY)+ng,  &
				     lo(DIRZ)-ng:hi(DIRZ)+ng,NCVRS)

	double precision :: residual(lo(DIRX)-ng:hi(DIRX)+ng,  &
			             lo(DIRY)-ng:hi(DIRY)+ng,  &
				     lo(DIRZ)-ng:hi(DIRZ)+ng,NCVRS)

	double precision :: vfraction(lo(DIRX)-ng:hi(DIRX)+ng,  &
			             lo(DIRY)-ng:hi(DIRY)+ng,  &
				     lo(DIRZ)-ng:hi(DIRZ)+ng)

	integer :: i,j,k,comp

	double precision :: vol_inv,vol

	vol = dx(DIRX)*dx(DIRY)*dx(DIRZ)
        vol_inv = 1.d0/vol
        
	!solving  V du/dt + R = 0

	do comp=1,NCVRS
	   do k=lo(DIRZ),hi(DIRZ)
	     do j=lo(DIRY),hi(DIRY)
	       do i=lo(DIRX),hi(DIRX)
			consvars(i,j,k,comp) = consvars_nm1(i,j,k,comp)    &
			- tfactor*residual(i,j,k,comp)*vol_inv*tstep*vfraction(i,j,k)
		enddo
	      enddo
	    enddo
	enddo

end subroutine advance_fwdeuler
!===============================================================================================
