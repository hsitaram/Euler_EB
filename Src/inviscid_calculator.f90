!===============================================================================================
subroutine update_primitive_vars(dens,velx,vely,velz,pres,consvar,lo,hi,ng,dx,prob_lo,prob_hi) bind(C,name="update_primitive_vars")
  
  use inviscidterms_module
  implicit none

  integer          :: lo(3), hi(3), ng

  double precision :: dens(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: velx(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: vely(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: velz(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: pres(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)

  double precision :: consvar(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,NCVRS)

  double precision :: dx(3) 
  double precision :: prob_lo(3) 
  double precision :: prob_hi(3)
  double precision :: primvars(5)

  integer :: i,j,k

  do k=lo(3),hi(3)
	do j=lo(2),hi(2)
		do i=lo(1),hi(1)

			call cons_to_primitive(consvar(i,j,k,:),GAMA,primvars)

			dens(i,j,k) = primvars(DENS_INDX)
			velx(i,j,k) = primvars(VELX_INDX)
			vely(i,j,k) = primvars(VELY_INDX)
			velz(i,j,k) = primvars(VELZ_INDX)
			pres(i,j,k) = primvars(PRES_INDX)

		enddo
	enddo
   enddo


end subroutine update_primitive_vars
!===============================================================================================
subroutine update_conservative_vars(dens,velx,vely,velz,pres,consvar,lo,hi,ng,dx,prob_lo,prob_hi) bind(C,name="update_conservative_vars")
  
  use inviscidterms_module
  implicit none

  integer          :: lo(3), hi(3), ng

  double precision :: dens(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: velx(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: vely(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: velz(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: pres(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)

  double precision :: consvar(lo(1)-ng :hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,NCVRS)

  double precision :: dx(3) 
  double precision :: prob_lo(3) 
  double precision :: prob_hi(3)
  double precision :: primvars(5)

  integer :: i,j,k

  do k=lo(3)-ng,hi(3)+ng
	do j=lo(2)-ng,hi(2)+ng
		do i=lo(1)-ng,hi(1)+ng

			primvars(DENS_INDX) = dens(i,j,k)
			primvars(VELX_INDX) = velx(i,j,k)
			primvars(VELY_INDX) = vely(i,j,k)
			primvars(VELZ_INDX) = velz(i,j,k)
			primvars(PRES_INDX) = pres(i,j,k)
			call prim_to_conservative(primvars,GAMA,consvar(i,j,k,:))
		enddo
	enddo
   enddo


end subroutine update_conservative_vars
!===============================================================================================
subroutine update_inviscidResidual(consvar,residual,lo,hi,ng,dx,prob_lo,prob_hi) bind(C,name="update_inviscidResidual")

  use inviscidterms_module
  implicit none

  integer          :: lo(3), hi(3), ng
  double precision :: consvar(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,NCVRS)
  double precision :: residual(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,NCVRS)
  double precision :: dx(3) 
  double precision :: prob_lo(3) 
  double precision :: prob_hi(3)

  integer :: i,j,k

  residual = 0.d0

  !x direction
  do k=lo(3),hi(3)
	do j=lo(2),hi(2)
		call oned_residualsweep(consvar(:,j,k,:),residual(:,j,k,:),lo(DIRX),hi(DIRX),ng,DIRX,dx(DIRY)*dx(DIRZ))
	enddo
  enddo
  
  !y direction
  do k=lo(3),hi(3)
	do i=lo(1),hi(1)
		call oned_residualsweep(consvar(i,:,k,:),residual(i,:,k,:),lo(DIRY),hi(DIRY),ng,DIRY,dx(DIRX)*dx(DIRZ))
	enddo
  enddo

  !z direction
  do j=lo(2),hi(2)
	do i=lo(1),hi(1)
		call oned_residualsweep(consvar(i,j,:,:),residual(i,j,:,:),lo(DIRZ),hi(DIRZ),ng,DIRZ,dx(DIRX)*dx(DIRY))
	enddo
  enddo
 
end subroutine update_inviscidResidual
!===============================================================================================
subroutine oned_residualsweep(consvars,res,lo,hi,ng,dir,area)

  use inviscidterms_module
  implicit none

  integer :: lo,hi,ng,dir
  double precision :: consvars(lo-ng:hi+ng,NCVRS)
  double precision :: res(lo-ng:hi+ng,NCVRS)

  double precision :: ul(NCVRS),ur(NCVRS)
  double precision :: fhalf(NCVRS)
  double precision :: area

  double precision :: normal(THREEDIM)
  integer :: i

  normal = 0.d0

  if(dir .eq. DIRX) normal(DIRX)=1.0
  if(dir .eq. DIRY) normal(DIRY)=1.0
  if(dir .eq. DIRZ) normal(DIRZ)=1.0

  do i=lo,hi+1
	
	ul=consvars(i-1,  :)
	ur=consvars(i,    :)

	!call centrdiff(ul,ur,GAMA,normal,fhalf)
	!call laxf(ul,ur,GAMA,normal,fhalf)
	call ausm(ul,ur,GAMA,normal,fhalf)

	res(i-1,:) = res(i-1,:) + fhalf(:)*area
	res(i  ,:) = res(i  ,:) - fhalf(:)*area
   enddo

	
	

end subroutine oned_residualsweep
!===============================================================================================
subroutine findcfltimestep(consvar,lo,hi,ng,dx,prob_lo,prob_hi,cfl,tstep,nanflag) bind(C,name="findcfltimestep")

    use inviscidterms_module
    implicit none

    integer :: lo(THREEDIM),hi(THREEDIM)
    integer :: ng
    double precision :: consvar(lo(DIRX)-ng:hi(DIRX)+ng,            &
	lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng,NCVRS)
    double precision :: dx(THREEDIM)
    integer :: prob_lo(THREEDIM),prob_hi(THREEDIM)
    integer :: nanflag

    double precision :: cfl
    double precision :: tstep

    double precision :: primvar(NCVRS),a,vel
  
    double precision :: dxmin,tstep_min,tstep_local;

    integer :: i,j,k

    dxmin = min(min(dx(DIRX),dx(DIRY)),dx(DIRZ))

    tstep_min = BIGVAL 
    nanflag = 0

    do k=lo(DIRZ),hi(DIRZ)
	do j=lo(DIRY),hi(DIRY)
	   do i=lo(DIRX),hi(DIRX)
		
		call cons_to_primitive(consvar(i,j,k,:),GAMA,primvar)

		vel = sqrt(primvar(VELX_INDX)**2+primvar(VELY_INDX)**2+primvar(VELZ_INDX)**2)
		a   = sqrt(GAMA*primvar(PRES_INDX)/primvar(DENS_INDX))

	        tstep_local = dxmin/(vel+a)

		if(isnan(tstep_local)) then
			print *,"tstep_local:",tstep_local
			print *,"pres,dens:",primvar(PRES_INDX),primvar(DENS_INDX)
			nanflag=1
			exit
		endif

		if(tstep_local .lt. tstep_min) then
			tstep_min = tstep_local
		endif

	     enddo
	   enddo
 	enddo

	tstep = cfl*tstep_min

end subroutine findcfltimestep 
!===============================================================================================
function isit_nan(x) result(isnan)

          real*8 :: x
          logical :: isnan

          isnan = .false.
          if(x .ne. x) then
              isnan = .true.
          endif

          if((x/x .ne. 1.0) .and. (x .ne. 0)) then
              isnan = .true.
          endif

  end function
!===================================================
