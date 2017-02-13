module inviscidterms_module

implicit none
	include "f_constants.h"

contains
!=====================================================
subroutine laxf(ul,ur,g,n,flux)

	implicit none	
	double precision :: g
	double precision :: n(THREEDIM)
	double precision :: ul(NCVRS),ur(NCVRS)
	double precision :: pl(NCVRS),pr(NCVRS)
	double precision :: fl(NCVRS),fr(NCVRS)
	double precision :: flux(NCVRS)

	double precision :: a,rho,p,vn,nx,ny,nz

	nx = n(DIRX)
	ny = n(DIRY)
	nz = n(DIRZ)

	call cons_to_primitive(ul,g,pl)
	call cons_to_primitive(ur,g,pr)

	call getflux(ul,g,n,fl) 
	call getflux(ur,g,n,fr) 
	
	rho = 0.5*(pl(DENS_INDX)+pr(DENS_INDX))
	p   = 0.5*(pl(PRES_INDX)+pr(PRES_INDX))

	vn = 0.5*(pl(VELX_INDX)+pr(VELX_INDX))*nx
	vn = vn + 0.5*(pl(VELY_INDX)+pr(VELY_INDX))*ny
	vn = vn + 0.5*(pl(VELZ_INDX)+pr(VELZ_INDX))*nz

	a = sqrt(g*p/rho)

	flux(:) = 0.5*(fl(:)+fr(:)) - 0.5*(abs(vn)+a)*(ur(:)-ul(:))

end subroutine laxf
!=====================================================
subroutine ausm(ul,ur,g,n,flux)
	
	implicit none	
	double precision :: g
	double precision :: n(THREEDIM)
	double precision :: ul(NCVRS),ur(NCVRS)
	double precision :: pl(NCVRS),pr(NCVRS)
	double precision :: flux(NCVRS)

	double precision :: nx,ny,nz
	double precision :: aL,aR,vL,vR,ML,MR
	double precision :: Mvl,MLplus,MRminus
	double precision :: PLplus,PRminus
	
	nx = n(DIRX)
	ny = n(DIRY)
	nz = n(DIRZ)

	call cons_to_primitive(ul,g,pl)
	call cons_to_primitive(ur,g,pr)
	
	aL = sqrt(g*pl(PRES_INDX)/pl(DENS_INDX))
	aR = sqrt(g*pr(PRES_INDX)/pr(DENS_INDX))

	vL = pl(VELX_INDX)*nx + pl(VELY_INDX)*ny + pl(VELZ_INDX)*nz
	vR = pr(VELX_INDX)*nx + pr(VELY_INDX)*ny + pr(VELZ_INDX)*nz
	
	ML = vL/aL
	MR = vR/aR

	call vleersplit_Mach(ML, 1.d0 ,MLplus);
	call vleersplit_Mach(MR,-1.d0,MRminus);

	Mvl = MLplus + MRminus

	if(Mvl .gt. 0) then
		flux(RHO_INDX)  = ul(RHO_INDX)*aL	
		flux(RHOU_INDX) = ul(RHOU_INDX)*aL	
		flux(RHOV_INDX) = ul(RHOV_INDX)*aL	
		flux(RHOW_INDX) = ul(RHOW_INDX)*aL	
		flux(RHOE_INDX) = (ul(RHOE_INDX)+pl(PRES_INDX))*aL
	else
		flux(RHO_INDX)  = ur(RHO_INDX)*aR	
		flux(RHOU_INDX) = ur(RHOU_INDX)*aR	
		flux(RHOV_INDX) = ur(RHOV_INDX)*aR	
		flux(RHOW_INDX) = ur(RHOW_INDX)*aR	
		flux(RHOE_INDX) = (ur(RHOE_INDX)+pr(PRES_INDX))*aR
	endif

	flux = flux*Mvl;

	!pressure term
	call vleersplit_pres(pl(PRES_INDX),ML, 1.d0,PLplus)
	call vleersplit_pres(pr(PRES_INDX),MR,-1.d0,PRminus)

	flux(RHOU_INDX) = flux(RHOU_INDX) +  (PLplus+PRminus)*nx
	flux(RHOV_INDX) = flux(RHOV_INDX) +  (PLplus+PRminus)*ny
	flux(RHOW_INDX) = flux(RHOW_INDX) +  (PLplus+PRminus)*nz
			

end subroutine ausm
!=====================================================
subroutine vleersplit_Mach(M,plus_or_minus,Mvl)
	
	implicit none	
	double precision :: M,plus_or_minus,Mvl;

	!Mvl(+-) = (+-)0.25*(M (+-) 1)**2   !subsonic
	!Mvl(+-) = 0.5 * (M (+-) |M|)       !supersonic
	
	if(plus_or_minus .lt. 0) then

		if(abs(M) .le. 1) then
			Mvl = -0.25*(M-1)**2
		else
			Mvl = 0.5*(M-abs(M))
		endif
	else
		if(abs(M) .le. 1) then
			Mvl = 0.25*(M+1)**2
		else
			Mvl = 0.5*(M+abs(M))
		endif
	endif

end subroutine vleersplit_Mach
!=====================================================
subroutine vleersplit_pres(P,M,plus_or_minus,Pvl)
	
	implicit none	
	double precision :: P,M,plus_or_minus,Pvl;

	!Pvl(+-) = (+-)0.25P (M (+-) 1)**2 (2 (-+) M)  !subsonic
	!Mvl(+-) = 0.5* P * (M (+-) |M|)/M             !supersonic
	
	if(plus_or_minus .lt. 0) then

		if(abs(M) .le. 1) then
			Pvl = 0.25 * P * (M-1)**2 * (2+M)
		else
			Pvl =  0.5 * P * (M-abs(M))/M
		endif
	else
		if(abs(M) .le. 1) then
			Pvl = 0.25 * P * (M+1)**2 * (2-M)
		else
			Pvl = 0.5 * P * (M+abs(M))/M
		endif
	endif

end subroutine vleersplit_pres
!=====================================================
subroutine centrdiff(ul,ur,g,n,flux)
	
	implicit none	
	double precision :: g
	double precision :: n(THREEDIM)
	double precision :: ul(NCVRS),ur(NCVRS)
	double precision :: fl(NCVRS),fr(NCVRS)
	double precision :: flux(NCVRS)

	call getflux(ul,g,n,fl) 
	call getflux(ur,g,n,fr) 
	
	flux(:) = 0.5*(fl(:)+fr(:))

end subroutine centrdiff
!=====================================================
subroutine getflux(u,g,n,flux)

	implicit none	
	double precision :: g
	double precision :: u(NCVRS),n(THREEDIM)
	double precision :: flux(NCVRS)
	double precision :: p,vx,vy,vz,rho,rho_e
	double precision :: nx,ny,nz
	double precision :: pvars(NCVRS)

	call cons_to_primitive(u,g,pvars)

	rho = pvars(DENS_INDX)
	vx  = pvars(VELX_INDX)
	vy  = pvars(VELY_INDX)
	vz  = pvars(VELZ_INDX)
	p   = pvars(PRES_INDX)

	rho_e = u(RHOE_INDX)

	nx = n(DIRX)
	ny = n(DIRY)
	nz = n(DIRZ)

	flux(RHO_INDX)  = rho*(vx*nx + vy*ny + vz*nz)

	flux(RHOU_INDX) = (rho*vx**2+p)*nx + rho*vx*vy*ny     + rho*vx*vz*nz
	flux(RHOV_INDX) = rho*vy*vx*nx     + (rho*vy**2+p)*ny + rho*vy*vz*nz
	flux(RHOW_INDX) = rho*vz*vx*nx      + rho*vz*vy*ny    + (rho*vz**2+p)*nz

	flux(RHOE_INDX) = (rho_e+p)*(vx*nx + vy*ny + vz*nz)

end subroutine getflux
!=====================================================
subroutine cons_to_primitive(u,g,primvars)
	
	implicit none	
	double precision :: u(NCVRS),primvars(NCVRS)
	double precision :: g

	double precision :: p,vx,vy,vz,rho

	rho = u(RHO_INDX)
	primvars(DENS_INDX) = u(RHO_INDX)
	
	vx = u(RHOU_INDX)/u(RHO_INDX)
	vy = u(RHOV_INDX)/u(RHO_INDX)
	vz = u(RHOW_INDX)/u(RHO_INDX)

	primvars(VELX_INDX)=vx
	primvars(VELY_INDX)=vy
	primvars(VELZ_INDX)=vz

	primvars(PRES_INDX) = (g-1)*( u(RHOE_INDX) - 0.5*rho*(vx**2+vy**2+vz**2) )

end subroutine cons_to_primitive
!=====================================================
subroutine prim_to_conservative(pvars,g,u)

	implicit none	
	double precision :: u(NCVRS),pvars(NCVRS)
	double precision :: g

	double precision :: p,vx,vy,vz,rho

	u(RHO_INDX) = pvars(DENS_INDX)

	u(RHOU_INDX) = pvars(DENS_INDX)*pvars(VELX_INDX)	
	u(RHOV_INDX) = pvars(DENS_INDX)*pvars(VELY_INDX)	
	u(RHOW_INDX) = pvars(DENS_INDX)*pvars(VELZ_INDX)	

	u(RHOE_INDX) = pvars(PRES_INDX)/(g-1)
	u(RHOE_INDX) = u(RHOE_INDX) + 0.5*pvars(DENS_INDX)*( pvars(VELX_INDX)**2+pvars(VELY_INDX)**2+pvars(VELZ_INDX)**2 )

end subroutine prim_to_conservative
!=====================================================

end module inviscidterms_module
