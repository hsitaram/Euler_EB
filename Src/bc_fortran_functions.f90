!==================================================
subroutine update_currentcomponent(val) bind(C,name="update_currentcomponent")

	use dirichlet_bcs
	implicit none

	integer :: val
	currentcomponent=val

end subroutine update_currentcomponent
!==================================================
subroutine set_dirichlet_bcs(d,vx,vy,vz,p) bind(C,name="set_dirichlet_bcs")

	use dirichlet_bcs
	implicit none

	double precision :: d(6)
	double precision :: vx(6),vy(6),vz(6)
	double precision :: p(6)

	dirc_bcs(:,1)=d
	dirc_bcs(:,2)=vx
	dirc_bcs(:,3)=vy
	dirc_bcs(:,4)=vz
	dirc_bcs(:,5)=p

end subroutine set_dirichlet_bcs
!==================================================
