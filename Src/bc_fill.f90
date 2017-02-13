module bc_fill_module

! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
! for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  subroutine fillboundary(phi,lo,hi,domlo,domhi,delta,xlo,time,bc,comp) &
       bind(C, name="fillboundary")

    use bl_fort_module, only : bl_spacedim, c_real
    use dirichlet_bcs

    implicit none
    include 'bc_types.fi'

    integer      :: comp
    integer      :: lo(3),hi(3)
    integer      :: bc(bl_spacedim,2)
    integer      :: domlo(3), domhi(3)
    real(c_real) :: delta(3), xlo(3), time
    real(c_real) :: phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    integer    i, j, k
    
    call filcc(phi,lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( (bc(1,1).eq.EXT_DIR) .and. (lo(1) .lt. domlo(1)) ) then
       do i = lo(1), domlo(1)-1
          do j = lo(2), hi(2)
             do k = lo(3), hi(3)
                   phi(i,j,k) = dirc_bcs(1,currentcomponent)
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2).eq.EXT_DIR) .and. (hi(1).gt.domhi(1)) ) then
       do i = domhi(1)+1, hi(1)
          do j = lo(2), hi(2)
             do k = lo(3), hi(3)
                   phi(i,j,k) = dirc_bcs(2,currentcomponent)
             end do
          end do
       end do
    end if

    !     YLO
    if ( (bc(2,1).eq.EXT_DIR) .and. (lo(2).lt.domlo(2)) ) then
       do i = lo(1), hi(1)
          do j = lo(2), domlo(2)-1
             do k = lo(3), hi(3)
                   phi(i,j,k) = dirc_bcs(3,currentcomponent)
             end do
          end do
       end do
    end if

    !     YHI
    if ( (bc(2,2).eq.EXT_DIR) .and. (hi(2).gt.domhi(2)) ) then
       do i = lo(1), hi(1)
          do j = domhi(2)+1, hi(2)
             do k = lo(3), hi(3)
                   phi(i,j,k) = dirc_bcs(4,currentcomponent)
             end do
          end do
       end do
    end if

    !     ZLO
    if ( (bc(3,1).eq.EXT_DIR) .and. (lo(3).lt.domlo(3)) ) then
       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             do k = lo(3), domlo(3)-1
                   phi(i,j,k) = dirc_bcs(5,currentcomponent)
             end do
          end do
       end do
    end if

    !     ZHI
    if ( (bc(3,2).eq.EXT_DIR) .and. (hi(3).gt.domhi(3)) ) then
       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             do k = domhi(3)+1, hi(3)
                   phi(i,j,k) = dirc_bcs(6,currentcomponent)
             end do
          end do
       end do
    end if

  end subroutine fillboundary
  
end module bc_fill_module
