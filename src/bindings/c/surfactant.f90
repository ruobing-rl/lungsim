module surfactant_c

  !use surfactant


  implicit none

  private 


contains
!
!##############################################################################

  subroutine update_surface_tension_c(surf_concentration, surface_tension) &
          bind(C, name="update_surface_tension_c")

    use surfactant, only: update_surface_tension
    use arrays,only: dp
    implicit none

    real(dp), intent(out) :: surface_tension(:)

    real(dp) :: surf_concentration(:)


    call update_surface_tension(surf_concentration, surface_tension)


  end subroutine update_surface_tension_c

!##############################################################################

end module surfactant_c


