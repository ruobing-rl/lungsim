module surfactant_c

  !use surfactant


  implicit none

  private 


contains
!
!##############################################################################

  subroutine calculate_surface_tension_c(volume_mean, volume_change, frequency, volumes, radii, area, dA, &
    surf_concentration, gamma_star, gamma_max, bulk_c, k_a, k_d, m2, sigma, sigma_hat) bind(C, name="calculate_surface_tension_c")

    use surfactant, only: calculate_surface_tension
    implicit none


    real(dp) :: volume_mean, volume_change, frequency, gamma_star, gamma_max, bulk_c, k_a, k_d, m2, sigma_hat
    real(dp), allocatable :: volumes(:), radii(:), area(:), dA(:), surf_concentration(:), sigma(:), Pc_alv(:)
    integer :: num_steps, i
    real(dp) :: dt, endt

    call calculate_surface_tension(volume_mean, volume_change, frequency, volumes, radii, area, dA, surf_concentration, &
            gamma_star, gamma_max, bulk_c, k_a, k_d, m2, sigma, sigma_hat)


  end subroutine calculate_surface_tension_c

!##############################################################################

end module surfactant_c


