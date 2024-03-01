module surfactant
!*Brief Description:* This module handles all code specific to
! simulating surfactant in the acinus.
!
!*LICENSE:*
!
!
!
!*Full Description:*
!This module handles all code specific to simulating surfactant in the acinus.
  

  use ventilation

  implicit none

  public calculate_surface_tension


  !Module parameters

  real(dp), parameter :: pi = 3.14159265358979


  !Module types


  !Module variables

  !Interfaces
  !private list1
  !public list2
  !public list3
  !public 
  

    
    
contains
!
!##############################################################################
!
!*subname:* calculate surface tension for unit

  subroutine calculate_surface_tension(volume_mean, volume_change, frequency, volumes, radii, area, dA, &
    surf_concentration, gamma_star, gamma_max, bulk_c, k_a, k_d, m2, sigma, sigma_hat)

    real(dp), intent(in) :: volume_mean, volume_change, frequency, gamma_star, gamma_max, bulk_c, k_a, & 
      k_d, sigma_hat, m2

    real(dp), intent(out) :: volumes(:), radii(:), area(:), dA(:), surf_concentration(:), sigma(:)

    integer :: num_steps, i

    real(dp) :: dt, t,endt
        
    character(len=60) :: sub_name


    endt=60.0
    dt = 0.005  ! Time step
    num_steps = int(endt / dt) + 1
        
! --------------------------------------------------------------------------

  sub_name = 'calculate_surface_tension'
    call enter_exit(sub_name,1)
  
     
    do i = 1, num_steps !size(volumes)
      t = (i-1) * dt
      volumes(i) = volume_mean + volume_change * sin(2.0 * pi * frequency * t + (pi / 2.0))
      radii(i) = (3.0 * volumes(i) / (4.0 * pi)) ** (1.0 / 3.0)
      area(i)=4.0*pi*(radii(i)**2.0);
    end do
        
        
    dA(1)=0.0 !First element is 0
        
    do i = 1, num_steps
      dA(i+1)=(area(i+1)-area(i))/dt; 
    end do
        
        
        
    surf_concentration(1) = gamma_star/2

    do i = 1, num_steps
      if (surf_concentration(i) .le. gamma_star) then     
        surf_concentration(i+1)=surf_concentration(i)+dt*(k_a*bulk_c & 
          *(gamma_star-surf_concentration(i))-k_d*surf_concentration(i)-(surf_concentration(i)/area(i))*dA(i))

      elseif ((surf_concentration(i) .ge. gamma_max) .and. (dA(i) .lt. 0.0)) then
        surf_concentration(i+1)=gamma_max
                
      else
        surf_concentration(i+1)=surf_concentration(i)+ dt*((-surf_concentration(i)/area(i))*dA(i))
                
      end if
    end do
        
        
    do i = 1, num_steps
      if (surf_concentration(i) .le. gamma_star) then 
        sigma(i)=(sigma_hat-70)*(surf_concentration(i)/gamma_star)+70

      elseif ((surf_concentration(i) .ge. gamma_max)) then!.and. (dA(i) .lt. 0.0)) then
        sigma(i)=-m2*(gamma_max/gamma_star)+(m2+sigma_hat)

      else
        sigma(i)=-m2*(surf_concentration(i)/gamma_star)+(m2+sigma_hat)
      end if
    end do
        
  call enter_exit(sub_name,2)

  end subroutine calculate_surface_tension
!##############################################################################

end module surfactant