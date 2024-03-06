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
  

  !use ventilation
  use precision
  use diagnostics


  implicit none


  !Module parameters

  real(dp), parameter :: pi = 3.14159265358979_dp
  real(dp), parameter :: frequency = 0.25_dp ! frequency of breathing per minute


  !Module types

  !Module variables

  !Interfaces
  private
  public update_surface_tension
  !private
  !public 
  

    
    
contains
!
!##############################################################################
!
!*subname:* calculate the change of surface area for each alveolar unit

  subroutine update_surface_area(volumes, radii, area, dA)

    !real(dp), intent(in) :: volume_mean, volume_change

    real(dp), intent(out) :: volumes(:), radii(:), area(:), dA(:)

    integer :: num_steps, i

    real(dp) :: dt, t,endt

    real(dp), parameter :: volume_mean = 4.2e-6_dp
    real(dp), parameter :: volume_change = 0.75e-6_dp
        
    character(len=60) :: sub_name


    endt=60.0_dp
    dt = 0.005_dp  ! Time step
    num_steps = int(endt / dt) + 1
        
! --------------------------------------------------------------------------

  sub_name = 'update_surface_area'
    call enter_exit(sub_name,1)
  
     
    do i = 1, num_steps !size(volumes)
      t = (i-1) * dt
      volumes(i) = volume_mean + volume_change * sin(2* pi * frequency * t + (pi / 2))
      radii(i) = (3 * volumes(i) / (4 * pi)) ** (1 / 3)
      area(i)=4 * pi * (radii(i)**2);
    end do
        
        
    dA(1)=0.0_dp !First element is 0
        
    do i = 1, num_steps
      dA(i+1)=(area(i+1)-area(i))/dt; 
    end do
        

        
  call enter_exit(sub_name,2)

  end subroutine update_surface_area
!##############################################################################
!
!*subname:* calculate surfactant concentration in each alveolar unit

    subroutine update_surfactant_concentration(area, dA, surf_concentration)

    !real(dp), intent(in) :: volume_mean, volume_change, gamma_star, gamma_max, bulk_c, k_a, & k_d, sigma_hat, m2

    real(dp), intent(in) :: area(:), dA(:)

    real(dp), intent(out) :: surf_concentration(:)

    integer :: num_steps, i

    real(dp) :: dt, t,endt

    real(dp), parameter :: gamma_star = 0.3e-6_dp !
    real(dp), parameter :: gamma_max = 0.345e-6_dp !
    real(dp), parameter :: bulk_c = 1e-3_dp ! bulk concentration  g/ml
    real(dp), parameter :: k_a = 10**4 !6*10**5 ! adsorption coefficient  ml/(g*sec)
    real(dp), parameter :: k_d = k_a/(1.2_dp*(10**5)) !desorption coefficient sec^(-1)


    character(len=60) :: sub_name


    endt=60.0_dp
    dt = 0.005_dp  ! Time step
    num_steps = int(endt / dt) + 1


! --------------------------------------------------------------------------

  sub_name = 'surf_concentration'
    call enter_exit(sub_name,1)



    !k_d = k_a/(1.2*(10**5))
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



  call enter_exit(sub_name,2)

  end subroutine update_surfactant_concentration
!##############################################################################
!
!*subname:* calculate surface tension for each alveolar unit

    subroutine update_surface_tension(surf_concentration, surface_tension)

    !real(dp), intent(in) :: volume_mean, volume_change, gamma_star, gamma_max, bulk_c, k_a, & k_d, sigma_hat, m2

    real(dp), intent(out) :: surface_tension(:)

    integer :: num_steps, i

    real(dp) :: dt, t,endt
    real(dp) :: surf_concentration(:)

    real(dp), parameter :: gamma_star = 0.3e-6_dp !
    real(dp), parameter :: gamma_max = 0.345e-6_dp !
    real(dp), parameter :: surface_tension_hat = 22.02_dp !dyn/cm
    real(dp), parameter :: m2 = 140.0_dp !slope


    character(len=60) :: sub_name


    endt=60.0_dp
    dt = 0.005_dp  ! Time step
    num_steps = int(endt / dt) + 1

! --------------------------------------------------------------------------

  sub_name = 'surface_tension'
    call enter_exit(sub_name,1)


    !call update_surfactant_concentration

    do i = 1, num_steps
      if (surf_concentration(i) .le. gamma_star) then
        surface_tension(i)=(surface_tension_hat-70)*(surf_concentration(i)/gamma_star)+70

      elseif ((surf_concentration(i) .ge. gamma_max)) then!.and. (dA(i) .lt. 0.0)) then
        surface_tension(i)=-m2*(gamma_max/gamma_star)+(m2+surface_tension_hat)

      else
        surface_tension(i)=-m2*(surf_concentration(i)/gamma_star)+(m2+surface_tension_hat)
      end if
    end do

  call enter_exit(sub_name,2)

  end subroutine update_surface_tension
!##############################################################################


  subroutine write_results(t,volumes,radii,area,dA, &
       surf_concentration,surface_tension)

    real(dp), allocatable :: volumes(:), radii(:), area(:), dA(:), surf_concentration(:), surface_tension(:)
    integer :: num_steps, i
    real :: dt, endt,t


    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'write_results'
    call enter_exit(sub_name,1)

    !do i = 1, num_steps
    !   print '("At t = ", F8.4, "  Volume = ", F12.8, "  Radius = ", F12.8, "  Area = ", F12.8, "  dA = ", F12.8,  &
    !    " surf_concentration = ",F12.8, " surface_tension= ", F12.8)', &
    !    (i-1) * 0.005_dp , volumes(i), radii(i), area(i), dA(i), surf_concentration(i), surface_tension(i)

    !end do

    num_steps = int(60 / 0.00083 ) + 1
    !num_steps = int(endt / dt) + 1


    allocate(volumes(num_steps))
    allocate(radii(num_steps))
    allocate(area(num_steps))
    allocate(dA(num_steps))
    allocate(surf_concentration(num_steps))
    allocate(surface_tension(num_steps))

    write(*,'(2X,''Time'',3X,''Volume'',4X,''Radius'',5X,''Area'',5X,&
            &''dA'',4X,''Surf_concentration'',5X,''surface_tension'')')

    write(*,'(F7.3,2(F8.1),8(F8.2))') &
         (i-1) * 0.005_dp , &  !time, flow, tidal
         volumes , & !res (cmH2O/L.s)
         radii , & !total model compliance
         area, & !Ppl (cmH2O)
         dA, & !mean Ptp (cmH2O)
         surf_concentration, & !total model volume (L)
         surface_tension  !Pmuscle (cmH2O)


    call enter_exit(sub_name,2)

  end subroutine write_results
!##############################################################################


end module surfactant