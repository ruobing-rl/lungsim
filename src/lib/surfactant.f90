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

!  use ventilation
!  use geometry
  use precision
  use diagnostics
  use arrays
  use indices
  use other_consts




  implicit none


  !Module parameters

  !real(dp), parameter :: frequency = 0.25_dp ! frequency of breathing per minute

  !Module types

  !Module variables

  !Interfaces

  private
  public evaluate_surf
  public update_surfactant_concentration
  public update_surface_tension
  public update_Pc_compliance
!  public alv_collapse_pressure

contains
!
!##############################################################################
!
subroutine evaluate_surf !(num_steps, dt, t, volumes, radii, area, dA, surf_concentration, surface_tension)


    !integer :: i,ntime
    integer :: num_steps, num_time,num_units,nu_vol,nalv

    real(dp) :: dt !, t, endt,time

!    real(dp), allocatable :: volumes(:)
!    real(dp), allocatable :: radii(:)
!    real(dp), allocatable :: area(:)
!    real(dp), allocatable :: dA(:)
    real(dp), dimension(:,:), allocatable ::alv_area_current
    real(dp), dimension(:,:), allocatable ::alv_dA
    real(dp), dimension(:,:), allocatable ::alv_radii_current
    real(dp), dimension(:,:), allocatable :: surf_concentration
    real(dp), dimension(:,:), allocatable :: surface_tension
    real(dp), dimension(:,:), allocatable :: Pc
!    real(dp), dimension(:,:), allocatable :: alv_collapse_pressure_pre
!    real(dp), dimension(:,:), allocatable :: alv_collapse_pressure_current


    logical :: CONTINUE,converged
    character(len=60) :: sub_name

! --------------------------------------------------------------------------
    sub_name = 'evaluate_surf'
    call enter_exit(sub_name,1)

!      endt=12.0_dp
!      dt = 0.05_dp
!      num_steps = int(endt / dt ) + 1
!      num_steps = int(time / dt )
!      num_time=240.0_dp

!      allocate(volumes(num_steps))
!      allocate(radii(num_steps))
!      allocate(area(num_steps))
!      allocate(dA(num_steps))
!      allocate(surf_concentration(nu_vol,num_units))
!      allocate(surface_tension(nu_vol,num_units))


      !call evaluate_vent

      !print*, 'size', alv_area(1,1)
      !print*, 'row', nu_vol
      !print*, 'column', num_units


!      call update_alveolar_info(dt,time
     ! call update_unit_volume(dt)
      !call update_surface_area(volumes, radii, area, dA)

      call update_surfactant_concentration(dt, alv_area_current, alv_dA, surf_concentration)
      call update_surface_tension( surf_concentration, surface_tension, alv_radii_current, Pc)
!      call alv_collapse_pressure(alv_radii_current, surface_tension, alv_collapse_pressure)

      !call write_results(volumes,radii,area,dA, surf_concentration,surface_tension)

      !print *, "Volume, Radius, Surface Area, surfactant concentration and surface tention values at selected timesteps:"
!      do i = 1, num_steps
!        print '("At t = ", F8.4, "  Volume = ", F12.8, "  Radius = ", F12.8, "  Area = ", F12.8, "  dA = ", &
!        F12.8,"  surf_con = ", F12.8,"  surf_ten = ", F12.8)', (i-1) * dt , volumes(i), radii(i), &
!        area(i), dA(i), surf_concentration(i), surface_tension(i)
!
!      end do



!      deallocate(alv_area_current)
!      deallocate(alv_dA)
!      deallocate(alv_radii_current)
!      deallocate(alv_collapse_pressure)
!      deallocate(surf_concentration)
!      deallocate(surface_tension)
!       deallocate (alv_collapse_pressure_pre)
!       deallocate (alv_collapse_pressure_current)

    call enter_exit(sub_name,2)


  end subroutine evaluate_surf

!##############################################################################
!
!*subname:* calculate surfactant concentration in each alveolar unit

  subroutine update_surfactant_concentration(dt, alv_area_current, alv_dA, surf_concentration)

    real(dp), dimension(:,:), intent(in) :: alv_area_current, alv_dA

    real(dp), dimension(:,:), intent(inout) :: surf_concentration

    integer :: nalv

    real(dp), intent(in) :: dt

    real(dp), parameter :: gamma_star = 0.3e-6_dp !g/cm^2
    real(dp), parameter :: gamma_max = 0.345e-6_dp !g/cm^2
    real(dp), parameter :: bulk_c = 10e-3_dp ! bulk concentration  g/ml
    real(dp), parameter :: k_a = 10.0_dp**4 ! adsorption coefficient  ml/(g*sec)
    real(dp), parameter :: k_d = k_a/(1.2_dp*(10.0_dp**5)) !desorption coefficient sec^(-1)
!    real(dp), parameter :: surf_concentration= gamma_star/2.0_dp

    character(len=60) :: sub_name



 !--------------------------------------------------------------------------

  sub_name = 'surf_concentration'
    call enter_exit(sub_name,1)


    do nalv = 1,num_units



      if (surf_concentration(nu_vol,nalv) .le. gamma_star) then
        surf_concentration(nu_vol,nalv)=surf_concentration(nu_vol,nalv)+dt*(k_a*bulk_c &
                *(gamma_star-surf_concentration(nu_vol,nalv))-k_d*surf_concentration(nu_vol,nalv) &
                -(surf_concentration(nu_vol,nalv)/alv_area_current(nu_vol,nalv))*alv_dA(nu_vol,nalv))

      elseif ((surf_concentration(nu_vol,nalv) .ge. gamma_max) .and. (alv_dA(nu_vol,nalv) .lt. 0.0)) then
        surf_concentration(nu_vol,nalv)=gamma_max

      else
        surf_concentration(nu_vol,nalv)=surf_concentration(nu_vol,nalv)+ &
                dt*((-surf_concentration(nu_vol,nalv)/alv_area_current(nu_vol,nalv))*alv_dA(nu_vol,nalv))

      end if
    end do
!    print*,'dt2',dt
!    print*, 'alveolar area2', alv_area(nu_vol,1)
!    print*, 'dA2', alv_dA(nu_vol,1)
!    print*,'surf_concentration2',surf_concentration(nu_vol,1)
    !deallocate(surf_concentration)
    call enter_exit(sub_name,2)

  end subroutine update_surfactant_concentration
!##############################################################################
!
!*subname:* calculate surface tension an collapse pressure for each alveolar unit

  subroutine update_surface_tension(surf_concentration, surface_tension, alv_radii_current, Pc)

    real(dp), dimension(:,:), intent(in) :: surf_concentration
    real(dp), dimension(:,:), intent(in) :: alv_radii_current
!    real(dp) :: alv_area_current
    real(dp), dimension(:,:), intent(inout) :: surface_tension
    real(dp), dimension(:,:), intent(out) :: Pc
!    real(dp), dimension(:,:), allocatable :: Pc_pre
!    real(dp), dimension(:,:), allocatable :: Pc_current

    integer :: nalv

    real(dp) :: dt


    real(dp), parameter :: gamma_star = 0.3e-6_dp !
    real(dp), parameter :: surface_tension_hat = 22.02_dp !dyn/cm
    real(dp), parameter :: surface_tension_min= 1.0_dp !dyn/cm
    real(dp), parameter :: m2 = 140.0_dp !slope
    real(dp), parameter :: gamma_max = 3.2e-7_dp!gamma_star*(1+(surface_tension_hat-surface_tension_min)/m2)!0.345e-6_dp !

    character(len=60) :: sub_name

 !--------------------------------------------------------------------------

  sub_name = 'surface_tension'
    call enter_exit(sub_name,1)



    do nalv = 1,num_units

      if (surf_concentration(nu_vol,nalv) .le. gamma_star) then
        surface_tension(nu_vol,nalv)=(surface_tension_hat-70)*(surf_concentration(nu_vol,nalv)/gamma_star)+70

      elseif ((surf_concentration(nu_vol,nalv) .ge. gamma_max)) then!.and. (dA(i) .lt. 0.0)) then
        surface_tension(nu_vol,nalv)=-m2*(gamma_max/gamma_star)+(m2+surface_tension_hat)

      else
        surface_tension(nu_vol,nalv)=-m2*(surf_concentration(nu_vol,nalv)/gamma_star)+(m2+surface_tension_hat)
      end if

!      Pc_pre(nu_vol,nalv) = alv_collapse_pressure(nu_vol,nalv)
      Pc(nu_vol,nalv)= ((2.0 *surface_tension(nu_vol,nalv)) /alv_radii_current(nu_vol,nalv))/10.0
!      Pc(nu_vol,nalv)= -((alv_area_current(nu_vol,nalv)*(alv_radii_current(nu_vol,nalv)**2)))/(2.0*surface_tension(nu_vol,nalv))
      !dyn/cm^2 to Pa (/10.0)
!      alv_collapse_pressure(nu_vol,nalv) = alv_collapse_pressure(nu_vol,nalv)
!      Pc_current(nu_vol,nalv) =alv_collapse_pressure(nu_vol,nalv)

    end do


!    print*,'surface tension2',Pc_pre(nu_vol,1)
!    print*,'surface tension2',Pc_current(nu_vol,1)

!    print*,'surface tension2',surface_tension(nu_vol,1)

    call enter_exit(sub_name,2)

  end subroutine update_surface_tension

!##############################################################################
!!
!!*subname:* calculate compliance calculated by surface tension for each acinar unit
!
  subroutine update_Pc_compliance(unit_field_current, Pc, Pc_com, smoothed_Pc_com)
!

    real(dp), dimension(:,:), intent(in) :: unit_field_current
    real(dp), dimension(:,:), intent(in) :: Pc
    real(dp), dimension(:,:), intent(out) :: Pc_com

    integer, parameter :: window_size = 5!101 ! Choose window size (example: 10% of the data length)

    real(dp), dimension(window_size) :: kernel
    real(dp), dimension(:,:), intent(out) :: smoothed_Pc_com
!    real(8) :: sigma, sum, gaussian_kernel(100)  ! Adjust array size if needed
!    real(8) :: p(n), smoothed_p(n)
    real(dp):: sum,sum_weights, value, weight !,mean_smoothed,mean_original

    integer :: nalv,j

    real(dp) :: dt


    real(dp), parameter :: sigma = 10.5!10.5_dp!window_size/4.0_dp


!    integer :: nalv, nu_vol, i, j
!    real(8), dimension(num_units) :: Pc_com
!    real(8), dimension(kernel_size) :: kernel
!    real(8), dimension(num_units) :: smoothed_Pc_com
!    real(8) :: sigma, sum_weights, value, weight

    ! Example values for testing


    ! Initialize example values (replace with actual data)


    character(len=60) :: sub_name
!
! !--------------------------------------------------------------------------
!
  sub_name = 'Pc_compliance'
    call enter_exit(sub_name,1)


    !Pc_compliance of acinus unit
    do nalv = 1,num_units

        Pc_com(nu_vol,nalv)= (3.0 *unit_field_current(nu_vol,nalv)) /Pc(nu_vol,nalv)

    end do

    ! Calculate mean of the original data
!    mean_original = sum(Pc_com) / num_units


!    moving average
        do nalv = 1, num_units
            sum_weights = 0.0_dp
            value = 0.0_dp
            do j =  max(1, nalv - window_size), min(num_units, nalv + window_size)
!                if (nalv + j > 0 .and. nalv + j <= num_units) then
                    value = value + Pc_com(nu_vol, nalv + j)
                    sum_weights = sum_weights + 1.0_dp
!                end if
            end do
            smoothed_Pc_com(nu_vol, nalv) = value / sum_weights
        end do

        ! Compute Gaussian kernel values
!    kernel = gaussian_kernel(window_size, sigma)

    ! Smooth the Pc_com values for this volume
!    do nalv = 1, num_units
!        sum_weights = 0.0_dp
!        value = 0.0_dp
!        do j = -window_size/2, window_size/2
!            if (nalv + j > 0 .and. nalv + j <= num_units) then
!                weight = kernel(j + window_size/2 + 1)
!                value= value + Pc_com(nu_vol, nalv + j) * weight
!                sum_weights = sum_weights + weight
!            end if
!        end do
!        smoothed_Pc_com(nu_vol,nalv) = value / sum_weights
!        if (sum_weights > 0.0_dp) then
!            smoothed_Pc_com(nu_vol, nalv) = value / sum_weights
!        else
!            smoothed_Pc_com(nu_vol, nalv) = Pc_com(nu_vol, nalv)
!        end if

!    end do

    ! Calculate mean of the smoothed data
!    mean_smoothed = sum(smoothed_Pc_com) / num_units

    ! Adjust the smoothed data by subtracting the mean shift
!    smoothed_Pc_com = smoothed_Pc_com - (mean_smoothed - mean_original)

     ! Initialize buffer and counters
!     data_buffer = 0.0
!     buffer_index = 0



    call enter_exit(sub_name,2)
!
  end subroutine update_Pc_compliance


!##############################################################################

function gaussian_kernel(size, sigma) result(kernel)

    integer, intent(in) :: size
    real(dp), intent(in) :: sigma
    real(dp), dimension(size):: kernel
    integer :: i
    real(dp) :: sum, exponent
!
! !--------------------------------------------------------------------------
!
!    sum = 0.0_dp
!
!    do i = -size/2, size/2
!        exponent = -0.5_dp * (i / sigma )**2
!        kernel(i + size/2 + 1) = exp(exponent) / (sigma * sqrt(2.0_dp * PI))
!        sum = sum + kernel(i + size/2 + 1)
!    end do
!
!    ! Normalize the kernel
!    do i = 1, size
!        kernel(i) = kernel(i) / sum
!
!    end do
!
!end function gaussian_kernel

        sum = 0.0_dp
        do i = 0, size-1
            exponent = -0.5_dp * ((i - size/2) / sigma)**2
            kernel(i + 1) = exp(exponent) / (sigma * sqrt(2.0_dp * PI))
            sum = sum + kernel(i + 1)
        end do

        ! Normalize the kernel
        do i = 1, size
            kernel(i) = kernel(i) / sum
        end do
end function gaussian_kernel
!##############################################################################
!
!function moving_average(y, n, size) result(smoothed)
!        integer, intent(in) :: n, size
!        real, intent(in) :: y(n)
!        real :: smoothed(n)
!        integer :: i, j
!        real :: sum
!
!        ! Loop through each data point
!        do i = 1, n
!            sum = 0.0
!            ! Calculate the moving average within the window
!            do j = max(1, i-size/2), min(n, i+size/2)
!                sum = sum + y(j)
!            end do
!            smoothed(i) = sum / (min(n, i+size/2) - max(1, i-size/2) + 1)
!        end do
!end function moving_average
!    moving average
!        do nalv = 1, num_units
!            sum_weights = 0.0_dp
!            value = 0.0_dp
!            do j =  max(1, nalv - window_size), min(num_units, nalv + window_size)
!!                if (nalv + j > 0 .and. nalv + j <= num_units) then
!                    value = value + Pc_com(nu_vol, nalv + j)
!                    sum_weights = sum_weights + 1.0_dp
!!                end if
!            end do
!            smoothed_Pc_com(nu_vol, nalv) = value / sum_weights
!        end do

!
!##############################################################################
!##############################################################################

end module surfactant