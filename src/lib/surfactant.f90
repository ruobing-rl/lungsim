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




  implicit none


  !Module parameters


  !real(dp), parameter :: pi = 3.14159265358979_dp
  !real(dp), parameter :: frequency = 0.25_dp ! frequency of breathing per minute


  !Module types

  !Module variables

  !real(dp),  allocatable :: volumes(:)
  !real(dp),  allocatable :: radii(:)
  !real(dp),  allocatable :: area(:)
  !real(dp),  allocatable :: dA(:)
  !real(dp),  allocatable :: surf_concentration(:)
  !real(dp),  allocatable :: surface_tension(:)







  !Interfaces

  !private
  !, volumes, radii, area, dA, surf_concentration, surface_tension
  
  private
  !public update_surface_tension
  public evaluate_surf
  public update_surfactant_concentration
  public update_surface_tension

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
    real(dp), dimension(:,:), allocatable :: alv_collapse_pressure
    !real(dp), dimension(:,:), allocatable :: alv_radii
    !real(dp), dimension(:,:), allocatable :: alv_area


    !real(dp), dimension(:,:), allocatable :: alv_area_over_breath

    logical :: CONTINUE,converged
    character(len=60) :: sub_name

! --------------------------------------------------------------------------
    sub_name = 'evaluate_surf'
    call enter_exit(sub_name,1)

!      endt=12.0_dp
!      dt = 0.05_dp
!      num_steps = int(endt / dt ) + 1
!      !num_steps = int(time / dt )
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


!      call update_alveolar_info(dt,time)
     ! call update_unit_volume(dt)
      !call update_surface_area(volumes, radii, area, dA)

      call update_surfactant_concentration(dt, alv_area_current, alv_dA, surf_concentration)
      call update_surface_tension( surf_concentration, surface_tension, alv_radii_current, alv_collapse_pressure)
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

    call enter_exit(sub_name,2)


  end subroutine evaluate_surf
!##############################################################################
!
!*subname:* calculate the change of surface area for each alveolar unit

!  subroutine update_surface_area(volumes, radii, area, dA)
!
!    !real(dp), intent(in) :: volume_mean, volume_change
!
!    real(dp), intent(out) :: volumes(:), radii(:), area(:), dA(:)
!
!    integer :: num_steps, i
!
!    real(dp) :: dt, t, endt
!
!    real(dp), parameter :: volume_mean = 4.2e-6_dp
!    real(dp), parameter :: volume_change = 0.75e-6_dp
!
!    character(len=60) :: sub_name
!
!
!    endt=12.0_dp
!    dt = 0.05_dp  ! Time step
!    num_steps = int(endt / dt) + 1
!
! !--------------------------------------------------------------------------
!
!  sub_name = 'update_surface_area'
!    call enter_exit(sub_name,1)
!
!
!    do i = 1, num_steps !size(volumes)
!      t = (i-1) * dt
!      volumes(i) = volume_mean + volume_change * sin(2.0_dp* pi * frequency * t + (pi / 2.0_dp))
!
!      radii(i) = (3.0_dp * volumes(i) / (4.0_dp * pi)) ** (1.0_dp / 3.0_dp)
!
!      area(i)=4.0_dp * pi * (radii(i)**2.0_dp)
!    end do
!
!
!    dA(1)=0.0_dp !First element is 0
!
!    do i = 2, num_steps
!      dA(i)=(area(i)-area(i-1))/dt;
!    end do
!
!
!    call enter_exit(sub_name,2)
!
!  end subroutine update_surface_area
!##############################################################################
!
!*subname:* calculate surfactant concentration in each alveolar unit

  subroutine update_surfactant_concentration(dt, alv_area_current, alv_dA, surf_concentration)

    real(dp), dimension(:,:), intent(in) :: alv_area_current, alv_dA

    real(dp), dimension(:,:), intent(inout) :: surf_concentration

    integer :: nalv

    real(dp), intent(in) :: dt

    real(dp), parameter :: gamma_star = 0.3e-6_dp !
    real(dp), parameter :: gamma_max = 0.345e-6_dp !
    real(dp), parameter :: bulk_c = 10e-3_dp ! bulk concentration  g/ml
    real(dp), parameter :: k_a = 10.0_dp**4 !6*10**5 ! adsorption coefficient  ml/(g*sec)
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
!*subname:* calculate surface tension for each alveolar unit

  subroutine update_surface_tension(surf_concentration, surface_tension, alv_radii_current, alv_collapse_pressure)

    real(dp), dimension(:,:), intent(in) :: surf_concentration
    real(dp), dimension(:,:), intent(in) :: alv_radii_current

    real(dp), dimension(:,:), intent(inout) :: surface_tension
    real(dp), dimension(:,:), intent(out) :: alv_collapse_pressure

    integer :: nalv

    real(dp) :: dt


    real(dp), parameter :: gamma_star = 0.3e-6_dp !
    real(dp), parameter :: gamma_max = 0.345e-6_dp !
    real(dp), parameter :: surface_tension_hat = 22.02_dp !dyn/cm
    real(dp), parameter :: m2 = 140.0_dp !slope


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


      alv_collapse_pressure(nu_vol,nalv)= ((2.0 *surface_tension(nu_vol,nalv)) /alv_radii_current(nu_vol,nalv))!/10.0

    end do




!    print*,'surface tension2',surface_tension(nu_vol,1)

    call enter_exit(sub_name,2)

  end subroutine update_surface_tension
!##############################################################################
!!
!  subroutine alv_collapse_pressure(alv_radii_current, surface_tension, alv_collapse_pressure)
!
!    !real(dp), intent(in) :: volume_mean, volume_change, gamma_star, gamma_max, bulk_c, k_a, &k_d, sigma_hat, m2
!    real(dp), dimension(:,:), intent(in) :: alv_radii_current, surface_tension
!
!    real(dp), dimension(:,:), intent(out) :: alv_collapse_pressure
!    !real(dp), allocatable :: surface_tension(:)
!    integer :: nalv
!
!    real(dp) :: dt
!
!
!
!    character(len=60) :: sub_name
!
!
!!    endt=12.0_dp
!!    dt = 0.05_dp  ! Time step
!!    num_steps = int(endt / dt) + 1
!
! !--------------------------------------------------------------------------
!
!  sub_name = 'alv_collapse_pressure'
!    call enter_exit(sub_name,1)
!
!
!    !call update_surfactant_concentration
!
!    do nalv = 1,num_units
!        alv_collapse_pressure(nu_vol,nalv)= ((2.0 *surface_tension(nu_vol,nalv)) /alv_radii_current(nu_vol,nalv))/10.0
!    end do
!
!    print*,'pressure2',alv_collapse_pressure(nu_vol,1)
!
!
!    call enter_exit(sub_name,2)
!
!  end subroutine alv_collapse_pressure
!##############################################################################

  !subroutine write_results(volumes,radii,area,dA, &
  !     surf_concentration,surface_tension)

  !  real(dp), allocatable :: volumes(:), radii(:), area(:), dA(:), surf_concentration(:), surface_tension(:)
  !  integer :: num_steps, i
  !  real :: dt, endt,t


  !  character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

  !  sub_name = 'write_results'
  !  call enter_exit(sub_name,1)

    !do i = 1, num_steps
    !   print '("At t = ", F8.4, "  Volume = ", F12.8, "  Radius = ", F12.8, "  Area = ", F12.8, "  dA = ", F12.8,  &
    !    " surf_concentration = ",F12.8, " surface_tension= ", F12.8)', &
    !    (i-1) * 0.005_dp , volumes(i), radii(i), area(i), dA(i), surf_concentration(i), surface_tension(i)

    !end do

  !  num_steps = int(60 / 0.00083 ) + 1
    !num_steps = int(endt / dt) + 1


    !do i = 1, num_steps
    !    print '("At t = ", F8.4, "  Volume = ", F12.8, "  Radius = ", F12.8, "  Area = ", F12.8, "  dA = ", &
    !    F12.8,"  surf_con = ", F12.8,"  surf_ten = ", F12.8)', (i-1) * 0.00083 , volumes(i), radii(i), &
    !    area(i), dA(i), surf_concentration(i), surface_tension(i)

    !end do

  !  write(*,'(2X,''Time'',3X,''Volume'',4X,''Radius'',5X,''Area'',5X,&
  !          &''dA'',4X,''Surf_concentration'',5X,''surface_tension'')')

  !  write(*,'(F7.3,2(F8.1),8(F8.2))') &
  !       (i-1) * 0.005_dp , &  !time, flow, tidal
  !       volumes , & !res (cmH2O/L.s)
  !       radii , & !total model compliance
  !       area, & !Ppl (cmH2O)
  !       dA, & !mean Ptp (cmH2O)
  !       surf_concentration, & !total model volume (L)
  !       surface_tension  !Pmuscle (cmH2O)


  !  call enter_exit(sub_name,2)

  !end subroutine write_results
!##############################################################################

!##############################################################################

end module surfactant