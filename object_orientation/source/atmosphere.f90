module atmosphere_m 
  use koch_m
  use jsonx_m 
  implicit none 
  
  type atmosphere_t 
    real, allocatable :: wind(:) 
    character(len=:), allocatable :: turb_model, turb_intensity 
    real :: wingspan, hstab_dist, vstab_dist 
    logical :: turb_repeatable

    real :: light_hag(3) = [2000.0, 8000.0, 17000.0]
    real :: light_sig(3) = [5.0, 5.0, 3.0]
    real :: moderate_hag(3) = [2000.0, 11000.0, 45000.0]
    real :: moderate_sig(3) = [10.0, 10.0, 3.0] 
    real :: severe_hag(4) = [2000.0, 4000.0, 20000.0, 80000.0]
    real :: severe_sig(4) = [15.0, 21.0, 21.0, 3.0]

    real, allocatable :: turb_hag(:), turb_sig(:)
    real :: prev_turb(6), prev_xyz(3), prev_f, prev_g 
    real :: Lu, Lv, Lw, Lb
    real, allocatable :: time_history(:,:)
    integer :: time_history_points
  end type atmosphere_t

  contains 
  
  subroutine atmosphere_init(t, j_atmosphere) 
    implicit none 
    type(atmosphere_t) :: t 
    type(json_value), pointer :: j_atmosphere, j_turb, j_sample 
    logical :: found 
    integer :: i, k, n, seed 
    integer, allocatable :: seed_array(:) 

    write(*,*) 'Initializing atmospheric model...' 
    call jsonx_get(j_atmosphere, 'constant_wind[ft/s]', t%wind, 0.0, 3) 

    write(*,*) '   Constant wind[ft/s] = ', t%wind(:) 
    call json_get(j_atmosphere, 'turbulence', j_turb, found) 
    if(found) then 
      call jsonx_get(j_turb, 'model', t%turb_model, 'none') 
      if (t%turb_model .ne. 'none') then 
        call jsonx_get(j_turb, 'wingspan[ft]',       t%wingspan)
        call jsonx_get(j_turb, 'hstab_distance[ft]', t%hstab_dist)
        call jsonx_get(j_turb, 'vstab_distance[ft]', t%vstab_dist)
        call jsonx_get(j_turb, 'intensity',          t%turb_intensity)
        call jsonx_get(j_turb, 'repeatable',         t%turb_repeatable)
        call jsonx_get(j_turb, 'time_history_points',t%time_history_points) 

        write(*,*) '    Turbulence intensity = ', t%turb_intensity
        write(*,*) '    Turbulence model = ', t%turb_model

        ! initialize time history
        allocate(t%time_history(t%time_history_points, 3))
        t%time_history = 0.0 
        do k = 1,t%time_history_points
          t%time_history(k,1) = (k-1) * max(t%hstab_dist, t%vstab_dist)
        end do         

        ! random number generator
        if(t%turb_repeatable) then 
          seed = 12345 
        else 
          call system_clock(count=seed) 
        end if 
        call random_seed(size=n) 
        allocate(seed_array(n)) 
        seed_array = seed + 7 * [(i-1, i=1, n)] ! given from hunsaker 
        call random_seed(put= seed_array) 
        deallocate(seed_array) 

        select case(trim(t%turb_intensity)) 
          case ('light') 
            allocate(t%turb_hag(3)) 
            allocate(t%turb_sig(3)) 
            t%turb_hag = t%light_hag
            t%turb_sig = t%light_sig
          case ('moderate') 
            allocate(t%turb_hag(3)) 
            allocate(t%turb_sig(3)) 
            t%turb_hag = t%moderate_hag
            t%turb_sig = t%moderate_sig
          case ('severe') 
            allocate(t%turb_hag(4)) 
            allocate(t%turb_sig(4)) 
            t%turb_hag = t%severe_hag
            t%turb_sig = t%severe_sig
        end select 

        select case(trim(t%turb_model)) 
        case('dryden_beal') 
          t%Lu = 1750.0 
          t%Lv = 875.0 
          t%Lw = 875.0 
          t%Lb = 4*t%wingspan / pi
        end select 
        
        t%prev_turb(:) = 0.0 
        t%prev_xyz(:) = 0.0 
        t%prev_f = 0.0 
        t%prev_g = 0.0 

        call json_get(j_turb, 'sample', j_sample, found) 
        if(found) call turbulence_sample(t, j_sample) 
      end if 
    end if 

  end subroutine atmosphere_init 

  subroutine turbulence_sample(t, j_sample) 
    implicit none 
    type(atmosphere_t) :: t 
    type(json_value), pointer :: j_sample 
    character(len=:), allocatable :: fn 
    integer :: i, j, k, n, n_psd, iunit, psd_mean_unit 
    real, allocatable :: vals(:,:), psd_mean(:,:), psd_temp(:,:) 
    real :: hag, sigma, dx, turb(6) 
    logical :: found 
    real :: temp_variance, mean_variance

    ! ! Test random number generator 
    ! call test_rand__normal() 

    write(*,*) '    Sampling Atmospheric Turbulence...' 
    call jsonx_get(j_sample, 'save_filename', fn) 
    open(newunit = iunit, file='output_files/' // fn, status='REPLACE') 
    write(*,*) '    Saving sample to ', fn 
    call jsonx_get(j_sample, 'number_of_points',    n) 
    call jsonx_get(j_sample, 'dx[ft]',              dx) 
    call jsonx_get(j_sample, 'height_above_ground[ft]', hag) 

    allocate(vals(n,4)) 

    sigma = interpolate_1D(t%turb_hag, t%turb_sig, hag) 

    write(*,*) '    Altitude[ft] = ', hag
    write(*,*) '    Turbulence Standard Deviation, sigma = ', sigma 

    write(iunit, *) 'distance[ft],uprime[ft/s],vprime[ft/s],wprime[ft/s],pprime[rad/s]'
    do i = 1,n 
      turb(:) = get_turbulence(t, dx, sigma, sigma, sigma) 
      write(iunit,*) dx*real(i-1),',',turb(1),',',turb(2),',',turb(3),',',turb(4)
      vals(i,:) = turb(:) 
    end do 
    close(iunit) 

    call json_get(j_sample, 'psd_analyses', n_psd, found) 
    if(found) then 
      call jsonx_get(j_sample, 'psd_analyses', n_psd) 
      allocate(psd_mean(n/2+1,5))
      allocate(psd_temp(n/2+1,2))
      psd_mean = 0.0 
      psd_temp = 0.0 

      write(*,*) 'Turbulence Normalized Mean PSD Analysis'
      open(newunit=psd_mean_unit, file = 'output_files/' // 'PSD_Mean_Analyses.csv', status='REPLACE')
      write(*,*) '  -saving normalized PSD analysis to PSD_Mean_Analyses.csv'
      do j=1,n_psd 
        write(*,*) 'PSD',j,'of',n_psd
        do i=1,n 
          turb(:) = get_turbulence(t, dx, sigma, sigma, sigma) 
          vals(i,:) = turb(:) 
        end do 

        ! frequency 
        call psd(vals(:,1),dx, psd_norm=psd_temp)
        if (j == 1) psd_mean(:,1) = psd_temp(:,1)

        ! u component
        psd_mean(:,2) = psd_mean(:,2) + psd_temp(:,2)/n_psd

        ! v component
        call psd(vals(:,2),dx, psd_norm=psd_temp)
        psd_mean(:,3) = psd_mean(:,3) + psd_temp(:,2)/n_psd

        ! w component
        call psd(vals(:,3),dx, psd_norm=psd_temp)
        psd_mean(:,4) = psd_mean(:,4) + psd_temp(:,2)/n_psd

        ! p component 
        call psd(vals(:,4),dx, psd_norm=psd_temp, variance=temp_variance)
        psd_mean(:,5) = psd_mean(:,5) + psd_temp(:,2)/n_psd
        mean_variance = mean_variance + temp_variance/n_psd

      end do 

        write(psd_mean_unit,*) 'omega,Su,Sv,Sw,Sp'

        do i = 1, n/2+1
          write(psd_mean_unit,*) psd_mean(i,1)*2*pi, ',',psd_mean(i,2)/(2*pi), ',',psd_mean(i,3)/(2*pi), ',',psd_mean(i,4)/(2*pi),',',psd_mean(i,5)/(2*pi)    
        end do
        close(psd_mean_unit) 
      end if 

      deallocate(vals) 

  end subroutine turbulence_sample 

  function atmosphere_get_turbulence(t, states) result(ans) 
    implicit none 
    type(atmosphere_t) :: t 
    real :: states(21) 
    real :: ans(6) 
    real :: dx, sigma 
     
    dx = sqrt((states(7) - t%prev_xyz(1))**2 + (states(8) - t%prev_xyz(2))**2 + (states(9) - t%prev_xyz(3))**2)
    sigma = interpolate_1D(t%turb_hag, t%turb_sig, -states(9))

    ans = get_turbulence(t, dx, sigma, sigma, sigma) 
    t%prev_xyz(:) = states(7:9) 
  end function atmosphere_get_turbulence

  function get_turbulence(t, dx, su, sv, sw) result(ans) 
    implicit none 
    type(atmosphere_t) :: t 
    real :: dx, su, sv, sw
    real :: ans(6) 

    select case(trim(t%turb_model)) 
    case('dryden_beal') 
      ans(:) = dryden_beal(t,dx,su,sv,sw) 
    end select 
  end function get_turbulence

  function dryden_beal(t, dx, su, sv, sw) result(ans) 
    implicit none
    type(atmosphere_t) :: t 
    real, intent(in) :: dx, su, sv, sw
    real :: ans(6) 
    real :: Au, Av, Aw, Ap, etau, etav, etaw, etap, f, g 
    real :: w_cg, v_cg, w_h, v_h
    integer :: n_hist 

    Au = 0.5  * dx/t%Lu
    Av = 0.25 * dx/t%Lv 
    Aw = 0.25 * dx/t%Lw 
    Ap = 0.5  * dx/t%Lb

    etau = rand_normal() * su * sqrt(2.0 * t%Lu/dx)
    etav = rand_normal() * sv * sqrt(2.0 * t%Lv/dx)
    etaw = rand_normal() * sw * sqrt(2.0 * t%Lw/dx)  
    etap = rand_normal() * sw * sqrt(0.8 * pi * (t%Lw/t%Lb)**(1.0/3.0) / (t%Lw * dx))

    f = ((1.0 - Av) * t%prev_f + 2.0 * Av * etav)/(1.0 + Av) 
    g = ((1.0 - Aw) * t%prev_g + 2.0 * Aw * etaw)/(1.0 + Aw) 

    ans(1) = ((1.0 - Au) * t%prev_turb(1) + 2.0 * Au * etau)/(1.0 + Au)
    ans(2) = ((1.0 - Av) * t%prev_turb(2) + Av*(f + t%prev_f) + sqrt(3.0)*(f - t%prev_f))/(1.0 + Av)
    ans(3) = ((1.0 - Aw) * t%prev_turb(3) + Aw*(g + t%prev_g) + sqrt(3.0)*(g - t%prev_g))/(1.0 + Aw)
    ans(4) = ((1.0 - Ap) * t%prev_turb(4) + 2.0 * Ap * etap)/(1.0 + Ap)

    ! disturbance at cg
    w_cg = ans(3)
    v_cg = ans(2)

    ! Update time history 
    n_hist = t%time_history_points

    write(*,*) 'new call '
    write(*,*) t%time_history(:,1)
    write(*,*) t%time_history(:,2)
    write(*,*) t%time_history(:,3)

    t%time_history(:,1) = t%time_history(:,1) + dx 
    t%time_history(2:n_hist,:) = t%time_history(1:n_hist-1,:)
    t%time_history(1,1) = 0.0
    t%time_history(1,2) = v_cg
    t%time_history(1,3) = w_cg  

    write(*,*)
    write(*,*) t%time_history(:,1)
    write(*,*) t%time_history(:,2)
    write(*,*) t%time_history(:,3)

    ! disturbance at tail 
    w_h = interpolate_1D(t%time_history(:,1), t%time_history(:,3), t%hstab_dist)
    v_h = interpolate_1D(t%time_history(:,1), t%time_history(:,2), t%vstab_dist)

    ans(5) =  (w_cg - w_h)/t%hstab_dist
    ans(6) = -(v_cg - v_h)/t%vstab_dist 

    if (t%hstab_dist > t%time_history(n_hist,1)) write(*,*) 'Time history too small to interpolate hstab.'
    if (t%vstab_dist > t%time_history(n_hist,1)) write(*,*) 'Time history too small to interpolate vstab.'

    t%prev_f = f 
    t%prev_g = g 
    t%prev_turb(:) = ans(:) 
  end function dryden_beal 

end module atmosphere_m 