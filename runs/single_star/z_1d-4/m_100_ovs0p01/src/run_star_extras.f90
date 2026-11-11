! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib
  use chem_def
  use num_lib
  use binary_def

  implicit none
  ! VMS wind state uses xtra(20:25), ixtra(20:22), and lxtra(20).
  ! These per-star slots are saved in photos and restored on retries.

  ! these routines are called by the standard run_star check_model
contains


  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


    ! the extras functions in this file will not be called
    ! unless you set their function pointers as done below.
    ! otherwise we use a null_ version which does nothing (except warn).

    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    s% how_many_extra_history_header_items => how_many_extra_history_header_items
    s% data_for_extra_history_header_items => data_for_extra_history_header_items
    s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
    s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
    s% other_wind => VMS_framework !my_wind

  end subroutine extras_controls


  subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

   if (.not. restart) then
   ! depletion event flags: s1 H, He-half, He, C; s2 H, He-half, He, C
   s% lxtra(1)  = .false. ! s1 H depletion  (center_h1  < 1d-4)
   s% lxtra(5)  = .false. ! s1 He-half      (center_he4 < 0.5)
   s% lxtra(6)  = .false. ! s1 He depletion (center_he4 < 1d-4)
   s% lxtra(7)  = .false. ! s1 C depletion  (center_c12 < 1d-4)
   s% xtra(20:25) = 0d0
   s% ixtra(20:22) = 1
   s% lxtra(20) = .true.
   else if (.not. s% lxtra(20)) then
      write(*,'(a)') 'Restart photo has no saved VMS wind state; use a photo made with this fix.'
      ierr = -1
      return
   end if

  end subroutine extras_startup


  integer function extras_start_step(id)
    integer, intent(in) :: id
    integer :: ierr, k
    type (star_info), pointer :: s

    ! depletion thresholds
    real(dp), parameter :: h_dep    = 1d-4
    real(dp), parameter :: he_half  = 0.5d0
    real(dp), parameter :: he_dep   = 1d-4
    real(dp), parameter :: c_dep    = 1d-4

   real(dp) :: gamma1_s1, integral_norm, gamma1_integral

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = 0

   ! Compute gamma integrals fresh every step
   gamma1_s1 = -1d0
   integral_norm = 0d0
   gamma1_integral = 0d0
   do k = 1, s% nz
      integral_norm = integral_norm + (s% Pgas(k) + s% Prad(k)) * s% dm(k) / s% rho(k)
      gamma1_integral = gamma1_integral + &
         (s% gamma1(k) - 4d0/3d0) * (s% Pgas(k) + s% Prad(k)) * s% dm(k) / s% rho(k)
   end do
   gamma1_s1 = gamma1_integral / max(1d-99, integral_norm)

    if (mod(s% model_number, 10) == 0) then
       write(*,*) 'gamma1_integral =', gamma1_s1
    end if

   ! Cache for history columns subroutine

    s% xtra(2) = gamma1_s1 ! average gamma1
    s% xtra(4) = s% gamma1(s% nz) - 4d0/3d0  ! center gamma1

   ! H depletion
   if (.not. s% lxtra(1)) then
      if (s% center_h1 < h_dep) then
         s% lxtra(1) = .true.
         s% xtra(1) = s% r(1)
         write(*,*) 'star H depletion'
         call star_write_model(id, &
            trim(s% log_directory)//'/H_dep.mod', ierr)
         call star_write_profile_info(id, &
            trim(s% log_directory)//'/H_dep.data', ierr)
      end if
   end if

   ! He-half and He depletion: only after H depletion reached
   if (s% lxtra(1)) then

      ! He-half
      if (.not. s% lxtra(5)) then
         if (s% center_he4 < he_half) then
            s% lxtra(5) = .true.
            write(*,*) 'star He half-depletion'
            call star_write_model(id, &
               trim(s% log_directory)//'/He_half.mod', ierr)
            call star_write_profile_info(id, &
               trim(s% log_directory)//'/He_half.data', ierr)
         end if
      end if

      ! He depletion
      if (.not. s% lxtra(6)) then
         if (s% center_he4 < he_dep) then
            s% lxtra(6) = .true.
            write(*,*) 'star He depletion'
            call star_write_model(id, &
               trim(s% log_directory)//'/He_dep.mod', ierr)
            call star_write_profile_info(id, &
               trim(s% log_directory)//'/He_dep.data', ierr)
         end if
      end if

   end if

   ! C depletion: only after H and He depletion reached
   if (s% lxtra(1) .and. s% lxtra(6)) then
      if (.not. s% lxtra(7)) then
         if (s% center_c12 < c_dep) then
            s% lxtra(7) = .true.
            write(*,*) 'star C depletion'
            call star_write_model(id, &
               trim(s% log_directory)//'/C_dep.mod', ierr)
            call star_write_profile_info(id, &
               trim(s% log_directory)//'/C_dep.data', ierr)
         end if
      end if
   end if

  end function extras_start_step


  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    ! Same criterion as star 2 in run_binary_extras: stop once the
    ! pressure-weighted <gamma1 - 4/3> (s% xtra(2), set in extras_start_step)
    ! drops below x_ctrl(1).
    if (s% xtra(2) < s% x_ctrl(1)) then
       write(*,*) 'gamma1_integral below hydro threshold - terminating'
       write(*,*) 'gamma1_integral =', s% xtra(2), ' threshold =', s% x_ctrl(1)
       write(*,*) 'gamma1_center =', s% xtra(4)
       call star_write_model(id, trim(s% log_directory)//'/prehydro.mod', ierr)
       call star_write_profile_info(id, trim(s% log_directory)//'/prehydro.data', ierr)
       extras_check_model = terminate
    end if

    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns =  3
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ! -------------------------------------
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

   ! Gamma integrals are computed fresh every step in extras_binary_finish_step
   ! and cached into b% xtra(2) and b% xtra(3)
    names(1) = 'gamma1_integral'
    vals(1) = s% xtra(2)

    names(2) = 'wind_flag'
    vals(2) = s% xtra(3)

    names(3) = 'gamma1_center'
    vals(3) = s% xtra(4)

  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'
    !names(1) = 'beta'
    !do k = 1, nz
    !   vals(k,1) = s% Pgas(k)/s% P(k)
    !end do

  end subroutine data_for_extra_profile_columns


  integer function how_many_extra_history_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_header_items = 0
  end function how_many_extra_history_header_items


  subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    !integer :: i
    !real(dp) :: Initial_X, Initial_Y, Initial_Z, initial_m

    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

  end subroutine data_for_extra_history_header_items


  integer function how_many_extra_profile_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_header_items = 0
  end function how_many_extra_profile_header_items


  subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    ! here is an example for adding an extra profile header item
    ! also set how_many_extra_profile_header_items
    ! names(1) = 'mixing_length_alpha'
    ! vals(1) = s% mixing_length_alpha

  end subroutine data_for_extra_profile_header_items


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going

    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
  end function extras_finish_step


  subroutine extras_after_evolve(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve

      subroutine VMS_framework(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)

   type (star_info), pointer :: s
   integer, intent(in) :: id
   real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z! surface values (cgs)
   ! NOTE: surface is outermost cell. not necessarily at photosphere.
   ! NOTE: don't assume that vars are set at this point.
   ! so if you want values other than those given as args,
   ! you should use values from s% xh(:,:) and s% xa(:,:) only.
   ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
   real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
   integer, intent(out) :: ierr

! -------------------------------------------------------- Variable decleration ------------------------------------------------------------

   integer :: use_V01, use_VMS 

   real(dp) ::  L1, M1, T1, R1, w_vink                     ! surface properties
   real(dp) ::  Z_init, X_c, X_c_init, Y_c_init, Zsolar, const_k          
   real(dp) ::  vesc, vesc_eff, vinf, vinf_fac
   real(dp) ::  M, Gamma_e_iter, R_iter, eta_vink, eta_switch, vcrit_iter, v_sq_iter, Boost_iter, Gamma_iter 
   real(dp) ::  diff_old, diff_new, w1, w2
   real(dp) ::  vcrit, vsurf, v_vcrit_sq, Boost, alpha, Teff_jump
   real(dp) ::  w_all, w_RSG, w_WR, alfa_WR, mass_loss
   real(dp) :: flag, mdot_check, L_iter, M_iter, Gamma_e_switch, Gamma_e_switch_old
   real(dp) :: Gamma_e, Gamma_e_old, Gamma, Omega_Gamma
   integer :: do_iter, switch_to_V01, switch_to_VMS

   integer :: nz, i
   w = 0
   call star_ptr(id,s,ierr)
   if (ierr /= 0) return

   flag = s% xtra(20)
   mdot_check = s% xtra(21)
   L_iter = s% xtra(22)
   M_iter = s% xtra(23)
   Gamma_e_switch = s% xtra(24)
   Gamma_e = s% xtra(25)
   do_iter = s% ixtra(20)
   switch_to_V01 = s% ixtra(21)
   switch_to_VMS = s% ixtra(22)

   use_V01 = 0
   use_VMS = 0

   Gamma_e_old = Gamma_e
   Gamma_e_switch_old = Gamma_e_switch

   L1 = Lsurf                                                    ! surface L, M, T and R
   M1 = Msurf
   T1 = Tsurf
   R1 = Rsurf

   nz = s% nz
   ! X = s% xa(s%net_iso(ih1),1)                                    ! surface X
   ! Y = s% xa(s%net_iso(ihe4),1)                                   ! surface Y
   ! Z = 1 - (Y + X)                                                ! surface Z
   
   Z_init = s% initial_z                                          ! mdot scaling with Z_fe for Vink recipe
   X_c = s% xa(s%net_iso(ih1),nz)                                 ! central X
   Y_c_init = 0.24 + 2*Z_init
   !Y_c_init = s% initial_y                                       ! if an initial y is given in the inlist file, use this line instead
   X_c_init = 1 - Y_c_init - Z_init
   Gamma_e = 10**(-4.813)*(1+X)*(L1/M1)*(Msun/Lsun)               ! electron scattering Eddington parameter

   Zsolar = 0.019
   const_k = (clight*Msun*1d5)/(Lsun*3600d0*24d0*365d0)           ! constant to evaluate eta


! ----------------------------------------- Store and retrieve switch properties for He burning --------------------------------------------

!   if (s% x_ctrl(2) == 0) then                                  
!         open(unit = 1,file = 'LM_switch.txt',status='replace',form='formatted')       
!         write(1,*) Gamma_e_switch, L_iter, M_iter, mdot_check, do_iter
!      close(1)
!   else if (s% x_ctrl(2) == 1 .and. check_once == 1) then          ! For helium burning, retrieve the L_switch, and do not perform the iteration at start of step
!         open (unit = 2, file = 'LM_switch.txt', status = 'old')
!         read(2,*) Gamma_e_switch, L_iter, M_iter, mdot_check, do_iter
!         close(2)
!         Gamma_e_switch_old = Gamma_e_switch
!         Gamma_e_old = Gamma_e
!         check_once = 0
!   end if



! ------------------------------------------- iteration procedure to estimate L_switch -----------------------------------------------------
            
   call eval_Vink01_wind(w, vinf_fac, Teff_jump, L1, M1)           ! evaluate mdot_vink
   w_vink = w


   if (X_c/X_c_init < 0.999 .and. do_iter == 1) then   

         do i = 1, 3000, 1
            L_iter = 5d0 + 0.001*i
            if (L_iter < 6.5) then
               call eval_homogeneous_mass_low(M)
            else
                  call eval_homogeneous_mass_high(M)
            end if
            L_iter = exp10(L_iter) * Lsun
            M_iter = M * Msun 
            Gamma_e_iter = 10**(-4.813)*(1+X)*(L_iter/M_iter)*(Msun/Lsun)              ! electron scattering gamma_e of iteration to correct M
            R_iter = pow((L_iter)/(4*pi*T1*T1*T1*T1*boltz_sigma), 0.5d0)
            
            vcrit_iter = pow((2d0/3d0)*standard_cgrav*M_iter/R_iter, 0.5d0)/1d5
            vsurf = s% v_rot_avg_surf/1d5 
            v_sq_iter = (vsurf/vcrit_iter)**2

            alpha = 0.52d0
            Gamma_iter = s% photosphere_opacity * L_iter/(4*pi*standard_cgrav*clight*M_iter) 
            Boost_iter = (pow((1-Gamma_iter), 1/alpha - 1))/(pow((1-Gamma_iter - (4d0/9d0)*v_sq_iter), 1/alpha - 1))

            call eval_Vink01_wind(w, vinf_fac, Teff_jump, L_iter, M_iter)
            !mdot_check = w * Boost_iter
            mdot_check = w
            
            vesc = pow(2d0*standard_cgrav*M_iter/R_iter, 0.5d0)/1d5
            vesc_eff = pow(2d0*standard_cgrav*M_iter*(1-Gamma_e_iter)/R_iter, 0.5d0)/1d5
            vinf =  vinf_fac * vesc_eff * pow(Z/Zsolar,0.20d0) 
            
            eta_vink = const_k*(mdot_check*vinf)/(L_iter/Lsun)
            eta_switch = 0.75/(1+(vesc*vesc)/(vinf*vinf)) 
                     
            if (eta_vink > eta_switch) then                                        ! Switch is reached when the etas are equal      
                  exit
            end if               
         end do

         Gamma_e_switch = (10**(-4.813))*(1+X)*(L_iter/M_iter)*(Msun/Lsun)           

   end if

   if (X_c/X_c_init < 0.999) then

         diff_old = Gamma_e_old - Gamma_e_switch_old
         diff_new = Gamma_e - Gamma_e_switch
         
         if (diff_new < 0) then            ! sets what mass loss recipe to use
            use_V01 = 1
         else if (diff_new > 0) then
            use_VMS = 1
            do_iter = 0
         end if
   

         if (diff_old < 0 .and. diff_new > 0) then              ! Additional switch condition to go from V01 to VMS
            switch_to_VMS = 1
            switch_to_V01 = 0
         else if (diff_old > 0 .and. diff_new < 0) then         ! Additional switch condition to go from VMS back to V01 
            switch_to_V01 = 1
            switch_to_VMS = 0
         end if  
   end if
   

! -------------------------------------------------------- Mass loss for different stages --------------------------------------------------

   if (use_V01 == 1 .and. switch_to_V01 == 1) then
         call eval_Vink01_wind(w, vinf_fac, Teff_jump, L1, M1)
         flag = 1
         w_all = w
   else if (use_VMS == 1 .and. switch_to_VMS == 1) then
         call eval_Vink11_wind(w, L1, M1)
         w_all = w
   flag = 2
   end if
   
   call eval_de_Jager_wind(w)                   ! de Jager recipe for below 4,000 K
   w_RSG = w
            
   call eval_sv2020_WR_mdot(w1)                  ! Sander and Vink 2021 recipe for temperatures above 100 kK, (Also has Vink 2017 now, but was not used in the paper)
   call eval_v2017_WR_mdot(w2)  
   w_WR = MAX(w1,w2)  

   


   if (X_c/X_c_init > 0.999) then               ! mass loss during PMS is set to zero
         mass_loss = 0
         flag = 0
   else if (T1 < 9d4 .and. T1 > 4d3) then       ! mass loss is the new VMS implementation in the temp range, regardless of the fuel burning in the core
         mass_loss = w_all
   else if (T1 < 4d3) then
         mass_loss = w_RSG
         flag = 4
   else if (T1 > 1d5) then
         mass_loss = w_WR
         flag = 3
   else if (T1 > 9d4 .and. T1 < 1d5) then
         alfa_WR = (T1 - 9d4) / (1d4)
         mass_loss = alfa_WR * w_WR + (1 - alfa_WR) * w_all
   end if

! -------------------------------------------------------- For rotation only ---------------------------------------------------------------

   vcrit = pow((2d0/3d0)*standard_cgrav*M1/R1, 0.5d0)/1d5
   vsurf = s% v_rot_avg_surf/1d5 
   v_vcrit_sq = (vsurf/vcrit)**2

   alpha = 0.52d0
   Gamma = s% photosphere_opacity * L1/(4*pi*standard_cgrav*clight*M1) 
   !Boost = (pow((1-Gamma), 1/alpha - 1))/(pow((1-Gamma - (4d0/9d0)*v_vcrit_sq), 1/alpha - 1))
   Boost = 1
   Omega_Gamma = Gamma/(1-(4d0/9d0)*v_vcrit_sq)
   mass_loss = Boost * mass_loss

   s% xtra(3) = flag

! -------------------------------------------------------- Outputs for terminal ------------------------------------------------------------

   write(*,*) 'flag = ', flag
   !write(*,*) 'use_V01 = ', use_V01
   !write(*,*) 'use_VMS = ', use_VMS
   !write(*,*) 'switch_to_V01 = ', switch_to_V01
   !write(*,*) 'switch_to_VMS = ', switch_to_VMS
   write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,*) 'L switch = ', log10(L_iter/Lsun)
   write(*,*) 'M switch = ', M_iter/Msun
   !write(*,*) 'mdot switch = ', log10(mdot_check)
   !write(*,*) 'vinf = ', vinf_fac 
   write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,*) 'L model = ', log10(L1/Lsun)
   write(*,*) 'M model = ', M1/Msun
   write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   !write(*,*) 'eta_vink = ', eta_vink
   !write(*,*) 'eta_switch = ', eta_switch
   write(*,*) 'Gamma_e = ', Gamma_e
   write(*,*) 'Gamma_e at switch = ', Gamma_e_switch
   write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,*) 'Rotation boost = ', Boost
   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,*) 'Mdot from V01 = ', log10(w_vink)
   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,*) 'Mdot used in evolution = ', log10(mass_loss)
   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


   ! enhance mass loss below a temperature threshold (x_ctrl(3) = fractional boost, x_ctrl(4) = T threshold in K)
   if (s% x_ctrl(3) /= 0d0 .and. T1 < s% x_ctrl(4)) then
      mass_loss = (1d0 + s% x_ctrl(3)) * mass_loss
      write(*,*) 'cool wind enhancement applied: factor = ', 1d0 + s% x_ctrl(3), ' T = ', T1
   end if

   w = s% Dutch_scaling_factor * mass_loss

   s% xtra(20:25) = [flag, mdot_check, L_iter, M_iter, Gamma_e_switch, Gamma_e]
   s% ixtra(20:22) = [do_iter, switch_to_V01, switch_to_VMS]

   contains


   subroutine eval_homogeneous_mass_low(M)               ! For logL < 6.5
         real(dp), intent(out) :: M
         real(dp) :: f, logL, g

         real(dp), parameter :: F1 = 4.026
         real(dp), parameter :: F2 = 4.277
         real(dp), parameter :: F3 = -1
         real(dp), parameter :: F4 = 25.48
         real(dp), parameter :: F5 = 36.93
         real(dp), parameter :: F6 = -2.792
         real(dp), parameter :: F7 = -3.226
         real(dp), parameter :: F8 = -5.317
         real(dp), parameter :: F9 = 1.648

         !logL = log10(L_iter/Lsun)
         f = F4 + F5*X + F6*pow2(X) + (F7 + F8*X)*L_iter
         g = pow(f, 0.5d0)
         M = (F1 + F2*X + F3*g)/(1+F9*X)
         M = pow(10d0, M)
         
   
   end subroutine eval_homogeneous_mass_low              

   subroutine eval_homogeneous_mass_high(M)              ! For logL > 6.5
         real(dp), intent(out) :: M
         real(dp) :: f, logL, g

         real(dp), parameter :: F1 = 10.05
         real(dp), parameter :: F2 = 8.204
         real(dp), parameter :: F3 = -1
         real(dp), parameter :: F4 = 151.7
         real(dp), parameter :: F5 = 254.5
         real(dp), parameter :: F6 = -11.46
         real(dp), parameter :: F7 = -13.16
         real(dp), parameter :: F8 = -31.68
         real(dp), parameter :: F9 = 2.408

         !logL = log10(L_iter/Lsun)
         f = F4 + F5*X + F6*pow2(X) + (F7 + F8*X)*L_iter
         g = pow(f, 0.5d0)
         M = (F1 + F2*X + F3*g)/(1+F9*X)
         M = pow(10d0, M)          
   
   end subroutine eval_homogeneous_mass_high




   subroutine eval_sv2020_WR_mdot(w)
      real(dp), intent(out) :: w
      real(dp) :: log10w, a, c, gamma_eb, log10w_off
      include 'formats'

      a = 2.932 
      c = -0.44*log10(Z_init/Zsolar) + 9.15
      gamma_eb = -0.324*log10(Z_init/Zsolar) + 0.244
      log10w_off = 0.23*log10(Z_init/Zsolar) - 2.61

      log10w = a*log10(-log10(1-Gamma_e)) - 0.3010*pow(Gamma_eb/Gamma_e, c) + log10w_off
      w = exp10(log10w)

   end subroutine eval_sv2020_WR_mdot


   subroutine eval_v2017_WR_mdot(w)           ! This part of the recipe was not used in the paper
      real(dp), intent(out) :: w
      real(dp) :: log10w
      
      include 'formats'
      
      log10w = -13.3 + 1.36 * log10(L1/Lsun) + 0.61 * log10(Z_init/Zsolar)
      w = exp10(log10w)

   end subroutine eval_v2017_WR_mdot


   subroutine eval_de_Jager_wind(w)
      ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
      real(dp), intent(out) :: w
      real(dp) :: log10w
      include 'formats'
      log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
      w = exp10(log10w)
   end subroutine eval_de_Jager_wind



   subroutine eval_Vink01_wind(w, vinf_fac, Teff_jump, L, M)
      real(dp), intent(inout) :: w, vinf_fac, Teff_jump
      real(dp), intent(in) :: L, M
      real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc_1, vinf_div_vesc_2, vinf_div_vesc, f_z

      f_z = 0.85d0


      ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
      Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_init/Zsolar)))
      dT = 2000d0
      if (T1 > Teff_jump + dT) then
            alfa = 1
      else if (T1 < Teff_jump - dT) then
            alfa = 0
      else
            alfa = (T1 - (Teff_jump - dT)) / (2*dT)
      end if
      

      if (alfa > 0) then ! eval hot side wind (eqn 24)
         vinf_div_vesc_1 = 2.6d0 ! this is the hot side galactic value
         vinf_div_vesc = vinf_div_vesc_1*pow(Z/Zsolar,0.20d0) ! corrected for Z
         logMdot = &
            - 6.697d0 &
            + 2.194d0*log10(L/Lsun/1d5) &
            - 1.313d0*log10(M/Msun/30) &
            - 1.226d0*log10(vinf_div_vesc/2d0) &
            + f_z*log10(Z_init/Zsolar)
         w1 = exp10(logMdot)
      else
         vinf_div_vesc_1 = 0
         w1 = 0
      end if

      if (alfa < 1) then ! eval cool side wind (eqn 25)
         vinf_div_vesc_2 = 1.3d0 ! this is the cool side galactic value
         vinf_div_vesc = vinf_div_vesc_2*pow(Z/Zsolar,0.20d0) ! corrected for Z
         logMdot = &
            - 6.688d0 &
            + 2.210d0*log10(L/Lsun/1d5) &
            - 1.339d0*log10(M/Msun/30) &
            - 1.601d0*log10(vinf_div_vesc/2d0) &
            + f_z*log10(Z_init/Zsolar)
         w2 = exp10(logMdot)
      else
         vinf_div_vesc_2 = 0
         w2 = 0
      end if

      w = alfa*w1 + (1 - alfa)*w2
      vinf_fac = alfa*vinf_div_vesc_1 + (1 - alfa)*vinf_div_vesc_2

   end subroutine eval_Vink01_wind


   subroutine eval_Vink11_wind(w, L, M)
      real(dp), intent(inout) :: w
      real(dp), intent(in) :: L, M
      real(dp) :: logMdot

      logMdot = &
            log10(mdot_check) &        
            + 4.77d0*log10(L/L_iter) & 
            - 3.99d0*log10(M/M_iter)
      w = exp10(logMdot)

   end subroutine eval_Vink11_wind

end subroutine VMS_framework


    subroutine my_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
   use star_def
   type (star_info), pointer :: s
   integer, intent(in) :: id
   real(dp), parameter :: Zbjork = 0.014d0, T_high = 11000, T_low = 10000
   real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
   real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
   integer, intent(out) :: ierr
   real(dp) :: M1, L1, T1, R1, w1, w2, alfa
   w = 0
   ierr = 0
   call star_ptr(id, s, ierr)

   M1 = Msurf
   L1 = Lsurf
   T1 = Tsurf

   if (T1 <= T_low) then
      call eval_lowT(w)
      !print*, "Decin"

   else if (T1 >= T_high) then
        call eval_highT(w)

   else ! transition
      call eval_lowT(w1)
      call eval_highT(w2)
      alfa = (T1 - T_low) / (T_high - T_low)
      w = (1-alfa)*w1 + alfa*w2
   end if

contains

    subroutine eval_de_Jager_wind(w4)
       real(dp), intent(out) :: w4
       real(dp) :: log10w

       log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
       w4 = exp10(log10w)
   end subroutine eval_de_Jager_wind

   subroutine eval_Decin_wind(w)
      ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
      real(dp), intent(out) :: w
      real(dp) :: log10w
      include 'formats'
      !log10w = -17.44 &
      !         +2.80*log10(L1/Lsun) &
      !         -0.13*s% initial_mass! &
      log10w = -20.71 &
               +3.50*log10(L1/Lsun) &
               -0.17*s% initial_mass+0.49! &
      !write(*,*) "Check Decin wind", log10w
      w = 10**(log10w)
   end subroutine eval_Decin_wind

  subroutine eval_lowT(w)
        real(dp), intent(out) :: w
        if (s% x_character_ctrl(2) == 'Decin') then
                call eval_Decin_wind(w)
        end if

        if (s% x_character_ctrl(2) == 'nowind') then
                call eval_Decin_wind(w1)
                w = min(w1, 1d-6)
        end if

   end subroutine eval_lowT

   subroutine eval_highT(w5)
      real(dp), intent(out) :: w5

      if (X < 0.4d0) then
         w5 = 1d-11 * pow(L1/Lsun, 1.29d0) * pow(Y, 1.7d0) * sqrt(Z)
         !print*, "Nugis & Lamers"
      else if (s% x_character_ctrl(1) == 'Bjorklund') then
         call eval_Bjorklund_wind(w5)
      else if (s% x_character_ctrl(1) == 'Vink') then
         call eval_Vink_wind(w5)
      else if (s% x_character_ctrl(1) == 'Sabhahit') then
        call eval_Sabhahit_wind(w5)
      end if
   end subroutine eval_highT

   subroutine eval_Bjorklund_wind(w6)
      real(dp), intent(out) :: w6
      real(dp) :: logw, Meff

      Meff = M1 * (1d0 - 0.34d0 * L1 / (pi4 * clight * s% cgrav(1) * M1)) ! effective mass

      logw = -5.52d0 &
           + 2.39d0 * log10(L1 / (1d6 * Lsun)) &
           - 1.48d0 * log10(Meff / (4.5d1 * Msun)) &
           + 2.12d0 * log10(T1 / 4.5d4) &
           + (0.75d0 - 1.87d0 * log10(T1 / 4.5d4)) * log10(Z / Zbjork)
      w6 = exp10(logw)
   end subroutine eval_Bjorklund_wind

    subroutine eval_Vink_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
    real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula

    ! alfa = 1 for hot side, = 0 for cool side
    if (T1 > 27500d0) then
       alfa = 1
    else if (T1 < 22500d0) then
       alfa = 0
    else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
       Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
       dT = 100d0
       if (T1 > Teff_jump + dT) then
          alfa = 1
       else if (T1 < Teff_jump - dT) then
          alfa = 0
       else
          alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
       end if
    end if

    if (alfa > 0) then ! eval hot side wind (eqn 24)
       vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
       vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
       logMdot = &
          - 6.697d0 &
          + 2.194d0*log10(L1/Lsun/1d5) &
          - 1.313d0*log10(M1/Msun/30d0) &
          - 1.226d0*log10(vinf_div_vesc/2d0) &
          + 0.933d0*log10(T1/4d4) &
          - 10.92d0*pow2(log10(T1/4d4)) &
          + 0.85d0*log10(Z/Zsolar)
       w1 = exp10(logMdot)
    else
       w1 = 0
    end if

    if (alfa < 1) then ! eval cool side wind (eqn 25)
       vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
       vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
       logMdot = &
          - 6.688d0 &
          + 2.210d0*log10(L1/Lsun/1d5) &
          - 1.339d0*log10(M1/Msun/30d0) &
          - 1.601d0*log10(vinf_div_vesc/2d0) &
          + 1.07d0*log10(T1/2d4) &
          + 0.85d0*log10(Z/Zsolar)
       w2 = exp10(logMdot)
    else
       w2 = 0
    end if

    w = alfa*w1 + (1 - alfa)*w2

   ! if (dbg) write(*,*) 'vink wind', w

 end subroutine eval_Vink_wind

subroutine eval_Sabhahit_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
            real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula
            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
               Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
               end if
            end if
            if (alfa > 0) then ! eval hot side wind (eqn 10 in Sabhahit)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z/Zsolar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 9.552 &
                  + 4.77d0*log10(1+X) &
                  + 4.77d0*log10(L1/Lsun/1d5) &
                  - 3.99d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.5d0*log10(Z/Zsolar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if
            w = alfa*w1 + (1 - alfa)*w2
      !      if (dbg) write(*,*) 'Sabhahit wind', w
         end subroutine eval_Sabhahit_wind

end subroutine my_wind
end module run_star_extras
