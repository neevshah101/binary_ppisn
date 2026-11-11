! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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
module run_binary_extras

  use star_lib
  use star_def
  use const_def
  use const_def
  use chem_def
  use num_lib
  use binary_def
  use math_lib

  implicit none

contains

  subroutine extras_binary_controls(binary_id, ierr)
    integer :: binary_id
    integer, intent(out) :: ierr
    type (binary_info), pointer :: b
    ierr = 0

    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in binary_ptr'
       return
    end if

    ! Set these function pointers to point to the functions you wish to use in
    ! your run_binary_extras. Any which are not set, default to a null_ version
    ! which does nothing.
    b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
    b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
    b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
    b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

    b% extras_binary_startup=> extras_binary_startup
    b% extras_binary_start_step=> extras_binary_start_step
    b% extras_binary_check_model=> extras_binary_check_model
    b% extras_binary_finish_step => extras_binary_finish_step
    b% extras_binary_after_evolve=> extras_binary_after_evolve

    ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
    ! to disable the printed warning message,
    b% warn_binary_extra =.false.

  end subroutine extras_binary_controls

  integer function how_many_extra_binary_history_header_items(binary_id)
    use binary_def, only: binary_info
    integer, intent(in) :: binary_id
    how_many_extra_binary_history_header_items = 0
  end function how_many_extra_binary_history_header_items


  subroutine data_for_extra_binary_history_header_items( &
       binary_id, n, names, vals, ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id, n
    character (len=maxlen_binary_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ierr = 0
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in binary_ptr'
       return
    end if
  end subroutine data_for_extra_binary_history_header_items


  integer function how_many_extra_binary_history_columns(binary_id)
    use binary_def, only: binary_info
    integer, intent(in) :: binary_id
    how_many_extra_binary_history_columns = 0
  end function how_many_extra_binary_history_columns


  subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(in) :: n
    character (len=maxlen_binary_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ierr = 0
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in binary_ptr'
       return
    end if

  end subroutine data_for_extra_binary_history_columns


  integer function extras_binary_startup(binary_id,restart,ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(out) :: ierr
    logical, intent(in) :: restart
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if

    ! ixtra(1): last accepted binary model checked for a stopping condition.
    b% ixtra(1) = -1
    b% ixtra_old(1) = -1

    if (.not. restart) then
       b% lxtra(1) = .false. ! flag for end of donor's main sequence
       b% lxtra(2) = .false. ! flag for beginning RLOF
       b% lxtra(3) = .false. ! flag for end of accretor's main sequence
       ! depletion event flags: s1 H, He-half, He, C; s2 H, He-half, He, C
       b% lxtra(4)  = .false. ! s1 H depletion  (center_h1  < 1d-4)
       b% lxtra(5)  = .false. ! s1 He-half      (center_he4 < 0.5)
       b% lxtra(6)  = .false. ! s1 He depletion (center_he4 < 1d-4)
       b% lxtra(7)  = .false. ! s1 C depletion  (center_c12 < 1d-4)
       b% lxtra(8)  = .false. ! s2 H depletion  (center_h1  < 1d-4)
       b% lxtra(9)  = .false. ! s2 He-half      (center_he4 < 0.5)
       b% lxtra(10) = .false. ! s2 He depletion (center_he4 < 1d-4)
       b% lxtra(11) = .false. ! s2 C depletion  (center_c12 < 1d-4)
       ! initialize some quantitites
       b% xtra(1) = -1d99      ! donor radius at TAMS
       ! extras are used to store the two tidal sychronization timescales (rad/conv) for each star.
       ! -1 if they are point masses
       if (b% point_mass_i /= 1) then
          b% s1% xtra(2) = -1.0d0 ! t_sync_conv_1
          b% s1% xtra(3) = -1.0d0 ! t_sync_rad_1
       end if
       if   (b% point_mass_i /= 2) then
          b% s2% xtra(2) = -1.0d0 ! t_sync_conv_2
          b% s2% xtra(3) = -1.0d0 ! t_sync_rad_2
       end if
       b% xtra(4) =  -1d99    ! radius at onset RLOF
    end if
    extras_binary_startup = keep_going
  end function  extras_binary_startup

  integer function extras_binary_start_step(binary_id,ierr)

    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(out) :: ierr

    extras_binary_start_step = keep_going
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) return

    ! Kolb never switches donor in set_donor_star, so do it here:
    ! the star closer to filling its Roche lobe is the donor.
    if (b% point_mass_i == 0) then
       if (b% d_i == 1 .and. b% rl_relative_gap(2) > b% rl_relative_gap(1)) then
          b% d_i = 2
          b% a_i = 1
          b% s_donor    => b% s2
          b% s_accretor => b% s1
          write(*,*) 'switching donor to star 2 at model', b% model_number
       else if (b% d_i == 2 .and. b% rl_relative_gap(1) > b% rl_relative_gap(2)) then
          b% d_i = 1
          b% a_i = 2
          b% s_donor    => b% s1
          b% s_accretor => b% s2
          write(*,*) 'switching donor to star 1 at model', b% model_number
       end if
    end if

  end function extras_binary_start_step

  !Return either keep_going, retry or terminate
  integer function extras_binary_check_model(binary_id)
    use binary_lib, only: binary_set_separation_eccentricity
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer:: i_don, i_acc
    real(dp) :: r_l2, d_l2
    real(dp) :: q
    integer :: ierr
    character (len=200) :: fname
    real(dp) :: gamma1_s1, gamma1_s2
    logical :: check_primary_prehydro

    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
    extras_binary_check_model = keep_going

    ! Check if s1 has reached prehydro threshold (only while still a live star)
    check_primary_prehydro = .false.
    if (b% have_star_1) check_primary_prehydro = b% s1% x_logical_ctrl(1)
    if (check_primary_prehydro) then
       gamma1_s1 = b% s1% xtra(2)
       if (gamma1_s1 < b% s1% x_ctrl(1)) then

          write(*,*) 'Star 1 gamma1_integral below hydro threshold — converting to point mass'
          write(*,*) 'gamma1_integral s1 =', gamma1_s1, ' threshold =', b% s1% x_ctrl(1)
          write(*,*) 'gamma1_center s1 =', b% s1% xtra(4)

          if (.not. b% job% evolve_both_stars) then
             print *, "Prehydro reached, but evolve_both_stars already false. Terminating."
             call do_saves_for_binary(b, ierr)

             print *, "save models at s1 prehydro"
             fname = trim(b% s1% log_directory)//'/prehydro.mod'
             call star_write_model(b% star_ids(1), fname, ierr)
             fname = trim(b% s1% log_directory)//'/prehydro.data'
             call star_write_profile_info(b% star_ids(1), fname, ierr)
             extras_binary_check_model = terminate
             return
          end if

          call do_saves_for_binary(b, ierr)

          print *, "save models at s1 prehydro"
          fname = trim(b% s1% log_directory)//'/prehydro.mod'
          call star_write_model(b% star_ids(1), fname, ierr)
          fname = trim(b% s1% log_directory)//'/prehydro.data'
          call star_write_profile_info(b% star_ids(1), fname, ierr)

          if (b% have_star_2) then
             fname = 's2_post_s1prehydro.mod'
             call star_write_model(b% star_ids(2), fname, ierr)
          end if

          b% s1% x_logical_ctrl(1) = .false.

          print *, "****************************************"
          print *, "* Switching s1 to point mass (prehydro). *"
          print *, "****************************************"

          b% job% evolve_both_stars = .false.

          b% d_i = 2
          b% a_i = 1
          b% s_donor    => b% s2
          b% s_accretor => b% s1

          b% point_mass_i  = 1
          b% have_star_1   = .false.

          b% s1% generations = 0
          b% s2% generations = 0
          b% generations     = 0

          b% eq_initial_bh_mass = b% m(1)

          b% limit_retention_by_mdot_edd          = .true.
          b% use_radiation_corrected_transfer_rate = .true.

          print *, "If this is a failed SN you can set new e and P here"
          call binary_set_separation_eccentricity(binary_id, 1000000d0*Rsun, 0d0, ierr)
          b% ignore_hard_limits_this_step = .true.

          print *, "----------------------------------------"
          print *, "period, separation after switch"
          print *, b% period, b% separation
          print *, "----------------------------------------"

          ierr = my_binary_finish_step(b)
          if (ierr == keep_going) then
             ierr = 0
          else
             ierr = -1
          end if

       end if
    end if

        ! Check if s2 has reached prehydro threshold
    if (b% have_star_2) then
       gamma1_s2 = b% s2% xtra(2)
       if (gamma1_s2 < b% s2% x_ctrl(1)) then

          write(*,*) 'Star 2 gamma1_integral below hydro threshold — terminating'
          write(*,*) 'gamma1_integral s2 =', gamma1_s2, ' threshold =', b% s2% x_ctrl(1)
          write(*,*) 'gamma1_center s2 =', b% s2% xtra(4)

          call do_saves_for_binary(b, ierr)

          print *, "save models at s2 prehydro"
          fname = trim(b% s2% log_directory)//'/prehydro.mod'
          call star_write_model(b% star_ids(2), fname, ierr)
          fname = trim(b% s2% log_directory)//'/prehydro.data'
          call star_write_profile_info(b% star_ids(2), fname, ierr)

          extras_binary_check_model = terminate
          return

       end if
    end if

    ! Neev: use the Misra et al. (2020) L2 checks throughout the evolution.
    ! Omit the separate Marchant et al. (2016) main-sequence check; see
    ! Misra et al. (2020) for its restriction to massive overcontact binaries.

    if (b% point_mass_i /= 1) then
       i_don = 1
       i_acc = 2
       if (b% m(i_don) > b% m(i_acc)) then !mdon>macc, q<1
          q = b% m(i_acc) / b% m(i_don)
          r_l2 = b% rl(i_don) * (0.784_dp * pow(q,1.05_dp) * exp(-0.188_dp*q) + 1.004_dp)
          d_l2 = b% rl(i_don) * (3.334_dp * pow(q, 0.514_dp) * exp(-0.052_dp*q) + 1.308_dp)
          !Condition to stop when star overflows L2
          if (b% r(i_don) >= (r_l2)) then
             extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1'
          end if
          if (b% r(i_don) >= (d_l2)) then
             extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1'
          end if

       else    !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
          q = b% m(i_acc) / b% m(i_don)
          r_l2 = b% rl(i_don) * (0.29066811_dp * pow(q, 0.82788069_dp) * exp(-0.01572339_dp*q) + 1.36176161_dp)
          d_l2 = b% rl(i_don) * (-0.04029713_dp * pow(q, 0.862143_dp) * exp(-0.04049814_dp*q) + 1.88325644_dp)
          if (b% r(i_don) >= (r_l2)) then
             extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1'
          end if
          if (b% r(i_don) >= (d_l2)) then
             extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1'
          end if
       end if
    end if

    if (b% point_mass_i /= 2) then
       i_don = 2
       i_acc = 1
       if (b% m(i_don) > b% m(i_acc)) then !mdon>macc, q<1
          q = b% m(i_acc) / b% m(i_don)
          r_l2 = b% rl(i_don) * (0.784_dp * pow(q, 1.05_dp) * exp(-0.188_dp * q) + 1.004_dp)
          d_l2 = b% rl(i_don) * (3.334_dp * pow(q,  0.514_dp) * exp(-0.052_dp * q) + 1.308_dp)
          !Condition to stop when star overflows L2
          if (b% r(i_don) >= (r_l2)) then
             extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2'
          end if
          if (b% r(i_don) >= (d_l2)) then
             extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 2'
          end if

       else             !mdonor<maccretor  Condition to stop when mass loss from L2 (previously it was L3) q>1
          q = b% m(i_acc) / b% m(i_don)
          r_l2 = b% rl(i_don) * (0.29066811_dp * pow(q, 0.82788069_dp) * exp(-0.01572339_dp*q) + 1.36176161_dp)
          d_l2 = b% rl(i_don) * (-0.04029713_dp * pow(q, 0.862143_dp) * exp(-0.04049814_dp*q) + 1.88325644_dp)
          if (b% r(i_don) >= (r_l2)) then
             !extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2'
          end if
          if (b% r(i_don) >= (d_l2)) then
             !extras_binary_check_model = terminate
             write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 2'
          end if
       end if
    end if

    if (b% point_mass_i/=0 .and. ((b% rl_relative_gap(1) >= 0.d0) &
         .or. (abs(b% mtransfer_rate/(Msun/secyer)) >= 1.0d-10))) then
       if (b% point_mass_i/=1) then
          i_don = 1
       else
          i_don = 2
       end if
    end if

  end function extras_binary_check_model


  ! returns either keep_going or terminate.
  ! note: cannot request retry; extras_check_model can do that.
  integer function extras_binary_finish_step(binary_id)
    use binary_lib, only: binary_set_separation_eccentricity
    type (binary_info), pointer :: b
    type (star_info), pointer :: s
    integer, intent(in) :: binary_id
    integer :: ierr, k
    character (len=200) :: fname
    real(dp) :: gamma1_s1, gamma1_s2, integral_norm, gamma1_integral

    ! depletion thresholds
    real(dp), parameter :: h_dep    = 1d-4
    real(dp), parameter :: he_half  = 0.5d0
    real(dp), parameter :: he_dep   = 1d-4
    real(dp), parameter :: c_dep    = 1d-4

    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
    extras_binary_finish_step = keep_going

    ! find beginning RLOF
    if (b% lxtra(2) .eqv. .false.) then
       ! RLOF has not started before
       if (b% rl_relative_gap(b% d_i) > 0) then
          if (b% have_star_1) then
             write(fname, fmt="(a19)") 's1_onset_RLOF.mod'
             call star_write_model(b% star_ids(1), fname, ierr)
          end if
          if (b% point_mass_i /= 2) then
             write(fname, fmt="(a19)") 's2_onset_RLOF.mod'
             call star_write_model(b% star_ids(2), fname, ierr)
          end if
          b% lxtra(2) = .true.
          b% xtra(4) = b% s_donor% r(1)
       end if
    end if

    if (b% point_mass_i == 0) then
       ! Check for simultaneous RLOF from both stars after TAMS of one star
       if (b% s2% center_h1 < h_dep .or. b% s1% center_h1 < h_dep) then
          if (b% rl_relative_gap(1) > 0.0_dp .and. b% rl_relative_gap(2) > 0.0_dp) then
             extras_binary_finish_step = terminate
             write(*,'(g0)') "termination code: Both stars fill their Roche Lobe and at least one of them is off MS"
          end if
       end if
    end if

    !check if mass transfer rate reached maximun, assume unstable regime if it happens
    if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
       extras_binary_finish_step = terminate
       write(*,'(g0)') "termination code: Reached maximum mass transfer rate: 1d-1"
    end if

    if (b% model_number == 1 ) then ! Saving initial_models
       write(*,*) "saving initial models"
       if (b% point_mass_i /= 1) then
          call star_write_model(b% s1% id, "initial_star1.mod",  ierr)
       end if
       if (ierr /= 0) return ! failure
       if (b% point_mass_i /= 2) then
          call star_write_model(b% s2% id, "initial_star2.mod",  ierr)
       end if
       if (ierr /= 0) return ! failure
    end if

  end function extras_binary_finish_step

  subroutine extras_binary_after_evolve(binary_id, ierr)
    type (binary_info), pointer :: b
    integer, intent(in) :: binary_id
    integer, intent(out) :: ierr
    call binary_ptr(binary_id, b, ierr)
    if (ierr /= 0) then ! failure in  binary_ptr
       return
    end if
    ! Save final profiles for stellar models that are present.
    ! MANOS reported that star 1 appeared to become a point mass after
    ! minimum-timestep failures, so it was previously saved unconditionally.
    ! TODO: verify these presence checks on that failure path.
    if (b% have_star_1) then
       call star_write_profile_info(b% s1% id, "LOGS1/final_profile.data", ierr)
       if (ierr /= 0) return
    end if
    if (b% have_star_2) then
       call star_write_profile_info(b% s2% id, "LOGS2/final_profile.data", ierr)
    end if
  end subroutine extras_binary_after_evolve

  subroutine do_saves_for_binary(b, ierr)
      type(binary_info), pointer :: b
      integer, intent(out) :: ierr
      integer :: iounit, id
      character (len = strlen) :: str_photo, filename, iomsg, report_str

      call string_for_model_number('x', b% model_number, b% photo_digits, str_photo)

      filename = trim(trim(b% photo_directory) // '/b_' // str_photo)
      report_str = trim('save ' // filename)
      open(newunit = iounit, file = trim(filename), action = 'write', &
         status = 'replace', iostat = ierr, iomsg = iomsg, form = 'unformatted')
      if (ierr /= 0) then
         write(*, *) 'failed in do_saves_for_binary', trim(filename)
         return
      end if
      call binary_photo_write(b% binary_id, iounit)
      close(iounit)

      if (b% have_star_1) then
         filename = trim(trim(b% s1% photo_directory) // '/1_' // str_photo)
         call star_save_for_restart(b% s1% id, filename, ierr)
         report_str = trim(trim(report_str) // ', ' // filename)
      end if
      if (b% have_star_2) then
         filename = trim(trim(b% s2% photo_directory) // '/2_' // str_photo)
         call star_save_for_restart(b% s2% id, filename, ierr)
         report_str = trim(trim(report_str) // ', ' // filename)
      end if
      if (ierr /= 0) then
         write(*, *) 'failed in do_saves_for_binary'
         return
      end if

      write(*, *) trim(trim(report_str) // ' for model'), b% model_number

   end subroutine do_saves_for_binary

   subroutine binary_photo_write(binary_id, iounit)
      integer, intent(in) :: binary_id, iounit
      type(binary_info), pointer :: b

      integer :: ierr, k, len_history_col_spec

      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      write(iounit) star_def_version

      write(iounit, iostat = ierr) &
         b% binary_age, b% binary_age_old, &
         b% model_number, b% model_number_old, &
         b% mtransfer_rate, b% mtransfer_rate_old, &
         b% angular_momentum_j, b% angular_momentum_j_old, &
         b% separation, b% separation_old, &
         b% eccentricity, b% eccentricity_old, &
         b% rl_relative_gap(1), b% rl_relative_gap_old(1), &
         b% rl_relative_gap(2), b% rl_relative_gap_old(2), &
         b% r(1), b% r_old(1), &
         b% r(2), b% r_old(2), &
         b% rl(1), b% rl_old(1), &
         b% rl(2), b% rl_old(2), &
         b% m(1), b% m_old(1), &
         b% m(2), b% m_old(2), &
         b% dt, b% dt_old, &
         b% env(1), b% env_old(1), &
         b% env(2), b% env_old(2), &
         b% eq_initial_bh_mass, &
         b% period, b% period_old, &
         b% max_timestep, b% max_timestep_old, &
         b% change_factor, b% change_factor_old, &
         b% min_binary_separation, &
         b% using_jdot_mb(1), b% using_jdot_mb_old(1), &
         b% using_jdot_mb(2), b% using_jdot_mb_old(2), &
         b% d_i, b% d_i_old, b% a_i, b% a_i_old, &
         b% point_mass_i, b% point_mass_i_old, &
         b% ignore_rlof_flag, b% ignore_rlof_flag_old, &
         b% model_twins_flag, b% model_twins_flag_old, &
         b% dt_why_reason, b% dt_why_reason_old, &
         b% have_star_1, b% have_star_2, &
         b% CE_flag, b% CE_flag_old, &
         b% CE_init, b% CE_init_old, &
         b% CE_nz, b% CE_initial_radius, b% CE_initial_separation, b% CE_initial_Mdonor, &
         b% CE_initial_Maccretor, b% CE_initial_age, b% CE_initial_model_number, &
         b% CE_b_initial_age, b% CE_b_initial_model_number, &
         b% CE_num1, b% CE_num1_old, &
         b% CE_num2, b% CE_num2_old, &
         b% CE_lambda1, b% CE_lambda1_old, &
         b% CE_lambda2, b% CE_lambda2_old, &
         b% CE_Ebind1, b% CE_Ebind1_old, &
         b% CE_Ebind2, b% CE_Ebind2_old, &
         b% CE_years_detached, b% CE_years_detached_old, &
         b% generations, &
         b% ixtra(:), b% ixtra_old(:), &
         b% xtra(:), b% xtra_old(:), &
         b% lxtra(:), b% lxtra_old(:)

      if (associated(b% binary_history_column_spec)) then
         len_history_col_spec = size(b% binary_history_column_spec)
         write(iounit) len_history_col_spec
         write(iounit) b% binary_history_column_spec(1:len_history_col_spec)
      else
         write(iounit) 0 ! len_log_col_spec
      end if
      write(iounit)  &
         b% number_of_binary_history_columns, b% model_number_of_binary_history_values, &
         b% need_to_set_binary_history_names_etc
      if (b% number_of_binary_history_columns > 0) then
         write(iounit) b% binary_history_value_is_integer(1:b% number_of_binary_history_columns)
         do k = 1, b% number_of_binary_history_columns
            write(iounit) b% binary_history_names(k)
         end do
      end if

      if (b% CE_init) then
         write(iounit, iostat = ierr) &
            b% CE_m(:), b% CE_entropy(:), b% CE_U_in(:), b% CE_U_out(:), b% CE_Omega_in(:), b% CE_Omega_out(:)
      end if

      call b% other_binary_photo_write(binary_id, iounit)

      if (ierr /= 0) stop "error in binary_photo_write"

   end subroutine binary_photo_write

      integer function my_binary_finish_step(b)
         type (binary_info), pointer :: b
         real(dp) :: spin_period

         my_binary_finish_step = keep_going
         ! update change factor in case mtransfer_rate has changed
         if(.not. b% doing_first_model_of_run) then
            if(b% mtransfer_rate_old /= b% mtransfer_rate .and. &
               b% mtransfer_rate /= 0 .and. b% mtransfer_rate_old /= 0) then
               if(b% mtransfer_rate < b% mtransfer_rate_old) then
                  b% change_factor = b% change_factor*(1d0-b% implicit_lambda) + b% implicit_lambda* &
                     (1+b% change_factor_fraction*(b% mtransfer_rate/b% mtransfer_rate_old-1))
               else
                  b% change_factor = b% change_factor*(1d0-b% implicit_lambda) + b% implicit_lambda* &
                     (1+b% change_factor_fraction*(b% mtransfer_rate_old/b% mtransfer_rate-1))
               end if
               if(b% change_factor > b% max_change_factor) b% change_factor = b% max_change_factor
               if(b% change_factor < b% min_change_factor) b% change_factor = b% min_change_factor
            end if
         end if

         ! store all variables into "old"

         b% model_number_old = b% model_number
         b% binary_age_old = b% binary_age
         b% mtransfer_rate_old = b% mtransfer_rate
         b% angular_momentum_j_old = b% angular_momentum_j
         b% separation_old = b% separation
         b% eccentricity_old = b% eccentricity
         b% dt_old = b% dt
         b% env_old(1) = b% env(1)
         b% env_old(2) = b% env(2)
         b% period_old = b% period
         b% rl_relative_gap_old(1) = b% rl_relative_gap(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap(2)
         b% r_old(1) = b% r(1)
         b% r_old(2) = b% r(2)
         b% rl_old(1) = b% rl(1)
         b% rl_old(2) = b% rl(2)
         b% m_old(1) = b% m(1)
         b% m_old(2) = b% m(2)
         b% using_jdot_mb_old = b% using_jdot_mb
         b% max_timestep_old = b% max_timestep
         b% change_factor_old = b% change_factor

         b% d_i_old = b% d_i
         b% a_i_old = b% a_i
         b% point_mass_i_old = b% point_mass_i

         b% ignore_rlof_flag_old = b% ignore_rlof_flag
         b% model_twins_flag_old = b% model_twins_flag

         b% CE_flag_old = b% CE_flag
         b% CE_init_old = b% CE_init

         b% CE_num1_old = b% CE_num1
         b% CE_num2_old = b% CE_num2
         b% CE_lambda1_old = b% CE_lambda1
         b% CE_lambda2_old = b% CE_lambda2
         b% CE_Ebind1_old = b% CE_Ebind1
         b% CE_Ebind2_old = b% CE_Ebind2
         b% CE_years_detached_old = b% CE_years_detached

         b% dt_why_reason_old = b% dt_why_reason

         !set all xtra variables
         b% ixtra_old(:) = b% ixtra(:)
         b% xtra_old(:) = b% xtra(:)
         b% lxtra_old(:) = b% lxtra(:)

      end function my_binary_finish_step


end module run_binary_extras
