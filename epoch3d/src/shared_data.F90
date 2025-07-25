! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! ****************************************************************
! All global variables defined here (cf F77 COMMON block).
! ****************************************************************

MODULE shared_data

  USE mpi
  USE constants
  USE sdf_job_info

  IMPLICIT NONE

#ifdef BOUND_HARMONIC
#ifndef CONSTEPS
#error "we now require BOUND_HARMONIC and CONSTEPS to be compiled together"
#endif
#endif

  ! This type represents parameters given to the parser.
  ! It can be extended by a developer freely
  ! It is the responsibility of the developer to ensure that a parameter is
  ! specified when needed

  ! If you set the use_grid_position parameter to .FALSE. then the deck parser
  ! will evaluate position x, y, z as being at the location pack_pos(1,2,3)
  ! rather than x(pack%ix), y(pack%iy), z(pack%iz). It is essential that the
  ! ix, iy, iz, parameters are still set to match, because other functions
  ! will still use them
  TYPE parameter_pack
    LOGICAL :: use_grid_position = .TRUE.
    INTEGER :: pack_ix = 1, pack_iy = 1, pack_iz = 1
    REAL(num), DIMENSION(c_ndims) :: pack_pos = 0.0_num
    REAL(num), DIMENSION(c_ndirs) :: pack_p = 0.0_num
  END TYPE parameter_pack

  TYPE stack_element
    INTEGER :: ptype
    INTEGER :: value
    REAL(num) :: numerical_data
#ifdef PARSER_DEBUG
    CHARACTER(LEN=string_length) :: text
#endif
  END TYPE stack_element

  TYPE primitive_stack
    TYPE(stack_element), POINTER :: entries(:)
    INTEGER :: stack_point, stack_size
    LOGICAL :: init = .FALSE.
    LOGICAL :: is_time_varying = .FALSE.
    LOGICAL :: should_simplify
  END TYPE primitive_stack

  TYPE deck_constant
    CHARACTER(LEN=string_length) :: name
    TYPE(primitive_stack) :: execution_stream
  END TYPE deck_constant

  INTEGER :: n_deck_constants = 0
  TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: deck_constant_list

  !----------------------------------------------------------------------------
  ! string handling
  !----------------------------------------------------------------------------

  CHARACTER(LEN=string_length) :: blank
  TYPE string_type
    CHARACTER(string_length) :: value
  END TYPE string_type
  CHARACTER(LEN=string_length) :: extended_error_string

  !----------------------------------------------------------------------------
  ! Particles
  !----------------------------------------------------------------------------
  
  ! Time to start the particle push - 0 by default, can be set in the control
  ! block of the deck using 'particle_tstart'.
  REAL(num) :: particle_push_start_time = 0.0_num

#ifdef BOUND_HARMONIC

  ! bound_safe_factor sets the allowance of large oscillations of bound species
  REAL(num) :: bound_safe_factor = 0.25 ! 4 cells ought to be enough for anyone
  INTEGER :: check_bound_interval = 0
  ! redo pointer list here due to it being in particles
  ! some day we do this not by ugly hack, consolidate in a type of particle.
  TYPE bp_item
    TYPE(bp_item), POINTER  :: next
    TYPE(particle), POINTER :: p
  END TYPE bp_item
  LOGICAL :: initial_association = .FALSE.
#endif
  ! Object representing a particle
  ! If you add or remove from this section then you *must* update the
  ! particle pack and unpack routines
  TYPE particle
    REAL(num), DIMENSION(3) :: part_p
    REAL(num), DIMENSION(c_ndims) :: part_pos
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    REAL(num) :: weight
#endif
#ifdef DELTAF_METHOD
    REAL(num) :: pvol
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    REAL(num) :: charge
    REAL(num) :: mass
#endif
    TYPE(particle), POINTER :: next, prev
#ifdef PARTICLE_DEBUG
    INTEGER :: processor
    INTEGER :: processor_at_t0
#endif
#ifdef PARTICLE_ID4
    INTEGER :: id
#elif PARTICLE_ID
    INTEGER(i8) :: id
#endif
#ifdef COLLISIONS_TEST
    INTEGER :: coll_count
#endif
#ifdef WORK_DONE_INTEGRATED
    REAL(num) :: work_x
    REAL(num) :: work_y
    REAL(num) :: work_z
    REAL(num) :: work_x_total
    REAL(num) :: work_y_total
    REAL(num) :: work_z_total
#endif
#ifdef PHOTONS
    REAL(num) :: optical_depth
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    REAL(num) :: particle_energy
#endif
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
    REAL(num) :: optical_depth_tri
#endif
#ifdef BREMSSTRAHLUNG
    REAL(num) :: optical_depth_bremsstrahlung
#endif
#ifdef BOUND_HARMONIC
    REAL(num), DIMENSION(3) :: part_ip

    TYPE(bp_item), POINTER :: partners_head
    INTEGER :: partner_count = -1 !number of partners bound to this particle
#endif
  END TYPE particle

  ! Data for migration between species
  TYPE particle_species_migration
    LOGICAL :: this_species, fluid, done
    LOGICAL :: promoteable, demoteable
    INTEGER :: promote_to_species, demote_to_species
    REAL(num) :: promotion_energy_factor, demotion_energy_factor
    REAL(num) :: promotion_density, demotion_density
    REAL(num), DIMENSION(:,:,:), POINTER :: fluid_energy, fluid_density
  END TYPE particle_species_migration

  LOGICAL :: use_particle_migration = .FALSE.
  INTEGER :: particle_migration_interval = 1

  ! Object representing a collection of particles
  ! Used internally by the MPI particle transfer code
  TYPE particle_list
    TYPE(particle), POINTER :: head
    TYPE(particle), POINTER :: tail
    INTEGER(i8) :: count
    INTEGER :: id_update
    ! Pointer is safe if the particles in it are all unambiguously linked
    LOGICAL :: safe

    ! Does this partlist hold copies of particles rather than originals
    LOGICAL :: holds_copies

    TYPE(particle_list), POINTER :: next, prev
  END TYPE particle_list

  ! Represents the initial conditions of a species
  TYPE initial_condition_block
    REAL(num), DIMENSION(:,:,:), POINTER :: density
    REAL(num), DIMENSION(:,:,:,:), POINTER :: temp
    REAL(num), DIMENSION(:,:,:,:), POINTER :: drift

    REAL(num) :: density_min
    REAL(num) :: density_max

    ! Delta-f settings
    REAL(num) :: temp_back(3)
    REAL(num) :: drift_back(3)
    REAL(num) :: density_back
  END TYPE initial_condition_block

  ! Object representing a particle species
  TYPE particle_species
    ! Core properties
    CHARACTER(string_length) :: name
    TYPE(particle_species), POINTER :: next, prev
    INTEGER :: id
    INTEGER :: dumpmask
    INTEGER :: count_update_step

    REAL(num) :: charge
    REAL(num) :: mass
    REAL(num) :: weight
    INTEGER(i8) :: count
    TYPE(particle_list) :: attached_list
    LOGICAL :: immobile
    LOGICAL :: fill_ghosts

    ! Parameters for relativistic and arbitrary particle loader
    INTEGER :: ic_df_type
    REAL(num) :: fractional_tail_cutoff

    TYPE(primitive_stack) :: dist_fn
    TYPE(primitive_stack) :: dist_fn_range(3)

#ifndef NO_TRACER_PARTICLES
    LOGICAL :: zero_current
#endif

#ifdef BREMSSTRAHLUNG
    INTEGER :: atomic_no
    LOGICAL :: atomic_no_set = .FALSE.
#endif
#ifdef BOUND_HARMONIC
    REAL(num), DIMENSION(3) :: harmonic_omega
    REAL(num), DIMENSION(3) :: harmonic_gamma
    REAL(num) :: bfield_sample_factor

    INTEGER, DIMENSION(:), POINTER :: bound_to
    LOGICAL :: dont_transfer_cpu
    REAL(num) :: diminish_factor
#ifdef NONLIN
    REAL(num) :: nl_alpha3 = 0.0_num
    REAL(num) :: linear_factor = 1.0_num
!end NONLIN
#endif
!end BOUND_HARMONIC
#endif
#ifdef CONSTEPS
    LOGICAL :: eps_off_on_ionise = .FALSE.
#endif
#if defined(BOUND_HARMONIC) || defined(MEDIUM)
    LOGICAL :: medium_species = .FALSE.
#endif
#ifdef MEDIUM
    INTEGER :: medium_index = -1
#endif
#ifdef MERGE_PARTICLES
    LOGICAL :: merge = .FALSE.
    INTEGER :: merge_max_particles = 200000000
    REAL(num) :: merge_max_energy_sig = 2.0_num
    REAL(num) :: merge_max_pcomp_sig =  2.0_num
    INTEGER :: merge_start = 0
    REAL(num) :: merge_energy_cut = HUGE(1.0_num)
#endif

    ! Specify if species is background species or not
    LOGICAL :: background_species = .FALSE.
    ! Background density
    REAL(num), DIMENSION(:,:,:), POINTER :: background_density

    ! ID code which identifies if a species is of a special type
    INTEGER :: species_type

    ! particle cell division
    INTEGER(i8) :: global_count
    LOGICAL :: split
    INTEGER(i8) :: npart_max
    ! Secondary list
    TYPE(particle_list), DIMENSION(:,:,:), POINTER :: secondary_list

    ! Loading of particles
    REAL(num) :: npart_per_cell
#ifdef BOUND_HARMONIC
    INTEGER(i8) :: dir_nparts(c_ndims)
#endif

    TYPE(primitive_stack) :: density_function, temperature_function(3)
    TYPE(primitive_stack) :: drift_function(3)

    ! Thermal boundaries
    REAL(num), DIMENSION(:,:,:), POINTER :: ext_temp_x_min, ext_temp_x_max
    REAL(num), DIMENSION(:,:,:), POINTER :: ext_temp_y_min, ext_temp_y_max
    REAL(num), DIMENSION(:,:,:), POINTER :: ext_temp_z_min, ext_temp_z_max

    ! Species_ionisation
    LOGICAL :: electron
    LOGICAL :: ionise
    INTEGER :: ionise_to_species
    INTEGER :: release_species
    INTEGER :: n
    INTEGER :: l
    REAL(num) :: ionisation_energy

    ! Attached probes for this species
#ifndef NO_PARTICLE_PROBES
    TYPE(particle_probe), POINTER :: attached_probes
#endif

    ! Particle migration
    TYPE(particle_species_migration) :: migrate

    ! Initial conditions
    TYPE(initial_condition_block) :: initial_conditions

    ! Per-species boundary conditions
    INTEGER, DIMENSION(2*c_ndims) :: bc_particle
  END TYPE particle_species

  REAL(num), ALLOCATABLE, TARGET :: global_species_density(:,:,:)
  REAL(num), ALLOCATABLE, TARGET :: global_species_temp(:,:,:,:)
  REAL(num), ALLOCATABLE, TARGET :: global_species_drift(:,:,:,:)

  !----------------------------------------------------------------------------
  ! file handling
  !----------------------------------------------------------------------------
  INTEGER :: deck_state
  INTEGER :: subtype_field, subtype_field_r4
  INTEGER :: subarray_field, subarray_field_r4
  INTEGER :: subarray_field_big, subarray_field_big_r4
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp
  INTEGER :: full_dump_every, restart_dump_every
  LOGICAL :: force_first_to_be_restartable
  LOGICAL :: force_final_to_be_restartable
  LOGICAL :: use_offset_grid
  INTEGER :: n_zeros_control, n_zeros = 4
  INTEGER, DIMENSION(num_vars_to_dump) :: dumpmask

  !----------------------------------------------------------------------------
  ! Time averaged IO
  !----------------------------------------------------------------------------
  TYPE averaged_data_block
    REAL(num), DIMENSION(:,:,:,:), POINTER :: array
    REAL(r4), DIMENSION(:,:,:,:), POINTER :: r4array
    REAL(num) :: real_time
    INTEGER :: species_sum, n_species
    LOGICAL :: started, dump_single
  END TYPE averaged_data_block
  LOGICAL :: any_average = .FALSE.

  TYPE io_block_type
    CHARACTER(LEN=string_length) :: name
    REAL(num) :: dt_snapshot, time_prev, time_first
    REAL(num) :: dt_average, dt_min_average, average_time, average_time_start
    REAL(num) :: time_start, time_stop
    REAL(num) :: walltime_interval, walltime_prev
    REAL(num) :: walltime_start, walltime_stop
    REAL(num), POINTER :: dump_at_times(:)
    REAL(num), POINTER :: dump_at_walltimes(:)
    INTEGER, POINTER :: dump_at_nsteps(:)
    INTEGER :: nstep_snapshot, nstep_prev, nstep_first, nstep_average
    INTEGER :: nstep_start, nstep_stop, dump_cycle, prefix_index
    INTEGER :: dump_cycle_first_index
    LOGICAL :: restart, dump, any_average, dump_first, dump_last
    LOGICAL :: dump_source_code, dump_input_decks, rolling_restart
    LOGICAL :: dump_first_after_restart
    LOGICAL :: disabled
    LOGICAL :: use_offset_grid
    INTEGER, DIMENSION(num_vars_to_dump) :: dumpmask
    TYPE(averaged_data_block), DIMENSION(num_vars_to_dump) :: averaged_data
  END TYPE io_block_type

  TYPE(io_block_type), POINTER :: io_block_list(:)
  INTEGER :: n_io_blocks
  LOGICAL :: track_ejected_particles, new_style_io_block
  INTEGER, DIMENSION(num_vars_to_dump) :: averaged_var_block
  INTEGER, DIMENSION(num_vars_to_dump) :: averaged_var_dims
  REAL(num) :: time_start, time_stop
  REAL(num) :: walltime_start, walltime_stop
  INTEGER :: nstep_start, nstep_stop
  CHARACTER(LEN=c_id_length), ALLOCATABLE :: file_prefixes(:)
  INTEGER, ALLOCATABLE :: file_numbers(:)
  INTEGER(i8) :: sdf_buffer_size

  !----------------------------------------------------------------------------
  ! Extended IO information
  !----------------------------------------------------------------------------

  ! Represents a 2 or 3D distribution
  TYPE distribution_function_block
    CHARACTER(LEN=string_length) :: name

    ! The number of dimensions left
    INTEGER :: ndims

    ! The dumpmask for the distribution function
    INTEGER :: dumpmask

    ! The variables which define the ranges and resolutions of the
    ! distribution function
    INTEGER, DIMENSION(c_df_maxdims) :: directions
    REAL(num), DIMENSION(2,c_df_maxdims) :: ranges
    INTEGER, DIMENSION(c_df_maxdims) :: resolution
    LOGICAL, DIMENSION(:), POINTER :: use_species
    REAL(num), DIMENSION(2,c_df_maxdirs) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs) :: use_restrictions

    ! Whether to output deltaf or totalf
    LOGICAL :: output_deltaf

    ! Pointer to next distribution function
    TYPE(distribution_function_block), POINTER :: next
  END TYPE distribution_function_block
  TYPE(distribution_function_block), POINTER :: dist_fns

  ! Represents a data subset
  TYPE subset
    CHARACTER(LEN=string_length) :: name

    ! The dumpmask for the subset
    INTEGER :: mask
    INTEGER, DIMENSION(:,:), POINTER :: dumpmask
    LOGICAL, DIMENSION(:), POINTER :: use_species
    LOGICAL :: use_gamma
    LOGICAL :: use_restriction(c_subset_max)
    LOGICAL :: use_restriction_function(c_subset_max)
    LOGICAL :: space_restrictions
    LOGICAL :: skip, dump_field_grid
    LOGICAL :: time_varying
    REAL(num) :: restriction(c_subset_max)
    TYPE(primitive_stack) :: restriction_function(c_subset_max)
    INTEGER :: subtype, subarray, subtype_r4, subarray_r4
    INTEGER, DIMENSION(c_ndims) :: skip_dir, n_local, n_global, n_start, starts
    ! Persistent subset
    LOGICAL :: persistent, locked
    REAL(num) :: persist_start_time
    INTEGER :: persist_start_step

    ! Pointer to next subset
    TYPE(subset), POINTER :: next
  END TYPE subset
  TYPE(subset), DIMENSION(:), POINTER :: subset_list
  INTEGER :: n_subsets
  LOGICAL :: any_persistent_subset

#ifndef NO_PARTICLE_PROBES
  TYPE particle_probe
    ! Arbitrary point on the plane
    REAL(num), DIMENSION(c_ndims) :: point
    ! The normal to the plane
    REAL(num), DIMENSION(c_ndims) :: normal
    REAL(num) :: ek_min, ek_max
    CHARACTER(LEN=string_length) :: name

    LOGICAL, DIMENSION(:), POINTER :: use_species
    TYPE(particle_list) :: sampled_particles
    TYPE(particle_probe), POINTER :: next
    INTEGER :: dumpmask
  END TYPE particle_probe
#endif

#ifndef NEWPML
  INTEGER :: cpml_thickness
#else
  INTEGER :: cpml_thicknesses(6)
#endif
  INTEGER :: cpml_x_min_start, cpml_x_min_end, cpml_x_min_offset
  INTEGER :: cpml_x_max_start, cpml_x_max_end, cpml_x_max_offset
  INTEGER :: cpml_y_min_start, cpml_y_min_end, cpml_y_min_offset
  INTEGER :: cpml_y_max_start, cpml_y_max_end, cpml_y_max_offset
  INTEGER :: cpml_z_min_start, cpml_z_min_end, cpml_z_min_offset
  INTEGER :: cpml_z_max_start, cpml_z_max_end, cpml_z_max_offset
  ! Indicate that we have a boundary on the current processor
  LOGICAL :: cpml_x_min = .FALSE., cpml_x_max = .FALSE.
  LOGICAL :: cpml_y_min = .FALSE., cpml_y_max = .FALSE.
  LOGICAL :: cpml_z_min = .FALSE., cpml_z_max = .FALSE.
  LOGICAL :: cpml_boundaries
  ! Indicate that the laser injection is located on the current processor
  INTEGER :: cpml_x_min_laser_idx, cpml_x_max_laser_idx
  INTEGER :: cpml_y_min_laser_idx, cpml_y_max_laser_idx
  INTEGER :: cpml_z_min_laser_idx, cpml_z_max_laser_idx
  REAL(num) :: cpml_kappa_max, cpml_a_max, cpml_sigma_max
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_kappa_ex, cpml_kappa_bx
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_kappa_ey, cpml_kappa_by
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_kappa_ez, cpml_kappa_bz
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_ex, cpml_sigma_ex
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_bx, cpml_sigma_bx
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_ey, cpml_sigma_ey
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_by, cpml_sigma_by
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_ez, cpml_sigma_ez
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_bz, cpml_sigma_bz
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: cpml_psi_eyx, cpml_psi_ezx
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: cpml_psi_byx, cpml_psi_bzx
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: cpml_psi_exy, cpml_psi_ezy
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: cpml_psi_bxy, cpml_psi_bzy
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: cpml_psi_exz, cpml_psi_eyz
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: cpml_psi_bxz, cpml_psi_byz
#ifdef CONSTEPS
#ifndef NEWPML
#error  "CONSTEPS must be compiled with NEWPML"
#endif

  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: iepsx, iepsy, iepsz
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: epsx, epsy, epsz

  TYPE(primitive_stack), SAVE :: epsx_func, epsy_func, epsz_func
  TYPE(primitive_stack), SAVE :: eps3_func

  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: eps0x, eps0y, eps0z

  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: eps_n1, eps_n2
  TYPE(primitive_stack), SAVE :: eps_n1_func, eps_n2_func
  LOGICAL :: use_eps_n1n2 = .FALSE., saturateable_n2 = .FALSE.

  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: eps3
  LOGICAL :: use_eps3 = .FALSE., saturateable_eps3 = .FALSE.

  
  LOGICAL :: use_eps_spatial_average = .FALSE.
  REAL(num) :: eps_a1 = 0.75_num
  LOGICAL :: eps_stored = .FALSE.
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: eps_temp
!end CONSTEPS
#endif
#ifdef NEWPML
  REAL(num) pml_thickness_real(6)
  REAL(num) :: newpml_coeff_a = 0.d0, newpml_coeff_m = 0.d0
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: pml_inv, pml_eye, pml_sig
  LOGICAL  use_newpml, use_manualpml
  LOGICAL :: floating_laser=.FALSE.
#endif!NEWPML
#ifdef GLOBALFIELD
  REAL(num) :: global_e(3) = 0.0_num
  REAL(num) :: global_b(3) = 0.0_num
#endif
#ifdef MEDIUM

  ! NEW SHIT: medium model, essentially, my particle-less attempt to model long
  !           pulse propagation
  INTEGER :: n_media = 0 ! number of ioniseable media

  TYPE medium
    CHARACTER(string_length) :: name
    INTEGER :: species = -1

    LOGICAL :: is_electron_species = .FALSE. ! this is an electron medium model

    ! minimum production density to produce a electron macroparticle
    REAL(num) :: particle_create_density = HUGE(1.0_num)
    REAL(num) :: next_create_min = HUGE(1.0_num)

    INTEGER :: nmax_create = 2000000

    LOGICAL :: use_field_ionisation = .FALSE.
    LOGICAL :: use_collisional_ionisation = .FALSE.
    LOGICAL :: quantised = .FALSE.
    LOGICAL :: per_coll = .FALSE.

    CHARACTER(string_length) :: ionisation_file
    CHARACTER(string_length) :: gamma_file
    LOGICAL :: dump_ionisation_rates = .FALSE.
    LOGICAL :: use_prob = .TRUE.

    LOGICAL :: compound = .FALSE.

    LOGICAL :: bound = .FALSE.
    INTEGER :: parent_index = -1

    LOGICAL :: contribute_n1 = .FALSE.
    REAL(num) :: mol_al1 = 0.0_num
    REAL(num) :: mol_al2 = 0.0_num
  END type medium

  TYPE(medium), DIMENSION(:), POINTER :: media_list
  REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: media_density

  INTEGER :: ielectron_medium = -1
  INTEGER :: media_prod_freq = 1
  CHARACTER(LEN=string_length) :: full_ionisation_file
  LOGICAL :: dump_ionisation_file = .FALSE., load_ionisation_file = .FALSE.
#endif

  !----------------------------------------------------------------------------
  ! Core code
  !----------------------------------------------------------------------------
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION
  INTEGER :: realsize

  ! fng is the number of ghost cells needed by the field solver
  ! sng is the number of ghost cells needed by the current smoother
  ! In constants:
  ! ng is the number of ghost cells allocated in the arrays
  ! png is the number of ghost cells needed by the particles
  ! jng is the number of ghost cells needed by the current arrays
  INTEGER :: sng = 1
  INTEGER :: fng, nx, ny, nz
  INTEGER :: nx_global, ny_global, nz_global
  INTEGER(i8) :: npart_global, particles_max_id
  INTEGER :: nsteps, n_species = -1
  LOGICAL :: smooth_currents
  INTEGER :: smooth_its = 1
  INTEGER :: smooth_comp_its = 0
  INTEGER, DIMENSION(:), ALLOCATABLE :: smooth_strides
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: ex, ey, ez, bx, by, bz, jx, jy, jz
  REAL(r4), ALLOCATABLE, DIMENSION(:,:,:) :: r4array
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: npart_per_cell_array
  LOGICAL :: pre_loading

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ex_x_min, ex_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ey_x_min, ey_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ez_x_min, ez_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bx_x_min, bx_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: by_x_min, by_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bz_x_min, bz_x_max

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ex_y_min, ex_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ey_y_min, ey_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ez_y_min, ez_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bx_y_min, bx_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: by_y_min, by_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bz_y_min, bz_y_max

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ex_z_min, ex_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ey_z_min, ey_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ez_z_min, ez_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bx_z_min, bx_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: by_z_min, by_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bz_z_min, bz_z_max

  REAL(num) :: initial_jx, initial_jy, initial_jz

  TYPE(particle_species), DIMENSION(:), POINTER :: species_list
  TYPE(particle_species), DIMENSION(:), POINTER :: ejected_list
  TYPE(particle_species), DIMENSION(:), POINTER :: io_list, io_list_data

  REAL(num), ALLOCATABLE, DIMENSION(:) :: x, xb, y, yb, z, zb

  ! c_max_path_length is the maximum path length on Linux machines
  CHARACTER(LEN=c_max_path_length) :: data_dir
  CHARACTER(LEN=c_max_prefix) :: filesystem

  LOGICAL :: neutral_background = .TRUE.
  LOGICAL :: use_random_seed = .FALSE.
  LOGICAL :: use_particle_lists = .FALSE.
  LOGICAL :: use_particle_count_update = .FALSE.
  LOGICAL :: use_accurate_n_zeros = .FALSE.
  LOGICAL :: use_injectors = .FALSE.
  LOGICAL :: use_more_setup_memory = .FALSE.
  LOGICAL :: all_deck_errcodes_fatal = .FALSE.

  REAL(num) :: dt, t_end, time, dt_multiplier, dt_laser, dt_plasma_frequency
  REAL(num) :: dt_from_restart
  REAL(num) :: dt_min_average, cfl
  ! x_min is the left-hand edge of the simulation domain as specified in
  ! the input deck.
  ! x_grid_min is the location of x(1). Since the grid is cell-centred,
  ! this is usually at x_min + dx/2.
  ! If CPML boundaries are used then the whole grid is shifted along by
  ! cpml_thicknes cells and then x(1) (and also x_grid_min) is at
  ! the location x_min + dx*(1/2-cpml_thickness)
  REAL(num) :: length_x, dx, x_grid_min, x_grid_max, x_min, x_max
  REAL(num) :: x_grid_min_local, x_grid_max_local, x_min_local, x_max_local
  REAL(num) :: length_y, dy, y_grid_min, y_grid_max, y_min, y_max
  REAL(num) :: y_grid_min_local, y_grid_max_local, y_min_local, y_max_local
  REAL(num) :: length_z, dz, z_grid_min, z_grid_max, z_min, z_max
  REAL(num) :: z_grid_min_local, z_grid_max_local, z_min_local, z_max_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_grid_mins, x_grid_maxs
  REAL(num), DIMENSION(:), ALLOCATABLE :: y_grid_mins, y_grid_maxs
  REAL(num), DIMENSION(:), ALLOCATABLE :: z_grid_mins, z_grid_maxs
  REAL(num) :: dir_d(c_ndims), dir_min(c_ndims), dir_max(c_ndims)
  REAL(num) :: dir_grid_min(c_ndims), dir_grid_max(c_ndims)
  REAL(num) :: dir_min_local(c_ndims), dir_max_local(c_ndims)

  LOGICAL :: ic_from_restart = .FALSE.
  LOGICAL :: need_random_state
  LOGICAL :: use_exact_restart, use_exact_restart_set
  LOGICAL :: allow_cpu_reduce
  LOGICAL :: simplify_deck
  LOGICAL :: print_deck_constants
  LOGICAL :: allow_missing_restart
  LOGICAL :: done_mpi_initialise = .FALSE.
  LOGICAL :: use_current_correction
  INTEGER, DIMENSION(2*c_ndims) :: bc_field, bc_particle, bc_allspecies
  INTEGER :: restart_number, step
  CHARACTER(LEN=c_max_path_length) :: full_restart_filename, restart_filename
  CHARACTER(LEN=c_max_path_length) :: status_filename

  TYPE particle_sort_element
    TYPE(particle), POINTER :: particle
  END TYPE particle_sort_element

  TYPE(particle_sort_element), POINTER, DIMENSION(:) :: coll_sort_array
  INTEGER :: coll_sort_array_size = 0

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: coll_pairs
  REAL(num) :: coulomb_log
  LOGICAL :: coulomb_log_auto, use_collisions
  LOGICAL :: use_nanbu = .TRUE.
#ifdef COLL_ELECCHECK
  LOGICAL :: quick_check_elec
  INTEGER :: quick_check_ispecies = -1
#endif

  LOGICAL :: use_field_ionisation, use_collisional_ionisation
  LOGICAL :: use_multiphoton, use_bsi
#ifdef BOUND_HARMONIC
  LOGICAL :: print_ionisation_counts = .FALSE.
  INTEGER, DIMENSION(:), ALLOCATABLE :: field_ionisation_counts
  INTEGER, DIMENSION(:), ALLOCATABLE :: coll_ionisation_counts
  INTEGER, DIMENSION(:), ALLOCATABLE :: last_field_ionisation_counts
  INTEGER, DIMENSION(:), ALLOCATABLE :: last_coll_ionisation_counts
#endif
#ifdef MERGE_PARTICLES
  INTEGER :: merge_nsteps = -1, merge_max_nstep = HUGE(1)
  INTEGER :: merge_scheme = 1
#endif

  INTEGER :: maxwell_solver = c_maxwell_solver_yee
  REAL(num) :: dt_custom

  !----------------------------------------------------------------------------
  ! Moving window
  !----------------------------------------------------------------------------
  LOGICAL :: move_window, inject_particles, window_started
  TYPE(primitive_stack), SAVE :: window_v_x_stack
  LOGICAL :: use_window_stack
  REAL(num) :: window_v_x
  REAL(num) :: window_start_time, window_stop_time
  INTEGER :: bc_x_min_after_move = c_bc_null
  INTEGER :: bc_x_max_after_move = c_bc_null
  INTEGER :: bc_y_min_after_move = c_bc_null
  INTEGER :: bc_y_max_after_move = c_bc_null
  INTEGER :: bc_z_min_after_move = c_bc_null
  INTEGER :: bc_z_max_after_move = c_bc_null
  REAL(num) :: window_offset

#ifdef PHOTONS
  !----------------------------------------------------------------------------
  ! QED - Written by C. P. Ridgers
  !----------------------------------------------------------------------------
  REAL(num), ALLOCATABLE :: log_chi2(:), epsilon_split(:), p_energy(:,:)
  REAL(num), ALLOCATABLE :: log_hsokolov(:,:), p_photon_energy(:,:)
  REAL(num), ALLOCATABLE :: log_eta(:), log_chi(:,:)
  REAL(num), ALLOCATABLE :: log_tpair(:,:), chimin_table(:), log_omegahat(:,:)
  INTEGER :: n_photon, n_pair
  INTEGER :: n_sample_epsilon, n_sample_chi2, n_sample_h
  INTEGER :: n_sample_eta, n_sample_chi, n_sample_t

  ! These track which species should be the species used by the QED routines
  INTEGER :: photon_species = -1, trident_electron_species = -1
  INTEGER :: breit_wheeler_electron_species = -1
  INTEGER :: trident_positron_species = -1, breit_wheeler_positron_species = -1

  REAL(num) :: photon_energy_min = EPSILON(1.0_num)
  REAL(num) :: qed_start_time = 0.0_num
  LOGICAL :: produce_pairs = .FALSE., use_radiation_reaction = .TRUE.
  LOGICAL :: produce_photons = .FALSE., photon_dynamics = .FALSE.
  CHARACTER(LEN=string_length) :: qed_table_location
#endif
  LOGICAL :: use_qed = .FALSE.

#ifdef BREMSSTRAHLUNG
  !----------------------------------------------------------------------------
  ! Bremsstrahlung
  !----------------------------------------------------------------------------
  TYPE interpolation_state
    REAL(num) :: x = HUGE(1.0_num), y = HUGE(1.0_num), val1d, val2d
    INTEGER :: ix1 = 1, ix2 = 1, iy1 = 1, iy2 = 1
  END TYPE interpolation_state
  ! Table declarations
  TYPE brem_tables
    REAL(num), ALLOCATABLE :: cdf_table(:,:), k_table(:,:)
    REAL(num), ALLOCATABLE :: cross_section(:), e_table(:)
    INTEGER :: size_k, size_t
    TYPE(interpolation_state) :: state
  END TYPE brem_tables
  TYPE(brem_tables), ALLOCATABLE :: brem_array(:)
  INTEGER, ALLOCATABLE :: z_values(:)
  INTEGER, ALLOCATABLE :: z_to_index(:)
  INTEGER :: size_brem_array

  ! Bremsstrahlung photon species flag
  INTEGER :: bremsstrahlung_photon_species = -1
#ifndef PHOTONS
  INTEGER :: photon_species = -1
#endif

  ! Deck variables
  REAL(num) :: photon_energy_min_bremsstrahlung = EPSILON(1.0_NUM)
  REAL(num) :: bremsstrahlung_start_time = 0.0_num
  REAL(num) :: photon_weight = 1.0_num
  LOGICAL :: use_bremsstrahlung_recoil = .TRUE.
  LOGICAL :: produce_bremsstrahlung_photons = .FALSE.
  LOGICAL :: bremsstrahlung_photon_dynamics = .FALSE.
  LOGICAL :: use_plasma_screening = .FALSE.
  CHARACTER(LEN=string_length) :: bremsstrahlung_table_location
#endif
  LOGICAL :: use_bremsstrahlung = .FALSE.

  !----------------------------------------------------------------------------
  ! MPI data
  !----------------------------------------------------------------------------
  INTEGER :: coordinates(c_ndims), neighbour(-1:1, -1:1, -1:1)
  INTEGER :: x_coords, proc_x_min, proc_x_max
  INTEGER :: y_coords, proc_y_min, proc_y_max
  INTEGER :: z_coords, proc_z_min, proc_z_max
  INTEGER :: errcode, comm, tag, rank
  INTEGER :: nproc, nprocx, nprocy, nprocz
  INTEGER :: nprocdir(c_ndims)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank, ny_each_rank, nz_each_rank
  INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank
  LOGICAL :: x_min_boundary, x_max_boundary
  LOGICAL :: y_min_boundary, y_max_boundary
  LOGICAL :: z_min_boundary, z_max_boundary
  LOGICAL :: any_open
#ifdef FASTBCS
  LOGICAL :: fastx_bcs = .FALSE.
#endif
  
  !----------------------------------------------------------------------------
  ! domain and loadbalancing
  !----------------------------------------------------------------------------
  LOGICAL :: use_balance, balance_first
  LOGICAL :: use_pre_balance
  LOGICAL :: use_optimal_layout
  REAL(num) :: dlb_threshold
  INTEGER :: dlb_maximum_interval, dlb_force_interval
  INTEGER(i8), PARAMETER :: npart_per_it = 1000000
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_global, y_global, z_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_offset_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: yb_offset_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: zb_offset_global
  ! The location of the processors
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_x_min, cell_x_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_y_min, cell_y_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_z_min, cell_z_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: old_x_max, old_y_max, old_z_max
  INTEGER :: nx_global_min, nx_global_max
  INTEGER :: ny_global_min, ny_global_max
  INTEGER :: nz_global_min, nz_global_max
  INTEGER :: n_global_min(c_ndims), n_global_max(c_ndims)
  LOGICAL :: debug_mode

  !----------------------------------------------------------------------------
  ! Particle injectors
  !----------------------------------------------------------------------------
  TYPE injector_block

    INTEGER :: boundary
    INTEGER :: id
    INTEGER :: species
    REAL(num) :: npart_per_cell
    REAL(num) :: density_min
    REAL(num) :: density_max
    LOGICAL :: use_flux_injector

    TYPE(primitive_stack) :: density_function
    TYPE(primitive_stack) :: temperature_function(3)
    TYPE(primitive_stack) :: drift_function(3)

    REAL(num) :: t_start, t_end
    LOGICAL :: has_t_end
    REAL(num), DIMENSION(:,:), POINTER :: depth

    TYPE(injector_block), POINTER :: next
  END TYPE injector_block

  TYPE(injector_block), POINTER :: injector_x_min, injector_x_max
  TYPE(injector_block), POINTER :: injector_y_min, injector_y_max
  TYPE(injector_block), POINTER :: injector_z_min, injector_z_max

  !----------------------------------------------------------------------------
  ! laser boundaries
  !----------------------------------------------------------------------------
  TYPE laser_block
    ! Boundary to which laser is attached
    INTEGER :: boundary
    ! A unique id number for the laser (not used directly by EPOCH)
    ! Only used if hard coding time profiles
    INTEGER :: id
    REAL(num), DIMENSION(:,:), POINTER :: profile
    REAL(num), DIMENSION(:,:), POINTER :: phase
    REAL(num) :: current_integral_phase

    LOGICAL :: use_time_function, use_phase_function, use_profile_function
    LOGICAL :: use_omega_function
    TYPE(primitive_stack) :: time_function, phase_function, profile_function
    TYPE(primitive_stack) :: omega_function

    REAL(num) :: amp, omega, pol_angle, t_start, t_end
    INTEGER :: omega_func_type

    TYPE(laser_block), POINTER :: next
  END TYPE laser_block

  TYPE(laser_block), POINTER :: laser_x_min, laser_x_max
  TYPE(laser_block), POINTER :: laser_y_min, laser_y_max
  TYPE(laser_block), POINTER :: laser_z_min, laser_z_max
  INTEGER :: n_laser_x_min = 0, n_laser_x_max = 0
  INTEGER :: n_laser_y_min = 0, n_laser_y_max = 0
  INTEGER :: n_laser_z_min = 0, n_laser_z_max = 0
  LOGICAL, DIMENSION(2*c_ndims) :: add_laser = .FALSE.

  TYPE(jobid_type) :: jobid

  INTEGER(i4) :: run_date
  INTEGER(i8) :: defines

  REAL(num) :: walltime_started, real_walltime_start
  REAL(num) :: stop_at_walltime
  REAL(num) :: elapsed_time = 0.0_num
  REAL(num) :: old_elapsed_time = 0.0_num
  INTEGER :: stdout_frequency, check_stop_frequency
  LOGICAL :: check_walltime, print_eta_string, reset_walltime

  LOGICAL, DIMENSION(c_dir_x:c_dir_z,0:c_stagger_max) :: stagger
  INTEGER(i8) :: push_per_field = 5

  ! Absorption diagnostic
  REAL(num) :: laser_inject_local = 0.0_num
  REAL(num) :: laser_absorb_local = 0.0_num
  REAL(num) :: laser_injected = 0.0_num
  REAL(num) :: laser_absorbed = 0.0_num
  LOGICAL :: dump_absorption = .FALSE.

  REAL(num) :: total_field_energy = 0.0_num
  REAL(num) :: total_particle_energy = 0.0_num
  REAL(num), ALLOCATABLE :: total_particle_energy_species(:)

  !----------------------------------------------------------------------------
  ! custom particle loading - written by MP Tooley
  !----------------------------------------------------------------------------

  INTEGER :: n_custom_loaders = 0
  INTEGER, PARAMETER :: c_loader_chunk_size = 131072 ! 1MB/variable data chunks

  TYPE custom_particle_loader
    INTEGER :: species_id

    ! Position Data
    CHARACTER(LEN=string_length) :: x_data
    CHARACTER(LEN=string_length) :: y_data
    CHARACTER(LEN=string_length) :: z_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: x_data_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: y_data_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: z_data_offset

    ! Momentum Data
    CHARACTER(LEN=string_length) :: px_data
    CHARACTER(LEN=string_length) :: py_data
    CHARACTER(LEN=string_length) :: pz_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: px_data_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: py_data_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: pz_data_offset
    LOGICAL :: px_data_given
    LOGICAL :: py_data_given
    LOGICAL :: pz_data_given
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    ! Weight data
    CHARACTER(LEN=string_length) :: w_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: w_data_offset
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    ! ID data
    CHARACTER(LEN=string_length) :: id_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: id_data_offset
    LOGICAL :: id_data_given
    LOGICAL :: id_data_4byte
#endif
  END TYPE custom_particle_loader

  TYPE(custom_particle_loader), DIMENSION(:), POINTER :: custom_loaders_list

END MODULE shared_data
