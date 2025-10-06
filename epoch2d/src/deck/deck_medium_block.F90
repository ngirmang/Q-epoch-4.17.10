! Copyright (C) 2025 Gregory K. Ngirmang
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

MODULE deck_medium_block

#ifdef MEDIUM
  USE strings_advanced

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: medium_deck_initialise, medium_deck_finalise, medium_block_start
  PUBLIC :: medium_block_end, medium_block_handle_element, medium_block_check

  INTEGER :: current_block
  LOGICAL :: dump_rates, dump_keldyshs, next_create_min_set

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: constituents

CONTAINS

  SUBROUTINE medium_deck_initialise
    
    current_block = 0
    IF (deck_state == c_ds_first) THEN
      n_media = 0
    ELSE
      ALLOCATE(constituents(n_media,n_media))
      constituents = -1
    END IF

  END SUBROUTINE medium_deck_initialise



  SUBROUTINE medium_deck_finalise

    INTEGER :: i, j, iu, io, ispecies, n

    IF (n_media == 0) RETURN

    IF (deck_state == c_ds_first) THEN
      ! not going to make a routine for what takes a few lines
      ALLOCATE(media_list(n_media))
      !eps_stored = .TRUE.
      !use_eps_n1n2 = .TRUE.
      RETURN
    END IF

    ! associate compound media
    DO j = 1,n_media
      IF (.NOT. media_list(j)%compound) CYCLE

      IF (ALL(constituents(:,i) == -1)) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'compound media with no consituents'
          END DO
        END IF
        CALL abort_code(c_err_bad_value)
      END IF

      ! count partners
      DO i=1,n_media
        IF(constituents(i,j) == -1) EXIT
        n = constituents(i,j)
        n = species_list(n)%medium_index
        media_list(n)%parent_index = j
        media_list(n)%bound = .TRUE.
      END DO
    END DO

    DEALLOCATE(constituents)

  END SUBROUTINE medium_deck_finalise



  SUBROUTINE medium_block_start
    
    current_block = current_block + 1
    IF (deck_state == c_ds_first) n_media = n_media + 1

    dump_rates = .FALSE. ; dump_keldyshs = .FALSE.
    next_create_min_set = .FALSE.
  END SUBROUTINE medium_block_start



  SUBROUTINE medium_block_end

    INTEGER :: i, species, iq
    

    IF (deck_state == c_ds_first) RETURN

    species = media_list(current_block)%species
    species_list(species)%medium_species = .TRUE.
    species_list(species)%medium_index = current_block
    media_list(current_block)%dump_ionisation_rates = dump_rates .AND. dump_keldyshs
    
    iq = NINT(species_list(species)%charge/q0)
    ! check if this is an electron species
    ! I think allocatables are faster than pointers...idk
        IF ( iq == -1 ) THEN
          media_list(current_block)%is_electron_species = .TRUE.
          ielectron_medium = current_block
          IF (rank == 0) PRINT '("electron medium is medium ", I4)', &
              current_block
        END IF

  END SUBROUTINE medium_block_end



  FUNCTION medium_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    TYPE(primitive_stack) :: stack
    REAL(num) :: v
    INTEGER :: i, species, err, n_constituents, n, j
    REAL(num), DIMENSION(:), POINTER :: r_media_indices

    errcode = c_err_none
    IF (str_cmp(element, 'molecular_alpha') &
        .OR. str_cmp(element, 'molecular_alpha2')) THEN
      eps_stored = .TRUE.
      use_eps_n1n2 = .TRUE.
      use_media_alpha = .TRUE.
      !IF (rank == 0) PRINT '("using molecular n1/n2")'
    END IF

    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN
    ! I think by now, most of the basic species information is here

    ! Since we have no real docs for this model yet, we put comments here instead
    ! Set the name of this model. I don't really use this yet.
    IF (str_cmp(element, 'name')) THEN
      media_list(current_block)%name = value
      RETURN
    END IF

    ! Set the particle species this is related to.
    IF (str_cmp(element, 'species')) THEN
      CALL initialise_stack(stack)
      CALL tokenize(value, stack, errcode)
      v = evaluate(stack, err)
      media_list(current_block)%species = NINT(v)
      CALL deallocate_stack(stack)
      RETURN
    END IF

    ! The maximum lumped density we can have before we produce particles
    ! In effect, we produce macroparticles of this density, essentially.
    IF (str_cmp(element, 'particle_creation_density')) THEN
      media_list(current_block)%particle_create_density = &
           as_real_print(value, element, errcode)
      RETURN
    END IF

    ! By default, we just assume the next ionisation level has a medium.
    ! Otherwise, the medium creates ions only if the delta N is greater
    ! Than this value. Each created macroparticle then has this density,
    ! Unless the created particle count exceeds nmax_create, see below.
    IF (str_cmp(element, 'next_creation_density')) THEN
      media_list(current_block)%next_create_min = &
           as_real_print(value, element, errcode)
      next_create_min_set = .TRUE.
      RETURN
    END IF

    ! To control creation of particles in the next ionisation model,
    ! We produce at maximum this number of macroparticles per time step.
    ! The model when the next ionisation level is not a medium, we would
    ! create macroparticles of next_create_density. However, if
    !          dN / next_create_density > nmax_create
    ! Then instead, each macroparticle gets the following density
    !          dN_per_macro = dN/nmax_create
    ! Thus, we always, at maximum, create nmax_create particles per
    ! ionisation time step.
    IF (str_cmp(element, 'nmax_create')) THEN
      media_list(current_block)%nmax_create = &
           as_integer_print(value, element, errcode)
      RETURN
    END IF

    ! use collisional ionisation for this medium
    IF (str_cmp(element, 'collisional_ionisation')) THEN
      media_list(current_block)%use_collisional_ionisation = &
           as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! use ppt field ionisation for this medium
    IF (str_cmp(element, 'field_ionisation')) THEN
      media_list(current_block)%use_field_ionisation = &
           as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! enable quantised ionisation, only ionise next species
    ! strictly in units of next_creation_density
    IF (str_cmp(element, 'quantised')) THEN
      media_list(current_block)%quantised = &
           as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! dump ionisation rates as a file for this
    IF (str_cmp(element, 'dump_ionisation_rates')) THEN
      media_list(current_block)%ionisation_file = value
      dump_rates = .TRUE.
      RETURN
    END IF

    ! dump ionisation rates as a file for this
    IF (str_cmp(element, 'dump_keldyshs')) THEN
      media_list(current_block)%gamma_file = value
      dump_keldyshs = .TRUE.
      RETURN
    END IF

    ! for predicted ionisation density changes below next_creation_density,
    ! we do probablistic monte carlo ionisation, as in the vanilla code.
    IF (str_cmp(element, 'probablistic_submin')) THEN
      media_list(current_block)%use_prob = &
           as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! this medium is a composite (molecule)
    IF (str_cmp(element, 'compound_medium')) THEN
      media_list(current_block)%compound = &
           as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! set the compound constituents
    IF (str_cmp(element, 'constituents')) THEN
      ALLOCATE(r_media_indices(n_media))
      CALL initialise_stack(stack)
      CALL tokenize(value, stack, errcode)
      CALL evaluate_and_return_all(stack, &
           n_constituents, r_media_indices, errcode)
      CALL deallocate_stack(stack)

      DO i = 1,n_constituents
        n = NINT(r_media_indices(i))
        constituents(i,current_block) = n
      END DO
      DEALLOCATE(r_media_indices)
      RETURN
    END IF

    ! this medium (meant for compound media really) has a "molecular dielectric"
    ! contribution
    IF (str_cmp(element, 'molecular_alpha')) THEN
      media_list(current_block)%contribute_n1 = .TRUE.
      media_list(current_block)%mol_al1 = &
           as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'molecular_alpha2')) THEN
      media_list(current_block)%contribute_n1 = .TRUE.
      media_list(current_block)%mol_al2 = &
           as_real_print(value, element, errcode)
      RETURN
    END IF

    ! here, we do restrict the numbers of ionisation counts on a per collision
    ! basis instead of a per time step basis. This generately creates many many
    ! more particles.
    IF (str_cmp(element, 'per_collision_treatment')) THEN
      media_list(current_block)%per_coll = &
           as_logical_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION medium_block_handle_element



  FUNCTION medium_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none
    IF ( current_block == 0) RETURN

    IF ( media_list(current_block)%quantised .AND. &
         .NOT. next_create_min_set ) THEN
      IF (rank == 0) THEN
        PRINT *, "*** ERROR ***"
        PRINT *, "in medium block, you must set a next_creation_density &
                 &to use quantisation!"
      END IF
      errcode = c_err_bad_value
    END IF

    !IF ( media_list(current_block
  END FUNCTION medium_block_check
#endif
END MODULE deck_medium_block
