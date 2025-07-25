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

MODULE deck_fields_block

  USE strings_advanced
  USE simple_io

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: fields_deck_initialise, fields_deck_finalise
  PUBLIC :: fields_block_start, fields_block_end
  PUBLIC :: fields_block_handle_element, fields_block_check

  INTEGER(KIND=MPI_OFFSET_KIND) :: offset

CONTAINS

  SUBROUTINE fields_deck_initialise

  END SUBROUTINE fields_deck_initialise



  SUBROUTINE fields_deck_finalise

  END SUBROUTINE fields_deck_finalise



  SUBROUTINE fields_block_start

    offset = 0

  END SUBROUTINE fields_block_start



  SUBROUTINE fields_block_end

  END SUBROUTINE fields_block_end



  FUNCTION fields_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    CHARACTER(LEN=string_length) :: filename
    INTEGER :: err
    LOGICAL :: got_file

    errcode = c_err_none
#ifndef CONSTEPS
    IF (deck_state == c_ds_first) RETURN
#else
    IF (str_cmp(element, 'eps3')) THEN
      use_eps3 = .TRUE.
      eps_stored = .TRUE.
    END IF
    IF (str_cmp(element, 'eps_n1') &
         .OR. str_cmp(element, 'eps_n2')) THEN
      use_eps_n1n2 = .TRUE.
      eps_stored = .TRUE.
    END IF
    IF (deck_state == c_ds_first) RETURN
#endif
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    CALL get_filename(value, filename, got_file, err)

    IF (str_cmp(element, 'ex')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ex, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_ex)
        CALL evaluate_string_in_space(value, ex, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'ey')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ey, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_ey)
        CALL evaluate_string_in_space(value, ey, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'ez')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ez, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_ez)
        CALL evaluate_string_in_space(value, ez, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'bx')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bx, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_bx)
        CALL evaluate_string_in_space(value, bx, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'by')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, by, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_by)
        CALL evaluate_string_in_space(value, by, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'bz')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bz, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_bz)
        CALL evaluate_string_in_space(value, bz, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF
#ifdef CONSTEPS
    
    IF (eps_stored) THEN
      IF (str_cmp(element, 'eps_x')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, eps0x, offset, errcode)
        ELSE
          CALL initialise_stack(epsx_func)
          CALL tokenize(value, epsx_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, eps0x, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        RETURN
      END IF

      IF (str_cmp(element, 'eps_y')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, eps0y, offset, errcode)
        ELSE
          CALL initialise_stack(epsy_func)
          CALL tokenize(value, epsy_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, eps0y, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        RETURN
      END IF

      IF (str_cmp(element, 'eps_z')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, eps0z, offset, errcode)
        ELSE
          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, eps0z, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        RETURN
      END IF

      IF (str_cmp(element, 'eps')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, eps0x, offset, errcode)
          eps0y = eps0x
          eps0z = eps0x
        ELSE
          CALL initialise_stack(epsx_func)
          CALL initialise_stack(epsy_func)
          CALL initialise_stack(epsz_func)
          CALL tokenize(value, epsx_func, errcode)
          CALL tokenize(value, epsy_func, errcode)
          CALL tokenize(value, epsz_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)

          CALL evaluate_string_in_space(value, eps0x, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
          CALL evaluate_string_in_space(value, eps0y, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
          CALL evaluate_string_in_space(value, eps0z, &
            1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)        
        END IF
        RETURN
      END IF
    ELSE ! continue using iepsx
      IF (str_cmp(element, 'eps_x')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, iepsx, offset, errcode)
        ELSE
          CALL initialise_stack(epsx_func)
          CALL tokenize(value, epsx_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, iepsx, &
               1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        iepsx = 1.0_num / iepsx
        RETURN
      END IF

      IF (str_cmp(element, 'eps_y')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, iepsy, offset, errcode)
        ELSE
          CALL initialise_stack(epsy_func)
          CALL tokenize(value, epsy_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, iepsy, &
               1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        iepsy = 1.0_num / iepsy
        RETURN
      END IF

      IF (str_cmp(element, 'eps_z')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, iepsz, offset, errcode)
        ELSE
          CALL initialise_stack(epsz_func)
          CALL tokenize(value, epsz_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, iepsz, &
               1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        iepsz = 1.0_num / iepsz
        RETURN
      END IF

      IF (str_cmp(element, 'eps')) THEN
        IF (got_file) THEN
          CALL load_single_array_from_file(filename, iepsx, offset, errcode)
          iepsy = iepsx
          iepsz = iepsx
        ELSE
          CALL initialise_stack(epsx_func)
          CALL initialise_stack(epsy_func)
          CALL initialise_stack(epsz_func)
          CALL tokenize(value, epsx_func, errcode)
          CALL tokenize(value, epsy_func, errcode)
          CALL tokenize(value, epsz_func, errcode)

          CALL set_tokenizer_stagger(c_stagger_centre)
          CALL evaluate_string_in_space(value, iepsx, &
               1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
          CALL evaluate_string_in_space(value, iepsy, &
               1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
          CALL evaluate_string_in_space(value, iepsz, &
               1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
        END IF
        iepsx = 1.0_num / iepsx
        iepsy = 1.0_num / iepsy
        iepsz = 1.0_num / iepsz
        RETURN
      END IF
    END IF  !if stored eps

    IF (str_cmp(element, 'eps_n1')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, eps_n1, offset, errcode)
      ELSE
        CALL initialise_stack(eps_n1_func)
        CALL tokenize(value, eps_n1_func, errcode)

        CALL set_tokenizer_stagger(c_stagger_centre)
        CALL evaluate_string_in_space(value, eps_n1, &
          1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'eps_n2')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, eps_n2, offset, errcode)
      ELSE
        CALL initialise_stack(eps_n2_func)
        CALL tokenize(value, eps_n2_func, errcode)

        CALL set_tokenizer_stagger(c_stagger_centre)
        CALL evaluate_string_in_space(value, eps_n2, &
          1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'eps3')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, eps3, offset, errcode)
      ELSE
        CALL initialise_stack(eps3_func)
        CALL tokenize(value, eps3_func, errcode)

        CALL set_tokenizer_stagger(c_stagger_centre)
        CALL evaluate_string_in_space(value, eps3, &
          1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'saturateable_eps3')) THEN
      saturateable_eps3 = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'saturateable_n2')) THEN
      saturateable_n2 = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_eps_spatial_average')) THEN
      use_eps_spatial_average = as_logical_print(value, element, errcode)
      RETURN
    END IF
!
!   IF (str_cmp(element, 'eps_medium_diminish_by_cell')) THEN
!     IF (as_logical_print(value, element, errcode)) &
!         medium_eps_mode = c_epsmode_bycell
!     RETURN
!   END IF
!
!   IF (str_cmp(element, 'eps_medium_diminish_explicit')) THEN
!     IF (as_logical_print(value, element, errcode)) &
!          medium_eps_mode = c_epsmode_explicit
!     RETURN
!   END IF
!end CONSTEPS
#endif
#ifdef GLOBALFIELD
    IF (str_cmp(element, 'global_ex')) THEN
      global_e(1) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'global_ey')) THEN
      global_e(2) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'global_ez')) THEN
      global_e(3) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'global_bx')) THEN
      global_b(1) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'global_by')) THEN
      global_b(2) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'global_bz')) THEN
      global_b(3) = as_real_print(value, element, errcode)
      RETURN
    END IF
#endif

  END FUNCTION fields_block_handle_element



  FUNCTION fields_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

#ifdef CONSTEPS
    IF (rank == 0) THEN
      IF (saturateable_eps3) PRINT "('using saturateable eps3')"
      IF (use_eps_spatial_average) PRINT "('spatially averaging eps')"
      IF (use_eps_n1n2) PRINT "('using n1 and/or n2')"
      IF (saturateable_n2) PRINT "('using saturateable n2')"
    END IF
#endif
  END FUNCTION fields_block_check

END MODULE deck_fields_block
