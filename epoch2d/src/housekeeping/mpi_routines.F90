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

MODULE mpi_routines

  USE helper

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init

    CALL MPI_INIT(errcode)
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    CALL MPI_COMM_SIZE(comm, nproc, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    done_mpi_initialise = .FALSE.
#ifdef MPI_DEBUG
    CALL mpi_set_error_handler
#endif

  END SUBROUTINE mpi_minimal_init



  SUBROUTINE split_domain

    LOGICAL :: reset
    INTEGER :: old_comm
    INTEGER :: ix, iy
    INTEGER :: nxsplit, nysplit
    INTEGER :: area, minarea
    INTEGER :: ranges(3,1), nproc_orig, oldgroup, newgroup
    CHARACTER(LEN=11) :: str
    CHARACTER(LEN=1) :: dir

    nproc_orig = nproc

    IF (nx_global < ncell_min .OR. ny_global < ncell_min) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(ncell_min, str)
        PRINT*,'*** ERROR ***'
        PRINT*,'Simulation domain is too small.'
        PRINT*,'There must be at least ' // TRIM(str) &
            // ' cells in each direction.'
      END IF
      CALL abort_code(c_err_bad_setup)
    END IF

    reset = .FALSE.
    IF (MAX(nprocx,1) * MAX(nprocy,1) > nproc) THEN
      reset = .TRUE.
      IF (rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Requested domain split exceeds CPUs. Ignoring'
      END IF
    ELSE IF (nprocx > 0 .OR. nprocy > 0) THEN
      ! Sanity check
      IF (nprocx == 0) nprocx = nproc / nprocy
      IF (nprocy == 0) nprocy = nproc / nprocx
      nproc = nprocx * nprocy
      nxsplit = nx_global / nprocx
      nysplit = ny_global / nprocy
      IF (nxsplit < ncell_min .OR. nysplit < ncell_min) THEN
        reset = .TRUE.
        IF (rank == 0) THEN
          IF (nxsplit < ncell_min) THEN
            dir = 'x'
          ELSE IF (nysplit < ncell_min) THEN
            dir = 'y'
          END IF
          PRINT*,'*** WARNING ***'
          PRINT'('' Requested domain split gives less than '', I1, &
              &  '' cells in the '', A, ''-direction. Ignoring'')', &
              ncell_min, dir
        END IF
      END IF
    END IF

    IF (reset) THEN
      IF (rank == 0) THEN
        PRINT *, 'Unable to use requested processor subdivision. Using ' &
            // 'default division.'
      END IF
      nprocx = 0
      nprocy = 0
    END IF

    IF (nprocx * nprocy == 0) THEN
      DO WHILE (nproc > 1)
        ! Find the processor split which minimizes surface area of
        ! the resulting domain

        minarea = nx_global + ny_global

        DO ix = 1, nproc
          iy = nproc / ix
          IF (ix * iy /= nproc) CYCLE

          nxsplit = nx_global / ix
          nysplit = ny_global / iy
          ! Actual domain must be bigger than the number of ghostcells
          IF (nxsplit < ncell_min .OR. nysplit < ncell_min) CYCLE

          area = nxsplit + nysplit
          IF (area < minarea) THEN
            nprocx = ix
            nprocy = iy
            minarea = area
          END IF
        END DO

        IF (nprocx > 0) EXIT

        ! If we get here then no suitable split could be found. Decrease the
        ! number of processors and try again.

        nproc = nproc - 1
      END DO
    END IF

    IF (nproc_orig /= nproc) THEN
      IF (.NOT.allow_cpu_reduce) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(nproc, str)
          PRINT*,'*** ERROR ***'
          PRINT*,'Cannot split the domain using the requested number of CPUs.'
          PRINT*,'Try reducing the number of CPUs to ',TRIM(str)
        END IF
        CALL abort_code(c_err_bad_setup)
        STOP
      END IF
      IF (rank == 0) THEN
        CALL integer_as_string(nproc, str)
        PRINT*,'*** WARNING ***'
        PRINT*,'Cannot split the domain using the requested number of CPUs.'
        PRINT*,'Reducing the number of CPUs to ',TRIM(str)
      END IF
      ranges(1,1) = nproc
      ranges(2,1) = nproc_orig - 1
      ranges(3,1) = 1
      old_comm = comm
      CALL MPI_COMM_GROUP(old_comm, oldgroup, errcode)
      CALL MPI_GROUP_RANGE_EXCL(oldgroup, 1, ranges, newgroup, errcode)
      CALL MPI_COMM_CREATE(old_comm, newgroup, comm, errcode)
      IF (comm == MPI_COMM_NULL) THEN
        CALL MPI_FINALIZE(errcode)
        STOP
      END IF
      CALL MPI_GROUP_FREE(oldgroup, errcode)
      CALL MPI_GROUP_FREE(newgroup, errcode)
      CALL MPI_COMM_FREE(old_comm, errcode)
    END IF

    CALL setup_communicator

  END SUBROUTINE split_domain



  SUBROUTINE setup_communicator

    INTEGER :: dims(c_ndims), idim, old_comm
    LOGICAL :: periods(c_ndims), reorder, op
    INTEGER :: test_coords(c_ndims)
    INTEGER :: ix, iy

    dims = (/nprocy, nprocx/)
    CALL MPI_DIMS_CREATE(nproc, c_ndims, dims, errcode)

    periods = .FALSE.
    reorder = .TRUE.

    ! Set boundary to be periodic if *any* boundary condition requires it.
    ! For per-species boundary conditions then this will be true
    ! if any of the species are periodic

    IF (bc_field(c_bd_x_min) == c_bc_periodic &
        .OR. bc_x_min_after_move == c_bc_periodic) THEN
      periods(c_ndims) = .TRUE.
    ELSE
      DO idim = 1, n_species
        IF (species_list(idim)%bc_particle(c_bd_x_min) == c_bc_periodic) THEN
          periods(c_ndims) = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    IF (bc_field(c_bd_y_min) == c_bc_periodic &
        .OR. bc_y_min_after_move == c_bc_periodic) THEN
      periods(c_ndims-1) = .TRUE.
    ELSE
      DO idim = 1, n_species
        IF (species_list(idim)%bc_particle(c_bd_y_min) == c_bc_periodic) THEN
          periods(c_ndims-1) = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    old_comm = comm
    CALL MPI_CART_CREATE(old_comm, c_ndims, dims, periods, reorder, comm, &
                         errcode)
    CALL MPI_COMM_FREE(old_comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, c_ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_y_min, proc_y_max, errcode)

    nprocx = dims(2)
    nprocy = dims(1)
    nprocdir = dims

    IF (rank == 0) THEN
      IF (.NOT.use_pre_balance .OR. .NOT.use_optimal_layout) THEN
        PRINT*, 'Processor subdivision is ', (/nprocx, nprocy/)
      END IF
    END IF

    x_coords = coordinates(c_ndims)
    x_min_boundary = .FALSE.
    x_max_boundary = .FALSE.
    IF (x_coords == 0) x_min_boundary = .TRUE.
    IF (x_coords == nprocx - 1) x_max_boundary = .TRUE.

    y_coords = coordinates(c_ndims-1)
    y_min_boundary = .FALSE.
    y_max_boundary = .FALSE.
    IF (y_coords == 0) y_min_boundary = .TRUE.
    IF (y_coords == nprocy - 1) y_max_boundary = .TRUE.

    neighbour = MPI_PROC_NULL
    DO iy = -1, 1
      DO ix = -1, 1
        test_coords = coordinates
        test_coords(1) = test_coords(1)+iy
        test_coords(2) = test_coords(2)+ix
        op = .TRUE.

        ! For some stupid reason MPI_CART_RANK returns an error rather than
        ! MPI_PROC_NULL if the coords are out of range.
        DO idim = 1, c_ndims
          IF ((test_coords(idim) < 0 &
              .OR. test_coords(idim) >= dims(idim)) &
              .AND. .NOT. periods(idim)) op = .FALSE.
        END DO
        IF (op) CALL MPI_CART_RANK(comm, test_coords, neighbour(ix,iy), errcode)
      END DO
    END DO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: ispecies, idim
    INTEGER :: nx0, nxp
    INTEGER :: ny0, nyp

#ifdef NEWPML
    IF (.NOT.(cpml_boundaries .OR. use_newpml)) &
         cpml_thicknesses = 0
#else
    IF (.NOT.cpml_boundaries) cpml_thickness = 0
#endif

    CALL split_domain

    ALLOCATE(npart_each_rank(nproc))
    ALLOCATE(x_grid_mins(0:nprocx-1), x_grid_maxs(0:nprocx-1))
    ALLOCATE(y_grid_mins(0:nprocy-1), y_grid_maxs(0:nprocy-1))
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))
    ALLOCATE(cell_y_min(nprocy), cell_y_max(nprocy))

#ifndef NEWPML
    nx_global = nx_global + 2 * cpml_thickness
    ny_global = ny_global + 2 * cpml_thickness
#else
    nx_global = nx_global + cpml_thicknesses(1) + cpml_thicknesses(2)
    ny_global = ny_global + cpml_thicknesses(3) + cpml_thicknesses(4)
#endif

    IF (use_exact_restart) THEN
      old_x_max(nprocx) = nx_global
      cell_x_max = old_x_max
      DEALLOCATE(old_x_max)

      old_y_max(nprocy) = ny_global
      cell_y_max = old_y_max
      DEALLOCATE(old_y_max)

      cell_x_min(1) = 1
      DO idim = 2, nprocx
        cell_x_min(idim) = cell_x_max(idim-1) + 1
      END DO

      cell_y_min(1) = 1
      DO idim = 2, nprocy
        cell_y_min(idim) = cell_y_max(idim-1) + 1
      END DO
    ELSE
      nx0 = nx_global / nprocx
      ny0 = ny_global / nprocy

      ! If the number of gridpoints cannot be exactly subdivided then fix
      ! The first nxp processors have nx0 grid points
      ! The remaining processors have nx0+1 grid points
      IF (nx0 * nprocx /= nx_global) THEN
        nxp = (nx0 + 1) * nprocx - nx_global
      ELSE
        nxp = nprocx
      END IF

      IF (ny0 * nprocy /= ny_global) THEN
        nyp = (ny0 + 1) * nprocy - ny_global
      ELSE
        nyp = nprocy
      END IF

      DO idim = 1, nxp
        cell_x_min(idim) = (idim - 1) * nx0 + 1
        cell_x_max(idim) = idim * nx0
      END DO
      DO idim = nxp + 1, nprocx
        cell_x_min(idim) = nxp * nx0 + (idim - nxp - 1) * (nx0 + 1) + 1
        cell_x_max(idim) = nxp * nx0 + (idim - nxp) * (nx0 + 1)
      END DO

      DO idim = 1, nyp
        cell_y_min(idim) = (idim - 1) * ny0 + 1
        cell_y_max(idim) = idim * ny0
      END DO
      DO idim = nyp + 1, nprocy
        cell_y_min(idim) = nyp * ny0 + (idim - nyp - 1) * (ny0 + 1) + 1
        cell_y_max(idim) = nyp * ny0 + (idim - nyp) * (ny0 + 1)
      END DO
    END IF

    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)
    n_global_min(1) = nx_global_min
    n_global_max(1) = nx_global_max

    ny_global_min = cell_y_min(y_coords+1)
    ny_global_max = cell_y_max(y_coords+1)
    n_global_min(2) = ny_global_min
    n_global_max(2) = ny_global_max

    nx = nx_global_max - nx_global_min + 1
    ny = ny_global_max - ny_global_min + 1

    subtype_field = 0

    DEALLOCATE(x, y)
    DEALLOCATE(xb, yb)
    ALLOCATE(x(1-ng:nx+ng), y(1-ng:ny+ng))
    ALLOCATE(xb(1-ng:nx+ng), yb(1-ng:ny+ng))
    ALLOCATE(x_global(1-ng:nx_global+ng))
    ALLOCATE(y_global(1-ng:ny_global+ng))
    ALLOCATE(xb_global(1-ng:nx_global+ng))
    ALLOCATE(yb_global(1-ng:ny_global+ng))
    ALLOCATE(xb_offset_global(1-ng:nx_global+ng))
    ALLOCATE(yb_offset_global(1-ng:ny_global+ng))
    ALLOCATE(ex(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(ey(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(ez(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(bx(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(by(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(bz(1-ng:nx+ng, 1-ng:ny+ng))
    ! Current may need an extra layer of ghostcells.
    ALLOCATE(jx(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jy(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jz(1-jng:nx+jng, 1-jng:ny+jng))
#if defined(CONSTEPS) || defined(MEDIUM)
    IF (.NOT. eps_stored) THEN
      ALLOCATE(iepsx(1-ng:nx+ng, 1-ng:ny+ng))
      ALLOCATE(iepsy(1-ng:nx+ng, 1-ng:ny+ng))
      ALLOCATE(iepsz(1-ng:nx+ng, 1-ng:ny+ng))
    ELSE
      ALLOCATE(epsx(1-ng:nx+ng, 1-ng:ny+ng))
      ALLOCATE(epsy(1-ng:nx+ng, 1-ng:ny+ng))
      ALLOCATE(epsz(1-ng:nx+ng, 1-ng:ny+ng))
      IF (use_eps_spatial_average) &
        ALLOCATE(eps_temp(1-ng:nx+ng, 1-ng:ny+ng))
      IF ( use_eps_n1n2 ) THEN
        ALLOCATE(eps_n1(1-ng:nx+ng, 1-ng:ny+ng))
        ALLOCATE(eps_n2(1-ng:nx+ng, 1-ng:ny+ng))
      ELSE IF ( use_eps3 ) THEN
        ALLOCATE(eps0x(1-ng:nx+ng, 1-ng:ny+ng))
        ALLOCATE(eps0y(1-ng:nx+ng, 1-ng:ny+ng))
        ALLOCATE(eps0z(1-ng:nx+ng, 1-ng:ny+ng))
        ALLOCATE(eps3(1-ng:nx+ng, 1-ng:ny+ng))
      END IF
    END IF
#ifdef MEDIUM
    IF (n_media > 0) &
      ALLOCATE(media_density(1-ng:nx+ng, 1-ng:ny+ng, n_media))
#endif
!end CONSTEPS or MEDIUM
#endif

    ! Setup the particle lists
    IF (n_species > 0) &
        NULLIFY(species_list(1)%prev, species_list(n_species)%next)
    DO ispecies = 1, n_species-1
      species_list(ispecies)%next => species_list(ispecies+1)
    END DO
    DO ispecies = 2, n_species
      species_list(ispecies)%prev => species_list(ispecies-1)
    END DO
    DO ispecies = 1, n_species
      species_list(ispecies)%id = ispecies
#ifndef NO_PARTICLE_PROBES
      NULLIFY(species_list(ispecies)%attached_probes)
#endif
      NULLIFY(species_list(ispecies)%attached_list%next)
      NULLIFY(species_list(ispecies)%attached_list%prev)
      CALL create_empty_partlist(species_list(ispecies)%attached_list)

      IF (species_list(ispecies)%bc_particle(c_bd_x_min) == c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_x_min(1-ng:ny+ng,1:3))
      END IF
      IF (species_list(ispecies)%bc_particle(c_bd_x_max) == c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_x_max(1-ng:ny+ng,1:3))
      END IF
      IF (species_list(ispecies)%bc_particle(c_bd_y_min) == c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_y_min(1-ng:nx+ng,1:3))
      END IF
      IF (species_list(ispecies)%bc_particle(c_bd_y_max) == c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_y_max(1-ng:nx+ng,1:3))
      END IF
    END DO

    ALLOCATE(total_particle_energy_species(n_species))

    CALL allocate_ic

    start_time = MPI_WTIME()
    done_mpi_initialise = .TRUE.

  END SUBROUTINE mpi_initialise



  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

    IF (rank == 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
#ifndef NO_IO
      WRITE(stat_unit, *)
      WRITE(stat_unit, '(''runtime = '', i4, ''h '', i2, ''m '', i2, &
          & ''s on '', i4, '' process elements.'')') hours, minutes, seconds, &
          nproc
#endif
    END IF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE mpi_close



#ifdef MPI_DEBUG
  SUBROUTINE mpi_set_error_handler

    INTEGER :: errhandler

    CALL MPI_COMM_CREATE_ERRHANDLER(mpi_error_handler, errhandler, errcode)
    CALL MPI_COMM_SET_ERRHANDLER(MPI_COMM_WORLD, errhandler, errcode)

  END SUBROUTINE mpi_set_error_handler



  SUBROUTINE mpi_error_handler(comm, error_code)

    INTEGER :: comm, error_code
    REAL :: tmp1, tmp2
    CHARACTER(LEN=29) :: errstring(0:MPI_ERR_LASTCODE)

    errstring(MPI_SUCCESS                  ) = 'MPI_SUCCESS                  '
    errstring(MPI_ERR_BUFFER               ) = 'MPI_ERR_BUFFER               '
    errstring(MPI_ERR_COUNT                ) = 'MPI_ERR_COUNT                '
    errstring(MPI_ERR_TYPE                 ) = 'MPI_ERR_TYPE                 '
    errstring(MPI_ERR_TAG                  ) = 'MPI_ERR_TAG                  '
    errstring(MPI_ERR_COMM                 ) = 'MPI_ERR_COMM                 '
    errstring(MPI_ERR_RANK                 ) = 'MPI_ERR_RANK                 '
    errstring(MPI_ERR_REQUEST              ) = 'MPI_ERR_REQUEST              '
    errstring(MPI_ERR_ROOT                 ) = 'MPI_ERR_ROOT                 '
    errstring(MPI_ERR_GROUP                ) = 'MPI_ERR_GROUP                '
    errstring(MPI_ERR_OP                   ) = 'MPI_ERR_OP                   '
    errstring(MPI_ERR_TOPOLOGY             ) = 'MPI_ERR_TOPOLOGY             '
    errstring(MPI_ERR_DIMS                 ) = 'MPI_ERR_DIMS                 '
    errstring(MPI_ERR_ARG                  ) = 'MPI_ERR_ARG                  '
    errstring(MPI_ERR_UNKNOWN              ) = 'MPI_ERR_UNKNOWN              '
    errstring(MPI_ERR_TRUNCATE             ) = 'MPI_ERR_TRUNCATE             '
    errstring(MPI_ERR_OTHER                ) = 'MPI_ERR_OTHER                '
    errstring(MPI_ERR_INTERN               ) = 'MPI_ERR_INTERN               '
    errstring(MPI_ERR_IN_STATUS            ) = 'MPI_ERR_IN_STATUS            '
    errstring(MPI_ERR_PENDING              ) = 'MPI_ERR_PENDING              '
    errstring(MPI_ERR_ACCESS               ) = 'MPI_ERR_ACCESS               '
    errstring(MPI_ERR_AMODE                ) = 'MPI_ERR_AMODE                '
    errstring(MPI_ERR_ASSERT               ) = 'MPI_ERR_ASSERT               '
    errstring(MPI_ERR_BAD_FILE             ) = 'MPI_ERR_BAD_FILE             '
    errstring(MPI_ERR_BASE                 ) = 'MPI_ERR_BASE                 '
    errstring(MPI_ERR_CONVERSION           ) = 'MPI_ERR_CONVERSION           '
    errstring(MPI_ERR_DISP                 ) = 'MPI_ERR_DISP                 '
    errstring(MPI_ERR_DUP_DATAREP          ) = 'MPI_ERR_DUP_DATAREP          '
    errstring(MPI_ERR_FILE_EXISTS          ) = 'MPI_ERR_FILE_EXISTS          '
    errstring(MPI_ERR_FILE_IN_USE          ) = 'MPI_ERR_FILE_IN_USE          '
    errstring(MPI_ERR_FILE                 ) = 'MPI_ERR_FILE                 '
    errstring(MPI_ERR_INFO_KEY             ) = 'MPI_ERR_INFO_KEY             '
    errstring(MPI_ERR_INFO_NOKEY           ) = 'MPI_ERR_INFO_NOKEY           '
    errstring(MPI_ERR_INFO_VALUE           ) = 'MPI_ERR_INFO_VALUE           '
    errstring(MPI_ERR_INFO                 ) = 'MPI_ERR_INFO                 '
    errstring(MPI_ERR_IO                   ) = 'MPI_ERR_IO                   '
    errstring(MPI_ERR_KEYVAL               ) = 'MPI_ERR_KEYVAL               '
    errstring(MPI_ERR_LOCKTYPE             ) = 'MPI_ERR_LOCKTYPE             '
    errstring(MPI_ERR_NAME                 ) = 'MPI_ERR_NAME                 '
    errstring(MPI_ERR_NO_MEM               ) = 'MPI_ERR_NO_MEM               '
    errstring(MPI_ERR_NOT_SAME             ) = 'MPI_ERR_NOT_SAME             '
    errstring(MPI_ERR_NO_SPACE             ) = 'MPI_ERR_NO_SPACE             '
    errstring(MPI_ERR_NO_SUCH_FILE         ) = 'MPI_ERR_NO_SUCH_FILE         '
    errstring(MPI_ERR_PORT                 ) = 'MPI_ERR_PORT                 '
    errstring(MPI_ERR_QUOTA                ) = 'MPI_ERR_QUOTA                '
    errstring(MPI_ERR_READ_ONLY            ) = 'MPI_ERR_READ_ONLY            '
    errstring(MPI_ERR_RMA_CONFLICT         ) = 'MPI_ERR_RMA_CONFLICT         '
    errstring(MPI_ERR_RMA_SYNC             ) = 'MPI_ERR_RMA_SYNC             '
    errstring(MPI_ERR_SERVICE              ) = 'MPI_ERR_SERVICE              '
    errstring(MPI_ERR_SIZE                 ) = 'MPI_ERR_SIZE                 '
    errstring(MPI_ERR_SPAWN                ) = 'MPI_ERR_SPAWN                '
    errstring(MPI_ERR_UNSUPPORTED_DATAREP  ) = 'MPI_ERR_UNSUPPORTED_DATAREP  '
    errstring(MPI_ERR_UNSUPPORTED_OPERATION) = 'MPI_ERR_UNSUPPORTED_OPERATION'
    errstring(MPI_ERR_WIN                  ) = 'MPI_ERR_WIN                  '
    errstring(MPI_ERR_LASTCODE             ) = 'MPI_ERR_LASTCODE             '

    PRINT*, 'Caught MPI error: ', TRIM(errstring(error_code))
    IF (comm == MPI_COMM_WORLD) THEN
      PRINT*, 'Communicator MPI_COMM_WORLD'
    ELSE
      PRINT*, 'Communicator ', comm, '(Not MPI_COMM_WORLD)'
    END IF

    ! Deliberately raise a divide-by-zero error
    tmp1 = 0.0
    tmp2 = 1.0 / tmp1

  END SUBROUTINE mpi_error_handler
#endif

END MODULE mpi_routines
