MODULE merge
#ifdef MERGE_PARTICLES

USE shared_data
USE partlist
USE utilities
USE constants, only: c

IMPLICIT NONE
PRIVATE
PUBLIC :: merge_particles

CONTAINS

  ! delete later
  SUBROUTINE dump_parts(ix, iy, ispecies)
    INTEGER, INTENT(IN) :: ix, iy, ispecies
    INTEGER :: ifile, i, n
    CHARACTER(len=128) :: name
    TYPE(particle_list), POINTER :: plist
    TYPE(particle), POINTER :: cur, next
    REAL(num) :: v

    ifile = 101

    plist => species_list(ispecies)%secondary_list(ix, iy)

    n = plist%count

    write(name, '("dump-", I2.2,".",I2.2,".dat")') ix, iy
    open(ifile, file=name, access='stream')

    write(ifile) n

    cur => plist%head
    DO i=1, n
      write(ifile) cur%weight
      cur => cur%next
    END DO

    cur => plist%head
    DO i=1, n
      write(ifile) cur%part_pos(1)
      cur => cur%next
    END DO

    cur => plist%head
    DO i=1, n
      write(ifile) cur%part_pos(2)
      cur => cur%next
    END DO

    cur => plist%head
    DO i=1, n
      write(ifile) cur%part_p(1)
      cur => cur%next
    END DO

    cur => plist%head
    DO i=1, n
      write(ifile) cur%part_p(2)
      cur => cur%next
    END DO

    cur => plist%head
    DO i=1, n
      write(ifile) cur%part_p(3)
      cur => cur%next
    END DO

    close(ifile)
  END SUBROUTINE dump_parts
  ! end delete

  SUBROUTINE merge_particles

    SELECT CASE (merge_scheme)
    CASE (1)
      IF (rank == 0) PRINT '("warning: merge scheme 1 is broken I think")'
      CALL merge_scheme1
    CASE (2)
      IF (rank == 0) PRINT '("warning: merge scheme 2 is broken I think")'
      CALL merge_scheme2
    CASE (3)
      CALL merge_scheme3
    CASE DEFAULT
      IF (rank == 0) PRINT '("ERROR, unknown merge scheme", I3.3)', merge_scheme
      CALL abort_code(c_err_bad_value)
    END SELECT

  END SUBROUTINE merge_particles

  SUBROUTINE merge_scheme1

    INTEGER :: ispecies, n, i, ix, iy, nmerges, max_ncount, nc, nkeep, &
         niter, ngroup, npre, nafter, ndiv, nrem, ncell_max, nmerge_start

    INTEGER :: npre_r, nafter_r, ierr, im
    REAL(num) :: avg_energy, avg_v(3), cur_v(3), cur_energy, &
         cell_v(3), cell_energy, cell_weight, &
         wid_energy, wid_v(3), cur_weight, &
         energy_fac, pcomp_fac, mass, weight, &
         cur_pos(3), avg_pos(3)
    LOGICAL :: keep, endprog
    TYPE(particle), POINTER ::  cur, next, repl
    TYPE(particle_list) :: temp_plist, keep_plist
    TYPE(particle_list), POINTER :: plist

isp:DO ispecies=1,n_species

      IF (.NOT. species_list(ispecies)%merge) CYCLE isp

      energy_fac = species_list(ispecies)%merge_max_energy_sig
      pcomp_fac = species_list(ispecies)%merge_max_pcomp_sig
      max_ncount = species_list(ispecies)%merge_max_particles
      mass = species_list(ispecies)%mass
      nmerge_start = species_list(ispecies)%merge_start
      npre = 0 ; nafter = 0

      DO iy=1,ny
ixlp: DO ix=1,nx

        plist => species_list(ispecies)%secondary_list(ix,iy)
        npre = npre + plist%count

        ! early filter to skip cells with too few particles
        IF (plist%count <= nmerge_start) THEN 
          nafter = nafter + plist%count
          CYCLE ixlp
        END IF

        ! calculate cell average
        cell_v = 0 ; cell_energy = 0
        cur => plist%head
        DO i = 1, plist%count
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          cell_energy = cell_energy + cur_energy*cur%weight
          cell_v = cell_v + cur_v*cur%weight

          cell_weight = cell_weight + cur%weight
          cur => cur%next
        END DO
        cell_v = cell_v / cell_weight
        cell_energy = cell_energy / cell_weight

        ! get first stddev (really average rms)
        wid_energy = 0 ; wid_v = 0
        cur => plist%head
        DO i=1,plist%count
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          wid_energy = wid_energy &
               + (cell_energy-cur_energy)**2*cur%weight

          wid_v = wid_v + (cell_v - cur_v)**2*cur%weight
          cur => cur%next
        END DO
        wid_energy = sqrt(wid_energy / cell_weight) * energy_fac
        wid_v = sqrt(wid_v / cell_weight) * pcomp_fac

        ! keep particles out of prescribed stddevs

        CALL create_empty_partlist(keep_plist)
        cur => plist%head
        DO i = 1, plist%count
          next => cur%next
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          !keep = ANY(ABS(cur_v - cell_v) > wid_v)
          !keep = keep .OR. (ABS(cell_energy - cur_energy) > wid_energy)
          keep = .FALSE.!(ABS(cell_energy - cur_energy) > wid_energy)

          IF(keep) THEN
            CALL remove_particle_from_partlist(plist, cur)
            CALL add_particle_to_partlist(keep_plist, cur)
          END IF
          cur => next
        END DO

        ncell_max = max_ncount - keep_plist%count

        ! cycle if none remain after filtering

        IF (ncell_max <= 0) goto 01100 !end of ixlp loop

        ndiv = plist%count / ncell_max
        nrem = MOD(INT(plist%count), ncell_max)

niters: DO niter=1,ncell_max

          ngroup = ndiv
          IF (niter .LE. nrem) ngroup = ngroup + 1
          IF (ngroup == 1) EXIT niters

          CALL create_empty_partlist(temp_plist)
          ! take up to ngroup from plist
ncloop:   DO i = 1,ngroup
            cur => plist%head
            IF (.NOT. ASSOCIATED(cur)) EXIT ncloop
            CALL remove_particle_from_partlist(plist, cur)
            CALL add_particle_to_partlist(temp_plist, cur)
          END DO ncloop
          IF (temp_plist%count == 0) CYCLE niters

          ! get group averages
          avg_v = 0 ; avg_pos = 0 ; weight = 0
          cur => temp_plist%head
          DO i=1,temp_plist%count
            cur_v = cur%part_p/(mass*c)
            cur_weight = cur%weight
            cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

            avg_v = avg_v + cur_v*cur_weight
            avg_pos = avg_pos + cur%part_pos*cur_weight
            weight = weight + cur%weight
            cur => cur%next
          END DO
          avg_pos = avg_pos / weight
          avg_v = avg_v / weight

          NULLIFY(repl)
          CALL create_particle(repl)
          repl%part_pos = avg_pos
          repl%part_ip = avg_pos
          repl%part_p = avg_v*mass
          repl%weight = weight
          ! add the lumped replacement particle
          CALL add_particle_to_partlist(keep_plist, repl)
          CALL destroy_partlist(temp_plist)

        END DO niters
        ! add keeps to plist
01100   CALL append_partlist(plist, keep_plist)
        CALL destroy_partlist(keep_plist)

        nafter = nafter + plist%count
      END DO ixlp
      END DO

      !IF (endprog) STOP

      CALL MPI_REDUCE(nafter, nafter_r, 1, MPI_INTEGER, MPI_SUM, 0, &
           comm, ierr)
      CALL MPI_REDUCE(npre, npre_r, 1, MPI_INTEGER, MPI_SUM, 0, &
           comm, ierr)

      IF (rank == 0) THEN
        PRINT '("particle_merge at step=",I7.7,": species_list(",I2, &
            & "), count change ", I9,"->",I9)', &
             step, ispecies, npre_r, nafter_r
      END IF
    END DO isp

  END SUBROUTINE merge_scheme1


  SUBROUTINE merge_scheme2

    INTEGER :: ispecies, n, i, ix, iy, nmerges, max_ncount, nc, nkeep, &
         niter, ngroup, npre, nafter, ndiv, nrem, ncell_max, nmerge_start

    INTEGER :: npre_r, nafter_r, ierr, im
    REAL(num) :: avg_energy, avg_v(3), cur_v(3), cur_energy, &
         cell_energy, energy_cut, cell_weight, wid_energy, cur_weight, &
         energy_fac, pcomp_fac, mass, weight, &
         cur_pos(3), avg_pos(3)
    LOGICAL :: keep !, endprog
    TYPE(particle), POINTER ::  cur, next, repl
    TYPE(particle_list) :: temp_plist, keep_plist
    TYPE(particle_list), POINTER :: plist

    INTEGER :: ncells, ncells_r
    REAL(num) :: avg_ppc_bef, avg_ppc_aft

    !endprog = .FALSE.

isp:DO ispecies=1,n_species

      IF (.NOT. species_list(ispecies)%merge) CYCLE isp

      energy_fac = species_list(ispecies)%merge_max_energy_sig
      pcomp_fac = species_list(ispecies)%merge_max_pcomp_sig
      max_ncount = species_list(ispecies)%merge_max_particles
      mass = species_list(ispecies)%mass
      nmerge_start = species_list(ispecies)%merge_start
      energy_cut = species_list(ispecies)%merge_energy_cut

      npre = 0 ; nafter = 0

      ncells = 0 ; avg_ppc_bef = 0 ; avg_ppc_aft = 0

      DO iy=1,ny
ixlp: DO ix=1,nx

        plist => species_list(ispecies)%secondary_list(ix,iy)
        npre = npre + plist%count

        IF (plist%count .GT. 0) ncells = ncells + 1

        ! early filter to skip cells with too few particles
        IF (plist%count <= nmerge_start) THEN
          nafter = nafter + plist%count
          CYCLE ixlp
        END IF

        ! calculate cell average
        cell_energy = 0
        cur => plist%head
        DO i = 1, plist%count
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          cell_energy = cell_energy + cur_energy*cur%weight

          cell_weight = cell_weight + cur%weight
          cur => cur%next
        END DO
        cell_energy = cell_energy / cell_weight

        ! get first stddev (really average rms)

        wid_energy = 0
        cur => plist%head
        DO i = 1, plist%count
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          wid_energy = wid_energy &
               + (cell_energy-cur_energy)**2*cur%weight

          cur => cur%next
        END DO
        wid_energy = sqrt(wid_energy / cell_weight) * energy_fac

        ! keep particles out of prescribed stddevs

        CALL create_empty_partlist(keep_plist)
        cur => plist%head
        DO i = 1, plist%count
          next => cur%next
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          keep = cur_energy .GT. cell_energy + wid_energy
          keep = keep .OR. cur_energy .GT. energy_cut

          IF(keep) THEN
            CALL remove_particle_from_partlist(plist, cur)
            CALL add_particle_to_partlist(keep_plist, cur)
          END IF
          cur => next
        END DO

        IF (plist%count .LE. 1) goto 01200

        ncell_max = max_ncount - keep_plist%count
        ! cycle if none remain after filtering
        IF (ncell_max .LE. 0) ncell_max = 1

        ndiv = plist%count / ncell_max
        nrem = MOD(INT(plist%count), ncell_max)

niters: DO niter=1,ncell_max

          ngroup = ndiv
          IF (niter .LE. nrem) ngroup = ngroup + 1
          IF (ngroup == 1) EXIT niters

          CALL create_empty_partlist(temp_plist)
          ! take up to ngroup from plist
ncloop:   DO i = 1,ngroup
            cur => plist%head
            IF (.NOT. ASSOCIATED(cur)) EXIT ncloop
            CALL remove_particle_from_partlist(plist, cur)
            CALL add_particle_to_partlist(temp_plist, cur)
          END DO ncloop
          IF (temp_plist%count == 0) CYCLE niters

          ! get group averages
          avg_v = 0 ; avg_pos = 0 ; weight = 0
          cur => temp_plist%head
          DO i=1,temp_plist%count
            cur_v = cur%part_p/(mass*c)
            cur_weight = cur%weight
            cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

            avg_v = avg_v + cur_v*cur_weight
            avg_pos = avg_pos + cur%part_pos*cur_weight
            weight = weight + cur%weight
            cur => cur%next
          END DO
          avg_pos = avg_pos / weight
          avg_v = avg_v / weight

          NULLIFY(repl)
          CALL create_particle(repl)
          repl%part_pos = avg_pos
          repl%part_ip = avg_pos
          repl%part_p = avg_v*mass
          repl%weight = weight
          ! add the lumped replacement particle
          CALL add_particle_to_partlist(keep_plist, repl)
          CALL destroy_partlist(temp_plist)

        END DO niters
        ! add keeps to plist
01200   CALL append_partlist(plist, keep_plist)
        CALL destroy_partlist(keep_plist)

        nafter = nafter + plist%count
      END DO ixlp
      END DO



      CALL MPI_REDUCE(nafter, nafter_r, 1, MPI_INTEGER, MPI_SUM, 0, &
           comm, ierr)
      CALL MPI_REDUCE(npre, npre_r, 1, MPI_INTEGER, MPI_SUM, 0, &
           comm, ierr)
      CALL MPI_REDUCE(ncells, ncells_r, 1, MPI_INTEGER, MPI_SUM, 0, &
           comm, ierr)


      IF (rank == 0) THEN
        avg_ppc_bef = REAL(npre,num) / REAL(ncells_r, num)
        avg_ppc_aft = REAL(nafter,num) / REAL(ncells_r, num)

        PRINT '("particle_merge: step=",I6.6,": species_list(",I2, &
            & "), count change ", I9,"->",I9)', &
            step, ispecies, npre_r, nafter_r
        PRINT '("          avgs: ", F6.3,"->",F6.3)', &
            avg_ppc_bef, avg_ppc_aft
      END IF

    END DO isp

  END SUBROUTINE merge_scheme2

  SUBROUTINE merge_scheme3

    INTEGER :: ispecies, n, i, ix, iy, nmerges, max_ncount, nc, nkeep, &
      niter, ngroup, npre, nafter, ndiv, nrem, ncell_max, nmerge_start, &
      nexcess

    INTEGER :: npre_r, nafter_r, ierr, im
    REAL(num) :: avg_energy, avg_v(3), cur_v(3), cur_energy, &
      cell_energy, energy_cut, cell_weight, wid_energy, cur_weight, &
      energy_fac, pcomp_fac, mass, weight, stddev, post_cell_weight, &
      cur_pos(3), avg_pos(3), max_energy
    LOGICAL :: keep, remove
    TYPE(particle), POINTER ::  cur, next, repl
    TYPE(particle_list) :: temp_plist
    TYPE(particle_list), POINTER :: plist

    INTEGER :: ncells, ncells_r, ndiffcells, ndiffcells_r, ndiffs
    INTEGER :: ncurdiff
    INTEGER :: merge_end_nstep, merge_start_nstep
    REAL(num) :: avg_ppc_bef, avg_ppc_aft, &
      sum_mean, sum_mean_r, avg_mean, &
      sum_mean_max_spread, sum_mean_max_spread_r, avg_mean_max_spread, &
      sum_stddev, sum_stddev_r, avg_stddev, &
      sum_spread2stddev, sum_spread2stddev_r, avg_spread2stddev
    REAL(num) :: rest_energy

isp:DO ispecies=1,n_species

      IF (.NOT. species_list(ispecies)%merge) CYCLE isp

      rest_energy = species_list(ispecies)%mass*m0 * c**2
      energy_fac = species_list(ispecies)%merge_max_energy_sig
      pcomp_fac = species_list(ispecies)%merge_max_pcomp_sig
      max_ncount = species_list(ispecies)%merge_max_particles
      mass = species_list(ispecies)%mass
      nmerge_start = species_list(ispecies)%merge_start
      energy_cut = species_list(ispecies)%merge_energy_cut / rest_energy
      merge_end_nstep = species_list(ispecies)%merge_end_nstep
      merge_start_nstep = species_list(ispecies)%merge_start_nstep

      IF (step .GT. merge_end_nstep .OR. step .LT. merge_start_nstep) THEN
        IF (rank == 0) PRINT '("skipping merge for ", A12)', &
          species_list(ispecies)%name
        CYCLE isp
      END IF

      npre = 0 ; nafter = 0

      ncells = 0 ; avg_ppc_bef = 0 ; avg_ppc_aft = 0 ; ndiffs = 0

      sum_mean = 0 ; sum_mean_max_spread = 0
      sum_stddev = 0 ; sum_spread2stddev = 0

      DO iy=1,ny
ixlp: DO ix=1,nx

        plist => species_list(ispecies)%secondary_list(ix,iy)
        ncurdiff = plist%count
        npre = npre + plist%count

        IF (plist%count .GT. 0) ncells = ncells + 1

        ! early filter to skip cells with too few particles
        IF (plist%count < nmerge_start) THEN
          nafter = nafter + plist%count
          CYCLE ixlp
        END IF

        ndiffcells = ndiffcells + 1
        nexcess = plist%count + 1 - nmerge_start

        ! calculate cell average and stddev
        cell_weight = 0 ; stddev = 0
        cell_energy = 0 ; wid_energy = 0 ; max_energy = 0
        cur => plist%head
        DO i = 1, plist%count
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

          IF (cur_energy .GT. max_energy) max_energy = cur_energy

          cell_energy = cell_energy + cur_energy*cur%weight

          stddev = stddev + cur_energy**2*cur%weight
          cell_weight = cell_weight + cur%weight
          cur => cur%next
        END DO
        cell_energy = cell_energy / cell_weight
        stddev = sqrt(stddev / cell_weight - cell_energy**2)
        wid_energy = stddev * energy_fac

        sum_mean = sum_mean + cell_energy
        sum_mean_max_spread = sum_mean_max_spread + &
          abs(max_energy - cell_energy)

        sum_stddev = sum_stddev + stddev
        sum_spread2stddev = sum_spread2stddev + &
          abs(max_energy - cell_energy) / stddev

        ! get particles below threshold
        CALL create_empty_partlist(temp_plist)
        cur => plist%head
        DO i = 1, plist%count
          next => cur%next
          cur_v = cur%part_p/(mass*c)
          cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1


          remove = cur_energy .LT. cell_energy + wid_energy
          remove = remove .AND. cur_energy .LT. energy_cut
          IF (remove) THEN
            CALL remove_particle_from_partlist(plist, cur)
            CALL add_particle_to_partlist(temp_plist, cur)
          END IF
          IF (temp_plist%count .GE. nexcess) EXIT
          cur => next
        END DO

        IF (temp_plist%count .GT. 0) THEN
          avg_v = 0 ; avg_pos = 0 ; weight = 0
          cur => temp_plist%head
          DO i=1,temp_plist%count
            cur_v = cur%part_p/(mass*c)
            cur_weight = cur%weight
            cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

            avg_v = avg_v + cur_v*cur_weight
            cur_pos = cur%part_pos*cur_weight
            !PRINT '("rank=",I3,": cur_pos -> [",3ES11.2,"] vs. [",3ES11.2,"]*",ES11.2)', &
            !  rank, cur_pos, cur%part_pos, cur_weight
            avg_pos = avg_pos + cur_pos
            weight = weight + cur_weight
            cur => cur%next
          END DO
          avg_pos = avg_pos / weight
          avg_v = avg_v / weight

          NULLIFY(repl)
          CALL create_particle(repl)
          repl%part_pos = avg_pos
          repl%part_ip = avg_pos
          repl%part_p = avg_v*mass
          repl%weight = weight
          ! add the lumped replacement particle
          CALL add_particle_to_partlist(plist, repl)
        END IF
        CALL destroy_partlist(temp_plist)

        nafter = nafter + plist%count
#ifdef DEBUG_MERGE_LV2
        post_cell_weight = 0
        cur => plist%head
        DO i = 1, plist%count
          post_cell_weight = post_cell_weight + cur%weight
          cur => cur%next
        END DO

        PRINT '("rank=",I3," (ix,iy)=(",I4.4,",",I4.4,"): ",&
          &"cell_weight ",ES8.1," -> ",ES8.1)', rank, ix, iy, &
          cell_weight, post_cell_weight
#endif
      END DO ixlp
      END DO

#ifdef DEBUG_MERGE
      CALL MPI_REDUCE(nafter, nafter_r, 1, MPI_INTEGER, MPI_SUM, 0, &
        comm, ierr)
      CALL MPI_REDUCE(npre, npre_r, 1, MPI_INTEGER, MPI_SUM, 0, &
        comm, ierr)
      CALL MPI_REDUCE(ncells, ncells_r, 1, MPI_INTEGER, MPI_SUM, 0, &
        comm, ierr)
      CALL MPI_REDUCE(sum_mean_max_spread, sum_mean_max_spread_r, 1, &
        mpireal, MPI_SUM, 0, comm, ierr)
      CALL MPI_REDUCE(sum_mean, sum_mean_r, 1, mpireal, MPI_SUM, 0, &
        comm, ierr)
      CALL MPI_REDUCE(sum_stddev, sum_stddev_r, 1, mpireal, &
        MPI_SUM, 0, comm, ierr)
      CALL MPI_REDUCE(sum_spread2stddev, sum_spread2stddev_r, 1, mpireal, &
        MPI_SUM, 0, comm, ierr)

      IF (rank == 0) THEN

        avg_ppc_bef = REAL(npre_r,num) / REAL(ncells_r, num)
        avg_ppc_aft = REAL(nafter_r,num) / REAL(ncells_r, num)

        avg_mean = sum_mean_r / REAL(ncells_r, num)
        avg_mean_max_spread = sum_mean_max_spread_r / REAL(ncells_r, num)
        avg_stddev = sum_stddev_r / REAL(ncells_r, num)
        avg_spread2stddev = sum_spread2stddev_r / REAL(ncells_r, num)

        IF (ncells_r > 0 ) THEN
          PRINT '(" particle_merge3: step=",I7.7,": species_list(",I2, &
            & "), count change ", I9,"->",I9)', &
            step, ispecies, npre_r, nafter_r
          PRINT '("            avgs: ", F6.3,"->",F6.3)', &
            avg_ppc_bef, avg_ppc_aft
          PRINT '("-------")'
          PRINT '("         avg(mean): ", ES8.2)', avg_mean
          PRINT '("     avg(mean2max): ", ES8.2)', avg_mean_max_spread
          PRINT '("       avg(stddev): ", ES8.2)', avg_stddev
          PRINT '("avg(spread2stddev): ", ES8.2)', avg_spread2stddev
        ELSE
          PRINT '(" particle_merge3: step=",I6.6,": no non-zero cells")', &
            step
        END IF
      END IF
#endif
    END DO isp

  END SUBROUTINE merge_scheme3


#endif
END MODULE merge
