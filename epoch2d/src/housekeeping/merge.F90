MODULE merge
#ifdef MERGE_PARTICLES

USE shared_data
USE partlist
USE utilities
USE constants, only: c

IMPLICIT NONE

CONTAINS
  SUBROUTINE merge_particles

    INTEGER :: ispecies, n, i, ix, iy, nmerges, max_ncount, nc, nkeep, &
         niter, ngroup, npre, nafter, ndiv, nrem, ncell_max, nmerge_start

    INTEGER :: npre_r, nafter_r, ierr, im
    REAL(num) :: avg_energy, avg_v(3), cur_v(3), cur_energy, &
         cell_v(3), cell_energy, cell_weight, &
         wid_energy, wid_v(3), cur_weight, &
         energy_fac, pcomp_fac, mass, weight, &
         cur_pos(3), avg_pos(3)
    LOGICAL :: keep
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

          keep = ANY(ABS(cur_v - cell_v) > wid_v)
          keep = keep .OR. (ABS(cell_energy - cur_energy) > wid_energy)

          IF(keep) THEN
            CALL remove_particle_from_partlist(plist, cur)
            CALL add_particle_to_partlist(keep_plist, cur)
          END IF
          cur => next
        END DO

        ncell_max = max_ncount - keep_plist%count

        ! cycle if none remain after filtering

        IF (ncell_max <= 0) THEN
          nafter = nafter + plist%count
          CYCLE ixlp
        END IF

!       PRINT '("cell_weight = ", ES14.8)', cell_weight/(dx*dy)
!
!       PRINT '("linear density = ", ES14.8)', &
!            media_density(ix,iy,1)
!       PRINT '("sum = ", ES14.8)', &
!            cell_weight/(dx*dy) + media_density(ix,iy,1)

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
        CALL append_partlist(plist, keep_plist)
        CALL destroy_partlist(keep_plist)

        nafter = nafter + plist%count
      END DO ixlp
      END DO

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

  END SUBROUTINE merge_particles
#endif
END MODULE merge
