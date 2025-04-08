MODULE merge_particles
CONTAINS
  SUBROUTINE merge_particles

    INTEGER :: ispecies, n, i, ix, iy, nmerges, max_ncount, nc, nremove, &
         ngroup, ngroups
    REAL(num) :: avg_energy, avg_v(3), cur_v(3), cur_energy, &
         wid_energy, wid_v(3),energy_fac, pcomp_fac, mass, weight, &
         cur_pos(3), avg_pos(3)
    LOGICAL :: remove
    TYPE(particle), POINTER :: last, last_ig, first_ig, cur, last, next
    TYPE(partlcle_list), POINTER :: cur_plist

isp:DO ispecies=1,n_species

      IF (.NOT. species_list(ispecies)%merge) CYCLE isp

      energy_fac = species_list(ispecies)%merge_max_energy_sig
      pcomp_fac = species_list(ispecies)%merge_max_pcomp_sig
      max_ncount = species_list(ispecies)%merge_max_particles
      mass = species_list(ispecies)%mass

      DO iy=1,ny
ixlp: DO ix=1,nx

        cur_plist => species_list(ispecies)%secondary_list(ix,iy)
        IF (cur_plist%count <= max_ncount) CYCLE ixlp

        last => cur_plist%tail
        first=> cur_plist%head
        first_in_g => first

        ngroups = cur_plis%count / max_ncount

ngrps:  DO ngroup=1,ngroups
          avg_v = 0 ; avg_pos = 0 ; cur_energy = 0
          weight = 0
          nc = 0
          DO WHILE(ASSOCIATED(cur) .AND. nc < max_ncount)
            nc = nc + 1
            IF (.NOT. ASSOCIATED(last_in_g)) THEN
              PRINT "(A)", "unexpected unassociated pointer in merge"
              abort_code(c_err_generic_error)
            END IF
          END DO

          ! get averages
          cur => first_in_g
          DO i=1,nc
            cur_v = cur%part_p/mass
            cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

            avg_energy = avg_energy + cur_energy*cur%weight
            avg_v = avg_v + cur_v*cur%weight

            weight = weight + cur%weight
            cur => cur%next
          END DO
          avg_energy = avg_energy / weight
          avg_v = avg_v / weight

          ! get stddev (really average rms)
          cur => first_in_g
          DO i=1,nc
            cur_v = cur%part_p/mass
            cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1

            wid_energy = wid_energy + (avg_energy-cur_energy)**2*cur%weight
            wid_v = wid_v + (avg_v - cur_v)**2*cur%weight
            cur => cur%next
          END DO
          wid_energy = sqrt(wid_energy / weight)
          wid_v = sqrt(wid_v / weight)
          wid_energy = wid_energy * energy_fac
          wid_v = wid_v * energy_fac

          !now remove if they lie outside of the prescribed stddevs
          nremove = 0
          cur => first_in_g
          DO n=1,nc
            next => cur%next
            cur_v = cur%part_p/mass
            cur_energy = sqrt(cur_v(1)**2+cur_v(2)**2+cur_v(3)**2+1)-1
            remove = .FALSE.
            IF ( ABS(cur_energy-avg_energy) > wid_energy ) remove=.TRUE.
            IF ( ANY( ABS(cur_v-avg_v) > wid_v ) ) remove = .TRUE.
            IF (remove) THEN
              IF (ASSOCIATED(cur%prev)) cur%prev%next => next
              IF (ASSOCIATED(cur, first))      first=>next
              IF (ASSOCIATED(cur, first_in_g)) first_in_g=>next

              cur%prev => last
              last%next=> cur
              last => cur
              NULLIFY(cur%next)
              nremove = nremove + 1
            END IF
            cur => next
          END DO
          nc = nc - nremove
          IF (nremove > 0) THEN
            avg_v = 0
            weight = 0
            ! one more time, calculate the average velocities
            cur => first_in_g
            DO n=1,nc
              cur_v = cur%part_p/mass
              avg_v = avg_v + cur_v*cur%weight
              weight = weight + cur%weight
              cur => cur%next
            END DO
            avg_v = avg_v / weight
          END IF
          ! get average position, summed weight
          cur => first_in_g
          DO n=1,nc
            avg_pos = avg_pos + cur%part_pos*cur%weight
            weight  = weight + cur%weight
          END DO
          avg_pos = avg_pos / weight
          NULLIFY(cur)
          CALL create_particle(cur)
          cur%part_pos = avg_pos
          cur%part_ip = avg_pos
          cur%part_p = avg_v*mass
          cur%weight = weight
          ! first, fix part list
          cur_plist%head => first
          cur_plist%tail => last
          cur => first_in_g
          ! remove parts
          DO n=1,nc
            next => cur%next
            CALL remove_particle_from_partlist(&
                 cur_plist, cur)
            cur => next
          END DO
          ! add the lumped replacement particle
          CALL add_particle_to_partlist(&
               cur_plist, cur)
        END DO ngrps
      END DO ixlp
      END DO
    END DO isp

  END SUBROUTINE merge_particles
END MODULE merge_particles
