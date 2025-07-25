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

MODULE particles

  USE boundary
  USE partlist
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE
#ifdef BOUND_HARMONIC
  
  PRIVATE
  
  PUBLIC :: push_particles, f0, check_bound_particles, clean_diminished_parts
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
  PUBLIC :: push_photons
#endif !photons

  REAL(num) :: eye(3,3) = reshape(   &
       (/ 1.0_num, 0.0_num, 0.0_num, &
          0.0_num, 1.0_num, 0.0_num, &
          0.0_num, 0.0_num, 1.0_num/), (/3,3/))
#endif !BOUND_HARMONIC
  
CONTAINS
#ifdef BOUND_HARMONIC
  PURE FUNCTION bsq(ibx,iby,ibz)
    REAL(num), INTENT(in) :: ibx, iby, ibz
    REAL(num) bsq(3,3)
    bsq(1,1) = ibx*ibx
    bsq(1,2) = ibx*iby
    bsq(1,3) = ibx*ibz
    bsq(2,2) = iby*iby
    bsq(2,3) = iby*ibz
    bsq(3,3) = ibz*ibz
    bsq(2,1) = bsq(1,2)
    bsq(3,1) = bsq(1,3)
    bsq(3,2) = bsq(2,3)
  END FUNCTION bsq
  PURE FUNCTION crossB(ibx,iby,ibz)
    REAL(num), INTENT(in) :: ibx, iby, ibz
    REAL(num) crossB(3,3)
    crossB = 0.0d0
    crossB(1,2) = ibz
    crossB(1,3) =-iby
    crossB(2,3) = ibx
    crossB(2,1) =-ibz
    crossB(3,1) = iby
    crossB(3,2) =-ibx
  END FUNCTION crossB
#endif

  SUBROUTINE push_particles

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3
    INTEGER :: cell_y1, cell_y2, cell_y3
    INTEGER :: cell_z1, cell_z2, cell_z3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num) :: jxh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jyh
    REAL(num), DIMENSION(sf0-1:sf1+1,sf0-1:sf1+1) :: jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight, part_m
#ifdef HC_PUSH
    REAL(num) :: beta_x, beta_y, beta_z, beta2, beta_dot_u, alpha, sigma
#endif

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, part_mc2
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: idx, idy, idz
    REAL(num) :: idtyz, idtxz, idtxy
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, dtfac, gamma_rel, part_u2, third
    REAL(num) :: delta_x, delta_y, delta_z
    REAL(num) :: xfac1, xfac2, yfac1, yfac2, zfac1, zfac2
    REAL(num) :: gz_iz, hz_iz, hygz, hyhz, hzyfac1, hzyfac2, yzfac
    INTEGER :: ispecies, ix, iy, iz, dcellx, dcelly, dcellz, cx, cy, cz
    INTEGER(i8) :: ipart
#ifdef WORK_DONE_INTEGRATED
#ifdef NONLIN
#error "WORK_DONE_INTEGRATED currently incompatible with NONLIN"
#endif
    REAL(num) :: tmp_x, tmp_y, tmp_z
    REAL(num) :: work_x, work_y, work_z
#endif
#ifndef NO_PARTICLE_PROBES
    LOGICAL :: probes_for_species
    REAL(num) :: gamma_rel_m1
#endif
#ifndef NO_TRACER_PARTICLES
    LOGICAL :: not_zero_current_species
#endif
    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif
#ifdef DELTAF_METHOD
    REAL(num) :: weight_back
#endif

    TYPE(particle), POINTER :: current, next
#ifdef BOUND_HARMONIC
    
    REAL(num), DIMENSION(3) :: lomegasq_dto2c, lgammadto2
    REAL(num), DIMENSION(3) :: one_p_lg, one_m_lg, upvec
    REAL(num) &
         Bdotprod,     &! B dot B*(1 + dt*gamma/2)
         detLm,        &! det(L_-)
         detLm_m_Bdot, &! determinant without Bdotprod
         dlx(3)         ! x, y, z of offset from binding centre
    REAL(num) :: bfield_factor
    REAL(num), DIMENSION(3,3) :: vm2vp, gminv
    !REAL(num) pre_part_pos(3),tmp(3)
#ifdef NONLIN
    REAL(num) :: fnl3, esq_part, elinfac
!end NONLIN
#endif
!end BOUND_HARMONIC
#endif

#ifdef PREFETCH
    CALL prefetch_particle(species_list(1)%attached_list%head)
#endif

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt * fac
    third = 1.0_num / 3.0_num

    idtyz = idt * idy * idz * fac
    idtxz = idt * idx * idz * fac
    idtxy = idt * idx * idy * fac

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%immobile) CYCLE
#ifdef SKIP_ZEROPARTS
      IF (species_list(ispecies)%count == 0) CYCLE
#endif
      IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
#ifdef BREMSSTRAHLUNG
        IF (ispecies == bremsstrahlung_photon_species) THEN
          IF (bremsstrahlung_photon_dynamics) THEN
            CALL push_photons(ispecies)
          ELSE
            CYCLE
          END IF
        END IF
#endif
#ifdef PHOTONS
        IF (photon_dynamics) CALL push_photons(ispecies)
#endif
        CYCLE
      END IF
#ifndef NO_PARTICLE_PROBES
      current_probe => species_list(ispecies)%attached_probes
      probes_for_species = ASSOCIATED(current_probe)
#endif
#ifndef NO_TRACER_PARTICLES
      not_zero_current_species = .NOT. species_list(ispecies)%zero_current
#endif

#ifdef PER_SPECIES_WEIGHT
      part_weight = species_list(ispecies)%weight
      fcx = idtyz * part_weight
      fcy = idtxz * part_weight
      fcz = idtxy * part_weight
#endif
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q   = species_list(ispecies)%charge
      part_m   = species_list(ispecies)%mass
      part_mc  = c * species_list(ispecies)%mass
      ipart_mc = 1.0_num / part_mc
      cmratio  = part_q * dtfac * ipart_mc
      ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
      part_mc2 = c * part_mc
#endif
#endif
#ifdef BOUND_HARMONIC
      lomegasq_dto2c= species_list(ispecies)%harmonic_omega**2 * dto2 / c
      lgammadto2    = species_list(ispecies)%harmonic_gamma * dto2
      one_p_lg = 1 + lgammadto2
      one_m_lg = 1 - lgammadto2
      detLm_m_Bdot = one_p_lg(1)*one_p_lg(2)*one_p_lg(3)
      bfield_factor = species_list(ispecies)%bfield_sample_factor

      gminv = 0.0_num
      ! gminv(1,1) = one_m_lg(1) / one_p_lg(1)
      ! gminv(2,2) = one_m_lg(2) / one_p_lg(2)
      ! gminv(3,3) = one_m_lg(3) / one_p_lg(3)
      ! vm2vp = gminv
      gminv(1,1) = one_m_lg(1)*one_p_lg(2)*one_p_lg(3)
      gminv(2,2) = one_p_lg(1)*one_m_lg(2)*one_p_lg(3)
      gminv(3,3) = one_p_lg(1)*one_p_lg(2)*one_m_lg(3)
#ifdef NONLIN
      fnl3 = species_list(ispecies)%nl_alpha3 * dtfac / c !dto2 or dtfac ???
      elinfac = species_list(ispecies)%linear_factor
!end NONLIN
#endif
!end BOUND_HARMONIC
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next
#ifdef PREFETCH
        CALL prefetch_particle(next)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
        fcx = idtyz * part_weight
        fcy = idtxz * part_weight
        fcz = idtxy * part_weight
#endif
#ifndef NO_PARTICLE_PROBES
        init_part_x = current%part_pos(1)
        init_part_y = current%part_pos(2)
        init_part_z = current%part_pos(3)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q   = current%charge
        part_m   = current%mass
        part_mc  = c * current%mass
        ipart_mc = 1.0_num / part_mc
        cmratio  = part_q * dtfac * ipart_mc
        ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
        part_mc2 = c * part_mc
#endif
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_grid_min_local
        part_y  = current%part_pos(2) - y_grid_min_local
        part_z  = current%part_pos(3) - z_grid_min_local
        part_ux = current%part_p(1) * ipart_mc
        part_uy = current%part_p(2) * ipart_mc
        part_uz = current%part_p(3) * ipart_mc

        ! Calculate v(t) from p(t)
        ! See PSC manual page (25-27)
        gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        root = dtco2 / gamma_rel

        ! Move particles to half timestep position to first order
        part_x = part_x + part_ux * root
        part_y = part_y + part_uy * root
        part_z = part_z + part_uz * root
#ifdef BOUND_HARMONIC
        ! should this also be at half step?
        dlx = current%part_pos - current%part_ip
        dlx = dlx + (/part_ux, part_uy, part_uz/)*root
#endif
        
#ifdef WORK_DONE_INTEGRATED
        ! This is the actual total work done by the fields: Results correspond
        ! with the electron's gamma factor
        root = cmratio / gamma_rel

        tmp_x = part_ux * root
        tmp_y = part_uy * root
        tmp_z = part_uz * root
#endif

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
        cell_z_r = part_z * idz - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
        cell_z_r = part_z * idz
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        cell_z1 = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z1, num) - cell_z_r
        cell_z1 = cell_z1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        cell_z2 = FLOOR(cell_z_r)
        cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
        cell_z2 = cell_z2 + 1

        dcellx = 0
        dcelly = 0
        dcellz = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

#ifdef GLOBALFIELD
        ex_part = ex_part + global_e(1)
        ey_part = ey_part + global_e(2)
        ez_part = ez_part + global_e(3)

        bx_part = bx_part + global_b(1)
        by_part = by_part + global_b(2)
        bz_part = bz_part + global_b(3)
#endif
        ! update particle momenta using weighted fields
#if defined(BOUND_HARMONIC) && defined(NONLIN)
        uxm = part_ux + cmratio * ex_part * elinfac
        uym = part_uy + cmratio * ey_part * elinfac
        uzm = part_uz + cmratio * ez_part * elinfac
#else
        uxm = part_ux + cmratio * ex_part
        uym = part_uy + cmratio * ey_part
        uzm = part_uz + cmratio * ez_part
#endif
#ifdef BOUND_HARMONIC
        ! add lorentzian omega terms
        uxm = uxm - lomegasq_dto2c(1)*dlx(1)
        uym = uym - lomegasq_dto2c(2)*dlx(2)
        uzm = uzm - lomegasq_dto2c(3)*dlx(3)
#ifdef NONLIN
        !third order nonlinearity for now
        esq_part = ex_part**2 + ey_part**2 + ez_part**2
        uxm = uxm + fnl3 * esq_part * ex_part
        uym = uym + fnl3 * esq_part * ey_part
        uzm = uzm + fnl3 * esq_part * ez_part
!end NONLIN
#endif
!end BOUND_HARMONIC
#endif
        
#ifdef HC_PUSH
        ! Half timestep, then use Higuera-Cary push see
        ! https://aip.scitation.org/doi/10.1063/1.4979989
        gamma_rel = uxm**2 + uym**2 + uzm**2 + 1.0_num
        alpha = 0.5_num * part_q * dt / part_m
        beta_x = alpha * bx_part
        beta_y = alpha * by_part
        beta_z = alpha * bz_part
        beta2 = beta_x**2 + beta_y**2 + beta_z**2
        sigma = gamma_rel - beta2
        beta_dot_u = beta_x * uxm + beta_y * uym + beta_z * uzm
        gamma_rel = sigma + SQRT(sigma**2 + 4.0_num * (beta2 + beta_dot_u**2))
        gamma_rel = SQRT(0.5_num * gamma_rel)
#else
        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
        gamma_rel = SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)
#endif
        root = ccmratio / gamma_rel
#ifdef BOUND_HARMONIC
        !print "(a,es11.2)", "multiplying b*_part with ", bfield_factor
        bx_part = bx_part*bfield_factor*root
        by_part = by_part*bfield_factor*root
        bz_part = bz_part*bfield_factor*root

        ! construct v- to v+
        vm2vp = 0
        ! B dot product
        Bdotprod = bx_part*bx_part*one_p_lg(1) &
                 + by_part*by_part*one_p_lg(2) &
                 + bz_part*bz_part*one_p_lg(3)
        ! diagonal
        vm2vp = gminv - Bdotprod*eye
        ! add cross
        vm2vp = vm2vp + 2.0_num*crossB(&
             bx_part*one_p_lg(1),&
             by_part*one_p_lg(2),&
             bz_part*one_p_lg(3))
        ! add B tensorprod B
        vm2vp = vm2vp + 2.0_num*Bsq(bx_part,by_part,bz_part)
        ! determinant
        detLm = detLm_m_Bdot + Bdotprod
        vm2vp = vm2vp / detLm
        ! obtain vp by matrix mul
        upvec = matmul(vm2vp, (/uxm,uym,uzm/))
        uxp = upvec(1)
        uyp = upvec(2)
        uzp = upvec(3)
#else

        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        taux2 = taux**2
        tauy2 = tauy**2
        tauz2 = tauz**2

        tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

        uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
            + 2.0_num * ((taux * tauy + tauz) * uym &
            + (taux * tauz - tauy) * uzm)) * tau
        uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
            + 2.0_num * ((tauy * tauz + taux) * uzm &
            + (tauy * taux - tauz) * uxm)) * tau
        uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
            + 2.0_num * ((tauz * taux + tauy) * uxm &
            + (tauz * tauy - taux) * uym)) * tau
#endif !BOUND_HARMONIC
        
        ! Rotation over, go to full timestep
#if defined(BOUND_HARMONIC) && defined(NONLIN)
        part_ux = uxp + cmratio * ex_part * elinfac
        part_uy = uyp + cmratio * ey_part * elinfac
        part_uz = uzp + cmratio * ez_part * elinfac
#else
        part_ux = uxp + cmratio * ex_part
        part_uy = uyp + cmratio * ey_part
        part_uz = uzp + cmratio * ez_part
#endif
#ifdef BOUND_HARMONIC
        ! add lorentzian omega terms
        part_ux = part_ux - lomegasq_dto2c(1)*dlx(1)
        part_uy = part_uy - lomegasq_dto2c(2)*dlx(2)
        part_uz = part_uz - lomegasq_dto2c(3)*dlx(3)
#ifdef NONLIN
        !third order nonlinearity for now
        esq_part = ex_part**2 + ey_part**2 + ez_part**2
        part_ux = part_ux + fnl3 * esq_part * ex_part
        part_uy = part_uy + fnl3 * esq_part * ey_part
        part_uz = part_uz + fnl3 * esq_part * ez_part
!end NONLIN
#endif
!end BOUND_HARMONIC
#endif

        ! Calculate particle velocity from particle momentum
        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel = SQRT(part_u2 + 1.0_num)
        root = dtco2 / gamma_rel

        delta_x = part_ux * root
        delta_y = part_uy * root
        delta_z = part_uz * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + delta_x
        part_y = part_y + delta_y
        part_z = part_z + delta_z

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = (/ part_x + x_grid_min_local, &
            part_y + y_grid_min_local, part_z + z_grid_min_local /)
        current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

#ifdef WORK_DONE_INTEGRATED
        ! This is the actual total work done by the fields: Results correspond
        ! with the electron's gamma factor
        root = cmratio / gamma_rel

        work_x = ex_part * (tmp_x + part_ux * root)
        work_y = ey_part * (tmp_y + part_uy * root)
        work_z = ez_part * (tmp_z + part_uz * root)

        current%work_x = work_x
        current%work_y = work_y
        current%work_z = work_z
        current%work_x_total = current%work_x_total + work_x
        current%work_y_total = current%work_y_total + work_y
        current%work_z_total = current%work_z_total + work_z
#endif

#ifndef NO_PARTICLE_PROBES
        final_part_x = current%part_pos(1)
        final_part_y = current%part_pos(2)
        final_part_z = current%part_pos(3)
#endif
        ! Original code calculates densities of electrons, ions and neutrals
        ! here. This has been removed to reduce memory footprint

        ! If the code is compiled with zero-current particle support then put in
        ! an IF statement so that the current is not calculated for this species
#ifndef NO_TRACER_PARTICLES
        IF (not_zero_current_species) THEN
#endif
          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + delta_x
          part_y = part_y + delta_y
          part_z = part_z + delta_z

          ! Delta-f calcuation: subtract background from
          ! calculated current.
#ifdef DELTAF_METHOD
          weight_back = current%pvol * f0(ispecies, part_mc / c, current%part_p)
          fcx = idtyz * (part_weight - weight_back)
          fcy = idtxz * (part_weight - weight_back)
          fcz = idtxy * (part_weight - weight_back)
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
          cell_x_r = part_x * idx - 0.5_num
          cell_y_r = part_y * idy - 0.5_num
          cell_z_r = part_z * idz - 0.5_num
#else
          cell_x_r = part_x * idx
          cell_y_r = part_y * idy
          cell_z_r = part_z * idz
#endif
          cell_x3 = FLOOR(cell_x_r + 0.5_num)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3 + 1

          cell_y3 = FLOOR(cell_y_r + 0.5_num)
          cell_frac_y = REAL(cell_y3, num) - cell_y_r
          cell_y3 = cell_y3 + 1

          cell_z3 = FLOOR(cell_z_r + 0.5_num)
          cell_frac_z = REAL(cell_z3, num) - cell_z_r
          cell_z3 = cell_z3 + 1

          hx = 0.0_num
          hy = 0.0_num
          hz = 0.0_num

          dcellx = cell_x3 - cell_x1
          dcelly = cell_y3 - cell_y1
          dcellz = cell_z3 - cell_z1
          ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          hx = hx - gx
          hy = hy - gy
          hz = hz - gz

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          xmin = sf_min + (dcellx - 1) / 2
          xmax = sf_max + (dcellx + 1) / 2

          ymin = sf_min + (dcelly - 1) / 2
          ymax = sf_max + (dcelly + 1) / 2

          zmin = sf_min + (dcellz - 1) / 2
          zmax = sf_max + (dcellz + 1) / 2

          fjx = fcx * part_q
          fjy = fcy * part_q
          fjz = fcz * part_q

          jzh = 0.0_num
          DO iz = zmin, zmax
            cz = cell_z1 + iz
            zfac1 =         gz(iz) + 0.5_num * hz(iz)
            zfac2 = third * hz(iz) + 0.5_num * gz(iz)

            gz_iz = gz(iz)
            hz_iz = hz(iz)

            jyh = 0.0_num
            DO iy = ymin, ymax
              cy = cell_y1 + iy
              yfac1 =         gy(iy) + 0.5_num * hy(iy)
              yfac2 = third * hy(iy) + 0.5_num * gy(iy)

              hygz = hy(iy) * gz_iz
              hyhz = hy(iy) * hz_iz
              yzfac = gy(iy) * zfac1 + hy(iy) * zfac2
              hzyfac1 = hz_iz * yfac1
              hzyfac2 = hz_iz * yfac2

              jxh = 0.0_num
              DO ix = xmin, xmax
                cx = cell_x1 + ix
                xfac1 =         gx(ix) + 0.5_num * hx(ix)
                xfac2 = third * hx(ix) + 0.5_num * gx(ix)

                wx = hx(ix) * yzfac
                wy = xfac1 * hygz + xfac2 * hyhz
                wz = gx(ix) * hzyfac1 + hx(ix) * hzyfac2

                ! This is the bit that actually solves d(rho)/dt = -div(J)
                jxh = jxh - fjx * wx
                jyh(ix) = jyh(ix) - fjy * wy
                jzh(ix, iy) = jzh(ix, iy) - fjz * wz

                jx(cx, cy, cz) = jx(cx, cy, cz) + jxh
                jy(cx, cy, cz) = jy(cx, cy, cz) + jyh(ix)
                jz(cx, cy, cz) = jz(cx, cy, cz) + jzh(ix, iy)
              END DO
            END DO
          END DO
#ifndef NO_TRACER_PARTICLES
        END IF
#endif
#if !defined(NO_PARTICLE_PROBES) && !defined(NO_IO)
        IF (probes_for_species) THEN
          ! Compare the current particle with the parameters of any probes in
          ! the system. These particles are copied into a separate part of the
          ! output file.

          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          current_probe => species_list(ispecies)%attached_probes

          ! Cycle through probes
          DO WHILE(ASSOCIATED(current_probe))
            ! Note that this is the energy of a single REAL particle in the
            ! pseudoparticle, NOT the energy of the pseudoparticle
            probe_energy = gamma_rel_m1 * part_mc2

            ! Unidirectional probe
            IF (probe_energy > current_probe%ek_min) THEN
              IF (probe_energy < current_probe%ek_max) THEN

                d_init  = SUM(current_probe%normal * (current_probe%point &
                    - (/init_part_x, init_part_y, init_part_z/)))
                d_final = SUM(current_probe%normal * (current_probe%point &
                    - (/final_part_x, final_part_y, final_part_z/)))
                IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                  ! this particle is wanted so copy it to the list associated
                  ! with this probe
                  ALLOCATE(particle_copy)
                  particle_copy = current
                  CALL add_particle_to_partlist(&
                      current_probe%sampled_particles, particle_copy)
                  NULLIFY(particle_copy)
                END IF

              END IF
            END IF
            current_probe => current_probe%next
          END DO
        END IF
#endif
        current => next
      END DO
      CALL current_bcs(species=ispecies)
    END DO

    CALL particle_bcs

  END SUBROUTINE push_particles

#ifdef BOUND_HARMONIC
  SUBROUTINE check_bound_particles
    INTEGER ispecies,ipart
    TYPE(particle), POINTER :: current
    REAL(num), DIMENSION(3) :: omega, gamma, dlx, part_u, maxd
    REAL(num) part_m
    
    IF (check_bound_interval == 0 ) RETURN
    IF (step == 0 .OR. MOD(step, check_bound_interval) > 0) RETURN
    IF (rank == 0) PRINT "(A,I6.6)", "checking bound at step = ", step
    
    maxd = (/dx, dy, dz/) / bound_safe_factor

    DO ispecies = 1, n_species
      omega = species_list(ispecies)%harmonic_omega
      gamma = species_list(ispecies)%harmonic_gamma
#ifndef PER_PARTICLE_CHARGE_MASS
      part_m   = species_list(ispecies)%mass
#endif
      IF (ALL(omega == 0.0_num)) CYCLE
      current => species_list(ispecies)%attached_list%head

      DO ipart = 1, species_list(ispecies)%attached_list%count
        IF (current%weight == 0.0_num) THEN
          current => current%next
          CYCLE
        END IF
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m = current%mass
#endif
        dlx = current%part_pos - current%part_ip
        
        IF (ANY(ABS(dlx) > maxd)) THEN
          part_u = current%part_p / (part_m * c)
          PRINT '(a)', '***ERROR***'
          PRINT '(a, 3ES13.5)', "dlx value large =", dlx
          PRINT '(A, F5.2, A, 3ES13.5)', "for dxs/", bound_safe_factor, &
               " = ", maxd
          PRINT '(A, I3)',"  rank = ", rank
          PRINT '(a, a)', "  species = ", TRIM(species_list(ispecies)%name)
          PRINT '(a, I4)', "  step = ", step
          PRINT '(a, 3ES13.5)', "  ux,uy,uz= ", part_u
          PRINT '(a, 3ES13.5)', "  part_pos=", current%part_pos
          PRINT '(a, 3ES13.5)', "  part_ip=", current%part_ip
          PRINT '(a, 3ES13.5)', "  omega=", omega
          PRINT '(a, 3ES13.5)', "  gamma=", gamma
          CALL abort_code(c_err_bad_value)
        END IF
        current => current%next
      END DO
    END DO
  END SUBROUTINE check_bound_particles
#endif
  

  ! Background distribution function used for delta-f calculations.
  ! Specialise to a drifting (tri)-Maxwellian to simplify and ensure
  ! zero density/current divergence.
  ! Can effectively switch off deltaf method by setting zero background density.

  FUNCTION f0(ispecies, mass, p)

    INTEGER, INTENT(IN) :: ispecies
    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(:), INTENT(IN) :: p
    REAL(num) :: f0
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz, density
    REAL(num) :: f0_exponent, norm, two_kb_mass, two_pi_kb_mass3
    TYPE(particle_species), POINTER :: species

    species => species_list(ispecies)

    IF (ABS(species%initial_conditions%density_back) > c_tiny) THEN
       two_kb_mass = 2.0_num * kb * mass
       two_pi_kb_mass3 = (pi * two_kb_mass)**3

       Tx = species%initial_conditions%temp_back(1)
       Ty = species%initial_conditions%temp_back(2)
       Tz = species%initial_conditions%temp_back(3)
       driftx  = species%initial_conditions%drift_back(1)
       drifty  = species%initial_conditions%drift_back(2)
       driftz  = species%initial_conditions%drift_back(3)
       density = species%initial_conditions%density_back
       f0_exponent = ((p(1) - driftx)**2 / Tx &
                    + (p(2) - drifty)**2 / Ty &
                    + (p(3) - driftz)**2 / Tz) / two_kb_mass
       norm = density / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
       f0 = norm * EXP(-f0_exponent)
    ELSE
       f0 = 0.0_num
    END IF

  END FUNCTION f0



#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
  SUBROUTINE push_photons(ispecies)

    ! Very simple photon pusher
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: delta_x, delta_y, delta_z
    INTEGER,INTENT(IN) :: ispecies
    TYPE(particle), POINTER :: current

    REAL(num) :: current_energy, dtfac, fac

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    LOGICAL :: probes_for_species
#endif

#ifndef NO_PARTICLE_PROBES
    current_probe => species_list(ispecies)%attached_probes
    probes_for_species = ASSOCIATED(current_probe)
#endif
    dtfac = dt * c**2

    ! set current to point to head of list
    current => species_list(ispecies)%attached_list%head
    ! loop over photons
    DO WHILE(ASSOCIATED(current))
      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
      current_energy = current%particle_energy

      fac = dtfac / current_energy
      delta_x = current%part_p(1) * fac
      delta_y = current%part_p(2) * fac
      delta_z = current%part_p(3) * fac
#ifndef NO_PARTICLE_PROBES
      init_part_x = current%part_pos(1)
      init_part_y = current%part_pos(2)
      init_part_z = current%part_pos(3)
#endif
      current%part_pos = current%part_pos + (/delta_x, delta_y, delta_z/)
#ifndef NO_PARTICLE_PROBES
      final_part_x = current%part_pos(1)
      final_part_y = current%part_pos(2)
      final_part_z = current%part_pos(3)
#endif

#ifndef NO_PARTICLE_PROBES
      IF (probes_for_species) THEN
        ! Compare the current particle with the parameters of any probes in
        ! the system. These particles are copied into a separate part of the
        ! output file.

        current_probe => species_list(ispecies)%attached_probes

        ! Cycle through probes
        DO WHILE(ASSOCIATED(current_probe))
          ! Unidirectional probe
          IF (current_energy > current_probe%ek_min) THEN
            IF (current_energy < current_probe%ek_max) THEN

              d_init  = SUM(current_probe%normal * (current_probe%point &
                  - (/init_part_x, init_part_y, init_part_z/)))
              d_final = SUM(current_probe%normal * (current_probe%point &
                  - (/final_part_x, final_part_y, final_part_z/)))
              IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                ! this particle is wanted so copy it to the list associated
                ! with this probe
                ALLOCATE(particle_copy)
                particle_copy = current
                CALL add_particle_to_partlist(&
                    current_probe%sampled_particles, particle_copy)
                NULLIFY(particle_copy)
              END IF

            END IF
          END IF
          current_probe => current_probe%next
        END DO
      END IF
#endif

      current => current%next
    END DO

  END SUBROUTINE push_photons
#endif
#ifdef BOUND_HARMONIC
  SUBROUTINE clean_diminished_parts
    TYPE(particle_species), POINTER :: species
    TYPE(particle), POINTER :: current, next
    REAL(num), PARAMETER :: epsmin = c_small_weight
    INTEGER i, nremoved

    DO i=1,n_species
      current => species_list(i)%attached_list%head
      nremoved = 0
      DO WHILE (ASSOCIATED(current))
        next => current%next
        IF (current%partner_count == 0 .AND. current%weight == 0) THEN
          ! no partners and zeroweight, physically diminished and
          ! all partners unlinked pointer wise, so we can remove from the simulation
          ! safely
          CALL remove_particle_from_partlist(species_list(i)%attached_list,&
               current)
          CALL destroy_particle(current)
          nremoved = nremoved + 1
        END IF
        current => next
      END DO
      IF (nremoved > 0) THEN
#ifdef VERBOSE_BNDHARM
        PRINT '(A,I2,A,I5.5,A,A)', &
             "rank=",rank," clean_diminished_parts: removed ",nremoved,&
             " from species ", TRIM(species_list(i)%name)
#endif
      END IF
    END DO ! species
  END SUBROUTINE clean_diminished_parts

  SUBROUTINE check_spec_cell_z(label)
    INTEGER, INTENT(IN) :: label
    TYPE(particle_species), POINTER :: species
    TYPE(particle), POINTER :: current, next
    REAL(num) d(3)
    INTEGER i, ncell(3)
    LOGICAL found
    found = .FALSE.
    current => species_list(1)%attached_list%head
    i = 1
    DO WHILE (ASSOCIATED(current))
      next => current%next
      d = current%part_pos - &
           [x_grid_min_local,y_grid_min_local,z_grid_min_local]
      d = d / [dx,dy,dz] + 1.5_num
      ncell = FLOOR(d)

      IF (ncell(3) > 50) THEN
        PRINT "(A,I2,A,I3,A,3I4,A,I4)", &
             "rank=",rank,": label=",label,", ncell = ",ncell,&
             ", i = ",i
        found = .TRUE.
      END IF
      current => next
      i = i + 1
    END DO
    IF (.NOT. found) PRINT "(A,I2,A,I3)", "rank=",rank,&
         ": check_spec_cell_z ok for label=",label
  END SUBROUTINE check_spec_cell_z
#endif

END MODULE particles
