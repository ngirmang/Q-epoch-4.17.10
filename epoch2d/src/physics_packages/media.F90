! Copyright (C) 2025 Gregory K. Ngirmang
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


MODULE media
#ifdef MEDIUM

  USE constants
  USE shared_data
  USE random_generator, only : random
  USE utilities
  USE partlist
  USE random_generator

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ionise_media, initialise_media, media_particle_production,&
       update_medium_eps, deallocate_media

  REAL(num), PARAMETER :: c_sq = c**2, e_mass = m0, e_restE = m0*c_sq

  ! relevant parameters for PPT
  INTEGER ii              ! 0.2  1.0    20   200
  INTEGER, PARAMETER :: ngm=1  + 400 + 380 + 360
  REAL(num), PARAMETER :: gms(ngm) = [&
         0.2_num, &
       ( 0.2_num + 0.002_num*ii, ii=1,400), &
       ( 1.0_num + 0.050_num*ii, ii=1,380), &
       (20.0_num + 0.500_num*ii, ii=1,360) ]

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: ppt_rates

  ! convience look up to save a few multiplies to easily calculate the
  ! keldysh gamma parameter
  REAL(num), DIMENSION(:), ALLOCATABLE :: keldysh_gmf

  ! relevant parameters for high field keldysh
  REAL(num), DIMENSION(:), ALLOCATABLE :: keldysh_prefs
  REAL(num), DIMENSION(:), ALLOCATABLE :: keldysh_Fi
  REAL(num), DIMENSION(:), ALLOCATABLE :: keldysh_fexp

  ! collision ionisation (taken from collisions.F90)
  REAL(num), DIMENSION(3,0:2), PARAMETER :: a_bell = RESHAPE( &
      (/ 0.5250_num, 0.5300_num, 0.1300_num, &
         0.0000_num, 0.6000_num, 0.3880_num, &
         0.0000_num, 0.0000_num, 0.3500_num /) * 1e-13_num, (/3,3/) )

  REAL(num), DIMENSION(3,0:2,7), PARAMETER :: b_bell = RESHAPE( &
      (/-0.5100_num, -0.4100_num,  0.2500_num, &
         0.0000_num, -0.4000_num, -0.2000_num, &
         0.0000_num,  0.0000_num,  1.6000_num, &
         0.2000_num,  0.1500_num, -1.5000_num, &
         0.0000_num, -0.7100_num, -0.2356_num, &
         0.0000_num,  0.0000_num, -3.0000_num, &
         0.0500_num,  0.1500_num,  2.4000_num, &
         0.0000_num,  0.6550_num,  0.5355_num, &
         0.0000_num,  0.0000_num,  4.0000_num, &
        -0.0250_num, -0.2000_num,  3.2200_num, &
         0.0000_num,  0.4250_num,  3.1500_num, &
         0.0000_num,  0.0000_num,  2.0000_num, &
        -0.1000_num, -0.1500_num, -3.6670_num, &
         0.0000_num, -0.7500_num, -8.5000_num, &
         0.0000_num,  0.0000_num, -5.0000_num, &
         0.0000_num,  0.0000_num,  0.0000_num, &
         0.0000_num,  0.0000_num,  5.0500_num, &
         0.0000_num,  0.0000_num, -1.5000_num, &
         0.0000_num,  0.0000_num,  0.0000_num, &
         0.0000_num,  0.0000_num,  0.3700_num, &
         0.0000_num,  0.0000_num,  3.5000_num /) * 1e-13_num, (/3,3,7/) )

  REAL(num), DIMENSION(0:2), PARAMETER :: &
      l_bell = (/ 1.27_num, 0.542_num, 0.95_num /) * 1e-13_num

  REAL(num) :: omega = -1.0_num

CONTAINS

  ! my factorial which is better, screw loops

  REAL(num) FUNCTION sfact(i)
    INTEGER, INTENT(IN) :: i

    REAL(num), PARAMETER :: facts(0:9) = (/ &
         1.0_num, 1.0_num, 2.0_num, 6.0_num, 24.0_num, &
         120.0_num, 720.0_num, 5040.0_num, 40320.0_num, &
         362880.0_num /)
    
    IF     (i .LT. 0) THEN
      sfact = 0.0_num
    ELSEIF (i .GT. 9) THEN
      sfact = gamma(REAL(i-1,num))
    ELSE
      sfact = facts(i)
    END IF

  END FUNCTION sfact

  ! calculate ppt rates, currently only the unpolarised variant

  SUBROUTINE calculate_ppt_rates(rates, &
       Qi, omega, U0, l, ma)

    REAL(num), INTENT(OUT), DIMENSION(ngm) :: rates
    REAL(num), INTENT(IN) :: Qi, omega, U0
    INTEGER, INTENT(IN) :: l, ma

    REAL(num), PARAMETER :: conv_thresh = 1d-8, vst=1d-5

    REAL(num), PARAMETER :: &
         beta(ngm) = 2*gms / sqrt(1 + gms**2), &
         alpha(ngm)= 2*asinh(gms) - beta, &
         gs(ngm) = 3.0_num / (2*gms) * ( &
           (1 + 1.0_num/(2*gms**2))*asinh(gms) - 1/beta)


    REAL(num) :: nus(ngm), dels(ngm), A0s(ngm), Es(ngm), rate_pref(ngm)
    REAL(num) :: eff_n, dkap, vm1, vm2, vacc, v, w0, pref, x, omi, Fi


    ! you should converge by 10000s
    INTEGER, PARAMETER :: nterm_max = 100000, nterm_min = 50, &
         nw = 9999 ! 10K ought to be enough for anyone
    REAL(num) :: w0_f(0:nw), ws(0:nw), dw, wst
    INTEGER :: i,j,n

    IF (ma .GT. 0 ) THEN
      IF (rank == 0) PRINT '("ppt for |m| > 0 not implemented yet!")'
      CALL abort_code(c_err_bad_setup)
    END IF

    ! calculate dkappa := U0 / (h_bar*omega)
    
    dkap = U0 / (h_bar*omega)
    eff_n = Qi/q0 / sqrt(2*U0 / hartree)

    IF (NINT(eff_n) - eff_n .lt. 1e-9) eff_n = NINT(eff_n)

    Fi = atomic_electric_field * (2*U0/hartree) ** 1.5_num
    Es = omega*sqrt(2*m0*U0)/(q0*gms)
    omi= U0 / h_bar

    rate_pref = omi * &
         sqrt(6/pi)*(2.0_num*l + 1)*sfact(l+ma)/sfact(l-ma)&
         / &
         (sfact(ma)*2.0_num**ma)

    ! ADK Cml^2 approximation
    rate_pref = rate_pref*&
         2.0_num**(2*eff_n) &
         / &                  ! why not just use
         eff_n*gamma(2*eff_n) ! intrinsic gamma??

    ! front field factor
    rate_pref = rate_pref* &
         (Es*sqrt(1+gms**2)/(2*Fi)) ** (ma + 1.5_num - 2 * eff_n)

    nus  = dkap*(1 + 0.5_num/gms**2)
    dels = nus - floor(nus)

    DO i=1,ngm
      vm2 = 2.0_num
      vm1 = 2.0_num
      vacc = 0

      pref = (4.0_num/sqrt(3*pi))*gms(i)**2/(1+gms(i)**2)

      DO n=1,nterm_max
        x  = sqrt(beta(i)*(n-dels(i)))
      
        ! for large x, we just use the asymptotic expansion for dawsn
        IF ( x .gt. 30) THEN
        ! asymptotic expansion for dawsn
          w0 = 1 / (2*x) + 1 / (4*x**3) + 3 / (8*x**5) + 15 / (16*x**7)

          ! ignore parts of the integrand that is well near zero so we
          ! can focus on the region near the peak which is always at
          ! the right edge of the integrand region
          
          ! calculating the w start
          !
          ! exp(wst**2 - x**2) == vst
          ! wst**2 - x**2 == log(vst)
          ! wst = sqrt(log(vst) + x**2)

          !  check that exp(-x**2) isn't greater than vst
        ELSEIF ( exp(-x**2) .gt. vst ) THEN
          wst = 0
          dw = (x - wst) / nw
          ! now allocate the integrand
          ! we do it in reverse order, knowing the convergence of
          ! the integral
          ws = [(x - dw*j, j=0, nw )]
          w0_f = exp(ws**2 - x**2)
          w0 = simps(w0_f, nw, dw)
          ! otherwise, use the vst calculation in the above comment
        ELSE
          wst = sqrt(log(vst)+x**2)
          dw = (x - wst) / nw
          ws = [(x - dw*j, j=0, nw )]
          w0_f = exp(ws**2 - x**2)
          w0 = simps(w0_f, nw, dw)
        END IF

        v  = w0*exp(-alpha(i)*(n-dels(i)))*pref
        ! accumulate the v's
        vacc = vacc + v

        ! check convergence
        IF ( 2*abs(vacc-vm2)/abs(vacc + vm2) .lt. conv_thresh &
             .AND. n .ge. nterm_min ) THEN! converged!
          EXIT
        !ELSEIF (n > 1000 .AND. MOD(n,1000) == 0) THEN
        !  PRINT '(" iteration n=",I7.7," x=", ES12.6,&
        !     " w0=",ES12.6, " v=",ES12.6)',&
        !     n, x, w0, v
        !  CONTINUE
        END IF
        vm2 = vm1
        vm1 = vacc

      END DO
      A0s(i) = vacc
      rates(i) = rate_pref(i)*A0s(i)*exp(-2.0_num/3.0_num*Fi/Es(i)*gs(i))
      PRINT '("gm=", ES10.4, ", A0=", ES10.4, ", rate=",ES10.4)',&
           gms(i), A0s(i), rates(i)
    END DO

    !rates = rate_pref*A0s*exp(-2.0_num/3.0_num*Fi/Es*gs)

  END SUBROUTINE calculate_ppt_rates



  ! simpsons' 3/8th rule with fast2sum
  ! assumes ny is divisible by 3
  REAL(num) FUNCTION simps(ds, ny, dy)

    INTEGER, INTENT(IN) :: ny
    REAL(num),INTENT(IN) :: ds(0:ny), dy

    REAL(num) s, z, err, rn, y
    INTEGER i

    simps = ds(0)
    err = 0.d0

    DO i = 1, ny-1
      rn = 3.0_num
      if (mod(i,3) == 0) rn = 2.0_num
      ! fast2sum
      y = rn*ds(i) + err
      s = simps + y
      z = s - simps
      err = y - z
      simps = s
    END DO
    ! add last term
    y = ds(ny) + err
    simps = simps + y
    simps = simps * dy * 0.375_num ! 3/8

  END FUNCTION simps



  SUBROUTINE calculate_keldysh_quantities(&
       pref, Fi, ofexp, Qi, omega, U0, l, ma)

    REAL(num), INTENT(OUT) :: pref, Fi, ofexp
    
    REAL(num), INTENT(IN) :: Qi, omega, U0
    INTEGER, INTENT(IN) :: l, ma

    REAL(num) :: eff_n, dkap, &
         vm1, vm2, vacc, v, &
         w0, x, omi

    INTEGER :: i,n

    eff_n = Qi/q0 / sqrt(2*U0 / hartree)

    omi = U0/h_bar
    pref = omi * &
         sqrt(6/pi)*(1.0_num + 2*l)*sfact(l+ma)/sfact(1-ma) &
         / &
         (sfact(ma)*2.0_num**ma)
    ! ADK Cml^2 approximation
    pref = pref*&
         2.0_num**(2*eff_n) &
         / &                  ! why not just use
         eff_n*gamma(2*eff_n) ! intrinsic gamma??
    Fi = atomic_electric_field * (2*U0/hartree) ** 1.5_num
    ofexp = ma + 1.5_num - 2*eff_n

  END SUBROUTINE calculate_keldysh_quantities



  SUBROUTINE dump_rates(im)

    INTEGER, INTENT(IN) :: im
    INTEGER ih

    ih = im + 100

    IF (rank .NE. 0) RETURN

    PRINT '("dumping ionisation rates")'

    OPEN(ih, file=media_list(im)%gamma_file, access='stream')

    WRITE(ih) ngm
    WRITE(ih) gms

    CLOSE(ih)

    ih = im + 200

    OPEN(ih, file=media_list(im)%ionisation_file, access='stream')

    WRITE(ih) omega
    WRITE(ih) ngm
    WRITE(ih) ppt_rates

    CLOSE(ih)

  END SUBROUTINE dump_rates

  ! Do most of the calculations on a per-species basis once at start of
  ! simulation to save on computational time

  SUBROUTINE initialise_media

    INTEGER :: i, iu, io, err_laser, im
    LOGICAL :: laser_set, any_field_ionise
    TYPE(laser_block), POINTER :: current_laser
    INTEGER :: ncur_species, next_species
    INTEGER :: l, ma
    REAL(num) :: U0, ion_charge

    ! look for a useable laser for an omega
    IF (n_media == 0) RETURN

    err_laser = 0
    laser_set = .FALSE.
    DO i=1,4
      SELECT CASE(i)
      CASE (1)
        current_laser => laser_x_min
      CASE (2)
        current_laser => laser_x_max
      CASE (3)
        current_laser => laser_y_min
      CASE (4)
        current_laser => laser_y_max
      END SELECT
      
      IF (.NOT. ASSOCIATED(current_laser)) CYCLE
      
      IF (laser_set .AND. ABS(current_laser%omega - omega) > c_tiny) THEN
        err_laser = 1
      ELSE
        omega = current_laser%omega
        laser_set = .TRUE.
      END IF
      DO WHILE (ASSOCIATED(current_laser%next))
        IF (ABS(current_laser%omega - omega) > c_tiny) THEN
          err_laser = 1
          EXIT
        END IF
        current_laser => current_laser%next
      END DO
    END DO
    IF (.NOT. laser_set) err_laser = 2

    IF (err_laser == 1) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Medium field ioniser model does not currently'
          WRITE(io,*) 'support lasers of differing frequencies attached to'
          WRITE(io,*) 'the boundaries. Please adjust your input deck.'
        END DO
      END IF
      CALL abort_code(c_err_bad_setup)
    ELSE IF (err_laser == 2) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Medium field ioniser found no laser, so we are assuming'
          WRITE(io,*) 'omega=0, that is, using Landau style field ionisation.'
          WRITE(io,*) 'If this is not what you intend, fix your deck.'
        END DO
      END IF
      omega = 0
    END IF

    !! load ionisation details
    !
    ! in the current iteration of this model for ppt field ionisation,
    ! the projection m quantum number is averaged over in the current subshell,
    ! because I'm going to assume I do NOT track quantum levels for each atom.

    ! Thus, with l=0, ofc we only have m = 0, for p-filled subshells, we have
    ! 
    ! w_0 1/3 + 2*w_|1|/3,
    ! and for d subshells,
    ! w_0 1/5 + 2*w_|1|/5 + 2*w_|2|/5

    ! ionised_species_present = .FALSE.

    IF (rank == 0) THEN
      PRINT '("initialising ",I2," media models")', n_media
    END IF
    ALLOCATE(keldysh_gmf(n_media))
    keldysh_gmf = 1.0_num
    ALLOCATE(ppt_rates(ngm, n_media))
    ppt_rates = 0.0_num

    ALLOCATE(keldysh_prefs(n_media),&
         keldysh_Fi(n_media),&
         keldysh_fexp(n_media))

    DO im=1, n_media
      ncur_species = media_list(im)%species

      IF (.NOT. species_list(ncur_species)%ionise) CYCLE

      next_species = species_list(ncur_species)%ionise_to_species
      ion_charge = species_list(next_species)%charge

      l = species_list(ncur_species)%l

      U0 = species_list(ncur_species)%ionisation_energy

      keldysh_gmf(im) = omega*sqrt(2*U0*m0)/q0

      ma = 0 !for now, we don't do m...

      CALL calculate_ppt_rates(&
           ppt_rates(:,im), ion_charge, omega, U0, l, ma)

      ! in the very low gamma limit (here, \gamma < 0.2
      ! we do "keldysh" ionisation, really Eq. 59 of PPT with
      ! the ADK approximation for C_nl and the extra factor on
      ! the field prefactor

      CALL calculate_keldysh_quantities(&
           keldysh_prefs(im), &
           keldysh_Fi(im), &
           keldysh_fexp(im), &
           ion_charge, omega, U0, l, ma)
      IF ( media_list(im)%dump_ionisation_rates ) CALL dump_rates(im)
    END DO

  END SUBROUTINE initialise_media



  SUBROUTINE deallocate_media

    INTEGER :: stat

    DEALLOCATE(keldysh_gmf, stat=stat)
    DEALLOCATE(ppt_rates, stat=stat)
    DEALLOCATE(keldysh_prefs, stat=stat)
    DEALLOCATE(keldysh_Fi, stat=stat)
    DEALLOCATE(keldysh_fexp, stat=stat)

  END SUBROUTINE deallocate_media


  ! ppt field ionisation
  FUNCTION ppt_ionise(E, im) result(rate)

    REAL(num) :: rate
    REAL(num), INTENT(IN) :: E
    INTEGER, INTENT(IN) :: im

    INTEGER :: ig, i, lowend, highend
    REAL(num) :: mygamma, U0, x

    mygamma = keldysh_gmf(im)/E

    IF (mygamma .GT. gms(ngm) ) THEN
      rate = 0.0_num
      RETURN
    ELSE IF (mygamma .LT. gms(1)) THEN
      rate = keldysh_ionise_species(E, mygamma, im)
      RETURN
    END IF
    ! gamma is GT small and LT 200, so in the gms array

    ! binary search on gms
    lowend = 1
    highend = ngm
    ig = highend / 2
    outer: DO WHILE (.TRUE.)
      IF (gms(ig) .le. mygamma .AND. mygamma .lt. gms(ig+1)) THEN
        EXIT outer
      ELSE IF (highend - lowend .LE. 4) THEN
        ! just linear search when it's just 4 items
        inner:  DO ig=lowend, highend-1
          IF (gms(ig) .le. mygamma .AND. mygamma .lt. gms(ig+1)) THEN
            EXIT outer
          END IF
        END DO inner
        IF (rank == 0) PRINT '("huh???")'
        CALL abort_code(c_err_generic_error)
      END IF
      IF (gms(ig) .gt. mygamma) THEN
        highend = ig
      ELSE
        lowend = ig
      END IF
      ig = (highend + lowend) / 2
      
    END DO outer

    ! linear interpolate
    x = (gms(ig+1) - mygamma) / (gms(ig+1) - gms(ig))
    rate = ppt_rates(ig,im)*x + ppt_rates(ig+1, im)*(1-x)

  END FUNCTION ppt_ionise

  FUNCTION keldysh_ionise_species(E, gm, species) result(rate)
    REAL(num) :: rate
    REAL(num), INTENT(IN) :: E, gm
    INTEGER, INTENT(IN) :: species

    REAL(num) :: Fi, fexp
    Fi = keldysh_Fi(species)

    rate = keldysh_prefs(species)*(E/(2*Fi))**fexp * &
         exp(-2.0_num*Fi/(3.0_num*E) * (1 - gm**2 / 10.0_num))
    
  END FUNCTION keldysh_ionise_species
  ! Redirect to correct routine; the only difference between these routines is
  ! in the rate calculated but it saves doing this logic for each ionisation

  SUBROUTINE ionise_media

    INTEGER :: im, ne, i, ix, iy, cur_spec, elec_spec, next_spec, nextm, &
         ncreate, nmax_create, ntmp, nl, nn, next_me, &
         ncreate_total
    INTEGER :: ncreate_totalr, ierr

    REAL(num) enorm, rate, dN, cur_N, U0, dN_field, &
         next_create_min, next_create_max, next_weight
    REAL(num) enorm_max, min_N, max_N, rate_max
    REAL(num) min_Nr, max_Nr, enorm_maxr, rate_maxr

    TYPE(particle), POINTER :: cur, newe, newi, next
    TYPE(medium), POINTER :: cur_medium, next_medium

    LOGICAL :: next_media, elec_media, use_probsm

    ! collisional ionisation
    INTEGER(8) :: nelec
    INTEGER :: npart
    REAL(num) :: e_p2_st, e_ke_st, e_v_st, e_ke_en, gr, red_inc, red_ion, &
         red, fion, eiics, ion_charge, full_ion_charge, my_dN
    REAL(num) :: e_p(3)

    DO im=1,n_media

      cur_medium => media_list(im)
      cur_spec = cur_medium%species

      ncreate_total = 0
      enorm_max = 0.0_num
      rate_max  = 0.0_num

      min_N = HUGE(1.0_num)
      max_N = 0.0_num

      IF (.NOT. species_list(cur_spec)%ionise) CYCLE

      next_spec = species_list(cur_spec)%ionise_to_species
      next_media = species_list(next_spec)%medium_species
      next_create_min = media_list(im)%next_create_min
      nmax_create = media_list(im)%nmax_create
      use_probsm = media_list(im)%use_prob

      IF (next_media) &
           nextm = species_list(next_spec)%medium_index

      elec_spec = species_list(cur_spec)%release_species
      elec_media = species_list(elec_spec)%medium_species

      IF (elec_media) next_me = ielectron_medium

      U0 = species_list(cur_spec)%ionisation_energy

      IF (cur_medium%use_collisional_ionisation) THEN
        ntmp = next_spec
        DO WHILE(species_list(ntmp)%ionise)
          ntmp = species_list(ntmp)%ionise_to_species
        END DO
        full_ion_charge = species_list(ntmp)%charge

        red_ion= e_restE / U0 ! we stay in SI
        ion_charge = species_list(next_spec)%charge
      END IF
      DO iy = 1,ny
ixlp: DO ix = 1,nx

        cur_N = media_density(ix,iy,im)
        ncreate = 0
        dN  = 0

        IF (cur_medium%use_field_ionisation) THEN

          ! need to set up temporal averaging
          enorm = SQRT( 0.25_num*(ex(ix,iy) + ex(ix-1,iy  ))**2 &
                     +  0.25_num*(ey(ix,iy) + ey(ix  ,iy-1))**2 &
                              +   ez(ix,iy)**2)

          IF (enorm > enorm_max) enorm_max = enorm

          rate  = ppt_ionise(enorm, im)

          ! don't touch dN if enorm is too small
          IF (enorm .GE. 1.0e2_num) THEN

            dN = rate * dt * cur_N

            IF (use_probsm .AND. dN .LT. next_create_min) THEN
              IF (random() .LT. 1.0_num - exp(-rate*dt)) THEN
                dN = next_create_min
              ELSE
                dN = 0.0_num
              END IF
            END IF

            IF (rate .GT. rate_max) rate_max = rate

          END IF ! enorm .lt. 1d2

        END IF ! use field ionisation

        IF (cur_medium%use_collisional_ionisation) THEN
          nelec = species_list(elec_spec)%secondary_list(ix,iy)%count
          cur => species_list(elec_spec)%secondary_list(ix,iy)%head

          DO ne=1,nelec
            e_p = cur%part_p
            e_p2_st = e_p(1)**2 + e_p(2)**2 + e_p(3)**2
            e_ke_st = SQRT(e_p2_st*c_sq + (e_restE)**2) - e_restE
            e_v_st  = SQRT(e_p2_st / (e_mass**2 + e_p2_st / c**2))

            IF (e_ke_st < U0) THEN
              cur => cur%next
              CYCLE ixlp
            END IF
            ! we only implement this for lower Zm and thus restrict to 
            ! A<36 and just use MBELL cross sections (see collisions.F90)
            red_inc = e_ke_st / U0
            red = red_inc + red_ion

            ! relativistic correction (just copied from collisions.F90,
            ! need to look this up)
            gr = &
(2*red_ion + 1) / (2*red_ion + red_inc) &
* (red/(1+red_ion))**2 * ( &
  (1+red_inc)*(red+red_ion)*(1+red_ion)**2 &
/ (red_ion**2*(2*red_ion+1) + red_inc*(red + red_ion)*(red_ion+1))&
  ) ** 1.5_num
            fion = 1 + &
              3.0_num*(red_inc*ion_charge/full_ion_charge)**l_bell(nl)
            eiics = 0
            DO i=1,7
              eiics = eiics + b_bell(nn,nl,i)*(1 - 1.0_num / red_inc)**i
            END DO
            
            eiics = (a_bell(nn,nl)*LOG(red_inc) + eiics) / (e_ke_st*U0)
            
            eiics = eiics * gr * fion

            rate = eiics*1e-4*e_v_st * cur%weight/dx/dy

            my_dN = rate*dt*cur_N
            dN = dN + my_dN

            !finally, slow electron
            e_ke_en = e_ke_st - my_dN*U0
            e_p = e_p*sqrt( ((e_ke_en/e_restE)**2 - 1) / e_p2_st) * m0*c

            cur%part_p = e_p

            cur => cur%next
          END DO
        END IF ! collisional ionisation


        IF (cur_medium%quantised) &
          dN = FLOOR(dN/next_create_min)*next_create_min

        ! make sure to subtract at most the remaining medium density
        dN = MIN(dN, cur_N)

        IF(.NOT. next_media .OR. .NOT. elec_media) THEN
          IF (dN .GE. next_create_min) THEN
            ncreate = NINT(dN / next_create_min)
            next_weight = next_create_min

            IF (ncreate .GT. nmax_create) THEN
              next_weight = dN / nmax_create
              ncreate = nmax_create
            END IF
          ELSE
            ncreate = 0
            dN = 0
          END IF
        END IF

        IF ( dN .GT. 0.0_num .AND. ncreate .GT. 0 ) THEN
          IF (.NOT. next_media .AND. .NOT. elec_media ) THEN
            ! create ions and electrons with same random positions
            CALL create_media_ions_electrons(ncreate, next_weight, &
                 next_spec, elec_spec, ix, iy)
          ELSE
            !ion species
            IF (next_media) THEN
              media_density(ix,iy,nextm) = media_density(ix,iy,nextm)+dN
            ELSE IF (ncreate > 0) THEN
              CALL create_media_particles(ncreate, next_weight,&
                   next_spec, ix, iy)
            END IF
            ! electron species
            IF (elec_media) THEN
              media_density(ix,iy,next_me) = &
                   media_density(ix,iy,next_me) + dN
            ELSE IF (ncreate > 0) THEN
              CALL create_media_particles(ncreate, next_weight,&
                   elec_spec, ix, iy)
            END IF
          END IF
        END IF
        media_density(ix,iy,im) = cur_N - dN

        ncreate_total = ncreate_total + ncreate

        IF      ( media_density(ix,iy,im) < min_N ) THEN
          min_N = media_density(ix,iy,im)
        ELSE IF ( media_density(ix,iy,im) > max_N ) THEN
          max_N = media_density(ix,iy,im)
        END IF

      END DO ixlp
      END DO

      CALL MPI_REDUCE(ncreate_total, ncreate_totalr, 1, MPI_INTEGER, &
           MPI_SUM, 0, comm, ierr)
      CALL MPI_REDUCE(min_N, min_Nr, 1, mpireal, MPI_MIN, 0, comm, ierr)
      CALL MPI_REDUCE(max_N, max_Nr, 1, mpireal, MPI_MAX, 0, comm, ierr)
      CALL MPI_REDUCE(enorm_max, enorm_maxr, 1, mpireal, MPI_MAX, 0, &
           comm, ierr)
      CALL MPI_REDUCE(rate_max, rate_maxr, 1, mpireal, MPI_MAX, 0, &
           comm, ierr)

      IF (rank == 0 .AND. stdout_frequency .GT. 0 .AND. &
        MOD(step, stdout_frequency) == 0) THEN

        PRINT '("media_ionise: step=",I7.7,", ncreate=",I5.5,", peak E=",&
             & ES9.3,", max_N=",ES9.3,", min_N=",ES13.7)',&
             step, ncreate_totalr, enorm_maxr, max_Nr, min_Nr
        PRINT '("                            rate_max=",ES9.3)', &
             rate_maxr

      END IF

    END DO ! each medium

  END SUBROUTINE ionise_media

  ! based on user input, break media into macroparticles
  SUBROUTINE media_particle_production
    
    INTEGER :: n, ix, iy, ncreate, cur_spec
    TYPE(medium), POINTER :: cur_medium
    REAL(num) :: create_min, cur_N

    DO n=1,n_media
      cur_medium => media_list(n)
      cur_spec = cur_medium%species
      create_min = cur_medium%particle_create_density

      DO iy = 1, ny
      DO ix = 1, nx

        cur_N = media_density(ix, iy, n)
        IF (cur_N .GT. create_min) THEN
          ncreate = FLOOR(cur_N/create_min)
          CALL create_media_particles(ncreate, create_min, cur_spec,&
                                      ix, iy)
          media_density(ix,iy,n) = cur_N - ncreate*create_min
        END IF

      END DO
      END DO

    END DO

  END SUBROUTINE media_particle_production



  SUBROUTINE create_media_particles(ncreate, weight, species, ix, iy)

    INTEGER, INTENT(IN) :: ncreate, species, ix, iy
    REAL(num), INTENT(IN) :: weight

    INTEGER :: n
    TYPE(particle), POINTER :: cur

    DO n=1,ncreate
      CALL create_particle(cur)
      cur%weight   = weight*dx*dy
      cur%part_pos = [ xb(ix) + random()*dx,&
                       yb(iy) + random()*dy,&
                       0.0_num ]
      cur%part_ip  = cur%part_pos
      CALL add_particle_to_partlist(&
           species_list(species)%secondary_list(ix,iy), cur)

    END DO

  END SUBROUTINE create_media_particles



  SUBROUTINE create_media_ions_electrons(ncreate, weight,&
       ion_species, elec_species, ix, iy)

    INTEGER, INTENT(IN) :: ncreate, ion_species, elec_species, ix, iy
    REAL(num), INTENT(IN) :: weight

    INTEGER :: n
    TYPE(particle), POINTER :: cur_ion, cur_elec

    DO n=1, ncreate
      CALL create_particle(cur_ion)
      cur_ion%weight   = weight*dx*dy
      cur_ion%part_pos = [xb(ix) + random()*dx,&
                          yb(iy) + random()*dy,&
                          0.0_num ]
      cur_ion%part_ip  = cur_ion%part_pos
      CALL add_particle_to_partlist(&
           species_list(ion_species)%secondary_list(ix,iy), cur_ion)
      CALL create_particle(cur_elec)
      cur_elec%weight   = cur_ion%weight
      cur_elec%part_pos = cur_ion%part_pos
      cur_elec%part_ip  = cur_ion%part_ip
      CALL add_particle_to_partlist(&
           species_list(elec_species)%secondary_list(ix,iy), cur_elec)

    END DO

  END SUBROUTINE create_media_ions_electrons



  SUBROUTINE update_medium_eps

    INTEGER ix,iy,n
    REAL(num) N_c

    IF (ielectron_medium == -1) RETURN

    n = ielectron_medium

    N_c = omega**2*epsilon0*m0/q0**2

    epsx = epsx - media_density(:,:,n) / N_c
    epsy = epsy - media_density(:,:,n) / N_c
    epsz = epsz - media_density(:,:,n) / N_c
    
  END SUBROUTINE update_medium_eps

#endif
END MODULE media
