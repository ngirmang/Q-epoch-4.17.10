MODULE perfmon
#ifdef PERFMON2D
#ifndef MEDIUM
#error  "PERFMON2D requires MEDIUM and CONSTEPS"
#endif

  USE constants

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: Nperfpts = c_perfpts, Nperfsamps = 10, Nperfstops=c_nperfstops
  INTEGER, PARAMETER :: perfstops(Nperfstops) = [3,5,10]
  INTEGER :: perfsamps(Nperfsamps+1, Nperfpts)
  REAL(num) :: &
    perfavg(Nperfstops, Nperfpts), &
    perfrms(Nperfstops, Nperfpts), &
    perfmed(Nperfstops, Nperfpts), &
    perfmin(Nperfstops, Nperfpts), &
    perfmax(Nperfstops, Nperfpts)
  INTEGER, PARAMETER :: Nstats = 5

  PUBLIC :: init_perfmon
  PUBLIC :: perfmon_update, perfavg, perfrms, perfmed, perfmin, perfmax
  PUBLIC :: perfmon_start, perfmon_stop, Nperfstops, Nperfpts

CONTAINS

  SUBROUTINE subsortlin(arr)

    INTEGER, DIMENSION(:) :: arr

    INTEGER i, j, nv

    IF (SIZE(arr) .LT. 2) RETURN

    DO i = 1, SIZE(arr)-1
      DO j = i + 1, SIZE(arr)
        IF (arr(i) .GT. arr(j)) THEN
          nv = arr(i)
          arr(i) = arr(j)
          arr(j) = nv
        END IF
      END DO
    END DO

  END SUBROUTINE subsortlin



  SUBROUTINE init_perfmon

    perfsamps = 0

  END SUBROUTINE init_perfmon



  SUBROUTINE perfmon_update

    INTEGER :: i, n, nsi, nsampstop, nsamps, nst, nmed
    INTEGER :: subsamp(Nperfsamps)
    REAL(num) :: v

    nperfptsloop: DO n = 1, Nperfpts
      !! update all statistics
      nsamstoploop: DO nsi = 1, Nperfstops
        nsampstop = perfstops(nsi)
        nst = Nperfsamps-nsampstop+1
        ! average
        v = 0.0_num
        DO i = nst, Nperfsamps
          v = v + REAL(perfsamps(i,n),num)
        END DO
        v = v / nsampstop
        perfavg(nsi, n) = v

        ! rms
        v = 0.0_num
        DO i = nst, Nperfsamps
          v = v + REAL(perfsamps(i,n),num)**2
        END DO
        perfrms(nsi, n) = SQRT(v / nsampstop)

        ! median
        subsamp(1:nsampstop) = perfsamps(nst:Nperfsamps, n)

        CALL subsortlin(subsamp(1:nsampstop))
        nmed = nsampstop / 2 + MOD(nsampstop, 2)

        perfmed(nsi, n) = subsamp(nmed)
        ! minimum
        perfmin(nsi, n) = subsamp(1)

        ! maximum
        perfmax(nsi, n) = subsamp(nsampstop)

      END DO nsamstoploop
      !! shift
      DO i = Nperfsamps, 1, -1
        perfsamps(i+1, n) = perfsamps(i, n)
      END DO

    END DO nperfptsloop
    
  END SUBROUTINE perfmon_update



  SUBROUTINE perfmon_start(ilabel)

    INTEGER, INTENT(IN) :: ilabel
    INTEGER :: ist

    CALL system_clock(ist)
    perfsamps(1, ilabel) = ist

  END SUBROUTINE perfmon_start



  SUBROUTINE perfmon_stop(ilabel)

    INTEGER, INTENT(IN) :: ilabel
    INTEGER :: i

    CALL system_clock(i)

    i = i - perfsamps(1, ilabel)
    perfsamps(1, ilabel) = i

  END SUBROUTINE perfmon_stop

#endif
END MODULE perfmon
