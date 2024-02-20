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

MODULE timer

  USE constants
  USE mpi
#ifdef PERFMON
  USE shared_data
#endif
  
  IMPLICIT NONE

  LOGICAL :: timer_collect
  REAL(num), PARAMETER :: avg_weight1 = 0.4_num
  REAL(num), PARAMETER :: avg_weight2 = 1.0_num - avg_weight1
  INTEGER, PARAMETER :: c_timer_step = 1
  INTEGER, PARAMETER :: c_timer_dt = 2
  INTEGER, PARAMETER :: c_timer_io = 3
  INTEGER, PARAMETER :: c_timer_balance = 4
  INTEGER, PARAMETER :: c_timer_max = 4

  REAL(num) :: timer_walltime
  REAL(num) :: timer_first(c_timer_max)
  REAL(num) :: timer_time(c_timer_max)
  REAL(num) :: timer_average(c_timer_max)
  LOGICAL :: timer_avg_first(c_timer_max)
#ifdef PERFMON

  INTEGER, PARAMETER :: nperf_avgsteps = 10
  INTEGER, PARAMETER :: &
       nperf_ebc_mpi  = 1, &
       nperf_ebc_zero = 2, &
       nperf_bbc_mpi  = 3, &
       nperf_bbc_zero = 4, nperf_num = 4
  
  REAL(num), PARAMETER :: c_perf_avgsteps = REAL(nperf_avgsteps, num)
  
  INTEGER nperf_steps
  INTEGER, PRIVATE :: nperf_curstep
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: timer_perf_times
  REAL(num), PRIVATE, DIMENSION(:,:), ALLOCATABLE :: perf_avg
  LOGICAL, PRIVATE :: perfticks(nperf_num) = .FALSE., perfpause = .TRUE.
#endif

CONTAINS

  SUBROUTINE timer_init

    timer_first = 0.0_num
    timer_time = 0.0_num
    timer_average = 0.0_num
    timer_avg_first = .TRUE.
    timer_collect = .FALSE.
    
  END SUBROUTINE timer_init



  SUBROUTINE timer_start(id, use_old)

    INTEGER, INTENT(IN) :: id
    LOGICAL, INTENT(IN), OPTIONAL :: use_old

    IF (PRESENT(use_old)) THEN
      IF (.NOT.use_old) timer_walltime = MPI_WTIME()
    ELSE
      timer_walltime = MPI_WTIME()
    END IF

    timer_first(id) = timer_walltime

  END SUBROUTINE timer_start



  SUBROUTINE timer_stop(id)

    INTEGER, INTENT(IN) :: id

    timer_walltime = MPI_WTIME()

    timer_time(id) = timer_walltime - timer_first(id)
    IF (timer_avg_first(id)) THEN
      timer_avg_first(id) = .FALSE.
      timer_average(id) = timer_time(id)
    ELSE
      timer_average(id) = avg_weight1 * timer_time(id) &
          + avg_weight2 * timer_average(id)
    END IF

  END SUBROUTINE timer_stop



  SUBROUTINE timer_reset

    INTEGER, PARAMETER :: id = c_timer_dt
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      first = .FALSE.
      timer_avg_first(c_timer_step) = .TRUE.
      RETURN
    END IF

    timer_time(id) = timer_time(c_timer_step) - timer_time(c_timer_io) &
        - timer_time(c_timer_balance)

    IF (timer_avg_first(id)) THEN
      timer_avg_first(id) = .FALSE.
      timer_average(id) = timer_time(id)
    ELSE
      timer_average(id) = avg_weight1 * timer_time(id) &
          + avg_weight2 * timer_average(id)
    END IF

    timer_first = 0.0_num
    timer_time = 0.0_num

  END SUBROUTINE timer_reset

#ifdef PERFMON
  SUBROUTINE timer_perf_init
    nperf_steps = 0
    IF (nsteps >= 0) THEN
      nperf_steps = nsteps / nperf_avgsteps
      IF ( MOD(nsteps, nperf_avgsteps) > 0 ) nperf_steps = nperf_steps + 1
    ELSE
      nperf_steps = NINT(t_end / dt)
    END IF
    ALLOCATE(&
         timer_perf_times(nperf_num,nperf_steps), &
         perf_avg(nperf_num, nperf_avgsteps))

    timer_perf_times  = 0.0_num
    perf_avg = 0.0_num
    
    IF (rank == 0) &
         PRINT '(A,I4.4,I5.4,A)', "initialising timer_perf arrays, size=[", &
         SHAPE(timer_perf_times),"]"
    nperf_curstep = 1
  END SUBROUTINE timer_perf_init

  SUBROUTINE timer_perf_tick(index)
    INTEGER i
    INTEGER, INTENT(IN) :: index
    LOGICAL tick

    tick = perfticks(index)
    IF (step == 0 .OR. perfpause) RETURN
    i = MOD(step, nperf_avgsteps)
    IF (i == 0) i = nperf_avgsteps
    perf_avg(index,i) = MPI_WTIME() - perf_avg(index,i)
    IF (rank==0 .AND. tick ) &
         PRINT "(A,I1,A,ES9.3)", "perf(index =",index,"), dt = ", perf_avg(index,i)
    IF (tick .AND. i == nperf_avgsteps) THEN !average
      timer_perf_times(index, nperf_curstep) = SUM(perf_avg(index,:)) / c_perf_avgsteps
      nperf_curstep = nperf_curstep + 1
      perf_avg(index,:) = 0.0_num
    END IF

    perfticks(index) = .NOT. tick
  END SUBROUTINE timer_perf_tick

  SUBROUTINE timer_perf_start
    perfpause = .FALSE.
  END SUBROUTINE timer_perf_start

  SUBROUTINE timer_perf_stop
    perfpause = .TRUE.
  END SUBROUTINE timer_perf_stop

  SUBROUTINE timer_perf_dump
    CHARACTER(len=512) fname
    ! 9999 ranks ought to be enough for anyone
    WRITE (fname, '(a, I4.4, a)') 'timerperf', rank, '.dat'
    
    OPEN(1337, file=TRIM(fname), access='stream')
    WRITE(1337) timer_perf_times
    CLOSE(1337)

    IF (rank /= 0) RETURN
    OPEN(1338, file=TRIM('timerperfhead.dat'), access='stream')
    WRITE(1338) nperf_num
    WRITE(1338) nperf_steps
    CLOSE(1338)
    
  END SUBROUTINE timer_perf_dump
#endif
END MODULE timer
