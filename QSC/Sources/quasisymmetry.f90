! Main program  10/12/18,3/1/19.wrapper for QSC().

program quasisymmetry

  use quasisymmetry_variables, only: start_time, total_time, general_option, general_option_single, general_option_scan, &
       N_procs, mpi_rank, proc0

  implicit none

  include 'mpif.h'

  call QSC()  !10/12/18,3/1/19.(6b)bulk of orig quasisymmetry.f90.

end program quasisymmetry
