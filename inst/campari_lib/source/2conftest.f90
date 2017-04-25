program blah
  implicit none
  integer tpi, OMP_GET_THREAD_NUM
  call OMP_SET_NUM_THREADS(10)
!$OMP parallel private(tpi)
  tpi = OMP_GET_THREAD_NUM()
!$OMP CRITICAL(PRINT)
  write(*,*) tpi(1000)
!$OMP END CRITICAL(PRINT)
!$OMP END parallel
end
