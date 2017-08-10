subroutine print_hi(X) bind(C)
    implicit real*8(a-h,o-z)

      write(*,*) "Hello from Fortran."
      write(*,*) "X: ", X
end subroutine print_hi
