! this is to test our LU solver

program main
  use mod_solve
  implicit none
  
  real, dimension(4,4):: A
  real, dimension(4):: b, x
  integer :: n=4
  A = reshape( (/9.,18.,9.,-27., &
              18., 45., 0., -45., &
              9., 0., 126., 9., &
              -27., -45., 9., 135./) , (/4,4/) )
  b = (/1., 2., 16., 8./)
  
  call lusolve(A,b,x,n)
  write(*,*) x
  ! correct answer is : 1/9, 1/9, 1/9, 1/9
endprogram
