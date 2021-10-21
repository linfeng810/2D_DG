! this is a direct LU decomposition solver

module mod_solve

  implicit none 
  
  contains
  subroutine lusolve(A,b,x,n)
    integer, intent(in)::n
    real, intent(in)::A(n,n),b(n)
    real, intent(out)::x(n)
    real :: l(n,n), u(n,n)
    integer :: i, j, k
    real :: y(n)
    
    ! do LU decomposition and stores in l & u respectively
    do j = 1,n
      u(1,j) = a(1,j)
    enddo
    do i = 2,n
      l(i,1) = a(i,1) / u(1,1)
    enddo
    do i = 2,n-1
      u(i,i) = a(i,i)
      do k = 1,i-1
        u(i,i) = u(i,i) - l(i,k) * u(k,i)
      enddo
      do j = i+1 , n
        u(i,j) = a(i,j)
        do k = 1,i-1
          u(i,j) = u(i,j) - l(i,k)*u(k,j)
        enddo
        l(j,i) = a(j,i)
        ! print*, 'lji,aji', l(j,i), a(j,i)
        do k = 1,i-1
          l(j,i) = l(j,i) - l(j,k) * u(k,i)
          ! print*, ' ljk, uki, lji' , l(j,k), u(k,i), l(j,i)
        enddo
        l(j,i) = l(j,i) / u(i,i)
        ! print*, 'j,i,lji', j,i, l(j,i)
      enddo
    enddo
    u(n,n) = a(n,n)
    do k = 1,n-1
      u(n,n) = u(n,n) - l(n,k) * u(k,n)
    enddo
    
    write(*,*) 'A is' 
    do i = 1,n
      write(*,*) (a(i,j), j=1,n)
    enddo
    write(*,*) 'u is' 
    do i = 1,n
      write(*,*) (u(i,j), j=1,n)
    enddo
    write(*,*) 'l is' 
    do i = 1,n
      write(*,*) (l(i,j), j=1,n)
    enddo
    
    ! solving Ax=b by Ux=y, Ly=b
    y(1) = b(1)
    do i = 2,n 
      y(i) = b(i)
      do j = 1,i-1 
        y(i) = y(i) - l(i,j) * y(j)
      enddo
    enddo
    x(n) = y(n) / u(n,n)
    do i = n-1,1,-1
      x(i) = y(i)
      do j = i+1,n
        x(i) = x(i) - u(i,j)*x(j)
      enddo
      x(i) = x(i) / u(i,i)
    enddo

  endsubroutine
endmodule
