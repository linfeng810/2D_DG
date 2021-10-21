program main
  ! we don't need input file at this moment as this is not to serve as a universal 
  ! solver. 
  ! we will solve a problem with a known analytical solution.
  use mesh_type
  use gmsh_input
  use basis_function
  use assemb_matrix_engine
  use mod_solve

  implicit none
  
  ! variable declaration
  type(element), allocatable::meshele(:)
  type(face), allocatable::meshface(:)
  type(vertex), allocatable::meshvertex(:)
  character (len=255)::gmsh_filename='example/square_fine.msh'
  integer::outfileno ! output file unit.
  
  real, dimension(:,:), allocatable::bigm 
  real, dimension(:), allocatable::rhs , phi
  integer :: i,j
  
  open(newunit=outfileno, file='tests/output.out', status='replace')
  
  ! read in mesh and create connectivity
  call gmsh_read(gmsh_filename, meshele, meshface, meshvertex)
  nele = size(meshele)
  nnod = 3*nele
  allocate( bigm(nele*3, nele*3) , rhs(nele*3) , phi(nele*3) )
  bigm= 0.
  rhs = 0.
  phi = 0.

  call assemb_matrix(meshele, meshface, meshvertex, bigm, rhs)

  call lusolve(bigm,rhs,phi,nnod)
  
  print*, nele

  ! do i = 1, nele*3
  !   write(outfileno, '(999E24.10)') (bigm(i,j),j=1,nele*3)
  ! enddo
  ! write(outfileno, *)
  ! do i = 1,nele*3
  !   write(outfileno, *) rhs(i)
  ! enddo
  write(outfileno, *) 'nnod=', nnod
  do i = 1,nnod 
    write(outfileno, *) meshvertex(i)%coor, phi(i)
  enddo
  close(outfileno)

end program