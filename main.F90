program main
  ! we don't need input file at this moment as this is not to serve as a universal 
  ! solver. 
  ! we will solve a problem with a known analytical solution.
  use mesh_type
  use gmsh_input
  implicit none

  ! variable declaration
  type(element), allocatable::meshele(:)
  type(face), allocatable::meshface(:)
  type(vertex), allocatable::meshvertex(:)
  character (len=255)::gmsh_filename='example/square.msh'

  ! read in mesh and create connectivity
  call gmsh_read(gmsh_filename, meshele, meshface, meshvertex)
end program