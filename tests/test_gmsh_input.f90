! this is to test if module gmsh_input works fine.
program main
  use mesh_type
  use gmsh_input

  implicit none
  
  type(element), allocatable::meshele(:)
  type(face), allocatable::meshface(:)
  type(vertex), allocatable::meshvertex(:)
  character (len=255)::gmsh_filename='example/square.msh'
  integer :: i

  call gmsh_read(gmsh_filename, meshele, meshface, meshvertex)
  open(unit=20, file='tests/output.out', status='replace')
  write(20,'(999A12)') 'ele', 'face1', 'face2', 'face3', &
    'node1', 'node2', 'node3'
  do i = 1,size(meshele)
    write(20,'(999I12)') i, meshele(i)
  enddo
  write(20,'(999A12)') 'face', 'node1', 'node2', &
    'nb_ele_1', 'nb_ele_2', 'bctype' , 'nb_node_1', 'nb_node_2', &
    'lc_node_1', 'lc_node_2', 'nb_iface'
  do i = 1,size(meshface)
    write(20,'(999I12)') i, meshface(i)
  enddo
  write(20,'(999A12)') 'vertex', 'x', 'y'
  do i = 1,size(meshvertex)
    write(20,'(I12,2F12.5)') i, meshvertex(i)
  enddo
  close(20)
  print*, 'complete'
endprogram