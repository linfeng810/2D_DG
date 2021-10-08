! this is to test if module gmsh_input works fine.
program main
    use mesh_type
    use gmsh_input

    implicit none
    
    type(element), allocatable::meshele(:)
    type(face), allocatable::meshface(:)
    type(vertex), allocatable::meshvertex(:)
    character (len=255)::gmsh_filename='example/square.msh'

    call gmsh_read(gmsh_filename, meshele, meshface, meshvertex)
    open(unit=1, file='output.txt', status='replace')
    write(1,*) meshele 
    write(1,*)
    write(1,*) meshface 
    write(1,*)
    write(1,*) meshvertex 
    close(1)
endprogram