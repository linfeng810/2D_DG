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
    open(unit=1, file='output.txt', status='replace')
    do i = 1,size(meshele)
        write(1,*) i, '--', meshele(i)
    enddo
    write(1,*)
    do i = 1,size(meshface)
        write(1,*) i, '--',  meshface(i)
    enddo
    write(1,*)
    do i = 1,size(meshvertex)
        write(1,*) i, '--', meshvertex(i)
    enddo
    close(1)
endprogram