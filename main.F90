program main
    ! we don't need input file at this moment as this is not to serve as a universal 
    ! solver. 
    ! we will solve a problem with a known analytical solution.
    use gmsh_io
    implicit none

    ! derived data types
    ! vertices of elements (global nodes)
    type :: vertex
        real, dimension(2) :: coor
    end type vertex
    ! faces of elements (global edges)
    type :: face
        integer, dimension(2) :: vertex ! index of vertices
        integer, dimension(2) :: neighbor ! global number of elements sharing the face
        integer :: bctype ! boundary type, 0 for interior, 1 Dirichlet, 2 Neumann
    end type face
    
    ! variable declaration
    real :: Lx(2), Ly(2), hx, hy 
    integer :: nnod, inod, jnod, itotal, jtotal, nedge
    type(vertex), dimension(:), allocatable :: meshvertex
    type(face), dimension(:), allocatable :: meshface
    integer :: node_num, node_dim, element_num, element_order 
    character :: gmsh_filename='example/square.msh'

    ! domain and subdivision (triangulation)
    Lx = (/0.,1./)    ! x range
    Ly = (/0.,1./)    ! y range
    hx = 0.1    ! x size
    hy = 0.1    ! y size
    itotal = int((Lx(2)-Lx(1)) / hx + 1)  ! x-direction nodes number
    jtotal = int((Ly(2)-Ly(1)) / hy + 1)  ! y-direction nodes number
    ! vertex coordinate
    nnod = itotal * jtotal    ! total number of nodes
    allocate(meshvertex(nnod))
    do inod = 1,itotal 
        do jnod = 1,jtotal 
            meshvertex( (jnod-1)*jtotal + inod )%coor(1) = (inod-1)*hx + Lx(1)
            meshvertex( (jnod-1)*jtotal + inod )%coor(2) = (jnod-1)*hy + Ly(1)
        enddo
    enddo
    ! element faces (edges)
    
    call gmsh_size_read(gmsh_filename,node_num, node_dim, element_num, element_order)
    ! call gmsh_data_read(gmsh_filename,2,)
    print*, node_num, node_dim, element_num, element_order

end program