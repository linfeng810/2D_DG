program main
    ! we don't need input file at this moment as this is not to serve as a universal 
    ! solver. 
    ! we will solve a problem with a known analytical solution.

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
    ! elements
    type :: element
        integer, dimension(3) :: face 
        integer, dimension(3) :: node
    end type element 
    
    ! variable declaration
    real :: Lx(2), Ly(2), hx, hy 
    integer :: nnod, inod, jnod, itotal, jtotal, nedge, i
    type(vertex), dimension(:), allocatable :: meshvertex
    type(face), dimension(:), allocatable :: meshface
    type(element), dimension(:), allocatable :: meshele
    integer :: node_num, node_dim, element_num, element_order 
    character (len=255) :: gmsh_filename='example/square.msh'
    real, allocatable:: node_x(:,:)
    integer, allocatable:: element_node(:,:)

    ! ! domain and subdivision (triangulation)
    ! Lx = (/0.,1./)    ! x range
    ! Ly = (/0.,1./)    ! y range
    ! hx = 0.1    ! x size
    ! hy = 0.1    ! y size
    ! itotal = int((Lx(2)-Lx(1)) / hx + 1)  ! x-direction nodes number
    ! jtotal = int((Ly(2)-Ly(1)) / hy + 1)  ! y-direction nodes number
    ! ! vertex coordinate
    ! nnod = itotal * jtotal    ! total number of nodes
    ! allocate(meshvertex(nnod))
    ! do inod = 1,itotal 
    !     do jnod = 1,jtotal 
    !         meshvertex( (jnod-1)*jtotal + inod )%coor(1) = (inod-1)*hx + Lx(1)
    !         meshvertex( (jnod-1)*jtotal + inod )%coor(2) = (jnod-1)*hy + Ly(1)
    !     enddo
    ! enddo
    ! ! element faces (edges)
    
    call gmsh_size_read(gmsh_filename,node_num, node_dim, element_num, element_order)
    
    allocate ( node_x(1:node_dim,1:node_num) )
    allocate ( element_node(1:element_order,1:element_num) )
    call gmsh_data_read (gmsh_filename, node_dim, node_num, node_x, &
    element_order, element_num, element_node)
    
    allocate ( meshvertex(node_num) )
    allocate ( meshele(element_num) )
    do i = 1, node_num
        meshvertex%coor(i) = node_x(:,i)
    enddo
    do i = 1, element_num 
        print *, element_node(:,i), i
    enddo

end program