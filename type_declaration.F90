! derived data types
module mesh_type
    implicit none

    ! vertices of elements (global nodes)
    type :: vertex
        real, dimension(2) :: coor
    end type vertex
    ! faces of elements (global edges)
    type :: face
        integer, dimension(2) :: vertex ! index of vertices
        integer, dimension(2) :: neighbor ! global number of elements sharing the face
        integer :: bctype ! boundary type, 0 for interior, 1 Dirichlet, 2 Neumann
        integer, dimension(2) :: nb_node ! global number of neighbouring nodes - corresponding to face%vertex(2).
        integer, dimension(2) :: lc_vertex ! local index of vertices
        integer :: nb_iface ! neighbouring face local index
    end type face
    ! elements
    type :: element
        integer, dimension(3) :: face 
        integer, dimension(3) :: node
    end type element 

    save
    ! dimensions
    integer,parameter :: ngi =3   ! volume quadrature points number
    integer,parameter :: nsgi=2   ! surface quadrature points number
    integer,parameter :: ndim=2   ! dimension of the problem
    integer,parameter :: nloc=3   ! local node number 
    integer :: nnod     ! total node number 
    integer :: nele     ! total element number

end module 