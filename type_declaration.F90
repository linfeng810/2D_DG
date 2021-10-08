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
    end type face
    ! elements
    type :: element
        integer, dimension(3) :: face 
        integer, dimension(3) :: node
    end type element 
end module 