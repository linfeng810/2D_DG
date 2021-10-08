! read in gmsh file into our element types

module gmsh_input
    use mesh_type

    implicit none

    public :: gmsh_read 
    
    contains 
    
    subroutine gmsh_read(file_name, meshele, meshface, meshvertex)
        character (len=*), intent(in)::file_name
        type(element), allocatable, intent(inout)::meshele(:)
        type(face), allocatable, intent(inout)::meshface(:)
        type(face), allocatable::bfaces(:) ! this is to store faces on boundaries
            ! actual mesh faces will have interior faces as well.
        type(vertex), allocatable, intent(inout)::meshvertex(:)
        type(vertex), allocatable::vertices(:)  ! this only stores one copy 
            ! of coordinates at same position. (but there can be multiple nodes)
        integer :: ierr , i, nnod, idx, ele, ele2, nface, inod
        integer :: edge_nodes(2)
        integer :: nelements=0, npoints=0, nbcface=0, nele=0
        real, dimension(3) :: coor
        integer, dimension(8):: element_info
        character(len=255) :: text, t, word
        
        open(unit=10, file=file_name, status='old', action='read') 
        do  
            ! this find string snippet is from:
            ! https://stackoverflow.com/questions/29125581/fortran-find-string-in-txt-file
            ! by Fortranner
            read (10, "(a)", iostat=ierr) text
            if (ierr /= 0) exit 
            read (text,*) word 
            if (word == '$Nodes') then ! start reading nodes information
                read (10, *) nnod
                ! print*, 'nnod=', nnod
                allocate(vertices(nnod))
                do 
                    read (10, *, iostat=ierr) idx, coor(1), coor(2), coor(3)
                    if (ierr /= 0) exit 
                    vertices(idx)%coor(1) = coor(1)
                    vertices(idx)%coor(2) = coor(2)
                    ! print*, 'coor', coor
                enddo
            endif
            if (word == '$Elements') then ! start reading element information 
                read (10, *) nelements
                ! print *, 'nelements=', nelements
                do  
                    read (10, "(A)", iostat=ierr) text
                    ! print *, text 
                    if (ierr /= 0) exit 
                    if (text(1:12)=="$EndElements") exit
                    i=0
                    ! print *, text
                    do while (text /= "")
                        i=i+1
                        call GetTag(t,text)
                        read(t,*) element_info(i)
                    enddo
                    select case (i)
                    case (6) ! 0D elements (points)
                        npoints = npoints + 1 
                    case (7) ! 1D elements (edges)
                        nbcface = nbcface + 1
                    case (8) ! 2D elements (faces)
                        nele = nele + 1
                    case default 
                        ! print *, 'something is wrong when parsing Elements data, what is being parsed is:'
                        ! print *, text
                    end select 
                    ! print*, 'npoints=',npoints, 'nbcface=', nbcface, 'nele=', nele
                enddo 
            endif
        enddo
        allocate(bfaces(nbcface), meshele(nele)) 
        rewind (10) ! go back and read again to store node info associated with 
                    ! points, edges, and faces.
        do  
            read (10, "(a)", iostat=ierr) text
            if (ierr /= 0) exit 
            read (text,*) word 
            if (word == '$Elements') then ! start reading element information
                ! jump 0D elements
                do i=1,npoints+1
                    read(10,*)
                enddo
                ! 1D elements
                do nface=1,nbcface 
                    read (10, "(A)", iostat=ierr) text 
                    i=0
                    do while (text /= "")
                        i=i+1
                        call GetTag(t,text)
                        read(t,*) element_info(i)
                    enddo
                    bfaces(nface)%vertex(1) = element_info(6)
                    bfaces(nface)%vertex(2) = element_info(7)
                    bfaces(nface)%bctype = 1
                    ! print *, 'bfaces=', bfaces(nface)
                enddo
                ! 2D elements
                do ele=1,nele
                    read (10, "(A)", iostat=ierr) text 
                    i=0
                    do while (text /= "")
                        i=i+1
                        call GetTag(t,text)
                        read(t,*) element_info(i)
                        ! print *, 'text=',text,'t=', t, 'element_info', element_info(i)
                    enddo
                    meshele(ele)%node(1) = element_info(6)
                    meshele(ele)%node(2) = element_info(7)
                    meshele(ele)%node(3) = element_info(8)
                    ! print *, 'meshele=', meshele(ele)
                enddo
            endif
        enddo
        close(10) ! gmsh file no longer has useful info for us.

        ! now we are going to build edges and ... vertices actually. (DG)
        allocate(meshvertex(nele*3), meshface(nele*3)) ! every 2D element 
            ! has 3 vertices and 3 faces
        do ele = 1,nele 
            ! brutely build meshvertex and mesh face
            do inod = 1,3
                ! global nodes
                meshvertex( (ele-1)*3+inod )%coor(:) = &
                    vertices( meshele(ele)%node(inod) )%coor(:)
                    ! this gives us a natural local to global
                    ! map. I guess IC-FERST uses the same way?
                print*, 'coordinate for global node', (ele-1)*3+inod, 'is' ,meshvertex( (ele-1)*3+inod )%coor
                
                ! global edges' nodes 
                    ! (local numbering)
                    !  3
                    !  | \
                    !  |  \
                    !  3   2
                    !  |    \
                    !  1--1--2
                if (inod.eq.3) edge_nodes=(/3,1/)
                if (inod.ne.3) edge_nodes=(/inod, inod+1/)
                meshface( (ele-1)*3+inod )%vertex(1) = (ele-1)*3+edge_nodes(1)
                meshface( (ele-1)*3+inod )%vertex(2) = (ele-1)*3+edge_nodes(2)
                print*, 'global nodes for edge', (ele-1)*3+inod, 'is ', meshface( (ele-1)*3+inod )%vertex

                ! is this meshface a boundary face or an interior face?
                meshface( (ele-1)*3+inod )%bctype=0 ! default: interior face
                do i = 1,nbcface
                    if ( any( meshele(ele)%node(inod) .eq. bfaces(i)%vertex ) .and. &
                        any( meshele(ele)%node(edge_nodes(2)) .eq. bfaces(i)%vertex)) then
                        meshface( (ele-1)*3+inod )%bctype=1 ! is boundary, assuming Dirichlet
                        exit
                    endif
                enddo
                print *, 'edge type (bc or not)', meshface( (ele-1)*3+inod )%bctype

                ! now we find and store neighbor elements to this meshface
                meshface( (ele-1)*3+inod )%neighbor(1) = ele ! first, this face belongs to 
                    ! its own element. Then we need to find its other neighbor
                if ( meshface( (ele-1)*3+inod )%bctype.eq.0) then
                    do ele2 = 1,nele
                        if (any(meshele(ele)%node(inod) .eq. meshele(ele2)%node) .and. &
                            any(meshele(ele)%node(edge_nodes(2)) .eq. meshele(ele2)%node)) then
                            meshface( (ele-1)*3+inod )%neighbor(1) = ele2
                            exit
                        endif
                    enddo
                endif
                print *, 'edge neighbor elements are', meshface( (ele-1)*3+inod )%neighbor

            enddo
        enddo
        
        contains
        subroutine GetTag(t, s)     ! remove first token, t, from s
            ! from http://computer-programming-forum.com/49-fortran/acb879484e8189d2.htm
            ! author: Jos Bergervo
            character(len=*), intent(in out) :: s
            character(len=*), intent(out)    :: t
            integer    :: ii
            s = adjustl(s)          ! Now the 1st token is directly at the start
            ii = index(s," ")
            t = s(1:ii)              ! returning s(1:0) gives just empty string
            s = s(ii:len(s))         ! Keep only the rest of string s
        end subroutine GetTag
    end subroutine gmsh_read
      

end module