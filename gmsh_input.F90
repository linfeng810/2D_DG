! read in gmsh file into our element types

module gmsh_input
  use mesh_type

  implicit none

  public :: gmsh_read 
  
  contains 
  
  subroutine gmsh_read(file_name, meshele, meshface, meshvertex)
    character (len=*), intent(in)::file_name  ! gmsh file name
    type(element), allocatable, intent(inout)::meshele(:) ! to store elements
    type(face), allocatable, intent(inout)::meshface(:) ! to store element faces, number = 3 * element number, including interior and boundary faces
    type(face), allocatable::bfaces(:) ! this is to store faces on boundaries
      ! actual mesh faces will have interior faces as well.
    type(vertex), allocatable, intent(inout)::meshvertex(:) ! to store element nodes, number = 3 * element number
    type(vertex), allocatable::vertices(:)  ! this only stores one copy 
      ! of coordinates at same position. (but there can be multiple nodes)
    integer :: ierr , i,  idx, ele, ele2, nface, inod, iface
    integer :: edge_nodes(2)
    integer :: nelements=0, npoints=0, nbcface=0
    real, dimension(3) :: coor
    integer, dimension(8):: element_info
    character(len=255) :: text, t, word
    integer :: mesh_file_unit
    
    open(newunit=mesh_file_unit, file=file_name, status='old', action='read') 

    ! first read in gmsh file raw data: 
    ! 1. nodes coordinate
    ! 2. boundary edge and nodes associated with it 
    !    (in the future we can use physical id to indentify different bc types)
    ! 3. element and nodes associated with it
    ! then we enrich mesh info:
    ! 4. create local nodes and faces, find neighbouring elements, faces, and nodes
    ! 5. reordering nodes numbering
    ! Step 1: ===================================================================
    do  
      ! this find string snippet is from:
      ! https://stackoverflow.com/questions/29125581/fortran-find-string-in-txt-file
      ! by Fortranner
      read (mesh_file_unit, "(a)", iostat=ierr) text
      if (ierr /= 0) exit 
      read (text,*) word 
      if (word == '$Nodes') then ! start reading nodes information
        read (mesh_file_unit, *) nnod ! total number of nodes
        ! print*, 'nnod=', nnod
        allocate(vertices(nnod))
        do 
          read (mesh_file_unit, *, iostat=ierr) idx, coor(1), coor(2), coor(3)
          if (ierr /= 0) exit ! if read in '$EndNodes', will exit this loop
          vertices(idx)%coor(1) = coor(1)
          vertices(idx)%coor(2) = coor(2)
          ! print*, 'idx', idx, 'coor', coor
        enddo
      endif
      if (word == '$Elements') then ! start reading element information 
        read (mesh_file_unit, *) nelements  ! total number of elements (1D, 2D, 3D)
        ! print *, 'nelements=', nelements
        do  
          read (mesh_file_unit, "(A)", iostat=ierr) text
          ! print *, text 
          if (ierr /= 0) exit 
          if (text(1:12)=="$EndElements") exit
          i=0
          ! print *, text
          ! parsing element information in one line
          !   according to gmsh document, this line may contain various (6-8) numbers
          do while (text /= "")
            i=i+1
            call GetTag(t,text)
            read(t,*) element_info(i)
          enddo
          select case (i)
          case (6) ! 0D elements (points)
            npoints = npoints + 1 ! useless
          case (7) ! 1D elements (edges)
            nbcface = nbcface + 1 ! number of boundary faces
          case (8) ! 2D elements (faces)
            nele = nele + 1 ! number of 2D elements (actual *DG* elements)
          case default 
            print *, 'something is wrong when parsing Elements data, what is being parsed is:'
            print *, text
          end select 
          ! print*, 'npoints=',npoints, 'nbcface=', nbcface, 'nele=', nele
        enddo 
      endif
    enddo

    allocate(bfaces(nbcface), meshele(nele)) 

    ! Step 2&3: ==========================================================================
    rewind (mesh_file_unit) ! go back and read again to store node info associated with 
          ! points, edges, and faces.
    do  
      read (mesh_file_unit, "(a)", iostat=ierr) text
      if (ierr /= 0) exit 
      read (text,*) word 
      if (word == '$Elements') then ! start reading element information
        ! jump 0D elements
        do i=1,npoints+1
          read(mesh_file_unit,*)
        enddo
        ! 1D elements
        do nface=1,nbcface 
          read (mesh_file_unit, "(A)", iostat=ierr) text 
          i=0
          do while (text /= "")
            i=i+1
            call GetTag(t,text)
            read(t,*) element_info(i)
          enddo
          bfaces(nface)%vertex(1) = element_info(6)
          bfaces(nface)%vertex(2) = element_info(7)
          bfaces(nface)%bctype = 1
          ! print *, 'bfaces=', bfaces(nface)%vertex
        enddo
        ! 2D elements
        do ele=1,nele
          read (mesh_file_unit, "(A)", iostat=ierr) text 
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
          ! print *, 'ele', ele, 'meshele=', meshele(ele)%node
        enddo
      endif
    enddo
    close(mesh_file_unit) ! gmsh file no longer has useful info for us.

    ! Step 4: =============================================================================
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
        ! meshele(ele)%node(inod) = (ele-1)*3+inod ! change element node numbering to DG global numbering
          ! can't do it here because we need this information to find neighbor element!
        ! print*, 'coordinate for global node', (ele-1)*3+inod, 'is' ,meshvertex( (ele-1)*3+inod )%coor
      enddo

      do iface=1,3
        ! global edges' nodes 
          ! (local numbering)
          !  3
          !  | \
          !  |  \
          !  3   2
          !  |    \
          !  1--1--2
        if (iface.eq.3) edge_nodes=(/3,1/)
        if (iface.ne.3) edge_nodes=(/iface, iface+1/)
        meshface( (ele-1)*3+iface )%vertex(1) = (ele-1)*3+edge_nodes(1)
        meshface( (ele-1)*3+iface )%vertex(2) = (ele-1)*3+edge_nodes(2)
        meshele( ele )%face(iface) = (ele-1)*3+iface
        meshface( (ele-1)*3+iface )%lc_vertex = edge_nodes

        ! is this meshface a boundary face or an interior face?
        meshface( (ele-1)*3+iface )%bctype=0 ! default: interior face
        do i = 1,nbcface
          if ( any( meshele(ele)%node(edge_nodes(1)) .eq. bfaces(i)%vertex ) .and. &
            any( meshele(ele)%node(edge_nodes(2)) .eq. bfaces(i)%vertex)) then
            meshface( (ele-1)*3+iface )%bctype=1 ! is boundary, assuming Dirichlet
            exit
          endif
        enddo

        ! now we find and store neighbor elements to this meshface
        meshface( (ele-1)*3+iface )%neighbor(1) = ele ! first, this face belongs to 
          ! its own element. Then we need to find its other neighbor
        if ( meshface( (ele-1)*3+iface )%bctype.eq.0) then
          do ele2 = 1,nele
            if (any(meshele(ele)%node(edge_nodes(1)) .eq. meshele(ele2)%node) .and. &
              any(meshele(ele)%node(edge_nodes(2)) .eq. meshele(ele2)%node) .and. & 
              ele .ne. ele2) then
              ! found neighbouring element
                meshface( (ele-1)*3+iface )%neighbor(2) = ele2
              ! find neighbouring node to first edge node
              do inod = 1,nloc 
                ! print*,  meshele(ele)%node(edge_nodes(1)), meshele(ele2)%node(inod)
                if ( meshele(ele)%node(edge_nodes(1)) .eq. meshele(ele2)%node(inod) ) then
                  ! print*, 'I got in!'
                  meshface( (ele-1)*3+iface )%nb_node(1) = (ele2-1)*3+inod 
                  exit 
                endif
              enddo
              ! find neighbouring node to 2nd edge node
              do inod = 1,nloc 
                if ( meshele(ele)%node(edge_nodes(2)) .eq. meshele(ele2)%node(inod) ) then 
                  meshface( (ele-1)*3+iface )%nb_node(2) = (ele2-1)*3+inod 
                  exit 
                endif
              enddo
              
              ! now find local index of this neighbouring face
              ! mod( nb_node(1) + nb_node(2) , 3) = 0 ..... local node 1&2 local face 1
              !                                   = 1 ..... local node 1&3 local face 3
              !                                   = 2 ..... local node 2&3 local face 2
              select case ( mod ( sum(meshface( (ele-1)*3+iface )%nb_node) , 3 ) )
              case(0)
                meshface( (ele-1)*3 + iface )%nb_iface=1
              case(1)
                meshface( (ele-1)*3 + iface )%nb_iface=3
              case(2)
                meshface( (ele-1)*3 + iface )%nb_iface=2
              case default 
                print*, 'in sbrt gmsh_read: something is wrong ', & 
                  'when finding local index of neighbouring face.'
              end select

              ! print*, 'ele', ele, 'ele2', ele2, 'node', meshface((ele-1)*3+iface)%vertex, &
              !   'nbnode', meshface( (ele-1)*3 + iface)%nb_node
              ! print*, 'ele', ele, 'ele2', ele2, 'iface', iface, 'iface2', meshface((ele-1)*3+iface)%nb_iface

              exit ! do ele2

            endif ! ele.node == ele2.node ??
          enddo ! do ele2
        endif ! bctype .eq. 0 ??
        ! print *, 'edge neighbor elements are', meshface( (ele-1)*3+inod )%neighbor

      enddo
    enddo

    ! Step 5: ===========================================================================
    ! after all, change element node numbering to DG global numbering
    do ele=1,nele 
      do inod = 1,3
        meshele(ele)%node(inod) = (ele-1)*3+inod 
      enddo
    enddo
    
    contains
    subroutine GetTag(t, s)   ! remove first token, t, from s
      ! from http://computer-programming-forum.com/49-fortran/acb879484e8189d2.htm
      ! author: Jos Bergervo
      character(len=*), intent(in out) :: s
      character(len=*), intent(out)  :: t
      integer  :: ii
      s = adjustl(s)      ! Now the 1st token is directly at the start
      ii = index(s," ")
      t = s(1:ii)        ! returning s(1:0) gives just empty string
      s = s(ii:len(s))     ! Keep only the rest of string s
    end subroutine GetTag
  end subroutine gmsh_read
    

end module