! this is to test the module "basis_function.f90"

program main
    use mesh_type
    use basis_function
    use gmsh_input

    implicit none 

    type(shape_func) :: ref_sf
    type(shape_func_dev) :: loc_sf 
    type(element) :: test_ele 
    type(vertex), dimension(3):: vertices

    type(element), allocatable::meshele(:)
    type(face), allocatable::meshface(:)
    type(vertex), allocatable::meshvertex(:)
    character (len=255)::gmsh_filename='example/square.msh'

    integer::ele, iface
    
    test_ele%node(:) = (/1,2,3/)
    vertices(1)%coor(:) = (/1.5,2.4/)
    vertices(2)%coor(:) = (/4.5,2.6/)
    vertices(3)%coor(:) = (/1.9,3.5/)

    call init_shape_func(ref_sf)
    call calc_local_shape_func(loc_sf, ref_sf, test_ele, vertices)

    open(unit=10, file='tests/output.out', status='replace')
    write(10,*) ' ======================= trivial test =================='
    write(10,*) '  basis functions at reference element:'
    write(10,*) '  (values on Gaussian quadrature point)'
    write(10,*) '  N1,              N2,              N3'
    write(10,*) ref_sf%fe_funs
    write(10,*)  
    write(10,*) '  derivatives of basis functions at local element:'
    write(10,*) '  (still, values on Gaussian quadrature point)'
    write(10,*) '   dN1/dx,         dN2/dx,          dN3/dx'
    write(10,*) loc_sf%dev_funs(1,:,1)
    write(10,*) '   dN1/dy,         dN2/dy,          dN3/dy'
    write(10,*) loc_sf%dev_funs(2,:,1)
    write(10,*) '  determinant x weight'
    write(10,*) loc_sf%detwei
    
    write(10,*) ' ======================= test on square.msh ==========='
    call gmsh_read(gmsh_filename, meshele, meshface, meshvertex)
    nele = size(meshele)
    write(10, *) 'nele=', nele
    write(10, '(2A12,999A26)') 'ele', 'iface', 'dev_funs111', 'dev_funs121', 'dev_funs131', &
        'dev_funs211', 'dev_funs221', 'dev_funs231', 'detwei', &
        'sdx111', 'sdy111', 'sdx211', 'sdy211', 'sdx311', 'sdy311',  &
        'sdx121', 'sdy121', 'sdx221', 'sdy221', 'sdx321', 'sdy321', &   ! surface derivative x/y of nodal shape function 1/2/3 on face 1/2/3
        '...', '...', '...', '...', '...', '...', &
        'sdetwei'
    do ele = 1,nele
        print*, 'input ele', ele
        call calc_local_shape_func(loc_sf, ref_sf, meshele(ele), meshvertex)
        do iface = 1,3
            ! print*, 'input outside face', meshface( meshele(ele)%face(iface) )
            call calc_local_surf_sf(loc_sf, ref_sf, meshface( meshele(ele)%face(iface) ), meshvertex)
            write(10,'(2I12,999E26.10)') ele, iface, loc_sf 
        enddo
    enddo
    write(10, *) 'end'
    close(10)
    
end program