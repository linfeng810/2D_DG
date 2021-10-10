! this is to test the module "basis_function.f90"

program main
    use mesh_type
    use basis_function

    implicit none 

    type(shape_func) :: ref_sf
    type(shape_func_dev) :: loc_sf 
    type(element) :: test_ele 
    type(vertex), dimension(3):: vertices
    
    test_ele%node(:) = (/1,2,3/)
    vertices(1)%coor(:) = (/1.5,2.4/)
    vertices(2)%coor(:) = (/4.5,2.6/)
    vertices(3)%coor(:) = (/1.9,3.5/)

    call init_shape_func(ref_sf)
    call calc_local_shape_func(loc_sf, ref_sf, test_ele, vertices)

    open(unit=10, file='tests/output.out', status='replace')
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
    

end program