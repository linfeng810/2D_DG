! this is to give basis function info
! 1. basis function at reference element
! 2. basis function at global element, its determinant, its derivatives

module basis_function
  use mesh_type

  implicit none 

  type :: shape_func
    real, dimension(3,1)::fe_funs=0.   ! assuming linear element, 3 local nodes, 1 Gaussian point
    real, dimension(1)::wei=0. ! Gaussian point weight
    real, dimension(3,3,1)::sfe_funs=0. ! 3 shape functions evaluated on 3 faces each at 1 Gauss point
    real, dimension(1)::swei=0. ! surface gaussian poitn weight
  end type

  type :: shape_func_dev
    real, dimension(2,3,1)::dev_funs=0. ! assuming linear element, 2 derivatives, 3 local nodes, 1 
      ! gaussian point
    real, dimension(1)::detwei=0. ! determinant (of Jacobian) * Gaussian point weight
    real, dimension(2,3,3,1)::sdev_funs=0. ! 2 derivatives of 3 local nodes on 3 edges evaluated @ 
      ! 1 gaussian point
    real, dimension(1)::sdetwei=0.  ! determinant (edge length) * Gaussian weight
  end type

  contains 
  subroutine init_shape_func(sf)
    type(shape_func) , intent(out):: sf 
    ! assuming 1st order linear element
    ! 1 Gaussian point
    ! gaussian point barycentric coordinate (1/3, 1/3, 1/3)
    ! shape function in terms of barycentric coorcinate 
    ! N1 = \zeta, N2 = \xi, N3 =\eta
    ! where \zeta = 1 - \xi - \eta
    ! therefore, shape function evaluated at gaussian point is:
    ! N1 = 1/3, N2 = 1/3, N3 = 1/3

    sf%fe_funs(:,1) = (/1./3.,1./3.,1./3./)
    sf%wei = (/1./)

    ! N1, N2, N3 evaluated at 3 edges
    ! (local numbering)
    !  3
    !  | \
    !  |  \
    !  3   2
    !  |    \
    !  1--1--2
    sf%sfe_funs(:,:,1) = reshape( (/0.5, 0., 0.5, 0.5, 0.5, 0., 0., 0.5, 0.5/) , (/3,3/))
    sf%wei = (/1./)

  end subroutine init_shape_func

  subroutine calc_local_shape_func(lsf, sf, lelement, meshvertex)
    type(shape_func_dev), intent(out):: lsf   ! local elements' shape function's derivatives
    type(shape_func), intent(in):: sf   ! a reference element's shape functions
    type(element), intent(in):: lelement  ! a local element
    ! type(face), dimension(:), intent(in):: meshface ! all mesh faces => seems useless
    type(vertex), dimension(:), intent(in):: meshvertex ! all mesh vertices
    ! real, dimension(2,3):: ref_coor=(/0.,0.,1.,0.,0.,1./) => seems this is useless.
    real, dimension(2,3):: x_all=0.
    real, dimension(2,2):: jacob, jacobit
    real :: a2
    integer :: inod

    ! print*, 'is input alright?', 'element' , lelement, ' mesh veretex' , meshvertex

    ! first we fetech Nodes coordinate in local element
    do inod = 1,3
      ! print*, 'lelement', lelement
      x_all(:,inod) = meshvertex( lelement%node(inod) )%coor(:)
      ! print*, 'x_all', x_all
    enddo

    ! assuming linear element
    ! physical coordinate (x,y), reference coordinate (\xi, \eta)
    ! x = x1 + (x2 - x1) \xi + (x3 - x1) \eta
    ! y = y1 + (y2 - y1) \xi + (y3 - y1) \eta   (1)
    ! we can easily calculate derivatives in reference element:
    ! dN1/d\xi = -1, dN1/d\eta = -1
    ! dN2/d\xi = 1,  dN2/d\eta = 0
    ! dN3/d\xi = 0,  dN3/d\eta = 1
    ! transform to local element, we would need to evaluate d(\xi, \eta)/d(x, y)
    ! d\xi/dx = (y3-y1)/2A,  d\xi/dy = (x1-x3)/2A
    ! d\eta/dx= (y1-y2)/2A,  d\eta/dy= (x2-x1)/2A
    ! (this is not obvious but can be derived from (1))
    ! then we can get derivatives in local element
    ! dN1/dx = dN1/d\xi * d\xi/dx + dN1/d\eta * d\eta/dx = (y2-y3)/2A
    ! dN1/dy = (x3-x2)/2A
    ! dN2/dx = (y3-y1)/2A
    ! dN2/dy = (x1-x3)/2A
    ! dN3/dx = (y1-y2)/2A
    ! dN3/dy = (x2-x1)/2A


    ! otherwise, we can use Gaussian quadrature

    ! 2 x triangle area
    a2  = (x_all(1,2) - x_all(1,1)) * (x_all(2,3) - x_all(2,1)) &
      - (x_all(1,3) - x_all(1,1)) * (x_all(2,2) - x_all(2,1)) ! = 2A
    ! det * weight (these two alway appears together)
    lsf%detwei(:) = sf%wei(:) * abs(a2)
    ! print*, 'a2', a2
    ! write(20, *) 'x', x_all

    ! calculate J^-T
    ! J
    jacob(1,1) = x_all(1,2) - x_all(1,1)
    ! jacob(1,2) = x_all(1,3) - x_all(1,1)
    ! jacob(2,1) = x_all(2,2) - x_all(2,1)
    jacob(2,2) = x_all(2,3) - x_all(2,1)
    ! transpose
    jacob(1,2) = x_all(2,2) - x_all(2,1)
    jacob(2,1) = x_all(1,3) - x_all(1,1)
    ! inverse
    jacobit(1,1) = jacob(2,2)/a2 
    jacobit(1,2) = -jacob(1,2)/a2 
    jacobit(2,1) = -jacob(2,1)/a2 
    jacobit(2,2) = jacob(1,1)/a2
    ! print*, 'jacob', jacob
    ! print*, 'jacobit', jacobit 

    ! calculate J^-T * \nabla N (->on reference element)
    ! \nabla N (->on reference element), for linear element, 
    ! derivatives on reference element is constant
    lsf%dev_funs(:,1,1) = (/-1.,-1./)
    lsf%dev_funs(:,2,1) = (/1.,0./)
    lsf%dev_funs(:,3,1) = (/0.,1./)
    ! print*, lsf%dev_funs(:,1,1)
    ! left multiply J^-T to transform to local element
    lsf%dev_funs(:,:,1) = matmul( jacobit , lsf%dev_funs(:,:,1) )
    ! print*, lsf%dev_funs
  end subroutine calc_local_shape_func

  subroutine calc_local_surf_sf(lsf, sf, lface, meshvertex)
    ! this is to calculate shape function value on a surface
    type(shape_func_dev), intent(inout)::lsf 
    type(shape_func), intent(in)::sf
    type(face), intent(in):: lface 
    type(vertex), dimension(:), intent(in):: meshvertex
    real, dimension(2,2):: x_all 
    integer :: inod
    
    ! print*, 'input inside face', lface
    ! first retrieve coordinates
    x_all(:,1) = meshvertex( (lface%vertex(1)) )%coor(:)
    x_all(:,2) = meshvertex( (lface%vertex(2)) )%coor(:)
    ! print*, x_all, 'sqrt', sqrt( (x_all(1,1)-x_all(1,2))**2 + (x_all(2,1)-x_all(2,2))**2 )
    ! edge length * weight
    lsf%sdetwei(1) = sqrt( (x_all(1,1)-x_all(1,2))**2 + (x_all(2,1)-x_all(2,2))**2 ) * 1.
    ! print*, 'edge length',  lsf%sdetwei(1)
    ! 2 derivatives of 3 local nodes shape functions evaluated on 3 edges, 1 gaussian point for each
    do inod = 1,3
      lsf%sdev_funs(:,:,inod,:) = lsf%dev_funs(:,:,:)
    enddo
  end subroutine calc_local_surf_sf
    

end module