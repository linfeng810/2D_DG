! this is to give basis function info
! 1. basis function at reference element
! 2. basis function at global element, its determinant, its derivatives

module basis_function
  use mesh_type

  implicit none 

  type :: shape_func
    real, dimension(nloc,ngi)::fe_funs=0.   ! assuming linear element, 3 local nodes, 3 Gaussian point
    real, dimension(ngi)::wei=0. ! Gaussian point weight
    real, dimension(nloc,3,nsgi)::sfe_funs=0. ! 3 shape functions evaluated on 3 faces each at 2 Gauss points
    real, dimension(nsgi)::swei=0. ! surface gaussian poitn weight
  end type

  type :: shape_func_dev
    real, dimension(ndim,nloc,ngi)::dev_funs=0. ! assuming linear element, 2 derivatives, 3 local nodes, 3 
      ! gaussian point
    real, dimension(ngi)::detwei=0. ! determinant (of Jacobian) * Gaussian point weight
    real, dimension(ndim,nloc,3,nsgi)::sdev_funs=0. ! 2 derivatives of 3 local nodes on 3 edges evaluated @ 
      ! 1 gaussian point
    real, dimension(nsgi)::sdetwei=0.  ! determinant (edge length) * Gaussian weight
  end type

  contains 
  subroutine init_shape_func(sf)
    type(shape_func) , intent(out):: sf 
    real :: a,b,c
    ! assuming 1st order linear element
    ! 3 Gaussian points (mid points of each edge)
    ! assuming GI poitns has same numbering as edge (shown below)
    !  3
    !  | \
    !  |  \
    !  3   2
    !  |    \
    !  1--1--2

    sf%fe_funs(1,:) = (/.5,0.,.5/); sf%fe_funs(2,:) = (/.5,.5,0./);
    sf%fe_funs(3,:) = (/0.,.5,.5/)
    sf%wei = 1./3.

    !======================SUSPECT==============================
    !=======DO WE REALLY NEED THIS? ++++++++++++++++++++++++++++
    !=====Yes we do! Only that we need to change to 1D quadrature (surface/face)
    ! N1, N2, N3 evaluated at 3 edges
    ! (local numbering) (with edge direction, Gaussian points should be destributed along this direction)
    !  3
    !  | ┌
    !  |  \
    !  3   2
    !  v    \
    !  1--1->2
    ! note that this time we need 1D quadrature, that is 2-point: +- 1/sqrt(3), weight 1
    ! can take value either (1-1/sqrt(3))/2 or (1+1/sqrt(3))/2 or 0
    a = (1.-1./sqrt(3.))/2.; b = (1.+1./sqrt(3.))/2.; c = 0.
    sf%sfe_funs(1,1,:) = (/b,a/);   sf%sfe_funs(1,2,:) = (/c,c/);   sf%sfe_funs(1,3,:) = (/a,b/);
    sf%sfe_funs(2,1,:) = (/a,b/);   sf%sfe_funs(2,2,:) = (/b,a/);   sf%sfe_funs(2,3,:) = (/c,c/);
    sf%sfe_funs(3,1,:) = (/c,c/);   sf%sfe_funs(3,2,:) = (/a,b/);   sf%sfe_funs(3,3,:) = (/b,a/);
    sf%wei = 1./3.

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
    real, dimension(ndim,nloc,ngi)::ref_sfdev
    real :: a2
    integer :: inod, gi

    ! print*, 'is input alright?', 'element' , lelement, ' mesh veretex' , meshvertex

    ! first we fetech Nodes coordinate in local element
    do inod = 1,nloc
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
    lsf%detwei(:) = sf%wei(:) * abs(a2)/2.
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
    ! \nabla N (->on reference element): for linear element, 
    ! derivatives on reference element is constant
    do gi = 1,ngi
      ref_sfdev(:,1,gi) = (/-1.,-1./)
      ref_sfdev(:,2,gi) = (/1.,0./)
      ref_sfdev(:,3,gi) = (/0.,1./)
    enddo
    ! print*, lsf%dev_funs(:,1,1)
    ! left multiply J^-T to transform to local element
    do gi = 1,ngi
      lsf%dev_funs(:,:,gi) = matmul( jacobit , ref_sfdev(:,:,gi) )
    enddo
    ! print*, lsf%dev_funs
  end subroutine calc_local_shape_func

  subroutine calc_local_surf_sf(lsf, sf, lface, meshvertex)
    ! this is to calculate shape function value on a surface
    type(shape_func_dev), intent(inout)::lsf 
    type(shape_func), intent(in)::sf
    type(face), intent(in):: lface 
    type(vertex), dimension(:), intent(in):: meshvertex
    real, dimension(2,2):: x_all 
    integer :: iface
    
    ! print*, 'input inside face', lface
    ! first retrieve coordinates
    x_all(:,1) = meshvertex( (lface%vertex(1)) )%coor(:)
    x_all(:,2) = meshvertex( (lface%vertex(2)) )%coor(:)
    ! print*, x_all, 'sqrt', sqrt( (x_all(1,1)-x_all(1,2))**2 + (x_all(2,1)-x_all(2,2))**2 )
    ! edge length * weight
    lsf%sdetwei = sqrt( (x_all(1,1)-x_all(1,2))**2 + (x_all(2,1)-x_all(2,2))**2 ) * 1. /2.
    ! print*, 'edge length',  lsf%sdetwei(1)
    ! 2 derivatives of 3 local nodes shape functions evaluated on 3 edges, 2 gaussian point for each
    do iface = 1,3
      ! print*, lsf%sdev_funs(:,:,iface,:)
      ! print*, lsf%dev_funs(:,:,:)
      lsf%sdev_funs(:,:,iface,:) = lsf%dev_funs(:,:,1:2)  ! because derivative is constant, 
          ! we take a short cut and applies its value at (2D) quadrature points 
          ! to (1D) quadrature points on edges
    enddo
  end subroutine calc_local_surf_sf
    

end module