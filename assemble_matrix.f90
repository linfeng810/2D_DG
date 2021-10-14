! this is to assemble left hand side and right hand side 

module assemb_matrix_engine
  use mesh_type
  use basis_function

  implicit none 

  contains
  subroutine assemb_matrix(meshele, meshface, meshvertex, bigm, rhs)
    type(element), dimension(:), intent(in) :: meshele 
    type(face), dimension(:), intent(in) :: meshface
    type(vertex), dimension(:), intent(in) :: meshvertex 
    real, dimension(:,:), intent(inout) :: bigm ! lhs matrix, declared out side 
      ! this scope and passed in here.
    real, dimension(:), intent(inout) :: rhs ! rhs vector, :
    integer :: ele, nele, iface, inod, jnod, nnod, &
               glob_i, glob_j, idim, glob_iface, ele2 
    real :: normal(2) ! normal vector
    real :: nn, nnx, nxn, nxnx, ml ! ml is local mass?
    real, dimension(2,3) :: x_all ! nodes coordinate
    real :: elength ! edge length
    real, parameter :: K(2,2) = reshape((/1.0,0.0,0.0,1.0/) , (/2,2/)), &  ! diffusion coefficient
        alpha = 0.  ,  &    ! reaction term in govering equation
        beta = 1.   ,  &    ! J0 term: 1 - w/o superpenalization; >1 - w/ superpenalization.
        sigma = 1.  ,  &    ! J0 term
        epsilon = 1.        ! symmetry/asymetry interior penalty

    ! shape/basis functions
    type(shape_func) :: sf 
    type(shape_func_dev) :: sf_dev, sf_dev2 ! neighbour element shape function derivatives

    ! total number of elements
    nele=size(meshele)

    call init_shape_func(sf)

    do ele = 1,nele 
      ! calculate shape functions'derivative in local element, determinant of Jacobian
      call calc_local_shape_func(sf_dev, sf, meshele(ele), meshvertex)
      ! to store local \nabla Ni \cdot \nabla Nj
      nxnx = 0.
      ! local node cycle
      do inod = 1,3
        ! global node number
        glob_i = (ele-1)*3 + inod 
        ! \sum_j (Ni Nj) -> this is used in rhs f
        ml = 0.
        ! local node cycle 2
        do jnod = 1,3 
          ! global node number
          glob_j = (ele-1)*3 + jnod 
          ! dot product of divergence
          do idim = 1,2
            nxnx = nxnx + sum( sf_dev%dev_funs(idim, inod, :) * sf_dev%dev_funs(idim, inod, :) & 
              * sf_dev%detwei(:) )
          enddo
          ! mass matrix Ni Nj
          nn = sum( sf%fe_funs(inod,:) * sf%fe_funs(jnod,:) * sf_dev%detwei(:) )
          ml = ml + nn
          ! contribution A = K \nabla Ni \cdot \nabla Nj + \alpha Ni Nj
          bigm(glob_i, glob_j) = bigm(glob_i, glob_j) + nxnx * K(inod, jnod) + nn * alpha
        enddo
        ! here we didn't use a direct local coordinate to evaluate right hand side term f
        ! (as used in Rivere's book), but instead we use a shape function interpolation.
        ! this is because practically we will only give right hand side f value on each 
        ! node.
        rhs(glob_i) = ml * f1( meshvertex(glob_i)%coor(1) , meshvertex(glob_i)%coor(1) )
      enddo
      ! enddo !nele

      ! ! element cycle
      ! do ele = 1, nele 
      ! face cycle
      do iface = 1,3
        glob_iface = meshele(ele)%face(iface)
        if (meshface(glob_iface)%bctype .eq. 0) then ! not boundary face
          ! element index on other sides of the face
          ele2 = meshface(glob_iface)%neighbor(2)
          call calc_local_shape_func(sf_dev2, sf, meshele(ele2), meshvertex)
          ! calculate derivatives @ edges (suspicious! watch out TODO)
          call calc_local_surf_sf(sf_dev, sf, meshface(glob_iface), meshvertex)
          call calc_local_surf_sf(sf_dev2, sf, meshface(glob_iface), meshvertex)
          ! nodes coordinate of the face
          x_all(:,1) = meshvertex( meshface(glob_iface)%vertex(1) )%coor(:)
          x_all(:,2) = meshvertex( meshface(glob_iface)%vertex(2) )%coor(:)
          ! this is the third node opposite to the face
          x_all(:,3) = meshvertex( 6-(meshface(glob_iface)%vertex(2))-(meshface(glob_iface)%vertex(1)) )%coor(:) 
          ! normal direction of the face n_e
          elength = sqrt( (x_all(2,1)-x_all(1,1))**2 + (x_all(2,2)-x_all(1,2)**2) )
          normal(1) = (x_all(1,2) - x_all(2,2)) / elength ! (y1-y2)/e
          normal(2) = (x_all(2,1) - x_all(1,1)) / elength ! (x2-x1)/e
          ! is this the outer normal?
          ! n * (V31) < 0 then outer; or inner
          if ( normal(1)*(x_all(1,3)-x_all(1,1)) + normal(2)*(x_all(2,3)-x_all(2,1)) .lt. 0. ) then 
            ! do nothing
          else
            normal(1) = - normal(1)
            normal(2) = - normal(2)
          endif

          ! start assembling
          ! contribution m11
          do inod = 1,nnod
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nnod
              glob_j = (ele-1)*3 + jnod 
              do idim = 1,2
                ! \nabla P1.n_e.n1
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev%sdetwei(:) ) * normal(idim) 
                ! \nabla n1.n_e.P1
                nxn = nxn + sum( sf_dev%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev%sdetwei(:) ) * normal(idim)
              enddo
              ! P1.n1
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                - 0.5 * K(inod,jnod)*nnx + 0.5 * epsilon * K(inod, jnod)*nxn & 
                + sigma/elength**beta * nn 
            enddo
          enddo
          ! contribution m22
          do inod = 1,nnod
            glob_i = (ele2-1)*3 + inod 
            do jnod = 1,nnod
              glob_j = (ele2-1)*3 + jnod 
              do idim = 1,2
                ! \nabla P1.n_e.n1
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev2%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev2%sdetwei(:) ) * normal(idim) 
                ! \nabla n1.n_e.P1
                nxn = nxn + sum( sf_dev2%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev2%sdetwei(:) ) * normal(idim)
              enddo
              ! P1.n1
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev2%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                + 0.5 * K(inod,jnod)*nnx - 0.5 * epsilon * K(inod, jnod)*nxn & 
                + sigma/elength**beta * nn 
            enddo
          enddo
          ! contribution m12
          do inod = 1,nnod
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nnod
              glob_j = (ele2-1)*3 + jnod 
              do idim = 1,2
                ! \nabla P1.n_e.n1
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev2%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev2%sdetwei(:) ) * normal(idim) 
                ! \nabla n1.n_e.P1
                nxn = nxn + sum( sf_dev%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev%sdetwei(:) ) * normal(idim) ! it doesn't really matter if 
                      ! we use sf_dev%sdetwei or sf_dev2%sdetwei because they are the same
                      ! as long as nodes are same on both sides
              enddo
              ! P1.n1
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev2%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                - 0.5 * K(inod,jnod)*nnx - 0.5 * epsilon * K(inod, jnod)*nxn & 
                - sigma/elength**beta * nn 
            enddo
          enddo
          ! contribution m21
          do inod = 1,nnod
            glob_i = (ele2-1)*3 + inod 
            do jnod = 1,nnod
              glob_j = (ele-1)*3 + jnod 
              do idim = 1,2
                ! \nabla P1.n_e.n1
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev%sdetwei(:) ) * normal(idim) 
                ! \nabla n1.n_e.P1
                nxn = nxn + sum( sf_dev2%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev2%sdetwei(:) ) * normal(idim)
              enddo
              ! P1.n1
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev2%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                + 0.5 * K(inod,jnod)*nnx + 0.5 * epsilon * K(inod, jnod)*nxn & 
                - sigma/elength**beta * nn 
            enddo
          enddo
        elseif( meshface(glob_iface)%bctype .eq. 1) then ! Dirichlet boundary face
          ! calculate derivatives @ edges (suspicious! watch out TODO)
          call calc_local_surf_sf(sf_dev, sf, meshface(glob_iface), meshvertex)
          ! only m11 contribution
          do inod = 1,nnod
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nnod
              glob_j = (ele-1)*3 + jnod 
              do idim = 1,2
                ! \nabla P1.n_e.n1
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev%sdetwei(:) ) * normal(idim) 
                ! \nabla n1.n_e.P1
                nxn = nxn + sum( sf_dev%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev%sdetwei(:) ) * normal(idim)
              enddo
              ! P1.n1
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                - K(inod,jnod)*nnx + epsilon * K(inod, jnod)*nxn & 
                - sigma/elength**beta * nn 
              ! right hand side -- boundary condition
              rhs(glob_j) = rhs(glob_j) + epsilon * K(inod, jnod)*nnx + sigma/elength**beta*nn &
                * g_D( meshvertex(glob_j)%coor(1), meshvertex(glob_j)%coor(2) )
            enddo
          enddo
          
        elseif( meshface(glob_iface)%bctype .eq. 2) then ! Neumann boundary face
          ! to be filled...
        endif ! face type (0 interior 1 Dirichlet 2 Neumann)
      enddo ! iface = 1,3
    enddo ! ele=1,nele
    

  end subroutine

  real function f1(x,y) result(f)
    real, intent(in):: x,y  ! coordinate to be evaluated
    ! real, intent(out) :: f
    f = exp(-x-y**2.) + (4.*y**2 - 2.) * exp(-x-y**2)
    ! analytical solution : p = exp(-x-y**2)
  end function

  real function f2(x,y) result(f)
    real, intent(in):: x,y  ! coordinate to be evaluated
    ! real, intent(out) :: f
    f = exp(-x**2 - y**2) * ( 2.*(x - 1.)*x*(2.*y**4 - 2.*y**3 - 5.*y**2 + 3.*y + 1.) & 
      + 2.*(2.*x**4 - 2.*x**3 - 5.**2 + 3.*x + 1.) )
    ! analytical solution : p = x(x-1)y(y-1)exp(-x^2-y^2)
  end function

  real function g_D(x,y) result(g)
    real, intent(in):: x, y
    g = exp(-x-y**2) ! for f1
    ! g = x*(x-1)*y*(y-1)*exp(-x**2-y**2) ! for f2
  end function
end module