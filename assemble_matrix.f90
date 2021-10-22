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
    integer :: ele, iface, iface2, inod, jnod, &
               glob_i, glob_j, idim, glob_iface,glob_iface2, &
               ele2,  third_nod, gi 
    logical :: sgi_order_same
    real :: normal(2) ! normal vector
    real :: nn, nnx, nxn, nxnx, ml ! ml is local mass?
    real, dimension(2,3) :: x_all ! nodes coordinate
    real :: elength ! edge length
    real, parameter :: K(2,2) = reshape((/1.0,0.0,0.0,1.0/) , (/ndim,ndim/)), &  ! diffusion coefficient
        alpha = 0.  ,  &    ! reaction term in govering equation
        beta = 1.   ,  &    ! J0 term: 1 - w/o superpenalization; >1 - w/ superpenalization.
        sigma = 1.  ,  &    ! J0 term
        epsilon = 1.        ! symmetry/asymetry interior penalty

    ! shape/basis functions
    type(shape_func) :: sf 
    type(shape_func_dev) :: sf_dev, sf_dev2 ! neighbour element shape function derivatives

    ! total number of elements ===> no need to do this since nele is stored in mesh_type module and defined after reading mesh
    ! nele=size(meshele)
    ! print*, 'nele', nele

    call init_shape_func(sf)

    do ele = 1,nele 
      ! calculate shape functions'derivative in local element, determinant of Jacobian
      ! write(20, *) '==============ele:', ele, '============'
      call calc_local_shape_func(sf_dev, sf, meshele(ele), meshvertex)
      ! write(20, *) 'ele=', ele, 'sf_dev = ', sf_dev
      ! local node cycle
      do inod = 1,nloc
        ! global node number
        glob_i = (ele-1)*3 + inod 
        ! \sum_j (Ni Nj) -> this is used in rhs f
        ml = 0.
        ! local node cycle 2
        do jnod = 1,nloc
          ! to store local \nabla Ni \cdot \nabla Nj
          nxnx = 0.
          ! global node number
          glob_j = (ele-1)*3 + jnod 
          ! dot product of divergence
          do idim = 1,ndim
            nxnx = nxnx + sum( sf_dev%dev_funs(idim, inod, :) * sf_dev%dev_funs(idim, jnod, :) & 
              * sf_dev%detwei(:) )
              ! write(20,*) 'nxnx', nxnx
              ! write(20,*) sf_dev%dev_funs(idim, inod, :) , sf_dev%dev_funs(idim, inod, :) ,                sf_dev%detwei(:)
          enddo
          ! mass matrix Ni Nj
          nn = sum( sf%fe_funs(inod,:) * sf%fe_funs(jnod,:) * sf_dev%detwei(:) )
          ! print*, inod, jnod, glob_i, glob_j, nn
          ml = ml + nn
          ! contribution A = K \nabla Ni \cdot \nabla Nj + \alpha Ni Nj
          bigm(glob_i, glob_j) = bigm(glob_i, glob_j) + nxnx + nn * alpha
          ! print*,  'ele', ele, 'glob_i', glob_i, 'glob_j', glob_j, 'bigm', bigm(glob_i, glob_j), 'nxnx' , nxnx
          ! print*, 'nxnx, nn, ml', nxnx, nn, ml

          ! here we didn't use a direct local coordinate to evaluate right hand side term f
          ! (as used in Rivere's book), but instead we use a shape function interpolation.
          ! this is because practically we will only give right hand side f value on each 
          ! node.
          rhs(glob_i) = rhs(glob_i) + nn * f1( meshvertex(glob_j)%coor(1) , meshvertex(glob_j)%coor(2) )
          ! print*, 'ele, inod, rhs, f1', ele, inod, rhs(glob_i), f1( meshvertex(glob_j)%coor(1) , meshvertex(glob_j)%coor(2) )
        enddo
      enddo
      ! enddo !nele

      ! ! element cycle
      ! do ele = 1, nele 
      ! face cycle
      do iface = 1,nloc
        glob_iface = meshele(ele)%face(iface)
        if (meshface(glob_iface)%bctype .eq. 0) then ! not boundary face
          ! element index on other sides of the face
          ele2 = meshface(glob_iface)%neighbor(2)
          ! face info on the other side
          iface2 = meshface(glob_iface)%nb_iface
          glob_iface2 = meshele(ele)%face(iface2)
          ! are gaussian point distribution in same order for two sides?
          sgi_order_same = all(meshvertex(meshface(glob_iface)%vertex(1))%coor & 
              .eq. meshvertex(meshface(glob_iface2)%vertex(1))%coor) 
          ! print*, 'ele', ele, 'glob_iface', glob_iface, 'ele2',  ele2
          call calc_local_shape_func(sf_dev2, sf, meshele(ele2), meshvertex)
          ! calculate derivatives @ edges (suspicious! watch out TODO)
          call calc_local_surf_sf(sf_dev, sf, meshface(glob_iface), meshvertex)
          call calc_local_surf_sf(sf_dev2, sf, meshface(glob_iface), meshvertex)
          ! nodes coordinate of the face
          x_all(:,1) = meshvertex( meshface(glob_iface)%vertex(1) )%coor(:)
          x_all(:,2) = meshvertex( meshface(glob_iface)%vertex(2) )%coor(:)
          ! this is the third node opposite to the face
          third_nod = iface - 1
          if (third_nod .eq. 0) third_nod=3
          x_all(:,3) = meshvertex( (ele-1)*3 + third_nod )%coor(:) 
          ! normal direction of the face n_e
          elength = sqrt( (x_all(1,2)-x_all(1,1))**2 + (x_all(2,2)-x_all(2,1))**2 ) ! sqrt( (x2-x1)^2 + (y2-y1)^2 )
          normal(1) = (x_all(2,1) - x_all(2,2)) / elength ! (y1-y2)/e
          normal(2) = (x_all(1,2) - x_all(1,1)) / elength ! (x2-x1)/e
          ! is this the outer normal?
          ! n * (V31) < 0 then outer; otherwise inner
          if ( normal(1)*(x_all(1,3)-x_all(1,1)) + normal(2)*(x_all(2,3)-x_all(2,1)) .lt. 0. ) then 
            ! do nothing
          else
            normal(1) = - normal(1)
            normal(2) = - normal(2)
          endif
          ! print*,  'ele fac enormal ', ele , iface, normal(1), normal(2) ! pass review!
          ! print*, elength, sf_dev%sdetwei(:), sf_dev2%sdetwei(:)

          ! start assembling
          ! contribution m11
          do inod = 1,nloc
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nloc
              glob_j = (ele-1)*3 + jnod 
              nnx=0.
              nxn=0.
              do idim = 1,ndim
                ! \nabla P1.n_e.n1
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev%sdetwei(:) ) * normal(idim) 
                ! print*, 'sdev_funs', sf_dev%sdev_funs(idim, jnod, iface, :)
                ! \nabla n1.n_e.P1
                nxn = nxn + sum( sf_dev%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev%sdetwei(:) ) * normal(idim)
              enddo
              ! P1.n1
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                - 0.5 *nnx + 0.5 * epsilon * nxn & 
                + sigma/(elength**beta) * nn 
              ! print*, 'nnx', nnx, 'nxn', nxn, 'elenght', elength, 'nn', nn
              ! print*,  'ele', ele, 'glob_i', glob_i, 'glob_j', glob_j, 'bigm', bigm(glob_i, glob_j)
            enddo
          enddo
          ! contribution m22
          do inod = 1,nloc
            glob_i = (ele2-1)*3 + inod 
            do jnod = 1,nloc
              glob_j = (ele2-1)*3 + jnod 
              nnx=0.
              nxn=0.
              do idim = 1,ndim
                ! \nabla P2.n_e.n2
                nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev2%sdev_funs(idim, jnod, iface, :) &
                     * sf_dev2%sdetwei(:) ) * normal(idim) 
                ! \nabla n2.n_e.P2
                nxn = nxn + sum( sf_dev2%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                     * sf_dev2%sdetwei(:) ) * normal(idim)
              enddo
              ! P2.n2
              nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev2%sdetwei(:) )
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                + 0.5 *nnx - 0.5 * epsilon * nxn & 
                + sigma/(elength**beta) * nn 
            enddo
          enddo

          ! now surface term across elements.
          ! in m12 and m21 we need neighbouring face local index (iface2) to correctly use sf%sfe_funs.
          ! we also need to know if gaussian points distribution is the same on both sides (sgi_order_same)
          ! contribution m12
          do inod = 1,nloc
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nloc
              glob_j = (ele2-1)*3 + jnod 
              nnx=0.
              nxn=0.
              do idim = 1,ndim
                if (sgi_order_same) then 
                ! \nabla P2.n_e.n1
                  nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev2%sdev_funs(idim, jnod, iface, :) &
                      * sf_dev2%sdetwei(:) ) * normal(idim) 
                  ! \nabla n1.n_e.P2
                  nxn = nxn + sum( sf_dev%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                      * sf_dev%sdetwei(:) ) * normal(idim) ! it doesn't really matter if 
                        ! we use sf_dev%sdetwei or sf_dev2%sdetwei because they are the same
                        ! as long as nodes are same on both sides
                else
                  do gi = 1,nsgi
                    ! \nabla P2.n_e.n1
                    nnx = nnx + sf%sfe_funs(inod,iface,gi) * sf_dev2%sdev_funs(idim, jnod, iface2, nsgi+1 - gi) &
                        * sf_dev2%sdetwei(nsgi+1-gi) * normal(idim)
                    ! \nabla n1.n_e.P2
                    nxn = nxn + sf_dev%sdev_funs(idim, inod, iface, gi) * sf%sfe_funs(jnod,iface2, nsgi+1 - gi) & 
                        * sf_dev%sdetwei(gi)  * normal(idim)
                  enddo
                endif
              enddo
              ! P2.n1
              if (sgi_order_same) then
                nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface2, :) * sf_dev2%sdetwei(:) )
                  ! TODO - suspicious. Do we need to identify which nodes correspond to which on the other side?
                  ! Seems we do! Or the golbi and globj won't be right!.
                  ! unless we don't use this nn surface term!
              else
                nn = 0.
                do gi = 1,nsgi 
                  nn = nn + sf%sfe_funs(inod,iface,gi) * sf%sfe_funs(jnod,iface2,nsgi+1-gi) * sf_dev2%sdetwei(nsgi+1-gi)
                enddo
              endif
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                - 0.5 *nnx - 0.5 * epsilon * nxn & 
                - sigma/(elength**beta) * nn 
            enddo
          enddo
          ! contribution m21
          do inod = 1,nloc
            glob_i = (ele2-1)*3 + inod 
            do jnod = 1,nloc
              glob_j = (ele-1)*3 + jnod 
              nnx=0.
              nxn=0.
              do idim = 1,ndim
                if (sgi_order_same) then 
                  ! \nabla P1.n_e.n2
                  nnx = nnx + sum( sf%sfe_funs(inod,iface, :) * sf_dev%sdev_funs(idim, jnod, iface, :) &
                      * sf_dev%sdetwei(:) ) * normal(idim) 
                  ! \nabla n2.n_e.P1
                  nxn = nxn + sum( sf_dev2%sdev_funs(idim, inod, iface, :) * sf%sfe_funs(jnod,iface, :) & 
                      * sf_dev2%sdetwei(:) ) * normal(idim)
                else
                  do gi = 1,nsgi
                    ! \nabla P1.n_e.n2
                    nnx = nnx + sf%sfe_funs(inod,iface2,nsgi+1-gi) * sf_dev%sdev_funs(idim,jnod,iface,gi) &
                        * sf_dev%sdetwei(gi) * normal(idim)
                    ! nabla n2.n_e.P1
                    nxn = nxn + sf_dev2%sdev_funs(idim,inod,iface2,nsgi+1-gi) * sf%sfe_funs(jnod,iface,gi) &
                        * sf_dev2%sdetwei(nsgi+1-gi) * normal(idim)
                  enddo 
                endif
              enddo
              if (sgi_order_same) then 
                ! P1.n2
                nn = sum( sf%sfe_funs(inod,iface, :) * sf%sfe_funs(jnod,iface, :) * sf_dev2%sdetwei(:) )
              else
                nn = 0.
                do gi = 1,nsgi 
                  ! P1.n2
                  nn = nn + sf%sfe_funs(inod,iface2,nsgi+1-gi) * sf%sfe_funs(jnod,iface,gi) * sf_dev2%sdetwei(nsgi+1-gi)
                enddo
              endif
              bigm(glob_i, glob_j) = bigm(glob_i, glob_j) &
                + 0.5 *nnx + 0.5 * epsilon * nxn & 
                - sigma/(elength**beta) * nn 
            enddo
          enddo
        elseif( meshface(glob_iface)%bctype .eq. 1) then ! Dirichlet boundary face
          ! calculate derivatives @ edges (suspicious! watch out TODO)
          call calc_local_surf_sf(sf_dev, sf, meshface(glob_iface), meshvertex)
          ! only m11 contribution
          do inod = 1,nloc
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nloc
              glob_j = (ele-1)*3 + jnod 
              nnx=0.
              nxn=0.
              do idim = 1,ndim
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
                - nnx + epsilon * nxn & 
                - sigma/(elength**beta) * nn 
              ! right hand side -- boundary condition
              rhs(glob_i) = rhs(glob_i) + epsilon *nxn + sigma/elength**beta*nn &
                * g_D( meshvertex(glob_j)%coor(1), meshvertex(glob_j)%coor(2) )
            enddo
          enddo
          
        elseif( meshface(glob_iface)%bctype .eq. 2) then ! Neumann boundary face
          print*, 'Caution: You are using Neumann boundary conditions now. Please check if the code is correctly implemented.' 
          ! to be tested... (because we don't mark meshface%bctype with 2 now)
          call calc_local_surf_sf(sf_dev, sf, meshface(glob_iface), meshvertex)
          do inod = 1,nloc 
            glob_i = (ele-1)*3 + inod 
            do jnod = 1,nloc 
              glob_j = (ele-1)*3 + jnod 
              nn = sum( sf%sfe_funs(inod, iface, :) * sf%sfe_funs(jnod, iface, :) * sf_dev%sdetwei(:) )
              ! when actually use this, change g_N to the desired BC function.
              ! rhs(glob_i) = rhs(glob_i) + nn*g_N( meshvertex(glob_j)%coor(1), meshvertex(glob_j)%coor(2) )
            enddo
          enddo
        endif ! face type (0 interior 1 Dirichlet 2 Neumann)
      enddo ! iface = 1,3
    enddo ! ele=1,nele
    
    return
  end subroutine

  real function f1(x,y) result(f)
    real, intent(in):: x,y  ! coordinate to be evaluated
    ! real, intent(out) :: f
    f = -exp(-x-y**2.) - (4.*y**2 - 2.) * exp(-x-y**2)
    ! analytical solution : p = exp(-x-y**2)
  end function

  real function f2(x,y) result(f)
    real, intent(in):: x,y  ! coordinate to be evaluated
    ! real, intent(out) :: f
    f = -2.*exp(-x**2 - y**2) * (2*x**4-2*x**3-5*x**2+3*x+1 ) * (y-1) & 
      - -2.*exp(-x**2 - y**2) * (2*y**4-2*y**3-5*y**2+3*y+1 ) * (x-1)
    ! analytical solution : p = x(x-1)y(y-1)exp(-x^2-y^2)
  end function

  real function f3(x,y) result(f)
    real, intent(in):: x,y
    f = 0.
  end function

  real function g_D(x,y) result(g)
    real, intent(in):: x, y
    g = exp(-x-y**2) ! for f1
    ! g = x*(x-1)*y*(y-1)*exp(-x**2-y**2) ! for f2
    ! g =1. ! for f3
  end function
end module