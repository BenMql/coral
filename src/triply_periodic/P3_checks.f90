
   ! perform checks
   if (verbose) print*, '=========================='
   if (verbose) print*, '... Now performing some checks'
   if (verbose) print*, '=========================='
   if (verbose) print*, '... '                                                                   
   if (verbose) print*, '... First checking the transforms'
   if (verbose) print*, '... Method: initialise one physical fields F. Transform to spectral and back'
   if (verbose) print*, '...         G = iDFT( DFT( F )) '
   if (verbose) print*, '... For each process, we compute the error of the first three moments '
   if (verbose) print*, '... of G and F normalized by the moments of F. '
   if (verbose) print*, '... More explicitly, for n=1,2,3, this quantity is:'
   if (verbose) print*, '... '                                                                   
   if (verbose) print*, '=========================='
   if (verbose) print*, '          NX,NY,NZ                   NX,NY,NZ              '
   if (verbose) print*, '            ===                        ===                 '
   if (verbose) print*, '            \                 n        \                 n '
   if (verbose) print*, '            /   [F(iz,iy,ix)]    __    /   [G(iz,iy,ix)]   '
   if (verbose) print*, '            ===                        ===                 '
   if (verbose) print*, '          ix,iy,iz                   ix,iy,iz              '
   if (verbose) print*, ' E  = ______________________________________________________'
   if (verbose) print*, '  n                NX,NY,NZ '
   if (verbose) print*, '                     ===     '
   if (verbose) print*, '                     \                 n '
   if (verbose) print*, '                     /     F(iz,iy,ix)   '
   if (verbose) print*, '                     ===                             '
   if (verbose) print*, '                   ix,iy,iz '
   if (verbose) print*, '=========================='
   
   do i3 = lbound(main%linear_variables(1)%phys, dim=3),&
           ubound(main%linear_variables(1)%phys, dim=3)
   do i2 = lbound(main%linear_variables(1)%phys, dim=2),&
           ubound(main%linear_variables(1)%phys, dim=2)
   do i1 = lbound(main%linear_variables(1)%phys, dim=1),&
           ubound(main%linear_variables(1)%phys, dim=1)
   main%linear_variables(1)%phys(i1,i2,i3) = cos(my_rank+cos(real(i1))+sin(real(i2))+cos(cos(real(i3))))
   end do
   end do
   end do
   !print *, main%linear_variables(1)%phys
   first_moment = sum(main%linear_variables(1)%phys)
   second_moment = sum((main%linear_variables(1)%phys)**2)
   third_moment = sum((main%linear_variables(1)%phys)**3)
   call main%linear_variables(1)%phys_to_spec()
   call main%linear_variables(1)%spec_to_phys()
   !print*, my_rank, first_moment, &
   !                 sum(main%linear_variables(1)%phys),&
   !                 abs(first_moment- sum(main%linear_variables(1)%phys)),&
   !                 second_moment,&
   !                 sum((main%linear_variables(1)%phys)**2),&
   !                 abs(second_moment- sum((main%linear_variables(1)%phys)**2)),&
   !                 third_moment,&
   !                 sum((main%linear_variables(1)%phys)**3),&
   !                 abs(third_moment- sum((main%linear_variables(1)%phys)**2))
   
   4262 Format ('core ',I4,' | E1=',ES10.3,' | E2=',ES10.3,' | E3=',ES10.3)
   write(*,4262)  my_rank,& 
                    abs((first_moment- sum(main%linear_variables(1)%phys))/first_moment),&
                    abs((second_moment- sum((main%linear_variables(1)%phys)**2))/second_moment),&
                    abs((third_moment- sum((main%linear_variables(1)%phys)**3))/third_moment)
