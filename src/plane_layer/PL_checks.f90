
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

   call wclock%init()
   call main%linear_variables(1)%phys_to_spec()
   call main%linear_variables(1)%spec_to_phys()
   call wclock%check(.false.)

   first_moment = sum(main%linear_variables(1)%phys)
   second_moment = sum((main%linear_variables(1)%phys)**2)
   third_moment = sum((main%linear_variables(1)%phys)**3)


   if (my_rank==0) print *, " >>> BEFORE THE TRANSFORMS:"

   4253 Format ('core ',I4,', moments before transforms | m1=',ES10.3,' | m2=',ES10.3,' | m3=',ES10.3)
   do irank = 0, world_size
   if (my_rank == irank) then
   write(*,4253)  my_rank,& 
                    first_moment,&
                    second_moment,&
                    third_moment
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   end do

   call wclock%init()
   call main%linear_variables(1)%phys_to_spec()
   call main%linear_variables(1)%spec_to_phys()
   call wclock%check(.false.)
   
   if (my_rank==0) print *, " >>> AFTER THE TRANSFORMS:"

   4254 Format ('core ',I4,', moments after transforms | m1=',ES10.3,' | m2=',ES10.3,' | m3=',ES10.3)
   do irank = 0, world_size
   if (my_rank == irank) then
   write(*,4254)  my_rank,& 
                    sum(main%linear_variables(1)%phys),&
                    sum((main%linear_variables(1)%phys)**2),&
                    sum((main%linear_variables(1)%phys)**3)
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   end do

   if (my_rank==0) print *, " >>> RELATIVE ERRORS (should be approx. machine precision):"

   4262 Format ('core ',I4,' | E1=',ES10.3,' | E2=',ES10.3,' | E3=',ES10.3)
   do irank = 0, world_size
   if (my_rank == irank) then
   write(*,4262)  my_rank,& 
                    abs((first_moment- sum(main%linear_variables(1)%phys))/first_moment),&
                    abs((second_moment- sum((main%linear_variables(1)%phys)**2))/second_moment),&
                    abs((third_moment- sum((main%linear_variables(1)%phys)**3))/third_moment)
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   end do

   if (my_rank==0) print *, " >>> TIMING:"
   if (my_rank.eq.0) then
      print *, "One forward/backward transform takes: ", wclock% time_now - wclock% start_time, " seconds"
   end if
   call wclock%init()
   do i1 = 1,10
   call main%linear_variables(1)%phys_to_spec()
   call main%linear_variables(1)%spec_to_phys()
   end do
   call wclock%check(.false.)
   if (my_rank.eq.0) then
   print *, "Ten forward/backward transforms take: ", wclock% time_now - wclock% start_time, " seconds"
   end if

