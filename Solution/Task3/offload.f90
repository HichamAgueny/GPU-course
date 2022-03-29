      program offload

      implicit none

       integer :: i,j,k
       integer :: iter,count_rate, count_max,count
       integer :: t_start,t_final
       integer, parameter :: nx=8192,ny=nx,max_iter=2000
       double precision, parameter    :: pi=4d0*dtan(1d0) 
       real, parameter    :: error=0.001
       double precision               :: max_err,time_s,&
                                         d2fx,d2fy
       double precision, allocatable  :: f(:,:), f_k(:,:)

       allocate(f(nx,ny)); allocate(f_k(nx,ny))

!generate a random vector 
      CALL RANDOM_NUMBER(f)

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)
     
        write(*,*)'--sum before loop:', sum(f(:,:))/(nx*ny)

!offload the data to the device and copy the result back to the host
!$acc data copyin(f) copyout(f_k)
      do iter=1,5

!$acc parallel loop
         do j=1,ny
            do i=1,nx

               f_k(i,j) = 2*f(i,j)

             enddo
          enddo
!$acc end parallel

         write(*,*)'--sum:No update', sum(f_k(:,:))/(nx*ny)

!copy from the device to host
!$acc update host(f_k)
          print*,"--Sum: update self",sum(f_k(:,:))/(nx*ny)

          f_k = f_k/2.
!copy from the host to device
!$acc update device(f_k)
       enddo
!$acc end data

       write(*,*)'--sum after end data:', sum(f_k(:,:))/(nx*ny)
       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       print*, '--Time it takes (s)', time_s
       end


       subroutine arrayf_new(nx,ny,f,f_k)

         implicit none

         integer :: i,j,nx,ny
         double precision :: f(nx,ny),f_k(nx,ny)

         do j=1,ny
            do i=1,nx

               f_k(i,j) = 2*f(i,j) 

             enddo
          enddo

          end subroutine arrayf_new

