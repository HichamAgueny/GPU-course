      program enterdata

      implicit none

       integer :: i,j,k
       integer :: iter,count_rate, count_max,count
       integer :: t_start,t_final
       integer, parameter :: nx=8192,ny=nx,max_iter=2000
       double precision, parameter    :: pi=4d0*dtan(1d0) 
       real, parameter    :: error=0.001
       double precision               :: max_err,time_s,&
                                         d2fx,d2fy
       double precision, allocatable  :: f(:,:), f_k(:,:),f0(:,:)

       allocate(f(nx,ny)); allocate(f_k(nx,ny)); allocate(f0(nx,ny))

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)
     
!generate a random vector to define the initial conditions
!       call system_clock(count, count_rate, count_max)

       do i=1,nx
          f(i,1)  = dsin(pi*(i-1)/(nx-1))
          f(i,ny) =0d0
       enddo
       do j=1,ny
          f(1,j) = dsin(pi*(j-1)/(ny-1))
          f(nx,j) = 0d0
       enddo

        write(*,*)'--sum before loop:', sum(f(:,:))

!offload data to the device
!$acc enter data copyin(f) create(f_k)

!$acc parallel loop
         do j=1,ny
            do i=1,nx

               f_k(i,j) = 2*f(i,j)

             enddo
          enddo
!$acc end parallel

!offload data to the device
!$acc enter data create(f0)
!$acc parallel loop
         do j=1,ny
            do i=1,nx

               f0(i,j) = 2*f(i,j) + f_k(i,j)

             enddo
          enddo
!$acc end parallel

!$acc exit data copyout(f_k) delete(f)
!$acc exit data copyout(f0)
 
       write(*,*)'--sum after end data fk:', sum(f_k(:,:))
       write(*,*)'--sum after end data f0:',sum(f0(:,:)) 

       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       !print*, '--Time it takes (s)', time_s
       end


       subroutine arrayf_new(nx,ny,f,f_k)

         implicit none

         integer :: i,j,nx,ny
         double precision :: f(nx,ny),f_k(nx,ny)

         do j=2,ny-1
            do i=2,nx-1

               f_k(i,j) = 2*f(i,j) 

             enddo
          enddo

          end subroutine arrayf_new
