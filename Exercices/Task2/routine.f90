      program routine

      implicit none

       integer :: i,j,k
       integer :: iter,count_rate, count_max,count
       integer :: t_start,t_final
       integer, parameter :: nx=8192,ny=nx
       double precision, parameter    :: pi=4d0*dtan(1d0) 
       double precision               :: max_err,time_s,&
                                         sum_x,total
                                         
       double precision, allocatable  :: f(:,:), f_k(:,:)

       allocate(f(nx,ny)); allocate(f_k(nx,ny))

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)
     
!generate a random vector

       CALL RANDOM_NUMBER(f)

       print*, "--Start the iteration"
       max_err=0.0

!parallelise the loop and include a directive to offload the function
!arrayf_k
       do iter=1,5
          call arrayf_k(nx,ny,f,f_k)
         do j=1,ny
            do i=1,nx

               max_err = max(dabs(f_k(i,j) - f(i,j)),max_err)
             enddo
          enddo

          write(*,*)'--max:', max_err 

       enddo
       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       print*, '--Time it takes (s)', time_s
       end

       subroutine arrayf_k(nx,ny,f,f_k) 

       implicit none

       integer :: i,j
       integer :: nx,ny
       double precision :: f_k(nx,ny),f(nx,ny)

       do j=1,ny
            do i=1,nx

               f_k(i,j) = 0.25*f(i,j)

             enddo
          enddo

       end 

