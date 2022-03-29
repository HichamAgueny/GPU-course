      program reduction

      implicit none

       integer :: i,j,k
       integer :: iter,count_rate, count_max,count
       integer :: t_start,t_final
       integer, parameter :: nx=8192,ny=nx
       double precision, parameter    :: pi=4d0*dtan(1d0) 
       double precision               :: max_err,time_s,&
                                         d2fx,d2fy,sum_x,&
                                         total
       double precision, allocatable  :: f(:,:), f_k(:,:)

       allocate(f(nx,ny)); allocate(f_k(nx,ny))

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)
     
!generate a random vector

       CALL RANDOM_NUMBER(f)

       max_err=0.0

!parallelise the loop and compute the max using an openacc directive
!$acc parallel loop reduction(max:max_err) collapse(2)
         do j=1,ny
            do i=1,nx

               f_k(i,j) = 0.25*f(i,j)

               max_err = max(dabs(f_k(i,j) - f(i,j)),max_err)
             enddo
          enddo
!$acc end parallel

          write(*,*)'--max:', max_err 

!parallelise the loop and compute the sum using an openacc directive
!1D
       total = 0d0
!$acc parallel loop reduction(+:total)
       do i=1,nx
          total = total + f(i,1)
       enddo
!$acc end parallel
       write(*,*)'--1D sum:', total/nx

!2D
       total = 0d0
!$acc parallel loop reduction(+:total,sum_x)
      do i=1,nx
        sum_x = 0d0
        do j=1,ny
           sum_x = sum_x + f(i,j)
        enddo
        total = total + sum_x
      enddo
!$acc end parallel

      write(*,*)'--2D sum:', total/(nx*ny)   
  
       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       print*, '--Time it takes (s)', time_s
       end
