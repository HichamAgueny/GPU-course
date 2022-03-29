      program laplace_acc

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

!generate a random vector to define the initial conditions

       do i=1,nx
          f(i,1)  = dsin(pi*(i-1)/(nx-1))
          f(i,ny) =0d0
       enddo
       do j=1,ny
          f(1,j) = dsin(pi*(j-1)/(ny-1))
          f(nx,j) = 0d0
       enddo

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)
      
       print*, "--Start the iteration"
       iter = 0; max_err=1.0

!$acc data copyin(f) copyout(f_k)
       do while (max_err.gt.error.and.iter.le.max_iter)

!$acc parallel loop collapse(2)
         do j=2,ny-1
            do i=2,nx-1
               d2fx = f(i+1,j) + f(i-1,j)
               d2fy = f(i,j+1) + f(i,j-1)

               f_k(i,j) = 0.25*(d2fx + d2fy)
             enddo
          enddo
!$acc end parallel

          max_err=0.

!$acc parallel loop reduction(max:max_err) collapse(2)
          do j=2,ny-1
            do i=2,nx-1
               max_err = max(dabs(f_k(i,j) - f(i,j)),max_err)
               f(i,j) = f_k(i,j)
            enddo
          enddo
!$acc end parallel

          if(mod(iter,20).eq.0 ) write(*,'(i5,f10.6)') iter,max_err 
          iter = iter +1 

        enddo
!$acc end data
       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       print*, '--Time it takes (s)', time_s
       print*,"--Sum",sum(f_k(:,:))
       end

