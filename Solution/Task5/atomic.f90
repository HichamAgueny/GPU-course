      program atomic

      implicit none

       integer :: i,j,k
       integer :: iter,count_rate, count_max,count
       integer :: t_start,t_final
       integer, parameter :: nx=8192
       double precision, parameter    :: pi=4d0*dtan(1d0) 
       double precision               :: total,time_s,sum_x
       double precision, allocatable  :: u(:),v(:),w(:)

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)
   
       allocate(u(nx)); allocate(v(nx)); allocate(w(nx))

       do i=1,nx
          u(i) = dsin(i*1d0)**2
          v(i) = dcos(i*1d0)**2 
       enddo
     
!offload data and eliminate the dependencies 
!$acc parallel loop
      do i=1,nx
!$acc atomic update
         u(i) = u(i) + v(i)
      enddo
!$acc end parallel

!sum up the vector elements of w
      sum_x = 0d0
!$acc parallel loop reduction(+:sum_x)
      do i=1,nx
         sum_x = sum_x + u(i)
      enddo
!$acc end parallel

      write(*,*)'--sum w/nx',sum_x/nx
                   
       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       print*, '--Time it takes (s)', time_s
       end

