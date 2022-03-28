       program loop

       integer :: i,j
       integer :: iter,count_rate, count_max,count
       integer :: t_start,t_final

       integer, parameter :: n=1024,m=1024
       double precision, allocatable   :: a(:,:),b(:,:),c(:,:)

       allocate(a(n,m)); allocate(b(n,m)); allocate(c(n,m))

!generate random vectors
       CALL RANDOM_NUMBER(b)
       CALL RANDOM_NUMBER(c)

       call system_clock(count_max=count_max, count_rate=count_rate)

       call system_clock(t_start)

!parallelise the following loop
       do i=1,n
         do j=1,m
            a(i,j) = b(i,j) + c(i,j)
         enddo
       enddo

       call system_clock(t_final)

       time_s = real(t_final - t_start)/real(count_rate)

       print*, '--Time it takes (s)', time_s
        write(*,*)'--it is done: sum', sum(a)
       end
