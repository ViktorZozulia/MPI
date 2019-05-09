module Task

    use mpi
    implicit none
    contains

        subroutine GetMaxCoordinates(A, x1, y1, x2, y2)
        implicit none
        real(8), dimension(:,:), intent(in) :: A
        integer(4), intent(out) :: x1, y1, x2, y2

        integer(4) :: n, L, R, Up, Down, m, tmp
        real(8), allocatable :: current_column(:), B(:,:)
        real(8) :: current_sum, max_sum, max_sum_for_all
        logical :: transpos

        integer(4) :: mpiErr,mpiRank, mpiSize, max_mpiRank, test_max_mpiRank

        m = size(A, dim=1) 
        n = size(A, dim=2) 
        transpos = .FALSE.

        if (m < n) then 
            transpos = .TRUE.   
            B = transpose(A)
            m = size(B, dim=1) 
            n = size(B, dim=2) 
        else
            B = A     
        endif

        allocate(current_column(m))


        max_sum=B(1,1)
        x1=1
        y1=1
        x2=1
        y2=1

        
        

        call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
     
         
        L=mpiRank+1        

        do while (L<=n)
          
        
            current_column = B(:, L)       
            do R=L,n
 
                if (R > L) then 
                    current_column = current_column + B(:, R)
                endif
                
                call FindMaxInArray(current_column, current_sum, Up, Down) 

                if (current_sum > max_sum) then
                    max_sum = current_sum
                    x1 = Up
                    x2 = Down
                    y1 = L
                    y2 = R
                endif
       
            end do

            L=L+mpiSize
 
        end do
  

        call mpi_reduce(max_sum, max_sum_for_all, 1, MPI_INTEGER4, MPI_MAX, 0, MPI_COMM_WORLD, mpiErr)  
        call mpi_bcast(max_sum_for_all, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr) 

        test_max_mpiRank = 0
        if (max_sum_for_all == max_sum) then 
           test_max_mpiRank = mpiRank
        end if
  
        call mpi_reduce(test_max_mpiRank, max_mpiRank, 1, MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(max_mpiRank, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)          

        call mpi_bcast(x1, 1, MPI_INTEGER4, max_mpiRank, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(x2, 1, MPI_INTEGER4, max_mpiRank, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(y1, 1, MPI_INTEGER4, max_mpiRank, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(y2, 1, MPI_INTEGER4, max_mpiRank, MPI_COMM_WORLD, mpiErr)

        
        deallocate(current_column)


        if (transpos) then  
            tmp = x1
            x1 = y1
            y1 = tmp
    
            tmp = y2
            y2 = x2
            x2 = tmp
        endif

        end subroutine


        subroutine FindMaxInArray(a, Sum, Up, Down)
            real(8), intent(in), dimension(:) :: a
            integer(4), intent(out) :: Up, Down
            real(8), intent(out) :: Sum
            real(8) :: cur_sum
            integer(4) :: minus_pos, i

            Sum = a(1)
            Up = 1
            Down = 1
            cur_sum = 0
            minus_pos = 0


            do i=1, size(a)
                
                cur_sum = cur_sum + a(i)
               
                if (cur_sum > Sum) then
                    Sum = cur_sum
                    Up = minus_pos + 1
                    Down = i
                endif
         
                if (cur_sum < 0) then
                    cur_sum = 0
                    minus_pos = i
                endif

            enddo

        end subroutine FindMaxInArray


end module Task


