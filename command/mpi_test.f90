! mpirun -np 2 ./test
program mpi_test
    use mpi, only: mpi_init, mpi_comm_rank, mpi_comm_size, mpi_finalize, MPI_COMM_WORLD
    implicit none

    integer :: ierror, rank, np

    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)
    call mpi_comm_size(MPI_COMM_WORLD, np, ierror)
    write(*,*) "Process and np:", rank, np
!    call mpi_finalize(ierror)

end program mpi_test
