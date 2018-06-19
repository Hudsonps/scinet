#include <mpi.h>

void teste(double *x, int dim)
{
  for(int i = 0; i < dim; i++)
    x[i] = i;
}


int main(int argc, char **argv) {
int rank, size;
int ierr;

/*double a[10];
 teste(a, 10);

 teste(a, 10);
  
 for(int i = 0; i <=9; i++)
 std::cout << a[i] << std::endl;*/

 
ierr = MPI_Init(&argc, &argv);
ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
int p = 10/size;

 int c = 1;

 if(rank == 2)
   c = 2;
 
 std::cout<< c << " " << rank << std::endl;


 
ierr = MPI_Finalize();


 return 0;
}
