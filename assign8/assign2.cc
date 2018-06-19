//Simulates a one-dimensional damped wave equation
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "assign2functions.h"
#include <rarray>
#include <mpi.h>
#include <string>
//include </Users/hudsonps/Dropbox/scinet/assign7/rarray/rarray>
//#include <omp.h>

int main(int argc, char **argv)
{

  int rank, size;
  int ierr;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);

  //The first part of this function should be run in parallel. We want all the processors to have access to the parameters of our problem
  
    std::ifstream infile("parameters.txt"); //where to read the data from to fill the parameters below
    // Physical parameters, 
    double  c;  // wave speed
    double  tau;  // damping time
    double  x1;  // left most x value
    double  x2;  // right most x value

    // Simulation parameters
    double  runtime;  // how long should the simulation try to compute?
    double  dx;   // spatial grid size

    // Output parameters
    double  outtime;   // how often should a snapshot of the wave be written out?
    std::string outfilename; // name of the file with the output data
    
    int ngrid; //number of x points
    int npnts; //number of x points including boundary points
    double dt; //time step size
    int nsteps; //number of steps of that size to reach runtime
    int nper;   //how many step s between snapshots

    //All the variables above will get inputs from the function read_data, which accesses the file indicated by infile
    read_data(infile, c, tau, x1, x2, runtime, dx, outtime, outfilename);
    
    //Attributes proper values to these parameters, which are not independent from the previous ones
    derived_parameters(ngrid, npnts, dt, nsteps, nper, c, x1, x2, runtime, dx, outtime);       

    //Open output file
    std::ofstream fout("myresults" + std::to_string(rank) + ".txt");
    //With the line above, there will be one file, myresult[rank].txt for each processor


    // Report all the values on fout by writting them on the outputfile (indicated by fout)
    report_values(fout, c, tau, x1, x2, runtime, dx, outtime, ngrid, dt, nsteps, nper);

    
    // Define and allocate arrays
    //double* rho_prev = new double[npnts]; // time step t-1
    //double* rho      = new double[npnts]; // time step t
    //double* rho_next = new double[npnts]; // time step t+1
    //double* x        = new double[npnts]; // x values
    //Our functions will have these arrays as inputs to solve the problem

    //previously these vectors worked belonged to a single processor, but now we must split them.
    //Here is how we split them.

    //size is the number of processors being used.
    //npnts is the number of points including the boundary conditions
    //ngrid is the number of points without the boundary conditions
    //In general, ngrid/npnts is not an integer. The function below takes that into consideration to assign sites to each core

    int npnts_percore = number_elements_percore(npnts, rank, size);
    std::cout << npnts_percore << " " << ngrid << std::endl;
    
    rarray<double,1> rho_prev(npnts_percore);
    rarray<double,1> rho(npnts_percore);
    rarray<double,1> rho_next(npnts_percore);
    rarray<double, 1> x(npnts_percore);
    //Now we have rarrays split between our different cores. 

    

    //initialize(x, rho_prev, rho, rho_next, x1, x2, ngrid, npnts);
    //initialize simply initializes our vectors.
    //x represents the position axis. This function discretizes the axis by creating a grid determined by x1,x2,ngrid,npnts
    //All the rho's are simply set to 0.

    // excite(x, rho_prev, rho, rho_next, x1, x2, npnts);
    //This function determines the initial shape of the pulse
    //rho_prev is made equal to rho at this level

    //shapevstime(fout, x, dx, rho_prev, rho, rho_next, c, x1, x2, tau, dt, ngrid, nsteps, nper, npnts);
    //This function writes lists with rho as a function of x for different instants of time
    //It does so by solving the PDE of a one-dimensional wave with damping
    
    // Close file
  
    //fout.close();
    //std::cout << "Results written to "<< outfilename << std::endl;
    
    // Deallocate memory
    //  delete[] rho_prev;
    // delete[] rho;
    //delete[] rho_next;
    //delete[] x;
  
    ierr = MPI_Finalize();

    return 0;
}
