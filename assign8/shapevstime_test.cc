#include <iostream>
#include <fstream>
#include <cmath>
#include "assign2functions.h"
#include <rarray>

int main()
{
    std::ifstream infile("parameters.txt"); //where to read the data from to fill the parameters below
// Physical parameters
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

        // Define and allocate arrays
    rarray<double,1> rho_prev(npnts);
    rarray<double,1> rho(npnts);
    rarray<double,1> rho_next(npnts);
    rarray<double, 1> x(npnts);

    std::ofstream fout(outfilename);

    // Initialize
    for (int i = 0; i < npnts; i++) {
        x[i] = x1 + ((i-1)*(x2-x1))/ngrid; 
        rho[i] = 0.0;
        rho_prev[i] = 0.0;
        rho_next[i] = 0.0;
    }

    shapevstime(fout, x, dx, rho_prev, rho, rho_next, c, x1, x2, tau, dt, ngrid, nsteps, nper);

    double eps = 1e-10;

    int flag = 0;
    for (int i = 0; i < npnts; i++) {
      if(fabs(rho[i] - 0) > eps)
	{
	  flag = 1;
	  break;
	}
    }

    if(flag == 0){
      std::cout << "The test of the function shapevstime was successful." << std::endl;
    }
    else{
      std::cout << "The test of the function shapevstime failed." << std::endl;
    }
    
    
}
