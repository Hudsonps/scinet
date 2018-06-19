//assign2functions
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <rarray>
#include <mpi.h>
//include </Users/hudsonps/Dropbox/scinet/assign7/rarray/rarray>
//#include <omp.h>

int number_elements_percore(int npnts, int rank, int size)
{
  int remainder = npnts%size;
  int p;
  p = npnts/size;

  if(rank < remainder)
    p = p+1; //if npnts/size is not an integer, we must distribute the remaining elements to the cores
  
  if(rank == 0)
    p = p+1;
  if(rank == size-1)
    p = p+1; //These if's introduce the sites that account for the boundary conditions
  return p;
}			     
			    
void read_data(std::ifstream& infile, double &c, double &tau, double &x1, double &x2, double &runtime, double &dx, double &outtime, std::string &outfilename){
  //Physical parameters 
  infile >> c;  //wave speed 
  infile >> tau; //damping time
  infile >> x1; //left most x value
  infile >> x2; //right most x value

  //Simulation parameters
  infile >> runtime; //how long should the simulation try to compute
  infile >> dx; //spatial grid size

  //Output parameters
  infile >> outtime; //how often snapshots of the wave are written out
  infile >> outfilename; //name of file with output data
  
  infile.close(); //close the file where the parameters are stored
}


void derived_parameters(int& ngrid, int& npnts, double& dt, int& nsteps, int& nper, double c, double x1, double x2, double runtime, double dx, double outtime){
  //This funcion defines a set of additional parameters 
    ngrid   = (x2-x1)/dx;  // number of x points
    npnts   = ngrid + 2;   // number of x points including boundary points
    dt      = 0.5*dx/c;    // time step size
    nsteps  = runtime/dt;  // number of steps of that size to reach runtime
    nper    = outtime/dt;  // how many step s between snapshots
}

void report_values(std::ofstream& fout, double c, double tau, double x1, double x2, double runtime, double dx, double outtime, int ngrid, double dt, int nsteps, int nper){
      // Report all the values on the outputfile indicated by fout
    fout << "#c       " << c       << std::endl;
    fout << "#tau     " << tau     << std::endl;
    fout << "#x1      " << x1      << std::endl;
    fout << "#x2      " << x2      << std::endl;
    fout << "#runtime " << runtime << std::endl;
    fout << "#dx      " << dx      << std::endl;
    fout << "#outtime " << outtime << std::endl; 
    fout << "#ngrid   " << ngrid   << std::endl;
    fout << "#dt      " << dt      << std::endl;
    fout << "#nsteps  " << nsteps  << std::endl;    
    fout << "#nper    " << nper    << std::endl;
}

void initialize(rarray<double, 1> x, rarray<double, 1> rho_prev, rarray<double, 1> rho, rarray<double, 1> rho_next, double x1, double x2, int ngrid, int npnts){
    //Initialize simply initializes our vectors.
    //x represents the position axis. This function discretizes the axis by creating a grid determined by x1,x2,ngrid,npnts
    //All the rho's are simply set to 0.
  //The three vectors rho represent the shape of the pulse for three consecutive time steps, which will be necessary to solve a PDE
  for (int i = 0; i < npnts; i++) {
    x[i] = x1 + ((i-1)*(x2-x1))/ngrid; 
    rho[i] = 0.0;
    rho_prev[i] = 0.0;
    rho_next[i] = 0.0;
    }
}

void excite(rarray<double, 1> x, rarray<double, 1> rho_prev, rarray<double, 1> rho, rarray<double, 1> rho_next, double x1, double x2, int npnts){
    //This function determines the initial shape of rho
    //rho_prev is made equal to rho at this level
    double x14 = 0.25*(x2-x1) + x1;
    double x12 = 0.5*(x2+x1);
    double x34 = 0.75*(x2-x1) + x1;
    for (int i = npnts/4 + 1; i < 3*npnts/4; i++) {
        if (x[i]<x14 || x[i]>x34)
            rho[i] = 0.0;
        else
            rho[i] = 0.25 - fabs(x[i]-x12)/(x2-x1);
        rho_prev[i] = rho[i];
    }
}

void write_snapshots(std::ofstream& fout, int s, int nper, int ngrid, rarray <double, 1> x, rarray<double , 1> rho, double dt)
{
  //nper = after how many steps a snapshot should be printed
  //s = tracks the current iteration number
  if ((s+1)%nper == 0) {
            fout << "\n\n# t = " << s*dt << "\n";
            for (int i = 1; i <= ngrid; i++) 
                fout << x[i] << " " << rho[i] << "\n";
  }
}

void write_snapshots_binary(FILE *fout, int s, int nper, int npnts, rarray <double, 1> x, rarray<double , 1> rho, double dt)
{
  //nper = after how many steps a snapshot should be printed
  //s = tracks the current iteration number
  double xaux[npnts];
  double rhoaux[npnts];
  for (int i = 0; i < npnts; i++) {
    xaux[i] = x[i];
    rhoaux[i] = rho[i];
  }
  
  if ((s+1)%nper == 0) {
    fwrite(xaux, sizeof(double), npnts, fout);
    fwrite(rhoaux, sizeof(double), npnts, fout);
  }
}

void evolve(rarray<double, 1> x, double dx, rarray<double, 1> rho_prev, rarray<double, 1> rho, rarray<double, 1> rho_next, double c, double tau, double dt, int ngrid, int nsteps, int nper, int npnts){
      // Evolve
    for (int i = 1; i <= ngrid; i++) {
      double laplacian = pow(c/dx,2)*(rho[i+1] + rho[i-1] - 2*rho[i]);
      double friction = (rho[i] - rho_prev[i])/tau;
      rho_next[i] = 2*rho[i] - rho_prev[i] + dt*(laplacian*dt-friction);
    }
}
  
void rotate_rho(rarray <double, 1> rho_prev, rarray <double, 1> rho, rarray <double, 1> rho_next, int ngrid){
      for (int i = 1; i <= ngrid; i++) {
      double temp = rho_prev[i];
      rho_prev[i] = rho[i];
      rho[i]      = rho_next[i];
      rho_next[i] = temp;
    }
}


void shapevstime(std::ofstream& fout, rarray<double, 1> x, double dx, rarray<double, 1> rho_prev, rarray<double, 1> rho, rarray<double, 1> rho_next, double c, double x1, double x2, double tau, double dt, int ngrid, int nsteps, int nper, int npnts){
  
  //Output in binary
  FILE *out;
  out = fopen("binary.txt", "wb");  //will be used for the output in binary.
  
    // Output initial signal
  //fout << "\n#t = " << 0 << "\n";
  //for (int i = 1; i <= ngrid; i++) 
  //fout << x[i] << " " << rho[i] << "\n";

  for (int s = 0; s < nsteps; s++) {
    // Set zero dirichlet boundary conditions
    rho[0] = 0.0;
    rho[ngrid+1] = 0.0;

    //evolve_parallel(x, dx, rho_prev, rho, rho_next, c, tau, dt, ngrid, nsteps, nper, npnts);

    //rotate_rho_parallel(rho_prev, rho, rho_next, ngrid);
    
        // Output density
    //write_snapshots(fout, s, nper, ngrid, x, rho, dt);
    //write_snapshots_binary(out, s, nper, npnts, x, rho, dt);
  }
  fclose(out);
}


		   

			
