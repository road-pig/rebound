/**
 * Self-gravitating disc with MPI
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator. Collisions are not resolved.
 * This program makes use of MPI. Note that you need 
 * to have MPI compilers (mpicc) installed. The code is using 
 * four root boxes to distribute to particles to one, two 
 * or four MPI nodes. How to efficiently run this code on
 * large clusters goes beyond this simple example and
 * almost certainly requires experimentation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"


//load in the polynomials for the force
int load_force_polynomial(long double coeffs[], char* filename, int n_terms){
	FILE *fp;
	char buffer[255];
	fp = fopen(filename, "r");
	char* endptr;
	for(int i = 0; i < n_terms; i++){
		fscanf(fp, "%s", buffer);
		coeffs[i] = strtold(buffer, &endptr);
		printf("%Lf", coeffs[i]);
	}
	fclose(fp);
	return 0;
}

//force due to the magnet
long double force(long double r, const long double coeffs[], int terms){
    long double sum = 0;
    for (int i = 0; i < terms; i++){
        sum += coeffs[i] * powl(r, i);
    }
    return sum;
}

//load initial distributions
int load_initial_data(char* filename, int length, double positions[]){
	FILE *fp;
	char buffer[255];
	fp = fopen(filename, "r");
	char* endptr;
	for(int i = 0
}	

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
	const char* FORCE_COEFFS_FILENAME = "forcecoeffs.csv"; //coefficients for force polynomial
    const int TERMS = 26; //number of terms in force polynomial
    const int INITIAL_DATA_LENGTH = 399460; //number of initial r values
    const char* INITIAL_DATA_FILENAME = "initial_data.csv"; //initial r values
    const int PARTICLES = 70; //number of particles to simulate
    const int MESH_FINENESS = 3000; //dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
    const int M = 100000; //number of timesteps
    const int N_THREADS = 6;

    const long double VISCOSITY = 0.0010518; //dynamic viscosity of water
    const long double RADIUS = 2e-3; //radius of particle
    const long double DENSITY = 2260; //density of particles
    const long double MASS = (4.0 / 3) * DENSITY * M_PI * pow(RADIUS, 3);
    const long double CD = 6 * M_PI * VISCOSITY * RADIUS; //stokes drag
    const long double TEMPERATURE = 28 + 273.15; //temperature
    const long double KB = 1.38064852e-23; //boltzmann's constant
    const long double RANDOM_COEFFICIENT = sqrt(2 * KB * TEMPERATURE / CD); //coefficient in front of the dW term
    const int T_START = 0;
    const int T_END = 1000;

	long double *coeffs = (long double*) malloc(TERMS * sizeof(long double));
	load_force_polynomial(coeffs, FORCE_COEFFS_FILENAME, 26);

    struct reb_simulation* const r = reb_create_simulation();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->gravity    = REB_GRAVITY_TREE;
    r->boundary    = REB_BOUNDARY_OPEN;
    r->opening_angle2    = 1.5;        // This constant determines the accuracy of the tree code gravity estimate.
    r->G         = 1;        
    r->softening     = 0.02;        // Gravitational softening length
    r->dt         = 3e-2;        // Timestep
    const double boxsize = 10.2;
    // Setup root boxes for gravity tree.
    // Here, we use 2x2=4 root boxes (each with length 'boxsize')
    // This allows you to use up to 4 MPI nodes.
    reb_configure_box(r,boxsize,2,2,1);

    // Initialize MPI
    // This can only be done after reb_configure_box.
    reb_mpi_init(r);

    // Setup particles only on master node
    // In the first timestep, the master node will 
    // distribute particles to other nodes. 
    // Note that this is not the most efficient method
    // for very large particle numbers.
    double disc_mass = 2e-1/r->mpi_num;    // Total disc mass
    int N = 10000/r->mpi_num;            // Number of particles
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
    if (r->mpi_id==0){
        reb_add(r, star);
    }
    for (int i=0;i<N;i++){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(r, boxsize/10.,boxsize/2./1.2,-1.5);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.001);
        double mu     = star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
        double vkep     = sqrt(r->G*mu/a);
        pt.vx         =  vkep * sin(phi);
        pt.vy         = -vkep * cos(phi);
        pt.vz         = 0;
        pt.m         = disc_mass/(double)N;
        reb_add(r, pt);
    }
    r->heartbeat = heartbeat;

#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->allocatedN *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->allocatedN);
#endif // OPENGL
    
    // Start the integration
    reb_integrate(r, INFINITY);

    // Cleanup
    reb_mpi_finalize(r);
    reb_free_simulation(r); 
}

void heartbeat(struct reb_simulation* const r){
    if (reb_output_check(r,10.0*r->dt)){
        reb_output_timing(r,0);
    }
}
