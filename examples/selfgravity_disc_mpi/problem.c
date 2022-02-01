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

#define DT 0.005

long double* coeffs;

long double RANDOM_COEFFICIENT;

//load in the polynomials for the force
int load_force_polynomial(long double coeffs[], char* filename, int n_terms){
	FILE *fp;
	char buffer[255];
	fp = fopen(filename, "r");
	char* endptr;
	for(int i = 0; i < n_terms; i++){
		fscanf(fp, "%s", buffer);
		coeffs[i] = strtold(buffer, &endptr);
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
int load_initial_data(char* filename, long double positions[], int length){
	FILE *fp;
    char buffer[255];
    long double r;
    char* endptr;
	fp = fopen(filename, "r");
    if(fp == NULL){
        printf("die");
        return 0;
    }
	for(int i = 0; i < length; i++){
        fscanf(fp, "%s", buffer);
        r = strtold(buffer, &endptr);
        positions[i] = r;
    }
    fclose(fp);
    return 1;
}	

//initialize positions of particles
int generate_initial_positions(long double initial_data[], long double* positions[], int particles, int** mesh, int initial_data_length, long double radius){
    long double r, theta, x, y;
    int index;
    for(int i = 0; i < particles; i++){
        index = rand() % initial_data_length;
        r = initial_data[index];
        theta = ((long double) rand()) / RAND_MAX * 2 * M_PI;
        x = r * cos(theta);
        y = r * sin(theta);
        do {
            int a = 1500 + ((int) (x / radius));
            int b = 1500 + ((int) (y / radius)) ;
            if (a>=3000 || a<0 || b>=3000 || b<0)
            {
                printf ("segmentation fault %d %d %Lf %Lf\n", a, b, r, theta) ;
                break;
            }
            if (mesh[a][b])
            {}
            else{
                break;
            } 
            theta = rand() / RAND_MAX * 2 * M_PI;
            x = r * cos(theta);
            y = r * sin(theta);
        } while (1) ;
        
        positions[i][0] = x;
        positions[i][1] = y;
    }
    return 0;
}

//rng for brownian motion
long double get_random(long double mu, long double sigma, long double coefficient)
{
    long double U1, U2, W, mult;
    long double X1;

    do
    {
        U1 = -1 + ((long double)rand() / RAND_MAX) * 2;
        U2 = -1 + ((long double)rand() / RAND_MAX) * 2;
        W = powl(U1, 2) + powl(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    return coefficient * (mu + sigma * (long double)X1);
}

void additional_forces(struct reb_simulation* const r){
    long double VISCOSITY = 0.0010518; //dynamic viscosity of water
    long double RADIUS = 75e-6; //radius of particle
    long double CD = 6 * M_PI * VISCOSITY * RADIUS; //stokes drag
    long double DENSITY = 2260.0; //density of particles
    long double MASS = (4.0 / 3) * DENSITY * M_PI * pow(RADIUS, 3);

    // Simplest velocity dependent drag force.
    struct reb_particle* const particles = r->particles;
    const int N = r->N;
    for (int i=0;i<N;i++){

        long double r = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
        particles[i].ax += (double) ((get_random(0, sqrt(DT), RANDOM_COEFFICIENT) + (force(r, coeffs, 26) * particles[i].x / r) - (CD*particles[i].vx)) / MASS);
        particles[i].ay += (double) ((get_random(0, sqrt(DT), RANDOM_COEFFICIENT) + (force(r, coeffs, 26) * particles[i].y / r) - (CD*particles[i].vy)) / MASS);
        particles[i].az += (double) (-CD*particles[i].vz / MASS);
    }
}


void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
	char* FORCE_COEFFS_FILENAME = "forcecoeffs.csv"; //coefficients for force polynomial
    int TERMS = 26; //number of terms in force polynomial
    int INITIAL_DATA_LENGTH = 399460; //number of initial r values
    char* INITIAL_DATA_FILENAME = "initial_data.csv"; //initial r values
    int PARTICLES = 20000; //number of particles to simulate
    int MESH_FINENESS = 3000; //dimensions of mesh (MESH_FINENESS * MESH_FINENESS)

    long double VISCOSITY = 0.0010518; //dynamic viscosity of water
    long double RADIUS = 75e-6; //radius of particle
    long double DENSITY = 2260.0; //density of particles
    long double MASS = (4.0 / 3) * DENSITY * M_PI * pow(RADIUS, 3);
    long double CD = 6 * M_PI * VISCOSITY * RADIUS; //stokes drag
    long double TEMPERATURE = 28 + 273.15; //temperature
    long double KB = 1.38064852e-23; //boltzmann's constant
    RANDOM_COEFFICIENT = sqrt(2 * KB * TEMPERATURE / CD); //coefficient in front of the dW term
    double COEFF = 17589.3;
    int T_END = 1000;
    printf("test\n");

    coeffs = (long double*) malloc(TERMS * sizeof(long double));
	load_force_polynomial(coeffs, FORCE_COEFFS_FILENAME, 26);

    long double *initial_positions = (long double*) malloc(INITIAL_DATA_LENGTH * sizeof(long double));
    load_initial_data(INITIAL_DATA_FILENAME, initial_positions, INITIAL_DATA_LENGTH);

    //initialize 2d array for positions
    long double** positions = (long double**) malloc(PARTICLES * sizeof(long double*));
    for(int i = 0; i < PARTICLES; i++){
        positions[i] = (long double*) malloc(2 * sizeof(long double));
    }
    //initialize 2d array for collision detection
    int** mesh = (int**) malloc(MESH_FINENESS * sizeof(int*));
    for(int i = 0; i < MESH_FINENESS; i++){
        mesh[i] = (int*) malloc(MESH_FINENESS * sizeof(int));
    }
    //initialize to false
    for (int i = 0; i < MESH_FINENESS; i++) {
        for (int j = 0; j < MESH_FINENESS; j++) {
            mesh[i][j] = 0;
        }
    }
    //initialize positions
    generate_initial_positions(initial_positions, positions, PARTICLES, mesh, INITIAL_DATA_LENGTH, RADIUS);
    free(initial_positions); // free up memory since its not needed any more
    free(mesh);

    struct reb_simulation* const r = reb_create_simulation();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->gravity    = REB_GRAVITY_TREE;
    r->boundary    = REB_BOUNDARY_OPEN;
    r->opening_angle2    = 1.5;        // This constant determines the accuracy of the tree code gravity estimate.
    r->G         = COEFF;        
    r->softening     = 0.0075163;        // Gravitational softening length
    r->dt         = DT;        // Timestep
    r->collision = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere; 
    r->coefficient_of_restitution = 0;
    r->additional_forces = additional_forces;
    const double boxsize = 0.10;
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
    double disc_mass = MASS * PARTICLES/r->mpi_num;    // Total disc mass
    int N = PARTICLES/r->mpi_num;            // Number of particles
    // Initial conditions
    
    for (int i=0;i<N;i++){
        struct reb_particle pt = {0};
        pt.x         = (double) positions[i][0];
        pt.y         = (double) positions[i][1];
        pt.z         = 0.0;
        pt.vx         =  0.0;
        pt.vy         =  0.0;
        pt.vz         = 0.0;
        pt.ax = 0;
        pt.ay = 0;
        pt.az = 0;
        pt.r = 0.001;
        pt.m         = disc_mass/(double)N;
        reb_add(r, pt);
    }
    printf("death\n");
    r->heartbeat = heartbeat;


#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->allocatedN *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->allocatedN);
#endif // OPENGL
    
    // Start the integration
    reb_integrate(r, T_END);

    // Cleanup
    reb_mpi_finalize(r);
    reb_free_simulation(r); 
}

void heartbeat(struct reb_simulation* const r){
    if (reb_output_check(r,10.0*r->dt)){
        reb_output_timing(r,0);
    }
}
