/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

const double k = 1.;
const double m = 1.;

void derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    const double omega = sqrt(k/m);
    yDot[0] = y[1]; 
    yDot[1] = -omega*omega*y[0];
}

double energy_ho(struct reb_ode* const ode){
    double x = ode->y[0];
    double v = ode->y[1];
    return 1./2.*k*x*x + 1./2.*m*v*v;
}

double run(double e, int integrator, double* energyerror, double* energyerror_ho){
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m", 1.);                // Central object
    reb_add_fmt(r, "m a e", 1e-3, 1., e);   // Jupiter mass planet
    reb_add_fmt(r, "m a e omega", 1e-3, 4.312, e, 4.6);   // Jupiter mass planet
    reb_move_to_com(r);
    double E0 = reb_tools_energy(r);

    r->integrator = integrator;
    r->ri_bs.eps_rel = 1e-8;
    r->ri_bs.eps_abs = 1e-8;
    r->dt = 1e-2;

    struct reb_ode* ho = reb_create_ode(r,2);   // Add an ODE with 2 dimensions
    ho->derivatives = derivatives;              // Right hand side of the ODE
    ho->y[0] = 1;                               // Initial conditions
    ho->y[1] = 0;
    double E0_ho = energy_ho(ho);

    struct timeval start; gettimeofday(&start, NULL);

    reb_integrate(r, 10.);
    
    struct timeval stop; gettimeofday(&stop, NULL);
    
    double E1 = reb_tools_energy(r);
    *energyerror = fabs((E1-E0)/E0);

    double E1_ho = energy_ho(ho);
    *energyerror_ho = fabs((E1_ho-E0_ho)/E0_ho);
    printf("\nt=%.19f\n",r->t);

    reb_free_ode(ho);
    reb_free_simulation(r);

    return stop.tv_sec - start.tv_sec + (stop.tv_usec-start.tv_usec)/1e6;

}
int main(int argc, char* argv[]) {
    double energy_BS, energy_IAS15, energy_WHFAST;
    double energy_ho_BS, energy_ho_IAS15, energy_ho_WHFAST;
    double e = 0.92;
    double runtime_BS = run(e, REB_INTEGRATOR_BS, &energy_BS, &energy_ho_BS);
    double runtime_IAS15 = run(e, REB_INTEGRATOR_IAS15, &energy_IAS15, &energy_ho_IAS15);
    double runtime_WHFAST = run(e, REB_INTEGRATOR_WHFAST, &energy_WHFAST, &energy_ho_WHFAST);
    printf("BS (MP):   %.5fs \t %.5e (nbody) \t %.5e (ho)\n",runtime_BS,   energy_BS,   energy_ho_BS  );
    printf("IAS15:     %.5fs \t %.5e (nbody) \t %.5e (ho)\n",runtime_IAS15, energy_IAS15, energy_ho_IAS15);
    printf("WHFAST:    %.5fs \t %.5e (nbody) \t %.5e (ho)\n",runtime_WHFAST, energy_WHFAST, energy_ho_WHFAST);

}

