/**
 * @file 	integrator.c
 * @brief 	BS integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the Gragg-Bulirsch-Stoer integration scheme.  
 *          It is a reimplementation of the fortran code by E. Hairer and G. Wanner.
 *          The starting point was the JAVA implementation in hipparchus:
 *          https://github.com/Hipparchus-Math/hipparchus/blob/master/hipparchus-ode/src/main/java/org/hipparchus/ode/nonstiff/GraggBulirschStoerIntegrator.java
 *
 * @section 	LICENSE
 * Copyright (c) 2021 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (c) 2004, Ernst Hairer
 *
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following
 * conditions are met:
 * 
 *  - Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <float.h> // for DBL_MAX
#include "rebound.h"
#include "integrator_bs.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
    
// Default configuration parameter. 
// They are hard coded here because it
// is unlikely that these need to be changed by the user.
static const int maxOrder = 18;// was 18 
static const int sequence_length = maxOrder / 2; 
static const double stepControl1 = 0.65;
static const double stepControl2 = 0.94;
static const double stepControl3 = 0.02;
static const double stepControl4 = 4.0;
static const double orderControl1 = 0.8;
static const double orderControl2 = 0.9;
static const double stabilityReduction = 0.5;
static const int maxIter = 2; // maximal number of iterations for which checks are performed
static const int maxChecks = 1; // maximal number of checks for each iteration


static int tryStep(struct reb_ode_state* state, const int k, const int n, const double t0, const double step, const int method) {
    const double subStep  = step / n;
    const int length = state->length;
    double t = t0;
    const double* y0 = state->y;
    const double* scale = state->scale;
    double* const y0Dot = state->y0Dot;
    double* const y1 = state->y1;

    switch (method) {
        case 0: // LeapFrog
            {
                // first substep
                for (int i = 0; i < length; ++i) {
                    if (i%6<3){ // Drift
                        y1[i] = y0[i] + 0.5*subStep * y0[i+3];
                    }
                }
                t += 0.5*subStep;
                state->derivatives(state, y0Dot, y1, t);
                for (int i = 0; i < length; ++i) {
                    if (i%6>2){ // Kick
                        y1[i] = y0[i] + subStep * y0Dot[i];
                    }
                }


                // other substeps
                for (int j = 1; j < n; ++j) {
                    t += subStep;
                    for (int i = 0; i < length; ++i) {
                        if (i%6<3){ // Drift
                            y1[i] = y1[i] + subStep * y1[i+3];
                        }
                    }
                    state->derivatives(state, y0Dot, y1, t);
                    for (int i = 0; i < length; ++i) {
                        if (i%6>2){ // Kick
                            y1[i] = y1[i] + subStep * y0Dot[i];
                        }
                    }

                    // stability checki // TODO
                    //if (performStabilityCheck && (j <= maxChecks) && (k < maxIter)) {
                    //    double initialNorm = 0.0;
                    //    for (int l = 0; l < length; ++l) {
                    //        const double ratio = y0Dot[l] / scale[l];
                    //        initialNorm += ratio * ratio;
                    //    }
                    //    double deltaNorm = 0.0;
                    //    for (int l = 0; l < length; ++l) {
                    //        const double ratio = (yDot[l] - y0Dot[l]) / scale[l];
                    //        deltaNorm += ratio * ratio;
                    //    }
                    //    //printf("iii   %e %e\n",initialNorm, deltaNorm);
                    //    if (deltaNorm > 4 * MAX(1.0e-15, initialNorm)) {
                    //        return 0;
                    //    }
                    //}
                }

                // correction of the last substep (at t0 + step)
                for (int i = 0; i < length; ++i) {
                    if (i%6<3){ // Drift
                        y1[i] = y1[i] + 0.5 * subStep * y1[i+3];
                    }
                }

                return 1;
            }
        case 1: // Modified Midpoint
            {
                // first substep
                t += subStep;
                for (int i = 0; i < length; ++i) {
                    y1[i] = y0[i] + subStep * y0Dot[i];
                }

                // other substeps
                double* const yTmp = state->yTmp; 
                double* const yDot = state->yDot; 

                state->derivatives(state, yDot, y1, t);
                for (int i = 0; i < length; ++i) {
                    yTmp[i] = y0[i];
                }

                for (int j = 1; j < n; ++j) {  // Note: iterating n substeps, not 2n substeps as in Eq. (9.13)
                    t += subStep;
                    for (int i = 0; i < length; ++i) {
                        const double middle = y1[i];
                        y1[i]       = yTmp[i] + 2.* subStep * yDot[i];
                        yTmp[i]       = middle;
                    }

                    state->derivatives(state, yDot, y1, t);

                    // stability check
                    if (j <= maxChecks && k < maxIter) {
                        double initialNorm = 0.0;
                        for (int l = 0; l < length; ++l) {
                            const double ratio = y0Dot[l] / scale[l];
                            initialNorm += ratio * ratio;
                        }
                        double deltaNorm = 0.0;
                        for (int l = 0; l < length; ++l) {
                            const double ratio = (yDot[l] - y0Dot[l]) / scale[l];
                            deltaNorm += ratio * ratio;
                        }
                        if (deltaNorm > 4 * MAX(1.0e-15, initialNorm)) {
                            return 0;
                        }
                    }

                }

                // correction of the last substep (at t0 + step)
                for (int i = 0; i < length; ++i) {
                    y1[i] = 0.5 * (yTmp[i] + y1[i] + subStep * yDot[i]); // = 0.25*(y_(2n-1) + 2*y_n(2) + y_(2n+1))     Eq (9.13c)
                }

                return 1;
            }
            return 0;
            break;
        default:
            printf("Error. method not implemented in BS\n");
            exit(1);
    }
}

static void extrapolate(const struct reb_ode_state* state, double * const coeff, const int k) {
    double* const y1 = state->y1;
    double* const C = state->C;
    double** const D =  state->D;
    double const length = state->length;
        for (int j = 0; j < k; ++j) {
        double xi = coeff[k-j-1];
        double xim1 = coeff[k];
        double facC = xi/(xi-xim1);
        double facD = xim1/(xi-xim1);
        for (int i = 0; i < length; ++i) {
            double CD = C[i] - D[k - j -1][i];
            C[i] = facC * CD; // Only need to keep one C value
            D[k - j - 1][i] = facD * CD; // Keep all D values for recursion
        }
    }
    for (int i = 0; i < length; ++i) {
        y1[i] = D[0][i];
    }
    for (int j = 1; j <= k; ++j) {
        for (int i = 0; i < length; ++i) {
        y1[i] += D[j][i];
        }
    }
}

static void rescale(struct reb_simulation_integrator_bs* ri_bs, double* const y1, double* const y2, double* const scale, int scale_length) {
    for (int i = 0; i < scale_length; ++i) {
        scale[i] = ri_bs->scalAbsoluteTolerance + ri_bs->scalRelativeTolerance *MAX(fabs(y1[i]), fabs(y2[i]));
        scale[i] = ri_bs->scalAbsoluteTolerance + ri_bs->scalRelativeTolerance * 1.0; //TODO. This sets scale to 1 manually
    }
} 

static void nbody_derivatives(struct reb_ode_state* state, double* const yDot, const double* const y, double const t){
    struct reb_simulation* const r = (struct reb_simulation* const)(state->ref);
    for (int i=0; i<r->N; i++){
         struct reb_particle* const p = &(r->particles[i]);
         p->x  = y[i*6+0];
         p->y  = y[i*6+1];
         p->z  = y[i*6+2];
         p->vx = y[i*6+3];
         p->vy = y[i*6+4];
         p->vz = y[i*6+5];
    }
    reb_update_acceleration(r);

    for (int i=0; i<r->N; i++){
        const struct reb_particle p = r->particles[i];
        yDot[i*6+0] = p.vx;
        yDot[i*6+1] = p.vy;
        yDot[i*6+2] = p.vz;
        yDot[i*6+3] = p.ax;
        yDot[i*6+4] = p.ay;
        yDot[i*6+5] = p.az;
    }
}




void reb_integrator_bs_part1(struct reb_simulation* r){
}

static void allocate_sequence_arrays(struct reb_simulation_integrator_bs* ri_bs){
    ri_bs->sequence        = malloc(sizeof(int)*sequence_length);
    ri_bs->costPerStep     = malloc(sizeof(int)*sequence_length);
    ri_bs->coeff           = malloc(sizeof(double)*sequence_length);
    ri_bs->costPerTimeUnit = malloc(sizeof(double)*sequence_length);
    ri_bs->optimalStep     = malloc(sizeof(double)*sequence_length);

    // step size sequence: 2, 6, 10, 14, ...  // only needed for dense output
     for (int k = 0; k < sequence_length; ++k) {
        ri_bs->sequence[k] = 4 * k + 2;
    }
    
    // step size sequence: 1,2,3,4,5 ...
    //for (int k = 0; k < sequence_length; ++k) {
    //    ri_bs->sequence[k] = 2*( k+1);
    //}

    // initialize the order selection cost array
    // (number of function calls for each column of the extrapolation table)
    ri_bs->costPerStep[0] = ri_bs->sequence[0] + 1;
    for (int k = 1; k < sequence_length; ++k) {
        ri_bs->costPerStep[k] = ri_bs->costPerStep[k - 1] + ri_bs->sequence[k];
    }
    ri_bs->costPerTimeUnit[0]       = 0;

    // initialize the extrapolation tables
    for (int j = 0; j < sequence_length; ++j) {
        double r = 1./((double) ri_bs->sequence[j]);
        ri_bs->coeff[j] = r*r;
    }
}

static void allocate_data_arrays(struct reb_ode_state* state){
    const int length = state->length;
    // create some internal working arrays
    
    state->D   = malloc(sizeof(double*)*(sequence_length));
    for (int k = 0; k < sequence_length; ++k) {
        state->D[k]   = malloc(sizeof(double)*length);
    }

    state->C     = realloc(state->C, sizeof(double)*length);
    state->y1    = realloc(state->y1, sizeof(double)*length);
    state->y0Dot = realloc(state->y0Dot, sizeof(double)*length);
    state->yTmp  = realloc(state->yTmp, sizeof(double)*length);
    state->yDot  = realloc(state->yDot, sizeof(double)*length);

    state->scale = realloc(state->scale, sizeof(double)*length);

}


void reb_integrator_bs_step(struct reb_simulation_integrator_bs* ri_bs){
    if (ri_bs->sequence==NULL){
        allocate_sequence_arrays(ri_bs);
    }

    if (ri_bs->state.length>ri_bs->state.allocatedN){
        allocate_data_arrays(&(ri_bs->state));
        ri_bs->state.allocatedN = ri_bs->state.length;
        ri_bs->firstOrLastStep = 1;
    }

    rescale(ri_bs, ri_bs->state.y, ri_bs->state.y, ri_bs->state.scale, ri_bs->state.length); // initial scaling

    // initial order selection
    if (ri_bs->targetIter == 0){
        const double tol    = ri_bs->scalRelativeTolerance;
        const double log10R = log10(MAX(1.0e-10, tol));
        ri_bs->targetIter = MAX(1, MIN(sequence_length - 2, (int) floor(0.5 - 0.6 * log10R)));
    }

    double  maxError                 = DBL_MAX;

    int y_length = ri_bs->state.length;
    double error;
    int reject = 0;

    // first evaluation, at the beginning of the step
    if (ri_bs->method == 1){ // Note: only for midpoint. leapfrog calculates it itself
        ri_bs->state.derivatives(&(ri_bs->state), ri_bs->state.y0Dot, ri_bs->state.y, ri_bs->state.t);
    }

    const int forward = (ri_bs->hNew >= 0.);
    const double stepSize = ri_bs->hNew;
    printf("step = %.7e    order== %d\n",stepSize, ri_bs->targetIter);


    // iterate over several substep sizes
    int k = -1;
    for (int loop = 1; loop; ) {

        ++k;
        
        // modified midpoint integration with the current substep
        if ( ! tryStep(&ri_bs->state, k, ri_bs->sequence[k], ri_bs->state.t, stepSize, ri_bs->method)) {

            // the stability check failed, we reduce the global step
            printf("S"); //TODO
            ri_bs->hNew   = fabs(stepSize * stabilityReduction);
            reject = 1;
            loop   = 0;

        } else {
            for (int i = 0; i < y_length; ++i) {
                double CD = ri_bs->state.y1[i];
                ri_bs->state.C[i] = CD;
                ri_bs->state.D[k][i] = CD;
                //if (i==6){
                //    printf("k=%d      y = %.8e\n",k,CD);
                //}

            }

            // the substep was computed successfully
            if (k > 0) {

                // extrapolate the state at the end of the step
                // using last iteration data
                extrapolate(&ri_bs->state, ri_bs->coeff, k);
                rescale(ri_bs, ri_bs->state.y, ri_bs->state.y1, ri_bs->state.scale, y_length);

                // estimate the error at the end of the step.
                error = 0;
                for (int j = 0; j < y_length; ++j) {
                    const double e = (ri_bs->state.C[j]) / ri_bs->state.scale[j];
                    error = MAX(error, e * e);
                }
                //error = sqrt(error / y_length);
                if (isnan(error)) {
                    printf("Error. NaN appearing during integration.");
                    exit(0);
                }

                if ((error > 1.0e25)){ // TODO: Think about what to do when error increases: || ((k > 1) && (error > maxError))) {
                    // error is too big, we reduce the global step
                    printf("R (error= %.5e)",error);  // TODO
                    ri_bs->hNew   = fabs(stepSize * stabilityReduction);
                    reject = 1;
                    loop   = 0;
                } else {

                    maxError = MAX(4 * error, 1.0);

                    // compute optimal stepsize for this order
                    const double exp = 1.0 / (2 * k + 1);
                    double fac = stepControl2 / pow(error / stepControl1, exp);
                    const double power = pow(stepControl3, exp);
                    fac = MAX(power / stepControl4, MIN(1. / power, fac));
                    ri_bs->optimalStep[k]     = fabs(stepSize * fac);
                    ri_bs->costPerTimeUnit[k] = ri_bs->costPerStep[k] / ri_bs->optimalStep[k];

                    // check convergence
                    switch (k - ri_bs->targetIter) {

                        case -1 : // one before target
                            if ((ri_bs->targetIter > 1) && ! ri_bs->previousRejected) {

                                // check if we can stop iterations now
                                if (error <= 1.0) {
                                    // convergence have been reached just before targetIter
                                    loop = 0;
                                } else {
                                    // estimate if there is a chance convergence will
                                    // be reached on next iteration, using the
                                    // asymptotic evolution of error
                                    const double ratio = ((double) ri_bs->sequence[ri_bs->targetIter] * ri_bs->sequence[ri_bs->targetIter + 1]) / (ri_bs->sequence[0] * ri_bs->sequence[0]);
                                    if (error > ratio * ratio) {
                                        // we don't expect to converge on next iteration
                                        // we reject the step immediately and reduce order
                                        reject = 1;
                                        loop   = 0;
                                        ri_bs->targetIter = k;
                                        if ((ri_bs->targetIter > 1) &&
                                                (ri_bs->costPerTimeUnit[ri_bs->targetIter - 1] <
                                                 orderControl1 * ri_bs->costPerTimeUnit[ri_bs->targetIter])) {
                                            ri_bs->targetIter -= 1;
                                        }
                                        ri_bs->hNew = ri_bs->optimalStep[ri_bs->targetIter];
                                        printf("O"); // TODO
                                    }
                                }
                            }
                            break;

                        case 0: // exactly on target
                            if (error <= 1.0) {
                                // convergence has been reached exactly at targetIter
                                loop = 0;
                            } else {
                                // estimate if there is a chance convergence will
                                // be reached on next iteration, using the
                                // asymptotic evolution of error
                                const double ratio = ((double) ri_bs->sequence[k + 1]) / ri_bs->sequence[0];
                                if (error > ratio * ratio) {
                                    // we don't expect to converge on next iteration
                                    // we reject the step immediately
                                    printf("o"); // TODO
                                    reject = 1;
                                    loop = 0;
                                    if ((ri_bs->targetIter > 1) &&
                                            (ri_bs->costPerTimeUnit[ri_bs->targetIter - 1] <
                                             orderControl1 * ri_bs->costPerTimeUnit[ri_bs->targetIter])) {
                                        --ri_bs->targetIter;
                                    }
                                    ri_bs->hNew = ri_bs->optimalStep[ri_bs->targetIter];
                                }
                            }
                            break;

                        case 1 : // one past target
                            if (error > 1.0) {
                                printf("e"); // TODO
                                reject = 1;
                                if ((ri_bs->targetIter > 1) &&
                                        (ri_bs->costPerTimeUnit[ri_bs->targetIter - 1] <
                                         orderControl1 * ri_bs->costPerTimeUnit[ri_bs->targetIter])) {
                                    --ri_bs->targetIter;
                                }
                                ri_bs->hNew = ri_bs->optimalStep[ri_bs->targetIter];
                            }
                            loop = 0;
                            break;

                        default :
                            if (ri_bs->firstOrLastStep && (error <= 1.0)) {
                                loop = 0;
                            }
                            break;

                    }
                }
            }
        }
    }


    if (! reject) {
        printf("."); // TODO
        ri_bs->state.t += stepSize;
        // Swap arrays
        double* y_tmp = ri_bs->state.y;
        ri_bs->state.y = ri_bs->state.y1; 
        ri_bs->state.y1 = y_tmp; 

        int optimalIter;
        if (k == 1) {
            optimalIter = 2;
            if (ri_bs->previousRejected) {
                optimalIter = 1;
            }
        } else if (k <= ri_bs->targetIter) { // Converged before or on target
            optimalIter = k;
            if (ri_bs->costPerTimeUnit[k - 1] < orderControl1 * ri_bs->costPerTimeUnit[k]) {
                optimalIter = k - 1;
            } else if (ri_bs->costPerTimeUnit[k] < orderControl2 * ri_bs->costPerTimeUnit[k - 1]) {
                optimalIter = MIN(k + 1, sequence_length - 2);
            }
        } else {                            // converged after target
            optimalIter = k - 1;
            if ((k > 2) && (ri_bs->costPerTimeUnit[k - 2] < orderControl1 * ri_bs->costPerTimeUnit[k - 1])) {
                optimalIter = k - 2;
            }
            if (ri_bs->costPerTimeUnit[k] < orderControl2 * ri_bs->costPerTimeUnit[optimalIter]) {
                optimalIter = MIN(k, sequence_length - 2);
            }
        }

        if (ri_bs->previousRejected) {
            // after a rejected step neither order nor stepsize
            // should increase
            ri_bs->targetIter = MIN(optimalIter, k);
            ri_bs->hNew = MIN(fabs(stepSize), ri_bs->optimalStep[ri_bs->targetIter]);
        } else {
            // stepsize control
            if (optimalIter <= k) {
                ri_bs->hNew = ri_bs->optimalStep[optimalIter];
            } else {
                if ((k < ri_bs->targetIter) &&
                        (ri_bs->costPerTimeUnit[k] < orderControl2 * ri_bs->costPerTimeUnit[k - 1])) {
                    ri_bs->hNew = ri_bs->optimalStep[k] * ri_bs->costPerStep[optimalIter + 1] / ri_bs->costPerStep[k];
                } else {
                    ri_bs->hNew = ri_bs->optimalStep[k] * ri_bs->costPerStep[optimalIter] / ri_bs->costPerStep[k];
                }
            }

            ri_bs->targetIter = optimalIter;

        }
    }

    ri_bs->hNew = fabs(ri_bs->hNew);

    if (ri_bs->hNew < ri_bs->minStep) {
        ri_bs->hNew = ri_bs->minStep;
        printf("Error. Minimal stepsize reached during integration."); // TODO
        exit(0);
    }

    if (ri_bs->hNew > ri_bs->maxStep && ri_bs->maxStep>0.) {
        ri_bs->hNew = ri_bs->maxStep;
        printf("Error. Maximum stepsize reached during integration."); // TODO
        exit(0);
    }

    if (! forward) {
        ri_bs->hNew = -ri_bs->hNew;
    }

    if (reject) {
        ri_bs->previousRejected = 1;
    } else {
        ri_bs->previousRejected = 0;
        ri_bs->firstOrLastStep = 0;
    }
}

void reb_integrator_bs_part2(struct reb_simulation* r){
    double t_initial = r->t;

    struct reb_simulation_integrator_bs* ri_bs = &(r->ri_bs);

    // Setup state for combined N-body + user states
    ri_bs->state.t = r->t;
    int nbody_length = r->N*3*2;
    ri_bs->state.length = nbody_length;
    if (!ri_bs->state.y){
        ri_bs->state.y = malloc(sizeof(double)*ri_bs->state.length);
    }

    {
        double* const y = ri_bs->state.y;
        for (int i=0; i<r->N; i++){
            const struct reb_particle p = r->particles[i];
            y[i*6+0] = p.x;
            y[i*6+1] = p.y;
            y[i*6+2] = p.z;
            y[i*6+3] = p.vx;
            y[i*6+4] = p.vy;
            y[i*6+5] = p.vz;
        }
    }

    
    ri_bs->state.derivatives  = nbody_derivatives;
    ri_bs->state.ref    = r;
    ri_bs->hNew   = r->dt;
    if (r->status==REB_RUNNING_LAST_STEP){
        ri_bs->firstOrLastStep = 1;
    }

    // Generic integrator stuff
    reb_integrator_bs_step(ri_bs);

    // N-body specific:
    {
        double* const y = ri_bs->state.y; // y might have changed
        for (int i=0; i<r->N; i++){
            struct reb_particle* const p = &(r->particles[i]);
            p->x  = y[i*6+0];
            p->y  = y[i*6+1];
            p->z  = y[i*6+2];
            p->vx = y[i*6+3];
            p->vy = y[i*6+4];
            p->vz = y[i*6+5];
        }
    }
    r->t = ri_bs->state.t;
    r->dt = ri_bs->hNew;
    r->dt_last_done = t_initial - r->t;
}

void reb_integrator_bs_synchronize(struct reb_simulation* r){
    // Do nothing.
}


void reb_integrator_bs_reset_struct(struct reb_simulation_integrator_bs* ri_bs){

    // Free data array
    free(ri_bs->state.y1);
    ri_bs->state.y1 = NULL;
    free(ri_bs->state.C);
    ri_bs->state.C = NULL;
    free(ri_bs->state.scale);
    ri_bs->state.scale = NULL;
    
    if (ri_bs->state.D){
        for (int k = 0; k < sequence_length; ++k) {
            ri_bs->state.D[k] = NULL;
        }
        free(ri_bs->state.D);
        ri_bs->state.D = NULL;
    }
    if (ri_bs->state.y0Dot){
        free(ri_bs->state.y0Dot);
        ri_bs->state.y0Dot = NULL;
    }
    if (ri_bs->state.yTmp){
        free(ri_bs->state.yTmp);
        ri_bs->state.yTmp = NULL;
    }
    if (ri_bs->state.yDot){
        free(ri_bs->state.yDot);
        ri_bs->state.yDot = NULL;
    }

    // Free sequence arrays
    free(ri_bs->sequence);
    ri_bs->sequence = NULL;
    
    free(ri_bs->coeff);
    ri_bs->coeff = NULL;
    free(ri_bs->costPerStep);
    ri_bs->costPerStep = NULL;
    free(ri_bs->costPerTimeUnit);
    ri_bs->costPerTimeUnit = NULL;
    free(ri_bs->optimalStep);
    ri_bs->optimalStep = NULL;
    
    
    // Default settings
    ri_bs->scalAbsoluteTolerance= 1e-5;
    ri_bs->scalRelativeTolerance= 1e-5;
    ri_bs->maxStep              = 10; // Note: always positive
    ri_bs->minStep              = 1e-8; // Note: always positive
    ri_bs->firstOrLastStep      = 1;
    ri_bs->previousRejected     = 0;
    ri_bs->method               = 1;  // 1== midpoint
        
}

void reb_integrator_bs_reset(struct reb_simulation* r){
    struct reb_simulation_integrator_bs* ri_bs = &(r->ri_bs);
    reb_integrator_bs_reset_struct(ri_bs);
}
