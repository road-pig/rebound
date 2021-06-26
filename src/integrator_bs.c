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
#include "rebound.h"
#include "integrator_bs.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void setStabilityCheck(struct reb_simulation_integrator_bs* ri_bs, int performStabilityCheck, int maxNumIter, int maxNumChecks, double stepsizeReductionFactor) {
    ri_bs->performTest = performStabilityCheck;
    ri_bs->maxIter     = (maxNumIter   <= 0) ? 2 : maxNumIter;
    ri_bs->maxChecks   = (maxNumChecks <= 0) ? 1 : maxNumChecks;

    if ((stepsizeReductionFactor < 0.0001) || (stepsizeReductionFactor > 0.9999)) {
        ri_bs->stabilityReduction = 0.5;
    } else {
        ri_bs->stabilityReduction = stepsizeReductionFactor;
    }
}

void setControlFactors(struct reb_simulation_integrator_bs* ri_bs, double control1, double control2, double control3, double control4) {

    if ((control1 < 0.0001) || (control1 > 0.9999)) {
        ri_bs->stepControl1 = 0.65;
    } else {
        ri_bs->stepControl1 = control1;
    }

    if ((control2 < 0.0001) || (control2 > 0.9999)) {
        ri_bs->stepControl2 = 0.94;
    } else {
        ri_bs->stepControl2 = control2;
    }

    if ((control3 < 0.0001) || (control3 > 0.9999)) {
        ri_bs->stepControl3 = 0.02;
    } else {
        ri_bs->stepControl3 = control3;
    }

    if ((control4 < 1.0001) || (control4 > 999.9)) {
        ri_bs->stepControl4 = 4.0;
    } else {
        ri_bs->stepControl4 = control4;
    }

}

void initializeArrays(struct reb_simulation_integrator_bs* ri_bs) {

    int size = ri_bs->maxOrder / 2;

    if ((ri_bs->sequence == NULL) || (ri_bs->sequence_length != size)) {
        // all arrays should be reallocated with the right size
        ri_bs->sequence        = realloc(ri_bs->sequence,sizeof(int)*size);
        ri_bs->costPerStep     = realloc(ri_bs->costPerStep,sizeof(int)*size);
        ri_bs->coeff           = realloc(ri_bs->coeff,sizeof(double*)*size);
        for (int k = ri_bs->sequence_length; k < size; ++k) {
            ri_bs->coeff[k] = NULL;
        }
        ri_bs->costPerTimeUnit = realloc(ri_bs->costPerTimeUnit,sizeof(double)*size);
        ri_bs->optimalStep     = realloc(ri_bs->optimalStep,sizeof(double)*size);
    }

    // step size sequence: 2, 6, 10, 14, ...
    for (int k = 0; k < size; ++k) {
        ri_bs->sequence[k] = 4 * k + 2;
    }

    // initialize the order selection cost array
    // (number of function calls for each column of the extrapolation table)
    ri_bs->costPerStep[0] = ri_bs->sequence[0] + 1;
    for (int k = 1; k < size; ++k) {
        ri_bs->costPerStep[k] = ri_bs->costPerStep[k - 1] + ri_bs->sequence[k];
    }

    // initialize the extrapolation tables
    for (int k = 1; k < size; ++k) {
        ri_bs->coeff[k] = realloc(ri_bs->coeff[k],sizeof(double)*k);
        for (int l = 0; l < k; ++l) {
            double ratio = ((double) ri_bs->sequence[k]) / ri_bs->sequence[k - l - 1];
            ri_bs->coeff[k][l] = 1.0 / (ratio * ratio - 1.0);
        }
    }

}

void setInterpolationControl(struct reb_simulation_integrator_bs* ri_bs, int useInterpolationErrorForControl, int mudifControlParameter) {

    ri_bs->useInterpolationError = useInterpolationErrorForControl;

    if ((mudifControlParameter <= 0) || (mudifControlParameter >= 7)) {
        ri_bs->mudif = 4;
    } else {
        ri_bs->mudif = mudifControlParameter;
    }

}

void setOrderControl(struct reb_simulation_integrator_bs* ri_bs, int maximalOrder, double control1, double control2) {

    if (maximalOrder > 6 && maximalOrder % 2 == 0) {
        ri_bs->maxOrder = maximalOrder;
    } else {
        ri_bs->maxOrder = 18;
    }

    if ((control1 < 0.0001) || (control1 > 0.9999)) {
        ri_bs->orderControl1 = 0.8;
    } else {
        ri_bs->orderControl1 = control1;
    }

    if ((control2 < 0.0001) || (control2 > 0.9999)) {
        ri_bs->orderControl2 = 0.9;
    } else {
        ri_bs->orderControl2 = control2;
    }

    // reinitialize the arrays
    initializeArrays(ri_bs);

}

void bs_constructor(struct reb_simulation* r){
    struct reb_simulation_integrator_bs* ri_bs = &(r->ri_bs);
    ri_bs->minStep = fabs(ri_bs->minStep);
    ri_bs->maxStep = fabs(ri_bs->maxStep);
    ri_bs->initialStep = -1;
    setStabilityCheck(ri_bs, 1, -1, -1, -1);
    setControlFactors(ri_bs, -1, -1, -1, -1);
    setOrderControl(ri_bs, -1, -1, -1);
    setInterpolationControl(ri_bs, 1, -1);
}

double* computeDerivatives(double t, double* const yEnd){
    // DUMMY
    return 0; // TODO 
}


int tryStep(struct reb_simulation_integrator_bs* ri_bs, const double t0, const double* y0, const int y0_length, const double step, const int k, const double* scale, const int scale_length, const double** f, double* const yMiddle, double* const yEnd) {

    const int    n        = ri_bs->sequence[k];
    const double subStep  = step / n;
    const double subStep2 = 2 * subStep;

    // first substep
    double t = t0 + subStep;
    for (int i = 0; i < y0_length; ++i) {
        yEnd[i] = y0[i] + subStep * f[0][i];
    }
    f[1] = computeDerivatives(t, yEnd);

    // other substeps
    double* const yTmp = malloc(sizeof(double)*y0_length); // IMPROVE: should allocate this only once
    for (int j = 1; j < n; ++j) {

        if (2 * j == n) {
            // save the point at the middle of the step
            for (int i = 0; i < y0_length; ++i) {
                yMiddle[i] = yEnd[i];
            }
        }

        t += subStep;
        for (int i = 0; i < y0_length; ++i) {
            const double middle = yEnd[i];
            yEnd[i]       = y0[i] + subStep2 * f[j][i];
            yTmp[i]       = middle;
        }

        f[j + 1] = computeDerivatives(t, yEnd);

        // stability check
        if (ri_bs->performTest && (j <= ri_bs->maxChecks) && (k < ri_bs->maxIter)) {
            double initialNorm = 0.0;
            for (int l = 0; l < scale_length; ++l) {
                const double ratio = f[0][l] / scale[l];
                initialNorm += ratio * ratio;
            }
            double deltaNorm = 0.0;
            for (int l = 0; l < scale_length; ++l) {
                const double ratio = (f[j + 1][l] - f[0][l]) / scale[l];
                deltaNorm += ratio * ratio;
            }
            if (deltaNorm > 4 * MAX(1.0e-15, initialNorm)) {
                return 0;
            }
        }

    }

    // correction of the last substep (at t0 + step)
    for (int i = 0; i < y0_length; ++i) {
        yEnd[i] = 0.5 * (yTmp[i] + yEnd[i] + subStep * f[n][i]);
    }

    free(yTmp);
    return 1;

}

void extrapolate(struct reb_simulation_integrator_bs* ri_bs, const int offset, const int k, double** const diag, double* const last, const int last_length) {
    // update the diagonal
    for (int j = 1; j < k; ++j) {
        for (int i = 0; i < last_length; ++i) {
            // Aitken-Neville's recursive formula
            diag[k - j - 1][i] = diag[k - j][i] +
                ri_bs->coeff[k + offset][j - 1] * (diag[k - j][i] - diag[k - j - 1][i]);
        }
    }

    // update the last element
    for (int i = 0; i < last_length; ++i) {
        // Aitken-Neville's recursive formula
        last[i] = diag[0][i] + ri_bs->coeff[k + offset][k - 1] * (diag[0][i] - last[i]);
    }
}



void reb_integrator_bs_part1(struct reb_simulation* r){
	r->t+=r->dt/2.;
}
void reb_integrator_bs_part2(struct reb_simulation* r){
	r->t+=r->dt/2.;
	r->dt_last_done = r->dt;
}
	
void reb_integrator_bs_synchronize(struct reb_simulation* r){
	// Do nothing.
}

void reb_integrator_bs_reset(struct reb_simulation* r){
	// Do nothing.
}
