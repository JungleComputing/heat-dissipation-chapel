/*
 * Copyright 2019 Vrije Universiteit Amsterdam, The Netherlands
 *                University of Amsterdam, The Netherlands
 *                Per Fuchs, Pieter Hijma, Clemens Grelck
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>

struct results {
    size_t niter; /* effective number of iterations */
    double tmin; /* minimum temperature in last state */
    double tmax; /* maximum temperature in last state*/
    double maxdiff; /* maximum difference during last update */
    double tavg; /* average temperature */
    double time; /* compute time in seconds */
};



/* input parameters for the program */
struct parameters {
    /* matrix size: N rows, M columns */
    size_t N, M;

    /* maximum number of iterations */
    size_t maxiter;

    /* convergence threshold */
    double threshold; 

    /* initial temperature, in row-major order */
    const double *tinit;

    /* conductivity values for the cylinder, in row-major order */
    const double *conductivity;

    /* temperature range in input and output images */
    double io_tmin;
    double io_tmax;
};
 
void report_results(const struct parameters *p, const struct results *r);
void read_parameters(struct parameters *p, int argc, char **argv);

#endif
