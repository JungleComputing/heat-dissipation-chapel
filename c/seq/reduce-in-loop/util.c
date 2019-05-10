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


#include <getopt.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

static double *tinit = 0;
static double *conductivity = 0;

static void cleanup(void) {
  free(tinit);
  free(conductivity);
}

static void die(const char *msg) {
  if (errno != 0) 
    perror(msg);
  else
    fprintf(stderr, "error: %s\n", msg);
  exit(1);
}   


static void readpgm_float(const char *fname,
                          size_t height, size_t width, double *data,
                          double dmin, double dmax) {
  char format[3];
  FILE *f;
  unsigned imgw, imgh, maxv, v;
  size_t i;

  format[2] = '\0';

  printf("Reading PGM data from %s...\n", fname);

  if (!(f = fopen(fname, "r"))) die("fopen");

  fscanf(f, "%2s", format);
  if (format[0] != 'P' || format[1] != '2') die("only ASCII PGM input is supported");
    
  if (fscanf(f, "%u", &imgw) != 1) die("invalid input");
  if (fscanf(f, "%u", &imgh) != 1) die("invalid input");
  if (fscanf(f, "%u", &maxv) != 1) die("invalid input");

  if (imgw != width || imgh != height) {
    fprintf(stderr, "input data size (%ux%u) does not match cylinder size (%zux%zu)\n",
	    imgw, imgh, width, height);
    die("invalid input");
  }

  for (i = 0; i < width * height; ++i)
    {
      if (fscanf(f, "%u", &v) != 1) die("invalid data");
      data[i] = dmin + (double)v * (dmax - dmin) / maxv;
    }

  fclose(f);
}


static void usage(const char *pname) {
  printf("Usage: %s [OPTION]...\n"
	 "\n"
	 "  -N NUM     Set cylinder height to ROWS.\n"
	 "  -M NUM     Set cylinder width to COLUMNS.\n"
	 "  -I NUM     Set the maximum number of iterations to NUM.\n"
	 "  -E NUM     The the convergence threshold to NUM.\n"
	 "  -C FILE    Read conductivity values from FILE.\n"
	 "  -T FILE    Read initial temperature values from FILE.\n"
	 "  -L NUM     Coldest temperature in input/output images.\n"
	 "  -H NUM     Warmest temperature in input/output images.\n"
	 "  -h         Print this help.\n",
	 pname);
  exit(0);
}


void read_parameters(struct parameters* p, int argc, char **argv) {
  const char *conductivity_fname = 0;
  const char *tinit_fname = 0;
  int ch;

  /* set defaults */
  p->N = 150;
  p->M = 100;
  p->maxiter = 42;
  p->threshold = 0.1;
  p->io_tmin = -100.0;
  p->io_tmax = 100.0;
  conductivity_fname = "pattern_100x150.pgm";
  tinit_fname = "pattern_100x150.pgm";

  /* struct option longopt; */
  /* int longindex; */

  //  while ((ch = getopt_long_only(argc, argv, "C:E:h:H:I:L:M:N:T:", &longopt, &longindex)) != -1)
  while ((ch = getopt(argc, argv, "C:E:h:H:I:L:M:N:T:")) != -1)
    {
      switch(ch) {
      case 'C': conductivity_fname = optarg; break;
      case 'T': tinit_fname = optarg; break;
      case 'I': p->maxiter = strtol(optarg, 0, 10); break;
      case 'M': p->M = strtol(optarg, 0, 10); break;
      case 'N': p->N = strtol(optarg, 0, 10); break;
      case 'E': p->threshold = strtod(optarg, 0); break;
      case 'L': p->io_tmin = strtod(optarg, 0); break;
      case 'H': p->io_tmax = strtod(optarg, 0); break;
      case 'h': default: usage(argv[0]);
      }
    }

  printf("Parameters:\n"
	 "  -N %zu # number of rows\n"
	 "  -M %zu # number of columns\n"
	 "  -I %zu # maximum number of iterations\n"
	 "  -E %e # convergence threshold\n"
	 "  -C %s # input file for conductivity\n"
	 "  -T %s # input file for initial temperatures\n"
	 "  -L %e # coolest temperature in input/output\n"
	 "  -H %e # highest temperature in input/output\n"
	 "  -h help\n",
	 p->N, p->M, p->maxiter, p->threshold,
	 conductivity_fname ? conductivity_fname : "(none)",
	 tinit_fname ? tinit_fname : "(none)",
	 p->io_tmin, p->io_tmax);

  if (!p->N || !p->M) die("empty grid");

  atexit(cleanup);

  if (!(tinit = calloc(p->N * p->M, sizeof(double)))) die("calloc");
  if (tinit_fname) 
    readpgm_float(tinit_fname, p->N, p->M, tinit, p->io_tmin, p->io_tmax);
  p->tinit = tinit;

  if (!(conductivity = calloc(p->N * p->M, sizeof(double)))) die("calloc");
  if (conductivity_fname) 
    readpgm_float(conductivity_fname, p->N, p->M, conductivity, 0.0, 1.0);
  p->conductivity = conductivity;
}


#define FPOPS_PER_POINT_PER_ITERATION (14 + 2 + 1)


 
void report_results(const struct parameters *p, const struct results *r) {
  printf("Output:\n\n"
	 "%13s %13s %13s %13s %13s %13s %13s\n",
               "Iterations",
               "T(min)", "T(max)", "T(diff)", "T(avg)", "Time", "FLOP/s");

  printf("%-13zu % 13.6f % 13.6f % 13.6f % 13.6f % 13.6f % 13.6f\n",
	 r->niter,
	 r->tmin,
	 r->tmax,
	 r->maxdiff,
	 r->tavg,
	 r->time,
	 (double)p->N * (double)p->M * 
	 (double)(r->niter * FPOPS_PER_POINT_PER_ITERATION) / r->time);
}

