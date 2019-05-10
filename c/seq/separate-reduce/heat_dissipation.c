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


#include "util.h"

#include <sys/time.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>


static void exchange_columns(size_t h, size_t w,
			     double (*restrict g)[h][w]) {
  size_t i;

  /* copy left and right column to opposite border */
  for (i = 0; i < h; ++i) {
    (*g)[i][w-1] = (*g)[i][1];
    (*g)[i][0] = (*g)[i][w-2];
  }
}


/* Does the reduction step and return if the convergence has setteled */
static void fill_result(const struct parameters *p, struct results *r,
			size_t h, size_t w,
			double (*restrict a)[h][w],
			double (*restrict b)[h][w],
			double iter,
			struct timeval *before) {
  /* compute min/max/avg */
  double tmin = INFINITY, tmax = -INFINITY;
  double sum = 0.0;
  double maxdiff = 0.0;

  /* We have said that the final reduction does not need to be included. */
  struct timeval after;
  gettimeofday(&after, NULL);

  for (size_t i = 1; i < h - 1; ++i)
    for (size_t j = 1; j < w - 1; ++j)
      {
	double v = (*a)[i][j];
	double v_old = (*b)[i][j];
	double diff = fabs(v - v_old);
	sum += v;
	if (tmin > v) tmin = v;
	if (tmax < v) tmax = v;
	if (diff > maxdiff) maxdiff = diff;
      }

  r->niter = iter;
  r->maxdiff = maxdiff;
  r->tmin = tmin;
  r->tmax = tmax;
  r->tavg = sum / (p->N * p->M);

  r->time = (double)(after.tv_sec - before->tv_sec) +
    (double)(after.tv_usec - before->tv_usec) / 1e6;
}

void print_src(int h, int w, double (* restrict src)[h][w]) {
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      printf("%f ", (*src)[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}
			       

void do_compute(const struct parameters* p, struct results *r) {
  size_t i, j;

  /* alias input parameters */
  const double (*restrict tinit)[p->N][p->M] = (const double (*)[(size_t) p->N][(size_t) p->M])p->tinit;
  const double (*restrict cinit)[p->N][p->M] = (const double (*)[(size_t) p->N][(size_t) p->M])p->conductivity;

  /* allocate grid data */
  const size_t h = p->N + 2;
  const size_t w = p->M + 2;
  double (*restrict g1)[h][w] = malloc(h * w * sizeof(double));
  double (*restrict g2)[h][w] = malloc(h * w * sizeof(double));

  /* allocate halo for conductivities */
  double (*restrict c)[h][w] = malloc(h * w * sizeof(double));

  struct timeval before;

  static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
  static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

  double max_difference = -DBL_MAX;

  /* set initial temperatures and conductivities */
  for (i = 1; i < h - 1; ++i)
    for (j = 1; j < w - 1; ++j)
      {
	(*g1)[i][j] = (*tinit)[i-1][j-1];
	(*g2)[i][j] = (*tinit)[i-1][j-1];
	(*c)[i][j] = (*cinit)[i-1][j-1];
      }

  /* smear outermost row to border */
  for (j = 1; j < w-1; ++j) {
    (*g1)[0][j] = (*g2)[0][j] = (*g1)[1][j];
    (*g1)[h-1][j] = (*g2)[h-1][j] = (*g1)[h-2][j];
  }
  exchange_columns(h, w, g1);
  exchange_columns(h, w, g2);

  /* compute */
  size_t iter;
  double (*restrict src)[h][w] = g1;
  double (*restrict dst)[h][w] = g2;

  /*
   * If initialization should be included in the timings
   * could be a point of discussion.
   */
  gettimeofday(&before, NULL);

  for (iter = 1; iter <= p->maxiter; ++iter) {
    /* printf("iter %d\n", iter); */

    /* print_src(h, w, src); */
    
    /* initialize halo on source */
    exchange_columns(h, w, src);

    max_difference = -DBL_MAX;

    /* compute */
    for (i = 1; i < h - 1; ++i) {
      for (j = 1; j < w - 1; ++j) {
	double w = (*c)[i][j];
	double restw = 1.0 - w;

	(*dst)[i][j] = w * (*src)[i][j] + 
	  ((*src)[i+1][j  ] + (*src)[i-1][j  ] +
	   (*src)[i  ][j+1] + (*src)[i  ][j-1]) * (restw * c_cdir) +

	  ((*src)[i-1][j-1] + (*src)[i-1][j+1] +
	   (*src)[i+1][j-1] + (*src)[i+1][j+1]) * (restw * c_cdiag);
      }
    }

    for (i = 1; i < h - 1; ++i) {
      for (j = 1; j < w - 1; ++j) {
    	double diff = fabs((*dst)[i][j] - (*src)[i][j]);
    	if (diff > max_difference) {
    	  max_difference = diff;
    	}
      }
    }

    /* printf("diff: %f\n", max_difference); */

    { void *tmp = src; src = dst; dst = tmp; }

    if (max_difference < p->threshold) {
      iter++;
      break;
    }
  }

  iter--;
  { void *tmp = src; src = dst; dst = tmp; }

  /* report at end in all cases */
  fill_result(p, r, h, w, dst, src, iter, &before);

  free(c);
  free(g2);
  free(g1);
}


int main(int argc, char **argv) {
  struct parameters p;
  struct results r;

  read_parameters(&p, argc, argv);

  do_compute(&p, &r);

  report_results(&p, &r);

  return 0;
}
