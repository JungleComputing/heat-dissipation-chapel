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


use util;
use Time;

config const N = 2000;
config const M = 2000;
config const I = 500;
config const E = 0.01;
config const L = -100.0;
config const H = 100.0;
config const C = "images/pat2_2000x2000.pgm";
config const T = "images/plasma_2000x2000.pgm";
config const help_params = false;

print_parameters();

use StencilDist;

var MyLocaleView = {0..numLocales-1, 0..0};
var MyLocales = reshape(Locales[0..numLocales-1], MyLocaleView);

const HaloDomain: domain(2) dmapped Stencil(boundingBox = {0..N+1, 0..M+1},
					    targetLocales = MyLocales,
					    fluff=(1,1)) =
  {0..N+1, 0..M+1};
const CylinderDomain: subdomain(HaloDomain) = {1..N, 1..M};

const LeftHalo: subdomain(HaloDomain) = {1..N, 0..0};
const RightHalo: subdomain(HaloDomain) = {1..N, M+1..M+1};
const UpperHalo: subdomain(HaloDomain) = {0..0, 0..M+1};
const LowerHalo: subdomain(HaloDomain) = {N+1..N+1, 0..M+1};

const LeftColumn: subdomain(HaloDomain) = {1..N, 1..1};
const RightColumn: subdomain(HaloDomain) = {1..N, M..M};
const UpperRow: subdomain(HaloDomain) = {1..1, 0..M+1};
const LowerRow: subdomain(HaloDomain) = {N..N, 0..M+1};

class Cylinder {
  var temperature : [HaloDomain] real;
}

const factor_direct_neighbors : real = sqrt(2.0) / (sqrt(2.0) + 1) / 4;
const factor_diagonal_neighbors : real = 1 / (sqrt(2.0) + 1) / 4;

var tinit: [CylinderDomain] real;
readpgm(T, N, M, CylinderDomain, tinit, L, H);
var tcond: [CylinderDomain] real;
readpgm(C, N, M, CylinderDomain, tcond, 0.0, 1.0);


proc exchangeColumns(temperature : [HaloDomain] real) {
  forall (i,j) in zip(LeftHalo, RightColumn) {
    temperature[i] = temperature[j];
  }
  forall (i,j) in zip(RightHalo, LeftColumn) {
    temperature[i] = temperature[j];
  }
}

proc initializeHalos(temperature : [HaloDomain] real) {
  exchangeColumns(temperature);
  temperature[UpperHalo] = temperature[UpperRow];
  temperature[LowerHalo] = temperature[LowerRow];
}

proc do_compute() {
  var src : Cylinder;
  src = new unmanaged Cylinder();
  src.temperature[CylinderDomain] = tinit;
  initializeHalos(src.temperature);
  src.temperature.updateFluff();

  var dst : Cylinder;
  dst = new unmanaged Cylinder();
  dst.temperature[LeftHalo] = src.temperature[LeftHalo];
  dst.temperature[RightHalo] = src.temperature[RightHalo];
  dst.temperature[UpperHalo] = src.temperature[UpperHalo];
  dst.temperature[LowerHalo] = src.temperature[LowerHalo];

  var conductivity : [CylinderDomain] real = tcond;

  var nrIterations = 0;
  var difference = 0.0;
  var max_difference : real = min(real);

  var t : Timer;
  var exchangeTimer : Timer;
  var convolutionTimer : Timer;
  var fluffTimer : Timer;
  var reductionTimer : Timer;
  var swapTimer : Timer;

  t.start();

  for iteration in 1..I {

    ref s = src.temperature;
    ref d = dst.temperature;

    exchangeTimer.start();
    exchangeColumns(s);
    exchangeTimer.stop();

    convolutionTimer.start();
    forall (i, j) in CylinderDomain {
      var weight = conductivity.localAccess[i, j];
      var remaining_weight = 1 - weight;

      d.localAccess[i, j] =
	weight * s.localAccess[i, j] +
	// four direct neighbors
	remaining_weight * factor_direct_neighbors *
	(s.localAccess[i-1, j] +
	 s.localAccess[i, j+1] +
	 s.localAccess[i+1, j] +
	 s.localAccess[i, j-1]) +
	// four diagonal neighbors
	remaining_weight * factor_diagonal_neighbors *
	 (s.localAccess[i-1, j-1] +
	  s.localAccess[i-1, j+1] +
	  s.localAccess[i+1, j+1] +
	  s.localAccess[i+1, j-1]);
    }
    convolutionTimer.stop();
    fluffTimer.start();
    d.updateFluff();
    fluffTimer.stop();

    reductionTimer.start();
    max_difference = max reduce [ij in CylinderDomain] abs(d.localAccess[ij] - s.localAccess[ij]);
    reductionTimer.stop();

    swapTimer.start();
    src <=> dst;
    swapTimer.stop();
    nrIterations = iteration;

    if max_difference < E {
      break;
    }
  }

  dst <=> src;
  
  t.stop();

  var r : results;

  r.tmin = min reduce dst.temperature[CylinderDomain];
  r.tmax = max reduce dst.temperature[CylinderDomain];
  r.maxdiff = max_difference;
  r.niter = nrIterations;
  r.tavg = + reduce dst.temperature[CylinderDomain] / CylinderDomain.size;
  r.time = t.elapsed();
  
  writeln("Exchange: ", exchangeTimer.elapsed(), ", s");
  writeln("Convolution: ", convolutionTimer.elapsed(), ", s");
  writeln("Fluff: ", fluffTimer.elapsed(), ", s");
  writeln("Reduction: ", reductionTimer.elapsed(), ", s");
  writeln("Swap: ", swapTimer.elapsed(), ", s");
  writeln("Total: ", t.elapsed(), ", s");

  return r;
}

util.main();
