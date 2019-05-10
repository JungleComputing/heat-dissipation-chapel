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
use BlockDist;
use VisualDebug;
use Barriers;

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

var MyLocaleView = {0..numLocales-1, 0..0};
var MyLocales = reshape(Locales[0..numLocales-1], MyLocaleView);

const HaloDomain: domain(2) dmapped Block(boundingBox = {0..N+1, 0..M+1},
					  targetLocales = MyLocales) = {0..N+1, 0..M+1};
const CylinderDomain: subdomain(HaloDomain)  = {1..N, 1..M};

const LeftHalo: subdomain(HaloDomain) = {1..N, 0..0};
const RightHalo: subdomain(HaloDomain) = {1..N, M+1..M+1};
const UpperHalo: subdomain(HaloDomain) = {0..0, 0..M+1};
const LowerHalo: subdomain(HaloDomain) = {N+1..N+1, 0..M+1};

const LeftColumn: subdomain(HaloDomain) = {1..N, 1..1};
const RightColumn: subdomain(HaloDomain) = {1..N, M..M};
const UpperRow: subdomain(HaloDomain) = {1..1, 0..M+1};
const LowerRow: subdomain(HaloDomain) = {N..N, 0..M+1};

var temperatures: [HaloDomain] real;
readpgm(T, N, M, CylinderDomain, temperatures[CylinderDomain], L, H);
var conductivity: [CylinderDomain] real;
readpgm(C, N, M, CylinderDomain, conductivity, 0.0, 1.0);


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




const Row : domain(1) = {0..M+1};

class Communicator {
  var firstRow: [Row] real;
  var lastRow: [Row] real;
}

var communicators: [LocaleSpace] Communicator;

var max_diffs: [LocaleSpace] real;
var totalTimeCommunicator: [LocaleSpace] real;

const allLocalesBarrier = new Barrier(Locales.size);

proc do_compute() {
  initializeHalos(temperatures);
  var r : results;

  var globalIterations = 0;
  var globalMaxDiff: real = min(real);
  
  var globalTotalTime: real;
  coforall l in Locales with (ref globalIterations, ref globalMaxDiff, ref globalTotalTime) {
    on l {
      const myCommunicator = new unmanaged Communicator();
      communicators[here.id] = myCommunicator;

      allLocalesBarrier.barrier();
      const isLast = here.id == LocaleSpace.last;
      const isFirst = here.id == LocaleSpace.first;

      const beforeCommunicator = if !isFirst then communicators[here.id - 1] else new unmanaged Communicator();
      const nextCommunicator = if !isLast then communicators[here.id + 1] else new unmanaged Communicator();

      const LocalCylinderDomain = CylinderDomain.localSubdomain();
      const LocalHaloDomain = LocalCylinderDomain.expand(1, 1);

      const LocalUpperRow = LocalHaloDomain.dim(1).first;
      const LocalLowerRow = LocalHaloDomain.dim(1).last;

      class Cylinder {
	var temperature : [LocalHaloDomain] real;
      }

      var src: Cylinder = new unmanaged Cylinder();
      var dst: Cylinder = new unmanaged Cylinder();
      forall (i, j) in LocalHaloDomain {
        src.temperature[i, j] = temperatures[i, j];
      }

      dst.temperature = src.temperature;

      const factor_direct_neighbors = (sqrt(2) / (sqrt(2) + 1)) / 4;
      const factor_diagonal_neighbors = (1 / (sqrt(2) + 1)) / 4;

      const localConductivity: [LocalCylinderDomain] real = conductivity[LocalCylinderDomain];

      var totalTimer: Timer;
      var swapTimer: Timer;
      var localTimer: Timer;
      var sendTimer: Timer;
      var barrierTimer: Timer;
      var receiveTimer: Timer;
      var gatherTimer: Timer;
      var rowExchangeTimer: Timer;
      var reduceTimer: Timer;

      var local_max_difference = min(real);
      var local_iteration = 0;
      totalTimer.start();
      do {
	ref s = src.temperature;
	ref d = dst.temperature;

	localTimer.start();
	local {
	  local_max_difference = min(real);
	  
	  forall (i, j) in LocalCylinderDomain with (max reduce local_max_difference) {
	    const weight = localConductivity[i, j];
	    const remaining_weight = 1 - weight;
	    const oldTemp = s[i, j];
	    
	    const newTemp = weight * oldTemp +
	      remaining_weight * factor_direct_neighbors
	      * (s[i-1, j] +
		 s[i, j+1] +
		 s[i+1, j] +
		 s[i, j-1])
	      + remaining_weight * factor_diagonal_neighbors
	      * (s[i-1, j-1] +
		 s[i-1, j+1] +
		 s[i+1, j+1] +
		 s[i+1, j-1]);
	    d[i, j] = newTemp;
	    local_max_difference = max(local_max_difference, abs(newTemp - oldTemp));
	  }
	  forall i in {LocalUpperRow..LocalLowerRow} {
	    d[i, 0] = d[i, M];
	    d[i, M+1] = d[i, 1];
	  }
        }
        localTimer.stop();

        //sendTimer.start();
        rowExchangeTimer.start();

        nextCommunicator.firstRow = d[LocalLowerRow - 1, ..];
        beforeCommunicator.lastRow = d[LocalUpperRow + 1, ..];

        rowExchangeTimer.stop();

	swapTimer.start();
	src <=> dst;
	swapTimer.stop();

        reduceTimer.start();
        max_diffs[here.id] = local_max_difference;
        reduceTimer.stop();

        //sendTimer.stop();

	barrierTimer.start();
        allLocalesBarrier.barrier();
	barrierTimer.stop();

        receiveTimer.start();
        if (!isLast) {
          forall col in Row {
            d[LocalLowerRow, col] = myCommunicator.lastRow[col];
          }
        }
        if (!isFirst) {
          forall col in Row {
            d[LocalUpperRow, col] = myCommunicator.firstRow[col];
          }
        }
        receiveTimer.stop();

        reduceTimer.start();
        local_max_difference = max reduce max_diffs;
        reduceTimer.stop();

        local_iteration += 1;
      } while (local_max_difference > E && local_iteration < I);
      totalTimer.stop();

      totalTimeCommunicator[here.id] = totalTimer.elapsed();

      barrierTimer.start();
      allLocalesBarrier.barrier();
      barrierTimer.stop();

      gatherTimer.start();
      forall (i, j) in LocalCylinderDomain {
        temperatures[i, j] = src.temperature[i, j];
      }
      if (here.id==0) {
        globalTotalTime = max reduce totalTimeCommunicator;
        globalIterations = local_iteration;
        globalMaxDiff = local_max_difference;
      }
      gatherTimer.stop();
      writeln("Locale ", here.id, " Local Timer ", localTimer.elapsed(), " sendTimer ", sendTimer.elapsed(), " barrierTimer ", barrierTimer.elapsed(), 
	      " receiveTimer ", receiveTimer.elapsed(), " gatherTimer ", gatherTimer.elapsed(),
	      " reduceTimer ", reduceTimer.elapsed(), " rowExchangeTimer ", rowExchangeTimer.elapsed());
    }
  }

  r.tmin = min reduce temperatures[CylinderDomain];
  r.tmax = max reduce temperatures[CylinderDomain];
  r.maxdiff = globalMaxDiff;
  r.niter = globalIterations;
  r.tavg = (+ reduce temperatures[CylinderDomain]) / (N * M);
  r.time = globalTotalTime;

  return r;
}

util.main();
