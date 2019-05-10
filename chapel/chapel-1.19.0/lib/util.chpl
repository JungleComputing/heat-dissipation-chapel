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


use heat_dissipation;

const FLO_PER_POINT_STENCIL = 14;
const FLO_PER_POINT_DIFFERENCE = 2;
const FLO_PER_POINT_REDUCTION = 1;
const FLO_PER_POINT_PER_ITERATION =
  FLO_PER_POINT_STENCIL +
  FLO_PER_POINT_DIFFERENCE +
  FLO_PER_POINT_REDUCTION;


record results {
  var niter : int;
  var tmin : real;
  var tmax : real;
  var maxdiff : real;
  var tavg : real;
  var time : real;
}

proc print_parameters() {
  if help_params {
    usage();
    exit();
  }

  writeln("Parameters:\n",
	  "  --N ", N, " # number of rows\n",
	  "  --M ", M, " # number of columns\n",
	  "  --I ", I, " # maximum number of iterations\n",
	  "  --E ", E, " # convergence threshold\n",
	  "  --C ", C, " # input file for conductivity\n",
	  "  --T ", T, " # input file for initial temperatures\n",
	  "  --L ", L, " # coolest temperature in input/output\n",
	  "  --H ", H, " # highest temperature in input/output\n");
}

proc print_header(){
  var dataParTasks = if dataParTasksPerLocale == 0 
    then here.numPUs(false)
      else dataParTasksPerLocale;
							  
  writeln("Output from heat-dissipation with ", numLocales, " locale",
	  if numLocales == 1 then "" else "s", 
	  " and ",
	  dataParTasks, " data parallel task",
	  if dataParTasks == 1 then "" else "s", ":\n\n",
	  "   Iterations",
	  "        T(min)", 
	  "        T(max)", 
	  "       T(diff)", 
	  "        T(avg)", 
	  "          Time",
	  "        FLOP/s");
}

proc report_results(r : results) 
{
  writeln("% 13i ".format(r.niter), 
	  "% 13.6r ".format(r.tmin),
	  "% 13.6r ".format(r.tmax),
	  "% 13.6r ".format(r.maxdiff),
	  "% 13.6r ".format(r.tavg),
	  "% 13.6r ".format(r.time),
	  "% 13.6er".format((N:real * M:real * 
			     r.niter:real * FLO_PER_POINT_PER_ITERATION) 
			    / r.time));
}

proc usage()
{
  writeln("Usage:  [OPTION]...\n",
	  "\n",
	  "  --N NUM         Set cylinder height to ROWS.\n",
	  "  --M NUM         Set cylinder width to COLUMNS.\n",
	  "  --I NUM         Set the maximum number of iterations to NUM.\n",
	  "  --E NUM         The the convergence threshold to NUM.\n",
	  "  --C FILE        Read conductivity values from FILE.\n",
	  "  --T FILE        Read initial temperature values from FILE.\n",
	  "  --L NUM         Lowest temperature in input/output images.\n",
	  "  --H NUM         Highest temperature in input/output images.\n",
	  "  --help_params   Print this help.\n"
	  );
  exit(0);
}


proc readpgm(fname, height, width, D, array : [D] real, min, max) {
  writeln("Reading PGM image.");
  var f = open(fname, iomode.r).reader();

  var header = f.readln(string);
  var (w,h) = f.readln(int, int);
  var maxv = f.readln(int);

  assert(w == width && h == height);

  for i in 1..height {
    for j in 1..width {
      array[i,j] = min + f.read(int):real * (max - min) / maxv:real;
    }
  }
}

proc main() {

  print_header();

  var r = do_compute();

  report_results(r);
}


