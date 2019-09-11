//
//
// Copyright 2017 Tom Deakin, University of Bristol
//
// This file is part of mega-stream.
//
// mega-stream is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mega-stream is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mega-stream.  If not, see <http://www.gnu.org/licenses/>.
//
// This aims to investigate the limiting factor for a simple kernel, in particular
// where bandwidth limits not to be reached, and latency becomes a dominating factor.
//

// Collect and print PAPI performance counters for the main kernel
// If not built with PAPI, these functions do nothing to save complicating the Fortran.
// This use of PAPI is inspired by XSBench

#ifdef PAPI
#include <papi.h>

static int event_set = PAPI_NULL;
static int num_events = 0;
static long_long *values = NULL;

#endif

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

// Initialise the PAPI library with the low-level API
void papi_init() {
#ifdef PAPI
  printf("Initilising PAPI...\n");
  int err;
  err = PAPI_library_init(PAPI_VER_CURRENT);
  if (err != PAPI_VER_CURRENT) {
    fprintf(stderr, "PAPI init error\n");
    exit(EXIT_FAILURE);
  }
  if (PAPI_thread_init(omp_get_thread_num) != PAPI_OK) {
    fprintf(stderr, "PAPI thread init error\n");
    exit(EXIT_FAILURE);
  }

#endif
}

// Start recording the chosen events
void papi_start() {
#ifdef PAPI
  printf("Starting PAPI counters...\n");  

  // List of events to collect
  int[] events = {PAPI_TOT_INS,PAPI_LD_INS,PAPI_FP_INS};

  num_events = sizeof(events) / sizeof(int);

  // Allocate data now to store the values of the counters when read at end
  values = malloc(sizeof(long_long) * num_events);

  // Create event set
  if (PAPI_create_eventset(&event_set) != PAPI_OK) {
    fprintf(stderr, "Error: creating PAPI event set\n");
    exit(EXIT_FAILURE);
  }

  // Add events to the set
  for (int i = 0; i < num_events; ++i) {
    if (PAPI_add_event(event_set, events[i]) != PAPI_OK) {
      fprintf(stderr, "Error: adding PAPI event %d\n", events[i]);
      exit(EXIT_FAILURE);
    }
  }


  // Start recording!
  if (PAPI_start(event_set) != PAPI_OK) {
    fprintf(stderr, "Error: starting PAPI\n");
    exit(EXIT_FAILURE);
  }
#endif
}

// Stop recording the events and print out their values
// TODO? Should call this in parallel?
void papi_stop() {
#ifdef PAPI

  // Stop counting!
  if (PAPI_stop(event_set, values) != PAPI_OK) {
    fprintf(stderr, "Error: stopping PAPI\n");
    exit(EXIT_FAILURE);
  }

  // Recover events array from event_set
  int * events = malloc(sizeof(int) * num_events);
  int n = num_events;
  PAPI_list_events(event_set, events, &n);

  // Print the values
  for (int i = 0; i < num_events; ++i) {
    PAPI_event_info_t info;
    PAPI_get_event_info(events[i], &info);
    printf("%-15lld\t%s\t%s\n", values[i], info.symbol, info.long_descr);
  }

  free(events);
  free(values);

#endif
}

