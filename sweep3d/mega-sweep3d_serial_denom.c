//
// Copyright 2018 Tom Deakin, University of Bristol
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


//
// Serial version of mega-sweep3d designed for use in an architectural simulator
// Data arrays are allocated statically and will likely *exceed* the default stack limit.
// Increase by first running `ulimit -s unlimited`
//


//
// Problem dimensions
//
#define NANG    48
#define NG      1
#define NX      80
#define NY      8
#define NZ      8
#define CHUNK   16
#define NTIMES  5
#define NSWEEPS 8


void sweep(const int nang, const int nx, const int ny, const int nz, const int ng, const int nsweeps, const int chunk,
           double * restrict aflux0, double * restrict aflux1, double * restrict sflux,
           double * restrict psii, double * restrict psij, double * restrict psik,
           double * restrict mu, double * restrict eta, double * restrict xi, double * restrict w,
           double * restrict denom,
           const double v, const double dx, const double dy, const double dz,
           double * restrict y_buf, double * restrict z_buf);


double _start()
{
  // Allocate arrays
  double aflux0[NANG*NX*NY*NZ*NSWEEPS*NG];
  double aflux1[NANG*NX*NY*NZ*NSWEEPS*NG];
  double sflux[NX*NY*NZ*NG];
  double mu[NANG];
  double eta[NANG];
  double xi[NANG];
  double psii[NANG*NY*NZ*NG];
  double psij[NANG*CHUNK*NZ*NG];
  double psik[NANG*CHUNK*NY*NG];
  double w[NANG];
  double denom[NANG*NX*NY*NZ*NG];
  double y_buf[NANG*CHUNK*NZ*NG];
  double z_buf[NANG*CHUNK*NY*NG];

  double v = 0.1;
  double dx = 1.0/(double)NX;
  double dy = 1.0/(double)NY;
  double dz = 1.0/(double)NZ;

  // Initialise data
  for (long i = 0; i < NANG*NX*NY*NZ*NSWEEPS*NG; ++i) {
    aflux0[i] = 1.0;
    aflux1[i] = 0.0;
  }

  for (int i = 0; i < NX*NY*NZ*NG; ++i) {
    sflux[i] = 0.0;
  }

  for (int i = 0; i < NANG; ++i) {
    mu[i] = 0.33;
    eta[i] = 0.66;
    xi[i] = 0.674907401055878; // sqrt(1 - mu*mu - eta*eta)
    w[i] = 0.4;
  }

  for (int g = 0; g < NG; ++g) {
    for (int i = 0; i < NX*NY*NZ; ++i) {
      for (int a = 0; a < NANG; ++a) {
        denom[a+NANG*i+NANG*NX*NY*NZ*g] = 1.0 / (
          0.07 + 2.0*mu[a]/dx + 2.0*eta[a]/dy + 2.0*xi[a]/dz + v);
      }
    }
  }

  for (int i = 0; i < NANG*NY*NZ*NG; ++i) {
    psii[i] = 0.0;
  }

  for (int i = 0; i < NANG*CHUNK*NZ*NG; ++i) {
    psij[i] = 0.0;
    y_buf[i] = 0.0;
  }

  for (int i = 0; i < NANG*CHUNK*NY*NG; ++i) {
    psik[i] = 0.0;
    z_buf[i] = 0.0;
  }

  double * aflux0_ptr = aflux0;
  double * aflux1_ptr = aflux1;

  // Call the sweep routine
  for (int t = 0; t < NTIMES; ++t) {

    sweep(NANG, NX, NY, NZ, NG, NSWEEPS, CHUNK,
          aflux0_ptr, aflux1_ptr, sflux, psii, psij, psik,
          mu, eta, xi, w, denom, v, dx, dy, dz,
          y_buf, z_buf);

    // Swap pointers
    double * tmp = aflux0_ptr;
    aflux0_ptr = aflux1_ptr;
    aflux1_ptr = tmp;

  }

  // Calculate total population
  double pop[NG];
  double population = 0.0;
  for (int g = 0; g < NG; ++g) {
    pop[g] = 0.0;
    for (int i = 0; i < NX*NY*NZ; ++i) {
      pop[g] += sflux[i+NX*NY*NZ*g];
    }
    pop[g] *= dx * dy * dz;
    pop[g] /= (double)(NG-g);
    population += pop[g];
  }

  //printf("Population: %E\n", population);

  return population;
}


void sweep(const int nang, const int nx, const int ny, const int nz, const int ng, const int nsweeps, const int chunk,
           double * restrict aflux0, double * restrict aflux1, double * restrict sflux,
           double * restrict psii, double * restrict psij, double * restrict psik,
           double * restrict mu, double * restrict eta, double * restrict xi, double * restrict w,
           double * restrict denom,
           const double v, const double dx, const double dy, const double dz,
           double * restrict y_buf, double * restrict z_buf) {


  // Calculate number of chunks in x-direction
  int nchunks = nx / chunk;

  // Calculate inverse of spatial width to manually hoist divide operations
  double idx = 1.0 / dx;
  double idy = 1.0 / dy;
  double idz = 1.0 / dz;

  int istep, jstep, kstep;
  int xmin, xmax, ymin, ymax, zmin, zmax, cmin, cmax;

  // Just do ONE sweep direction to reduce the runtime by 8
  //for (int sweep = 0; sweep < nsweeps; ++sweep) {
  int sweep = 0; {

    // Set sweep directions
    // Minimum is inclusive, maximum is exclusive
    // Loop bounds will for (for i = min; i != max; i+=step)
    switch (sweep) {
      case 0:
        istep = -1;
        xmin = chunk-1;
        xmax = -1;
        cmin = nchunks-1;
        cmax = -1;
        jstep = -1;
        ymin = ny-1;
        ymax = -1;
        kstep = -1;
        zmin = nz-1;
        zmax = -1;
        break;

      case 1:
        istep = 1;
        xmin = 0;
        xmax = chunk;
        cmin = 0;
        cmax = nchunks;
        jstep = -1;
        ymin = ny-1;
        ymax = -1;
        kstep = -1;
        zmin = nz-1;
        zmax = -1;
        break;

      case 2:
        istep = -1;
        xmin = chunk-1;
        xmax = -1;
        cmin = nchunks-1;
        cmax = -1;
        jstep = 1;
        ymin = 0;
        ymax = ny;
        kstep = -1;
        zmin = nz-1;
        zmax = -1;
        break;

      case 3:
        istep = 1;
        xmin = 0;
        xmax = chunk;
        cmin = 0;
        cmax = nchunks;
        jstep = 1;
        ymin = 0;
        ymax = ny;
        kstep = -1;
        zmin = nz-1;
        zmax = -1;
        break;


      case 4:
        istep = -1;
        xmin = chunk-1;
        xmax = -1;
        cmin = nchunks-1;
        cmax = -1;
        jstep = -1;
        ymin = ny-1;
        ymax = -1;
        kstep = 1;
        zmin = 0;
        zmax = nz;
        break;

      case 5:
        istep = 1;
        xmin = 0;
        xmax = chunk;
        cmin = 0;
        cmax = nchunks;
        jstep = -1;
        ymin = ny-1;
        ymax = -1;
        kstep = 1;
        zmin = 0;
        zmax = nz;
        break;

      case 6:
        istep = -1;
        xmin = chunk-1;
        xmax = -1;
        cmin = nchunks-1;
        cmax = -1;
        jstep = 1;
        ymin = 0;
        ymax = ny;
        kstep = 1;
        zmin = 0;
        zmax = nz;
        break;

      case 7:
        istep = 1;
        xmin = 0;
        xmax = chunk;
        cmin = 0;
        cmax = nchunks;
        jstep = 1;
        ymin = 0;
        ymax = ny;
        kstep = 1;
        zmin = 0;
        zmax = nz;
        break;

    } // end switch

    // Zero boundary data every sweep
    for (int i = 0; i < nang*ny*nz*ng; ++i)
      psii[i] = 0.0;
    for (int i = 0; i < nang*chunk*nz*ng; ++i)
      psij[i] = 0.0;
    for (int i = 0; i < nang*chunk*ny*ng; ++i)
      psik[i] = 0.0;

    // Loop over chunks
    for (int c = cmin; c != cmax; c += istep) {

      // Fake MPI Recv
      for (int i = 0; i < nang*chunk*nz*ng; ++i)
        psij[i] = 0.0;
      for (int i = 0; i < nang*chunk*ny*ng; ++i)
        psik[i] = 0.0;

      for (int g = 0; g < ng; ++g) { // Loop over groups
        for (int k = zmin; k != zmax; k += kstep) { // Loop over z-dimension
          for (int j = ymin; j != ymax; j += jstep) { // Loop over y-dimension
            for (int ci = xmin; ci != xmax; ci += istep) { // Loop over cells in chunk (x-direction)

              // Calculate i index with respect to nx
              int i = c * chunk + ci;

              for (int a = 0; a < nang; ++a) { // Loop over angles

                // Calculate angular flux
                double psi = (
                  mu[a]  * psii[a + nang*j  + nang*ny*k    + nang*ny*nz*g] +
                  eta[a] * psij[a + nang*ci + nang*chunk*k + nang*chunk*nz*g] +
                  xi[a]  * psik[a + nang*ci + nang*chunk*j + nang*chunk*ny*g] +
                  v * aflux0[a + nang*i + nang*nx*j + nang*nx*ny*k + nang*nx*ny*nz*sweep + nang*nx*ny*nz*nsweeps*g])
                  * denom[a + nang*i + nang*nx*j + nang*nx*ny*k + nang*nx*ny*nz*g];

                // Outgoing diamond difference
                psii[a + nang*j  + nang*ny*k    + nang*ny*nz*g]    = 2.0*psi - psii[a + nang*j  + nang*ny*k    + nang*ny*nz*g];
                psij[a + nang*ci + nang*chunk*k + nang*chunk*nz*g] = 2.0*psi - psij[a + nang*ci + nang*chunk*k + nang*chunk*nz*g];
                psik[a + nang*ci + nang*chunk*j + nang*chunk*ny*g] = 2.0*psi - psik[a + nang*ci + nang*chunk*j + nang*chunk*ny*g];
                aflux1[a + nang*i + nang*nx*j + nang*nx*ny*k + nang*nx*ny*nz*sweep + nang*nx*ny*nz*nsweeps*g] = 2.0*psi - aflux0[a + nang*i + nang*nx*j + nang*nx*ny*k + nang*nx*ny*nz*sweep + nang*nx*ny*nz*nsweeps*g];

                // Reduction
                sflux[i + nx*j + nx*ny*k + nx*ny*nz*g] += psi*w[a];

              } // end angle loop
            } // end x-loop
          } // end y-loop
        } // end z-loop
      } // end group loop

      // Fake MPI Send
      for (int i = 0; i < nang*chunk*nz*ng; ++i)
        y_buf[i] = psij[i];
      for (int i = 0; i < nang*chunk*ny*ng; ++i)
        z_buf[i] = psik[i];

    } // end chunk loop

  } // end sweep loop

}

