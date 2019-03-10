/**
 * Implementation of the Particles functions. Most of your edits should be in here.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "particles.h"
#include "util.h"

using namespace std;

void Particles::
add_particle(const Vec2f &px, const Vec2f &pu)
{
   Particle p = Particle(px, pu);
   P.push_back(p);
   ++np;
}

template<class T>
void Particles::
accumulate(T &accum, float q, int i, int j, float fx, float fy, Array2f &sum)
{
   float weight;

   weight=(1-fx)*(1-fy);
   accum(i,j)+=weight*q;
   sum(i,j)+=weight;

   weight=fx*(1-fy);
   accum(i+1,j)+=weight*q;
   sum(i+1,j)+=weight;

   weight=(1-fx)*fy;
   accum(i,j+1)+=weight*q;
   sum(i,j+1)+=weight;

   weight=fx*fy;
   accum(i+1,j+1)+=weight*q;
   sum(i+1,j+1)+=weight;
}

float Particles::get_mass(float px, float py) {
   int i, j;
   float fx, fy;
   grid.bary_x(px, i, fx);
   std::cout << "fx = " << fx << std::endl;
   grid.bary_y_centre(py, j, fy);
   float mu = sum_u.bilerp(i, j, fx, fy);
   std::cout << "mu = " << mu << std::endl;

   grid.bary_x_centre(px, i, fx);
   grid.bary_y(py, j, fy);
   std::cout << "fy = " << fy << std::endl;
   float mv = sum_v.bilerp(i, j, fx, fy);
   std::cout << "mv = " << mv << std::endl;

   Vec2f mass = Vec2f(mu, mv);
   std::cout << mag(mass) << endl;
   return mag(mass);
}

void Particles::update_vol_dens(void) {
   Vec2f midx, gu;
   float xmin=1.001*grid.h, xmax=grid.lx-1.001*grid.h;
   float ymin=1.001*grid.h, ymax=grid.ly-1.001*grid.h;
   for(int p=0; p<np; ++p){
      P[p].d = get_mass(P[p].x[0], P[p].x[1]) / (grid.h * grid.h);
      P[p].V = (P[p].d == 0)? 0: 1.0/P[p].d;
   }
}


void Particles::
transfer_to_grid(void)
{
   int p, i, ui, j, vj;
   float fx, ufx, fy, vfy;

   grid.u.zero();
   sum_u.zero();

   for(int p=0; p<np; ++p) {
      grid.bary_x(P[p].x[0], ui, ufx);
      grid.bary_y_centre(P[p].x[1], j, fy);
      accumulate(grid.u, P[p].u[0], ui, j, ufx, fy, sum_u);
   }

   for(j=0; j<grid.u.ny; ++j) for(i=0; i<grid.u.nx; ++i){
      if(sum_u(i,j)!=0) grid.u(i,j)/=sum_u(i,j);
   }

   grid.v.zero();
   sum_v.zero();

   for(int p=0; p<np; ++p) {
      grid.bary_y(P[p].x[1], vj, vfy);
      grid.bary_x_centre(P[p].x[0], i, fx);
      accumulate(grid.v, P[p].u[1], i, vj, fx, vfy, sum_v);
   }

   for(j=0; j<grid.v.ny; ++j) for(i=0; i<grid.v.nx; ++i){
      if(sum_v(i,j)!=0) grid.v(i,j)/=sum_v(i,j);
   }

   // identify where particles are in grid
   grid.marker.zero();
   for(int p=0; p<np; ++p) {
      grid.bary_x(P[p].x[0], i, fx);
      grid.bary_y(P[p].x[1], j, fy);
      grid.marker(i,j)=FLUIDCELL;
      P[p].i = i;
      P[p].j = j;
   }
}

/* this function computes c from the gradient of w and the velocity field from the grid. */
Vec2f Particles::
computeC(Array2f &ufield, int i, int j, float fx, float fy)
{
   /* TODO: fill this in */
   return Vec2f(0.f,0.f);
}

void Particles::
update_from_grid(void)
{
   int p;
   int i, ui, j, vj;
   float fx, ufx, fy, vfy;
   for(int p=0; p<np; ++p) {
      grid.bary_x(P[p].x[0], ui, ufx);
      grid.bary_x_centre(P[p].x[0], i, fx);
      grid.bary_y(P[p].x[1], vj, vfy);
      grid.bary_y_centre(P[p].x[1], j, fy);
      P[p].u += Vec2f(grid.du.bilerp(ui, j, ufx, fy), grid.dv.bilerp(i, vj, fx, vfy)); // FLIP
   }
}

void Particles::
move_particles_in_grid(float dt)
{
   Vec2f midx, gu;
   float xmin=1.001*grid.h, xmax=grid.lx-1.001*grid.h;
   float ymin=1.001*grid.h, ymax=grid.ly-1.001*grid.h;

   for(int p=0; p<np; ++p) {
      // first stage of Runge-Kutta 2 (do a half Euler step)
      grid.bilerp_uv(P[p].x[0], P[p].x[1], gu[0], gu[1]);
      midx = P[p].x + 0.5*dt*gu;
      clamp(midx[0], xmin, xmax);
      clamp(midx[1], ymin, ymax);
      // second stage of Runge-Kutta 2
      grid.bilerp_uv(midx[0], midx[1], gu[0], gu[1]);
      P[p].x += dt*gu;
      clamp(P[p].x[0], xmin, xmax);
      clamp(P[p].x[1], ymin, ymax);
   }
}

void Particles::
write_to_file(const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
   char *filename;
   vasprintf(&filename, filename_format, ap);
   FILE *fp=fopen(filename, "wt");
   free(filename);
   va_end(ap);

   fprintf(fp, "%d\n", np);

   for(int p=0; p<np; ++p) {
      fprintf(fp, "%.5g %.5g\n", P[p].x[0], P[p].x[1]);
   }

   fclose(fp);
}

