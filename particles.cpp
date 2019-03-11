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
#include <Eigen/Dense>
#include <iostream>
#include "particles.h"
#include "particle.h"
#include "util.h"

using namespace std;

void Particles::
add_particle(const Vec2f &px, const Vec2f &pu)
{
   Particle p = Particle(px, pu);
   P.push_back(p);
   ++np;
}

// STEP 1
void Particles::transfer_mass_to_grid(void) {
   grid.marker.zero();
   float fx, fy;
   float weight;
   grid.mass.zero();
   grid.nodes_grad_weights_x.zero();
   grid.nodes_grad_weights_y.zero();

   for (int n = 0; n < np; n++) {
      P[n].i = (int) (P[n].x(0) / grid.h);
      P[n].j = (int) (P[n].x(1) / grid.h);

      int low_i = max((int)((P[n].x(0) - 2) / grid.h), 0);
      int low_j = max((int)((P[n].x(1) - 2) / grid.h), 0);
      int high_i = min((int)((P[n].x(0) + 2) / grid.h), grid.mass.nx - 1);
      int high_j = min((int)((P[n].x(1) + 2) / grid.h), grid.mass.ny - 1);

      // iterate for all grid cells within distance 2
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {

            float wx = grid.bspline_weight((P[n].x[0] - i * grid.h)/ grid.h);
            float wy = grid.bspline_weight((P[n].x[1] - j * grid.h) / grid.h);
            grid.mass(i, j) += (wx * wy);
            P[n].weights.push_back(wx * wy);

            Eigen::Vector2d grad_weight = Eigen::Vector2d(wx * grid.bspline_gradweight((P[n].x[0] - i * grid.h) / grid.h),
                                       wy * grid.bspline_gradweight((P[n].x[1] - j * grid.h) / grid.h));
            P[n].grad_weights.push_back(grad_weight);
            grid.grad_weights_x(i, j) += grad_weight(0);
            grid.grad_weights_y(i, j) += grad_weight(1);
            
         }
      }

      //adding to marker 
      grid.marker(i, j) = FLUIDCELL;
   }
}

//STEP 1
void Particles::transfer_v_to_grid(void) {
   grid.v_x.zero();
   grid.v_y.zero();

   for (int n = 0; n < np; n++) {
      for (int i = P[n].i ; i <= min(P[n].i + 1, grid.v_x.nx - 1); i++) {
         for (int j = P[n].j; j <= min(P[n].j + 1, grid.v_x.ny - 1); j++) {
            double dx = abs(P[n].x(0) - (i * grid.h)) / grid.h;
            double dy = abs(P[n].y(0) - (j * grid.h)) / grid.h;

            grid.v_x(i, j) += dx * dy * P[n].v(0);
            grid.v_y(i, j) += dx * dy * P[n].v(1);
         } 
      }
   }

   for (int i = 0; i < grid.v_x.nx; i++) {
      for (int j = 0; j < grid.v_y.ny; j++) {
         grid.v_x(i, j) /= grid.mass(i, j);
         grid.v_y(i, j) /= grid.mass(i, j);
      }
   }
}

//STEP 2
void Particles::compute_vol_dens(void) {
   for (int n = 0; n < np; n++) {
      int i = P[n].i;
      int j = P[n].j;

      //computing density
      int low_i = max((int)((P[n].x(0) - 2) / grid.h), 0);
      int low_j = max((int)((P[n].x(1) - 2) / grid.h), 0);
      int high_i = min((int)((P[n].x(0) + 2) / grid.h), grid.mass.nx - 1);
      int high_j = min((int)((P[n].x(1) + 2) / grid.h), grid.mass.ny - 1);

      P[n].d = 0.;
      // iterate for all grid cells within distance 2
      int index = 0;
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            P[n].d += (P[n].weights[index] / (grid.h * grid.h));
            index++;
         }
      }

      P[n].V = (P[n].d > 0.)? (1 / P[n].d) : 0.;
   }
}

void Particles::update_defgrad() {
   for (int n = 0; n < np; n++) {
      //compute velocity gradient
      Eigen::Matrix2d v_grad = Eigen::Matrix2d::Zero();

      int low_i = max((int)((P[n].x(0) - 2) / grid.h), 0);
      int low_j = max((int)((P[n].x(1) - 2) / grid.h), 0);
      int high_i = min((int)((P[n].x(0) + 2) / grid.h), grid.mass.nx - 1);
      int high_j = min((int)((P[n].x(1) + 2) / grid.h), grid.mass.ny - 1);

      // iterate for all grid cells within distance 2
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            v_grad += Eigen::Vector2d(v_star_x(i, j), v_star_y(i, j)) * P[n].grad_weights[index].transpose();
         }
      }

      Eigen::Matrix2d def_star_elas = (Eigen::Matrix2d::Identity() + dt * v_grad) * (P[n].def_elas);
      Eigen::Matrix2d def_star_plas = P[n].def_plas;

      P[n].def_grad = def_star_elas * def_star_plas;
      P[n].update_def_grads();
   }
}

void Particles::transfer_v_to_p(void) {
   for (int n = 0; n < np; n++) {
      Eigen::Vector2d v_pic = Eigen::Vector2d::Zero();
      Eigen::Vector2d v_flip = P[p].v;

      int low_i = max((int)((P[n].x(0) - 2) / grid.h), 0);
      int low_j = max((int)((P[n].x(1) - 2) / grid.h), 0);
      int high_i = min((int)((P[n].x(0) + 2) / grid.h), grid.mass.nx - 1);
      int high_j = min((int)((P[n].x(1) + 2) / grid.h), grid.mass.ny - 1);

      // iterate for all grid cells within distance 2
      index = 0;
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            Eigen::Vector2d v_star = Eigen::Vector2d(grid.v_star_x(i, j), grid.v_star_y(i, j));
            Eigen::Vector2d v = Eigen::Vector2d(grid.v_x(i, j), gird.v_y(i, j));
            v_pic +=  P[n].weights[index] * v_star;
            v_flip += P[n].weights[index] * (v_star - v);
            index++;
         }
      }

      P[n].v = (1 - ALPHA) * v_pic + ALPHA * v_flip;
   }
}

//STEP 9
void Particle::resolve_collisions(void) {

}

void Particle::update_x(void) { 
   for (int n = 0; n < np; n++) {
      P[n].x = P[n].x + dt * P[n].v;
   }
}


void Particles::
update_from_grid(void)
{
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



