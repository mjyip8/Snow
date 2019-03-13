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
add_particle(const Eigen::Vector2d &px, const Eigen::Vector2d &pu)
{
   Particle p = Particle(px, pu);
   P.push_back(p);
   ++np;
}

// STEP 1
void Particles::transfer_mass_to_grid(void) {
   grid.marker.zero();
   grid.mass = Eigen::ArrayXXd::Zero(grid.mass.rows(), grid.mass.cols());
   grid.grad_weights_x = Eigen::ArrayXXd::Zero(grid.grad_weights_x.rows(), grid.grad_weights_x.cols());
   grid.grad_weights_y = Eigen::ArrayXXd::Zero(grid.grad_weights_y.rows(), grid.grad_weights_y.cols());

   for (int n = 0; n < np; n++) {
      P[n].i = (int) (P[n].x(0) / grid.h);
      P[n].j = (int) (P[n].x(1) / grid.h);

      int low_i = max(P[n].i - 2, 0);
      int low_j = max(P[n].j - 2, 0);
      int high_i = min(P[n].i + 2, (int) grid.mass.rows() - 1);
      int high_j = min(P[n].j + 2, (int) grid.mass.cols() - 1);

      // iterate for all grid cells within distance 2
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {

            float wx = grid.bspline_weight((P[n].x[0] - i * grid.h)/ grid.h);
            float wy = grid.bspline_weight((P[n].x[1] - j * grid.h) / grid.h);
            grid.mass(i, j) += (wx * wy);
            P[n].weights.push_back(wx * wy);

            Eigen::Vector2d grad_weight = Eigen::Vector2d(wx * grid.bspline_gradweight((P[n].x[0] - i * grid.h) / grid.h),
                                       wy * grid.bspline_gradweight((P[n].x[1] - j * grid.h) / grid.h));
            cout << grad_weight << endl;
            P[n].grad_weights.push_back(grad_weight);
            grid.grad_weights_x(i, j) += grad_weight(0);
            grid.grad_weights_y(i, j) += grad_weight(1);
            grid.marker(i, j) = FLUIDCELL;
            
         }
      }
   }
}

//STEP 1
void Particles::transfer_v_to_grid(void) {
   grid.v_x = Eigen::ArrayXXd::Zero(grid.v_x.rows(), grid.v_x.cols());
   grid.v_y = Eigen::ArrayXXd::Zero(grid.v_y.rows(), grid.v_y.cols());

   for (int n = 0; n < np; n++) {
      int low_i = max(P[n].i - 2, 0);
      int low_j = max(P[n].j - 2, 0);
      int high_i = min(P[n].i + 2, (int) grid.mass.rows() - 1);
      int high_j = min(P[n].j + 2, (int) grid.mass.cols() - 1);

      for (int i = low_i ; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            double dx = abs(P[n].x(0) - (i * grid.h)) / grid.h;
            double dy = abs(P[n].x(1) - (j * grid.h)) / grid.h;

            grid.v_x(i, j) += dx * dy * P[n].v(0);
            grid.v_y(i, j) += dx * dy * P[n].v(1);
         } 
      }
   }

   for (int i = 0; i < grid.v_x.rows(); i++) {
      for (int j = 0; j < grid.v_y.cols(); j++) {
         grid.v_x(i, j) = (grid.mass(i, j) == 0) ? 0 : grid.v_x(i, j) / grid.mass(i, j);
         grid.v_y(i, j) = (grid.mass(i, j) == 0) ? 0 : grid.v_y(i, j) / grid.mass(i, j);
         //cout << "grid m = " << grid.mass(i, j) << endl;
         //cout << "grad_weights = (" << grid.grad_weights_x(i, j) << ", " << grid.grad_weights_y(i, j) << ")" << endl; 
         //cout << "grid v = (" << grid.v_x(i, j) << ", " << grid.v_y(i, j) << ")" << endl;
      }
   }
}

//STEP 2
void Particles::compute_vol_dens(void) {
   for (int n = 0; n < np; n++) {
      //computing density
      int low_i = max(P[n].i - 2, 0);
      int low_j = max(P[n].j - 2, 0);
      int high_i = min(P[n].i + 2, (int) grid.mass.rows() - 1);
      int high_j = min(P[n].j + 2, (int) grid.mass.cols() - 1);

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

//STEP 3
void Particles::compute_grid_forces(void) {
   cout << "entering compute_grid_forces" << endl;
   grid.f_x = Eigen::ArrayXXd::Zero(grid.f_x.rows(), grid.f_x.cols());
   grid.f_y = Eigen::ArrayXXd::Zero(grid.f_y.rows(), grid.f_y.cols());

   int np = P.size();

   for (int n = 0; n < np; n++) {
      int low_i = max(P[n].i - 2, 0);
      int low_j = max(P[n].j - 2, 0);
      int high_i = min(P[n].i + 2, (int) grid.mass.rows() - 1);
      int high_j = min(P[n].j + 2, (int) grid.mass.cols() - 1);

      int index = 0;
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            double over_Jp = 0.;
            if (P[n].def_plas.determinant() != 0.) {
               over_Jp = 1.0 / P[n].def_plas.determinant();
            }            
            Eigen::Matrix2d d_energy = P[n].get_energy_deriv();
            Eigen::Vector2d gw = P[n].grad_weights[index];
            Eigen::Vector2d df = P[n].V * over_Jp * d_energy * P[n].def_elas.transpose() * gw;
            grid.f_x(i, j) += df(0);
            grid.f_y(i, j) += df(1);

            index++;
         }
      }
   }
   cout << "leaving compute_grid_forces" << endl;

}

void Particles::update_defgrad() {
   for (int n = 0; n < np; n++) {
      //compute velocity gradient
      Eigen::Matrix2d v_grad = Eigen::Matrix2d::Zero();

      int low_i = max(P[n].i - 2, 0);
      int low_j = max(P[n].j - 2, 0);
      int high_i = min(P[n].i + 2, (int) grid.mass.rows() - 1);
      int high_j = min(P[n].j + 2, (int) grid.mass.cols() - 1);

      int index = 0;
      // iterate for all grid cells within distance 2
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            v_grad += Eigen::Vector2d(grid.v_star_x(i, j), grid.v_star_y(i, j)) * P[n].grad_weights[index].transpose();
            index++;
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
      Eigen::Vector2d v_flip = P[n].v;

      int low_i = max(P[n].i - 2, 0);
      int low_j = max(P[n].j - 2, 0);
      int high_i = min(P[n].i + 2, (int) grid.mass.rows() - 1);
      int high_j = min(P[n].j + 2, (int) grid.mass.cols() - 1);

      // iterate for all grid cells within distance 2
      int index = 0;
      for (int i = low_i; i <= high_i; i++) {
         for (int j = low_j; j <= high_j; j++) {
            Eigen::Vector2d v_star = Eigen::Vector2d(grid.v_star_x(i, j), grid.v_star_y(i, j));
            Eigen::Vector2d v = Eigen::Vector2d(grid.v_x(i, j), grid.v_y(i, j));
            v_pic +=  P[n].weights[index] * v_star;
            v_flip += P[n].weights[index] * (v_star - v);
            index++;
         }
      }

      P[n].v = (1 - ALPHA) * v_pic + ALPHA * v_flip;
   }
}

//STEP 9
void Particles::resolve_collisions(void) {
   for (int n = 0; n < np; n++) {
      Eigen::Vector2d x_star = P[n].x + dt * P[n].v;

      if (x_star(0) < 2 * grid.h || x_star(0) > 1 - 2 * grid.h) {
         P[n].v(0) = FRICTION * P[n].v(0);
      }

      if (x_star(1) < 2 * grid.h || x_star(1) > 1 - 2 * grid.h) {
         P[n].v(1) = FRICTION * P[n].v(1);
      }
   }
}

void Particles::update_x(void) { 
   for (int n = 0; n < np; n++) {
      P[n].x = P[n].x + dt * P[n].v;
      P[n].x(0) = clamp(P[n].x(0), 0.0, 1.0);
      P[n].x(1) = clamp(P[n].x(1), 0.0, 1.0);
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



