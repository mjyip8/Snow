/**
 * Implementation of the MAC grid. Most of your edits should be in particles.cpp instead.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cmath>
#include <vector>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include "grid.h"

using namespace std;

void Grid::
init(float gravity_, int cell_nx, int cell_ny, float lx_)
{
   gravity=gravity_;
   lx=lx_;
   ly=cell_ny*lx/cell_nx;
   h=lx/cell_nx;
   overh=cell_nx/lx;
   // allocate all the grid variables
   grad_weights_x.resize(cell_nx, cell_ny);
   grad_weights_y.resize(cell_nx, cell_ny);
   f_x.resize(cell_nx, cell_ny);
   f_y.resize(cell_nx, cell_ny);
   v_x.resize(cell_nx, cell_ny);
   v_y.resize(cell_nx, cell_ny);
   v_x_n1.resize(cell_nx, cell_ny);
   v_y_n1.resize(cell_nx, cell_ny);
   v_star_x.resize(cell_nx, cell_ny);
   v_star_y.resize(cell_nx, cell_ny);
   mass.resize(cell_nx, cell_ny);

   marker.init(cell_nx, cell_ny);
}

//====================================== private helper functions ============================

void Grid::
init_phi(void)
{
   int i, j;
   // start off with indicator inside the fluid and overestimates of distance outside
   float large_distance=phi.nx+phi.ny+2;
   for(i=0; i<phi.size; ++i)
      phi.data[i]=large_distance;
   for(j=1; j<phi.ny-1; ++j) for(i=1; i<phi.nx-1; ++i){
      if(marker(i,j)==FLUIDCELL){
         phi(i,j)=-0.5;
      }
   }
}

float Grid::get_max(void) {
   float r1 = 0;
   float r2 = 0;
   for (int i = 0; i < v_x.rows(); i++) {
      for (int j = 0; j < v_y.cols(); j++) {
         if (!(abs(v_x(i, j)) <= r1)) r1 = abs(v_x(i, j));
         if (!(abs(v_y(i, j)) <= r2)) r2 = abs(v_y(i, j));
      }
   }

   return r1 * r1 + r2 * r2;
}

float Grid::CFL(void)
{
   double maxv2=max(h*gravity, get_max());
   if(maxv2<1e-16) maxv2=1e-16;
   return h/sqrt(maxv2);
}

float Grid::bspline_weight(float n) {
   float x = abs(n);
   if (x >= 0. && x < 1.) {
      return 0.5 * x * x * x - x * x + (2./3.);
   } else if (x >= 1. && x < 2.) {
      return -(1./6.)*x*x*x + x*x - 2*x + (4./3.);
   } else {
      return 0;
   }
}

float Grid::bspline_gradweight(float n) {
   float x = abs(n);
   if (x >= 0. && x < 1.) {
      return 1.5*x*n - 2*n;
   } else if (x >= 1. && x < 2.) {
      return -.5*x*n + 2*n - 2*n/abs(n);
   } else {
      return 0;
   }
}

//STEP 4
void Grid::update_v(void) {
   v_star_x = Eigen::ArrayXXd::Zero(v_star_x.rows(), v_star_x.cols());
   v_star_y = Eigen::ArrayXXd::Zero(v_star_y.rows(), v_star_y.cols());
   v_x_n1 = Eigen::ArrayXXd::Zero(v_x_n1.rows(), v_x_n1.cols());
   v_y_n1 = Eigen::ArrayXXd::Zero(v_y_n1.rows(), v_y_n1.cols());

   for (int i = 0; i < v_star_x.rows(); i++) {
      for (int j = 0; j < v_star_y.cols(); j++) {
         v_star_x(i, j) = (mass(i, j) == 0)? 0 : v_x(i, j) + dt * (f_x(i, j) / mass(i, j));
         v_star_y(i, j) = (mass(i, j) == 0)? 0 : v_y(i, j) + dt * (-gravity + f_y(i, j) / mass(i, j));
      }
   }
}

//STEP 5
void Grid::resolve_collisions(void) {
   int i, j;
   // first mark where solid is
   for(j=0; j<marker.ny; ++j)
      marker(0,j)=marker(marker.nx-1,j)=SOLIDCELL;
   for(i=0; i<marker.nx; ++i)
      marker(i,0)=marker(i,marker.ny-1)=SOLIDCELL;
      // now make sure nothing leaves the domain

   for (int i = 0; i < marker.nx; ++i) {
      for (int j = 0; j < marker.ny; ++j) {
         if (marker(i, j) == FLUIDCELL) {
            double x_star_x = i * h + v_star_x(i, j) * dt; 
            double x_star_y = j * h + v_star_y(i, j) * dt; 

            if (x_star_x < 2 * h || x_star_x > 1 - 2 * h) {
               v_star_x(i, j) = 0;
               v_star_y(i, j) *= FRICTION;
            }

            if (x_star_y < 2 * h || x_star_y > 1 - 2 * h) {
               v_star_y(i, j) = 0;
               v_star_x(i, j) *= FRICTION;
            }
         }
      }
   }
}


