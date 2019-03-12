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
   v_star_x.resize(cell_nx, cell_ny);
   v_star_y.resize(cell_nx, cell_ny);
   mass.resize(cell_nx, cell_ny);

   marker.init(cell_nx, cell_ny);
}

float Grid::
CFL(void)
{
   float maxv2=max(h*gravity, sqr(u.infnorm())+sqr(v.infnorm()));
   if(maxv2<1e-16) maxv2=1e-16;
   return h/sqrt(maxv2);
}

void Grid::
save_velocities(void)
{
   u.copy_to(du);
   v.copy_to(dv);
}

/* centered gravity is the spherical gravity I added. */
/* for the normal uniform gravity, see lines 107-108 */
void Grid::
add_gravity(float dt, bool centered, float cx, float cy)
{
   float dtg=dt*gravity;
   if( centered )
   {
      for(int i = 0; i < u.nx; ++i)
      {
         for( int j = 0; j < v.ny; ++j)
         {
            float x, y, dx, dy, dr2, dr;
            if ( j < u.ny) 
            {
               x = ( i ) * h;
               y = ( 0.5 + j ) * h;
               dx = x - cx;
               dy = y - cy;
               dr2 = dx * dx + dy * dy;
               if (dr2 < 0.0001 * h * h){
                  printf("dr2 too small: %f \n", dr2);
                  dr2 = 0.0001 * h * h;
               }
               dr = sqrt(dr2);
               dx /= dr;
               dy /= dr;
               u(i,j) -= dtg * dx / dr2;
            }
            if( i < v.nx )
            {
               x = ( 0.5 + i ) * h;
               y = ( j ) * h;
               dx = x - cx;
               dy = y - cy;
               dr2 = dx * dx + dy * dy;
               if (dr2 < 0.0001 * h * h){
                  printf("dr2 too small: %f \n", dr2);
                  dr2 = 0.0001 * h * h;
               }
               dr = sqrt(dr2);
               dx /= dr;
               dy /= dr;
               v(i,j) -= dtg * dy / dr2;
            }

         }
      }
   }
   else
   {
      for(int i=0; i<v.size; ++i)
         v.data[i]-=dtg;
   }
}

void Grid::
compute_distance_to_fluid(void)
{
   init_phi();
   for(int i=0; i<2; ++i)
      sweep_phi();
}

void Grid::
extend_velocity(void)
{
   for(int i=0; i<4; ++i)
      sweep_velocity();
}

void Grid::
apply_boundary_conditions(void)
{
   int i, j;
   // first mark where solid is
   for(j=0; j<marker.ny; ++j)
      marker(0,j)=marker(marker.nx-1,j)=SOLIDCELL;
   for(i=0; i<marker.nx; ++i)
      marker(i,0)=marker(i,marker.ny-1)=SOLIDCELL;
   // now makre sure nothing leaves the domain
   for(j=0; j<u.ny; ++j)
      u(0,j)=u(1,j)=u(u.nx-1,j)=u(u.nx-2,j)=0;
   for(i=0; i<v.nx; ++i)
      v(i,0)=v(i,1)=v(i,v.ny-1)=v(i,v.ny-2)=0;
}

void Grid::
make_incompressible(void)
{
   find_divergence();
   form_poisson();
   form_preconditioner();
   solve_pressure(100, 1e-5);
   add_gradient();
}

void Grid::
get_velocity_update(void)
{
   int i;
   for(i=0; i<u.size; ++i)
      du.data[i]=u.data[i]-du.data[i];
   for(i=0; i<v.size; ++i)
      dv.data[i]=v.data[i]-dv.data[i];
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

float Grid::bspline_weight(float x) {
   if (x >= 0. && x < h * 1.) {
      return 0.5 * x * x * x - x * x + (2./3.);
   } else if (x >= h * 1. && x < h * 2) {
      return -(1./6.)*x*x*x + x*x - 2*x + (4./3.);
   } else {
      return 0;
   }
}

float Grid::bspline_gradweight(float x) {
   if (x >= 0. && x < h * 1.) {
      return 1.5*x*x - 2*x;
   } else if (x >= h * 1. && x < h * 2) {
      return -.5*x*x + 2*x - 2;
   } else {
      return 0;
   }
}

//STEP 3
void Grid::compute_grid_forces(vector<Particle>& P) {
   int np = P.size();

   for (int n = 0; n < np; n++) {
      int low_i = max((int)((P[n].x(0) - 2) / h), 0);
      int low_j = max((int)((P[n].x(1) - 2) / h), 0);
      int high_i = min((int)((P[n].x(0) + 2) / h), mass.nx - 1);
      int high_j = min((int)((P[n].x(1) + 2) / h), mass.ny - 1);

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

            f_x(i, j) += df(0);
            f_y(i, j) += df(1);

            index++;
         }
      }   
   }
}

//STEP 4
void update_v(void) {
   v_star_x = v_x + dt * (1/mass) * f_x;
   v_star_y = v_y + dt * (1/mass) * f_y;
}

//STEP 5
void resolve_collisions(void) {
   int i, j;
   // first mark where solid is
   for(j=0; j<marker.ny; ++j)
      marker(0,j)=marker(marker.nx-1,j)=SOLIDCELL;
   for(i=0; i<marker.nx; ++i)
      marker(i,0)=marker(i,marker.ny-1)=SOLIDCELL;
      // now makre sure nothing leaves the domain

   for (int i = 0; i < marker.nx; ++i) {
      for (int j = 0; j < marker.ny; ++j) {
         if (marker(i, j) == FLUIDCELL) {
            double x_star_x = v_star_x(i. j) * dt / h; 
            double x_star_y = v_star_y(i. j) * dt / h; 

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

