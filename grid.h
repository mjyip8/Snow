/**
 * Definition of the Grid struct, which holds the MAC grid.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#ifndef GRID_H
#define GRID_H

#include "array2.h"
#include "vec2.h"
#include "util.h"
#include "particles.h"
#include <map>
#include <Eigen/Dense>

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2
#define FRICTION 0.8

using namespace std;

struct Grid{
   float gravity;
   float lx, ly;
   float h, overh;
   double dt;

   Array2c marker; // identifies what sort of cell we have

   Eigen::ArrayXXd mass;

   //stuff for the elastoplastic solve
   Eigen::ArrayXXd grad_weights_x;
   Eigen::ArrayXXd grad_weights_y;

   //for computing weights
   Eigen::ArrayXXd f_x;
   Eigen::ArrayXXd f_y;

   //for computing velocity
   Eigen::ArrayXXd v_x;
   Eigen::ArrayXXd v_y;

   Eigen::ArrayXXd v_star_x;
   Eigen::ArrayXXd v_star_y;

   Grid(void)
   {}

   Grid(float gravity_, int cell_nx, int cell_ny, float lx_)
   { 
      init(gravity_, cell_nx, cell_ny, lx_); 
   }

   void init(float gravity_, int cell_nx, int cell_ny, float lx_);
   float CFL(void);
   void save_velocities(void);
   void add_gravity(float dt, bool centered, float cx, float cy);
   void compute_distance_to_fluid(void);
   void extend_velocity(void);
   void apply_boundary_conditions(void);
   void make_incompressible(void);
   void get_velocity_update(void);

   void bary_x(float x, int &i, float &fx)
   {
      float sx=x*overh;
      i=(int)sx;
      fx=sx-floor(sx);
   }

   void bary_x_centre(float x, int &i, float &fx)
   {
      float sx=x*overh-0.5;
      i=(int)sx;
      if(i<0){ i=0; fx=0.0; }
      else if(i>pressure.nx-2){ i=pressure.nx-2; fx=1.0; }
      else{ fx=sx-floor(sx); }
   }

   void bary_y(float y, int &j, float &fy)
   {
      float sy=y*overh;
      j=(int)sy;
      fy=sy-floor(sy);
   }

   void bary_y_centre(float y, int &j, float &fy)
   {
      float sy=y*overh-0.5;
      j=(int)sy;
      if(j<0){ j=0; fy=0.0; }
      else if(j>pressure.ny-2){ j=pressure.ny-2; fy=1.0; }
      else{ fy=sy-floor(sy); }
   }

   void bilerp_uv(float px, float py, float &pu, float &pv)
   {
      int i, j;
      float fx, fy;
      bary_x(px, i, fx);
      bary_y_centre(py, j, fy);
      pu=u.bilerp(i, j, fx, fy);
      bary_x_centre(px, i, fx);
      bary_y(py, j, fy);
      pv=v.bilerp(i, j, fx, fy);
   }


   float bspline_weight(float x);
   float bspline_gradweight(float x);

   void update_v(void);
   void resolve_collisions(void);

   private:

};

#endif

