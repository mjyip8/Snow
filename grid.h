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

   // active variables
   //Array2f u, v; // staggered MAC grid of velocities
   //Array2f du, dv; // saved velocities and differences for particle update
   Array2c marker; // identifies what sort of cell we have
   Array2f phi; // decays away from water into air (used for extrapolating velocity)
   Array2d pressure;

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

   float bspline_weight(float n);
   float bspline_gradweight(float n);

   void update_v(void);
   void resolve_collisions(void);

   private:
      float get_max(void);
      void init_phi(void);
};

#endif

