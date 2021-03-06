/**
 * Declaration for the Particles struct, which holds the particles.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#ifndef PARTICLES_H
#define PARTICLES_H
#define ALPHA 0.95

#include <vector>
#include "grid.h"
#include "vec2.h"
#include "particle.h"
#include "array2.h"

typedef enum SimulationTypeEnum { PIC = 0, FLIP = 1, APIC = 2 } SimulationType;

struct Particles{
   Grid &grid;
   int np; // number of particles
   double dt;

   std::vector<Particle> P;

   // transfer stuff
   Array2f sum_u;
   Array2f sum_v;
   Array2f mass_weight;
   SimulationType simType;

   Particles(Grid &grid_, SimulationType simType_)
      :grid(grid_), np(0), mass_weight(grid_.pressure.nx+1, grid_.pressure.ny+1),
       sum_u(grid_.pressure.nx+1, grid_.pressure.ny+1), sum_v(grid_.pressure.nx+1, grid_.pressure.ny+1), simType( simType_ )
   {}

   void add_particle(const Eigen::Vector2d &px, const Eigen::Vector2d &pu);
   void write_to_file(const char *filename_format, ...);

   void compute_grid_forces();
   void compute_vol_dens(void);
   void transfer_mass_to_grid(void);
   void transfer_v_to_grid(void);
   void update_vn1(void);
   void update_defgrad(void);
   void transfer_v_to_p(void);
   void resolve_collisions(void);
   void update_x(void);


   private:
   template<class T> void accumulate(T &accum, float q, int i, int j, float fx, float fy, Array2f &sum);
   Vec2f computeC(Array2f &ufield, int i, int j, float fx, float fy);
};

#endif
