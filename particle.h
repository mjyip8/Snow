#include "array2.h"
#include "vec2.h"

struct Particle {
	Vec2f x;
	Vec2f u; // positions and velocities
   
	Vec2f cx;
	Vec2f cy;

	double d;
	double V;
	Array2f def_grad;
	Array2f def_elas;
	Array2f def_plas;
	Array2f rot_mat;
	Array2f scal_mat;
	int i;
	int j;

	Particle(Vec2f X, Vec2f U): x(X), u(U), cx(Vec2f(0.f, 0.f)), cy(Vec2f(0.f, 0.f)),
	d(0.), V(0.), def_grad(Array2f()), def_elas(Array2f()), def_plas(Array2f()), 
	rot_mat(Array2f()), scal_mat(Array2f()), i(0), j(0) {

	}
};
