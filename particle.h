#include "array2.h"
#include "vec2.h"
#include <Eigen/Dense>

#ifndef PARTICLE
#define PARTICLE

struct Particle {
	Vec2f x;
	Vec2f u; // positions and velocities

	Vec2f x_star; //deformed position
   
	Vec2f cx;
	Vec2f cy;

	double d;
	double V;
	Eigen::Vector2d grad_w;
	Eigen::Matrix2d def_grad;
	Eigen::Matrix2d def_elas;
	Eigen::Matrix2d def_plas;
	Eigen::Matrix2d rot_mat;
	Eigen::Matrix2d scal_mat;

	Eigen::Matrix2d svd_v;
	Eigen::Matrix2d svd_s;
	Eigen::Matrix2d svd_u;

	std::vector<Eigen::Vector2i> weight_nodes;
	std::vector<Eigen::Vector2d> grad_weights;

	int i;
	int j;

	Particle(Vec2f X, Vec2f U): x(X), u(U), cx(Vec2f(0.f, 0.f)), cy(Vec2f(0.f, 0.f)),
	d(0.), V(0.), i(0), j(0) {
		weight_nodes.clear();
		grad_weights.clear();
		//setting as identity matrix
		def_grad = Eigen::Matrix2d::Zero();
		rot_mat = Eigen::Matrix2d::Identity();
		scal_mat = Eigen::Matrix2d::Identity();
		def_elas = Eigen::Matrix2d::Identity();
		def_plas = Eigen::Matrix2d::Identity();
		grad_w = Eigen::Vector2d::Zero();

		svd_v = Eigen::Matrix2d::Identity();
		svd_s = Eigen::Matrix2d::Identity();
		svd_u = Eigen::Matrix2d::Identity();
	}

	void updateDeformGradients(void);
};
#endif
