#include "array2.h"
#include <Eigen/Dense>

#ifndef PARTICLE
#define PARTICLE
#define YOUNG_MOD (0.000014)
#define POSS_R (0.2)
#define CRIT_STRETCH (1 + 0.0075)
#define CRIT_COMPRESS (1 - 0.025)
#define HARDEN (5.0)

struct Particle {
	Eigen::Vector2d x;
	Eigen::Vector2d v; // positions and velocities
	Eigen::Vector2d x_star; //deformed position
 
	double d;
	double V;
	Eigen::Matrix2d grad_v_n1;

	Eigen::Matrix2d def_grad;
	Eigen::Matrix2d def_elas;
	Eigen::Matrix2d def_plas;
	Eigen::Matrix2d rot_mat;
	Eigen::Matrix2d str_mat;

	Eigen::Matrix2d svd_v;
	Eigen::Matrix2d svd_s;
	Eigen::Matrix2d svd_u;

	std::vector<Eigen::Vector2d> grad_weights;
	std::vector<double> weights;


	double mu;
	double lambda;

	int i;
	int j;

	Particle(Eigen::Vector2d X, Eigen::Vector2d U): x(X), v(U), 
	d(0.), V(0.), i(0), j(0) {
		grad_weights.clear();
		//setting as identity matrix
		grad_v_n1 = Eigen::Matrix2d::Identity();
		def_grad = Eigen::Matrix2d::Identity();
		def_elas = Eigen::Matrix2d::Identity();
		def_plas = Eigen::Matrix2d::Identity();
		rot_mat = Eigen::Matrix2d::Identity();
		str_mat = Eigen::Matrix2d::Identity();

		svd_v = Eigen::Matrix2d::Identity();
		svd_s = Eigen::Matrix2d::Identity();
		svd_u = Eigen::Matrix2d::Identity();

		lambda = (YOUNG_MOD * POSS_R)/((1 + POSS_R) * (1 - 2 * POSS_R));
		mu = (YOUNG_MOD)/ (2 * (1 + POSS_R));
	}

	void update_def_grads(double dt);
	Eigen::Matrix2d get_energy_deriv(void);

	private:
		void clamp_def_elas(void);
};
#endif
