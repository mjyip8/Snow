#include "array2.h"
#include "vec2.h"
#include "particle.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;

void Particle::update_def_grads(double dt) {
    def_elas = (Eigen::Matrix2d::Identity() + dt * grad_v_n1) * def_elas;
	def_grad = def_elas * def_plas;

	Eigen::JacobiSVD<Eigen::Matrix2d> svd(def_elas, Eigen::ComputeFullV | Eigen::ComputeFullU);
	svd_v = svd.matrixV();
	Eigen::Vector2d singVal = svd.singularValues();

	singVal(0) = clamp(singVal(0), CRIT_COMPRESS, CRIT_STRETCH);
	singVal(1) = clamp(singVal(1), CRIT_COMPRESS, CRIT_STRETCH);
	svd_s = singVal.asDiagonal();
	svd_u = svd.matrixU();

	def_elas = svd_u * svd_s * svd_v.transpose();
	def_plas = svd_v * svd_s.inverse() * svd_u.transpose() * def_grad;

	rot_mat = svd_u * svd_v.transpose();
	str_mat = svd_v * svd_s * svd_v.transpose();
}

/*Eigen::Matrix2d Particle::get_energy_deriv(void) {
	double Je = def_elas.determinant();
	//double harden = exp(HARDEN * (1 - def_plas.determinant()));
	Eigen::Matrix2d result = 2 * mu * (def_elas - rot_mat) * def_elas.transpose();
	result += lambda * (Je - 1) * Je * Eigen::Matrix2d::Identity();
	//return harden * result;
	return result;
}*/

Eigen::Matrix2d Particle::get_energy_deriv(void) {
	double Je = def_elas.determinant();
	double harden = exp(HARDEN * (1 - def_plas.determinant()));
	Eigen::Matrix2d result = 2 * mu * (def_elas - svd_u * svd_v.transpose()) * def_elas.transpose();
	result += lambda * (Je - 1) * Je * Eigen::Matrix2d::Identity();
	return harden * result;
	//return result;
}
