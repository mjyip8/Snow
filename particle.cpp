#include "array2.h"
#include "vec2.h"
#include "particle.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;

void Particle::update_def_grads(void) {
	Eigen::JacobiSVD<Eigen::Matrix2d> svd(def_elas, Eigen::ComputeFullV | Eigen::ComputeFullU);
	svd_v = svd.matrixV();
	Eigen::Vector2d singVal = svd.singularValues();
	singVal(0) = clamp(singVal(0), CRIT_COMPRESS, CRIT_STRETCH);
	singVal(1) = clamp(singVal(1), CRIT_COMPRESS, CRIT_STRETCH);
	svd_s = singVal.asDiagonal();
	svd_u = svd.matrixU();

	def_elas = svd_u * svd_s * svd_v.transpose();
	def_plas = svd_v * svd_s.inverse() * svd_u.transpose() * def_plas;

	str_mat = svd_v * svd_s * svd_v.transpose();
	rot_mat = svd_u * svd_v.transpose();
}

Eigen::Matrix2d Particle::get_energy_deriv(void) {
	Eigen::Matrix2d result = 2 * mu * (def_elas - rot_mat);
	result += lambda * (def_elas.determinant() - 1) * def_elas.determinant() * def_elas.transpose().inverse();
	return result;
}
