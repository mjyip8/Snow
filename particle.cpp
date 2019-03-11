#include "array2.h"
#include "vec2.h"
#include "particle.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;

void Particle::updateDeformGradients(void) {
	Eigen::Jacobi<Matrix2d> svd(def_elastic, ComputeFullV | ComputeFullU);
	svd_v = svd.matrixV();
	svd_s = svd.singularValues().asDiagonal();
	svd_u = svd.matrixU();
}

Eigen::Matrix2d Particle::get_energy_deriv(void) {
	Matrix2d result = 2 * mu * (def_elas - rot_mat);
	result += lambda * (def_elas.determinant() - 1) * def_elas.determinant() * def_elas.transpose();
	return result;
}
