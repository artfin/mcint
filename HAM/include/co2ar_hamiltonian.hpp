#include <math.h>
#include <Eigen/Dense>

using Eigen::Matrix3d;
using Eigen::Matrix2d;
using Eigen::Matrix;

void fill_inertia_tensor(Matrix3d &inertia_tensor, double &R, double &theta);
void fill_a_matrix(Matrix2d &a, double &R, double &theta);
void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &theta);

double kinetic_energy(double R, double theta, double pR, double pT, double Jx, double Jy, double Jz);
void lagrange_vars( Eigen::VectorXf ham_vars, Eigen::Vector3d &omega, Eigen::Vector2d &q_dot );
