#pragma once

#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

void inertia_tensor(Matrix<double, 3, 3> &inertia_tensor, const double &R, const double &Theta);
void a_matrix(Matrix<double, 2, 2> &a, const double &R, const double &Theta);
void A_matrix(Matrix<double, 3, 2> &A, const double &R, const double &Theta);
void W_matrix( Matrix<double, 3, 3> &W, const double &theta, const double &psi );

double kinetic_energy( const double &R, const double &Theta, const double &pR, const double &pT, const double &Jx, const double &Jy, const double &Jz );
double kinetic_energy_euler( const double &R, const double &Theta, const double &pR, const double &pT, const double &phi, const double &theta, const double &psi, const double &p_phi, const double &p_theta, const double &p_psi );

