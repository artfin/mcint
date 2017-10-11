#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0;
const double mu2 = 36440.0;
const double l = 4.398;

void fill_inertia_tensor(Matrix3d &inertia_tensor, double &R, double &theta)
{
    double sin_t = sin(theta);
    double cos_t = cos(theta);
    double l2 = l * l;
    double R2 = R * R;

    inertia_tensor(0, 0) = mu1 * l2 * cos_t * cos_t + mu2 * R2; 
    inertia_tensor(1, 0) = 0;
    inertia_tensor(2, 0) = -mu1 * l2 * sin_t * cos_t;

    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = 0;
    inertia_tensor(2, 2) = mu1 * l2 * sin_t * sin_t;

    inertia_tensor(0, 1) = 0;
    inertia_tensor(1, 1) = inertia_tensor(0, 0) + inertia_tensor(2, 2);
    inertia_tensor(2, 1) = 0;
}

void fill_a_matrix(Matrix2d &a, double &R, double &theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;
}

void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;
}

// ham_vars = [R, theta, pR, pT, Jx, Jy, Jz]
// omega = [omega_x, omega_y, omega_z]
// q_dot = [ dot{r}, dot{theta} ]
void lagrange_vars( VectorXf ham_vars, Vector3d &omega, Vector2d &q_dot )
{
	Vector3d j_vector( ham_vars(4), ham_vars(5), ham_vars(6) );
	Vector2d p_vector( ham_vars(2), ham_vars(3) );

	for ( int i = 0; i < 3; i++ )
	{
		cout << "j_vector[" << i << "] = " << j_vector(i) << endl;
	}

	//double R = ham_vars(0);
	//double Theta = ham_vars(1);

	//Matrix<double, 3, 3> I;
	//fill_inertia_tensor( I, R, Theta );
	//Matrix<double, 3, 3> I_inv = I.inverse();

	//Matrix<double, 2, 2> a;
	//fill_a_matrix( a, R, Theta );
	//Matrix<double, 2, 2> a_inv = a.inverse();

	//Matrix<double, 3, 2> A;
	//fill_A_matrix( A, R, Theta );

	//Matrix<double, 3, 3> t1 = I;
	//t1.noalias() -= A * a_inv * A.transpose();
	//Matrix<double, 3, 3> G11 = t1.inverse();

	//Matrix<double, 2, 2> t2 = a;
	//t2.noalias() -= A.transpose() * I_inv * A;
	//Matrix<double, 2, 2> G22 = t2.inverse();

	//Matrix<double, 3, 2> G12 = - G11 * A * a.inverse();
	//Matrix<double, 2, 3> G21 = G12.transpose();

	//omega = G11 * j_vector + G12 * p_vector;
	//q_dot = G21 * j_vector + G22 * p_vector;
}

double kinetic_energy(double R, double theta, double pR, double pT, double Jx, double Jy, double Jz)
{
	Vector3d j_vector(Jx, Jy, Jz);
	Vector2d p_vector(pR, pT);

	Matrix<double, 3, 3> I;
	Matrix<double, 2, 2> a;
	Matrix<double, 3, 2> A;

	fill_inertia_tensor(I, R, theta);
	fill_a_matrix(a, R, theta);
	fill_A_matrix(A, R, theta);

	Matrix<double, 3, 3> I_inv = I.inverse();
	Matrix<double, 2, 2> a_inv = a.inverse();

	Matrix<double, 3, 3> G11;
	Matrix<double, 2, 2> G22;
	Matrix<double, 3, 2> G12;

	Matrix<double, 3, 3> t1;
	Matrix<double, 2, 2> t2;

	t1 = I;
	t1.noalias() -= A * a_inv * A.transpose();
	G11 = t1.inverse();

	t2 = a;
	t2.noalias() -= A.transpose() * I_inv * A;
	G22 = t2.inverse();

	G12.noalias() = - G11 * A * a.inverse();

	double ang_term = 0.5 * j_vector.transpose() * G11 * j_vector;
	double kin_term = 0.5 * p_vector.transpose() * G22 * p_vector;
	double cor_term = j_vector.transpose() * G12 * p_vector;

	return ang_term + kin_term + cor_term; 
}






