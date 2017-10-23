#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double mu2 = 0.5;
const double mu3 = 0.5;
const double mu4 = 0.25;
const double mu1 = 1960/209;
const double l1 = 1.0;
const double l2 = sqrt(2);
//=============== q[1] = R, q[2] = theta, q[3] = phi==============//

void fill_inertia_tensor(Matrix<double, 3, 3> &inertia_tensor, double &q1, double &q2, double &q3)
{
	double cos_q2 = cos(q2);
	double cos_q3 = cos(q3);
	double sin_q2 = sin(q2);
	double sin_q3 = sin(q3);
	double q1_sq = q1 * q1;
	double c2sq = cos_q2 * cos_q2;
	double c3sq = cos_q3 * cos_q3;

	//inertia_tensor(0, 0) = (3/2)*cos(q2)*cos(q2)*cos(q3)*cos(q3)-cos(q3)*cos(q3)-cos(q2)*cos(q2)+5/2+(1960/209)*q1*q1;
    inertia_tensor(0, 0) = 1.5 * c2sq * c3sq - c3sq - c2sq + 2.5 + mu1 * q1_sq;

	inertia_tensor(1, 0) = 0.5*cos_q3*sin_q3*( c2sq*sqrt(2.0) + c2sq - sqrt(2.0) );
    inertia_tensor(2, 0) = -0.5*cos_q2*sin_q2*( cos_q3*sqrt(2.0) - cos_q3 + 2.0 * sin_q3 );

    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = 0.5*cos_q2*sin_q2*( -sin_q3 + 2.0*cos_q3 );
    inertia_tensor(2, 2) = 0.5*c2sq*c3sq - 0.5*c3sq + 0.5*c2sq +1.0;

    inertia_tensor(0, 1) = inertia_tensor(1,0);
    inertia_tensor(1, 1) = -c2sq*c3sq + 0.5*c3sq + 0.5*c2sq + 1.5 + (1960/209)*q1_sq;
    inertia_tensor(2, 1) = inertia_tensor(1,2);
}
void fill_a_matrix(Matrix3d &a, double &q1, double &q2, double &q3)
{
    double cos_q2 = cos(q2);
	double cos_q3 = cos(q3);
	double sin_q2 = sin(q2);
	double sin_q3 = sin(q3);
	double c2sq = cos_q2 * cos_q2;
	double c3sq = cos_q3 * cos_q3;

    a(0, 0) = mu1;
    a(0, 1) = 0;
    a(0, 2) = 0;
    a(1, 0) = 0;
    a(2, 0) = 0;
    a(1, 1) = -0.5*c2sq*c3sq + 2.5;
    a(1, 2) = 0.5*cos_q3*cos_q2*sin_q3*sin_q2;
    a(2, 1) = a(1, 2);
    a(2, 2) = -0.5*c2sq*c3sq + 0.5*c3sq + c2sq + 0.5;
}

void fill_A_matrix(Matrix3d &A, double &q1, double &q2, double &q3)
{
    double cos_q2 = cos(q2);
	double cos_q3 = cos(q3);
	double sin_q2 = sin(q2);
	double sin_q3 = sin(q3);
	double c2sq = cos_q2 * cos_q2;

    A(0, 0) = A(1, 0) = A(2, 0) = A(2, 1) = 0;
    A(0, 1) = -cos_q3 - 1.5*sin_q3;
    A(1, 1) = 0.5 * cos_q3 * sqrt(2) + 0.5 * cos_q3 - sin_q3;
    A(0, 2) = -0.5 * cos_q2 * sin_q2 * (cos_q3 + 2.0 * sin_q3);
    A(1, 2) = 0.5 * cos_q2 * sin_q2 * (- sin_q3 * sqrt(2.0) + 2.0 * cos_q3 + sin_q3);
    A(2, 2) = -0.5 * c2sq * sqrt(2.0) + 0.5 * sqrt(2.0) + 1.5 * c2sq;
}

double kinetic_energy(double q1, double q2, double q3, double p1, double p2, double p3, double Jx, double Jy, double Jz)
{
	Vector3d j_vector(Jx, Jy, Jz);
	Vector3d p_vector(p1, p2, p3);

	Matrix<double, 3, 3> I;
	Matrix<double, 3, 3> a;
	Matrix<double, 3, 3> A;

	fill_inertia_tensor(I, q1, q2, q3);
	cout << "I: " << endl << I << endl;

	fill_a_matrix(a, q1, q2, q3);
	cout << "a: " << endl << a << endl;

	fill_A_matrix(A, q1, q2, q3);
	cout << "A: " << endl << A << endl;

	Matrix<double, 3, 3> I_inv = I.inverse();
	Matrix<double, 3, 3> a_inv = a.inverse();

	Matrix<double, 3, 3> G11;
	Matrix<double, 3, 3> G22;
	Matrix<double, 3, 3> G12;

	Matrix<double, 3, 3> t1;
	Matrix<double, 3, 3> t2;

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


int main()
{
	double ke = kinetic_energy( 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 2.0, 2.0, 2.0 );
	//fill_inertia_tensor( I, 5.0, 0.25, 0.30 );
	cout << "ke: " << endl << ke << endl;
	return 0;
}



