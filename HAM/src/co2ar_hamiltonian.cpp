#include "co2ar_hamiltonian.hpp"

const double mu1 = 14579.0;
const double mu2 = 36440.0;
const double l = 4.398;

void inertia_tensor(Matrix<double, 3, 3> &inertia_tensor, const double &R, const double &Theta)
{
    double sin_t = sin(Theta);
    double cos_t = cos(Theta);
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

void a_matrix(Matrix<double, 2, 2> &a, const double &R, const double &Theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;
}

void A_matrix(Matrix<double, 3, 2> &A, const double &R, const double &Theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;
}

double kinetic_energy( const double &R, const double &Theta, const double &pR, const double &pT, const double &Jx, const double &Jy, const double &Jz)
{
	Vector3d j_vector(Jx, Jy, Jz);
	Vector2d p_vector(pR, pT);

	Matrix<double, 3, 3> I;
	Matrix<double, 2, 2> a;
	Matrix<double, 3, 2> A;

	inertia_tensor( I, R, Theta );
	a_matrix( a, R, Theta );
	A_matrix( A, R, Theta );

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

void W_matrix( Matrix<double, 3, 3> &W, const double &theta, const double &psi )
{
	double sin_psi = sin( psi );
	double cos_psi = cos( psi );

	double sin_theta = sin( theta );
	double cos_theta = cos( theta );

	W(0, 0) = sin_psi / sin_theta;
	W(0, 1) = cos_psi;
	W(0, 2) = - sin_psi * cos_theta / sin_theta;

	W(1, 0) = cos_psi / sin_theta;
	W(1, 1) = - sin_psi;
	W(1, 2) = - cos_psi * cos_theta / sin_theta;

	W(2, 0) = 0;
	W(2, 1) = 0;
	W(2, 2) = 1;
}

double kinetic_energy_euler( const double &R, const double &Theta, const double &pR, const double &pT, const double &phi, const double &theta, const double &psi, const double &p_phi, const double &p_theta, const double &p_psi )
{
	Vector3d pe( phi, theta, psi );
	Vector2d p( pR, pT );

	// ##################################################################
    // filling matrices I, a, A, W
	//
	Matrix<double, 3, 3> W;
	W_matrix( W, theta, psi ); 
	
	Matrix<double, 3, 3> I;
	inertia_tensor( I, R, theta );
	
	Matrix<double, 2, 2> a;
	a_matrix( a, R, theta );
	
	Matrix<double, 3, 2> A;
	A_matrix( A, R, theta );
	// ##################################################################

	// ##################################################################
	// filling matrices G
	//
	Matrix<double, 3, 3> I_inv = I.inverse();
	Matrix<double, 2, 2> a_inv = a.inverse();

	Matrix<double, 3, 3> t1 = I;
	t1.noalias() -= A * a_inv * A.transpose();
	Matrix<double, 3, 3> G11 = t1.inverse();

	Matrix<double, 2, 2> t2 = a;
	t2.noalias() -= A.transpose() * I_inv * A;
	Matrix<double, 2, 2> G22 = t2.inverse();

	Matrix<double, 3, 2> G12;
	G12.noalias() = - G11 * A * a.inverse();
	// ##################################################################

	double ang_term = 0.5 * pe.transpose() * W.transpose() * G11 * W * pe;
	double kin_term = 0.5 * p.transpose() * G22 * p;
	double cor_term = pe.transpose() * W.transpose() * G12 * p;

	return ang_term + kin_term + cor_term;
}



