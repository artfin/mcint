#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// step size
double STEP = 0.5;

// uniform random in [0, 1]
double RAND01()
{
	return (double) rand() / (double) RAND_MAX;
}

// probability profile
double profile( double x, double y )
{
	double r2 = x*x + y*y;
	double p = r2 * exp( -r2 );

	return p;
}

//take a Metropolis step
void metropolis_step( double *x, double *y, int *rejected )
{
	double angle;
	double xt, yt, xx, yy;
	double old_probability, new_probability;
	double r, ratio;

	xx = *x;
	yy = *y;
	old_probability = profile( xx, yy );

	// trial step
	angle = 2 * M_PI * RAND01();
	xt = *x + STEP * cos( angle );
	yt = *y + STEP * sin( angle );

	// new probability
	new_probability = profile( xt, yt );

	// keep trial point if ratio > 1
	if ( new_probability > old_probability )
	{
		*x = xt;
		*y = yt;
		return;
	}

	// keep trial point with ratio probability
	ratio = new_probability / old_probability;
	r = RAND01();
	
	if ( r < ratio )
	{
		*x = xt;
		*y = yt;
		return;
	}

	*x = xx;
	*y = yy;
	*rejected++;
	return;

}

// set (reset) parameters from argument line
void set_parameters( int argc, char *argv[], int *print_walk, int *nsteps, int *warmup_steps, int *exit_parameter )
{
	*warmup_steps = 0;
	*nsteps = 100;
	*print_walk = 0;
	
	int help = 0;

	for ( int n = 1; n < argc; n++ )
	{
		if ( argv[n][0] == '-' )
		{
			switch ( argv[n][1] )
			{
				case 'p':
					*print_walk = 1;
					break;
				case 'w':
					*warmup_steps = atoi( argv[++n] );
					break;
				case 's':
					*nsteps = atoi( argv[++n] );
					break;
				default:
					help = 1;
			}
		}
	}

	if ( help == 1 )
	{
		cout << ">> HELP: ./guided walk -p print_walk [0] -w warmup_steps[0] -s nsteps[0] " << endl << endl;
		cout << "Exiting program..." << endl;

		*exit_parameter = 1;
	}
}

int main( int argc, char *argv[] )
{
	double x, y, xt, yt;
	int nsteps, rejected, warmup_steps, exit_parameter;
	int print_walk;

	// seed the random numbers
	srand( 171751 );

	// controls -- command line arguments
	set_parameters( argc, argv, &print_walk, &nsteps, &warmup_steps, &exit_parameter );
	
	if ( exit_parameter == 1 )
	{
		exit( 1 );
	}
	
	cout << "-------------------------------------" << endl;
	cout << "Parsed command-line arguments: " << endl;
	cout << "print_walk: " << print_walk << endl;
    cout << "nsteps: " << nsteps << endl;
	cout << "warmup_steps: " << warmup_steps << endl;	
	cout << "-------------------------------------" << endl;

	// warmup Metropolis steps
	for ( int istep = 0; istep < warmup_steps; istep++ )
	{
		metropolis_step( &x, &y, &rejected );
	}

	// generate the walk
	rejected = 0;
	for ( int istep = 0; istep < nsteps; istep++ )
	{
		metropolis_step( &x, &y, &rejected );

		if ( print_walk == 1 )
		{
			cout << "step: " << istep << "; x: " << x << "; y: " << y << endl;
		}
	}
	
	cout << "-----------------------" << endl << endl;
	cout << "NSTEPS: " << nsteps << "; STEP: " << STEP << endl;
	cout << "Accepted: " << nsteps - rejected << "; percent: " << 100.0 * (double) ( nsteps - rejected ) / (double) nsteps << endl;

	return 0;
}


