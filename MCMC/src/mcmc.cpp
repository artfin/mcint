#include "mcmc.hpp"

MCMC::MCMC ( function<double(VectorXf)> f, const int& MOVES, const int& DIM, const double& alpha ) : DIM (DIM), MOVES(MOVES), alpha(alpha), f(f)
{ 
}

MCMC::~MCMC()
{
	if ( set_histograms )
	{
		for ( int i = 0; i < DIM; i++ )
		{
			gsl_histogram_free( histograms[i] );
		}
	}
}

void MCMC::gsl_histogram_normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
	double min = gsl_histogram_min( h );
	double step = (max - min) / NBINS;

	double sum = gsl_histogram_sum( h ) * step;
	gsl_histogram_scale( h, 1.0 / sum );
}

void MCMC::save_histogram( gsl_histogram *histogram, string filename )
{
	std::ofstream file( filename );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
}

void MCMC::initialize_histograms( vector<gsl_histogram*> histograms,
			   					  vector<string> names )
{
	assert( histograms.size() == this->DIM );

	this->set_histograms = true;
	this->histograms = histograms;
	this->names = names;
}

double MCMC::nextDouble( const double& min, const double& max )
{
	std::uniform_real_distribution<double> distribution( min, max );
	return distribution( this->generator );
}

void MCMC::nextGaussianVec( VectorXf &v, VectorXf &mean )
{
	for ( int i = 0; i < this->DIM; i++ )
	{
		std::normal_distribution<double> d( mean(i), this->alpha );
		v(i) = d( this->generator );
	}
}

// wrap x -> [0, max)
double MCMC::wrapMax( const double& x, const double& max )
{
	return fmod( max + fmod(x, max), max );
}

VectorXf MCMC::metro_step( VectorXf& x )
{
	VectorXf prop( this->DIM );
	nextGaussianVec( prop, x );

	if ( nextDouble(0.0, 1.0) < std::min( 1.0, this->f(prop) / this->f(x) ))
	{
		return prop;
	}
	
	return x;	
}

void MCMC::burnin( VectorXf& initial_point, const int& burnin_length )
{
	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

	VectorXf x = metro_step( initial_point );

	for ( size_t i = 0; i < burnin_length; i++ )
	{
		x = metro_step( x );
	}

	cout << "Burnin finished. Time elapsed: " << 
		std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s " << endl;

	this->point_after_burnin = x;
}

void MCMC::run_chain( vector<pair<int, double>>& to_wrap, int& block_size )
{
	int moves = 0;
	int attempted_steps = 0;	

	VectorXf x = this->point_after_burnin;
	VectorXf xnew( this->DIM );

	std::chrono::milliseconds time_for_block;
	std::chrono::milliseconds time_for_blocks;

	vector<std::chrono::high_resolution_clock::time_point> block_times;
	block_times.push_back( std::chrono::high_resolution_clock::now() );

	int block_counter = 0;

	while ( moves < this->MOVES )
	{
		if ( attempted_steps % block_size == 0 && attempted_steps != 0 )
		{
			block_counter++;

			block_times.push_back( std::chrono::high_resolution_clock::now() );

			time_for_block = std::chrono::duration_cast<std::chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = std::chrono::duration_cast<std::chrono::milliseconds>( block_times.end()[-1] - block_times[0] );

			cout << endl << "Block " << block_counter << " finished." << endl;
			cout << "Attempted steps: " << attempted_steps << "; moves made: " << moves << endl;
			cout << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cout << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;
		}
		
		attempted_steps++;

		xnew = metro_step( x );

		if ( to_wrap != DEFAULT_VECTOR )
		{
			int curr_var; // number of current variable
			double curr_max; // maximum of current variable
			double before;

			for ( size_t i = 0; i < to_wrap.size(); i++ )
			{
				curr_var = to_wrap[i].first;
				curr_max = to_wrap[i].second;
				
				before = xnew(curr_var);
				xnew(curr_var) = wrapMax( xnew(curr_var), curr_max );
				//cout << "(BEFORE): " << before << "; (AFTER): " << xnew(curr_var) << endl;
			}
		}

		if ( xnew != x )
		{
			x = xnew;
			moves++;
		}

		if ( set_histograms )
		{
			for ( size_t i = 0; i < DIM; i++ )
			{
				gsl_histogram_increment( this->histograms[i], xnew(i) );
			}
		}
	}
	
	if ( set_histograms )
	{
		for ( size_t i = 0; i < DIM; i++ )
		{
			gsl_histogram_normalize( this->histograms[i] );
			save_histogram( this->histograms[i], this->names[i] );
		}
	}

    cout << "-----------------------------------" << endl;
   	cout << "Attempted steps: " << attempted_steps << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100.0 << "%" << endl;
   	cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - block_times[0]).count() / 1000.0 << " s" << endl; 
   	cout << "-----------------------------------" << endl;
}

