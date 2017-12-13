#pragma once

#include <iostream>
#include <assert.h>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <functional>
#include <iomanip>
#include <algorithm> // std::min

#include <gsl/gsl_histogram.h>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::function;
using std::vector;
using std::pair;
using std::string;

using Eigen::VectorXf;

const int NBINS = 500;

static vector<pair<int,double>> DEFAULT_VECTOR;
static int BLOCK_SIZE = 1e6;

class MCMC
{
	public:
		int DIM;
		int MOVES;
		double alpha;
		vector<gsl_histogram*> histograms;
		vector<string> names;

		bool set_histograms = false;

		VectorXf point_after_burnin;
		function<double(VectorXf)> f;
		std::mt19937 generator;
	
		double nextDouble( const double &min, const double& max );
		void nextGaussianVec( VectorXf &v, VectorXf &mean );
		double wrapMax( const double& x, const double& max );

		void initialize_histograms( vector<gsl_histogram*> histograms,
					   				vector<string> names
								  );
		void save_histogram( gsl_histogram *histogram, string filename );
		void gsl_histogram_normalize( gsl_histogram* h );

		VectorXf metro_step( VectorXf& x );
		void burnin( VectorXf& initial_point, const int& burnin_length );
		void run_chain( 
			vector<pair<int, double>>& to_wrap = DEFAULT_VECTOR, 
			int& block_size = BLOCK_SIZE 
					  );

		MCMC( function<double(VectorXf)> f, const int &MOVES, const int &DIM, const double& alpha );
		~MCMC();
};
