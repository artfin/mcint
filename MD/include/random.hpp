#pragma once

#include <random>
#include <ctime>

// thread_local is not supported in clang 
#if defined(__clang__)
static std::mt19937 generator;
#else
static thread_local std::mt19937 generator;
#endif

class Random
{
    public:
        static void seed( uint64_t seed = std::mt19937_64::default_seed )
        {
            generator.seed( seed );
        }

        static void randomSeed() 
        {
            generator.seed( std::time(0) );
        }

        // random float between min and max
        static float nextFloat( const float &min, const float &max )
        {
            std::uniform_real_distribution<float> distribution( min, max );
            return distribution( generator );
        }

        // random float between 0 and 1
        static float nextFloat() 
        {
            std::uniform_real_distribution<float> distribution( 0, 1 );
            return distribution( generator );
        }

        // random double between min and max
        static double nextDouble( const double &min, const double &max )
        {
            std::uniform_real_distribution<double> distribution( min, max );
            return distribution( generator );
        }
        
        // random double between 0 and 1
        static double nextDouble()
        {
            std::uniform_real_distribution<double> distribution( 0, 1 );
            return distribution( generator );
        }

        // normally distributed random number
        static double nextGaussian( const double &mean, const double &sigma )
        {
            std::normal_distribution<double> distribution( mean, sigma );
            return distribution( generator );
        }

        // random int between min and max
        static int nextInt( const int &min, const int &max )
        {
            std::uniform_int_distribution<int> distribution( min, max );
            return distribution( generator );
        }

        // random long between min and max
        static long nextLong( const long &min, const long &max )
        {
            std::uniform_int_distribution<long> distribution( min, max );
            return distribution( generator );
        }

        // random bool
        static bool nextBool()
        {
            std::uniform_int_distribution<int> distribution( 0, 1 );
            return distribution( generator );
        }
};
