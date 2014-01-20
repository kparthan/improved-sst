#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <memory>
#include <cstdlib>
#include <vector>
#include <array>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <ctime>
#include <cassert>
#include <thread>
#include <chrono>
#include <omp.h>
#include <liblcb/liblcb.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

#define HOME_DIRECTORY "/home/pkas7/"
#define CURRENT_DIRECTORY "/home/pkas7/Research/Work/improved-sst/"

// numeric constants
#define AOM 0.001
#define AOM_angle 0.2
#define LARGE_NUMBER 1000000000
#define PI boost::math::constants::pi<double>()
#define LOG_PI log(PI)
#define ZERO std::numeric_limits<double>::epsilon()
#define TOLERANCE 1e-6

#define SET 1 
#define UNSET 0

#define PRINT_NON_DETAIL 0
#define PRINT_DETAIL 1

#define DEFAULT_RESOLUTION 0.1
#define DEFAULT_COMPONENTS 2
#define MAX_COMPONENTS 10
#define MAX_KAPPA ceil(1/(0.09*AOM_angle*AOM_angle))
//#define MAX_KAPPA 100 

// for mixture components visualization
#define VISUALIZE_2D 0
#define VISUALIZE_3D 1
#define USING_MIXTURE_WEIGHTS 0
#define RANDOM_SAMPLE_SIZE 1
#define DEFAULT_SAMPLE_SIZE 2000

// for sst
#define NUM_IDEAL_MODELS 7
#define MIN_SEGMENT_SIZE 4
#define MAX_SEGMENT_SIZE 49

using namespace std;
using namespace std::chrono;
using namespace lcb;
using namespace lcb::geometry;
using namespace boost::program_options;
using namespace boost::filesystem;

#endif

