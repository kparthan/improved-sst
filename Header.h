#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <pwd.h>
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
#include <liblcb/liblcb.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

// numeric constants
#define AOM 0.001
#define AOM_angle 0.2
#define LARGE_NUMBER 1000000000
//#define PI boost::math::constants::pi<double>()
#define PI M_PI 
#define LOG_PI log(PI)
#define LOG2_PI log2(PI)
#define ZERO std::numeric_limits<double>::epsilon()
#define TOLERANCE 1e-6
#define NORMAL_MEAN 3.8
#define NORMAL_SIGMA 0.2

#define SET 1 
#define UNSET 0

#define PRINT_NON_DETAIL 0
#define PRINT_DETAIL 1

#define DEFAULT_RESOLUTION 0.1
#define DEFAULT_COMPONENTS 2
#define MAX_COMPONENTS 10
#define MAX_KAPPA ceil(1/(0.09*AOM_angle*AOM_angle))
//#define MAX_KAPPA 50 

// for mixture components visualization
#define VISUALIZE_2D 0
#define VISUALIZE_3D 1
#define USING_MIXTURE_WEIGHTS 0
#define RANDOM_SAMPLE_SIZE 1
#define DEFAULT_SAMPLE_SIZE 2000

// for sst
#define NUM_IDEAL_MODELS 7 
#define MIN_SIZE_HELIX  4 
#define MIN_SIZE_STRAND 4
#define MAX_SEGMENT_SIZE 40 
#define DEFAULT_ORIENTATION 1
#define FIT_ENTIRE_STRUCTURE 0
#define FIT_SINGLE_SEGMENT 1

#define ONE_COMPONENT_ADAPTIVE 1
#define MIXTURE_ADAPTIVE 2
#define NON_ADAPTIVE 3

using namespace std;
using namespace std::chrono;
using namespace lcb;
using namespace lcb::geometry;
using namespace boost::program_options;
using namespace boost::filesystem;

#endif

