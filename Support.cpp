#include "Support.h"
#include "Geometry3D.h"
#include "VonMises3D.h"
#include "Normal.h"

int initialize_components_from_file;
string HOME_DIRECTORY,CURRENT_DIRECTORY,STRUCTURE;
vector<double> ORIGIN = {0,0,0};
vector<double> XAXIS = {1,0,0};
vector<double> NEGATIVE_XAXIS = {-1,0,0};
vector<double> YAXIS = {0,1,0};
vector<double> ZAXIS = {0,0,1};
int DSSP_DATA_COLLECT;
int DEBUG;
//ofstream debug("debug");

//////////////////////// GENERAL PURPOSE FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This function sets the home and current working directory variables.
 */
void getHomeAndCurrentDirectory()
{
  struct passwd *pw = getpwuid(getuid());
  HOME_DIRECTORY = pw->pw_dir;
  //cout << "home_dir: " << HOME_DIRECTORY << endl;

  char current_dir[256];
  CURRENT_DIRECTORY = getcwd(current_dir,255);
  //cout << "current_dir: " << CURRENT_DIRECTORY << endl;
}

/*!
 *  \brief This function checks to see if valid arguments are given to the 
 *  command line output.
 *  \param argc an integer
 *  \param argv an vector of strings
 *  \return the parameters of the model
 */
struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string pdb_id,scop_id,generation,constrain,sst_method;

  bool noargs = 1;

  cout << "Checking command-line input ..." << endl;
  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("verbose","print some details")
       ("file",value<string>(&parameters.file),"pdb file")
       ("pdbid",value<string>(&pdb_id),"PDB ID")
       ("scopid",value<string>(&scop_id),"SCOP ID")
       ("force","force build the profile")
       ("directory",value<string>(&parameters.profiles_dir),"path to all profiles")
       ("constrain",value<string>(&constrain),"to constrain kappa")
       ("bins","to compute the frequencies of angles")
       ("res",value<double>(&parameters.res),"heat map resolution")
       ("mixture","flag to do mixture modelling")
       ("infer_components","flag to infer the number of components")
       ("k",value<int>(&parameters.fit_num_components),"number of components")
       ("update_weights_new","flag to update weights using modified rule")
       // parameters to aid in visualization/modelling of mixture components
       ("load",value<string>(&parameters.mixture_file),"mixture file")
       ("generate",value<string>(&generation),"method to generate samples")
       ("samples",value<int>(&parameters.num_samples),"sample size")
       ("simulate","flag to run the simulation")
       ("components",value<int>(&parameters.simulate_num_components),
                              "# of mixture components used in the simulation")
       ("initialize_components_from_file","to read the components")
       // parameters to run sst
       ("sst","flag to execute sst script")
       ("orientation",value<int>(&parameters.orientation),
          "orientation of the mean/direction used in the adaptive encoding")
       ("segment",value<vector<string>>(&parameters.end_points)->multitoken(),
                                  "segment to be fit")
       ("debug","flag to print out values to assist in debugging")
       ("method",value<string>(&sst_method),"the sst method to be used")
       // parameters to run DSSP
       ("dssp","flag to use DSSP output")
       ("type",value<string>(&parameters.dssp_sst_type),"model type")
       ("parse",value<string>(&parameters.dssp_file),"file with the DSSP assignment")
       ("collect",value<int>(&parameters.dssp_data_collect),"method to collect DSSP data")
  ;
  variables_map vm;
  store(parse_command_line(argc,argv,desc),vm);
  notify(vm);

  if (vm.count("help")) {
    Usage(argv[0],desc);
  }

  if (vm.count("force")) {
    parameters.force = SET;
  } else {
    parameters.force = UNSET;
  }

  if (vm.count("debug")) {
    DEBUG = SET;
  } else {
    DEBUG = UNSET;
  }

  if (vm.count("initialize_components_from_file")) {
    initialize_components_from_file = SET;
  } else {
    initialize_components_from_file = UNSET;
  }
  
  if (vm.count("directory")) {
    if (vm.count("file") || vm.count("pdbid") || vm.count("scopid")) {
      cout << "conflicting options: directory/file [USE ONLY ONE] ..." << endl;
      Usage(argv[0],desc);
    }
    parameters.read_profiles = SET;
  } else {
    parameters.read_profiles = UNSET;
  }

  if (constrain.compare("kappa") == 0) {
    parameters.constrain_kappa = SET;
  } else {
    parameters.constrain_kappa = UNSET;
  }

  // to JUST build the angular profile
  if (!vm.count("directory") && !vm.count("mixture")) {
    if (vm.count("file") && vm.count("pdbid")) {
      cout << "Please use one of --file or --pdbid ..." << endl;
      Usage(argv[0],desc);
    } else if (vm.count("file")) {
      cout << "Using protein file: " << parameters.file << endl;
    } else if (vm.count("pdbid")) {
      cout << "Using PDB ID: " << pdb_id << endl;
      parameters.file = getPDBFilePath(pdb_id);
    } else if (vm.count("scopid")) {
      cout << "Using SCOP ID: " << scop_id << endl;
      parameters.file = getSCOPFilePath(scop_id);
    } 
  }
  noargs = 0;

  if (vm.count("mixture") && vm.count("sst")) {
    cout << "Not allowed to use \"mixture\" and \"sst\" options together\n";
    Usage(argv[0],desc);
  }

  if (vm.count("mixture")) {  // run mixture modelling
    parameters.mixture_model = SET;
    if (vm.count("load")) {  // for visualization/simulation
      parameters.load_mixture = SET;
      if (generation.compare("using_mixture_weights") == 0) {
        parameters.sample_generation = USING_MIXTURE_WEIGHTS;
        if (!vm.count("samples")) {
          parameters.num_samples = DEFAULT_SAMPLE_SIZE; 
        }
      } else if (generation.compare("random_sample_size") == 0) {
        parameters.sample_generation = RANDOM_SAMPLE_SIZE;
      } else {
        cout << "Unrecognized generation flag ..." << endl;
        Usage(argv[0],desc);
      }
    } else {
      parameters.load_mixture = UNSET;
    }
    if (vm.count("simulate")) {
      parameters.simulation = SET;
      cout << "Simulating a mixture model ..." << endl;
      if (!vm.count("components") && !vm.count("load")) {
        parameters.simulate_num_components = DEFAULT_COMPONENTS;
        cout << "# of components used in simulation: " 
             << parameters.simulate_num_components << endl;
      }
      if (!vm.count("samples")) {
        parameters.num_samples = DEFAULT_SAMPLE_SIZE; 
      }
    } else {
      parameters.simulation = UNSET;
    }
    if (vm.count("simulate") || !vm.count("load")) {
      if (vm.count("infer_components")) {
        parameters.infer_num_components = SET;
      } else {
        parameters.infer_num_components = UNSET;
        if (!vm.count("k")) {
          parameters.fit_num_components = DEFAULT_COMPONENTS; 
        } else if (parameters.fit_num_components == 0) {
          cout << "# of inferred components should be non-zero ..." << endl;
          Usage(argv[0],desc);
        }
        cout << "# of inferred components: " << parameters.fit_num_components 
             << endl;
      }
      if (vm.count("update_weights_new")) {
        parameters.update_weights_new = SET;
      } else {
        parameters.update_weights_new = UNSET;
      }
    }
  } else {  // run for single von mises
    parameters.mixture_model = UNSET;
  }

  if (vm.count("bins")) {
    parameters.heat_map = SET;
    if (!vm.count("res")) {
      parameters.res = DEFAULT_RESOLUTION;
    }
  } else {
    parameters.heat_map = UNSET;
  }

  if (vm.count("sst")) {
    parameters.sst = SET;
    if (vm.count("load")) { 
      parameters.load_mixture = SET;
    } else if (!vm.count("dssp")) {
      cout << "Please supply a mixture file ...\n";
      Usage(argv[0],desc);
    }
    if (vm.count("method")) {
      if (sst_method.compare("mixture_adaptive") == 0) {
        parameters.method = MIXTURE_ADAPTIVE;
      } else if (sst_method.compare("one_component_adaptive") == 0) {
        parameters.method = ONE_COMPONENT_ADAPTIVE;
      } else if (sst_method.compare("non_adaptive") == 0) {
        parameters.method = NON_ADAPTIVE;
      } else if (sst_method.compare("dssp_non_adaptive") == 0) {
        parameters.method = DSSP_NON_ADAPTIVE;
      }
    } else {
      parameters.method = MIXTURE_ADAPTIVE;
    }
    if (vm.count("pdbid")) {
      parameters.file = getPDBFilePath(pdb_id);
    } else if (vm.count("scopid")) {
      parameters.file = getSCOPFilePath(scop_id);
    }
    if (!vm.count("orientation")) {
      parameters.orientation = DEFAULT_ORIENTATION;
    }
    if (vm.count("segment")) {
      cout << "Fitting a single segment between the residues "
           << "[" << parameters.end_points[1] << ", " << parameters.end_points[2]
           << "] of chain " << parameters.end_points[0] << endl;
      parameters.portion_to_fit = FIT_SINGLE_SEGMENT;
    } else {
      parameters.portion_to_fit = FIT_ENTIRE_STRUCTURE;
    }
    parameters.simulation = UNSET;
  } else {
    parameters.sst = UNSET;
  }

  if (vm.count("dssp")) {
    parameters.dssp = SET;
    if (vm.count("parse")) {
      if(checkFile(parameters.dssp_file) == 0) {
        cout << "DSSP file: " << parameters.dssp_file << " doesn't exist.\n";
        Usage(argv[0],desc);
      }
      parameters.parse_dssp = SET;
      if (vm.count("pdbid")) {
        parameters.file = getPDBFilePath(pdb_id);
      } else if (vm.count("scopid")) {
        parameters.file = getSCOPFilePath(scop_id);
      }
    } else {
      if (!vm.count("type")) {
        parameters.dssp_sst_type = "all";
      }
    }
    DSSP_DATA_COLLECT = parameters.dssp_data_collect;
  } else {
    parameters.dssp = UNSET;
  }

  if (noargs) {
    cout << "Not enough arguments supplied..." << endl;
    Usage(argv[0],desc);
  }

  return parameters;
}

/*!
 *  \brief This module prints the acceptable input format to the program
 *  \param exe a reference to a const char
 *  \param desc a reference to a options_description object
 */
void Usage(const char *exe, options_description &desc)
{
  cout << "Usage: " << exe << " [options]" << endl;
  cout << desc << endl;
  exit(1);
}

/*!
 *  \brief This module checks whether the input file exists or not.
 *  \param file_name a reference to a string
 *  \return true or false depending on whether the file exists or not.
 */
bool checkFile(string &file_name)
{
  ifstream file(file_name);
  return file;
}

/*!
 *  \brief This module prints the elements of a vector<vector<>> to a file
 *  \param v a reference to vector<vector<double>>
 *  \param file_name a pointer to a const char
 */
void writeToFile(vector<vector<double>> &v, const char *file_name)
{
  ofstream file(file_name);
  for (int i=0; i<v.size(); i++) {
    file << "(";
    for (int j=0; j<v[i].size()-1; j++) {
      file << v[i][j] << ", ";
    }
    file << v[i][v[i].size()-1] << ")" << endl;
  }
  file.close(); 
}

/*!
 *  \brief This module extracts the file name from the path
 *  \param file a reference to a string
 *  \return the extracted portion of the file name
 */
string extractName(string &file)
{
  unsigned pos1 = file.find_last_of("/");
  unsigned pos2 = file.find(".");
  int length = pos2 - pos1 - 1;
  string sub = file.substr(pos1+1,length);
  return sub;
}

/*!
 *  \brief This function initializes a matrix.
 *  \param matrix a reference to a vector<vector<double>>
 *  \param rows an integers
 *  \param cols an integer
 */
void initializeMatrix(vector<vector<double>> &matrix, int rows, int cols)
{
  vector<double> tmp(cols,0);
  for (int i=0; i<rows; i++) {
    matrix.push_back(tmp);
  }
}

/*!
 *  \brief This function converts the cartesian coordinates into spherical.
 *  \param cartesian a reference to a vector<double> 
 *  \param spherical a reference to a vector<double> 
 */
void cartesian2spherical(vector<double> &cartesian, vector<double> &spherical)
{
  double r = vectorNorm(cartesian);

  double x = cartesian[0] / r;
  double y = cartesian[1] / r;
  double z = cartesian[2] / r;

  // theta \in [0,PI]: angle with Z-axis
  double theta = acos(z);

  // phi \in[0,2 PI]: angle with positive X-axis
  double ratio = x/sin(theta);
  if (ratio > 1) {
    ratio = 1;
  } else if (ratio < -1) {
    ratio = -1;
  }
  double angle = acos(ratio);
  double phi = 0;
  if (x == 0 && y == 0) {
    phi = 0;
  } else if (x == 0) {
    if (y > 0) {
      phi = angle;
    } else {
      phi = 2 * PI - angle;
    }
  } else if (y >= 0) {
    phi = angle;
  } else if (y < 0) {
    phi = 2 * PI - angle;
  }

  spherical[0] = r;
  spherical[1] = theta;
  spherical[2] = phi;
}

/*!
 *  \brief This function computes the unit spherical coordinates of a vector.
 *  Any given vector in Cartesian coordinates is first transformed to its 
 *  spherical equivalent, which is then converted to a vector on a unit sphere
 *  with direction (theta,phi) being the same as the original vector. 
 *  \param cartesian a reference to a vector<double> 
 *  \param unit_spherical a reference to a vector<double> 
 */
void cartesian2unitspherical(vector<double> &cartesian, 
                             vector<double>&unit_spherical)
{
  vector<double> spherical(3,0);
  cartesian2spherical(cartesian,spherical);
  spherical[0] = 1;
  spherical2cartesian(spherical,unit_spherical);
}

/*!
 *  \brief This function converts the spherical coordinates into cartesian.
 *  \param spherical a reference to a vector<double> 
 *  \param cartesian a reference to a vector<double> 
 */
void spherical2cartesian(vector<double> &spherical, vector<double> &cartesian)
{
  cartesian[0] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
  cartesian[1] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
  cartesian[2] = spherical[0] * cos(spherical[1]);
}

/*!
 *  \brief This function converts a liblcb::Point object to a std::vector object
 *  \param p a reference to a Point<double>
 *  \param v a reference to a vector<double>
 */
void point2vector(Point<double> &p, vector<double> &v)
{
  v[0] = p.x();
  v[1] = p.y();
  v[2] = p.z();
}

/*!
 *  \brief This function scales the value to AOM.
 *  \param angle_rad (in radians) a reference to a double
 */
void scaleToAOM(double &angle_rad)
{
  int scale = 1 / (double) AOM_angle;
  int angle = (angle_rad) * scale;
  angle_rad = angle / (double) scale;
  cout << angle_rad << endl;
}

/*!
 *  \brief This function finds the minimum of the two elements 
 *  \param a an element of type Realtype 
 *  \param b an element of type Realtype 
 *  \return the minimum value 
 */
template <typename RealType>
RealType minimum(RealType a, RealType b)
{
  if (a <= b) {
    return a;
  } else {
    return b;
  }
}
template int minimum(int,int);
template float minimum(float,float);
template double minimum(double,double);
template long double minimum(long double,long double);

/*!
 *  \brief This function gets the index of the maximum element.
 *  \param list a reference to a vector<RealType>
 *  \return the index
 */
int getIndexOfMaximumElement(vector<double> &list)
{
  int max_index = 0;
  double max_val = list[max_index];
  for (int i=1; i<list.size(); i++) {
    if (list[i] > max_val) {
      max_val = list[i];
      max_index = i;
    }
  }
  return max_index;
}

/*!
 *  \brief This function prints the elements of an vector.
 *  \param os a reference to a ostream
 *  \param v a reference to a vector<double>
 */
void print(ostream &os, vector<double> &v)
{
  os << fixed << setprecision(4) << "(" << v[0] << "," << v[1] << "," << v[2] <<
  //os << fixed << setprecision(4) << "(" << v[0] << "," << v[1]*180/PI << "," << v[2]*180/PI <<
")\t";
}
/*!
 *  \brief This function prints the elements of an array.
 *  \param os a reference to a ostream
 *  \param a a reference to a array<double,3>
 */
void print(ostream &os, array<double,3> &a)
{
  os << fixed << setprecision(4) << "(" << a[0] << "," << a[1] << "," << a[2] <<
")\t";
}

/*!
 *  \brief This function computes the ratio of Bessel functions -- A3(k)
 *  \param k a double
 *  \return the value of the ratio at a given k
 */
double ratioBesselFunction(double k)
{
  if (k < TOLERANCE) {
    return 0;
  } else {
    double k_inv = 1 / (double)k;
    double cothk = 1 / (double) tanh(k);
    return cothk - k_inv;
  }
}

/*!
 *  \brief This function computes the first derivative of the 
 *  ratio of Bessel functions -- A3'(k)
 *  \param k a double
 *  \return the value of the first derivative of the ratio at a given k
 */
double ratioBesselFunction_firstDerivative(double k)
{
  if (k < TOLERANCE) {
    return 0;
  } else {
    double value = 0;
    double cschk = 1 / (double) sinh(k);
    value -= cschk * cschk;
    double k_inv = 1 / (double)k;
    value += k_inv * k_inv; 
    return value ;
  }
}

/*!
 *  \brief This function computes the second derivative of the 
 *  ratio of Bessel functions -- A3''(k)
 *  \param k a double
 *  \return the value of the second derivative of the ratio at a given k
 */
double ratioBesselFunction_secondDerivative(double k)
{
  if (k < TOLERANCE) {
    return 0;
  } else {
    double value = 0;
    double cothk = 1 / (double) tanh(k);
    double cschk = 1 / (double) sinh(k);
    value += cothk * cschk * cschk;
    double k_inv = 1 / (double)k;
    value -= k_inv * k_inv * k_inv;
    return 2 * value ;
  }
}

/*!
 *  \brief This function computes the third derivative of the 
 *  ratio of Bessel functions -- A3'''(k)
 *  \param k a double
 *  \return the value of the third derivative of the ratio at a given k
 */
double ratioBesselFunction_thirdDerivative(double k)
{
  double value = 0;
  double cothk = 1 / (double) tanh(k);
  double cothksq = cothk * cothk;
  double cschk = 1 / (double) sinh(k);
  double cschksq = cschk * cschk;

  value += -2 * cschksq * (cothksq + cschksq);
  double k_inv = 1 / (double)k;
  value += 6 * k_inv * k_inv * k_inv * k_inv;
  return value ;
}

/*!
 *  \brief This function computes the approximation of the constant term for
 *  the constant term in the message length expression (pg. 257 Wallace)
 *  \param d an integer
 *  \return the constant term
 */
double computeConstantTerm(int d)
{
  double ad = 0;
  ad -= 0.5 * d * log(2 * PI);
  ad += 0.5 * log(d * PI);
  return ad;
}

/*!
 *  \brief This function computes the lattice constant for higher dimensions.
 *  \param d an integer
 *  \return the d-dimensional lattice constant
 */
double getLatticeConstant(int d)
{
  double ad = computeConstantTerm(d);
  double tmp = ((2.0/d) * ad) - 1;
  return exp(tmp);
}

/*!
 *  \brief This function converts the angle from degrees to radians.
 *  \param theta (measured in degrees) a reference to a double
 */
void angleInRadians(double &theta)
{
  theta *= PI / 180;
}

/*!
 *  \brief This function sorts the elements in the list
 *  \param list a reference to a vector<double>
 *  \return the sorted list
 */
template <typename RealType>
vector<RealType> sort(vector<RealType> &list)
{
  int num_samples = list.size();
	vector<RealType> sortedList(list);
  vector<int> index(num_samples,0);
	for(int i=0; i<num_samples; i++) {
			index[i] = i;
  }
	quicksort(sortedList,index,0,num_samples-1);
  return sortedList;
}
template vector<int> sort(vector<int> &);
template vector<float> sort(vector<float> &);
template vector<double> sort(vector<double> &);
template vector<long double> sort(vector<long double> &);

/*!
 *  \brief This function sorts the elements in the list
 *  \param list a reference to a vector<double>
 *  \return the sorted list
 */
template <typename RealType>
vector<int> sortedListIndex(vector<RealType> &list)
{
  int num_samples = list.size();
	vector<RealType> sortedList(list);
  vector<int> index(num_samples,0);
	for(int i=0; i<num_samples; i++) {
			index[i] = i;
  }
	quicksort(sortedList,index,0,num_samples-1);
  return index;
}
template vector<int> sortedListIndex(vector<int> &);
template vector<int> sortedListIndex(vector<float> &);
template vector<int> sortedListIndex(vector<double> &);
template vector<int> sortedListIndex(vector<long double> &);

/*!
 *  This is an implementation of the classic quicksort() algorithm to sort a
 *  list of data values. The module uses the overloading operator(<) to 
 *  compare two Point<T> objects. 
 *  Pivot is chosen as the right most element in the list(default)
 *  This function is called recursively.
 *  \param list a reference to a vector<double>
 *	\param index a reference to a vector<int>
 *  \param left an integer
 *  \param right an integer
 */
template <typename RealType>
void quicksort(vector<RealType> &list, vector<int> &index, int left, int right)
{
	if(left < right)
	{
		int pivotNewIndex = partition(list,index,left,right);
		quicksort(list,index,left,pivotNewIndex-1);
		quicksort(list,index,pivotNewIndex+1,right);
	}
}
template void quicksort(vector<float> &, vector<int> &, int, int);
template void quicksort(vector<double> &, vector<int> &, int, int);
template void quicksort(vector<long double> &, vector<int> &, int, int);

/*!
 *  This function is called from the quicksort() routine to compute the new
 *  pivot index.
 *  \param list a reference to a vector<double>
 *	\param index a reference to a vector<int>
 *  \param left an integer
 *  \param right an integer
 *  \return the new pivot index
 */
template <typename RealType>
int partition(vector<RealType> &list, vector<int> &index, int left, int right)
{
	RealType temp,pivotPoint = list[right];
	int storeIndex = left,temp_i;
	for(int i=left; i<right; i++) {
		if(list[i] < pivotPoint) {
			temp = list[i];
			list[i] = list[storeIndex];
			list[storeIndex] = temp;
			temp_i = index[i];
			index[i] = index[storeIndex];
			index[storeIndex] = temp_i;
			storeIndex += 1;	
		}
	}
	temp = list[storeIndex];
	list[storeIndex] = list[right];
	list[right] = temp;
	temp_i = index[storeIndex];
	index[storeIndex] = index[right];
	index[right] = temp_i;
	return storeIndex;
}
template int partition(vector<float> &, vector<int> &, int, int);
template int partition(vector<double> &, vector<int> &, int, int);
template int partition(vector<long double> &, vector<int> &, int, int);


//////////////////////// PROTEIN FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This module returns the file path associated with a PDB ID.
 *  \param pdb_id a reference to a string
 *  \return the file path
 */
string getPDBFilePath(string &pdb_id)
{
  boost::algorithm::to_lower(pdb_id);
  string path = string(HOME_DIRECTORY) + "/Research/PDB/" ;
  string directory(pdb_id,1,2);
  path += directory + "/pdb" + pdb_id + ".ent.gz";
  return path;
}

/*!
 *  \brief This module returns the file path associated with a PDB ID.
 *  \param scop_id a reference to a string
 *  \return the file path
 */
string getSCOPFilePath(string &scop_id)
{
  string path = string(HOME_DIRECTORY) + "/Research/SCOP/pdbstyle-1.75B/" ;
  string directory(scop_id,2,2);
  path += directory + "/" + scop_id + ".ent";
  return path;
}

/*!
 *  \brief This function is used to parse the DSSP assignment file and save
 *  the directional data for different types of secondary structure motifs.
 *  \param parameters a reference to a struct Parameters
 */
void parseDSSP(struct Parameters &parameters)
{
  ifstream file(parameters.dssp_file.c_str());
  int count = 0;
  string line;
  vector<string> chain_lines;
  vector<vector<string>> all_lines;

  int HASH = 3 - 1;
  char prev_chain_id,current_chain_id;
  bool begin = 0;
  bool ignore_current_chain = 1;
  bool new_chain = 1;
  int CHAIN_ID = 12 - 1;
  int CHAIN_BREAK = 14 - 1;
  int CHAIN_END = 15 - 1;

  vector<char> suitable_chains;
  int num_chains = 0;
  int num_suitable_chains = 0;
  ofstream log("dssp_parser.log",ios::app);
  log << STRUCTURE;
  while (getline(file,line)) {
    count++;
    if (line[HASH] == '#') {
      begin = 1;
      log << "\tParsing begins at line # " << count+1 << endl;
    }
    if (begin == 1 && line[HASH] != '#') {
      if (line[CHAIN_BREAK] == '!' && line[CHAIN_END] == '*') { // chain ends
        if (!ignore_current_chain) {
          all_lines.push_back(chain_lines);
          num_suitable_chains++;
          suitable_chains.push_back(current_chain_id);
        }
        new_chain = 1;
        chain_lines.clear();
      } else if (line[CHAIN_BREAK] == '!') {  // ignore the current chain
        ignore_current_chain = 1;
        log << "\tChain break in chain " << current_chain_id << " at line " 
            << count << endl;
      } else {
        current_chain_id = line[CHAIN_ID];
        if (!new_chain) {
          assert(current_chain_id == prev_chain_id);
        } else {
          new_chain = 0;
          ignore_current_chain = 0;
          num_chains++;
        }
        if (!ignore_current_chain) {
          chain_lines.push_back(line);
        }
        prev_chain_id = current_chain_id;
      }
    }
  }
  if (!ignore_current_chain) {
    all_lines.push_back(chain_lines);
    num_suitable_chains++;
    suitable_chains.push_back(current_chain_id);
  }
  log << "\tnum_chains: " << num_chains << ";"
      << "\tnum_suitable: " << num_suitable_chains << " | ";
  for (int i=0; i<all_lines.size(); i++) {
    log << "\t [" << suitable_chains[i] << "," << all_lines[i].size() << "]";
  }
  log << endl;
    
  file.close();

  /*ofstream check("check");
  for (int i=0; i<all_lines.size(); i++) {
    for (int j=0; j<all_lines[i].size(); j++) {
      check << all_lines[i][j] << endl;
    }
  }
  check.close();*/

  if (checkFile(parameters.file)) {
    ProteinStructure *p = parsePDBFile(parameters.file);
    Protein protein(p,STRUCTURE);
    bool success = checkParsedDSSPFile(protein,p,all_lines,log);
    if (success) {
      if (parameters.dssp_data_collect == 1) {
        collectData(protein,all_lines,log);
      } else if (parameters.dssp_data_collect == 2) {
        collectData2(protein,all_lines,log);
      }
    } 
  }
  log << endl;
  log.close();
}

/*!
 *  \brief This function performs some consistency checks on the parsed DSSP 
 *  file by testing it against the PDB file.
 *  \param protein a reference to a Protein object
 *  \param p a pointer to a ProteinStructure object
 *  \param all_lines a reference to a vector<vector<string>>
 *  \param log a reference to a ostream
 *  \return whether the sequences in the parsed matches with the original one
 */
bool checkParsedDSSPFile(Protein &protein, ProteinStructure *p,
                         vector<vector<string>> &all_lines, ostream &log)
{
  int CHAIN_ID = 12 - 1;
  int RESIDUE_ID = 14 - 1;

  string chain_id;
  for (int i=0; i<all_lines.size(); i++) {
    chain_id = all_lines[i][0][CHAIN_ID];
    cout << chain_id << endl;
    // get chain sequence from the parsed file
    string residue_ids;
    for (int j=0; j<all_lines[i].size(); j++) {
      assert(chain_id[0] == all_lines[i][j][CHAIN_ID]);
      residue_ids += all_lines[i][j][RESIDUE_ID];
    }
    //cout << "parsed: " << residue_ids << endl; 
    // get original sequence
    Chain chain = p->getDefaultModel()[chain_id];
    string seq = chain.sequence();
    //cout << "actual: " << seq << endl;
    if (residue_ids.compare(seq) != 0) {
      log << "\tSequences of chain " << chain_id << " don't match.\n";
      return 0;
    }
  }

  vector<string> chains = protein.getChainIds();
  if (chains.size() != all_lines.size()) {
    log << "\t# of suitable chains don't match.\n";
    return 0;
  }
  for (int i=0; i<chains.size(); i++) {
    chain_id = all_lines[i][0][CHAIN_ID];
    if(chain_id != chains[i]) {
      log << "\tCorresponding chains don't match.\n";
      return 0;
    }
  }
  log << "\tThis is a suitable structure." << endl;
  return 1;
}

/*!
 *  \brief This function is used to collect the data from the DSSP assignment file.
 *  \param protein a reference to a Protein object
 *  \param all_lines a reference to a vector<vector<string>>
 *  \param log a reference to a ostream
 */
void collectData(Protein &protein, vector<vector<string>> &all_lines,
                 ostream &log)
{
  int CHAIN_ID = 12 - 1;
  int ASSIGN_ID = 17 - 1;
  int k;
  string line;
  char type;
  string DSSP_DIR = CURRENT_DIRECTORY + "/dssp/parsed1/models/";
  string type_dir;
  type_dir = DSSP_DIR + "coil/" + STRUCTURE + ".profile";
  ofstream coil(type_dir.c_str());
  type_dir = DSSP_DIR + "strand/" + STRUCTURE + ".profile";
  ofstream strand(type_dir.c_str());
  type_dir = DSSP_DIR + "helix_310/" + STRUCTURE + ".profile";
  ofstream helix_310(type_dir.c_str());
  type_dir = DSSP_DIR + "helix_alpha/" + STRUCTURE + ".profile";
  ofstream helix_alpha(type_dir.c_str());
  type_dir = DSSP_DIR + "helix_pi/" + STRUCTURE + ".profile";
  ofstream helix_pi(type_dir.c_str());

  protein.computeSphericalTransformation();
  vector<vector<vector<double>>>
  spherical_coordinates = protein.getSphericalCoordinatesList();
  for (int i=0; i<all_lines.size(); i++) {
    for (int j=3; j<all_lines[i].size(); j++) {
      line = all_lines[i][j];
      type = line[ASSIGN_ID];
      vector<double> spherical = spherical_coordinates[i][j-3];

      //cout << type << endl;
      //if (type == ' ') cout << "yes\n";
      switch(type) {
        case 'E': // extended beta sheet
          strand << all_lines[i][j][CHAIN_ID];
          for (k=0; k<3; k++) {
            strand << fixed << setw(10) << setprecision(4) << spherical[k];
          }
          strand << endl;
          break;

        case 'G': //3-10 helix
          helix_310 << all_lines[i][j][CHAIN_ID];
          for (k=0; k<3; k++) {
            helix_310 << fixed << setw(10) << setprecision(4) << spherical[k];
          }
          helix_310 << endl;
          break;

        case 'H': // alpha helix
          helix_alpha << all_lines[i][j][CHAIN_ID];
          for (k=0; k<3; k++) {
            helix_alpha << fixed << setw(10) << setprecision(4) << spherical[k];
          }
          helix_alpha << endl;
          break;

        case 'I': // pi helix
          helix_pi << all_lines[i][j][CHAIN_ID];
          for (k=0; k<3; k++) {
            helix_pi << fixed << setw(10) << setprecision(4) << spherical[k];
          }
          helix_pi << endl;
          break;

        default:  // coil
          coil << all_lines[i][j][CHAIN_ID];
          for (k=0; k<3; k++) {
            coil << fixed << setw(10) << setprecision(4) << spherical[k];
          }
          coil << endl;
          break;
      }
    }
  }
  coil.close();
  strand.close();
  helix_310.close();
  helix_alpha.close();
  helix_pi.close();
}

/*!
 *  \brief This function is used to collect the data from the DSSP assignment file
 *  using the second kind of parsing scheme.
 *  \param protein a reference to a Protein object
 *  \param all_lines a reference to a vector<vector<string>>
 *  \param log a reference to a ostream
 */
void collectData2(Protein &protein, vector<vector<string>> &all_lines, ostream &log)
{
  int CHAIN_ID = 12 - 1;
  int ASSIGN_ID = 17 - 1;
  int k;
  string line;
  char type[4];
  string DSSP_DIR = CURRENT_DIRECTORY + "/dssp/parsed2/models/";
  string type_dir;
  type_dir = DSSP_DIR + "coil/profiles/" + STRUCTURE + ".profile";
  ofstream coil(type_dir.c_str());
  type_dir = DSSP_DIR + "strand/profiles/" + STRUCTURE + ".profile";
  ofstream strand(type_dir.c_str());
  type_dir = DSSP_DIR + "helix_310/profiles/" + STRUCTURE + ".profile";
  ofstream helix_310(type_dir.c_str());
  type_dir = DSSP_DIR + "helix_alpha/profiles/" + STRUCTURE + ".profile";
  ofstream helix_alpha(type_dir.c_str());
  type_dir = DSSP_DIR + "helix_pi/profiles/" + STRUCTURE + ".profile";
  ofstream helix_pi(type_dir.c_str());

  protein.load(STRUCTURE);
  vector<vector<vector<double>>>
  spherical_coordinates = protein.getSphericalCoordinatesList();
  for (int i=0; i<all_lines.size(); i++) {
    vector<char> sequence = getSequence(all_lines[i]);
    for (int j=3; j<sequence.size(); j++) {
    //for (int j=3; j<all_lines[i].size(); j++) {
      for (int k=0; k<4; k++) {
        type[k] = sequence[j+k-3];
        //type[k] = all_lines[i][j+k-3][ASSIGN_ID];
      }
      if (type[0] == type[1] && type[1] == type[2] && type[2] == type[3]) {
        vector<double> spherical = spherical_coordinates[i][j-3];
        switch(type[0]) {
          case 'E': // extended beta sheet
            strand << all_lines[i][j][CHAIN_ID];
            for (k=0; k<3; k++) {
              strand << fixed << setw(10) << setprecision(4) << spherical[k];
            }
            strand << endl;
            break;

          case 'G': //3-10 helix
            helix_310 << all_lines[i][j][CHAIN_ID];
            for (k=0; k<3; k++) {
              helix_310 << fixed << setw(10) << setprecision(4) << spherical[k];
            }
            helix_310 << endl;
            break;

          case 'H': // alpha helix
            helix_alpha << all_lines[i][j][CHAIN_ID];
            for (k=0; k<3; k++) {
              helix_alpha << fixed << setw(10) << setprecision(4) << spherical[k];
            }
            helix_alpha << endl;
            break;

          case 'I': // pi helix
            helix_pi << all_lines[i][j][CHAIN_ID];
            for (k=0; k<3; k++) {
              helix_pi << fixed << setw(10) << setprecision(4) << spherical[k];
            }
            helix_pi << endl;
            break;

          case 'C':  // coil
            coil << all_lines[i][j][CHAIN_ID];
            for (k=0; k<3; k++) {
              coil << fixed << setw(10) << setprecision(4) << spherical[k];
            }
            coil << endl;
            break;
        }
      }
    }
  }
  coil.close();
  strand.close();
  helix_310.close();
  helix_alpha.close();
  helix_pi.close();
}

/*!
 *  \brief This function converts the sequence of DSSP output assignment string to a
 *  sequence which is made up of G,H,I,E,C.
 *  \param lines a reference to a vector<string>
 *  \return the sequence of symbols
 */
vector<char> getSequence(vector<string> &lines)
{
  int ASSIGN_ID = 17 - 1;
  vector<char> sequence;
  char symbol;
  for (int i=0; i<lines.size(); i++) {
    symbol = lines[i][ASSIGN_ID];
    if (!(symbol == 'E' || symbol == 'G' || symbol == 'H' || symbol == 'I')) {
      symbol = 'C';
    }
    sequence.push_back(symbol);
  }
  assert(lines.size() == sequence.size());
  for (int i=0; i<sequence.size(); i++) {
    cout << sequence[i];
  } cout << endl;
  return sequence;
}

/*!
 *  \brief This method constructs the angular profiles of the protein structure.
 *  \param parameters a reference to a struct Parameters
 */
void buildAngularProfile(struct Parameters &parameters)
{
  bool exist = checkIfSphericalProfileExists(STRUCTURE);
  if (!exist || parameters.force == SET) {
    /* Obtain protein coordinates */
    ProteinStructure *p = parsePDBFile(parameters.file);
    Protein protein(p,STRUCTURE);
    protein.computeSphericalTransformation();
    protein.save();
    updateLogFile(STRUCTURE,protein.getCPUTime(),protein.getNumberOfChains());
  } else {  
    Protein protein;
    cout << "Profile of " << STRUCTURE << " exists ..." << endl;
    protein.load(STRUCTURE);
  }
}

/*!
 *  \brief This module checks for the existence of the spherical coordinate
 *  profile.
 *  \param name a reference to a string
 *  \return present or not 
 */
bool checkIfSphericalProfileExists(string &name)
{
  string spherical_profile = string(CURRENT_DIRECTORY) 
                             + "/spherical_system/profiles/" + name + ".profile";
  return checkFile(spherical_profile); 
}

/*!
 *  \brief This module parses the input PDB file.
 *  \param pdb_file a reference to a string 
 *  \return a pointer to a ProteinStructure object
 */
ProteinStructure *parsePDBFile(string &pdb_file)
{
  if(!checkFile(pdb_file)){
    cout << "\nFile \"" << pdb_file << "\" does not exist ..." << endl;
    ofstream log("files_not_present",ios::app);
    log << extractName(pdb_file) << endl;
    log.close();
    exit(1);
  } else {
    cout << "Parsing PDB file ...";
    BrookhavenPDBParser parser;
    ProteinStructure
    *structure = parser.getStructure(pdb_file.c_str())->select(CASelector());
    ProteinStructure 
    *one_model = new ProteinStructure(structure->getIdentifier());
    one_model->select(CASelector());
    std::shared_ptr<lcb::Model>
    newmodel = std::make_shared<lcb::Model>(structure->getDefaultModel());
    one_model->addModel(newmodel);
    delete structure;
    cout << " [OK]" << endl;
    return one_model;
  }
}

/*!
 *  \brief This function transforms a set of four points (p0,p1,p2,p3) to a 
 *  canonical form such that p2 is the origin, p1 lies on the -ve X-axis and 
 *  p0 lies on the XY plane.
 *  \param four_mer a reference to a vector<vector<double>>
 *  \param transformed_four_mer a reference to a vector<vector<double>>
 *  \param rotation_matrix a reference to a vector<vector<double>>
 */
void convertToCanonicalForm(vector<vector<double>> &four_mer,
                            vector<vector<double>> &transformed_four_mer, 
                            vector<vector<double>> &rotation_matrix)
{
  assert(four_mer.size() == 4);
  vector<double> translate(four_mer[2]);
  vector<vector<double>> interim_points;
  
  // translate p2 to origin
  for (int i=0; i<4; i++) {
    computeDirectionRatios(four_mer[i],translate,transformed_four_mer[i]);
  }

  vector<vector<double>> rotation1,rotation2;
  initializeMatrix(rotation1,3,3);
  initializeMatrix(rotation2,3,3);
  double angle;
  vector<double> normal(3,0);  // unit vector

  // rotate so that p1 is on -ve x-axis
  angle = computeAngle(transformed_four_mer[1],NEGATIVE_XAXIS);
  computeNormal(transformed_four_mer[1],NEGATIVE_XAXIS,normal);
  computeRotationMatrix(normal,angle,rotation1);
  interim_points = transformed_four_mer;
  for (int i=0; i<4; i++) {
    rotateVector(rotation1,interim_points[i],transformed_four_mer[i]);
  }

  // rotate so that p0 is on XY plane
  computeNormal(transformed_four_mer[0],NEGATIVE_XAXIS,normal);
  angle = computeAngle(normal,ZAXIS);
  if (transformed_four_mer[0][2] < 0) {
    angle *= -1;
  }
  computeRotationMatrix(NEGATIVE_XAXIS,angle,rotation2);
  interim_points = transformed_four_mer;
  for (int i=0; i<4; i++) {
    rotateVector(rotation2,interim_points[i],transformed_four_mer[i]);
  }
  multiply(rotation2,rotation1,rotation_matrix);
}

/*!
 *  \brief This function computes the transformation to align a given vector
 *  with the Z-axis.
 *  \param sp a reference to a vector<double>
 *  \param ep a reference to a vector<double>
 *  \param rotation_matrix a reference to a vector<vector<double>>
 */
void alignWithZAxis(vector<double> &sp, vector<double> &ep, 
                    vector<vector<double>> &rotation_matrix)
{
  // find drs corresponding to sp-->ep vector
  vector<double> dratios(3,0);
  computeDirectionRatios(ep,sp,dratios);
  // angle with Z-axis
  double theta = computeAngle(dratios,ZAXIS);
  // find normal to the plane (axis of rotation)
  vector<double> normal(3,0);
  computeCrossProduct(dratios,ZAXIS,normal);
  vector<double> unit_normal(3,0);
  computeDirectionCosines(normal,unit_normal);
  computeRotationMatrix(unit_normal,theta,rotation_matrix);
}

/*!
 *  \brief This function applies the transformation of the ideal model and returns 
 *  the corresponding unit vector.
 *  \param matrix a reference to a Matrix<double>
 *  \param sp a reference to a vector<double>
 *  \param ep a reference to a vector<double>
 *  \param suffstats a reference to a vector<double>
 */
void applyIdealModelTransformation(vector<vector<double>> &rotation_matrix, 
                                   vector<double> &sp, vector<double> &ep,
                                   vector<double> &suffstats)
{
  // find drs corresponding to sp-->ep vector
  vector<double> dratios(3,0);
  computeDirectionRatios(ep,sp,dratios);
  // find dcs of the vector
  vector<double> dcosines(3,0);
  computeDirectionCosines(dratios,dcosines);
  // rotate this unit vector
  vector<double> rotated(3,0);
  rotateVector(rotation_matrix,dcosines,rotated);
  for (int i=0; i<3; i++) {
    suffstats[i] += rotated[i];
  }
}

/*!
 *  \brief This function is used to read the angular profiles and use this data
 *  to estimate parameters of a Von Mises distribution.
 *  \param parameters a reference to a struct Parameters
 */
void computeEstimators(struct Parameters &parameters)
{
  if (parameters.mixture_model == UNSET) {  // no mixture modelling
    vector<double> resultant_direction(3,0);
    double n = 0;
    bool success = readProfiles(parameters,resultant_direction,n);
    if (success) {
      modelOneComponent(parameters,resultant_direction,n);
    }
  } else if (parameters.mixture_model == SET) { // mixture modelling
    vector<vector<double>> unit_coordinates;
    bool success = gatherData(parameters,unit_coordinates);
    if (success) {
      modelMixture(parameters,unit_coordinates);
    }
  }
}

/*!
 *  \brief This function models a single component.
 *  \param parameters a reference to a struct Parameters
 *  \param direction a reference to a vector<double>
 *  \param num_samples a reference to a double
 */
void modelOneComponent(struct Parameters &parameters, vector<double> &direction,
                       double &num_samples)
{
  cout << "sample size: " << num_samples << endl;
  Component component(direction,num_samples,parameters.constrain_kappa);
  component.minimizeMessageLength();
}

/*!
 *  \brief This function models a mixture of several components.
 *  \param parameters a reference to a struct Parameters
 *  \param data a reference to a vector<vector<double,3>>
 */
void modelMixture(struct Parameters &parameters, vector<vector<double>> &data)
{
  // if the optimal number of components need to be determined
  if (parameters.infer_num_components == SET) {
    vector<double> msglens;
    vector<int> components;
    for (int i=1; i<=20; i++) {
      //if (i % 2 == 1) {
        int k = 5 * i;
        cout << "Running for K: " << k << endl;
        components.push_back(k);
        Mixture mixture(k,data,parameters.update_weights_new,
                        parameters.constrain_kappa,parameters.simulation);
        if (parameters.dssp == SET) {
          mixture.setDSSPFlag(parameters.dssp_sst_type);
        }
        double msg = mixture.estimateParameters();
        msglens.push_back(msg);
      //}
    }
    plotMessageLengthAgainstComponents(components,msglens,parameters);
  } else if (parameters.infer_num_components == UNSET) {
    // for a given value of number of components
    // do the mixture modelling
    Mixture mixture(parameters.fit_num_components,data,
                    parameters.update_weights_new,parameters.constrain_kappa,
                    parameters.simulation);
    if (parameters.dssp == SET) {
      mixture.setDSSPFlag(parameters.dssp_sst_type);
    }
    mixture.estimateParameters();
  }
}

/*!
 *  \brief This function reads through the profiles from a given directory.
 *  \param parameters a reference to a struct Parameters
 *  \param direction a reference to a vector<double>
 *  \param n a reference to a double
 *  \return success or failure of browsing through the profiles
 */
bool readProfiles(struct Parameters &parameters, vector<double> &direction,
                  double &n)
{
  path p(parameters.profiles_dir);
  cout << "path: " << p.string() << endl;
  if (exists(p)) { 
    if (is_directory(p)) { 
      vector<path> files; // store paths,
      copy(directory_iterator(p), directory_iterator(), back_inserter(files));
      cout << "# of profiles: " << files.size() << endl;
      vector<vector<int>> bins;
      if (parameters.heat_map == SET) {
        int num_rows = 180 / parameters.res;
        int num_cols = 360 / parameters.res;
        cout << "rows: " << num_rows << endl;
        cout << "cols: " << num_cols << endl;
        vector<int> tmp(num_cols,0);
        for (int i=0; i<num_rows; i++) {
          bins.push_back(tmp);
        }
      }
      //vector<double> direction(3,0);
      //double n = 0;
      /*ofstream log("directions");
      for (int j=0; j<3; j++) {
        log << scientific << direction[j] << "\t";
      }
      log << endl;*/
      for (int i=0; i<files.size(); i++) {
        Protein protein;
        protein.load(files[i]);
        updateMeanDirection(direction,n,protein);
        if (parameters.heat_map == SET) {
          updateBins(bins,parameters.res,protein);
        }
        /*log << files[i].string() << "\t";
        for (int j=0; j<3; j++) {
          log << scientific << direction[j] << "\t";
        }
        log << endl;*/
      }
      //log.close();
      if (parameters.heat_map == SET) {
        outputBins(bins,parameters);
      }
      return 1;
    } else {
      cout << p << " exists, but is neither a regular file nor a directory\n";
    }
  } else {
    cout << p << " does not exist\n";
  }
  return 0;
}

/*!
 *  \brief This function reads through the profiles from a given directory
 *  and collects the data to do mixture modelling.
 *  \param parameters a reference to a struct Parameters
 *  \param unit_coordinates a reference to a vector<vector<double>>
 */
bool gatherData(struct Parameters &parameters, 
                vector<vector<double>> &unit_coordinates)
{
  path p(parameters.profiles_dir);
  cout << "path: " << p.string() << endl;
  if (exists(p)) { 
    if (is_directory(p)) { 
      vector<path> files; // store paths,
      copy(directory_iterator(p), directory_iterator(), back_inserter(files));
      cout << "# of profiles: " << files.size() << endl;
      for (int i=0; i<files.size(); i++) {
        Protein protein;
        protein.load(files[i]);
        vector<vector<double>> coords = protein.getUnitCoordinatesList();
        for (int j=0; j<coords.size(); j++) {
          unit_coordinates.push_back(coords[j]);
        }
      }
      return 1;
    } else {
      cout << p << " exists, but is neither a regular file nor a directory\n";
    }
  } else {
    cout << p << " does not exist\n";
  }
  return 0;
}

/*!
 *  \brief This function updates the run time
 *  \param name a reference to a string
 *  \param time a double
 */
void updateLogFile(string &name, double time, int num_chains)
{
  ofstream log("runtime_statistics",ios::app);
  log << name;
  log << fixed << setw(10) << setprecision(4) << time;
  log << fixed << setw(10) << num_chains;
  log << endl;
  log.close();
}

/*!
 *  \brief This function updates the mean estimator of the Von Mises
 *  distribution.
 *  \param direction a reference to a vector<double,3>
 *  \param n a pointer to an integer 
 *  \param protein a reference to a Protein object.
 */
void updateMeanDirection(vector<double> &direction, double &n, Protein &protein)
{
  vector<double> mean = protein.computeMeanDirection();
  for (int i=0; i<3; i++) {
    direction[i] += mean[i];
  }
  n += protein.getNumberOfSphericalCoordinates();
  //cout << *n << endl;
}

/*!
 *  \brief This function updates the frequencies of angles.
 *  \param bins a reference to a vector<vector<int>>
 *  \param res a double
 *  \param protein a reference to a Protein object.
 */
void updateBins(vector<vector<int>> &bins, double res, Protein &protein)
{
  vector<vector<double>> unit_coordinates = protein.getUnitCoordinatesList();
  double theta,phi;
  int row,col;
  vector<double> spherical(3,0);
  for (int i=0; i<unit_coordinates.size(); i++) {
    //cout << "i: " << i << endl; 
    cartesian2spherical(unit_coordinates[i],spherical);
    theta = spherical[1] * 180 / PI;
    if (fabs(theta) <= ZERO) {
      row = 0;
    } else {
      row = (int)(ceil(theta/res) - 1);
    }
    phi = spherical[2] * 180 / PI;
    if (fabs(phi) <= ZERO) {
      col = 0;
    } else {
      col = (int)(ceil(phi/res) - 1);
    }
    if (row >= bins.size() || col >= bins[0].size()) {
      cout << "outside bounds: " << row << " " << col << "\n";
      cout << "theta: " << theta << " phi: " << phi << endl;
      cout << "spherical_1: " << spherical[1] << " spherical_2: " << spherical[2] << endl;
      cout << "unit_coordinates[i]_1: " << unit_coordinates[i][1] << " unit_coordinates[i]_2: " << unit_coordinates[i][2] << endl;
      fflush(stdout);
    }
    bins[row][col]++;
    //cout << "row,col: " << row << "," << col << endl;
  }
}

/*!
 *  \brief This function outputs the bin data.
 *  \param bins a reference to a vector<vector<int>>
 *  \param parameters a reference to a struct Parameters
 */
void outputBins(vector<vector<int>> &bins, struct Parameters &parameters)
{
  double theta=0,phi;
  string fbins2D_file,fbins3D_file;
  if (parameters.dssp == UNSET) {
    fbins2D_file = CURRENT_DIRECTORY + "/matlab/bins_2D";
    fbins3D_file = CURRENT_DIRECTORY + "/matlab/bins_3D";
  } else if (parameters.dssp == SET) {
    string tmp;
    if (parameters.dssp_data_collect == 1) {
      tmp = CURRENT_DIRECTORY + "/dssp/parsed1/models/" + parameters.dssp_sst_type
            + "/matlab/";
    } else if (parameters.dssp_data_collect == 2) {
      tmp = CURRENT_DIRECTORY + "/dssp/parsed2/models/" + parameters.dssp_sst_type
            + "/matlab/";
    }
    fbins2D_file = tmp + "bins_2D";
    fbins3D_file = tmp +"bins_3D";
  }
  ofstream fbins2D(fbins2D_file.c_str());
  ofstream fbins3D(fbins3D_file.c_str());
  vector<double> cartesian(3,0);
  vector<double> spherical(3,1);
  for (int i=0; i<bins.size(); i++) {
    phi = 0;
    spherical[1] = theta * PI / 180;
    for (int j=0; j<bins[i].size(); j++) {
      fbins2D << fixed << setw(10) << bins[i][j];
      phi += parameters.res;
      spherical[2] = phi * PI / 180;
      spherical2cartesian(spherical,cartesian);
      for (int k=0; k<3; k++) {
        fbins3D << fixed << setw(10) << setprecision(4) << cartesian[k];
      }
      fbins3D << fixed << setw(10) << bins[i][j] << endl;
    }
    theta += parameters.res;
    fbins2D << endl;
  }
  fbins2D.close();
  fbins3D.close();
}

/*!
 *  \brief This function generates the data to visualize the mixture components.
 *  \param parameters a reference to a struct Parameters
 */
void visualizeMixtureComponents(struct Parameters &parameters)
{
  Mixture mixture;
  mixture.load(parameters.mixture_file);
  bool save = 1;
  if (parameters.sample_generation == USING_MIXTURE_WEIGHTS) {
    mixture.generateProportionally(parameters.num_samples,save);
  } else if (parameters.sample_generation == RANDOM_SAMPLE_SIZE) {
    mixture.generateRandomSampleSize(save);
  }
  if (parameters.heat_map == SET) {
    mixture.generateHeatmapData(parameters.res);
  }
}

/*!
 *  \brief This function is used to simulate the mixture model.
 *  \param parameters a reference to a struct Parameters
 */
void simulateMixtureModel(struct Parameters &parameters)
{
  vector<vector<double>> data;
  if (parameters.load_mixture == SET) {
    Mixture original;
    original.load(parameters.mixture_file);
    bool save = 0;
    if (parameters.sample_generation == USING_MIXTURE_WEIGHTS) {
      data = original.generateProportionally(parameters.num_samples,save);
    } else if (parameters.sample_generation == RANDOM_SAMPLE_SIZE) {
      data = original.generateRandomSampleSize(save);
    }
  } else {
    // # of components
    int K = parameters.simulate_num_components;

    // generate random weights
    vector<double> weights = generateRandomWeights(K,1);

    // generate random components
    vector<Component> 
    components = generateRandomComponents(K,parameters.constrain_kappa);

    // original mixture model
    Mixture original(K,components,weights);
    original.setSimulationFlag();
    data = original.generateProportionally(parameters.num_samples,0);
    original.printParameters();
  }

  // model a mixture using the original data
  modelMixture(parameters,data);
}

/*!
 *  \brief This function is used to generate a list of random weights.
 *  \param num_weights an integer
 *  \param range a double
 *  \return the list of weights
 */
vector<double> generateRandomWeights(int num_weights, double range)
{
  auto ts = high_resolution_clock::now();
  usleep(10);
  auto te = high_resolution_clock::now();
  double t = duration_cast<nanoseconds>(ts-te).count();
  srand(t);
  vector<double> weights;
  //double range = 1;
  for (int i=0; i<num_weights-1; i++) {
    double w = (rand() / (double) RAND_MAX) * range;
    assert(w > 0 && w < range);
    range -= w;
    weights.push_back(w);
  }
  weights.push_back(range);
  return weights;
}

/*!
 *  \brief This function is used to generate random components.
 *  \param num_components an integer
 *  \return the list of components
 *  \param constrain_kappa an integer
 */
vector<Component>
generateRandomComponents(int num_components, int constrain_kappa)
{
  auto ts = high_resolution_clock::now();
  usleep(10);
  auto te = high_resolution_clock::now();
  double t = duration_cast<nanoseconds>(ts-te).count();
  srand(t);
  vector<double> spherical_mean(3,1);
  vector<double> unit_mean(3,0);
  vector<Component> components;
  for (int i=0; i<num_components; i++) {
    // initialize component parameters
    spherical_mean[1] = (rand()/(double)RAND_MAX)*PI;
    spherical_mean[2] = (rand()/(double)RAND_MAX)*2*PI;
    spherical2cartesian(spherical_mean,unit_mean);
    double kappa = (rand() / (double) RAND_MAX) * MAX_KAPPA;
    Component component(unit_mean,kappa);
    component.printParameters(cout);
    components.push_back(component);
  }
  return components;
}

/*!
 *  \brief This function is used to plot the message lengths for different
 *  number of components.
 *  \param components a reference to a vector<int>
 *  \param msglens a reference to a vector<double>
 *  \param parameters a reference to a struct Parameters 
 */
void plotMessageLengthAgainstComponents(vector<int> &components,
                                        vector<double> &msglens, 
                                        struct Parameters &parameters)
{
  assert(components.size() == msglens.size());

  string data_file,plot_file,parsed;

  if (parameters.dssp == UNSET) {
    data_file = string(CURRENT_DIRECTORY) + "/mixture/";
    plot_file = string(CURRENT_DIRECTORY) + "/mixture/";
    if (parameters.simulation == SET) {
      data_file += "simulation/";
      plot_file += "simulation/";
    }
  } else if (parameters.dssp == SET) {
    if (DSSP_DATA_COLLECT == 1) {
      parsed = "parsed1";
    } else if (DSSP_DATA_COLLECT == 2) {
      parsed = "parsed2";
    }
    data_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/" + parameters.dssp_sst_type
                + "/mixture/";
    plot_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/" + parameters.dssp_sst_type
                + "/mixture/";
  }
  data_file += "msglens-infer.dat";
  plot_file += "msglens-infer.eps";
  ofstream file(data_file.c_str());
  for (int i=0; i<msglens.size(); i++) {
    file << components[i] << "\t" << msglens[i] << endl;
  }
  file.close();

  ofstream script("script.p");
	script << "# Gnuplot script file for plotting data in file \"data\"\n\n" ;
	script << "set terminal post eps" << endl ;
	script << "set autoscale\t" ;
	script << "# scale axes automatically" << endl ;
	script << "set xtic auto\t" ;
	script << "# set xtics automatically" << endl ;
	script << "set ytic auto\t" ;
	script << "# set ytics automatically" << endl ;
  script << "set xr [0:]" << endl;
	//script << "set title \"# of components: " << K << "\"" << endl ;
	script << "set xlabel \"# of components\"" << endl ;
	script << "set ylabel \"message length (in bits)\"" << endl ;
	script << "set output \"" << plot_file << "\"" << endl ;
	script << "plot \"" << data_file << "\" using 1:2 notitle " 
         << "with linespoints lc rgb \"red\"" << endl ;
  script.close();
  system("gnuplot -persist script.p");	
}

////////////////////////// SST FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void assignSecondaryStructure(string mixture_file, string structure_file,
                              int orientation, int portion_to_fit,
                              vector<string> &end_points, int sst_method)
{
  cout << "Assigning secondary structure to " << structure_file << endl;

  // read protein coordinate data
  string name = extractName(structure_file);
  ProteinStructure *p = parsePDBFile(structure_file);
  vector<ProteinStructure *> partitions = createPartitions(p);
  Protein protein(p,name);
  int num_residues = p->getNumberOfResidues();
  cout << "Number of residues: " << num_residues << endl;

  // compute the message length to transmit using the sphere approach
  protein.computeSuccessiveDistances();
  double sphere_msglen = protein.computeMessageLengthUsingSphereModel();

  // compute the message length to tansmit using the null model
  // null model: mixture model
  // read mixture data
  Mixture mixture;
  mixture.load(mixture_file);
  protein.computeSphericalTransformation();
  protein.getUnitCoordinatesList();
  double null_msglen = protein.computeMessageLengthUsingNullModel(mixture);

  // compute the message length using the compression model
  // using ideal models
  clock_t c_start = clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  protein.compressUsingIdealModels(mixture,orientation,portion_to_fit,end_points,sst_method);
  cout << "Sphere model message length: " << sphere_msglen << " bits. (" 
       << sphere_msglen / (double)num_residues << " bpr)" << endl;
  cout << "Null model message length: " << null_msglen << " bits. (" 
       << null_msglen / (double)num_residues << " bpr)" << endl;
  clock_t c_end = clock();
  auto t_end = std::chrono::high_resolution_clock::now();
  double cpu_time = double(c_end-c_start)/(double)(CLOCKS_PER_SEC);
  double wall_time = std::chrono::duration_cast<std::chrono::seconds>(t_end-t_start).count();
  cout << "CPU time: " << cpu_time << " secs." << endl;
  cout << "Wall time: " << wall_time << " secs." << endl;
}

/*!
 *  \brief This function loads the ideal models to be used in the compression
 *  model.
 *  \return the vector of ideal models
 */
vector<IdealModel> loadIdealModels()
{
  vector<IdealModel> ideal_models;
  string name,path;
  int num_residues;

  // load idealAlphaHelix_LH
  name = "AlphaHelix_LH";
  path = "./ideal_models/idealAlphaHelix_LH.pdb";
  ProteinStructure *alpha_lh = parsePDBFile(path);
  num_residues = alpha_lh->getNumberOfResidues();
  //cout << "# residues (alpha_lh): " << num_residues << endl;
  IdealModel m0(alpha_lh,name);
  ideal_models.push_back(m0);

  // load idealAlphaHelix_RH
  name = "AlphaHelix_RH";
  path = "./ideal_models/idealAlphaHelix_RH.pdb";
  ProteinStructure *alpha_rh = parsePDBFile(path);
  num_residues = alpha_rh->getNumberOfResidues();
  //cout << "# residues (alpha_rh): " << num_residues << endl;
  IdealModel m1(alpha_rh,name);
  ideal_models.push_back(m1);

  // load idealPiHelix_LH
  name = "PiHelix_LH";
  path = "./ideal_models/idealPiHelix_LH.pdb";
  ProteinStructure *pi_lh = parsePDBFile(path);
  num_residues = pi_lh->getNumberOfResidues();
  //cout << "# residues (pi_lh): " << num_residues << endl;
  IdealModel m2(pi_lh,name);
  ideal_models.push_back(m2);

  // load idealPiHelix_RH
  name = "PiHelix_RH";
  path = "./ideal_models/idealPiHelix_RH.pdb";
  ProteinStructure *pi_rh = parsePDBFile(path);
  num_residues = pi_rh->getNumberOfResidues();
  //cout << "# residues (pi_rh): " << num_residues << endl;
  IdealModel m3(pi_rh,name);
  ideal_models.push_back(m3);

  // load ideal310Helix_LH
  name = "310Helix_LH";
  path = "./ideal_models/ideal310Helix_LH.pdb";
  ProteinStructure *three10_lh = parsePDBFile(path);
  num_residues = three10_lh->getNumberOfResidues();
  //cout << "# residues (three10_lh): " << num_residues << endl;
  IdealModel m4(three10_lh,name);
  ideal_models.push_back(m4);

  // load ideal310Helix_RH
  name = "310Helix_RH";
  path = "./ideal_models/ideal310Helix_RH.pdb";
  ProteinStructure *three10_rh = parsePDBFile(path);
  num_residues = three10_rh->getNumberOfResidues();
  //cout << "# residues (three10_rh): " << num_residues << endl;
  IdealModel m5(three10_rh,name);
  ideal_models.push_back(m5);

  // load idealParallelBetaStrand
  name = "ParallelBetaStrand";
  path = "./ideal_models/idealParallelBetaStrand.pdb";
  ProteinStructure *beta_strand = parsePDBFile(path);
  num_residues = beta_strand->getNumberOfResidues();
  //cout << "# residues (beta_strand): " << num_residues << endl;
  IdealModel m6(beta_strand,name);
  ideal_models.push_back(m6);

  // load idealAntiParallelBetaStrand
  /*name = "AntiParallelBetaStrand";
  path = "./ideal_models/idealAntiParallelBetaStrand.pdb";
  ProteinStructure *beta_strand_anti = parsePDBFile(path);
  num_residues = beta_strand_anti->getNumberOfResidues();
  //cout << "# residues (beta_strand): " << num_residues << endl;
  IdealModel m7(beta_strand_anti,name);
  ideal_models.push_back(m7);*/

  return ideal_models;
}

/*!
 *  \brief This function assigns a representative mixture component to each of
 *  the ideal models.
 *  \param ideal_models a reference to a vector<IdealModel>
 *  \param components a reference to a vector<Component>
 *  \param weights a reference to a vector<double>
 *  \return the assignment list 
 */
vector<int> assignMixtureComponents(vector<IdealModel> &ideal_models, 
                                    vector<Component> &components, 
                                    vector<double> &weights)
{
  vector<double> mean(3,0),spherical(3,0),unit_mean(3,0);
  vector<double> density(components.size(),0);
  vector<int> assignment(ideal_models.size(),0);
  int max_index;
  double cumulative_weight = 0;

  for (int i=0; i<ideal_models.size(); i++) {
    ProteinStructure *p = ideal_models[i].getStructure();
    string name = ideal_models[i].getName();
    Protein protein(p,name);
    protein.computeSphericalTransformation();
    //protein.save();
    cout << name << endl;
    mean = protein.computeMeanDirection();
    //cout << "mean: "; print(cout,mean);
    cartesian2spherical(mean,spherical);
    //cout << "\tspherical mean of ideal model: "; print(cout,spherical);
    cout << "\t(theta,phi): (" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")";
    cartesian2unitspherical(mean,unit_mean);
    //cout << "unit mean: "; print(cout,unit_mean);
    string file = "density_" + name;
    ofstream out(file.c_str());
    double pr = 0;
    for (int j=0; j<components.size(); j++) {
      density[j] = weights[j] * components[j].likelihood(unit_mean);
      pr += density[j];
    }
    vector<int> sorted = sortedListIndex(density);
    for (int j=0; j<components.size(); j++) {
      out << sorted[j]+1 << "\t\t" << weights[j] << "\t\t" << density[sorted[j]]/pr << endl;
    }
    out.close();
    max_index = getIndexOfMaximumElement(density);
    assignment[i] = max_index;
    cumulative_weight += weights[max_index];
    cout << "\n\tassigned component: " << max_index+1;
    cout << "\tprob density: " << density[max_index]/pr << "\t" ; 
    cout << "\tweight: " << weights[max_index] << "\t";
    components[max_index].printParameters(cout);
  }
  cout << "Cumulative weight: " << cumulative_weight << endl;
  return assignment;
}

/*!
 *  \brief This function is used to load the ideal mixture models.
 *  \return the list of ideal mixtures
 */
vector<Mixture> loadIdealMixtureModels()
{
  vector<Mixture> ideal_mixtures;
  string mixture_file,dssp_sst_type,parsed;

  if (DSSP_DATA_COLLECT == 1) {
    parsed = "parsed1";
  } else if (DSSP_DATA_COLLECT == 2) {
    parsed = "parsed2";
  }

  // load ideal alpha-helix
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/50_helix_alpha.mixture";
  mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/helix_alpha.mixture";
  Mixture m0;
  m0.load(mixture_file);
  dssp_sst_type = "helix_alpha";
  m0.setDSSPFlag(dssp_sst_type);
  ideal_mixtures.push_back(m0);

  // load ideal pi-helix
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/50_helix_pi.mixture";
  mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/helix_pi.mixture";
  Mixture m1;
  m1.load(mixture_file);
  dssp_sst_type = "helix_pi";
  m1.setDSSPFlag(dssp_sst_type);
  ideal_mixtures.push_back(m1);

  // load ideal 310-helix
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed +"/models/ideal_mixture_models/50_helix_310.mixture";
  mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed +"/models/ideal_mixture_models/helix_310.mixture";
  Mixture m2;
  m2.load(mixture_file);
  dssp_sst_type = "helix_310";
  m2.setDSSPFlag(dssp_sst_type);
  ideal_mixtures.push_back(m2);

  // load ideal strand 
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/50_strand.mixture";
  mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/strand.mixture";
  Mixture m3;
  m3.load(mixture_file);
  dssp_sst_type = "strand";
  m3.setDSSPFlag(dssp_sst_type);
  ideal_mixtures.push_back(m3);

  // load ideal coil 
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/coil_blanks.mixture";
  mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/coil_binary_bins.mixture";
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/50_coil.mixture";
  //mixture_file = CURRENT_DIRECTORY + "/dssp/" + parsed + "/models/ideal_mixture_models/coil.mixture";
  Mixture m4;
  m4.load(mixture_file);
  dssp_sst_type = "coil";
  m4.setDSSPFlag(dssp_sst_type);
  ideal_mixtures.push_back(m4);

  return ideal_mixtures;
}

/*!
 *  \brief This function computes the relative weights of the mixtures based
 *  on their sample size.
 *  \param mixtures a reference to a vector<Mixture>
 *  \return the relative weights
 */
vector<double> computeRelativeWeights(vector<Mixture> &mixtures)
{
  vector<double> weights(mixtures.size(),0);
  double sum = 0;
  for (int i=0; i<mixtures.size(); i++) {
    weights[i] = mixtures[i].getSampleSize();
    sum += weights[i];
  }
  for (int i=0; i<mixtures.size(); i++) {
    //weights[i] /= sum;
    weights[i] = 1/(double)mixtures.size();
  }
  return weights;
}

/*!
 *  \brief This function returns the mean length of an ideal model.
 *  \param name a string
 *  \return the mean length
 */
double getMeanLength(string name)
{
  /*double MEAN_HELIX_ALPHA = 6;
  double MEAN_STRAND = 5;
  double MEAN_COIL = 2;*/

  double MEAN_COIL = 5.18;
  double MEAN_STRAND = 5.40;
  double MEAN_HELIX_310 = 3.36;
  double MEAN_HELIX_ALPHA = 11.21;
  double MEAN_HELIX_PI = 5.24;

  if (name.compare("strand") == 0) {
    return MEAN_STRAND;
  } else if (name.compare("coil") == 0) {
    return MEAN_COIL;
  } else if (name.compare("helix_alpha") == 0) {
    return MEAN_HELIX_ALPHA;
  } else if (name.compare("helix_310") == 0) {
    return MEAN_HELIX_310;
  } else if (name.compare("helix_pi") == 0) {
    return MEAN_HELIX_PI;
  }
}

/*!
 *  \brief This function collects all the continuous stretches within a 
 *  protein chain. 
 *  \param id a reference to a string
 *  \param atoms a reference to a vector<Atoms> 
 *  \return whether there are chain breks or not
 */
vector<array<int,2>> identifyChainBreaks(string &id, vector<Atom> &atoms)
{
  vector<array<int,2>> stretches;
  array<int,2> one_stretch;
  int start = 1,end;
  for (int i = 1; i < atoms.size(); ++i) {
    double dist = distance<double>(atoms[i-1], atoms[i]);
    if (dist > 4) {
      end = i - 1;
      one_stretch[0] = start;
      one_stretch[1] = end;
      stretches.push_back(one_stretch);
      cout << "Break in chain " << id << " ...\n";
    }
  }
  return stretches;
}

/*!
 *  \brief This function is used to create partitions within the protein.
 *  \param original a pointer to a ProteinStructure
 *  \return the list of partitions
 */
vector<ProteinStructure *> createPartitions(ProteinStructure *original)
{
  vector<ProteinStructure *> partitions;
  // get all chains
  vector<string> chain_ids = original->getChainIdentifiers();
  for (int i=0; i<chain_ids.size(); i++) {
    string id = chain_ids[i];
    Chain chain = original->getDefaultModel()[id];
    vector<Atom> atoms = chain.getAtoms();
    vector<array<int,2>> stretches = identifyChainBreaks(id,atoms);
  }
  return partitions;
}

