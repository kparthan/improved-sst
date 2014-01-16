#include "Support.h"
#include "VonMises3D.h"
#include "Mixture.h"

int initialize_components_from_file;

//////////////////////// GENERAL PURPOSE FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This function checks to see if valid arguments are given to the 
 *  command line output.
 *  \param argc an integer
 *  \param argv an array of strings
 *  \return the parameters of the model
 */
struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string pdb_id,scop_id,generation,constrain;

  bool noargs = 1;

  cout << "Checking command-line input ..." << endl;
  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("test","perform a demo")
       ("verbose","print some details")
       ("file",value<string>(&parameters.file),"pdb file")
       ("pdbid",value<string>(&pdb_id),"PDB ID")
       ("scopid",value<string>(&scop_id),"SCOP ID")
       ("force","force build the profile")
       ("directory",value<string>(&parameters.profiles_dir),
                                      "path to all profiles")
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

  if (vm.count("initialize_components_from_file")) {
    initialize_components_from_file = SET;
  } else {
    initialize_components_from_file = UNSET;
  }
  
  if (vm.count("directory")) {
    if (vm.count("file") || vm.count("pdbid") || vm.count("scopid")) {
      cout << "conflicting options: directory/file type ..." << endl;
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
  /*ifstream file(fileName);
  return file;*/
  return boost::filesystem::exists(file_name);
}

/*!
 *  \brief This module prints the list of coordinates to a file
 *  \param coordinates a reference to vector<Point<double>>
 *  \param file_name a pointer to a const char
 */
void writeToFile(vector<Point<double>> &coordinates, const char *file_name)
{
  ofstream file(file_name);
  for (int i=0; i<coordinates.size(); i++){
    file << coordinates[i] << endl;
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
 *  \brief This function computes the spherical coordinate values
 *  \param point a reference to a Point<double>
 *  \return the theta-phi values (in degrees) corresponding to a point on the
 *  sphere
 */
array<double,3> convertToSpherical(Point<double> &point)
{
  array<double,3> values;
  Point<double> origin(0,0,0);
  values[0] = distance<double>(origin,point);
  Vector<double> i_plus_1 = point.positionVector();
  array<double,2> angles = angleWithAxes<double>(i_plus_1); // angles in radians
  values[1] = angles[0] * (180/PI);
  values[2] = angles[1] * (180/PI);
  return values;
}

/*!
 *  \brief This function converts spherical coordinates into Cartesian
 *  coordinates.
 *  \param r a double
 *  \param theta (in degrees) a double
 *  \param phi (in degrees) a double
 *  \return the Cartesian coordinates
 */
array<double,3> convertToCartesian(double r, double theta, double phi)
{
  array<double,3> x;
  double theta_rad = theta * PI / 180;
  //scaleToAOM(&theta_rad);
  double phi_rad = phi * PI / 180;
  //scaleToAOM(&phi_rad);
  x[0] = r * sin(theta_rad) * cos(phi_rad);
  x[1] = r * sin(theta_rad) * sin(phi_rad);
  x[2] = r * cos(theta_rad);
  return x;
}

/*!
 *  \brief This function scales the value to AOM.
 *  \param angle_rad (in radians) a double
 *  \return angle to AOM radians
 */
void scaleToAOM(double *angle_rad)
{
  int scale = 1 / (double) AOM;
  int angle = (*angle_rad) * scale;
  *angle_rad = angle / (double) scale;
  cout << *angle_rad << endl;
}

/*!
 *  \brief This function prints the elements of an array.
 *  \param os a reference to a ostream
 *  \param a a reference to an array<double,3>
 */
void print(ostream &os, array<double,3> &a)
{
  os << "(" << a[0] << "," << a[1] << "," << a[2] << ")\n";
}

/*!
 *  \brief This function generates the heat map for theta vs phi
 *  \param estimates a reference to a array<double,3>
 *  \param res a double
 */
void vonMisesDistribution_2DPlot(array<double,3> &estimates, double res)
{
  double kappa = estimates[0];
  array<double,2> mu({estimates[1],estimates[2]});
  VonMises3D von_mises(mu,kappa);
  ofstream log("matlab/distribution_heat_map_2D.data");
  double density;
  for (double theta=0; theta<=180; theta+=res) {
    for (double phi=0; phi<=360; phi+=res) {
      //log << fixed << setw(10) << setprecision(2) << theta;
      //log << fixed << setw(10) << setprecision(2) << phi;
      density = von_mises.density(theta,phi);
      log << fixed << setw(10) << setprecision(4) << density;
    }
    log << endl;
  }
  log.close();
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
 *  \param theta (measured in degrees) a double
 *  \return the value in radians
 */
double angleInRadians(double theta)
{
  return theta * PI / 180;
}

//////////////////////// PROTEIN FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This module returns the file path associated with a PDB ID.
 *  \param pdb_id a reference to a string
 *  \return the file path
 */
string getPDBFilePath(string &pdb_id)
{
  boost::algorithm::to_lower(pdb_id);
  string path = string(HOME_DIRECTORY) + "Research/PDB/" ;
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
  string path = string(HOME_DIRECTORY) + "Research/SCOP/pdbstyle-1.75B/" ;
  string directory(scop_id,2,2);
  path += directory + "/" + scop_id + ".ent";
  return path;
}

/*!
 *  \brief This method constructs the angular profiles of the protein structure.
 *  \param parameters a reference to a struct Parameters
 */
void buildAngularProfile(struct Parameters &parameters)
{
  string name = extractName(parameters.file);
  // check if spherical profile exists
  bool exist = checkIfSphericalProfileExists(name);
  Protein protein;
  if (!exist || parameters.force == SET) {
    /* Obtain protein coordinates */
    ProteinStructure *p = parsePDBFile(parameters.file);
    Protein protein(p,name);
    protein.computeSphericalTransformation();
    cout << "Saving profile of " << name << endl;
    protein.save();
    //updateLogFile(name,protein.getCPUTime(),protein.getNumberOfChains());
  } else {  
    cout << "Profile of " << name << " exists ..." << endl;
    protein.load(name);
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
                             + "spherical_system/profiles/" + name + ".profile";
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
    ProteinStructure *structure = 
        parser.getStructure(pdb_file.c_str())->select(CASelector());
    ProteinStructure *one_model = 
        new ProteinStructure(structure->getIdentifier());
    one_model->select(CASelector());
    std::shared_ptr<lcb::Model> newmodel = 
        std::make_shared<lcb::Model>(structure->getDefaultModel());
    one_model->addModel(newmodel);
    delete structure;
    cout << " [OK]" << endl;
    return one_model;
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
    pair<array<double,3>,double> data = readProfiles(parameters);
    modelOneComponent(parameters,data);
  } else if (parameters.mixture_model == SET) { // mixture modelling
    vector<array<double,3>> data = gatherData(parameters);
    modelMixture(parameters,data);
  }
}

/*!
 *  \brief This function models a single component.
 *  \param parameters a reference to a struct Parameters
 *  \param data a reference to a pair<array<double,3>,double>
 */
void modelOneComponent(struct Parameters &parameters,
                       pair<array<double,3>,double> &data)
{
  array<double,3> direction = data.first;
  double num_samples = data.second;
  //vonMisesDistribution_2DPlot(direction,parameters.res);
  Component component(direction,num_samples,parameters.constrain_kappa);
  component.minimizeMessageLength();
}

/*!
 *  \brief This function models a mixture of several components.
 *  \param parameters a reference to a struct Parameters
 *  \param data a reference to a vector<array<double,3>>
 */
void modelMixture(struct Parameters &parameters, vector<array<double,3>> &data)
{
  // if the optimal number of components need to be determined
  if (parameters.infer_num_components == SET) {
    vector<double> msglens;
    vector<int> components;
    for (int i=2; i<=16; i++) {
      if (i % 3 == 2) {
        int k = 10 * i;
        cout << "Running for K: " << k << endl;
        components.push_back(k);
        Mixture mixture(k,data,parameters.update_weights_new,
                        parameters.constrain_kappa,parameters.simulation);
        double msg = mixture.estimateParameters();
        msglens.push_back(msg);
      }
    }
    plotMessageLengthAgainstComponents(components,msglens,parameters.simulation);
  } else if (parameters.infer_num_components == UNSET) {
    // for a given value of number of components
    // do the mixture modelling
    Mixture mixture(parameters.fit_num_components,data,
                    parameters.update_weights_new,parameters.constrain_kappa,
                    parameters.simulation);
    mixture.estimateParameters();
  }
}

/*!
 *  \brief This function reads through the profiles from a given directory.
 *  \param parameters a reference to a struct Parameters
 *  \return a pair of mean direction and the sample size
 */
pair<array<double,3>,double> readProfiles(struct Parameters &parameters)
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
        for (int i=0; i<num_rows; i++) {
          vector<int> tmp(num_cols,0);
          bins.push_back(tmp);
        }
      }
      array<double,3> direction({0,0,0});
      double n = 0;
      /*ofstream log("directions");
      for (int j=0; j<3; j++) {
        log << scientific << direction[j] << "\t";
      }
      log << endl;*/
      for (int i=0; i<files.size(); i++) {
        Protein protein;
        protein.load(files[i]);
        updateMeanDirection(direction,&n,protein);
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
        outputBins(bins,parameters.res);
      }
      return pair<array<double,3>,double>(direction,n);
    } else {
      cout << p << " exists, but is neither a regular file nor a directory\n";
    }
  } else {
    cout << p << " does not exist\n";
  }
  exit(1);
}

/*!
 *  \brief This function reads through the profiles from a given directory
 *  and collects the data to do mixture modelling.
 *  \param parameters a reference to a struct Parameters
 *  \return a pair of mean direction and the sample size
 */
vector<array<double,3>> gatherData(struct Parameters &parameters)
{
  path p(parameters.profiles_dir);
  cout << "path: " << p.string() << endl;
  if (exists(p)) { 
    if (is_directory(p)) { 
      vector<path> files; // store paths,
      copy(directory_iterator(p), directory_iterator(), back_inserter(files));
      cout << "# of profiles: " << files.size() << endl;
      vector<array<double,3>> coordinates;
      for (int i=0; i<files.size(); i++) {
        Protein protein;
        protein.load(files[i]);
        vector<array<double,3>> coords = protein.getSphericalCoordinatesList();
        for (int j=0; j<coords.size(); j++) {
          double theta = coords[j][1];
          double phi = coords[j][2];
          array<double,3> x = convertToCartesian(1,theta,phi);
          coordinates.push_back(x);
        }
      }
      return coordinates;
    } else {
      cout << p << " exists, but is neither a regular file nor a directory\n";
    }
  } else {
    cout << p << " does not exist\n";
  }
  exit(1);
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
 *  \param direction a reference to a array<double,3>
 *  \param n a pointer to an integer 
 *  \param protein a reference to a Protein object.
 */
void updateMeanDirection(array<double,3> &direction, double *n, Protein &protein)
{
  array<double,3> mean = protein.computeMeanDirection();
  for (int i=0; i<3; i++) {
    direction[i] += mean[i];
  }
  *n += protein.getNumberOfSphericalCoordinates();
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
  vector<array<double,3>> spherical_coordinates = protein.getSphericalCoordinatesList();
  double theta,phi;
  int row,col;
  for (int i=0; i<spherical_coordinates.size(); i++) {
    theta = spherical_coordinates[i][1];
    if (fabs(theta) <= ZERO) {
      row = 0;
    } else {
      row = (int)(ceil(theta/res) - 1);
    }
    phi = spherical_coordinates[i][2];
    if (fabs(phi) <= ZERO) {
      col = 0;
    } else {
      col = (int)(ceil(phi/res) - 1);
    }
    if (row >= bins.size() || col >= bins[0].size()) {
      cout << "outside bounds: " << row << " " << col << "\n";
      cout << "theta: " << theta << " phi: " << phi << endl;
      cout.flush();
    }
    bins[row][col]++;
    //cout << "row,col: " << row << "," << col << endl;
  }
}

/*!
 *  \brief This function outputs the bin data.
 *  \param bins a reference to a vector<vector<int>>
 *  \param res a double
 */
void outputBins(vector<vector<int>> &bins, double res)
{
  double theta=0,phi;
  ofstream fbins("matlab/bins_2D");
  ofstream fbins3D("matlab/bins_3D");
  for (int i=0; i<bins.size(); i++) {
    phi = 0;
    for (int j=0; j<bins[i].size(); j++) {
      fbins << fixed << setw(10) << bins[i][j];
      phi += res;
      array<double,3> point = convertToCartesian(1,theta,phi);
      for (int k=0; k<3; k++) {
        fbins3D << fixed << setw(10) << setprecision(4) << point[k];
      }
      fbins3D << fixed << setw(10) << bins[i][j] << endl;
    }
    theta += res;
    fbins << endl;
  }
  fbins.close();
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
  vector<array<double,3>> data;
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
    Mixture original(K,weights,components);
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
  vector<Component> components;
  for (int i=0; i<num_components; i++) {
    // initialize component parameters
    array<double,2> mu;
    mu[0] = (rand()/(double)RAND_MAX)*180;
    mu[1] = (rand()/(double)RAND_MAX)*360;
    double kappa = (rand() / (double) RAND_MAX) * MAX_KAPPA;
    Component component(mu,kappa,constrain_kappa);
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
 *  \param simulation an integer
 */
void plotMessageLengthAgainstComponents(vector<int> &components,
                                        vector<double> &msglens, int simulation)
{
  assert(components.size() == msglens.size());
  // output the data to a file
  string data_file = string(CURRENT_DIRECTORY) + "mixture/";
  if (simulation == SET) {
    data_file += "simulation/";
  }
  data_file += "msglens-infer2.dat";
  ofstream file(data_file.c_str());
  for (int i=0; i<msglens.size(); i++) {
    file << components[i] << "\t" << msglens[i] << endl;
  }
  file.close();

  // prepare gnuplot script file
  string output_file = string(CURRENT_DIRECTORY) + "mixture/";
  if (simulation == SET) {
    output_file += "simulation/";
  }
  output_file += "msglens-infer.eps";
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
	script << "set output \"" << output_file << "\"" << endl ;
	script << "plot \"" << data_file << "\" using 1:2 notitle " 
         << "with linespoints lc rgb \"red\"" << endl ;
  script.close();
  system("gnuplot -persist script.p") ;	
}

