#include "Support.h"
#include "VonMises3D.h"
#include "Message.h"

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
  string pdb_id,scop_id;

  bool noargs = 1;

  cout << "Checking command-line input ..." << endl;
  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help message")
       ("test","perform a demo")
       ("verbose","print some details")
       ("file",value<string>(&parameters.file),"pdb file")
       ("pdbid",value<string>(&pdb_id),"PDB ID")
       ("scopid",value<string>(&scop_id),"SCOP ID")
       ("force","force build the profile")
       ("profiles_dir",value<string>(&parameters.profiles_dir),
                                      "path to all profiles")
       ("bins","to compute the frequencies of angles")
       ("res",value<double>(&parameters.res),"heat map resolution")

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
  
  if (vm.count("profiles_dir")) {
    parameters.read_profiles = SET;
    if (vm.count("bins")) {
      parameters.update_bins = SET;
      if (!vm.count("res")) {
        parameters.res = DEFAULT_RESOLUTION;
      }
    } else {
      parameters.update_bins = UNSET;
    }
  } else {
    parameters.read_profiles = UNSET;
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
    } else {
      cout << "Input protein structure file not provided ..." << endl;
      Usage(argv[0],desc);
    }
  }
  noargs = 0;

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
 *  \return the spherical coordinates at a given index
 */
array<double,3> convertToSpherical(Point<double> &point)
{
  array<double,3> values;
  Point<double> origin(0,0,0);
  values[0] = distance<double>(origin,point);
  Vector<double> i_plus_1 = point.positionVector();
  array<double,2> angles = angleWithAxes<double>(i_plus_1);
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
  double phi_rad = phi * PI / 180;
  x[0] = r * sin(theta_rad) * cos(phi_rad);
  x[1] = r * sin(theta_rad) * sin(phi_rad);
  x[2] = r * cos(theta_rad);
  return x;
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
  double k_inv = 1 / (double)k;
  double cothk = 1 / (double) tanh(k);
  return cothk - k_inv;
}

/*!
 *  \brief This function computes the first derivative of the 
 *  ratio of Bessel functions -- A3'(k)
 *  \param k a double
 *  \return the value of the first derivative of the ratio at a given k
 */
double ratioBesselFunction_firstDerivative(double k)
{
  double value = 0;
  double cschk = 1 / (double) sinh(k);
  value -= cschk * cschk;
  double k_inv = 1 / (double)k;
  value += k_inv * k_inv; 
  return value ;
}

/*!
 *  \brief This function computes the second derivative of the 
 *  ratio of Bessel functions -- A3''(k)
 *  \param k a double
 *  \return the value of the second derivative of the ratio at a given k
 */
double ratioBesselFunction_secondDerivative(double k)
{
  double value = 0;
  double cothk = 1 / (double) tanh(k);
  double cschk = 1 / (double) sinh(k);
  value += cothk * cschk * cschk;
  double k_inv = 1 / (double)k;
  value -= k_inv * k_inv * k_inv;
  return 2 * value ;
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
  pair<array<double,3>,double> data = readProfiles(parameters);
  array<double,3> direction = data.first;
  double num_samples = data.second;
  //vonMisesDistribution_2DPlot(direction,parameters.res);
  Message message(direction,num_samples);
  message.minimize();
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
      if (parameters.update_bins == SET) {
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
        if (parameters.update_bins == SET) {
          updateBins(bins,parameters.res,protein);
        }
        /*log << files[i].string() << "\t";
        for (int j=0; j<3; j++) {
          log << scientific << direction[j] << "\t";
        }
        log << endl;*/
      }
      //log.close();
      if (parameters.update_bins == SET) {
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
  ofstream fbins("matlab/bins");
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

