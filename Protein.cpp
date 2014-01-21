#include "Support.h"
#include "Segment.h"

/*!
 *  \brief Null constructor module.
 */
Protein::Protein()
{}

/*!
 *  \brief This is a constructor function used to instantiate the Protein
 *  object from a ProteinStructure
 *  \param name a reference to a string
 *  \param structure a reference to a ProteinStructure
 */
Protein::Protein(ProteinStructure *structure, string &name) : 
                 structure(structure), name(name)
{
  vector<string> chain_ids = structure->getChainIdentifiers();
  for (int i=0; i<chain_ids.size(); i++) {
    string id = chain_ids[i];
    Chain chain = structure->getDefaultModel()[id];
    vector<Atom> atoms = chain.getAtoms();
    bool chain_break = checkChainBreak(id,atoms);
    if (!chain_break) {
      chains.push_back(id);
      vector<Point<double>> chain_coordinates;
      for (int j=0; j<atoms.size(); j++) {
        Point<double> p = atoms[j].point<double>();
        chain_coordinates.push_back(p);
      }
      translateProteinToOrigin(chain_coordinates);
      coordinates.push_back(chain_coordinates);
    }
  }
  cout << "# of suitable chains: " << coordinates.size() << endl;
  if (coordinates.size() == 0) {
    cout << name << " is an unsuitable structure ..." << endl;
    ofstream log("errors.log",ios::app);
    log << name << endl;
    log.close();
    exit(1);
  }
}

/*!
 *  \brief This function translates the protein so that its first point
 *  coincides with the origin
 *  \param chain_coordinates a reference to a vector<double>
 */
void Protein::translateProteinToOrigin(vector<Point<double>> &chain_coordinates)
{
  // before translation
  writeToFile(chain_coordinates,"before_translation");
  translation_vector.x(-chain_coordinates[0].x());
  translation_vector.y(-chain_coordinates[0].y());
  translation_vector.z(-chain_coordinates[0].z());
  cout << "Translation vector: " << translation_vector << endl;

  // translate the protein
  for (int i=0; i<chain_coordinates.size(); i++) {
    chain_coordinates[i] += translation_vector;
  }
  // after translation
  writeToFile(chain_coordinates,"after_translation");
}

/*!
 *  \brief This function checks if the pdb file has any chain breaks.
 *  \param id a reference to a string
 *  \param atoms a reference to a vector<Atoms> 
 *  \return whether there are chain breks or not
 */
bool Protein::checkChainBreak(string &id, vector<Atom> &atoms)
{
  for (unsigned i = 1; i < atoms.size(); ++i) {
    double dist = distance<double>(atoms[i-1], atoms[i]);
    if (dist > 4) {
      cout << "Break in chain " << id << " ...\n";
      return 1;
    }
  }
  return 0;
}

/*!
 *  \brief This function returns the number of chain breaks. 
 *  \param id a reference to a string
 *  \param atoms a reference to a vector<Atoms> 
 *  \return the number of chain breaks 
 */
int Protein::getNumberOfChainBreaks(string &id, vector<Atom> &atoms)
{
  int num_breaks = 0;
  for (unsigned i = 1; i < atoms.size(); ++i) {
    double dist = distance<double>(atoms[i-1], atoms[i]);
    if (dist > 4) {
      ++num_breaks;
    }
  }
  cout << "There are " << num_breaks << " breaks in chain " 
       << id << " of structure " << name  << endl;
  return num_breaks;
}

/*!
 *  \brief This function transforms the protein to a soherical coordinate
 *  system.
 */
void Protein::computeSphericalTransformation()
{
  clock_t c_start = clock();
  auto t_start = high_resolution_clock::now();

  for (int i=0; i<coordinates.size(); i++) {
    vector<array<double,3>> chain_coordinates;
    for (int j=2; j<coordinates[i].size()-1; j++) {
      vector<Point<double>> transformed_coordinates = computeTransformation(i,j);
      //string file_index = boost::lexical_cast<string>(j);
      //writeToFile(transformed_coordinates,file_index.c_str());
      array<double,3> values = convertToSpherical(transformed_coordinates[3]);
      chain_coordinates.push_back(values);
    }
    spherical_coordinates.push_back(chain_coordinates);
  }

  clock_t c_end = clock();
  auto t_end = high_resolution_clock::now();
  cpu_time = double(c_end-c_start)/(double)(CLOCKS_PER_SEC);
  wall_time = duration_cast<seconds>(t_end-t_start).count();
}

/*!
 *  \brief This function computes the spherical coordinate transformation to the
 *  canonical form at a coordinate index
 *  \param chain_index an integer
 *  \param index an integer
 *  \return the transformed list of four coordinates
 */
vector<Point<double>> Protein::computeTransformation(int chain_index, int index)
{
  vector<Point<double>> four_mer(4,Point<double>());
  four_mer[0] = coordinates[chain_index][index-2];
  four_mer[1] = coordinates[chain_index][index-1];
  four_mer[2] = coordinates[chain_index][index];
  four_mer[3] = coordinates[chain_index][index+1];
  pair<vector<Point<double>>,Matrix<double>> 
  transformation = convertToCanonicalForm(four_mer);
  return transformation.first; 
}

/*!
 *  \brief Loads the profile from an existing file.
 *  \param identifier a reference to a string
 */
void Protein::load(string &identifier)
{
  name = identifier;
  string file_name = string(CURRENT_DIRECTORY) + "spherical_system/profiles/"
                     + name + ".profile";
  read_profile(file_name);
}

/*!
 *  \brief Loads the profile from an existing file.
 *  \param path_to_file a reference to path object 
 */
void Protein::load(path &path_to_file)
{
  name = "example";
  string file_name = path_to_file.string();
  read_profile(file_name);
}

/*!
 *  \brief This functions reads the profile from a regular file.
 *  (theta,phi angles read are in degrees)
 *  \param file_name a reference to a string
 */
void Protein::read_profile(string &file_name)
{
  cout << "Reading " << file_name << " ..." << endl;
  ifstream profile(file_name.c_str());
  string line;
  spherical_coordinates.clear();
  all_spherical_coordinates.clear();

  vector<array<double,3>> chain_coordinates;
  while(getline(profile,line)) {
    boost::char_separator<char> sep(",() ");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    int i = -1;
    array<double,3> values; 
    BOOST_FOREACH (const string& t, tokens) {
      if (i == -1) {
        if (chains.size() == 0) {
          chains.push_back(t);
        } else {
          if (chains[chains.size()-1].compare(t) != 0) {
            chains.push_back(t);
            spherical_coordinates.push_back(chain_coordinates);
            chain_coordinates.clear();
          }
        }
        i++;
      } else {
        istringstream iss(t);
        double x;
        iss >> x;
        values[i++] = x;
      }
    }
    chain_coordinates.push_back(values);
  }
  spherical_coordinates.push_back(chain_coordinates);
  profile.close();
}

/*!
 *  \brief Saves the spherical system profile.
 *  (theta,phi angles saved are in degrees)
 */
void Protein::save()
{
  string file_name = string(CURRENT_DIRECTORY) + "spherical_system/profiles/"
                     + name + ".profile";
  ofstream profile(file_name.c_str());
  for (int i=0; i<spherical_coordinates.size(); i++) {
    for (int j=0; j<spherical_coordinates[i].size(); j++) {
      profile << chains[i];
      for (int k=0; k<3; k++) {
        profile << fixed << setw(10) << setprecision(4) 
                << spherical_coordinates[i][j][k];
      }
      profile << endl;
    }
  }
  profile.close();
}

/*!
 *  \brief This function gets the list of all spherical coordinates.
 *  \return the lsit of all spherical coordinates.
 */
vector<array<double,3>> Protein::getSphericalCoordinatesList()
{
  if (all_spherical_coordinates.size() == 0) {
    for (int i=0; i<spherical_coordinates.size(); i++) {
      for (int j=0; j<spherical_coordinates[i].size(); j++) {
        all_spherical_coordinates.push_back(spherical_coordinates[i][j]);
      }
    }
  }
  return all_spherical_coordinates;
}

/*!
 *  \brief This function returns the size of transformed spherical coordinates.
 *  \return the size
 */
int Protein::getNumberOfSphericalCoordinates()
{
  return all_spherical_coordinates.size();
}

/*!
 *  \brief This function gets the CPU time.
 *  \return the CPU time
 */
double Protein::getCPUTime()
{
  return cpu_time;
}

/*!
 *  \brief This function returns the number of usable chains in the structure.
 *  \return the number of chains
 */
int Protein::getNumberOfChains()
{
  return chains.size();
}

/*!
 *  \brief This function computes the mean direction vector as per the
 *  von Mises distribution.
 *  \return the mean direction
 */
array<double,3> Protein::computeMeanDirection()
{
  if (all_spherical_coordinates.size() == 0) {
    getSphericalCoordinatesList();
  }
  double r,theta,phi;
  array<double,3> estimate({0,0,0}),x;
  
  for (int i=0; i<all_spherical_coordinates.size(); i++) {
    r = all_spherical_coordinates[i][0];
    theta = all_spherical_coordinates[i][1];
    phi = all_spherical_coordinates[i][2];

    x = convertToCartesian(1,theta,phi);

    for (int j=0; j<3; j++) {
      estimate[j] += x[j];
    }
  }
  return estimate;
}

/*!
 *  \brief This function is used to compute the distance between successive
 *  residues of the protein.
 */
void Protein::computeSuccessiveDistances()
{
  for (int i=0; i<coordinates.size(); i++) {
    vector<double> dist;
    for (int j=0; j<coordinates[i].size()-1; j++) {
      double d = lcb::geometry::distance<double>(coordinates[i][j],coordinates[i][j+1]);
      dist.push_back(d);
    }
    distances.push_back(dist);
  }
}

/*!
 *  \brief This function computes the message length to communicate the protein
 *  coordinates using the sphere model.
 *  \return the message length
 */
double Protein::computeMessageLengthUsingSphereModel()
{
  double msglen = 0;

  // message length to state the number of chains
  int num_chains = chains.size(); // alternately spherical_coordinates.size()
  msglen += encodeUsingLogStarModel(num_chains);

  double constant = log2(4*PI) - 2*log2(AOM);
  for (int i=0; i<distances.size(); i++) {
    // for each chain state the number of residues
    int num_residues = distances[i].size();
    msglen += encodeUsingLogStarModel(num_residues);

    // state the residues
    vector<double> radii;
    for (int j=0; j<distances[i].size(); j++) {
      double r = distances[i][j];
      radii.push_back(r);
      msglen += 2 * log2(r);  // state the points on the surface of sphere
    }
    // state the points on the surface of sphere
    msglen += num_residues * constant;
    // collect the radii and send them together
    // state the radii
    msglen += encodeUsingNormalModel(radii);
  }  
  return msglen;
}

/*!
 *  \brief This function computes the message length to communicate the protein
 *  coordinates using the null model.
 *  \param mixture a reference to a Mixture
 *  \return the message length
 */
double Protein::computeMessageLengthUsingNullModel(Mixture &mixture)
{
  double msglen = 0;

  // message length to state the number of chains
  int num_chains = chains.size(); // alternately spherical_coordinates.size()
  msglen += encodeUsingLogStarModel(num_chains);

  double constant = log2(4*PI) - 2*log2(AOM);
  for (int i=0; i<spherical_coordinates.size(); i++) {
    // for each chain state the number of residues
    int num_residues = spherical_coordinates[i].size();
    msglen += encodeUsingLogStarModel(num_residues);

    // first point is origin
    // state the second & third points using the sphere model
    vector<double> radii(2,0);
    radii[0] = distances[i][0];
    radii[1] = distances[i][1];
    // state the two points on the surface of sphere
    msglen += 2 * (log2(radii[0]) + log2(radii[1]));
    msglen += 2 * constant;
    // state the radii of the first two points
    msglen += encodeUsingNormalModel(radii);
    radii.clear();
    vector<array<double,2>> points;
    for (int j=0; j<spherical_coordinates[i].size(); j++) {
      double r = spherical_coordinates[i][j][0];
      radii.push_back(r);
      array<double,2> x;
      x[0] = spherical_coordinates[i][j][1];  // theta
      x[1] = spherical_coordinates[i][j][2];  // phi
      points.push_back(x);
    }
    // state the radii
    msglen += encodeUsingNormalModel(radii);
    // state the theta,phi on unit spheres
    msglen += encodeUsingMixtureModel(points,mixture);
  }
  return msglen;
}

/*!
 *  \brief This function is used to initialize code length matrices.
 *  \param chain_index an integer
 */
void Protein::initializeCodeLengthMatrices(int chain_index)
{
  for (int i=0; i<optimal_model.size(); i++) {
    optimal_model[i].clear();
    optimal_code_length[i].clear();
  }
  optimal_model.clear();
  optimal_code_length.clear();
  int n = coordinates[chain_index].size();
  vector<OptimalFit> optimal(n,OptimalFit());
  vector<double> code_length(n,LARGE_NUMBER);
  for (int i=0; i<n; i++) {
    optimal_model.push_back(optimal);
    optimal_code_length.push_back(code_length);
  }
}

/*!
 *  \brief This function computes the code length matrix of individual
 *  pairs in the protein structure.
 *  \param mixture a reference to a Mixture
 */
void Protein::computeCodeLengthMatrix(Mixture &mixture)
{
  vector<IdealModel> ideal_models = loadIdealModels();
  for (int i=0; i<coordinates.size(); i++) {
    initializeCodeLengthMatrices(i);
    int j = 0;
    cout << "cartesian coordinates size: " << coordinates[i].size() << endl;
    cout << "distances size: " << distances[i].size() << endl;
    cout << "spherical coordinates size: " << spherical_coordinates[i].size() << endl;
    while (j < coordinates[i].size()) {
      int range = minimum((int)coordinates[i].size(),j+MAX_SEGMENT_SIZE);
      for (int k=j+MIN_SEGMENT_SIZE-1; k<range; k++) {
        //cout << j << ":" << k << endl;
        Segment segment(j,k,coordinates[i],spherical_coordinates[i]);
        if (j == 0) {
          segment.setInitialDistances(distances[i][0],distances[i][1]);
        }
        OptimalFit fit,ideal_fit;
        // fit null model to the segment
        ideal_fit = segment.fitNullModel(mixture);
        for (int m=0; m<NUM_IDEAL_MODELS; m++) {
          fit = segment.fitIdealModel(ideal_models[m],mixture);
          if (fit < ideal_fit) {
            ideal_fit = fit;
          }
        }
        optimal_model[j][k] = ideal_fit;
        optimal_code_length[j][k] = ideal_fit.getMessageLength();
      }
      if (j == 0) {
        j += MIN_SEGMENT_SIZE - 1;
      } else {
        j++;
      }
    }
    printCodeLengthMatrix(i);
  }
}

/*!
 *  \brief This function outputs the code length matrix to a file.
 *  \param chain_index an integer
 */
void Protein::printCodeLengthMatrix(int chain_index)
{
  string name = "code_length_matrix_chain_";
  name += boost::lexical_cast<string>(chain_index+1);
  ofstream file(name.c_str());
  for (int i=0; i<optimal_code_length.size(); i++) {
    for (int j=0; j<optimal_code_length[i].size(); j++) {
      file << fixed << scientific << optimal_code_length[i][j] << "\t";
    }
    file << endl;
  }
  file.close();
}

