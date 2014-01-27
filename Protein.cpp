#include "Support.h"
//#include "Segment.h"
//#include "Message.h"

extern string CURRENT_DIRECTORY,STRUCTURE;

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
  cout << "# of residues: " << structure->getNumberOfResidues() << endl;
  vector<string> chain_ids = structure->getChainIdentifiers();
  for (int i=0; i<chain_ids.size(); i++) {
    string id = chain_ids[i];
    Chain chain = structure->getDefaultModel()[id];
    vector<Atom> atoms = chain.getAtoms();
    bool chain_break = checkChainBreak(id,atoms);
    if (!chain_break) {
      chains.push_back(id);
      vector<vector<double>> chain_coordinates;
      vector<double> v(3,0);
      for (int j=0; j<atoms.size(); j++) {
        Point<double> p = atoms[j].point<double>();
        point2vector(p,v);
        chain_coordinates.push_back(v);
      }
      //translateProteinToOrigin(chain_coordinates);
      coordinates.push_back(chain_coordinates);
    }
  }
  cout << "# of suitable chains: " << coordinates.size() << endl;
  if (coordinates.size() == 0) {
    cout << name << " is an unsuitable structure ..." << endl;
    ofstream log("unsuitable_structures.log",ios::app);
    log << name << endl;
    log.close();
    exit(1);
  }
}

/*!
 *  \brief This function translates the protein so that its first point
 *  coincides with the origin
 *  \param chain_coordinates a reference to a vector<vector<double>>
 */
void Protein::translateProteinToOrigin(vector<vector<double>> &chain_coordinates)
{
  // before translation
  writeToFile(chain_coordinates,"before_translation");
  initial_translation_vector = vector<double>(3,0);
  for (int i=0; i<3; i++) {
    initial_translation_vector[i] = -chain_coordinates[0][i];
  }
  cout << "Translation vector: ";
  print(cout,initial_translation_vector);

  // translate the protein
  for (int i=0; i<chain_coordinates.size(); i++) {
    for (int j=0; j<3; j++) {
      chain_coordinates[i][j] += initial_translation_vector[j];
    }
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

  // initialize 4-mer and transformarion variables
  vector<vector<double>> four_mer,transformed_four_mer,rotation_matrix;
  initializeMatrix(four_mer,4,3);
  initializeMatrix(transformed_four_mer,4,3);
  initializeMatrix(rotation_matrix,3,3);

  // transform all coordinates
  for (int i=0; i<coordinates.size(); i++) {
    vector<vector<double>> spherical_chain_coordinates;
    vector<double> spherical(3,0);
    for (int j=2; j<coordinates[i].size()-1; j++) {
      computeTransformation(i,j,four_mer,transformed_four_mer,rotation_matrix);
      string file_index = "tmp/" + boost::lexical_cast<string>(j);
      writeToFile(transformed_four_mer,file_index.c_str());
      cartesian2spherical(transformed_four_mer[3],spherical);
      spherical_chain_coordinates.push_back(spherical);
    }
    spherical_coordinates.push_back(spherical_chain_coordinates);
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
void Protein::computeTransformation(int chain_index, int index,
                                    vector<vector<double>> &four_mer, 
                                    vector<vector<double>> &transformed_four_mer,
                                    vector<vector<double>> &rotation_matrix)
{
  // update four-mer
  four_mer[0] = coordinates[chain_index][index-2];
  four_mer[1] = coordinates[chain_index][index-1];
  four_mer[2] = coordinates[chain_index][index];
  four_mer[3] = coordinates[chain_index][index+1];
  // transform
  convertToCanonicalForm(four_mer,transformed_four_mer,rotation_matrix);
}

/*!
 *  \brief Loads the profile from an existing file.
 *  \param identifier a reference to a string
 */
void Protein::load(string &identifier)
{
  name = identifier;
  string file_name = string(CURRENT_DIRECTORY) + "/spherical_system/profiles/"
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
 *  (theta,phi angles read are in radians)
 *  \param file_name a reference to a string
 */
void Protein::read_profile(string &file_name)
{
  cout << "Reading " << file_name << " ..." << endl;
  ifstream profile(file_name.c_str());
  string line;
  spherical_coordinates.clear();
  all_spherical_coordinates.clear();

  vector<vector<double>> chain_coordinates;
  while(getline(profile,line)) {
    boost::char_separator<char> sep(",() ");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    int i = -1;
    vector<double> values(3,0); 
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
 *  (theta,phi angles saved are in radians)
 */
void Protein::save()
{
  cout << "Saving profile of " << STRUCTURE << endl;
  string file_name = string(CURRENT_DIRECTORY) + "/spherical_system/profiles/"
                     + name + ".profile";
  cout << "New profile created: " << file_name << endl;
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
vector<vector<double>> Protein::getSphericalCoordinatesList()
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
vector<double> Protein::computeMeanDirection()
{
  if (all_spherical_coordinates.size() == 0) {
    getSphericalCoordinatesList();
  }
  vector<double> estimate(3,0);
  vector<double> x(3,0);
  
  for (int i=0; i<all_spherical_coordinates.size(); i++) {
    all_spherical_coordinates[i][0] = 1;
    spherical2cartesian(all_spherical_coordinates[i],x);
    for (int j=0; j<3; j++) {
      estimate[j] += x[j];
    }
  }
  return estimate;
}

///*!
// *  \brief This function is used to compute the distance between successive
// *  residues of the protein.
// */
//void Protein::computeSuccessiveDistances()
//{
//  for (int i=0; i<coordinates.size(); i++) {
//    vector<double> dist;
//    for (int j=0; j<coordinates[i].size()-1; j++) {
//      double d = lcb::geometry::distance<double>(coordinates[i][j],coordinates[i][j+1]);
//      dist.push_back(d);
//    }
//    distances.push_back(dist);
//  }
//}
//
///*!
// *  \brief This function computes the message length to communicate the protein
// *  coordinates using the sphere model.
// *  \return the message length
// */
//double Protein::computeMessageLengthUsingSphereModel()
//{
//  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
//  Message message;
//  double msglen = 0;
//
//  // message length to state the number of chains
//  int num_chains = chains.size(); // alternately spherical_coordinates.size()
//  msglen += message.encodeUsingLogStarModel(num_chains);
//
//  for (int i=0; i<distances.size(); i++) {
//    // for each chain state the number of residues
//    int num_residues = distances[i].size();
//    msglen += message.encodeUsingLogStarModel(num_residues);
//
//    for (int j=0; j<distances[i].size(); j++) {
//      msglen += message.encodeUsingSphereModel(distances[i][j],normal);
//    }
//  }  
//  return msglen;
//}
//
///*!
// *  \brief This function computes the message length to communicate the protein
// *  coordinates using the null model.
// *  \param mixture a reference to a Mixture
// *  \return the message length
// */
//double Protein::computeMessageLengthUsingNullModel(Mixture &mixture)
//{
//  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
//  Message message;
//  double msglen = 0;
//
//  // message length to state the number of chains
//  int num_chains = chains.size(); // alternately spherical_coordinates.size()
//  msglen += message.encodeUsingLogStarModel(num_chains);
//
//  for (int i=0; i<spherical_coordinates.size(); i++) {
//    // for each chain state the number of residues
//    int num_residues = spherical_coordinates[i].size();
//    msglen += message.encodeUsingLogStarModel(num_residues);
//
//    // first point is origin
//    // state the second & third points using the sphere model
//    msglen += message.encodeUsingSphereModel(distances[i][0],normal);
//    msglen += message.encodeUsingSphereModel(distances[i][1],normal);
//
//    // state the remaining points using the mixture model
//    double r;
//    array<double,2> x;
//    for (int j=0; j<spherical_coordinates[i].size(); j++) {
//      // state radius
//      r = spherical_coordinates[i][j][0];
//      msglen += message.encodeUsingNormalModel(r,normal);
//      // state theta,phi
//      x[0] = spherical_coordinates[i][j][1];  // theta
//      x[1] = spherical_coordinates[i][j][2];  // phi
//      msglen += message.encodeUsingMixtureModel(x,mixture);
//    }
//  }
//  return msglen;
//}
//
///*!
// *  \brief This function is used to initialize code length matrices.
// *  \param chain_index an integer
// */
//void Protein::initializeCodeLengthMatrices(int chain_index)
//{
//  for (int i=0; i<optimal_model.size(); i++) {
//    optimal_model[i].clear();
//    optimal_code_length[i].clear();
//  }
//  optimal_model.clear();
//  optimal_code_length.clear();
//  int n = coordinates[chain_index].size();
//  vector<OptimalFit> optimal(n,OptimalFit());
//  vector<double> code_length(n,LARGE_NUMBER);
//  for (int i=0; i<n; i++) {
//    optimal_model.push_back(optimal);
//    optimal_code_length.push_back(code_length);
//  }
//}
//
///*!
// *  \brief This function compresses the protein using the expert ideal models.
// *  \param mixture a reference to a Mixture
// *  \param orientation an integer
// */
//void Protein::compressUsingIdealModels(Mixture &mixture, int orientation)
//{
//  vector<IdealModel> ideal_models = loadIdealModels();
//  for (int i=0; i<coordinates.size(); i++) {
//    /*cout << "cartesian coordinates size: " << coordinates[i].size() << endl;
//    cout << "distances size: " << distances[i].size() << endl;
//    cout << "spherical coordinates size: " << spherical_coordinates[i].size() << endl;*/
//    computeCodeLengthMatrix(ideal_models,mixture,orientation,i);
//    pair<double,vector<int>> segmentation = computeOptimalSegmentation(i);
//    vector<int> segments = segmentation.second;
//    cout << "Compression fit: " << segmentation.first << " bits." << endl;
//    cout << "Bits per residue: " << segmentation.first/coordinates[i].size() 
//             << endl << endl; 
//    cout << "# of segments: " << segments.size()-1 << endl << endl;
//    cout << "Internal segmentation:" << endl;
//    int j;
//    for (j=0; j<segments.size()-1; j++) {
//      cout << segments[j] << "-->";
//    }
//    cout << segments[j] << endl << endl;
//  }
//}
//
///*!
// *  \brief This function computes the code length matrix of individual
// *  pairs in the protein structure.
// *  \param ideal_models a reference to a vector<IdealModel> 
// *  \param mixture a reference to a Mixture
// *  \param orientation an integer
// *  \param chain an integer
// */
//void Protein::computeCodeLengthMatrix(vector<IdealModel> &ideal_models,
//                                      Mixture &mixture, int orientation, 
//                                      int chain)
//{
//  initializeCodeLengthMatrices(chain);
//  int chain_size = coordinates[chain].size();
//  for (int i=0; i<chain_size-1; i++) {
//    int bound = minimum(chain_size,i+MAX_SEGMENT_SIZE);
//    for (int j=i+1; j<chain_size; j++) {
//      cout << i << ":" << j << endl;
//      Segment segment(i,j,coordinates[chain],spherical_coordinates[chain]);
//      if (i == 0) {
//        segment.setInitialDistances(distances[chain][0],distances[chain][1]);
//      }
//      int segment_length = j - i + 1; 
//      OptimalFit fit,ideal_fit;
//      // fit null model to the segment
//      ideal_fit = segment.fitNullModel(mixture);
//      if (j < bound) {
//        for (int m=0; m<NUM_IDEAL_MODELS; m++) {
//          if ((m != NUM_IDEAL_MODELS-1 && segment_length >= MIN_SIZE_HELIX) ||
//              (m == NUM_IDEAL_MODELS-1 && segment_length >= MIN_SIZE_STRAND)) {
//            fit = segment.fitIdealModel(ideal_models[m],mixture,orientation);
//            if (fit < ideal_fit) {
//              ideal_fit = fit;
//            }
//          }
//        }
//      }
//      optimal_model[i][j] = ideal_fit;
//      optimal_code_length[i][j] = ideal_fit.getMessageLength();
//    }
//  }
//  printCodeLengthMatrix(chain);
//}
//
///*!
// *  \brief This module computes the optimal segmentation using
// *  dynamic programming
// *  \param chain an integer
// *  \return the indices of the segments
// */
//pair<double,vector<int>> Protein::computeOptimalSegmentation(int chain)
//{
//  pair <double,vector<int>> segmentation;
//  int chain_size = coordinates[chain].size();
//  vector<double> optimal_msglen(chain_size,100000);
//  vector<int> optimal_index(chain_size,-1);
//
//  for (int i=0; i<chain_size; i++){
//    optimal_msglen[i] = optimal_code_length[0][i];
//    optimal_index[i] = i;
//    for (int j=1; j<i; j++){
//      if (optimal_code_length[j][i] + optimal_msglen[j] < optimal_msglen[i]){
//        optimal_msglen[i] = optimal_code_length[j][i] + optimal_msglen[j];
//        optimal_index[i] = j;
//      }
//    }
//  }
//  segmentation.first = optimal_msglen[chain_size-1];
//  int index = chain_size - 1;
//  vector<int> backtrack; 
//  backtrack.push_back(chain_size-1);
//  while (1){
//    if (index == optimal_index[index]){
//      break;
//    }
//    index = optimal_index[index];
//    backtrack.push_back(index);
//  }
//  backtrack.push_back(0);
//  vector<int> segments;
//  for (int i=backtrack.size()-1; i>=0; i--){
//    segments.push_back(backtrack[i]);
//  }
//  segmentation.second = segments;
//  return segmentation;
//}
///*!
// *  \brief This function outputs the code length matrix to a file.
// *  \param chain_index an integer
// */
//void Protein::printCodeLengthMatrix(int chain_index)
//{
//  string name = "code_length_matrix_chain_";
//  name += boost::lexical_cast<string>(chain_index+1);
//  ofstream file(name.c_str());
//  for (int i=0; i<optimal_code_length.size(); i++) {
//    for (int j=0; j<optimal_code_length[i].size(); j++) {
//      file << fixed << scientific << optimal_code_length[i][j] << "\t";
//    }
//    file << endl;
//  }
//  file.close();
//}
//
