#include "Support.h"
#include "Message.h"
#include "Geometry3D.h"
#include "Segmentation.h"

#define MIN_SIZE_HELIX_ALPHA 6
#define MIN_SIZE_HELIX_PI 4 
#define MIN_SIZE_HELIX_310 4
#define MIN_SIZE_STRAND 5 
#define MAX_SEGMENT_SIZE 40 

extern string CURRENT_DIRECTORY,STRUCTURE;
int PORTION_TO_FIT;

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
  //cout << "# of residues: " << structure->getNumberOfResidues() << endl;
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
      translateProteinToOrigin(chain_coordinates);
      cartesian_coordinates.push_back(chain_coordinates);
    }
  }
  //cout << "# of suitable chains: " << cartesian_coordinates.size() << endl;
  if (cartesian_coordinates.size() == 0) {
    cout << name << " is an unsuitable structure ..." << endl;
    ofstream log("unsuitable_structures.log",ios::app);
    log << name << endl;
    log.close();
    exit(1);
  }
}

/*!
 *  \brief This function returns the suitable chain ids used.
 *  \return the list of chain ids
 */
vector<string> Protein::getChainIds()
{
  return chains;
}

/*!
 *  \brief This function translates the protein so that its first point
 *  coincides with the origin
 *  \param chain_coordinates a reference to a vector<vector<double>>
 */
void Protein::translateProteinToOrigin(vector<vector<double>> &chain_coordinates)
{
  // before translation
  //writeToFile(chain_coordinates,"before_translation");
  initial_translation_vector = vector<double>(3,0);
  for (int i=0; i<3; i++) {
    initial_translation_vector[i] = -chain_coordinates[0][i];
  }
  //cout << "Translation vector: ";
  //print(cout,initial_translation_vector);

  // translate the protein
  for (int i=0; i<chain_coordinates.size(); i++) {
    for (int j=0; j<3; j++) {
      chain_coordinates[i][j] += initial_translation_vector[j];
    }
  }
  // after translation
  //writeToFile(chain_coordinates,"after_translation");
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
  for (int i=0; i<cartesian_coordinates.size(); i++) {
    vector<vector<double>> chain_coordinates;
    vector<double> spherical(3,0);
    for (int j=2; j<cartesian_coordinates[i].size()-1; j++) {
      computeTransformation(i,j,four_mer,transformed_four_mer,rotation_matrix);
      string file_index = "tmp/" + boost::lexical_cast<string>(j);
      writeToFile(transformed_four_mer,file_index.c_str());
      cartesian2spherical(transformed_four_mer[3],spherical);
      chain_coordinates.push_back(spherical);
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
void Protein::computeTransformation(int chain_index, int index,
                                    vector<vector<double>> &four_mer, 
                                    vector<vector<double>> &transformed_four_mer,
                                    vector<vector<double>> &rotation_matrix)
{
  // update four-mer
  four_mer[0] = cartesian_coordinates[chain_index][index-2];
  four_mer[1] = cartesian_coordinates[chain_index][index-1];
  four_mer[2] = cartesian_coordinates[chain_index][index];
  four_mer[3] = cartesian_coordinates[chain_index][index+1];
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
 *  Each line in the file should contain the chain index followed by
 *  radius (\neq 1), theta, phi (measured in radians) values.
 *  \param file_name a reference to a string
 */
void Protein::read_profile(string &file_name)
{
  cout << "Reading " << file_name << " ..." << endl;
  ifstream profile(file_name.c_str());
  string line;
  cartesian_coordinates.clear();
  spherical_coordinates.clear();
  unit_coordinates.clear();
  all_unit_coordinates.clear();

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
 *  Each line in the file contains the chain index followed by
 *  radius (\neq 1), theta, phi (measured in radians) values.
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
 *  \brief This function gets the list of all coordinates on the unit sphere.
 *  \return the list of all unit Cartesian coordinates.
 */
vector<vector<double>> Protein::getUnitCoordinatesList()
{
  if (all_unit_coordinates.size() == 0) {
    vector<double> x(3,1);
    vector<double> cartesian(3,0);
    for (int i=0; i<spherical_coordinates.size(); i++) {
      vector<vector<double>> chain_coordinates;
      for (int j=0; j<spherical_coordinates[i].size(); j++) {
        x[1] = spherical_coordinates[i][j][1];
        x[2] = spherical_coordinates[i][j][2];
        spherical2cartesian(x,cartesian);
        chain_coordinates.push_back(cartesian);
        all_unit_coordinates.push_back(cartesian);
      }
      unit_coordinates.push_back(chain_coordinates);
    }
  }
  return all_unit_coordinates;
}

/*!
 *  \brief This function gets the list of all spherical coordinates
 */
vector<vector<vector<double>>> Protein::getSphericalCoordinatesList()
{
  return spherical_coordinates;
}

/*!
 *  \brief This function returns the size of transformed spherical coordinates.
 *  \return the size
 */
int Protein::getNumberOfSphericalCoordinates()
{
  return all_unit_coordinates.size();
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
  if (all_unit_coordinates.size() == 0) {
    getUnitCoordinatesList();
  }
  vector<double> estimate(3,0);
  
  for (int i=0; i<all_unit_coordinates.size(); i++) {
    for (int j=0; j<3; j++) {
      estimate[j] += all_unit_coordinates[i][j];
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
  for (int i=0; i<cartesian_coordinates.size(); i++) {
    vector<double> dist;
    for (int j=0; j<cartesian_coordinates[i].size()-1; j++) {
      double d = computeEuclideanDistance(cartesian_coordinates[i][j],cartesian_coordinates[i][j+1]);
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
  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
  Message message;
  double msglen = 0;

  // message length to state the number of chains
  int num_chains = chains.size(); // alternately spherical_coordinates.size()
  //msglen += message.encodeUsingLogStarModel(num_chains);

  for (int i=0; i<spherical_coordinates.size(); i++) {
    // for each chain state the number of residues
    //int num_residues = distances[i].size();
    //msglen += message.encodeUsingLogStarModel(num_residues);

    for (int j=0; j<spherical_coordinates[i].size(); j++) {
      double radius = spherical_coordinates[i][j][0];
      double theta = spherical_coordinates[i][j][1];
      msglen += message.encodeUsingSphereModel(radius,normal);
      //msglen += message.encodeUsingSphereModel(radius,theta,normal);
    }
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
  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
  Message message;
  double msglen = 0;

  // message length to state the number of chains
  int num_chains = chains.size(); // alternately spherical_coordinates.size()
  //msglen += message.encodeUsingLogStarModel(num_chains);

  for (int i=0; i<spherical_coordinates.size(); i++) {
    // for each chain state the number of residues
    //int num_residues = spherical_coordinates[i].size();
    //msglen += message.encodeUsingLogStarModel(num_residues);

    // first point is origin
    // state the second & third points using the sphere model
    //msglen += message.encodeUsingSphereModel(distances[i][0],normal);
    //msglen += message.encodeUsingSphereModel(distances[i][1],normal);

    // state the remaining points using the mixture model
    double r;
    for (int j=0; j<spherical_coordinates[i].size(); j++) {
      // state radius
      r = spherical_coordinates[i][j][0];
      //msglen += message.encodeUsingNormalModel(r,normal);
      // state direction 
      msglen += message.encodeUsingMixtureModel(unit_coordinates[i][j],mixture,r);
    }
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
  int n = cartesian_coordinates[chain_index].size();
  vector<OptimalFit> optimal(n,OptimalFit());
  vector<double> code_length(n,LARGE_NUMBER);
  for (int i=0; i<n; i++) {
    optimal_model.push_back(optimal);
    optimal_code_length.push_back(code_length);
  }
}

/*!
 *  \brief This function compresses the protein using the expert ideal models.
 *  \param mixture a reference to a Mixture
 *  \param orientation an integer
 *  \param portion_to_fit an integer
 *  \param end_points a reference to a vector<string>
 *  \param sst_method an integer
 */
void Protein::compressUsingIdealModels(Mixture &mixture, int orientation, 
                                       int portion_to_fit, vector<string> &end_points,
                                       int sst_method)
{
  PORTION_TO_FIT = portion_to_fit;
  vector<IdealModel> ideal_models = loadIdealModels();
  ofstream log("compression.log");
  if (portion_to_fit == FIT_ENTIRE_STRUCTURE) {
    for (int i=0; i<cartesian_coordinates.size(); i++) {
      /*cout << "cartesian coordinates size: " << cartesian_coordinates[i].size() << endl;
      cout << "distances size: " << distances[i].size() << endl;
      cout << "spherical coordinates size: " << spherical_coordinates[i].size() << endl;*/
      computeCodeLengthMatrix(ideal_models,mixture,orientation,i,sst_method,log);
      pair<double,vector<int>> segmentation = computeOptimalSegmentation(i);
      vector<int> segments = segmentation.second;
      cout << "\nCompression fit: " << segmentation.first << " bits. ("
           << segmentation.first/cartesian_coordinates[i].size() << " bpr)\n\n";
      cout << "# of segments: " << segments.size()-1 << endl << endl;
      cout << "Internal segmentation:" << endl;
      int j;
      for (j=0; j<segments.size()-1; j++) {
        cout << segments[j]+1 << "-->";
      }
      cout << segments[j]+1 << endl << endl;
      int a,b;
      vector<string> model_names;
      string name;
      for (j=0; j<segments.size()-2; j++) {
        a = segments[j];
        b = segments[j+1];
        name = optimal_model[a][b].getName();
        cout << name << "-->";
        model_names.push_back(name);
      }
      a = segments[j]; b = segments[j+1];
      name = optimal_model[a][b].getName();
      cout << name << endl;
      model_names.push_back(name);
      Segmentation assignment(structure,chains[i],segments,model_names);
      assignment.postprocess();
    }
  } else if (portion_to_fit == FIT_SINGLE_SEGMENT) {
    int start = boost::lexical_cast<int>(end_points[1]) - 1;
    int end = boost::lexical_cast<int>(end_points[2]) - 1;
    int chain = -1;
    for (int i=0; i<chains.size(); i++) {
      if (chains[i].compare(end_points[0]) == 0) {
        chain = i;
        cout << "chain_index: " << chain << endl;
      }
    }
    if (chain == -1) {
      cout << "Chain " << end_points[0] << " doesn't exist!\n";
      exit(1);
    }
    Segment segment(start,end,cartesian_coordinates[chain],
                    spherical_coordinates[chain],unit_coordinates[chain]);
    log << "Segment: chain " << chains[chain] << "\t[" << start+1 << "," << end+1 << "]\n";
    if (start == 0 || start == 1) {
      segment.setInitialDistances(distances[chain][0],distances[chain][1]);
    }
    fitOneSegment(ideal_models,segment,mixture,sst_method,orientation,log);
  }
  log.close();
}

/*!
 *  \brief This function fits the ideal models to a single segment of the
 *  protein structure using the adaptive superposition method.
 *  \param ideal_models a reference to a vector<IdealModel>
 *  \param segment a reference to a Segment 
 *  \param mixture a reference to a Mixture
 *  \param orientation an integer
 *  \param log a reference to a ostream
 *  \return the ideal fit for the segment
 */
OptimalFit Protein::fitSegment ( // adaptive superposition
  vector<IdealModel> &ideal_models,
  Segment &segment,
  Mixture &mixture, 
  int orientation, 
  ostream &log
) {
  int segment_length = segment.length(); 
  vector<OptimalFit> fit;
  OptimalFit ideal_fit,current_fit;

  ideal_fit = segment.fitNullModel(mixture,log);
  fit.push_back(ideal_fit);
  for (int m=0; m<NUM_IDEAL_MODELS; m++) {
    if ((m != NUM_IDEAL_MODELS-1 && segment_length >= MIN_SIZE_HELIX_ALPHA) ||
        (m == NUM_IDEAL_MODELS-1 && segment_length >= MIN_SIZE_STRAND)) {
      current_fit = segment.fitIdealModel(ideal_models[m],mixture,orientation,log);
      fit.push_back(current_fit);
      if (current_fit < ideal_fit) {
        ideal_fit = current_fit;
      }
    }
  }
  if (PORTION_TO_FIT == FIT_SINGLE_SEGMENT) {
    cout << "\n\nPrinting fit info:\n";
    for (int i=0; i<fit.size(); i++) {
      fit[i].printFitInfo();
      cout << endl;
    }
    cout << "\nBest fit: ";
    ideal_fit.printFitInfo();
  }
  return ideal_fit;
}

/*!
 *  \brief This function is used to fit a single segment in a non-adaptive encoding
 *  scheme using the residual mixture.
 *  \param ideal_models a reference to a vector<IdealModel>
 *  \param segment a reference to a Segment 
 *  \param residual_mixture a reference to a Mixture
 *  \param components a reference to a vector<Component>
 *  \param weights a reference to a vector<double>
 *  \param assignment a reference to a vector<int>
 *  \param log a reference to a ostream
 *  \return the ideal fit for the segment
 */
OptimalFit Protein::fitSegment (  // non-adaptive method
  vector<IdealModel> &ideal_models,
  Segment &segment,
  Mixture &residual_mixture,
  double &sum_residual_weights,
  vector<Component> &components,
  vector<double> &weights,
  vector<int> &assignment,
  ostream &log
) {
  int segment_length = segment.length(); 
  vector<OptimalFit> fit;
  OptimalFit current_fit,ideal_fit;

  ideal_fit = segment.fitNullModel(residual_mixture,sum_residual_weights,log);
  fit.push_back(ideal_fit);
  for (int m=0; m<NUM_IDEAL_MODELS; m++) {
    if ((m != NUM_IDEAL_MODELS-1 && segment_length >= MIN_SIZE_HELIX_ALPHA) ||
        (m == NUM_IDEAL_MODELS-1 && segment_length >= MIN_SIZE_STRAND)) {
      Component assigned_component = components[assignment[m]];
      double weight = weights[assignment[m]];
      current_fit = segment.fitIdealModel(ideal_models[m],residual_mixture,
                                          assigned_component,weight,log);
      fit.push_back(current_fit);
      if (current_fit < ideal_fit) {
        ideal_fit = current_fit;
      }
    }
  }
  if (PORTION_TO_FIT == FIT_SINGLE_SEGMENT) {
    cout << "\n\nPrinting fit info:\n";
    for (int i=0; i<fit.size(); i++) {
      fit[i].printFitInfo();
      cout << endl;
    }
    cout << "\nBest fit: ";
    ideal_fit.printFitInfo();
  }
  return ideal_fit;
}

/*!
 *  \brief This function is used to fit a single segment using the ideal mixture 
 *  models generated using DSSP
 *  \param ideal_mixture_models a reference to a vector<Mixture>
 *  \param model_weights a reference to a vector<double>
 *  \param segment a reference to a Segment 
 *  \param log a reference to a ostream
 *  \return the ideal fit for the segment
 */
OptimalFit Protein::fitSegment (  // dssp non-adaptive
  vector<Mixture> &ideal_mixture_models,
  vector<double> &model_weights,
  Segment &segment,
  ostream &log
) {
  int segment_length = segment.length(); 
  vector<OptimalFit> fit;
  OptimalFit current_fit,ideal_fit;
  vector<double> min_sizes(ideal_mixture_models.size()-1,0);
  for (int i=0; i<min_sizes.size(); i++) {
    if (i == 0) {
      min_sizes[i] = MIN_SIZE_HELIX_ALPHA;
    } else if (i == 1) {
      min_sizes[i] = MIN_SIZE_HELIX_PI;
    } else if (i == 2) {
      min_sizes[i] = MIN_SIZE_HELIX_310;
    } else if (i == 3) {
      min_sizes[i] = MIN_SIZE_STRAND;
    }
  }

  ideal_fit = segment.fitIdealModel(ideal_mixture_models[4],ideal_mixture_models[4],model_weights[4],log);
  fit.push_back(ideal_fit);
  for (int m=0; m<4; m++) {
    if (segment_length >= min_sizes[m]) { 
      current_fit = segment.fitIdealModel(ideal_mixture_models[4],ideal_mixture_models[m],model_weights[m],log);
      fit.push_back(current_fit);
      if (current_fit < ideal_fit) {
        ideal_fit = current_fit;
      }
    }
  }
  if (PORTION_TO_FIT == FIT_SINGLE_SEGMENT) {
    cout << "\n\nPrinting fit info:\n";
    for (int i=0; i<fit.size(); i++) {
      fit[i].printFitInfo();
      cout << endl;
    }
    cout << "\nBest fit: ";
    ideal_fit.printFitInfo();
  }
  return ideal_fit;
}

/*!
 *  \brief This function fits the ideal models to a JUST one segment of the
 *  protein structure using the adaptive superposition method.
 *  \param ideal_models a reference to a vector<IdealModel>
 *  \param segment a reference to a Segment 
 *  \param mixture a reference to a Mixture
 *  \param sst_method an integer
 *  \param orientation an integer
 *  \param log a reference to a ostream
 */
void Protein::fitOneSegment(
  vector<IdealModel> &ideal_models,
  Segment &segment,
  Mixture &mixture, 
  int sst_method,
  int orientation, 
  ostream &log
) {
  OptimalFit ideal_fit;
  switch(sst_method) {
    case ONE_COMPONENT_ADAPTIVE:
      break;

    case MIXTURE_ADAPTIVE:
      ideal_fit = fitSegment(ideal_models,segment,mixture,orientation,log);
      break;

    case NON_ADAPTIVE:
    {
      vector<Component> components = mixture.getComponents();
      vector<double> weights = mixture.getWeights();
      vector<int> assignment = assignMixtureComponents(ideal_models,components,weights);
      double sum_residual_weights = 1;
      Mixture residual_mixture = mixture.getResidualMixture(assignment,sum_residual_weights);
      ideal_fit = fitSegment(ideal_models,segment,residual_mixture,
                  sum_residual_weights,components,weights,assignment,log);
      break;
    }

    case DSSP_NON_ADAPTIVE:
    {
      vector<Mixture> ideal_mixture_models = loadIdealMixtureModels();
      vector<double> model_weights = computeRelativeWeights(ideal_mixture_models);
      ideal_fit = fitSegment(ideal_mixture_models,model_weights,segment,log);
      break;
    }
  }
}

/*!
 *  \brief This function computes the code length matrix of individual
 *  pairs in the protein structure.
 *  \param ideal_models a reference to a vector<IdealModel> 
 *  \param mixture a reference to a Mixture
 *  \param orientation an integer
 *  \param chain an integer
 *  \param sst_method an integer
 *  \param log a reference to a ostream
 */
void Protein::computeCodeLengthMatrix(vector<IdealModel> &ideal_models,
                                      Mixture &mixture, int orientation, 
                                      int chain, int sst_method, ostream &log)
{
  initializeCodeLengthMatrices(chain);
  int chain_size = cartesian_coordinates[chain].size();
  int i,j,m,bound,segment_length;
  OptimalFit ideal_fit;

  switch(sst_method) {
    case ONE_COMPONENT_ADAPTIVE:
      break;

    case MIXTURE_ADAPTIVE:
    {
      for (i=0; i<chain_size-1; i++) {
        bound = minimum(chain_size,i+MAX_SEGMENT_SIZE);
        for (j=i+1; j<chain_size; j++) {
          log << "Segment " << i+1 << ":" << j+1 << endl;
          cout << "Segment " << i+1 << ":" << j+1 << endl;
          Segment segment(i,j,cartesian_coordinates[chain],
                          spherical_coordinates[chain],unit_coordinates[chain]);
          if (i == 0 || i == 1) {
            segment.setInitialDistances(distances[chain][0],distances[chain][1]);
          }
          if (j >= bound) {
            ideal_fit = segment.fitNullModel(mixture,log);
          } else {
            ideal_fit = fitSegment(ideal_models,segment,mixture,orientation,log);
          }
          optimal_model[i][j] = ideal_fit;
          optimal_code_length[i][j] = ideal_fit.getMessageLength();
        }
      }
      break;
    }

    case NON_ADAPTIVE:
    {
      vector<Component> components = mixture.getComponents();
      vector<double> weights = mixture.getWeights();
      vector<int> assignment = assignMixtureComponents(ideal_models,components,weights);
      double sum_residual_weights = 1;
      Mixture residual_mixture = mixture.getResidualMixture(assignment,sum_residual_weights);

      for (i=0; i<chain_size-1; i++) {
        bound = minimum(chain_size,i+MAX_SEGMENT_SIZE);
        for (j=i+1; j<chain_size; j++) {
          log << "Segment " << i+1 << ":" << j+1 << endl;
          cout << "Segment " << i+1 << ":" << j+1 << endl;
          Segment segment(i,j,cartesian_coordinates[chain],
                          spherical_coordinates[chain],unit_coordinates[chain]);
          if (i == 0 || i == 1) {
            segment.setInitialDistances(distances[chain][0],distances[chain][1]);
          }
          if (j >= bound) {
            ideal_fit = segment.fitNullModel(residual_mixture,sum_residual_weights,log);
          } else {
            ideal_fit = fitSegment(ideal_models,segment,residual_mixture,
                        sum_residual_weights,components,weights,assignment,log);
          }
          optimal_model[i][j] = ideal_fit;
          optimal_code_length[i][j] = ideal_fit.getMessageLength();
        }
      }
      break;
    }

    case DSSP_NON_ADAPTIVE:
    {
      vector<Mixture> ideal_mixture_models = loadIdealMixtureModels();
      vector<double> model_weights = computeRelativeWeights(ideal_mixture_models);
      for (i=0; i<chain_size-1; i++) {
        bound = minimum(chain_size,i+MAX_SEGMENT_SIZE);
        for (j=i+1; j<chain_size; j++) {
          log << "Segment " << i+1 << ":" << j+1 << endl;
          cout << "Segment " << i+1 << ":" << j+1 << endl;
          Segment segment(i,j,cartesian_coordinates[chain],
                          spherical_coordinates[chain],unit_coordinates[chain]);
          if (i == 0 || i == 1) {
            segment.setInitialDistances(distances[chain][0],distances[chain][1]);
          }
          if (j >= bound) {
            ideal_fit = segment.fitIdealModel(ideal_mixture_models[4],ideal_mixture_models[4],model_weights[4],log);
          } else {
            ideal_fit = fitSegment(ideal_mixture_models,model_weights,segment,log);
          }
          optimal_model[i][j] = ideal_fit;
          optimal_code_length[i][j] = ideal_fit.getMessageLength();
        }
      }
      break;
    } 
  }
  printCodeLengthMatrix(chain);
}

/*!
 *  \brief This module computes the optimal segmentation using
 *  dynamic programming
 *  \param chain an integer
 *  \return the indices of the segments
 */
pair<double,vector<int>> Protein::computeOptimalSegmentation(int chain)
{
  pair <double,vector<int>> segmentation;
  int chain_size = cartesian_coordinates[chain].size();
  vector<double> optimal_msglen(chain_size,LARGE_NUMBER);
  vector<int> optimal_index(chain_size,-1);

  for (int i=0; i<chain_size; i++){
    optimal_msglen[i] = optimal_code_length[0][i];
    optimal_index[i] = i;
    for (int j=1; j<i; j++){
      if (optimal_code_length[j][i] + optimal_msglen[j] < optimal_msglen[i]){
        optimal_msglen[i] = optimal_code_length[j][i] + optimal_msglen[j];
        optimal_index[i] = j;
      }
    }
  }
  segmentation.first = optimal_msglen[chain_size-1];
  //cout << "Net optimal msglen: " << segmentation.first << " bits. ("
  //     << segmentation.first/chain_size << " bpr)\n";
  int index = chain_size - 1;
  vector<int> backtrack; 
  backtrack.push_back(chain_size-1);
  while(1) {
    if (index == optimal_index[index]){
      break;
    }
    index = optimal_index[index];
    backtrack.push_back(index);
  }
  backtrack.push_back(0);
  vector<int> segments;
  for (int i=backtrack.size()-1; i>=0; i--){
    segments.push_back(backtrack[i]);
  }
  segmentation.second = segments;
  return segmentation;
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
      file << "Segment " << i+1 << ":" << j+1 << "\t";
      //file << fixed << scientific << optimal_code_length[i][j] << "\t";
      file << fixed << setw(10) << setprecision(4) << optimal_code_length[i][j] << "\n";
    }
    file << endl;
  }
  file.close();
}

