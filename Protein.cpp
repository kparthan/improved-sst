#include "Protein.h"

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
      coordinates.push_back(chain_coordinates);
    }
  }
  //cout << "# of suitable chains: " << coordinates.size() << endl;
  if (coordinates.size() == 0) {
    cout << name << " is an unsuitable structure ..." << endl;
    ofstream log("errors.log",ios::app);
    log << name << endl;
    log.close();
    exit(1);
  }
  //writeToFile(coordinates,"coordinates");
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
      array<double,3> values = computeSphericalValues(transformed_coordinates);
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
  vector<Point<double>> transformed_coordinates(4,Point<double>());
  // translate index to origin
  for (int i=0; i<4; i++) {
    transformed_coordinates[i] = coordinates[chain_index][index+i-2] -
                                 coordinates[chain_index][index];
  }
  // rotate so that i-1 is on -ve x-axis
  double angle;
  Vector<double> vi_minus_1,negative_xaxis,rotation_axis;
  vector<double> vec1 = {-1,0,0};
  negative_xaxis = vec1;
  vi_minus_1 = transformed_coordinates[1].positionVector();
  rotation_axis = Vector<double>::crossProduct(vi_minus_1,negative_xaxis);
  angle = Vector<double>::angleBetween(vi_minus_1,negative_xaxis);
  for (int i=0; i<4; i++) {
    transformed_coordinates[i] = 
          rotate<double>(transformed_coordinates[i],rotation_axis,angle);
  }
  // rotate so that i-2 is on XY plane
  Vector<double> vi_minus_2,normal,positive_zaxis;
  vector<double> vec2 = {0,0,1};
  positive_zaxis = vec2;
  vi_minus_2 = transformed_coordinates[0].positionVector();
  normal = Vector<double>::crossProduct(vi_minus_2,negative_xaxis);
  angle = Vector<double>::angleBetween(normal,positive_zaxis);
  rotation_axis = negative_xaxis;
  if (vi_minus_2[2] < 0) {
    angle *= -1;
  }
  for (int i=0; i<4; i++) {
    transformed_coordinates[i] =
          rotate<double>(transformed_coordinates[i],rotation_axis,angle);
  }
  return transformed_coordinates;
}

/*!
 *  \brief This function computes the spherical coordinate values
 *  \param transformed_coordinates a reference to a vector<Point<double>>
 *  \return the spherical coordinates at a given index
 */
array<double,3>
Protein::computeSphericalValues(vector<Point<double>> &transformed_coordinates)
{
  array<double,3> values;
  Point<double> origin(0,0,0);
  values[0] = distance<double>(origin,transformed_coordinates[3]);
  Vector<double> i_plus_1 = transformed_coordinates[3].positionVector();
  array<double,2> angles = angleWithAxes<double>(i_plus_1);
  values[1] = angles[0] * (180/PI);
  values[2] = angles[1] * (180/PI);
  return values;
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
  ifstream profile(file_name.c_str());
  string line;
  spherical_coordinates.clear();

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
  if (all_spherical_coordinates.size() != 0) {
    for (int i=0; i<spherical_coordinates.size(); i++) {
      for (int j=0; j<spherical_coordinates[i].size(); j++) {
        all_spherical_coordinates.push_back(spherical_coordinates[i][j]);
      }
    }
  }
  return all_spherical_coordinates;
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

