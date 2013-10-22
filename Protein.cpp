#include "Protein.h"
//#include "Support.h"

/*!
 *  \brief This is a constructor function used to instantiate the Protein
 *  object from a ProteinStructure
 *  \param protein_structure a reference to a ProteinStructure
 */
Protein::Protein(ProteinStructure *protein_structure) : 
                 protein_structure(protein_structure)
{
  vector<Atom> atoms = protein_structure->getAtoms();
  Point<double> p;
  for (int i=0; i<atoms.size(); i++) {
    p = atoms[i].point<double>();
    coordinates.push_back(p);
  } 
  cout << "Size: " << coordinates.size() << endl;
  //writeToFile(coordinates,"coordinates");
}

/*!
 *  \brief This function transforms the protein to a soherical coordinate
 *  system.
 */
void Protein::getSphericalCoordinatesList()
{
  for (int i=2; i<coordinates.size()-1; i++) {
    vector<Point<double>> transformed_coordinates = computeTransformation(i);
    string file_index = boost::lexical_cast<string>(i);
    //writeToFile(transformed_coordinates,file_index.c_str());
    array<double,3> values = computeSphericalValues(transformed_coordinates);
    spherical_coordinates.push_back(values);
  }
}

/*!
 *  \brief This function computes the spherical coordinate transformation to the
 *  canonical form at a coordinate index
 *  \param index an integer
 *  \return the transformed list of four coordinates
 */
vector<Point<double>> Protein::computeTransformation(int index)
{
  vector<Point<double>> transformed_coordinates(4,Point<double>());
  // translate index to origin
  for (int i=0; i<4; i++) {
    transformed_coordinates[i] = coordinates[index+i-2] - coordinates[index];
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
  values[1] = angles[0];
  values[2] = angles[1];
  return values;
}

