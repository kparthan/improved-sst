#include "IdealModel.h"

/*!
 *  \brief This is a null constructor.
 */
IdealModel::IdealModel()
{}

/*!
 *  \brief This is a constructor function.
 *  \param length an integer
 *  \param name a string
 */
IdealModel::IdealModel(int length, string name): length(length), name(name)
{}

/*!
 *  \brief This is a constructor function.
 *  \param structure a reference to a ProteinStructure
 *  \param name a string
 */
IdealModel::IdealModel(ProteinStructure *structure, string name): 
                       model(structure), name(name)
{
  // get length
  length = structure->getNumberOfResidues();
  // get residues
  vector<string> chain_ids = structure->getChainIdentifiers();
  //assert(chain_ids.size() == 1);
  Chain chain = structure->getDefaultModel()[chain_ids[0]];
  vector<Atom> atoms = chain.getAtoms();
  for (int j=0; j<atoms.size(); j++) {
    Point<double> p = atoms[j].point<double>();
    residues.push_back(p);
  }
  //assert(residues.size() == length);
}

/*!
 *  \brief This function sets the length of the ideal model.
 *  \param new_length an integer
 */
void IdealModel::setLength(int new_length)
{
  length = new_length;
}

/*!
 *  \brief This module assigns an IdealModel object on the rhs to one
 *  on the lhs
 *  \param source a reference to an IdealModel object
 */
IdealModel IdealModel::operator=(const IdealModel &source)
{
  if (this != &source) {
    model = source.model;
    residues = source.residues;
    length = source.length;
    name = source.name;
  }
  return *this;
}

/*!
 *  \brief This function returns a condensed model of a given length
 *  \param num_residues an integer
 *  \return the list of coordinates of the condensed ideal model
 */
vector<Point<double>> IdealModel::getResidues(int num_residues)
{
  vector<Point<double>> points;
  for (int i=0; i<num_residues; i++) {
    points.push_back(residues[i]);
  }
  return points;
}

/*
 *  \brief This function returns the name of the ideal model.
 *  \return the name of the model
 */
string IdealModel::getName()
{
  return name;
}

