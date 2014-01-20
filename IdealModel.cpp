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
 *  \param p a reference to a ProteinStructure
 *  \param name a string
 */
IdealModel::IdealModel(ProteinStructure *p, string name): model(p), name(name)
{
  length = p->getNumberOfResidues();
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
    length = source.length;
    name = source.name;
  }
  return *this;
}

