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
 *  \brief This module assigns an IdealModel object on the rhs to one
 *  on the lhs
 *  \param source a reference to an IdealModel object
 */
IdealModel IdealModel::operator=(const IdealModel &source)
{
  if (this != &source) {
    length = source.length;
    name = source.name;
  }
  return *this;
}

