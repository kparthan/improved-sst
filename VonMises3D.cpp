#include "VonMises3D.h"

/*!
 *  \brief This is a null constructor module
 *  sets default values of mean vector as 0 and kappa as 1
 */
VonMises3D::VonMises3D() : mean(array<double,2>()),kappa(1)
{}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  standard deviation of the distribution
 *  \param mean a reference to a array<double,2>
 *  \param kappa a double
 */
VonMises3D::VonMises3D(array<double,2> &mean, double kappa) : 
                       mean(mean), kappa(kappa)
{}

/*!
 *  \brief This function assigns a source VonMises3D distribution.
 *  \param source a reference to a VonMises3D
 */
VonMises3D VonMises3D::operator=(const VonMises3D &source)
{
  if (this != &source) {
    mean = source.mean;
    kappa = source.kappa;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
const array<double,2> VonMises3D::mean(void)
{
	return mean;
}

/*!
 *  \brief This function returns the standard deviation of the distribution
 *  \return the scale parameter of the distribution
 */
const double VonMises3D::scale(void)
{
	return kappa;
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a reference to a double
 *  \return value of the function given x
 */
double VonMises3D::value(array<double,2> &x)
{
}

