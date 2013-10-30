#include "VonMises3D.h"
#include "Support.h"

/*!
 *  \brief This is a null constructor module
 *  sets default values of mean vector as 0 and kappa as 1
 */
VonMises3D::VonMises3D() : mu(array<double,2>()),kappa(1)
{
  computeConstants();
}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  standard deviation of the distribution
 *  \param mean a reference to a array<double,2>
 *  \param kappa a double
 */
VonMises3D::VonMises3D(array<double,2> &mu, double kappa) : 
                       mu(mu), kappa(kappa)
{
  computeConstants();
}

/*!
 *  \brief This function pre computes the constants.
 */
void VonMises3D::computeConstants()
{
  c3k = kappa / (4 * PI * sinh(k));
  array<double,3> x = convertToCartesian(1,mu[0],mu[1]);
  for (int i=0; i<3; i++) {
    kmu[i] = kappa * x[i];
  }
}

/*!
 *  \brief This function assigns a source VonMises3D distribution.
 *  \param source a reference to a VonMises3D
 */
VonMises3D VonMises3D::operator=(const VonMises3D &source)
{
  if (this != &source) {
    mu = source.mu;
    kappa = source.kappa;
    c3k = source.c3k;
    kmu = source.kmu;
  }
  return *this;
}

/*!
 *  \brief This function returns the mu of the distribution
 *  \return the mu of the distribution
 */
const array<double,2> VonMises3D::mean(void)
{
	return mu;
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
 *  \param theta a double
 *  \param phi a double
 *  \return the probability density
 */
double VonMises3D::density(double theta, double phi)
{
  array<double,3> x = convertToCartesian(1,theta,phi);
  double exponent = 0;
  for (int i=0; i<3; i++) {
    exponent += kmu[i] * x[i];
  }
  return c3k * exp(exponent);
}

