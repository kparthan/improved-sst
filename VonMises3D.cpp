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
  double tmp = exp(-2 * kappa);
  double denom = 2 * PI * (1-tmp);
  constant = kappa / denom;
  unit_mean = convertToCartesian(1,mu[0],mu[1]);
}

/*!
 *  \brief This function assigns a source VonMises3D distribution.
 *  \param source a reference to a VonMises3D
 */
VonMises3D VonMises3D::operator=(const VonMises3D &source)
{
  if (this != &source) {
    mu = source.mu;
    unit_mean = source.unit_mean;
    kappa = source.kappa;
    constant = source.constant;
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
  if (kappa <= TOLERANCE) {
    return 1.0/(4*PI);
  } else {
    array<double,3> x = convertToCartesian(1,theta,phi);
    double tmp = 0;
    for (int i=0; i<3; i++) {
      tmp += unit_mean[i] * x[i];
    }
    double exponent = kappa * (tmp - 1);
    return constant * exp(exponent);
  }
}

