#include "Poisson.h"

/*!
 *  \brief This is a constructor module
 *  sets default values of mean as 1
 */
Poisson::Poisson() : mu(1)
{}

/*!
 *  \brief constructor function which sets the value of mean  
 *  \param mu a double
 */
Poisson::Poisson(double mu) : mu(mu)
{
  assert(mu > 0);
}

/*!
 *  \brief This function assigns a source Poisson distribution.
 *  \param source a reference to a Poisson
 */
Poisson Poisson::operator=(const Poisson &source)
{
  if (this != &source) {
    mu = source.mu;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
const double Poisson::mean(void)
{
	return mu;
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a double
 *  \return density of the function given x
 */
double Poisson::density(int x)
{
  double tmp1 = exp(-mu);
  double tmp2 = 1;
  for (int i=0; i<x; i++) {
    tmp2 *= mu;
  } 
  long double tmp3 = boost::math::factorial<long double>(x);
  cout << tmp3<<endl;
  return tmp1*tmp2/tmp3;
}

