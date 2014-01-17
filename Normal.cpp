#include "Normal.h"

/*!
 *  \brief This is a constructor module
 *  sets default values of mean as 0 and standard deviation as 1
 */
Normal::Normal() : mu(0), sigma(1)
{}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  standard deviation of the distribution
 *  \param mu a double
 *  \param sigma a double
 */
Normal::Normal(double mu, double sigma) : mu(mu), sigma(sigma)
{}

/*!
 *  \brief This function assigns a source Normal distribution.
 *  \param source a reference to a Normal
 */
Normal Normal::operator=(const Normal &source)
{
  if (this != &source) {
    mu = source.mu;
    sigma = source.sigma;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
const double Normal::mean(void)
{
	return mu;
}

/*!
 *  \brief This function returns the standard deviation of the distribution
 *  \return the standard deviation of the distribution
 */
const double Normal::standardDeviation(void)
{
	return sigma;
}

/*!
 *  \brief This function computes the value of the distribution at a given x
 *  \param x a double
 *  \return value of the function given x
 */
double Normal::value(double x)
{
	double expNumerator = (-1) * (x-mu) * (x-mu);
	double expDenominator = (double)2 * sigma * sigma;
	double val = (exp (expNumerator/expDenominator)) / ((sqrt(2*PI)) * sigma);
  return val;
}

/*!
 *  \brief This function computes the negative log likelihood of given data.
 *  \param sample a reference to a vector<double>
 *  \return the negative log likelihood (base e)
 */
double Normal::negativeLogLikelihood(vector<double> &sample)
{
  double value = 0;
  double tmp = 0.5 * log(2*PI) + log(sigma);
  value += (sample.size() * tmp);

  double num = 0;
  for (int i=0; i<sample.size(); i++) {
    num += (sample[i]-mu) * (sample[i]-mu);
  }
  double denom = 2 * sigma * sigma;
  value += num / (double)denom;
  return value;
}

