#include "Component.h"
#include "Support.h"

/*!
 *  \brief The null constructor module.
 */
Component::Component()
{}

/*!
 *  \brief This is a constructor function.
 *  \param mean_direction a reference to an array<double,3>
 *  \param N a double
 */
Component::Component(array<double,3> &mean_direction, double N) :
                     mean_direction(mean_direction), N(N)
{}

/*!
 *  \brief This is a constructor function.
 *  \param mean_direction a reference to an array<double,2>
 *  \param kappa a double
 */
Component::Component(array<double,2> &mean_direction, double kappa) :
                     mu(mean_direction), kappa_mml(kappa)
{}

/*!
 *  \brief This function is used to minimize the message length expression
 *  by finding the suitable model parameters.
 */
void Component::minimizeMessageLength()
{
  estimateVonMisesMean();
  kappa_ml = estimateKappa_ML();
  kappa_mml = estimateKappa_MML(kappa_ml);
  von_mises = VonMises3D(mu,kappa_mml);
}

/*!
 *  \brief This function is used to estimate the Von Mises mean
 */
void Component::estimateVonMisesMean()
{
  cout << "\nCartesian coordinates of mean direction vector: ";
  print(cout,mean_direction);
  Point<double> point(mean_direction);
  array<double,3> spherical_coordinates = convertToSpherical(point);
  cout << "Spherical coordinates of mean direction vector: ";
  print(cout,spherical_coordinates);

  // magnitude of the mean direction vector
  R = spherical_coordinates[0];

  // normalize the mean direction vector to obtain an unit vector
  for (int i=0; i<3; i++) {
    mean_direction[i] /= R;
  }
  cout << "Cartesian coordinates of unit mean direction vector: ";
  print(cout,mean_direction);
  point = Point<double>(mean_direction);
  spherical_coordinates = convertToSpherical(point);
  cout << "Spherical coordinates of unit mean direction vector: ";
  print(cout,spherical_coordinates);

  // Estimates of theta, phi (in degrees)
  mu[0] = spherical_coordinates[1];
  mu[1] = spherical_coordinates[2];

  rbar = R / (double)N ;
  cout << "R of direction vector: " << R << endl;
  cout << "N: " << N << endl;
  cout << "rbar: " << rbar << endl;

  // print the estimates of theta, phi
  cout << "\nEstimates:\n";
  cout << "(theta,phi) = (" << spherical_coordinates[1] << ", " 
       << spherical_coordinates[2] << ")" << endl;
}

/*!
 *  \brief This function computes the ML estimate (approximate) of kappa.
 *  (ref: http://en.wikipedia.org/wiki/Von_Mises%E2%80%93Fisher_distribution)
 *  \return the ML estimate of kappa
 */
double Component::estimateKappa_ML()
{
  double kappa = (rbar * (3 - (rbar * rbar))) / (1 - (rbar * rbar));
  cout << "Kappa (ML): " << kappa << endl;
  return kappa;
}

/*!
 *  \brief This function is used to compute the MML estimate of kappa.
 *  \param initial a double
 *  \return the ML estimate of kappa
 */
double Component::estimateKappa_MML(double initial)
{
  double prev = initial;
  double current;
  int num_iterations = 0;
  while(1) {
    num_iterations++;
    if (prev < 0) {
      prev = fabs(prev);
    }
    double fx = computeFirstDerivative(prev);
    double fx_der = computeSecondDerivative(prev);
    if (fabs(fx_der) > TOLERANCE) {
      current = prev - (fx/(double)fx_der);
      cout << "Iteration " << num_iterations << ": [" << prev << ", " 
           << current <<  ", " << fx << ", " << fx_der << "]" << endl;
      if (fabs(current - prev) > TOLERANCE) {
        prev = current;
      } else {
        cout << "No significant change ..." << endl;
        cout << "current: " << current << endl;
        return current;
      }
    } else {
      cout << "Derivative is zero ..." << endl;
      return prev;
    }
  }
}

/*!
 *  \brief This function computes the minimum message length value for an
 *  assumed kappa.
 *  \param kappa a double
 *  \return the message length for a given kappa value
 */
double Component::computeComponentLength(double kappa)
{
  double part1 = computeFirstPart(kappa);
  double part2 = computeSecondPart(kappa);
  double msglen = part1 + part2;
  return msglen;
}

/*!
 *  \brief This function computes the first part of the message
 *  \param kappa a double
 *  \return the the first message length for a given kappa value
 */
double Component::computeFirstPart(double kappa)
{
  double lattice_constant = 0.08;
  double part1 = 1.5 * (log(lattice_constant) + log(N)) + 2 * LOG_PI
                 - log(kappa) + 2 * log(1 + kappa * kappa) 
                 + log(ratioBesselFunction(kappa)) 
                 + 0.5 * log(ratioBesselFunction_firstDerivative(kappa));
  return part1;
}

/*!
 *  \brief This function computes the second part of the message
 *  \param kappa a double
 *  \return the the second message length for a given kappa value
 */
double Component::computeSecondPart(double kappa)
{
  double part2 = 1.5 - kappa * R -  N * log(kappa) + N * log(4) + N * LOG_PI
                 + N * log(sinh(kappa));
  return part2;
}

/*!
 *  \brief This function computes the first derivative of the message length
 *  expression at a given kappa.
 *  \param kappa a double
 *  \return the derivative value
 */
double Component::computeFirstDerivative(double kappa)
{
  double A = ratioBesselFunction(kappa);
  double A1 = ratioBesselFunction_firstDerivative(kappa);
  double A2 = ratioBesselFunction_secondDerivative(kappa);

  double derivative_msglen = -R + N * A - (1/(double)kappa);
  derivative_msglen += A1 / (double) A;
  derivative_msglen += (0.5 * A2) / (double) A1;
  derivative_msglen += (4 * kappa) / (double) (1 + kappa * kappa);

  return derivative_msglen;
}

/*!
 *  \brief This function computes the second derivative of the message length
 *  expression at a given kappa.
 *  \param kappa a double
 *  \return the derivative value
 */
double Component::computeSecondDerivative(double kappa)
{
  double A = ratioBesselFunction(kappa);
  double A1 = ratioBesselFunction_firstDerivative(kappa);
  double A2 = ratioBesselFunction_secondDerivative(kappa);
  double A3 = ratioBesselFunction_thirdDerivative(kappa);

  double kappa_sq = kappa * kappa;
  double second_derivative_msglen = 0;
  second_derivative_msglen += 4.0 * (1-kappa_sq) / ((1+kappa_sq)*(1+kappa_sq));
  second_derivative_msglen += 1 / (double) kappa_sq;
  second_derivative_msglen += N * A1;

  double tmp = (A * A2 - A1 * A1) / (double) (A*A);
  second_derivative_msglen += tmp;

  tmp = (A1 * A3 - A2 * A2) / (double) (2 * A1 * A1);
  second_derivative_msglen += tmp;

  return second_derivative_msglen;
}

/*!
 *  \brief This function computes the likelihood of a datum.
 *  \param x a reference to an array<double,2>
 *  \return the likelihood value
 */
double Component::likelihood(array<double,2> &x)
{
  return von_mises.density(x[0],x[1]);
}

/*!
 *  \brief This function computes the likelihood of a datum for a given kappa.
 *  \param x a reference to an array<double,2>
 *  \param kappa a double
 *  \return the likelihood value
 */
double Component::likelihood(array<double,2> &x, double kappa)
{
  return von_mises.density(x[0],x[1]);
}

/*!
 *  \brief This function computes the probability of the component parameters
 *  by discretizing the parameter space using the Fisher information.
 *  \return the probability value
 */
double Component::computeParametersProbability()
{
  double prior_density = computePriorDensity();
  double expected_fisher = computeFisherInformation();
  if (expected_fisher <= TOLERANCE) {
    return LARGE_NUMBER;
  } else {
    return prior_density/sqrt(expected_fisher);
  }
}

/*!
 *  \brief This function computes the prior density on the parameters.
 *  \return the joint prior density value
 */
double Component::computePriorDensity()
{
  double kappa_sq = kappa_mml * kappa_mml;
  double num = kappa_sq * fabs(sin(mu[0]));
  double denom = PI * PI * (1+kappa_sq) * (1+kappa_sq);
  return num/denom; 
}

/*!
 *  \brief This function computes the expected Fisher information using the
 *  component parameters.
 *  \return the expected Fisher value
 */
double Component::computeFisherInformation()
{
  if (kappa_mml <= TOLERANCE) {
    return 0;
  } else {
    double log_fisher = 0;
    log_fisher += 3 * log(N);
    log_fisher += 2 * log(kappa_mml);
    log_fisher += 2 * log(ratioBesselFunction(kappa_mml));
    log_fisher += log(ratioBesselFunction_firstDerivative(kappa_mml));
    log_fisher += 2 * log(fabs(sin(mu[0])));
    return exp(log_fisher);
  }
}

/*!
 *  \brief This function prints the optimal parameters of the component.
 *  \param os a reference to a ostream
 */
void Component::printParameters(ostream &os)
{
  os << "[mu,kappa]: [(" << mu[0] << "," << mu[1]<< ")," << kappa_mml << "]\n";
}

/*!
 *  \brief This function returns the optimum theta.
 *  \return theta that corresponds to the minimum message length
 */
double Component::getTheta()
{
  return mu[0];
}

/*!
 *  \brief This function returns the optimum phi.
 *  \return phi that corresponds to the minimum message length
 */
double Component::getPhi()
{
  return mu[1];
}

/*!
 *  \brief This function returns the optimum kappa.
 *  \return kappa that corresponds to the minimum message length
 */
double Component::getKappa()
{
  return kappa_mml;
}

