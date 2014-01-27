#include "Component.h"
#include "Support.h"

/*!
 *  \brief The null constructor module.
 */
Component::Component()
{}

/*!
 *  \brief This is a constructor function.
 *  \param mean_direction a reference to an vector<double>
 *  \param N a double
 *  \param constrain_kappa an integer
 */
Component::Component(vector<double> &mean_direction, double N, 
                     int constrain_kappa) : mean_direction(mean_direction),
                     N(N), constrain_kappa(constrain_kappa)
{
  updateMu(mean_direction);
}

/*!
 *  \brief This is a constructor function.
 *  \param unit_mean a reference to an vector<double>
 *  \param kappa a double
 */
Component::Component(vector<double> &unit_mean, double kappa):
                     unit_mean(unit_mean), kappa_ml(kappa),
                     kappa_mml(kappa)
{
  constrain_kappa = SET;
  von_mises = VonMises3D(unit_mean,kappa_mml);
  updateMu(unit_mean);
}

/*!
 *  \brief This function updates Mu the directions (theta/phi) in radians
 *  \param mean a reference to an vector<double>
 */
void Component::updateMu(vector<double> &mean)
{
  vector<double> spherical(3,0);
  cartesian2spherical(mean,spherical);
  mu[0] = spherical[1];
  mu[1] = spherical[2];
}

/*!
 *  \brief This module assigns an Component object on the rhs to one
 *  on the lhs
 *  \param source a reference to an Component object
 */
Component Component::operator=(const Component &source)
{
  if (this != &source) {
    mean_direction = source.mean_direction;
    unit_mean = source.unit_mean;
    mu = source.mu;
    N = source.N;
    rbar = source.rbar;
    R = source.R;
    constrain_kappa = source.constrain_kappa;
    kappa_ml = source.kappa_ml;
    kappa_mml = source.kappa_mml;
    von_mises = source.von_mises;
  }
  return *this;
}

/*!
 *  \brief This function is used to minimize the message length expression
 *  by finding the suitable model parameters.
 */
void Component::minimizeMessageLength()
{
  estimateVonMisesMean();
  kappa_ml = estimateKappa_ML();
  kappa_mml = estimateKappa_MML(kappa_ml);
  cout << "Kappa (MML): " << kappa_mml << endl;
  von_mises = VonMises3D(unit_mean,kappa_mml);
}

/*!
 *  \brief This function is used to estimate the Von Mises mean
 */
void Component::estimateVonMisesMean()
{
  cout << "\nCartesian coordinates of resultant mean direction vector: ";
  print(cout,mean_direction);
  vector<double> spherical_coordinates(3,0);
  cartesian2spherical(mean_direction,spherical_coordinates);

  cout << "Spherical coordinates of mean direction vector: ";
  print(cout,spherical_coordinates);

  // magnitude of the mean direction vector
  R = spherical_coordinates[0];

  // normalize the mean direction vector to obtain an unit vector
  unit_mean = vector<double>(3,0);
  for (int i=0; i<3; i++) {
    unit_mean[i] = mean_direction[i] / R;
  }
  cout << "Cartesian coordinates of unit mean direction vector: ";
  print(cout,unit_mean);

  cartesian2spherical(unit_mean,spherical_coordinates);
  mu[0] = spherical_coordinates[1];
  mu[1] = spherical_coordinates[2];
  cout << "Spherical coordinates of unit mean direction vector: ";
  print(cout,spherical_coordinates);

  rbar = R / (double)N ;
  cout << "R of direction vector: " << R << endl;
  cout << "N: " << N << endl;
  cout << "rbar: " << rbar << endl;

  // print the estimates of theta, phi
  cout << "\nEstimates:\n";
  cout << "(theta,phi) = (" << spherical_coordinates[1]*180/PI << ", " 
       << spherical_coordinates[2]*180/PI << ")" << endl;
}

/*!
 *  \brief This function computes the ML estimate (approximate) of kappa.
 *  (ref: http://en.wikipedia.org/wiki/Von_Mises%E2%80%93Fisher_distribution)
 *  \return the ML estimate of kappa
 */
double Component::estimateKappa_ML()
{
  double kappa = (rbar * (3 - (rbar * rbar))) / (1 - (rbar * rbar));
  if (constrain_kappa == SET) {
    if (fabs(kappa) >= MAX_KAPPA) {
      kappa = MAX_KAPPA;
    }
  }
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
    if (num_iterations > 1000) {
      if (constrain_kappa == SET && current >= MAX_KAPPA) {
        return MAX_KAPPA;
      } else {
        return prev;
      }
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
        if (constrain_kappa == SET && current >= MAX_KAPPA) {
          return MAX_KAPPA;
        } else {
          return current;
        }
      }
    } else {
      cout << "Derivative is zero ..." << endl;
      if (constrain_kappa == SET && prev >= MAX_KAPPA) {
        return MAX_KAPPA;
      } else {
        return prev;
      }
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
  if (constrain_kappa == SET) {
    double tmp = atan(MAX_KAPPA) - (MAX_KAPPA/(1+MAX_KAPPA*MAX_KAPPA));
    double constant = PI / (2 * tmp);
    part1 -= log(constant);
  }
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
                 + N * log(sinh(kappa)) - 2 * N * log(AOM);
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
 *  \param x a reference to an vector<double,2>
 *  \return the likelihood value
 */
double Component::likelihood(vector<double> &x)
{
  return von_mises.density(x);
}

/*!
 *  \brief This function computes the likelihood of a datum for a given kappa.
 *  \param x a reference to an vector<double,2>
 *  \param kappa a double
 *  \return the likelihood value
 */
double Component::likelihood(vector<double> &x, double kappa)
{
  return von_mises.density(x);
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
  double num = kappa_sq;// * fabs(sin(mu[0]));
  if (constrain_kappa == SET) {
    double tmp = atan(MAX_KAPPA) - (MAX_KAPPA/(1+MAX_KAPPA*MAX_KAPPA));
    double constant = PI / (2 * tmp);
    num *= constant;
  }
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
    //log_fisher += 2 * log(fabs(sin(mu[0])));
    return exp(log_fisher);
  }
}

/*!
 *  \brief This function prints the optimal parameters of the component.
 *  \param os a reference to a ostream
 */
void Component::printParameters(ostream &os)
{
  os << "[mu,kappa]: [(" << mu[0] * 180/PI << "," << mu[1] * 180/PI << ")," << kappa_mml << "]\n";
}

/*!
 *  \brief This function returns the mean direction of the distribution.
 *  \return the angular mean
 */
vector<double> Component::getMeanDirection()
{
  return unit_mean;
}

/*!
 *  \brief This function returns the optimum kappa.
 *  \return kappa that corresponds to the minimum message length
 */
double Component::getKappa()
{
  return kappa_mml;
}

/*!
 *  \brief This function is used to generate random samples.
 *  \param num_points an integer
 *  \return the list of generated points
 */
vector<vector<double>> Component::generate(int num_points)
{
  von_mises = VonMises3D(unit_mean,kappa_mml);
  return von_mises.generateCoordinates(num_points);
}

/*!
 *  \brief This function conflates two components.
 *  \param component a reference to a Component
 *  \return the conflated component
 */
Component Component::conflate(Component &component)
{
  vector<double> mean1 = unit_mean; 
  vector<double> mean2 = component.getMeanDirection();
  double k1 = kappa_mml;
  double k2 = component.getKappa();

  vector<double> resultant(3,0);
  for (int i=0; i<3; i++) {
    resultant[i] = k1 * mean1[i] + k2 * mean2[i];
  }
  vector<double> resultant_spherical(3,0);
  cartesian2spherical(resultant,resultant_spherical);
  double conflated_kappa = resultant_spherical[0];
  vector<double> spherical_mean(3,1);
  vector<double> conflated_mean(3,0);  // conflated
  spherical_mean[1] = resultant_spherical[1];
  spherical_mean[2] = resultant_spherical[2];
  spherical2cartesian(spherical_mean,conflated_mean);
  
  return Component(conflated_mean,conflated_kappa,constrain_kappa); 
}

