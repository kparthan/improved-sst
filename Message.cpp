#include "Message.h"
#include "Support.h"

/*!
 *  \brief The null constructor module.
 */
Message::Message()
{}

/*!
 *  \brief This is a constructor function.
 *  \param direction a reference to an array<double,3>
 *  \param N an integer
 */
Message::Message(array<double,3> &mean_direction, double N) :
                 mean_direction(mean_direction), N(N)
{
  //double lattice_constant = 0.08;
  //constant = 1.5 * (1 + log(lattice_constant) + log(N)) + (N+2) * log(PI)
  //           + N * log(4);
}

/*!
 *  \brief This function is used to minimize the message length expression
 *  by finding the suitable model parameters.
 */
void Message::minimize()
{
  estimateVonMisesMean();
  kappa_ml = estimateKappa_ML();
  kappa_mml = estimateKappa_MML(kappa_ml);
}

/*!
 *  \brief This function is used to estimate the Von Mises mean
 */
void Message::estimateVonMisesMean()
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
double Message::estimateKappa_ML()
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
double Message::estimateKappa_MML(double initial)
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
    if (fabs(fx_der) > ZERO) {
      current = prev - (fx/(double)fx_der);
      cout << "Iteration " << num_iterations << ": [" << prev << ", " 
           << current <<  ", " << fx << ", " << fx_der << "]" << endl;
      if (fabs(current - prev) > ZERO) {
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
double Message::computeMessageLength(double kappa)
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
double Message::computeFirstPart(double kappa)
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
double Message::computeSecondPart(double kappa)
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
double Message::computeFirstDerivative(double kappa)
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
double Message::computeSecondDerivative(double kappa)
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

