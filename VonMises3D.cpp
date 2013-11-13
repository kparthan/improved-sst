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
 *  \brief This function computes the value of the distribution at a given 
 *  theta and phi (in degrees)
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
    double pr = constant * exp(exponent);
    //assert(pr <= 1);
    return pr;
  }
}

/*!
 *  \brief This function is used to generate samples from a distribution whose
 *  unit mean is directed along the Z-axis (canonical form).
 *  \param N an integer
 *  \return the list of random samples
 */
vector<array<double,3>> VonMises3D::generateCanonical(int N)
{
  auto ts = high_resolution_clock::now();
  usleep(1);
  auto te = high_resolution_clock::now();
  double t = duration_cast<nanoseconds>(ts-te).count();
  srand(t);
  double exponent = exp(-2 * kappa);
  double k_inv = 1 / (double)kappa;
  vector<array<double,3>> coordinates;
  for (int i=0; i<N; i++) {
    array<double,3> x;

    // generate probability value p \in [0,1]
    double p = rand() / (double) RAND_MAX;

    // generate a number (z coordinate) \in [-1,1]
    x[2] = 1 + k_inv * log(p + ((1-p) * exponent));
    assert(fabs(x[2]) <= 1);
    double theta = acos(x[2]);
    //scaleToAOM(&theta);

    // generate phi \in [0,2 PI]
    double phi = (rand() / (double) RAND_MAX) * 2 * PI;
    //scaleToAOM(&phi);

    // compute x and y coordinates
    x[0] = sin(theta) * cos(phi);
    x[1] = sin(theta) * sin(phi);
    coordinates.push_back(x);
  }
  return coordinates;
}

/*!
 *  \brief This function constructs the rotation matrix to transform points
 *  from the canonical representation to the current form.
 *  \return the rotation matrix
 */
Matrix<double> VonMises3D::constructRotationMatrix()
{
  vector<double> axis(3,0);
  // rotation by theta degrees
  double theta_rad = mu[0] * (PI/180.0);
  axis[1] = 1;  // rotation around Y-axis
  Vector<double> rot_axis1(axis);
  Matrix<double> m1 = rotationMatrix(rot_axis1,theta_rad); 

  // rotation by phi degrees
  axis = vector<double>(3,0);
  axis[2] = 1;  // rotation around Z-axis
  double phi_rad = mu[1] * (PI/180.0);
  Vector<double> rot_axis2(axis);
  Matrix<double> m2 = rotationMatrix(rot_axis2,phi_rad);
  
  // final rotation matrix
  Matrix<double> m = m2 * m1;
  return m;
}

/*!
 *  \brief This function is used to generate random samples (coordinates) from 
 *  the distribution.
 *  \param N an integer
 *  \return the list of random samples
 */
vector<array<double,3>> VonMises3D::generateCoordinates(int N)
{
  Matrix<double> rotate = constructRotationMatrix();
  vector<array<double,3>> canonical = generateCanonical(N);
  vector<array<double,3>> coordinates;
  array<double,3> r,x;
  for (int i=0; i<N; i++) {
    Point<double> p(canonical[i]);
    Point<double> pr = lcb::geometry::transform<double>(p,rotate);
    x[0] = pr.x();
    x[1] = pr.y();
    x[2] = pr.z();
    coordinates.push_back(x);
  }
  return coordinates;
}

/*!
 *  \brief This function is used to generate random samples (angles) from the 
 *  distribution.
 *  \param N an integer
 *  \return the list of random samples
 */
vector<array<double,2>> VonMises3D::generateAngles(int N)
{
  vector<array<double,3>> coordinates = generateCoordinates(N);
  vector<array<double,2>> angles;
  for (int i=0; i<N; i++) {
    Point<double> point(coordinates[i]);
    array<double,3> spherical = convertToSpherical(point);
    array<double,2> angle_pair({spherical[1],spherical[2]});
    angles.push_back(angle_pair);
  }
  return angles;
}

