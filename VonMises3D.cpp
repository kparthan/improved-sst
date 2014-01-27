#include "VonMises3D.h"
#include "Support.h"
#include "Geometry3D.h"

extern vector<double> YAXIS;
extern vector<double> ZAXIS;

/*!
 *  \brief This is a null constructor module
 *  sets default values of mean vector as 0 and kappa as 1
 */
VonMises3D::VonMises3D() : unit_mean(vector<double>({0,0,1})),kappa(1)
{
  computeNormalizationConstant();
}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  standard deviation of the distribution
 *  \param unit_mean a reference to a vector<double>
 *  \param kappa a double
 */
VonMises3D::VonMises3D(vector<double> &unit_mean, double kappa) : 
                       unit_mean(unit_mean), kappa(kappa)
{
  computeNormalizationConstant();
}

/*!
 *  \brief This function pre computes the normalization constant.
 */
void VonMises3D::computeNormalizationConstant()
{
  double tmp = exp(-2 * kappa);
  double denom = 2 * PI * (1-tmp);
  norm_constant = kappa / denom;
}

/*!
 *  \brief This function assigns a source VonMises3D distribution.
 *  \param source a reference to a VonMises3D
 */
VonMises3D VonMises3D::operator=(const VonMises3D &source)
{
  if (this != &source) {
    unit_mean = source.unit_mean;
    kappa = source.kappa;
    norm_constant = source.norm_constant;
  }
  return *this;
}

/*!
 *  \brief This function returns the mu of the distribution
 *  \return the mean of the distribution
 */
const vector<double> VonMises3D::mean(void)
{
	return unit_mean;
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
 *  \brief This function computes the value of the distribution at a given x. 
 *  \param x a reference to a vector<double> 
 *  \return the probability density
 */
double VonMises3D::density(vector<double> &x)
{
  if (kappa <= TOLERANCE) {
    return 1.0/(4*PI);
  } else {
    double dot_product;
    computeDotProduct(unit_mean,x,dot_product);
    double exponent = kappa * (dot_product - 1);
    double pr = norm_constant * exp(exponent);
    /*if (pr < ZERO) {
      pr = ZERO;
    }*/
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
vector<vector<double>> VonMises3D::generateCanonical(int N)
{
  auto ts = high_resolution_clock::now();
  usleep(1);
  auto te = high_resolution_clock::now();
  double t = duration_cast<nanoseconds>(ts-te).count();
  srand(t);
  double exponent = exp(-2 * kappa);
  double k_inv = 1 / (double)kappa;
  vector<vector<double>> coordinates;
  for (int i=0; i<N; i++) {
    vector<double> x(3,0);

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
vector<vector<double>> VonMises3D::constructRotationMatrix()
{
  vector<double> spherical(3,0);
  cartesian2spherical(unit_mean,spherical);

  vector<vector<double>> rotate1,rotate2,rotation;
  initializeMatrix(rotate1,3,3);
  initializeMatrix(rotate2,3,3);
  initializeMatrix(rotation,3,3);

  // rotation around Y-axis (by theta radians)
  computeRotationMatrix(YAXIS,spherical[1],rotate1);

  // rotation around Z-axis (by phi radians)
  computeRotationMatrix(ZAXIS,spherical[2],rotate2);

  // final rotation matrix
  multiplyVectors(rotate2,rotate1,rotation);

  return rotation;
}

/*!
 *  \brief This function is used to generate random samples (coordinates) from 
 *  the distribution.
 *  \param N an integer
 *  \return the list of random samples
 */
vector<vector<double>> VonMises3D::generateCoordinates(int N)
{
  vector<vector<double>> rotation_matrix = constructRotationMatrix();
  vector<vector<double>> canonical = generateCanonical(N);
  vector<vector<double>> coordinates;
  vector<double> x(3,0);
  for (int i=0; i<N; i++) {
    rotateVector(rotation_matrix,canonical[i],x);
    coordinates.push_back(x);
  }
  return coordinates;
}

/*!
 *  \brief This function is used to generate random samples (angles) from the 
 *  distribution (angles generated are measured in degrees).
 *  \param N an integer
 *  \return the list of random samples
 */
vector<vector<double>> VonMises3D::generateAngles(int N)
{
  vector<vector<double>> coordinates = generateCoordinates(N);
  vector<vector<double>> angles;
  vector<double> spherical(3,0);
  vector<double> angles_pair(2,0);
  for (int i=0; i<N; i++) {
    cartesian2spherical(coordinates[i],spherical);
    angles_pair[0] = spherical[1];
    angles_pair[1] = spherical[2];
    angles.push_back(angles_pair);
  }
  return angles;
}

