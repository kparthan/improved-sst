#include "Message.h"
extern double MSGLEN_RADIUS,MSGLEN_CELL,MSGLEN_VMF_CELL;

/*!
 *  \brief Null constructor module.
 */
Message::Message()
{}

/*!
 *  \brief This module computes the statement cost to communicate an
 *  integer over a log star distribution. 
 *  \return the length of encoding (in bits)
 */
double Message::encodeUsingLogStarModel(double value)
{
  double result = 0;
  if( value < 1 ){
    throw range_error("Not a positive real integer");
  }
  double partial = log2(value);
  while(partial > 0){
    result += partial;
    partial = log2(partial);
  }
  return (result + log2(2.865));
}

/*!
 *  \brief This function computes the message length to encode using a Poisson
 *  distribution.
 *  \param x a int
 *  \param poisson a reference to a Poisson
 *  \return the message length
 */
double Message::encodeUsingPoissonModel(int x, Poisson &poisson)
{
  //double density = poisson.density(x);
  //return -log2(density);
  double msglen = 0;
  double lambda = poisson.mean();
  msglen += lambda - (x * log(lambda));
  long double fact = boost::math::factorial<long double>(x);
  msglen += log(fact);
  return msglen/log(2);
}

/*!
 *  \brief This function encodes the length assuming a Normal distribution.
 *  \param length a double
 *  \param normal a reference to a Normal
 *  \return the message length to encode the length (in bits)
 */
double Message::encodeUsingNormalModel(double length, Normal &normal)
{
  double density = normal.density(length);
  double msglen = -log2(AOM) - log2(density);
  return msglen;
}

/*!
 *  \brief This function computes the message radius to communicate a point
 *  using the sphere model.
 *  \param radius a double
 *  \param normal a reference to a Normal
 *  \return the message length (in bits)
 */
double Message::encodeUsingSphereModel(double radius, Normal &normal)
{
  double msglen = 0;
  // encode the radius
  msglen += encodeUsingNormalModel(radius,normal);
  MSGLEN_RADIUS += msglen;
  // encode the position on the sphere
  double cell = 2 + LOG2_PI + 2 * log2(radius) - 2*log2(AOM);
  MSGLEN_CELL += cell;
  msglen += cell;
  return msglen;
}

double Message::encodeUsingSphereModel(double radius, double theta, Normal &normal)
{
  // encode the radius
  double msglen = encodeUsingNormalModel(radius,normal);
  MSGLEN_RADIUS += msglen;
  // encode the position on the sphere
  double cell = 2 + LOG2_PI -log2(sin(theta)) - 2*log2(AOM);
  MSGLEN_CELL += cell;
  msglen += cell;
  return msglen;
}

/*!
 *  \brief This function is used to compute the encoding length of a direction
 *  using a mixture model.
 *  \param direction a reference to a vector<double>
 *  \param mixture a reference to a Mixture
 *  \return the message length (in bits)
 */
double Message::encodeUsingMixtureModel(vector<double> &direction, 
                                        Mixture &mixture)
{
  double density = mixture.probability(direction); 
  double msglen = -log2(density) - 2*log2(AOM);
  return msglen;
}

// scaled AOM
double Message::encodeUsingMixtureModel(vector<double> &direction, 
                                        Mixture &mixture, double radius)
{
  double density = mixture.probability(direction); 
  double msglen = -log2(density) - 2*log2(AOM) + 2 * log2(radius);
  MSGLEN_VMF_CELL += msglen;
  return msglen;
}

/*!
 *  \brief This function is used to compute the encoding length of a direction
 *  using a single component.
 *  \param direction a reference to a vector<double>
 *  \param component a reference to a Component
 *  \return the message length (in bits)
 */
double Message::encodeUsingComponent(vector<double> &direction,
                                     Component &component)
{
  double density = component.likelihood(direction);
  double msglen = -log2(density) - 2*log2(AOM);
  return msglen;
}

