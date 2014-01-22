#include "Message.h"

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
  // encode the radius
  double msglen = encodeUsingNormalModel(radius,normal);
  // encode the position on the sphere
  msglen += 2 + LOG2_PI + 2 * log2(radius) - 2*log2(AOM);
  return msglen;
}

/*!
 *  \brief This function is used to compute the encoding length of a direction
 *  using a mixture model.
 *  \param direction a reference to a array<double,2>
 *  \param mixture a reference to a Mixture
 *  \return the message length (in bits)
 */
double Message::encodeUsingMixtureModel(array<double,2> &direction, 
                                        Mixture &mixture)
{
  double density = mixture.probability(direction); 
  double msglen = -log2(density) - 2*log2(AOM);
  return msglen;
}

