#include "Segment.h"

/*!
 *  \brief This is a constructor function.
 *  \param start an integer 
 *  \param end an integer
 *  \param cartesian a reference to a vector<Point<double>>
 *  \param spherical a reference to a vector<array<double,2>>
 */
Segment::Segment(int start, int end, vector<Point<double>> &cartesian,
                 vector<array<double,3>> &spherical) : start(start), end(end),
                cartesian(cartesian) , spherical(spherical)
{}

/*!
 *  \brief This function is used to set the initial two radii if the
 *  start index of the protein is zero.
 *  \param d1 a double
 *  \param d2 a double
 */
void Segment::setInitialDistances(double d1, double d2)
{
  radii[0] = d1;
  radii[1] = d2;
}

/*!
 *  \brief This function fits a null model to the protein segment.
 *  \return the optimal fit using the null model
 */
OptimalFit Segment::fitNullModel()
{
}

/*!
 *  \brief This function fits a null model to the protein segment.
 *  \param model an integer
 *  \return the optimal fit using the ideal model
 */
OptimalFit Segment::fitIdealModel(int model)
{
  switch(model) {
    
  }
}

