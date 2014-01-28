#include "OptimalFit.h"

/*!
 *  \brief Null constructor
 */
OptimalFit::OptimalFit()
{
  message_length = LARGE_NUMBER; 
}

/*!
 *  \brief This module instantiates a new optimalInfo object
 *  \param ideal_model a reference to a IdealModel 
 *  \param message_length a double
 */
OptimalFit::OptimalFit(IdealModel &ideal_model, double message_length) : 
                       ideal_model(ideal_model), message_length(message_length)
{}

/*
 *  \brief This module is used to create a copy of an OptimalFit object.
 *  \param source a reference to an OptimalFit object
 */
OptimalFit::OptimalFit(const OptimalFit &source) :
                       ideal_model(source.ideal_model),
                       message_length(source.message_length)
{}

/*!
 *  \brief This module returns the stored message length.
 *  \return the message length
 */
double OptimalFit::getMessageLength() const
{
  return message_length;
}

/*!
 *  \brief This module assigns an OptimalFit object on the rhs to one
 *  on the lhs
 *  \param source a reference to an OptimalFit object
 */
OptimalFit OptimalFit::operator=(const OptimalFit &source)
{
  if (this != &source) {
    ideal_model = source.ideal_model;
    message_length = source.message_length;
  }
  return *this;
}

/*
 *  \brief This module compares the optimal message lengths of
 *  two segments.
 *  \param other a reference to an OptimalFit object
 *  \return true if the message length is less compared to that of
 *  the other segment
 */
bool OptimalFit::operator<(const OptimalFit &other)
{
  if (message_length < other.getMessageLength()) {
    return 1;
  } else {
    return 0;
  }
}

/*!
 *  \brief This function prints the details of the optimal fit.
 */
void OptimalFit::printFitInfo()
{
  ideal_model.printModelInfo();
  cout << "total msglen: " << message_length << endl;
}

/*!
 *  \brief This function returns the name of the representative ideal model.
 *  \return the name
 */
string OptimalFit::getName()
{
  return ideal_model.getName();
}

