#ifndef OPTIMAL_FIT_H
#define OPTIMAL_FIT_H

#include "IdealModel.h"

class OptimalFit
{
  private:
    //! The respresentative ideal model;
    IdealModel ideal_model;

    //! Message length
    double message_length;

  public:
    //! Null constructor
    OptimalFit();

    //! Constructor
    OptimalFit(IdealModel &, double);

    //! Copy constructor
    OptimalFit(const OptimalFit &);

    //! Gets the message length
    double getMessageLength() const;

    //! Assignment operator
    OptimalFit operator=(const OptimalFit &);

    //! Compares the optimal message length of two segments
    bool operator<(const OptimalFit &);

    //! Prints the details of the optimal fit
    void printFitInfo();
};

#endif

