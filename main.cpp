#include "Support.h"

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.read_profiles == SET) {
    array<double,3> estimates = readProfiles(parameters.profiles_dir,parameters.res);
    cout << "Estimates: ";
    print(cout,estimates);
    //vonMisesDistribution_2DPlot(estimates,parameters.res);
  } else if (parameters.read_profiles == UNSET) {
    buildAngularProfile(parameters);
  }

  return 0;
}

