#include "Support.h"

int main(int argc, char **argv)
{
  //cout << "max kappa: " << MAX_KAPPA << endl;
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.load_mixture == SET) {
    visualizeMixtureComponents(parameters);
  } else {
    if (parameters.read_profiles == SET) {
      computeEstimators(parameters);
    } else if (parameters.read_profiles == UNSET) {
      buildAngularProfile(parameters);
    }
  }
  return 0;
}

