#include "Support.h"

int main(int argc, char **argv)
{
  //cout << "max kappa: " << MAX_KAPPA << endl;
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.read_profiles == UNSET && parameters.mixture_model == UNSET) {
    buildAngularProfile(parameters);
  }

  if (parameters.read_profiles == SET && parameters.mixture_model == SET) {
    computeEstimators(parameters);
  } 

  if (parameters.simulation == SET) {
    simulateMixtureModel(parameters);
  }

  if (parameters.load_mixture == SET) {
    visualizeMixtureComponents(parameters);
  } 

  return 0;
}

