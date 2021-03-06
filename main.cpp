#include "Support.h"

int main(int argc, char **argv)
{
  //cout << "max kappa: " << MAX_KAPPA << endl;
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.sst == UNSET && parameters.read_profiles == UNSET 
      && parameters.mixture_model == UNSET) {
    buildAngularProfile(parameters);
  }

  if (parameters.read_profiles == SET) {
    computeEstimators(parameters);
  } 

  if (parameters.simulation == SET) {
    simulateMixtureModel(parameters);
  }

  if (parameters.load_mixture == SET && parameters.simulation == UNSET) {
    visualizeMixtureComponents(parameters);
  } 

  if (parameters.load_mixture == SET && parameters.sst == SET) {
    assignSecondaryStructure(parameters.mixture_file,parameters.structure,
                             parameters.orientation);
  }

  return 0;
}

