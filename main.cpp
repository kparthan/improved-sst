#include "Support.h"

extern string STRUCTURE;

int main(int argc, char **argv)
{
  getHomeAndCurrentDirectory();
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  STRUCTURE = extractName(parameters.file);

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
    assignSecondaryStructure(parameters.mixture_file,parameters.file,
                             parameters.orientation,parameters.portion_to_fit,
                             parameters.end_points);
  }

  return 0;
}

