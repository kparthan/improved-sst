#include "Support.h"

extern string STRUCTURE;

int main(int argc, char **argv)
{
  getHomeAndCurrentDirectory();
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  STRUCTURE = extractName(parameters.file);

  if (parameters.dssp == SET && parameters.parse_dssp == SET) {
    parseDSSP(parameters);
  }

  if (parameters.dssp == UNSET && parameters.sst == UNSET && parameters.read_profiles == UNSET 
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

  if (parameters.dssp == SET && parameters.sst == SET) {
    assignSecondaryStructure(parameters.mixture_file,parameters.file,
                             parameters.orientation,parameters.portion_to_fit,
                             parameters.end_points,parameters.method);
  }

  return 0;
}

