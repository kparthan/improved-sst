#include "Support.h"

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.comparison == SET) {
    if (parameters.standardize == SET) {
      if (parameters.standardize_status == UNSET) {
        standardizePremeasures(parameters);
      } 
    }
    if (parameters.comparison_files.size() > 0) {
      compareStructuresList(parameters);
    }
  } else if (parameters.comparison == UNSET) {
    build(parameters);
  }

  return 0;
}

