#include "Segmentation.h"

/*!
 *  \brief Null constructor module.
 */
Segmentation::Segmentation()
{}

/*!
 *  \brief Constructor
 *  \param structure a reference to a ProteinStructure
 *  \param chain_id a string
 *  \param segmentation a reference to a vector<int>
 *  \param names a reference to a vector<string>
 */
Segmentation::Segmentation(ProteinStructure *structure, string chain_id,
                           vector<int> &segmentation, vector<string> &names) : 
                           structure(structure), chain_id(chain_id)
{
  assert(segmentation.size() == names.size()+1);
  array<int,2> segment;
  segment[0] = segmentation[0];
  for (int i=1; i<segmentation.size(); i++) {
    segment[0] = segmentation[i-1];
    segment[1] = segmentation[i];
    if (i != 1 && names[i-2].compare("strand") == 0 && names[i-1].compare("strand") == 0) {
      segments[segments.size()-1][1] = segment[1];
    } else {
      segments.push_back(segment);
      model_names.push_back(names[i-1]);
    }
  }
  assert(segments.size() == model_names.size());
  for (int i=0; i<segments.size(); i++) {
    cout << segments[i][0]+1 << "\t" << segments[i][1]+1 << "\t" << model_names[i] << endl;
  }
}

/*!
 *  \brief This function is used to get the beta strands.
 *  \return the list of strands
 */
vector<array<int,2>> Segmentation::getStrandSegments()
{
  vector<array<int,2>> strand_segments;
  for (int i=0; i<model_names.size(); i++) {
    if (model_names[i].compare("strand") == 0) {
      array<int,2> segment = segments[i];
      strand_segments.push_back(segment);
    }
  }
  return strand_segments;
}

/*!
 *  \brief This function is used to trim the strand assignments as a
 *  postprocessing step.
 */
void Segmentation::postprocess()
{
  structure->undoLastSelection();
  structure->select(FullBackboneSelector());
  Chain chain = structure->getDefaultModel()[chain_id];
  /*vector<string> residue_ids = chain.getResidueIdentifiers();
  vector<Residue> residues;
  for (int i=0; i<residue_ids.size(); i++) {
    Residue res = chain[residue_ids[i]];
    residues.push_back(res);
  }*/

  vector<array<int,2>> strand_segments = getStrandSegments();
  vector<vector<Residue>> strands;
  for (int i=0; i<strand_segments.size(); i++) {
    int start = strand_segments[i][0]+1;
    int end = strand_segments[i][1]+1;
    vector<Residue> residues;
    for (int j=start; j<=end; j++) {
      string res_id = boost::lexical_cast<string>(j);
      Residue res = chain[res_id];
      residues.push_back(res);
    }
    strands.push_back(residues);
  }
  cout << "# of strands: " << strands.size() << endl;

  string CA = "CA";
  string C = "C";
  string N = "N";
  string O = "O";

  array<double,3> coords;
  vector<array<double,3>> coords_N,coords_O;
  vector<vector<array<double,3>>> all_coords_N,all_coords_O;

  for (int i=0; i<strands.size(); i++) {
    for (int j=0; j<strands[i].size(); j++) {
      Residue residue = strands[i][j];
      Atom atomN = residue[N];
      coords = atomN.getAtomicCoordinate<double>();
      coords_N.push_back(coords);
      Atom atomO = residue[O];
      coords = atomO.getAtomicCoordinate<double>();
      coords_O.push_back(coords);
    }
    all_coords_N.push_back(coords_N);
    all_coords_O.push_back(coords_O);
  }

  vector<array<double,3>> other_coords_N,other_coords_O;
  for (int i=0; i<strands.size(); i++) {
    coords_N = all_coords_N[i];
    coords_O = all_coords_O[i];
    for (int j=0; j<strands.size(); j++) {
      if (j != i) {
        other_coords_N = all_coords_N[j];
        other_coords_O = all_coords_O[j];
        
      }
    }
  }
}

