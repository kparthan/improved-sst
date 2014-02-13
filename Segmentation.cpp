#include "Segmentation.h"
#include "Support.h"

#define HYDROGEN_BOND_THRESHOLD 3.2 

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

  vector<vector<Atom>> all_atoms_N,all_atoms_O;
  for (int i=0; i<strands.size(); i++) {
    vector<Atom> atoms_N,atoms_O;
    for (int j=0; j<strands[i].size(); j++) {
      Residue residue = strands[i][j];
      Atom atomN = residue["N"];
      atoms_N.push_back(atomN);
      Atom atomO = residue["O"];
      atoms_O.push_back(atomO);
    }
    all_atoms_N.push_back(atoms_N);
    all_atoms_O.push_back(atoms_O);
  }
  assert(all_atoms_N.size() == all_atoms_O.size());
  for (int i=0; i<all_atoms_N.size(); i++) {
    assert(all_atoms_N[i].size() == all_atoms_O[i].size());
    //cout << all_atoms_N[i].size() << endl;
    // print coordinates
    /*array<double,3> coords;
    for (int j=0; j<all_atoms_N[i].size(); j++) {
      Atom N = all_atoms_N[i][j];
      coords = N.getAtomicCoordinate<double>();
      print(cout,coords);
      Atom O = all_atoms_O[i][j];
      coords = O.getAtomicCoordinate<double>();
      print(cout,coords);
      cout << endl;
    }
    cout << endl;*/
  }

  ofstream log("pruning.log");
  // initialize pruning matrix
  vector<bool> consider(strands.size(),1);
  vector<vector<bool>> pruned;
  for (int i=0; i<strands.size(); i++) {
    vector<bool> tmp(strands[i].size(),1);
    pruned.push_back(tmp);
  }
  Atom current_N,current_O,other_N,other_O;
  for (int i=0; i<strands.size(); i++) {
    log << "Strand #" << i+1 << endl;
    log << "\tPruning from the start ...\n";
    int j = 0;
    while (j < strands[i].size()) { // pruning from the segment start
      log << "\tCURRENT RESIDUE #" << j+1 << endl;
      current_N = all_atoms_N[i][j];
      current_O = all_atoms_O[i][j];
      // check whether these atoms are close to any other atoms
      for (int k=0; k<strands.size(); k++) {
        if (k != i && consider[k] == 0) {
          log << "\t\t" << k+1 << " is not a strand anymore ...\n";
        }
        if (k != i && consider[k] == 1) {
          log << "\t\tComparing with residues in strand #" << k+1 << endl;
          for (int m=0; m<strands[k].size(); m++) {
            if (pruned[k][m] == 0) {
              log << "\t\t\tResidue " << m+1 << " is pruned already ...\n";
            } else if (pruned[k][m] == 1) {
              log << "\t\t\tResidue " << m+1 << ":\n";
              other_N = all_atoms_N[k][m];
              other_O = all_atoms_O[k][m];
              bool check = checkProximity(current_N,current_O,other_N,other_O,log);
              if (check == 1) { // stop and begin pruning from the other end
                log << "\t\t\tCurrent residue is retained ...\n";
                pruned[i][j] = 1;
                goto prune_other_end;
              } else if (check == 0) {
                pruned[i][j] = 0;
                if (j == strands[i].size() - 1) {
                  log << "\t\tStrand #" << i+1 << " discarded completely ...\n";
                  consider[i] = 0;
                  goto finish;
                }
              }
            }
          }
        }
      }
      if (pruned[i][j] == 0) {
        log << "\t\t\tCurrent residue is pruned ...\n";
      }
      j++;
    }
    prune_other_end:
    j = strands[i].size() - 1;
    log << "\tPruning from the end ...\n";
    while (j >= 1) {  // pruning from the segment end
      log << "\tCURRENT RESIDUE #" << j+1 << endl;
      current_N = all_atoms_N[i][j];
      current_O = all_atoms_O[i][j];
      // check whether these atoms are close to any other atoms
      for (int k=0; k<strands.size(); k++) {
        if (k != i && consider[k] == 0) {
          log << "\t\t" << k+1 << " is not a strand anymore ...\n";
        }
        if (k != i && consider[k] == 1) {
          log << "\t\tComparing with residues in strand #" << k+1 << endl;
          for (int m=0; m<strands[k].size(); m++) {
            if (pruned[k][m] == 0) {
              log << "\t\t\tResidue " << m+1 << " is pruned already ...\n";
            } else if (pruned[k][m] == 1) {
              log << "\t\t\tResidue " << m+1 << ":\n";
              other_N = all_atoms_N[k][m];
              other_O = all_atoms_O[k][m];
              bool check = checkProximity(current_N,current_O,other_N,other_O,log);
              if (check == 1) { // stop and begin pruning from the other end
                log << "\t\t\tCurrent residue is retained ...\n";
                pruned[i][j] = 1;
                goto finish;
              } else if (check == 0) {
                pruned[i][j] = 0;
              }
            }
          }
        }
      }
      if (pruned[i][j] == 0) {
        log << "\t\t\tCurrent residue is pruned ...\n";
      }
      j--;
    }
    finish:
    for (int j=0; j<pruned[i].size(); j++) {
      cout << pruned[i][j] << " ";
    }
    cout << endl;
    log << endl;
  }
  log.close();
}

/*!
 *  \brief This function checks whether a residue in one strand is involved in 
 *  hydrogen bonding with any other strand.
 *  \param current_N a reference to an Atom
 *  \param current_O a reference to an Atom
 *  \param other_N a reference to an Atom
 *  \param other_N a reference to an Atom
 *  \param log a reference to a ostream
 *  \return a boolean
 */
bool Segmentation::checkProximity(Atom &current_N, Atom &current_O,
                                  Atom &other_N, Atom &other_O, ostream &log)
{
  double d1,d2;
  // check current_N with other_O
  d1 = distance<double>(current_N,other_O);
  log << "\t\t\t\tdist(current_N,other_O) : " << d1 << endl;
  // check current_O with other_N
  d2 = distance<double>(current_O,other_N); 
  log << "\t\t\t\tdist(current_O,other_N) : " << d2 << endl;
  if (d1 <= HYDROGEN_BOND_THRESHOLD || d2 <= HYDROGEN_BOND_THRESHOLD) {
    return 1;
  } else {
    return 0;
  }
}

