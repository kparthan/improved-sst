#include "Segment.h"
#include "Message.h"
#include "Superpose3D.h"
#include "Geometry3D.h"

extern vector<double> ZAXIS;
extern int DEBUG;

/*!
 *  \brief This is a constructor function.
 *  \param start an integer 
 *  \param end an integer
 *  \param cartesian a reference to a vector<vector<double>>
 *  \param spherical a reference to a vector<vector<double,2>>
 */
Segment::Segment(int start, int end, vector<vector<double>> &cartesian_coordinates,
                 vector<vector<double>> &spherical_coordinates, 
                 vector<vector<double>> &unit_coordinates) : start(start), end(end),
                 cartesian_coordinates(cartesian_coordinates), 
                 spherical_coordinates(spherical_coordinates),
                 unit_coordinates(unit_coordinates)
{
  num_residues = end - start + 1;
  for (int i=start; i<=end; i++) {
    observed_residues.push_back(cartesian_coordinates[i]);
  }
}

/*!
 *  \brief This function is used to set the initial two radii if the
 *  start index of the protein is zero.
 *  \param d1 a double
 *  \param d2 a double
 */
void Segment::setInitialDistances(double d1, double d2)
{
  radii = vector<double>(2,0);
  radii[0] = d1;
  radii[1] = d2;
}

/*!
 *  \brief This function fits a null model to the protein segment.
 *  \param mixture a reference to a Mixture
 *  \return the optimal fit using the null model
 */
OptimalFit Segment::fitNullModel(Mixture &mixture)
{
  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
  Message message;
  double msglen = 0;
  int begin_loop = start;
  
  // state the length of segment
  msglen += message.encodeUsingLogStarModel(num_residues);
  if (start == 0) {
    // if segment begins with the first point in the protein, the receiver
    // has to wait until the first point(origin) and the next two points are
    // transmitted, using the sphere model
    msglen += message.encodeUsingSphereModel(radii[0],normal);
    if (end > 1) {
      msglen += message.encodeUsingSphereModel(radii[1],normal);
      begin_loop = 2;
    } else {
      begin_loop = 1;
    }
  } else if (start == 1) {
    msglen += message.encodeUsingSphereModel(radii[1],normal);
    begin_loop = 2;
  }
  // state the remaining points using the mixture model
  double r;
  for (int i=begin_loop; i<end; i++) {
    // state radius
    r = spherical_coordinates[i-2][0];
    msglen += message.encodeUsingNormalModel(r,normal);
    // state direction
    msglen += message.encodeUsingMixtureModel(unit_coordinates[i-2],mixture);
  }

  string name = "NullModel";
  IdealModel null_model(num_residues,name);
  return OptimalFit(null_model,msglen);
}

/*!
 *  \brief This function fits a null model to the protein segment.
 *  \param model a reference to a IdealModel;
 *  \param mixture a reference to a Mixture
 *  \param orientation an integer
 *  \return the optimal fit using the ideal model
 */
OptimalFit Segment::fitIdealModel(IdealModel &model, Mixture &mixture,
                                  int orientation)
{
  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
  Message message;
  int begin_loop = start;

  double msglen = 0;
  // state the length of segment
  msglen += message.encodeUsingLogStarModel(num_residues);

  if (start == 0) { // the first segment
    // if segment begins with the first point in the protein, the receiver
    // has to wait until the first point(origin) and the next two points are
    // transmitted, using the sphere model
    msglen += message.encodeUsingSphereModel(radii[0],normal);
    if (end > 1) {
      msglen += message.encodeUsingSphereModel(radii[1],normal);
      begin_loop = 2;
    }
    begin_loop = 1;
  } else if (start == 1) {
    msglen += message.encodeUsingSphereModel(radii[1],normal);
    begin_loop = 2;
  } else {  // an intermediate segment
    // the start point of the segment is the last point of the previous segment
    // the next two points in the segment are stated using the null model
    double r;
    for (int i=start; i<start+2; i++) {
      // state radius
      r = spherical_coordinates[i-2][0];
      msglen += message.encodeUsingNormalModel(r,normal);
      // state direction
      msglen += message.encodeUsingMixtureModel(unit_coordinates[i-2],mixture);
    }
    begin_loop = start + 2;
  }

  if (num_residues > 3) {
    // get the ideal model of length equal to the segment length
    vector<vector<double>> ideal_residues = model.getResidues(num_residues);
    // [observed,MOVING] and [ideal,FIXED]
    vector<vector<double>> observed,ideal;  // observed, ideal are running lists
                                            // of the orginal untransformed
                                            // coordinates
    vector<vector<double>> observed_copy,ideal_copy;
    suffStatClass suff_stats;
    vector<vector<double>> four_mer,transformed_four_mer;
    vector<vector<double>> rotation_matrix,zaxis_transform;
    initializeMatrix(four_mer,4,3);
    initializeMatrix(transformed_four_mer,4,3);
    initializeMatrix(rotation_matrix,3,3);
    initializeMatrix(zaxis_transform,3,3);
    vector<double> unit_mean(3,0);
    vector<double> x(3,0);
    vector<double> vonmises_suffstats(3,0);
    double kappa = 5;
    Mixture conflated_mixture;

    // INITIAL SUPERPOSITION
    for (int i=0; i<3; i++) {
      observed.push_back(observed_residues[i]);
      ideal.push_back(ideal_residues[i]);
    }
    observed_copy = observed;
    ideal_copy = ideal;
    Superpose3DClass superpose(observed_copy,ideal_copy);
    cout << "RMSD:" << superpose.getRMSD() << endl;
    suff_stats = superpose.getSufficientStatistics();
    observed_copy.push_back(observed_residues[3]);
    superpose.transformVectors(observed_copy);
    for (int i=0; i<4; i++) {
      four_mer[i] = observed_copy[i];
    }
    convertToCanonicalForm(four_mer,transformed_four_mer,rotation_matrix);
    getCurrentMeanAndDirection(four_mer[2],rotation_matrix,transformed_four_mer[2],
                               transformed_four_mer[3],ideal_residues[2],
                               ideal_residues[3],orientation,unit_mean,x);
    ofstream debug("debug",ios::app);
    ofstream msg("msg",ios::app);
    if (DEBUG == SET) {
      print(debug,unit_mean);print(debug,x);debug<<endl;
    }
    double MSG;
    Component component(unit_mean,kappa);
    conflated_mixture = mixture.conflate(component);
    MSG = message.encodeUsingMixtureModel(x,conflated_mixture);
    msglen += MSG;//msg << MSG << endl;
    alignWithZAxis(ideal_residues[2],ideal_residues[3],zaxis_transform);
    applyIdealModelTransformation(zaxis_transform,observed_copy[2],
                                  observed_copy[3],vonmises_suffstats);
    int N = 1;

    // ADAPTIVE SUPERPOSITION
    for (int om=start+3; om<end; om++) {
      // start,...,om : points known to receiver
      //           om : most recent point communicated to the receiver
      //         om+1 : current point being transmitted
      // superpose (start,...,om) with the ideal model's (i1,...,im)
      observed.push_back(observed_residues[om-start]);
      ideal.push_back(ideal_residues[om-start]);
      // copy the contents of the observed,ideal into temporary containers so that
      // these temp lists can be transformed without affecting the original lists
      observed_copy = observed;
      ideal_copy = ideal;
      //writeToFile(observed_copy,"observed_copy_before");
      //writeToFile(ideal_copy,"ideal_copy_before");
      Superpose3DClass superpose(suff_stats,observed_copy[om-start],
                                 ideal_copy[om-start]);
      cout << "\"RMSD\":" << superpose.getRMSD() << endl;
      suff_stats = superpose.getSufficientStatistics();
      observed_copy.push_back(observed_residues[om-start+1]);
      superpose.transformVectors(observed_copy);
      //writeToFile(observed_copy,"observed_copy_after");
      //writeToFile(ideal_copy,"ideal_copy_after");
      // get the current four_mer
      for (int i=0; i<4; i++) {
        four_mer[i] = observed_copy[N+i];
      }
      // apply the canonical transformation on the ideal four_mer
      // it is sufficient to transform im and i_{m+1}
      convertToCanonicalForm(four_mer,transformed_four_mer,rotation_matrix);
      getCurrentMeanAndDirection(four_mer[2],rotation_matrix,transformed_four_mer[2],
                                 transformed_four_mer[3],ideal_residues[om-start],
                                 ideal_residues[om-start+1],orientation,unit_mean,x);
      if (N > 1) {
        Component adaptive_component(vonmises_suffstats,N+1,SET);
        adaptive_component.minimizeMessageLength(ZAXIS);
        kappa = adaptive_component.getKappa();
      }
      Component component(unit_mean,kappa);
      conflated_mixture = mixture.conflate(component);
      MSG = message.encodeUsingMixtureModel(x,conflated_mixture);
      msglen += MSG;
      if (DEBUG == SET) {
        print(debug,unit_mean);print(debug,x);debug<<endl;
        msg << MSG << endl;
      }
      // update von mises suff stats
      alignWithZAxis(ideal_residues[om-start],ideal_residues[om-start+1],zaxis_transform);
      applyIdealModelTransformation(zaxis_transform,observed_copy[om-start],
                                    observed_copy[om-start+1],vonmises_suffstats);
      N++;
    }
    if (DEBUG == SET) {
      debug << endl;debug.close();
      msg << endl;msg.close();
    }
  }
  ProteinStructure *p = model.getStructure();
  string name = model.getName();
  IdealModel ideal_model(p,name);
  ideal_model.setLength(num_residues);
  return OptimalFit(ideal_model,msglen);
}

/*!
 *  \brief This function computes the current mean of the VMF component and 
 *  computes the direction from the mean to be used in the density function.
 *  \param pre_origin a reference to a vector<double>
 *  \param rotation_matrix a reference to a vector<vector<double>>
 *  \param om_new a reference to a vector<double>
 *  \param om_plus_1_new a reference to a vector<double>
 *  \param im a reference to a vector<double>
 *  \param im_plus_1 a reference to a vector<double>
 *  \param orientation an integer
 *  \param unit_mean a reference to a vector<double>
 *  \param x a reference to a vector<double>
 */
void Segment::getCurrentMeanAndDirection(
  vector<double> &pre_origin, 
  vector<vector<double>> &rotation_matrix,
  vector<double> &om_new, 
  vector<double> &om_plus_1_new, 
  vector<double> &im, 
  vector<double> &im_plus_1, 
  int orientation, 
  vector<double> &unit_mean, 
  vector<double> &x
) {
  // transform im and i_{m+1}
  // translation
  vector<double> im_tmp(3,0);
  vector<double> im_plus_1_tmp(3,0);
  for (int i=0; i<3; i++) {
    im_tmp[i] = im[i] - pre_origin[i];
    im_plus_1_tmp[i] = im_plus_1[i] - pre_origin[i];
  }
  // rotation
  vector<double> im_new(3,0);
  vector<double> im_plus_1_new(3,0);
  rotateVector(rotation_matrix,im_tmp,im_new);
  rotateVector(rotation_matrix,im_plus_1_tmp,im_plus_1_new);
  
  vector<double> dratios(3,0);
  switch(orientation) {
    case 1:
      computeDirectionRatios(im_plus_1_new,im_new,dratios);
      cartesian2unitspherical(dratios,unit_mean);
      computeDirectionRatios(om_plus_1_new,om_new,dratios);
      cartesian2unitspherical(dratios,x);
      break;

    case 2:
      computeDirectionRatios(im_plus_1_new,om_new,dratios);
      cartesian2unitspherical(dratios,unit_mean);
      computeDirectionRatios(om_plus_1_new,om_new,dratios);
      cartesian2unitspherical(dratios,x);
      break;

    case 3:
      computeDirectionRatios(im_plus_1_new,im_new,dratios);
      cartesian2unitspherical(dratios,unit_mean);
      computeDirectionRatios(om_plus_1_new,im_new,dratios);
      cartesian2unitspherical(dratios,x);
      break;
  }
}

