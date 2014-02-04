#include "Segment.h"
#include "Message.h"
#include "Superpose3D.h"
#include "Geometry3D.h"

extern vector<double> ZAXIS;
extern int DEBUG;
extern ofstream debug;

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

  // state the null model
  msglen += log2(NUM_IDEAL_MODELS+1);

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
    double MSG;
    MSG = message.encodeUsingMixtureModel(unit_coordinates[i-2],mixture);
    msglen += MSG;
  }

  string name = "NullModel";
  IdealModel null_model(num_residues,name);
  return OptimalFit(null_model,msglen);
}

/*!
 *  \brief This function fits a null model to the protein segment.
 *  \param residual_mixture a reference to a Mixture
 *  \param sum_residual_weights a reference to a double
 *  \return the optimal fit using the null model
 */
OptimalFit Segment::fitNullModel(Mixture &residual_mixture, double &sum_residual_weights)
{
  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
  Message message;
  double msglen = 0;
  int begin_loop = start;
  
  // state the length of segment
  msglen += message.encodeUsingLogStarModel(num_residues);

  // state the null/coil model
  msglen += -log2(sum_residual_weights);

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
  // state the remaining points using the residual mixture 
  double r;
  for (int i=begin_loop; i<end; i++) {
    // state radius
    r = spherical_coordinates[i-2][0];
    msglen += message.encodeUsingNormalModel(r,normal);
    // state direction
    double MSG;
    MSG = message.encodeUsingMixtureModel(unit_coordinates[i-2],residual_mixture);
    msglen += MSG;
  }

  string name = "NullModel";
  IdealModel null_model(num_residues,name);
  return OptimalFit(null_model,msglen);
}

/*!
 *  \brief This function fits a ideal model to the protein segment using an adaptive
 *  encoding scheme.
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

  // state the ideal model
  msglen += log2(NUM_IDEAL_MODELS+1);

  if (start == 0) { // the first segment
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
    double kappa = 25;
    double rmsd;
    double MSG,MSG1;
    Mixture conflated_mixture;
    double r = spherical_coordinates[start][0];
    MSG1 = message.encodeUsingNormalModel(r,normal);
    msglen += MSG1;

    // INITIAL SUPERPOSITION
    for (int i=0; i<3; i++) {
      observed.push_back(observed_residues[i]);
      ideal.push_back(ideal_residues[i]);
    }
    observed_copy = observed;
    ideal_copy = ideal;
    Superpose3DClass superpose(observed_copy,ideal_copy);
    rmsd = superpose.getRMSD();
    suff_stats = superpose.getSufficientStatistics();
    observed_copy.push_back(observed_residues[3]);
    superpose.transformVectors(observed_copy);
    for (int i=0; i<4; i++) {
      four_mer[i] = observed_copy[i];
    }
    double theta;
    vector<double> diff1(3,0),diff2(3,0);
    computeDirectionRatios(four_mer[3],four_mer[2],diff1);
    computeDirectionRatios(ideal_residues[3],ideal_residues[2],diff2);
    theta = computeAngle(diff1,diff2);
    convertToCanonicalForm(four_mer,transformed_four_mer,rotation_matrix);
    getCurrentMeanAndDirection(four_mer[2],rotation_matrix,transformed_four_mer[2],
                               transformed_four_mer[3],ideal_residues[2],
                               ideal_residues[3],orientation,unit_mean,x);
    Component component(unit_mean,kappa);
    conflated_mixture = mixture.conflate(component);
    MSG = message.encodeUsingMixtureModel(x,conflated_mixture);
    //MSG = message.encodeUsingComponent(x,component);
    msglen += MSG;
    alignWithZAxis(ideal_residues[2],ideal_residues[3],zaxis_transform);
    applyIdealModelTransformation(zaxis_transform,observed_copy[2],
                                  observed_copy[3],vonmises_suffstats);
    int N = 1;

    // ADAPTIVE SUPERPOSITION
    for (int om=start+3; om<end; om++) {
      r = spherical_coordinates[om-2][0];
      MSG1 = message.encodeUsingNormalModel(r,normal);
      msglen += MSG1;
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
      Superpose3DClass superpose(suff_stats,observed_copy[om-start],
                                 ideal_copy[om-start]);
      rmsd = superpose.getRMSD();
      suff_stats = superpose.getSufficientStatistics();
      observed_copy.push_back(observed_residues[om-start+1]);
      superpose.transformVectors(observed_copy);
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
      //MSG = message.encodeUsingComponent(x,component);
      msglen += MSG;
      // update von mises suff stats
      alignWithZAxis(ideal_residues[om-start],ideal_residues[om-start+1],zaxis_transform);
      applyIdealModelTransformation(zaxis_transform,observed_copy[om-start],
                                    observed_copy[om-start+1],vonmises_suffstats);
      N++;
    }
  }
  ProteinStructure *p = model.getStructure();
  string name = model.getName();
  IdealModel ideal_model(p,name);
  ideal_model.setLength(num_residues);
  return OptimalFit(ideal_model,msglen);
}

/*!
 *  \brief This function fits a ideal model to the protein segment using a non-adaptive
 *  encoding scheme.
 *  \param model a reference to a IdealModel;
 *  \param mixture a reference to a Mixture
 *  \param assigned_component a reference to a Component
 *  \param weight a double
 *  \return the optimal fit using the ideal model
 */
OptimalFit Segment::fitIdealModel(IdealModel &model, Mixture &mixture,
                                     Component &assigned_component, double weight)
{
  Normal normal(NORMAL_MEAN,NORMAL_SIGMA);
  Message message;
  int begin_loop = start;

  double msglen = 0;
  // state the length of segment
  msglen += message.encodeUsingLogStarModel(num_residues);

  // state the ideal model
  msglen += -log2(weight);

  if (start == 0) { // the first segment
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
    for (int om=start+2; om<end; om++) {
      // state radius
      double r = spherical_coordinates[om-2][0];
      msglen += message.encodeUsingNormalModel(r,normal);
      msglen += message.encodeUsingComponent(unit_coordinates[om-2],assigned_component);
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
  vector<double> &pre_origin, // om 
  vector<vector<double>> &rotation_matrix,
  vector<double> &om_new, 
  vector<double> &om_plus_1_new, 
  vector<double> &im, 
  vector<double> &im_plus_1, 
  int orientation, 
  vector<double> &unit_mean, 
  vector<double> &x
) {
  vector<double> dratios_mean(3,0);
  vector<double> dratios_x(3,0);
  switch(orientation) {
    case 1:
      computeDirectionRatios(im_plus_1,im,dratios_mean);
      computeDirectionRatios(om_plus_1_new,om_new,dratios_x);
      break;

    case 2:
      computeDirectionRatios(im_plus_1,om_new,dratios_mean);
      computeDirectionRatios(om_plus_1_new,om_new,dratios_x);
      break;

    case 3:
      computeDirectionRatios(im_plus_1,im,dratios_mean);
      computeDirectionRatios(om_plus_1_new,im,dratios_x);
      break;
  }
  vector<double> unit_mean_tmp(3,0);
  cartesian2unitspherical(dratios_mean,unit_mean_tmp);
  rotateVector(rotation_matrix,unit_mean_tmp,unit_mean);
  cartesian2unitspherical(dratios_x,x);
}

