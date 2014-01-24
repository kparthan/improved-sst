#include "Segment.h"
#include "Message.h"
#include "Superpose3D.h"

/*!
 *  \brief This is a constructor function.
 *  \param start an integer 
 *  \param end an integer
 *  \param cartesian a reference to a vector<Point<double>>
 *  \param spherical a reference to a vector<array<double,2>>
 */
Segment::Segment(int start, int end, vector<Point<double>> &cartesian,
                 vector<array<double,3>> &spherical) : start(start), end(end),
                 cartesian(cartesian) , spherical(spherical)
{
  num_residues = end - start + 1;
  for (int i=start; i<=end; i++) {
    vector<double> point = point2vector(cartesian[i]);
    observed_residues.push_back(point);
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
    msglen += message.encodeUsingSphereModel(radii[1],normal);
    begin_loop = 2;
  }
  // state the remaining points using the mixture model
  double r;
  array<double,2> x;
  for (int i=begin_loop; i<end; i++) {
    // state radius
    double r = spherical[i-2][0];
    msglen += message.encodeUsingNormalModel(r,normal);
    // state theta.phi
    x[0] = spherical[i-2][1];
    x[1] = spherical[i-2][2];
    msglen += message.encodeUsingMixtureModel(x,mixture);
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

  double msglen = 0;
  // state the length of segment
  msglen += message.encodeUsingLogStarModel(num_residues);

  if (start == 0) { // the first segment
    // if segment begins with the first point in the protein, the receiver
    // has to wait until the first point(origin) and the next two points are
    // transmitted, using the sphere model
    msglen += message.encodeUsingSphereModel(radii[0],normal);
    msglen += message.encodeUsingSphereModel(radii[1],normal);
  } else {  // an intermediate segment
    // the start point of the segment is the last point of the previous segment
    // the next two points in the segment are stated using the null model
    double r;
    array<double,2> x;
    for (int i=start; i<start+2; i++) {
      // state radius
      double r = spherical[i-2][0];
      msglen += message.encodeUsingNormalModel(r,normal);
      // state theta.phi
      x[0] = spherical[i-2][1];
      x[1] = spherical[i-2][2];
      msglen += message.encodeUsingMixtureModel(x,mixture);
    }
  }

  // get the ideal model of length equal to the segment length
  vector<vector<double>> ideal_residues = model.getResidues(num_residues);
  // [observed,MOVING] and [ideal,FIXED]
  vector<vector<double>> observed,ideal;  // observed, ideal are running lists
                                          // of the orginal untransformed
                                          // coordinates
  vector<vector<double>> observed_copy,ideal_copy;
  suffStatClass suff_stats;
  pair<vector<Point<double>>,Matrix<double>> canonical_transformation;
  Matrix<double> canonical_transformation_matrix;
  vector<Point<double>> four_mer(4,Point<double>());
  pair<array<double,2>,array<double,2>> mu_x;
  array<double,2> mu,x;
  array<double,3> vonmises_suffstats;
  double kappa = 5;
  double density;
  Mixture conflated_mixture;

  // INITIAL SUPERPOSITION
  for (int i=0; i<3; i++) {
    observed.push_back(observed_residues[i]);
    ideal.push_back(ideal_residues[i]);
  }
  observed_copy = observed;
  ideal_copy = ideal;
  Superpose3DClass superpose(observed_copy,ideal_copy);
  suff_stats = superpose.getSufficientStatistics();
  observed_copy.push_back(observed_residues[3]);
  superpose.transformVectors(observed_copy);
  for (int i=0; i<4; i++) {
    four_mer[i] = Point<double>(observed_copy[i]);
  }
  canonical_transformation = convertToCanonicalForm(four_mer);
  canonical_transformation_matrix = canonical_transformation.second;
  mu_x = getCurrentMeanAndDirection(canonical_transformation,ideal_residues[2],
                                    ideal_residues[3],orientation);
  mu = mu_x.first;
  x = mu_x.second;
  Component component(mu,kappa,mixture.constrain_kappa);
  conflated_mixture = mixture.conflate(component);
  msglen += message.encodeUsingMixtureModel(x,conflated_mixture);
  Matrix<double> zaxis_transform = alignWithZAxis(ideal_residues[2],ideal_residues[3]);
  vonmises_suffstats = applyIdealModelTransformation(zaxis_transform,
                        observed_residues[2],observed_residues[3]);
  //vonmises_suffstats = convertToCartesian(1,x[0],x[1]);
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
    suff_stats = superpose.getSufficientStatistics();
    observed_copy.push_back(observed_residues[om-start+1]);
    superpose.transformVectors(observed_copy);
    //writeToFile(observed_copy,"observed_copy_after");
    //writeToFile(ideal_copy,"ideal_copy_after");
    // get the current four_mer
    for (int i=0; i<4; i++) {
      four_mer[i] = Point<double>(observed_copy[N+i]);
    }
    canonical_transformation = convertToCanonicalForm(four_mer);
    canonical_transformation_matrix = canonical_transformation.second;
    // apply the canonical transformation on the ideal four_mer
    // it is sufficient to transform im and i_{m+1}
    mu_x = getCurrentMeanAndDirection(canonical_transformation,
           ideal_residues[om-start],ideal_residues[om-start+1],orientation);
    if (N > 1) {
      Component adaptive_component(vonmises_suffstats,N,SET);
      adaptive_component.minimizeMessageLength();
      kappa = adaptive_component.getKappa();
    }
    mu = mu_x.first;
    x = mu_x.second;
    component = Component(mu,kappa,mixture.constrain_kappa);
    conflated_mixture = mixture.conflate(component);
    msglen += message.encodeUsingMixtureModel(x,conflated_mixture);
    // update von mises suff stats
    zaxis_transform = alignWithZAxis(ideal_residues[om-start],ideal_residues[om-start+1]);
    array<double,3> dir = applyIdealModelTransformation(zaxis_transform,
                        observed_residues[om-start],observed_residues[om-start+1]);
    N++;
    //array<double,3> tmp = convertToCartesian(1,x[0],x[1]);
    for (int i=0; i<3; i++) {
      vonmises_suffstats[i] += dir[i];
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
 *  \param canonical_transformation a reference to a pair<vector<Point<double>>,Matrix<double>>
 *  \param im a reference to a vector<double>
 *  \param im_plus_1 a reference to a vector<double>
 *  \param orientation an integer
 *  \return the pair of mean and direction
 */
pair<array<double,2>,array<double,2>> 
Segment::getCurrentMeanAndDirection
         (pair<vector<Point<double>>,Matrix<double>> &canonical_transformation,
          vector<double> &im, vector<double> &im_plus_1, int orientation)
{
  // transform im and i_{m+1}
  Point<double> im_new =
  lcb::geometry::transform<double>(im,canonical_transformation.second);
  Point<double> im_plus_1_new =
  lcb::geometry::transform<double>(im_plus_1,canonical_transformation.second);
  
  // get the transformed om and o_{m+1}
  Point<double> om_new = canonical_transformation.first[2];
  Point<double> om_plus_1_new = canonical_transformation.first[3];

  array<double,2> mu,x;
  array<double,3> difference;
  Point<double> diff;
  switch(orientation) {
    case 1:
      diff = im_plus_1_new - im_new;
      difference = convertToSpherical(diff);
      mu[0] = difference[1];
      mu[1] = difference[2];
      diff = om_plus_1_new - om_new;
      difference = convertToSpherical(diff);
      x[0] = difference[1];
      x[1] = difference[2];
      break;

    case 2:
      diff = im_plus_1_new - om_new;
      difference = convertToSpherical(diff);
      mu[0] = difference[1];
      mu[1] = difference[2];
      diff = om_plus_1_new - om_new;
      difference = convertToSpherical(diff);
      x[0] = difference[1];
      x[1] = difference[2];
      break;

    case 3:
      diff = im_plus_1_new - im_new;
      difference = convertToSpherical(diff);
      mu[0] = difference[1];
      mu[1] = difference[2];
      diff = om_plus_1_new - im_new;
      difference = convertToSpherical(diff);
      x[0] = difference[1];
      x[1] = difference[2];
      break;
  }
  return pair<array<double,2>,array<double,2>>(mu,x);
}

