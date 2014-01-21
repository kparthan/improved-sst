#include "Segment.h"
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
  double msglen = 0;
  int begin_loop = start;
  
  // state the length of segment
  msglen += encodeUsingLogStarModel(num_residues);
  if (start == 0) {
    // if segment begins with the first point in the protein, the receiver
    // has to wait until the first point(origin) and the next two points are
    // transmitted, using the sphere model
    double constant = log2(4*PI) - 2*log2(AOM);
    msglen += 2 * (log2(radii[0]) + log2(radii[1]));
    msglen += 2 * constant;
    msglen += encodeUsingNormalModel(radii);
    begin_loop = 2;
  }
  // collect the remaining radii and the directional angles
  vector<double> distances;
  vector<array<double,2>> points;
  for (int i=begin_loop; i<end; i++) {
    double r = spherical[i-2][0];
    distances.push_back(r);
    array<double,2> x;
    x[0] = spherical[i-2][1];
    x[1] = spherical[i-2][2];
    points.push_back(x);
  }
  // state the radii of remaining points in the segment
  msglen += encodeUsingNormalModel(distances);
  // state the directions
  msglen += encodeUsingMixtureModel(points,mixture);

  string name = "NullModel";
  IdealModel null_model(num_residues,name);
  return OptimalFit(null_model,msglen);
}

/*!
 *  \brief This function fits a null model to the protein segment.
 *  \param model a reference to a IdealModel;
 *  \param mixture a reference to a Mixture
 *  \return the optimal fit using the ideal model
 */
OptimalFit Segment::fitIdealModel(IdealModel &model, Mixture &mixture)
{
  // get the ideal model of length equal to the segment length
  vector<vector<double>> ideal_residues = model.getResidues(num_residues);

  double msglen = 0;
  // state the length of segment
  msglen += encodeUsingLogStarModel(num_residues);

  if (start == 0) { // the first segment
    // if segment begins with the first point in the protein, the receiver
    // has to wait until the first point(origin) and the next two points are
    // transmitted, using the sphere model
    double constant = log2(4*PI) - 2*log2(AOM);
    msglen += 2 * (log2(radii[0]) + log2(radii[1]));
    msglen += 2 * constant;
    msglen += encodeUsingNormalModel(radii);
  } else {  // an intermediate segment
    // the start point of the segment is the last point of the previous segment
    // the next two points in the segment are stated using the null model
    vector<double> distances;
    vector<array<double,2>> points;
    for (int i=start; i<start+2; i++) {
      double r = spherical[i-2][0];
      distances.push_back(r);
      array<double,2> x;
      x[0] = spherical[i-2][1];
      x[1] = spherical[i-2][2];
      points.push_back(x);
    }
    // state the radii of remaining points in the segment
    msglen += encodeUsingNormalModel(distances);
    // state the directions
    msglen += encodeUsingMixtureModel(points,mixture);
  }

  // [o: observed,MOVING] and [i: ideal,FIXED]
  // INITIAL SUPERPOSITION
  vector<vector<double>> observed,ideal;
  for (int i=0; i<3; i++) {
    observed.push_back(observed_residues[i]);
    ideal.push_back(ideal_residues[i]);
  }
  Superpose3DClass superpose(observed,ideal,3);
  suffStatClass stats = superpose.getSufficientStatistics();
  superpose.transformVectors(observed);
  //tranformToCanonicalForm(observed);
  vector<Point<double>> four_mer;
  /*for (int om=start+2; om<end-1; om++) {
    // start,...,om : points known to receiver
    //           om : most recent point communicated to the receiver
    //         om+1 : current point being transmitted
    // superpose (start,...,om) with the ideal model's (i1,...,im)
    cout << "here\n";
    vector<vector<double>> moving = ideal_residues;
    vector<vector<double>> fixed = observed_residues;
    Superpose3DClass superpose(moving,fixed,N);
    superpose.transformVectors(moving);
    superpose.transformVectors(fixed);
    N++; 
    exit(1);
  }*/

  ProteinStructure *p = model.getStructure();
  string name = model.getName();
  IdealModel ideal_model(p,name);
  ideal_model.setLength(num_residues);
  return OptimalFit(ideal_model,msglen);
}

/*!
 *
 */
void Segment::transform()
{
}

