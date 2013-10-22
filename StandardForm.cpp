#include "StandardForm.h"
#include "Message.h"

/*!
 *  \brief This is a constructor function which is used to instantiate the
 *  object
 *  \param parameters a reference to a struct Parameters
 *  \param s a reference to a Structure
 */
StandardForm::StandardForm(struct Parameters &parameters, Structure *s) :
                           parameters(parameters), structure(s)
{
  original_coordinates = structure->getCoordinates();
  numResidues = original_coordinates.size();
  rotation = Matrix<double>::identity(3);
  for (int i=0; i<numResidues; i++) {
    vector<double> tmp(numResidues,0);
    codeLength.push_back(tmp);
  }
  if (parameters.controls.size() > 0) {
    output_file = createOutputFile(1);
    ofstream log_file(output_file.c_str());
    log_file << "Using structure file: " << parameters.file << endl;
    log_file << "Number of residues: " << numResidues << endl;
    log_file << "Doing a BEZIER segmentation ..." << endl << endl;
    log_file.close();
    for (int i=0; i<numResidues; i++) {
      vector<OptimalFit> tmp(numResidues,OptimalFit());
      optimalBezierFit.push_back(tmp);
    }
  } else {
    output_file = createOutputFile(0);
    ofstream log_file(output_file.c_str(),ios::app);
    log_file << "Using structure file: " << parameters.file << endl;
    log_file << "Number of residues: " << numResidues << endl;
    log_file << "Doing a LINEAR segmentation ..." << endl << endl;
    log_file.close();
  }
}

/*!
 *  \brief This module creates a new output file to print the
 *  segmentation results to.
 *  \param status a boolean variable (0 - linear, 1 - bezier)
 *  \return the output file path
 */
string StandardForm::createOutputFile(bool status)
{
  string filtered = extractName(parameters.file);
  string output_file;
  string current_dir = string(CURRENT_DIRECTORY);
  if (status) {
    output_file = current_dir + "experiments/segmentations/logs/";
    for (int i=0; i<parameters.controls.size(); i++) {
      output_file += boost::lexical_cast<string>(parameters.controls[i]); 
    }
    output_file += "/" + filtered + ".log";
  } else {
    output_file = current_dir + "experiments/segmentations/logs/linear_";
    output_file += filtered;
  }
  return output_file;
}

/*!
 *  \brief This module gets the number of residues in the protein structure
 *  \return number of residues
 */
int StandardForm::getNumberOfResidues(void)
{
  return numResidues; 
}

/*!
 *  \brief This module returns the coordinates at the given index
 *  \param index an integer
 *  \return the coordinates at an index
 */
array<double,3> StandardForm::getCoordinates(int index)
{
  return coordinates[index];
}

/*!
 *  \brief This module computes the minimum coordinates (x,y or z) value
 *  \param index an unsigned integer
 *  \return the minimum coordinate value
 */
double StandardForm::getMinimum(unsigned index)
{
  int size = coordinates.size();
  if (size <= 0){
    cout << "Empty list of coordinates passed ..." << endl;
    cout << "exiting ..." << endl;
    exit(1);
  }
  if (index > coordinates[0].size()){
    cout << "Index exceeds std::array<> size ..." << endl;
    cout << "exiting ..." << endl;
    exit(1);
  }
  double minimum = coordinates[0][index];
  for(int i=1; i<size; i++){
    if(coordinates[i][index] < minimum){
      minimum = coordinates[i][index];
    }
  }
  return minimum; 
}

/*!
 *  \brief This module computes the maximum coordinate (x,y or z) value
 *  \param index an unsigned integer
 *  \return the maximum coordinate value
 */
double StandardForm::getMaximum(unsigned index)
{
  int size = coordinates.size();
  if (size <= 0){
    cout << "Empty list of coordinates passed ..." << endl;
    cout << "exiting ..." << endl;
    exit(1);
  }
  if (index > coordinates[0].size()){
    cout << "Index exceeds std::array<> size ..." << endl;
    cout << "exiting ..." << endl;
    exit(1);
  }
  double maximum = coordinates[0][index];
  for(int i=1; i<size; i++){
    if(coordinates[i][index] > maximum){
      maximum = coordinates[i][index];
    }
  }
  return maximum;
}

/*!
 *  \brief This module constructs a segment between the specified indices
 *  \param i an unsigned integer
 *  \param j an unsigned integer
 *  \return a segment
 */
Segment StandardForm::getSegment(unsigned i, unsigned j)
{
  if (j < i) {
    cout << "Index of segment's end point is less than the index "
         << " of the segment's start point..." << endl;
    exit(1);
  }
  int numPoints = j - i + 1;
  vector<array<double,3>> coordinates;
  for (int k=0; k<numPoints; k++){
    coordinates.push_back(getCoordinates(i+k));
  }
  return Segment(coordinates,parameters,volume);
}

/*!
 *  \brief This module updates the list of coordinates with respect to
 *  the current configuration of the protein
 */
void StandardForm::updateCoordinates(void)
{
  coordinates = structure->getCoordinates();  
}

/*!
 *  \brief This module is used to transform the protein to the standard 
 *  canonical form. The standard form is defined such that the start point
 *  of the protein coincides with the origin, the last point of the protein
 *  lies on the X-axis, and the second point of the protein lies on the XY
 *  plane.
 */
void StandardForm::transform(void)
{
  updateCoordinates();
  if(parameters.print == PRINT_DETAIL) {
    cout << "Transforming the protein to a standard canonical form ...\n";
    writeToFile(coordinates,"experiments/before_translation");
  }

  /* translate the protein so that first point is at origin */
  translateToOrigin();
  updateCoordinates();
  if(parameters.print == PRINT_DETAIL) {
    writeToFile(coordinates,"experiments/after_translation");
  }

  /* move the last point onto the X-axis */
  rotateLastPoint();
  updateCoordinates();
  if(parameters.print == PRINT_DETAIL) {
    writeToFile(coordinates,"experiments/rotate_last_point");
  }

  /* rotate second point of the protein onto the XY plane */
  rotateSecondPoint();
  updateCoordinates();
  if(parameters.print == PRINT_DETAIL) {
    writeToFile(coordinates,"experiments/rotate_second_point");
  }
  //writeToFile(coordinates,"experiments/histograms/rotate_second_point");

  /* overall transformation matrix */
  transformationMatrix();

  if(parameters.print == PRINT_DETAIL) {
    /* validate transformation */
    structure->validateTransformation(transformation);
    cout << "Transformation to standard form done ..." << endl;
    cout << "Number of residues: " << getNumberOfResidues() << endl;
  }
}

/*!
 *  This method constructs the overall transformation matrix
 */
void StandardForm::transformationMatrix()
{
  rotation.changeDimensions(4,4);
  rotation[3][3] = 1;
  transformation = rotation * translation;
  if(parameters.print == PRINT_DETAIL) {
    ofstream log_file("experiments/transformation_matrices");
    log_file << "Transformation matrix:" << endl;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
        log_file << transformation[i][j] << " ";
      }
      log_file << endl;
    }
  }
}

/*!
 *  \brief This module shifts the whole protein so that the first coordinate
 *  of the protein is at the origin.
 */
void StandardForm::translateToOrigin()
{
  if(parameters.print == PRINT_DETAIL) {
    cout << "Translation to origin ... ";
  }
  updateCoordinates();
  double offsetx,offsety,offsetz; 
  offsetx = -coordinates[0][0];
  offsety = -coordinates[0][1];
  offsetz = -coordinates[0][2];
  translation = translationMatrix(offsetx,offsety,offsetz);
  structure->transform(translation);
  if(parameters.print == PRINT_DETAIL) {
    cout << "[OK]" << endl;
  }
}

/*!
 *  \brief This module rotates the protein so that its last point lies on
 *  the X-axis.
 */
void StandardForm::rotateLastPoint()
{
  if(parameters.print == PRINT_DETAIL) {
    cout << "Rotating protein so that last point lies on X-axis ... ";
  }
  Point<double> end(coordinates[coordinates.size()-1]);

  /* brings the last point in the protein to lie on the XY plane */
  Matrix<double> rotate = projectAndRotateLast(end); 
  structure->transform(rotate);
  rotation = rotate* rotation;

  updateCoordinates();
  end = Point<double>(coordinates[coordinates.size()-1]);

  /* brings the last point on the protein to lie on the X-axis */
  rotate = rotateInXYPlane(end);
  structure->transform(rotate);
  rotation = rotate* rotation;

  if(parameters.print == PRINT_DETAIL) {
    cout << "[OK]" << endl;
  }
}

/*!
 *  \brief This module rotates the last point in the protein onto the
 *  XY plane
 *  \param p a reference to a Point<double>
 *  \return a rotation matrix to effect that transformation
 */
Matrix<double> StandardForm::projectAndRotateLast(Point<double> &p)
{
  /* project the point onto the XZ plane first */
  Point<double> projection = projectOnXZPlane(p);
  //cout << projection.x() << " " << projection.y() << " " << projection.z() << endl;

  /* now rotate the point onto the XY plane */
  return rotateLastOntoXYPlane(projection); 
}

/*!
 *  \brief This module projects a point onto the XZ plane
 *  \param p a reference to a Point<double>
 *  \return point of projection
 */
Point<double> StandardForm::projectOnXZPlane(Point<double> &p)
{
  Point<double> origin(0,0,0);
  Point<double> unitVectorX(1,0,0);
  Point<double> unitVectorZ(0,0,1);
  Plane<double> xzPlane(origin,unitVectorX,unitVectorZ);
  return project(p,xzPlane);
}

/*!
 *  \brief This module generates a rotation matrix such that the last point
 *  lies on the XY plane
 *  \param projection a reference to a Point<double>
 *  \return a transformation matrix
 */
Matrix<double> StandardForm::rotateLastOntoXYPlane(Point<double> &projection)
{
  Point<double> origin(0,0,0);
  Line<double> xaxis(Point<double> {0,0,0},Point<double> {1,0,0});
  Line<double> projectedLine(origin,projection);
  double angleWithX = angle(xaxis,projectedLine); 
  double theta;
  Vector<double> yaxis(vector<double>{0,1,0});
  if (projection.z() > 0) {
    theta = angleWithX;
  } else {
    theta = -angleWithX;
  }
  return rotationMatrix(yaxis,theta);
}

/*!
 *  \brief This module generates a rotation matrix to bring a point in XY
 *  plane onto the X axis
 *  \param p a reference to a Point<double>
 *  \return a transformation matrix
 */
Matrix<double> StandardForm::rotateInXYPlane(Point<double> &p)
{
  Point<double> origin(0,0,0);
  Line<double> xaxis(Point<double> {0,0,0},Point<double> {1,0,0});
  Line<double> projectedLine(origin,p); 
  double angleWithX = angle(xaxis,projectedLine);
  double theta;
  if (p.y() > 0){
    theta = -angleWithX;
  } else {
    theta = angleWithX;
  }
  Vector<double> zaxis(vector<double>{0,0,1});
  return rotationMatrix(zaxis,theta);
}

/*!
 *  \brief This module transforms the protein such that the second point 
 *  lies on the XY plane
 */
void StandardForm::rotateSecondPoint()
{
  if(parameters.print == PRINT_DETAIL) {
    cout << "Rotating protein so that second point lies on XY plane ... ";
  }
  updateCoordinates();
  Point<double> second = Point<double>(coordinates[1]);

  /* brings the second point in the protein to lie on the XY plane */
  Matrix<double> rotate = projectAndRotateSecond(second); 
  structure->transform(rotate);
  rotation = rotate * rotation;

  if(parameters.print == PRINT_DETAIL) {
    cout << "[OK]" << endl;
  }
}

/*!
 *  \brief This module rotates the second point in the protein onto the
 *  XY plane
 *  \param p a reference to a Point<double>
 *  \return a rotation matrix to effect that transformation
 */
Matrix<double> StandardForm::projectAndRotateSecond(Point<double> &p)
{
  /* project the point onto the YZ plane first */
  Point<double> projection = projectOnYZPlane(p);
  //cout << projection.x() << " " << projection.y() << " " << projection.z() << endl;

  /* now rotate the point onto the XY plane */
  return rotateSecondOntoXYPlane(projection); 
}

/*!
 *  \brief This module projects a point onto the YZ plane
 *  \param p a reference to a Point<double>
 *  \return point of projection
 */
Point<double> StandardForm::projectOnYZPlane(Point<double> &p)
{
  Point<double> origin(0,0,0);
  Point<double> unitVectorY(0,1,0);
  Point<double> unitVectorZ(0,0,1);
  Plane<double> yzPlane(origin,unitVectorZ,unitVectorY);
  return project(p,yzPlane);
}

/*!
 *  \brief This module generates a rotation matrix such that the second 
 *  point lies in the XY plane
 *  \param projection a reference to a Point<double>
 *  \return a transformation matrix
 */
Matrix<double> StandardForm::rotateSecondOntoXYPlane(Point<double> &projection)
{
  Point<double> origin(0,0,0);
  Line<double> yaxis(Point<double> {0,0,0},Point<double> {0,1,0});
  Line<double> projectedLine(origin,projection);
  double angleWithY = angle(yaxis,projectedLine); 
  double theta;
  if (projection.z() > 0) {
    theta = -angleWithY;
  } else {
    theta = angleWithY;
  }
  Vector<double> xaxis(vector<double>{1,0,0});
  return rotationMatrix(xaxis,theta);
}

/*!
 *  \brief This module fits the different models to the StandardForm object
 *  \return a Segmentation
 */
Segmentation StandardForm::fitModels()
{
  /* Transform the protein structure to the standard canonical form */
  transform();

  /* Construct the bounding box */
  boundingBox(); 

  /* Null model fit */
  fitNullModel();

  /* Sphere model fit */
  //fitSphereModel();

  if (parameters.controls.size() == 0) {
    /* Linear model fit */
    return fitLinearModel();
  } else {
    /* Bezier Curve fit */
    return fitBezierCurveModel();
  }
} 

/*!
 *  \brief This module computes the volume of the bounding box for
 *  the standard protein configuration
 *  \return volume of the bounding box
 */
void StandardForm::boundingBox()
{
  updateCoordinates();
  double xmin = getMinimum(0);
  double xmax = getMaximum(0);
  double ymin = getMinimum(1);
  double ymax = getMaximum(1);
  double zmin = getMinimum(2);
  double zmax = getMaximum(2);
  volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin); 
  if(parameters.print == PRINT_DETAIL) {
    cout << "boundary values:\n";
    cout << xmin << " " << xmax << endl;
    cout << ymin << " " << ymax << endl;
    cout << zmin << " " << zmax << endl;
    cout << "bounding box volume: " << volume << endl;
  }
  ofstream log_file(output_file.c_str(),ios::app);
  log_file << "Bounding box volume: " << volume << endl;
  log_file.close();
}

/*!
 *  \brief This module computes the sphere model fit to the structure.
 */
void StandardForm::fitSphereModel(void)
{
  double msglen = 0;
  for (int i=1; i<numResidues; i++){
    Point<double> previous = Point<double>(coordinates[i-1]);
    Point<double> current = Point<double>(coordinates[i]);
    double r = distance(previous,current);
    Message msg;
    msglen += msg.encodeUsingSphereModel(r);
  }
  if(parameters.print == PRINT_DETAIL) {
    cout << "*** SPHERE FIT ***" << endl;
    cout << "Sphere Model Fit: " << msglen << " bits." << endl;
    cout << "Bits per residue: " << msglen/(numResidues-1) << endl << endl;
  }
  ofstream log_file(output_file.c_str(),ios::app);
  log_file << "Sphere Model Fit: " << msglen << " bits." << endl;
  log_file << "Bits per residue: " << msglen/(numResidues-1) << endl << endl;
  log_file.close();
}

/*!
 *  \brief This module computes the null model fit to the structure.
 */
void StandardForm::fitNullModel(void)
{
  double msglen = 0;
  Message msg;
  msglen = (numResidues-1) * msg.encodeUsingNullModel(volume,AOM);
  if(parameters.print == PRINT_DETAIL) {
    cout << "*** NULL MODEL ***" << endl;
    cout << "Null Model Fit: " << msglen << " bits." << endl;
    cout << "Bits per residue: " << msglen/(numResidues-1) << endl << endl;
  }
  ofstream log_file(output_file.c_str(),ios::app);
  log_file << "Null Model Fit: " << msglen << " bits." << endl;
  null_bpr = msglen / (numResidues-1);
  log_file << "Bits per residue: " << null_bpr << endl << endl;
  log_file.close();
}

/*!
 *  \brief This module computes the linear model fit to the structure.
 */
Segmentation StandardForm::fitLinearModel(void)
{
  if(parameters.print == PRINT_DETAIL) {
    cout << "*** LINEAR FIT ***" << endl;
  }
  /* compute the code length matrix */
  computeCodeLengthMatrix();

  /* compute the optimal segmentation using dynamic programming */
  pair<double,vector<int>> segmentation = optimalSegmentation();
  printLinearSegmentation(segmentation);
  return Segmentation();
}

/*!
 *  \brief This module fits Bezier curve to the structure
 */
Segmentation StandardForm::fitBezierCurveModel()
{
  if(parameters.print == PRINT_DETAIL) {
    cout << "*** BEZIER CURVE FIT ***" << endl;
  }
  Segmentation segmentation_profile;
  if (parameters.portion_to_fit == FIT_ENTIRE_STRUCTURE) {
    clock_t c_start = clock();
    auto t_start = high_resolution_clock::now();

    /* compute the code length matrix for the Bezier curve fit */
    computeCodeLengthMatrixBezier();

    /* compute the optimal segmentation using dynamic programming */
    pair<double,vector<int>> segmentation = optimalSegmentation();
    clock_t c_end = clock();
    auto t_end = high_resolution_clock::now();
    double cpu_time = double(c_end-c_start)/(double)(CLOCKS_PER_SEC);
    double wall_time = duration_cast<seconds>(t_end-t_start).count();
    printBezierSegmentation(segmentation,cpu_time,wall_time);
    segmentation_profile = structure->reconstruct(parameters.file,output_file,
                           codeLength,optimalBezierFit,segmentation.second,
                           parameters.control_string,transformation);
    segmentation_profile.setBitsPerResidue(null_bpr,bezier_bpr);
    //segmentation_profile.setMaximumRadius(getMaximumDistance(original_coordinates));
    vector<Point<double>> coords = structure->getCoordinatesPoints();
    segmentation_profile.setCoordinates(coords);
    segmentation_profile.setTime(cpu_time,wall_time);
  } else if (parameters.portion_to_fit == FIT_SINGLE_SEGMENT) {
    fitOneSegment();
  }
  return segmentation_profile;
}

/*
 *  \brief This function fits Bezier curves to just this segment
 */
void StandardForm::fitOneSegment()
{
  array<int,2> indexes = structure->getEndPoints(parameters.end_points);
  Segment segment = getSegment(indexes[0],indexes[1]);
  segment.estimateFreeParameters();
  OptimalFit min_fit, current_fit;
  cout << "\nFit (" << parameters.controls[0] << " intermediate control points):- ";
  min_fit = segment.fitBezierCurve(parameters.controls[0]);
  min_fit.printFitInfo();
  cout << "------------------------------------------------" << endl;
  for (int k=1; k<parameters.controls.size(); k++) {
    cout << "\nFit (" << parameters.controls[k] << " intermediate control points):- "; 
    current_fit = segment.fitBezierCurve(parameters.controls[k]);
    current_fit.printFitInfo();
    cout << "------------------------------------------------" << endl;
    if (current_fit < min_fit) {
      min_fit = current_fit;
    }
  }

  cout << "\nOPTIMAL BEZIER FIT:" << endl;
  min_fit.printFitInfo();
}

/*!
 *  \brief This module constructs the code length matrix
 */
void StandardForm::computeCodeLengthMatrix(void)
{
  int procs = omp_get_num_procs();
  omp_set_num_threads(procs);
  int j;
  //#pragma omp parallel for private(j)
  for (int i=0; i<numResidues; i++){
    //#pragma omp parallel for 
    for (j=i+1; j<numResidues; j++){
      cout << "Segment: " << i << " " << j <<endl;
      Segment segment = getSegment(i,j);
      segment.fitLinear();
      codeLength[i][j] = segment.getLinearFit(); 
    }
  }
  /*ofstream codeLengthFile("codeLengthLinear");
  codeLengthFile << "# of residues: " << numResidues << endl; 
  for (int i=0; i<numResidues; i++){
    for (int j=0; j<numResidues; j++){
      codeLengthFile << fixed << setw(9) << setprecision(3) << codeLength[i][j];
    }
    codeLengthFile << endl;
  }*/

  /*
  ofstream log_file("test_linear");
  for (int i=0; i<numResidues; i++) {
    for (int j=0; j<numResidues; j++) {
      if (j > i) {
        Segment segment = getSegment(i,j);
        segment.fitLinear();
        log_file << "[" << i << ", " << j << "]: " << segment.getLinearFit() << endl;
      }
    }
  }

  Segment segment = getSegment(4,60);
  segment.fitLinear();
  cout << "linear fit: " << segment.getLinearFit() << endl;
  */
}

/*!
 *  \brief This module constructs the code length matrix for the
 *  Bezier curve fit
 */
void StandardForm::computeCodeLengthMatrixBezier(void)
{
  int procs = omp_get_num_procs();
  omp_set_num_threads(procs);
  int i,j,k,window_size,limit;
  if (parameters.constrain_segment_length == CONSTRAIN) {
    window_size = parameters.max_segment_length;
    //#pragma omp parallel for private(j)
    //cout << "#residues: " << numResidues << endl;
    for (i=0; i<numResidues; i++) {
      if (i + window_size <= numResidues) {
        limit = i + window_size - 1;
      } else {
        limit = numResidues - 1;
      }
      //#pragma omp parallel for
      for (j=i+1; j<=limit; j++) {
        Segment segment = getSegment(i,j);
        segment.estimateFreeParameters();
        OptimalFit min_fit,current_fit;
        min_fit = segment.fitBezierCurve(parameters.controls[0]);
        for (k=1; k<parameters.controls.size(); k++) {
          current_fit = segment.fitBezierCurve(parameters.controls[k]);
          if (current_fit < min_fit) {
            min_fit = current_fit;
          }
        }
        //OptimalFit fit = segment.stateUsingCurve(min_fit);
        //min_fit.setSegmentLength(segment.length());
        optimalBezierFit[i][j] = min_fit;
        codeLength[i][j] = optimalBezierFit[i][j].getMessageLength();
        if (parameters.print == PRINT_DETAIL) {
          cout << "Segment: " << i << " " << j << " ";
          cout << optimalBezierFit[i][j].getMessageLength() << endl;
        }
      }
      if (limit != numResidues - 1) {
        for (k=limit+1; k<numResidues; k++) {
          codeLength[i][k] = LARGE_NUMBER;
        }
      }
    }
  } else {
    for (i=0; i<numResidues; i++) {
      //#pragma omp parallel for
      for (j=i+1; j<numResidues; j++) {
        Segment segment = getSegment(i,j);
        segment.estimateFreeParameters();
        OptimalFit min_fit,current_fit;
        min_fit = segment.fitBezierCurve(parameters.controls[0]);
        for (k=1; k<parameters.controls.size(); k++) {
          current_fit = segment.fitBezierCurve(parameters.controls[k]);
          if (current_fit < min_fit) {
            min_fit = current_fit;
          }
        }
        optimalBezierFit[i][j] = min_fit;
        codeLength[i][j] = optimalBezierFit[i][j].getMessageLength();
        if (parameters.print == PRINT_DETAIL) {
          cout << "Segment: " << i << " " << j << " ";
          cout << optimalBezierFit[i][j].getMessageLength() << endl;
        }
      }
    }
  }
  if(parameters.print == PRINT_DETAIL) {
    ofstream codeLengthFile("codeLengthBezier");
    codeLengthFile << "# of residues: " << numResidues << endl; 
    for (int i=0; i<numResidues; i++){
      for (int j=0; j<numResidues; j++){
        codeLengthFile << fixed << setw(9) << setprecision(3) << codeLength[i][j];
      }
      codeLengthFile << endl;
    }
  }
}

/*!
 *  \brief This module computes the optimal segmentation using
 *  dynamic programming
 *  \return the indices of the segments
 */
pair<double,vector<int>> StandardForm::optimalSegmentation(void)
{
  pair <double,vector<int>> segmentation;
  vector<double> optimal(numResidues,100000);
  vector<int> optimalIndex(numResidues,-1);

  for (int i=0; i<numResidues; i++){
    optimal[i] = codeLength[0][i];
    optimalIndex[i] = i;
    for (int j=1; j<i; j++){
      if (codeLength[j][i] + optimal[j] < optimal[i]){
        optimal[i] = codeLength[j][i] + optimal[j];
        optimalIndex[i] = j;
      }
    }
  }
  segmentation.first = optimal[numResidues-1];
  bezier_bpr = segmentation.first / (numResidues - 1);
  if(parameters.print == PRINT_DETAIL) {
    cout << "Best fit: " << segmentation.first << " bits." << endl;
    cout << "Bits per residue: " << bezier_bpr << endl;
  }
  int index = numResidues-1;
  vector<int> backtrack; 
  backtrack.push_back(numResidues-1);
  while (1){
    if (index == optimalIndex[index]){
      break;
    }
    index = optimalIndex[index];
    backtrack.push_back(index);
  }
  backtrack.push_back(0);
  vector<int> segments;
  for (int i=backtrack.size()-1; i>=0; i--){
    segments.push_back(backtrack[i]);
  }
  segmentation.second = segments;
  return segmentation;
}

/*!
 *  \brief This module prints the details of the linear segmentation
 *  \param segmentation a reference to a pair<double,vector<int>>
 */
void StandardForm::printLinearSegmentation(pair<double,
                                           vector<int>> &segmentation)
{
  ofstream log_file(output_file.c_str(),ios::app);
  int i;
  vector<int> segments = segmentation.second;
  log_file << "Linear model fit: " << segmentation.first << " bits." << endl;
  log_file << "Bits per residue: " << segmentation.first / (numResidues-1) 
           << endl << endl; 
  log_file << "# of linear segments: " << segments.size()-1 << endl;
  log_file << "Internal segmentation:" << endl;
  for (i=0; i<segments.size()-1; i++) {
    log_file << segments[i] << "-->";
  }
  log_file << segments[i] << endl << endl;
  for (i=0; i<segments.size()-1; i++) {
    int segment_start = segments[i];
    int segment_end = segments[i+1];
    log_file << "Segment # " << i+1 << endl;
    log_file << "Residue stretch: [" <<  segment_start << ", " << segment_end << "]" 
    << endl;
    log_file << "Length of the segment: " << segment_end - segment_start + 1 << endl;
    log_file << "\t\tMessage length: " << codeLength[segment_start][segment_end] 
    << endl << endl;
  } 
  log_file.close();
}

/*!
 *  \brief This module prints the segmentation details for the 
 *  Bezier curve fit.
 *  \param segmentation a reference to a pair<double,vector<int>>
 *  \param cpu_time a double
 *  \param wall_time a double
 */
void
StandardForm::printBezierSegmentation(pair<double,vector<int>> &segmentation,
                                      double cpu_time, double wall_time)
{
  ofstream log_file(output_file.c_str(),ios::app);
  int i;
  vector<int> segments = segmentation.second;
  log_file << "Bezier model fit: " << segmentation.first << " bits." << endl;
  log_file << "Bits per residue: " << segmentation.first / (numResidues-1) 
           << endl << endl; 
  log_file << "# of segments: " << segments.size()-1 << endl;
  log_file << "Internal segmentation:" << endl;
  for (i=0; i<segments.size()-1; i++) {
    log_file << segments[i] << "-->";
  }
  log_file << segments[i] << endl << endl;
  log_file << "CPU time used: " << cpu_time << " secs." << endl;
  log_file << "Wall clock time elapsed: " << wall_time << " secs." << endl << endl;
  log_file << setw(15) << "SEGMENT"
           << setw(15) << "START"
           << setw(15) << "END"
           << setw(15) << "NUM_RESIDUES"
           << setw(15) << "MESSAGE_LENGTH"
           << setw(30) << "END_POINTS"
           << setw(35) << "INTERMEDIATE_CONTROL_POINTS"
           << endl;
  log_file.close();
}

