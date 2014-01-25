#include "Geometry3D.h"

using namespace std;

void print(double x[3])
{
  //x[1] = 188.4;
  for (int i=0; i<3; i++) {
    cout << x[i] << " ";
  }
  cout << endl;
}

void print(double x[3][3])
{
  for (int i=0; i<3; i++) {
    print(x[i]);
  }
}
double XAXIS[3];

main()
{
  print(XAXIS);
  //array<double,3> x = {1,0,0};
  // testing rotation matrix
  double xaxis[3] = {1,0,1};
  cout << "xaxis: ";print(xaxis);
  double *y = (double *)malloc(3*sizeof(double));
  cout << "y1: ";print(y);
  //y[1] = 188.4;
  //double rm[3][3];
  double **rm = (double **) malloc(3*sizeof(double *));
  for (int i=0; i<3; i++) {
    rm[i] = (double *)malloc(3*sizeof(double));
  }
  rm[2][2] = 1;
  cout << "rm:\n";
  print(rm[2]);
  cout << "norm: " << vectorNorm(xaxis) << endl;
  vector<double*>v;
  v.push_back(y);
  cout << "v[0]: "; print(v[0]);

  vector<double> v2 = {1,2,3};
  cout << "v2: ";
  //print(v2);
  vector<vector<double>> v3(4,vector<double>());
  v3[0] = v2;
}


