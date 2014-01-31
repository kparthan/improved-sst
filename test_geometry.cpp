#include "Geometry3D.h"

vector<double> ORIGIN(3,0);

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

void print(vector<double> &x)
{
  for (int i=0; i<3; i++) {
    cout << fixed << setw(5) << setprecision(3) << x[i] << " ";
  }
  cout << endl;
}

void print(vector<vector<double>> &x)
{
  for (int i=0; i<3; i++) {
    print(x[i]);
  }
  cout << endl;
}

main()
{
  /*print(XAXIS);
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
  //cout << "norm: " << vectorNorm(xaxis) << endl;
  vector<double*>v;
  v.push_back(y);
  cout << "v[0]: "; print(v[0]);

  vector<double> v2 = {1,2,3};
  cout << "v2: \n";
  //print(v2);
  vector<vector<double>> v3(4,vector<double>());
  v3[0] = v2;
  vector<double>v4({1,7,8}),v5;
  cout << v4[0] << " " << v4[1] << " " << v4[2] << endl;
  v5 = vector<double>({1,3,4,5});
  cout << v5.size() << endl;
  cout << v5[0] << " " << v5[1] << " " << v5[2] << " " << v5[3] << endl;*/

  vector<double> xaxis = {1,0,0};
  vector<double> tmp(3,0);
  vector<vector<double>> rm1;
  for (int i=0; i<3; i++) {
    rm1.push_back(tmp);
  }
  //print(rm1);
  computeRotationMatrix(xaxis,M_PI/2,rm1);
  print(rm1);

  vector<double> yaxis = {0,1,0};
  vector<vector<double>> rm2;
  for (int i=0; i<3; i++) {
    rm2.push_back(tmp);
  }
  computeRotationMatrix(yaxis,M_PI/3,rm2);
  print(rm2);

  vector<vector<double>> rm3;
  for (int i=0; i<3; i++) {
    rm3.push_back(tmp);
  }
  multiply(rm2,rm1,rm3);
  print(rm3);
}


