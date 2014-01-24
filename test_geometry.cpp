#include "Geometry3D.h"

using namespace std;

void print(double x[3])
{
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
extern double XAXIS[3];

main()
{
  print(XAXIS);
  // testing rotation matrix
  /*double xaxis[3] = {1,0,1};
  print(xaxis);
  double rm[3][3];
  rm[2][2] = 1;
  print(rm);
  cout << "norm: " << vectorNorm(xaxis) << endl;*/
}
