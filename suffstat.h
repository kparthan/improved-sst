#ifndef SUFFSTAT_H
#define SUFFSTAT_H

class suffStatClass {
  public: 
  double sum_xm, sum_ym, sum_zm;
  double sum_xp, sum_yp, sum_zp;

  double sum_xmxm, sum_xmym, sum_xmzm;
  double sum_ymym, sum_ymzm;
  double sum_zmzm;
  double sum_xpxp, sum_xpyp, sum_xpzp;
  double sum_ypyp, sum_ypzp;
  double sum_zpzp;
  double sum_xmyp, sum_xmzp;
  double sum_ymxp, sum_ymzp;
  double sum_zmxp, sum_zmyp;

  // these two store the suff statistics to compute the centers of mass
  // of the two sets
  double sum_Ax, sum_Ay, sum_Az;
  double sum_Bx, sum_By, sum_Bz;

  size_t nVecs;
};

#endif
