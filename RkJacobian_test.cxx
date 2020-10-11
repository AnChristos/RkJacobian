#include "RkJacobian.h"
#include <cmath>
#include <iostream>
int
main()
{
  double P[45] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
                   1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4,
                   2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
                   3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5 };
  double H0[3] = { 0.11, 0.24, 0.35 };
  double H1[3] = { 0.21, 0.40, 0.61 };
  double H2[3] = { 0.31, 0.61, 0.91 };
  double A[3] = { 0.15, 0.24, 0.39 };
  double A0[3] = { 0.22, 0.47, 0.63 };
  double A3[3] = { 0.31, 0.66, 0.90 };
  double A4[3] = { 0.41, 0.9, 1.23 };
  double A6[3] = { 0.61, 1.29, 1.57 };
  std::cout << "========== SCALAR UNROLL ==========" << '\n';
  JacProp(P, H0, H1, H2, A, A0, A3, A4, A6, 0.3);
  std::cout << P[21] << ", " << P[22] << ", " << P[23] << '\n'
            << P[24] << ", " << P[25] << ", " << P[26] << '\n'
            << P[28] << ", " << P[29] << ", " << P[30] << '\n'
            << P[31] << ", " << P[32] << ", " << P[33] << '\n'
            << P[35] << ", " << P[36] << ", " << P[37] << '\n'
            << P[38] << ", " << P[39] << ", " << P[40] << '\n';

  std::cout <<'\n' <<"============ SCALAR LOOP ============" << '\n';
  double PLoop[45] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                       1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
                       1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
                       2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
                       3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5 };
  JacPropLoop(PLoop, H0, H1, H2, A, A0, A3, A4, A6, 0.3);
  std::cout << PLoop[21] << ", " << PLoop[22] << ", " << PLoop[23] << '\n'
            << PLoop[24] << ", " << PLoop[25] << ", " << PLoop[26] << '\n'
            << PLoop[28] << ", " << PLoop[29] << ", " << PLoop[30] << '\n'
            << PLoop[31] << ", " << PLoop[32] << ", " << PLoop[33] << '\n'
            << PLoop[35] << ", " << PLoop[36] << ", " << PLoop[37] << '\n'
            << PLoop[38] << ", " << PLoop[39] << ", " << PLoop[40] << '\n';

  for (int i = 0; i < 45; ++i) {
    if (std::abs(P[i] - PLoop[i]) > 1e-14) {
      std::cout << "element  " << i << " differs " <<  std::abs(P[i] - PLoop[i]) <<'\n';
    }
  }

  std::cout<<'\n' << "============ SSE UNROLL ============" << '\n';
  double Pvec[45] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                      1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
                      1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
                      2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
                      3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5 };
  JacPropVec(Pvec, H0, H1, H2, A, A0, A3, A4, A6, 0.3);
  std::cout << Pvec[21] << ", " << Pvec[22] << ", " << Pvec[23] << '\n'
            << Pvec[24] << ", " << Pvec[25] << ", " << Pvec[26] << '\n'
            << Pvec[28] << ", " << Pvec[29] << ", " << Pvec[30] << '\n'
            << Pvec[31] << ", " << Pvec[32] << ", " << Pvec[33] << '\n'
            << Pvec[35] << ", " << Pvec[36] << ", " << Pvec[37] << '\n'
            << Pvec[38] << ", " << Pvec[39] << ", " << Pvec[40] << '\n';

  for (int i = 0; i < 45; ++i) {
    if (std::abs(P[i] - Pvec[i]) > 1e-14) {
      std::cout << "element  " << i << " differs " << std::abs(P[i] - Pvec[i]) <<'\n';
    }
  }

  std::cout<<'\n' << "============ SSE LOOP ============" << '\n';
  double PloopVec[45] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                      1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
                      1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
                      2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
                      3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5 };
  JacPropLoopVec(PloopVec, H0, H1, H2, A, A0, A3, A4, A6, 0.3);
  std::cout << PloopVec[21] << ", " << PloopVec[22] << ", " << PloopVec[23] << '\n'
            << PloopVec[24] << ", " << PloopVec[25] << ", " << PloopVec[26] << '\n'
            << PloopVec[28] << ", " << PloopVec[29] << ", " << PloopVec[30] << '\n'
            << PloopVec[31] << ", " << PloopVec[32] << ", " << PloopVec[33] << '\n'
            << PloopVec[35] << ", " << PloopVec[36] << ", " << PloopVec[37] << '\n'
            << PloopVec[38] << ", " << PloopVec[39] << ", " << PloopVec[40] << '\n';

  for (int i = 0; i < 45; ++i) {
    if (std::abs(P[i] - PloopVec[i]) > 1e-14) {
      std::cout << "element  " << i << " differs " << std::abs(P[i] - PloopVec[i]) <<'\n';
    }
  }


  return 0;
}

