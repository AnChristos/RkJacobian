#ifndef RKJACOBIAN_H
#define RKJACOBIAN_H
#include "vec.h"

inline void
JacProp(double* __restrict__ P,
        const double* __restrict__ H0,
        const double* __restrict__ H1,
        const double* __restrict__ H2,
        const double* __restrict__ A,
        const double* __restrict__ A0,
        const double* __restrict__ A3,
        const double* __restrict__ A4,
        const double* __restrict__ A6,
        const double S3)
{
  double* d2A = &P[24];
  double* d3A = &P[31];
  double* d4A = &P[38];
  double d2A0 = H0[2] * d2A[1] - H0[1] * d2A[2];
  double d2B0 = H0[0] * d2A[2] - H0[2] * d2A[0];
  double d2C0 = H0[1] * d2A[0] - H0[0] * d2A[1];
  double d3A0 = H0[2] * d3A[1] - H0[1] * d3A[2];
  double d3B0 = H0[0] * d3A[2] - H0[2] * d3A[0];
  double d3C0 = H0[1] * d3A[0] - H0[0] * d3A[1];
  double d4A0 = (A0[0] + H0[2] * d4A[1]) - H0[1] * d4A[2];
  double d4B0 = (A0[1] + H0[0] * d4A[2]) - H0[2] * d4A[0];
  double d4C0 = (A0[2] + H0[1] * d4A[0]) - H0[0] * d4A[1];
  double d2A2 = d2A0 + d2A[0];
  double d2B2 = d2B0 + d2A[1];
  double d2C2 = d2C0 + d2A[2];
  double d3A2 = d3A0 + d3A[0];
  double d3B2 = d3B0 + d3A[1];
  double d3C2 = d3C0 + d3A[2];
  double d4A2 = d4A0 + d4A[0];
  double d4B2 = d4B0 + d4A[1];
  double d4C2 = d4C0 + d4A[2];
  double d0 = d4A[0] - A[0];
  double d1 = d4A[1] - A[1];
  double d2 = d4A[2] - A[2];
  double d2A3 = (d2A[0] + d2B2 * H1[2]) - d2C2 * H1[1];
  double d2B3 = (d2A[1] + d2C2 * H1[0]) - d2A2 * H1[2];
  double d2C3 = (d2A[2] + d2A2 * H1[1]) - d2B2 * H1[0];
  double d3A3 = (d3A[0] + d3B2 * H1[2]) - d3C2 * H1[1];
  double d3B3 = (d3A[1] + d3C2 * H1[0]) - d3A2 * H1[2];
  double d3C3 = (d3A[2] + d3A2 * H1[1]) - d3B2 * H1[0];
  double d4A3 = ((A3[0] + d0) + d4B2 * H1[2]) - d4C2 * H1[1];
  double d4B3 = ((A3[1] + d1) + d4C2 * H1[0]) - d4A2 * H1[2];
  double d4C3 = ((A3[2] + d2) + d4A2 * H1[1]) - d4B2 * H1[0];
  double d2A4 = (d2A[0] + d2B3 * H1[2]) - d2C3 * H1[1];
  double d2B4 = (d2A[1] + d2C3 * H1[0]) - d2A3 * H1[2];
  double d2C4 = (d2A[2] + d2A3 * H1[1]) - d2B3 * H1[0];
  double d3A4 = (d3A[0] + d3B3 * H1[2]) - d3C3 * H1[1];
  double d3B4 = (d3A[1] + d3C3 * H1[0]) - d3A3 * H1[2];
  double d3C4 = (d3A[2] + d3A3 * H1[1]) - d3B3 * H1[0];
  double d4A4 = ((A4[0] + d0) + d4B3 * H1[2]) - d4C3 * H1[1];
  double d4B4 = ((A4[1] + d1) + d4C3 * H1[0]) - d4A3 * H1[2];
  double d4C4 = ((A4[2] + d2) + d4A3 * H1[1]) - d4B3 * H1[0];
  double d2A5 = 2. * d2A4 - d2A[0];
  double d2B5 = 2. * d2B4 - d2A[1];
  double d2C5 = 2. * d2C4 - d2A[2];
  double d3A5 = 2. * d3A4 - d3A[0];
  double d3B5 = 2. * d3B4 - d3A[1];
  double d3C5 = 2. * d3C4 - d3A[2];
  double d4A5 = 2. * d4A4 - d4A[0];
  double d4B5 = 2. * d4B4 - d4A[1];
  double d4C5 = 2. * d4C4 - d4A[2];
  double d2A6 = d2B5 * H2[2] - d2C5 * H2[1];
  double d2B6 = d2C5 * H2[0] - d2A5 * H2[2];
  double d2C6 = d2A5 * H2[1] - d2B5 * H2[0];
  double d3A6 = d3B5 * H2[2] - d3C5 * H2[1];
  double d3B6 = d3C5 * H2[0] - d3A5 * H2[2];
  double d3C6 = d3A5 * H2[1] - d3B5 * H2[0];
  double d4A6 = d4B5 * H2[2] - d4C5 * H2[1];
  double d4B6 = d4C5 * H2[0] - d4A5 * H2[2];
  double d4C6 = d4A5 * H2[1] - d4B5 * H2[0];

  //--->
  double* dR = &P[21];
  dR[0] += (d2A2 + d2A3 + d2A4) * S3;
  dR[1] += (d2B2 + d2B3 + d2B4) * S3;
  dR[2] += (d2C2 + d2C3 + d2C4) * S3;
  d2A[0] = ((d2A0 + 2. * d2A3) + (d2A5 + d2A6)) * (1. / 3.);
  d2A[1] = ((d2B0 + 2. * d2B3) + (d2B5 + d2B6)) * (1. / 3.);
  d2A[2] = ((d2C0 + 2. * d2C3) + (d2C5 + d2C6)) * (1. / 3.);

  dR = &P[28];
  dR[0] += (d3A2 + d3A3 + d3A4) * S3;
  dR[1] += (d3B2 + d3B3 + d3B4) * S3;
  dR[2] += (d3C2 + d3C3 + d3C4) * S3;
  d3A[0] = ((d3A0 + 2. * d3A3) + (d3A5 + d3A6)) * (1. / 3.);
  d3A[1] = ((d3B0 + 2. * d3B3) + (d3B5 + d3B6)) * (1. / 3.);
  d3A[2] = ((d3C0 + 2. * d3C3) + (d3C5 + d3C6)) * (1. / 3.);

  dR = &P[35];
  dR[0] += (d4A2 + d4A3 + d4A4) * S3;
  dR[1] += (d4B2 + d4B3 + d4B4) * S3;
  dR[2] += (d4C2 + d4C3 + d4C4) * S3;
  d4A[0] = ((d4A0 + 2. * d4A3) + (d4A5 + d4A6 + A6[0])) * (1. / 3.);
  d4A[1] = ((d4B0 + 2. * d4B3) + (d4B5 + d4B6 + A6[1])) * (1. / 3.);
  d4A[2] = ((d4C0 + 2. * d4C3) + (d4C5 + d4C6 + A6[2])) * (1. / 3.);
}
inline void
JacPropLoop(double* __restrict__ P,
            const double* __restrict__ H0,
            const double* __restrict__ H1,
            const double* __restrict__ H2,
            const double* __restrict__ A,
            const double* __restrict__ A0,
            const double* __restrict__ A3,
            const double* __restrict__ A4,
            const double* __restrict__ A6,
            const double S3)

{

  for (int i = 21; i < 42; i += 7) {
    double* dR = &P[i];
    double* dA = &P[i + 3];
    double dA0 = H0[2] * dA[1] - H0[1] * dA[2];
    double dB0 = H0[0] * dA[2] - H0[2] * dA[0];
    double dC0 = H0[1] * dA[0] - H0[0] * dA[1];
    if (i == 35) {
      dA0 += A0[0];
      dB0 += A0[1];
      dC0 += A0[2];
    }
    double dA2 = dA0 + dA[0];
    double dB2 = dB0 + dA[1];
    double dC2 = dC0 + dA[2];
    double dA3 = dA[0] + dB2 * H1[2] - dC2 * H1[1];
    double dB3 = dA[1] + dC2 * H1[0] - dA2 * H1[2];
    double dC3 = dA[2] + dA2 * H1[1] - dB2 * H1[0];
    if (i == 35) {
      dA3 += A3[0] - A[0];
      dB3 += A3[1] - A[1];
      dC3 += A3[2] - A[2];
    }
    double dA4 = dA[0] + dB3 * H1[2] - dC3 * H1[1];
    double dB4 = dA[1] + dC3 * H1[0] - dA3 * H1[2];
    double dC4 = dA[2] + dA3 * H1[1] - dB3 * H1[0];
    if (i == 35) {
      dA4 += A4[0] - A[0];
      dB4 += A4[1] - A[1];
      dC4 += A4[2] - A[2];
    }
    double dA5 = dA4 + dA4 - dA[0];
    double dB5 = dB4 + dB4 - dA[1];
    double dC5 = dC4 + dC4 - dA[2];
    double dA6 = dB5 * H2[2] - dC5 * H2[1];
    double dB6 = dC5 * H2[0] - dA5 * H2[2];
    double dC6 = dA5 * H2[1] - dB5 * H2[0];
    if (i == 35) {
      dA6 += A6[0];
      dB6 += A6[1];
      dC6 += A6[2];
    }
    dR[0] += (dA2 + dA3 + dA4) * S3;
    dA[0] = (dA0 + dA3 + dA3 + dA5 + dA6) * (1. / 3.);
    dR[1] += (dB2 + dB3 + dB4) * S3;
    dA[1] = (dB0 + dB3 + dB3 + dB5 + dB6) * (1. / 3.);
    dR[2] += (dC2 + dC3 + dC4) * S3;
    dA[2] = (dC0 + dC3 + dC3 + dC5 + dC6) * (1. / 3.);
  }
}

inline void
JacPropVec(double* __restrict__ P,
           const double* __restrict__ H0,
           const double* __restrict__ H1,
           const double* __restrict__ H2,
           const double* __restrict__ A,
           const double* __restrict__ A0,
           const double* __restrict__ A3,
           const double* __restrict__ A4,
           const double* __restrict__ A6,
           const double S3)
{

  using namespace CxxUtils;
  using vec2 = CxxUtils::vec<double, 2>;

  /***
   * d step 2  PART
   */
  vec2 d2R_01{};
  vload(d2R_01, &P[21]);
  vec2 d2R_2{ P[23] };
  vec2 d2_01;
  vload(d2_01, &P[24]);
  vec2 d2_12;
  vload(d2_12, &P[25]);
  vec2 d2_20{};
  vblend<1, 2>(d2_20, d2_12, d2_01);

  // Magnetic field
  // H0
  vec2 H0_12;
  vload(H0_12, &H0[1]);
  vec2 H0_0{ H0[0], 0 };
  vec2 H0_20{};
  vblend<1, 2>(H0_20, H0_12, H0_0);
  // double d2A0 = H0[2]*d2A[1]-H0[1]*d2A[2];
  // double d2B0 = H0[0]*d2A[2]-H0[2]*d2A[0];
  // double d2C0 = H0[1]*d2A[0]-H0[0]*d2A[1];
  vec2 d20_01 = H0_20 * d2_12 - H0_12 * d2_20;
  vec2 d20_2 = H0_12 * d2_01 - H0_0 * d2_12;
  // double d2A2 = d2A0+d2A[0];
  // double d2B2 = d2B0+d2A[1];
  // double d2C2 = d2C0+d2A[2];
  vec2 d22_01 = d20_01 + d2_01;
  vec2 d22_2 = d20_2 + d2_20;
  vec2 d22_12;
  vblend<1, 2>(d22_12, d22_01, d22_2);
  vec2 d22_20;
  vblend<0, 2>(d22_20, d22_2, d22_01);
  // double d2A3 = (d2A[0] + d2B2 * H1[2]) - d2C2 * H1[1];
  // double d2B3 = (d2A[1] + d2C2 * H1[0]) - d2A2 * H1[2];
  // double d2C3 = (d2A[2] + d2A2 * H1[1]) - d2B2 * H1[0];
  // H1
  vec2 H1_12;
  vload(H1_12, &H1[1]);
  vec2 H1_0{ H1[0], 0 };
  vec2 H1_20{};
  vblend<1, 2>(H1_20, H1_12, H1_0);

  vec2 d23_01 = (d2_01 + d22_12 * H1_20) - d22_20 * H1_12;
  vec2 d23_2 = (d2_20 + d22_01 * H1_12) - d22_12 * H1_0;
  vec2 d23_12;
  vblend<1, 2>(d23_12, d23_01, d23_2);
  vec2 d23_20;
  vblend<0, 2>(d23_20, d23_2, d23_01);
  // double d2A4 = (d2A[0] + d2B3 * H1[2]) - d2C3 * H1[1];
  // double d2B4 = (d2A[1] + d2C3 * H1[0]) - d2A3 * H1[2];
  // double d2C4 = (d2A[2] + d2A3 * H1[1]) - d2B3 * H1[0];
  vec2 d24_01 = (d2_01 + d23_12 * H1_20) - d23_20 * H1_12;
  vec2 d24_2 = (d2_20 + d23_01 * H1_12) - d23_12 * H1_0;
  // double d2A5 = 2. * d2A4 - d2A[0];
  // double d2B5 = 2. * d2B4 - d2A[1];
  // double d2C5 = 2. * d2C4 - d2A[2];
  vec2 d25_01 = 2 * d24_01 - d2_01;
  vec2 d25_2 = 2 * d24_2 - d2_20;
  vec2 d25_12;
  vblend<1, 2>(d25_12, d25_01, d25_2);
  vec2 d25_20;
  vblend<0, 2>(d25_20, d25_2, d25_01);
  // double d2A6 = d2B5 * H2[2] - d2C5 * H2[1];
  // double d2B6 = d2C5 * H2[0] - d2A5 * H2[2];
  // double d2C6 = d2A5 * H2[1] - d2B5 * H2[0];
  // H2
  vec2 H2_12;
  vload(H2_12, &H2[1]);
  vec2 H2_0{ H2[0], 0 };
  vec2 H2_20{};
  vblend<1, 2>(H2_20, H2_12, H2_0);
  vec2 d26_01 = d25_12 * H2_20 - d25_20 * H2_12;
  vec2 d26_2 = d25_01 * H2_12 - d25_12 * H2_0;
  // dR[0] += (d2A2 + d2A3 + d2A4) * S3;
  // dR[1] += (d2B2 + d2B3 + d2B4) * S3;
  // dR[2] += (d2C2 + d2C3 + d2C4) * S3;
  d2R_01 += (d22_01 + d23_01 + d24_01) * S3;
  d2R_2 += ((d22_2 + d23_2 + d24_2) * S3);
  /***
   * d step 3  PART
   * Same as d2 in principle
   */
  vec2 d3R_01{};
  vload(d3R_01, &P[28]);
  vec2 d3R_2{ P[30] };

  vec2 d3_01;
  vload(d3_01, &P[31]);
  vec2 d3_12;
  vload(d3_12, &P[32]);

  vec2 d3_20{};
  vblend<1, 2>(d3_20, d3_12, d3_01);
  // double d3A0 = H0[2]*d3A[1]-H0[1]*d3A[2];
  // double d3B0 = H0[0]*d3A[2]-H0[2]*d3A[0];
  // double d3C0 = H0[1]*d3A[0]-H0[0]*d3A[1];
  vec2 d30_01 = H0_20 * d3_12 - H0_12 * d3_20;
  vec2 d30_2 = H0_12 * d3_01 - H0_0 * d3_12;
  // double d3A2 = d3A0+d3A[0];
  // double d3B2 = d3B0+d3A[1];
  // double d3C2 = d3C0+d3A[2];
  vec2 d32_01 = d30_01 + d3_01;
  vec2 d32_2 = d30_2 + d3_20;
  vec2 d32_12;
  vblend<1, 2>(d32_12, d32_01, d32_2);
  vec2 d32_20;
  vblend<0, 2>(d32_20, d32_2, d32_01);
  // double d3A3 = (d3A[0] + d3B2 * H1[2]) - d3C2 * H1[1];
  // double d3B3 = (d3A[1] + d3C2 * H1[0]) - d3A2 * H1[2];
  // double d3C3 = (d3A[2] + d3A2 * H1[1]) - d3B2 * H1[0];
  vec2 d33_01 = (d3_01 + d32_12 * H1_20) - d32_20 * H1_12;
  vec2 d33_2 = (d3_20 + d32_01 * H1_12) - d32_12 * H1_0;
  vec2 d33_12;
  vblend<1, 2>(d33_12, d33_01, d33_2);
  vec2 d33_20;
  vblend<0, 2>(d33_20, d33_2, d33_01);
  // double d3A4 = (d3A[0] + d3B3 * H1[2]) - d3C3 * H1[1];
  // double d3B4 = (d3A[1] + d3C3 * H1[0]) - d3A3 * H1[2];
  // double d3C4 = (d3A[2] + d3A3 * H1[1]) - d3B3 * H1[0];
  vec2 d34_01 = (d3_01 + d33_12 * H1_20) - d33_20 * H1_12;
  vec2 d34_2 = (d3_20 + d33_01 * H1_12) - d33_12 * H1_0;
  // double d3A5 = 2. * d3A4 - d3A[0];
  // double d3B5 = 2. * d3B4 - d3A[1];
  // double d3C5 = 2. * d3C4 - d3A[2];
  vec2 d35_01 = 2 * d34_01 - d3_01;
  vec2 d35_2 = 2 * d34_2 - d3_20;
  vec2 d35_12;
  vblend<1, 2>(d35_12, d35_01, d35_2);
  vec2 d35_20;
  vblend<0, 2>(d35_20, d35_2, d35_01);
  // double d3A6 = d3B5 * H2[2] - d3C5 * H2[1];
  // double d3B6 = d3C5 * H2[0] - d3A5 * H2[2];
  // double d3C6 = d3A5 * H2[1] - d3B5 * H2[0];
  vec2 d36_01 = d35_12 * H2_20 - d35_20 * H2_12;
  vec2 d36_2 = d35_01 * H2_12 - d35_12 * H2_0;
  // dR[0] += (d3A2 + d3A3 + d3A4) * S3;
  // dR[1] += (d3B2 + d3B3 + d3B4) * S3;
  // dR[2] += (d3C2 + d3C3 + d3C4) * S3;
  d3R_01 += (d32_01 + d33_01 + d34_01) * S3;
  d3R_2 += ((d32_2 + d33_2 + d34_2) * S3);
  /***
   * d step 4  PART
   */
  // double* d4A = &P[38];
  vec2 d4R_01{};
  vload(d4R_01, &P[35]);
  vec2 d4R_2{ P[37] };

  vec2 d4_01;
  vload(d4_01, &P[38]);
  vec2 d4_12;
  vload(d4_12, &P[39]);

  vec2 d4_20{};
  vblend<1, 2>(d4_20, d4_12, d4_01);
  // A0
  vec2 A0_01;
  vload(A0_01, &A0[0]);
  vec2 A0_2{ A0[2], 0 };
  // double d4A0 =(A0[0]+H0[2]*d4A[1])-H0[1]*d4A[2];
  // double d4B0 =(A0[1]+H0[0]*d4A[2])-H0[2]*d4A[0];
  // double d4C0 =(A0[2]+H0[1]*d4A[0])-H0[0]*d4A[1];
  vec2 d40_01 = (A0_01 + H0_20 * d4_12) - H0_12 * d4_20;
  vec2 d40_2 = ((A0_2 + H0_12 * d4_01) - H0_0 * d4_12);
  // double d4A2 = d4A0 + d4A[0];
  // double d4B2 = d4B0 + d4A[1];
  // double d4C2 = d4C0 + d4A[2];
  vec2 d42_01 = d40_01 + d4_01;
  vec2 d42_2 = d40_2 + d4_20;
  vec2 d42_12;
  vblend<1, 2>(d42_12, d42_01, d42_2);
  vec2 d42_20;
  vblend<0, 2>(d42_20, d42_2, d42_01);
  // A
  vec2 A_01;
  vload(A_01, &A[0]);
  vec2 A_2{ A[2], 0 };
  // double d0 = d4A[0] - A[0];
  // double d1 = d4A[1] - A[1];
  // double d2 = d4A[2] - A[2];
  vec2 d_01 = d4_01 - A_01;
  vec2 d_2 = d4_20 - A_2;
  // A3
  vec2 A3_01;
  vload(A3_01, &A3[0]);
  vec2 A3_2{ A3[2], 0 };
  // double d4A3 = ((A3[0] + d0) + d4B2 * H1[2]) - d4C2 * H1[1];
  // double d4B3 = ((A3[1] + d1) + d4C2 * H1[0]) - d4A2 * H1[2];
  // double d4C3 = ((A3[2] + d2) + d4A2 * H1[1]) - d4B2 * H1[0];
  vec2 d43_01 = ((A3_01 + d_01) + d42_12 * H1_20) - d42_20 * H1_12;
  vec2 d43_2 = ((A3_2 + d_2) + d42_01 * H1_12) - d42_12 * H1_0;
  vec2 d43_12;
  vblend<1, 2>(d43_12, d43_01, d43_2);
  vec2 d43_20;
  vblend<0, 2>(d43_20, d43_2, d43_01);
  // A4
  vec2 A4_01;
  vload(A4_01, &A4[0]);
  vec2 A4_2{ A4[2], 0 };
  // double d4A4 = ((A4[0] + d0) + d4B3 * H1[2]) - d4C3 * H1[1];
  // double d4B4 = ((A4[1] + d1) + d4C3 * H1[0]) - d4A3 * H1[2];
  // double d4C4 = ((A4[2] + d2) + d4A3 * H1[1]) - d4B3 * H1[0];
  vec2 d44_01 = ((A4_01 + d_01) + d43_12 * H1_20) - d43_20 * H1_12;
  vec2 d44_2 = ((A4_2 + d_2) + d43_01 * H1_12) - d43_12 * H1_0;
  // double d4A5 = 2. * d4A4 - d4A[0];
  // double d4B5 = 2. * d4B4 - d4A[1];
  // double d4C5 = 2. * d4C4 - d4A[2];
  vec2 d45_01 = 2 * d44_01 - d4_01;
  vec2 d45_2 = 2 * d44_2 - d4_20;
  vec2 d45_12;
  vblend<1, 2>(d45_12, d45_01, d45_2);
  vec2 d45_20;
  vblend<0, 2>(d45_20, d45_2, d45_01);
  // double d4A6 = d4B5 * H2[2] - d4C5 * H2[1];
  // double d4B6 = d4C5 * H2[0] - d4A5 * H2[2];
  // double d4C6 = d4A5 * H2[1] - d4B5 * H2[0];
  vec2 d46_01 = d45_12 * H2_20 - d45_20 * H2_12;
  vec2 d46_2 = d45_01 * H2_12 - d45_12 * H2_0;
  // A6
  vec2 A6_01;
  vload(A6_01, &A6[0]);
  vec2 A6_2{ A6[2], 0 };
  // dR[0] += (d4A2 + d4A3 + d4A4) * S3;
  // dR[1] += (d4B2 + d4B3 + d4B4) * S3;
  // dR[2] += (d4C2 + d4C3 + d4C4) * S3;
  d4R_01 += (d42_01 + d43_01 + d44_01) * S3;
  d4R_2 += ((d42_2 + d43_2 + d44_2) * S3);

  vstore(&P[21], d2R_01);
  P[23] = d2R_2[0];
  // d2A[0] = ((d2A0 + 2. * d2A3) + (d2A5 + d2A6)) * (1. / 3.);
  // d2A[1] = ((d2B0 + 2. * d2B3) + (d2B5 + d2B6)) * (1. / 3.);
  // d2A[2] = ((d2C0 + 2. * d2C3) + (d2C5 + d2C6)) * (1. / 3.);
  vstore(&P[24], ((d20_01 + 2 * d23_01) + (d25_01 + d26_01)) * (1. / 3.));
  P[26] = (((d20_2 + 2 * d23_2) + (d25_2 + d26_2)) * (1. / 3.))[0];

  vstore(&P[28], d3R_01);
  P[30] = d3R_2[0];
  // d3A[0] = ((d3A0 + 2. * d3A3) + (d3A5 + d3A6)) * (1. / 3.);
  // d3A[1] = ((d3B0 + 2. * d3B3) + (d3B5 + d3B6)) * (1. / 3.);
  // d3A[2] = ((d3C0 + 2. * d3C3) + (d3C5 + d3C6)) * (1. / 3.);
  vstore(&P[31], ((d30_01 + 2 * d33_01) + (d35_01 + d36_01)) * (1. / 3.));
  P[33] = (((d30_2 + 2 * d33_2) + (d35_2 + d36_2)) * (1. / 3.))[0];

  vstore(&P[35], d4R_01);
  P[37] = d4R_2[0];
  // d4A[0] = ((d4A0 + 2. * d4A3) + (d4A5 + d4A6 + A6[0])) * (1. / 3.);
  // d4A[1] = ((d4B0 + 2. * d4B3) + (d4B5 + d4B6 + A6[1])) * (1. / 3.);
  // d4A[2] = ((d4C0 + 2. * d4C3) + (d4C5 + d4C6 + A6[2])) * (1. / 3.);
  vstore(&P[38],
         ((d40_01 + 2 * d43_01) + (d45_01 + d46_01 + A6_01)) * (1. / 3.));
  P[40] = (((d40_2 + 2 * d43_2) + (d45_2 + d46_2 + A6_2)) * (1. / 3.))[0];
}
#endif
