#include "RkJacobian.h"
#include "benchmark/benchmark.h"
#include <algorithm>
#include <cstring>
#include <random>
/*
 * A bit hacky way to create random inputs
 */
double P[45];
double Pvec[45];
double Pvec2[45];
double H0[3];
double H1[3];
double H2[3];
double A[3] ;
double A0[3];
double A3[3];
double A4[3];
double A6[3];

class InitArray
{
public:
  InitArray()
  {
    std::mt19937 gen;
    std::uniform_real_distribution<> dis(1.0, 10.0);
    for (size_t i = 0; i < 45; ++i) {
      double in= dis(gen);
      P[i]=in;
      Pvec[i]=in;
      Pvec2[i]=in;
    }
    for (size_t i = 0; i < 3; ++i) {
      H0[i] = dis(gen);
      H1[i] = dis(gen);
      H2[i] = dis(gen);
      A[i] = dis(gen);
      A0[i] = dis(gen);
      A3[i] = dis(gen);
      A4[i] = dis(gen);
      A6[i] = dis(gen);
    }
  }
};
InitArray initArray;

static void
RkJacobian_bench(benchmark::State& state)
{
  for (auto _ : state) {
    const int n = state.range(0);
    for (int i = 0; i < n; ++i) {
      JacProp(P, H0, H1, H2, A, A0, A3, A4, A6, 1.4);
    }
  }
}
BENCHMARK(RkJacobian_bench)->RangeMultiplier(2)->Range(1024, 4096);

static void
RkJacobianVec_bench(benchmark::State& state)
{
  for (auto _ : state) {
    const int n = state.range(0);
    for (int i = 0; i < n; ++i) {
      JacPropVec(Pvec, H0, H1, H2, A, A0, A3, A4, A6, 1.4);
    }
  }
}
BENCHMARK(RkJacobianVec_bench)->RangeMultiplier(2)->Range(1024, 4096);




BENCHMARK_MAIN();
