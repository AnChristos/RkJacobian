# If you do not have already google benchmark installed

``git clone https://github.com/google/benchmark.git``

The following script from this repo will install the benchmark under ``${HOME}/.local``

``source buildGoogleBenchmark.sh``

# Build RkJacobian

``git clone https://github.com/AnChristos/RkJacobian.git``

``mkdir build; cd build``

``cmake ../RkJacobian``

``make``

#Example

``./RkJacobian_bench --benchmark_report_aggregates_only=true --benchmark_repetitions=20``
