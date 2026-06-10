# OFEC FreePeaks-MPMMO Experiment Commands

## Build the executable with CMake

```powershell
cmake -S D:\code\multi_party_multi_modal\hamaza_code `
  -B D:\code\multi_party_multi_modal\hamaza_code\build\cmake_root `
  -G "Visual Studio 17 2022" -A x64

cmake --build D:\code\multi_party_multi_modal\hamaza_code\build\cmake_root `
  --config Debug --target FreePeaks_Benchmark --parallel
```

Executable:

```text
D:\code\multi_party_multi_modal\hamaza_code\build\freepeaks_benchmark\bin\FreePeaks_Benchmark.exe
```

## Run NBN-CMEAS at MaxFE = 1e5

```powershell
cd D:\code\multi_party_multi_modal\hamaza_code
.\run_nbn_cmaes_1e6.ps1 -OutputRoot Visualization/multiparty_nbn_cmaes_1e5 `
  -Iterations 100000000 -PartyPop 180 -CoordPop 220 -NbnCap 1800 -MaxFE 100000
```

## Run / resume PlatEMO parallel comparison at MaxFE = 1e5

```matlab
cd('D:\code\multi_party_multi_modal\hamaza_code');
addpath('Matlab');
T = run_ofec_mpmmo_platemo_parallel_1e6( ...
    'Runs',1, ...
    'MaxFE',100000, ...
    'NumWorkers',4, ...
    'OutputDir','Results/ofec_mpmmo_platemo_1e5', ...
    'Resume',true);
```

The script saves one `.mat` file per algorithm/suite/dimension/run under:

```text
D:\code\multi_party_multi_modal\hamaza_code\Results\ofec_mpmmo_platemo_1e5\tasks
```

It can be safely interrupted and resumed.

## Summarize results

```matlab
cd('D:\code\multi_party_multi_modal\hamaza_code');
addpath('Matlab');
[S,A] = summarize_ofec_mpmmo_1e5_results();
```

Summary CSV files:

```text
D:\code\multi_party_multi_modal\hamaza_code\Results\ofec_mpmmo_1e5_summary\summary_by_algorithm_dimension.csv
D:\code\multi_party_multi_modal\hamaza_code\Results\ofec_mpmmo_1e5_summary\summary_by_suite_dimension.csv
D:\code\multi_party_multi_modal\hamaza_code\Results\ofec_mpmmo_1e5_summary\all_metrics.csv
```
