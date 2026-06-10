param(
    [string]$Config = "Debug",
    [string]$OutputRoot = "Visualization/multiparty_nbn_cmaes_1e6",
    [int]$Iterations = 100000000,
    [int]$PartyPop = 180,
    [int]$CoordPop = 220,
    [int]$NbnCap = 1800,
    [int]$MaxFE = 1000000
)

$ErrorActionPreference = "Stop"
$Repo = Split-Path -Parent $MyInvocation.MyCommand.Path
$CMake = "D:\software\cmake\bin\cmake.exe"
$BuildDir = Join-Path $Repo "build\ofec_mpmmo_cmake"
$SourceDir = Join-Path $Repo "Sourcecode\freepeak_ofec\dll"
$Exe = Join-Path $Repo "build\freepeaks_benchmark\bin\FreePeaks_Benchmark.exe"

& $CMake -S $SourceDir -B $BuildDir -G "Visual Studio 17 2022" -A x64
& $CMake --build $BuildDir --config $Config --target FreePeaks_Benchmark --parallel

if (!(Test-Path $Exe)) {
    throw "Executable not found: $Exe"
}

Push-Location $Repo
try {
    & $Exe --run-multiparty-nbn-cmaes-all-dims $OutputRoot $Iterations $PartyPop $CoordPop $NbnCap $MaxFE
}
finally {
    Pop-Location
}
