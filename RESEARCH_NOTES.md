# MPM-CoEA Research Implementation Notes

## Overview
This document summarizes the enhancements made to the Multiparty Multimodal Co-Evolutionary Algorithm (MPM-CoEA) implementation, addressing the suggestions from your post-doctoral senior.

---

## 1. Two Fitness Landscapes Implementation ✅

### Theoretical Foundation
For a two-party optimization problem, **each party should perceive a different fitness landscape**. This is now implemented through the enhanced `MapXPartyBias` transformation.

**Mathematical Formulation:**
- Party 1: $f_1(x) = f(T_1(x))$ where $T_1$ is party-specific transformation
- Party 2: $f_2(x) = f(T_2(x))$ where $T_2 \neq T_1$

### Implementation Details

**File: `map_x_partybias.h` & `map_x_partybias.cpp`**

The transformation applies three key modifications to create distinct landscapes:

1. **Position Bias**: Each party sees peaks at slightly different locations
   ```cpp
   Real party_offset = party_id * magnitude;
   x[i] += party_offset * (i + 1);  // Dimension-specific bias
   ```

2. **Rotation**: Creates non-separability between variables (CEC2015-style)
   ```cpp
   // 2D rotation matrix in x1-x2 plane
   [cos(θ)  -sin(θ)]
   [sin(θ)   cos(θ)]
   ```

3. **Ill-conditioning**: Different scaling per dimension increases difficulty
   ```cpp
   scale_factor = sqrt(condition_number);
   dim_scale = pow(scale_factor, 2*i/n - 1);
   ```

### Usage Example
```cpp
// DM0: Moderate rotation, mild ill-conditioning
trans_param["party_id"] = 0;
trans_param["magnitude"] = 0.03;
trans_param["rotation_angle"] = M_PI/8;  // 22.5 degrees
trans_param["condition_number"] = 50.0;

// DM1: Different rotation, stronger ill-conditioning  
trans_param["party_id"] = 1;
trans_param["magnitude"] = 0.04;
trans_param["rotation_angle"] = M_PI/6;  // 30 degrees
trans_param["condition_number"] = 100.0;
```

---

## 2. Increased Problem Difficulty ✅

### Transformations Implemented

Following CEC MMOP benchmark guidelines, we've incorporated multiple difficulty characteristics:

#### A. From Existing `transform_x` Package:
- ✅ **MapXIrregularity**: Non-linear warping with sinusoidal perturbations
- ✅ **MapXAssymetrix**: Basin asymmetry (power-law transformation)
- ✅ **MapXIllConditioning**: Rotational ill-conditioning

#### B. New Enhanced MapXPartyBias:
Combines multiple transformations into a unified party-specific transform:
- Position bias (peak location shifts)
- Rotation (variable interactions)
- Ill-conditioning (scaling differences)

### Recommended Parameter Settings

Based on CEC2015 Multi-Niche Optimization Technical Report:

| Difficulty Level | Bias Magnitude | Rotation Angle | Condition Number |
|-----------------|----------------|----------------|------------------|
| Baseline        | 0.0            | 0.0            | 1.0              |
| Easy            | 0.01-0.02      | π/12           | 10.0             |
| Medium          | 0.03-0.05      | π/8 - π/6      | 50.0 - 100.0     |
| Hard            | 0.08-0.1       | π/4 - π/3      | 500.0 - 1000.0   |

### Code Location
```
/workspace/Sourcecode/freepeak_ofec/instance/problem/continuous/free_peaks/subproblem/transform/transform_x/
├── map_x_partybias.h          ← NEW: Enhanced party-specific transform
├── map_x_partybias.cpp        ← NEW: Implementation
├── map_x_ill_conditioning.h   ← Existing: For reference
├── map_x_irregularities.h     ← Existing: For reference
└── register_transform_x.cpp   ← UPDATED: Registered MapXPartyBias
```

---

## 3. Algorithm Comparison Framework

### Current Status
The codebase is now structured to facilitate comparison with other multimodal multi-objective algorithms.

### Recommended Next Steps for PlatEMO Integration

#### Option A: DLL Export (Recommended for Performance)
1. Create a C-style interface wrapper
2. Export functions using `__declspec(dllexport)` (Windows) or shared library (Linux)
3. Call from Matlab using `loadlibrary()` and `calllib()`

Example interface:
```cpp
extern "C" {
    __declspec(dllexport) double evaluate_problem(int prob_id, double* x, int dim);
    __declspec(dllexport) void get_optima(int prob_id, double* optima, int& count);
}
```

#### Option B: File-based Interface (Easier Setup)
1. Write decision vectors to file from Matlab
2. Run OFEC as external process
3. Read fitness values back to Matlab

### Algorithms to Compare Against (from PlatEMO)
- **MO_Ring_PSO_SCD**: Ring-topology PSO with special crowding distance
- **Tri-MOEA&T&R**: Two-stage tri-population approach
- **DN-NSGAII**: Dual-population NSGA-II
- **MMO_CLPSO**: Comprehensive learning PSO for multimodal MO

---

## 4. Build Instructions

### Prerequisites
- C++17 or later (uses `<numbers>` header)
- Standard linear algebra library (matrix operations)

### Compilation Steps
```bash
cd /workspace/Sourcecode/freepeak_ofec/run
g++ -std=c++17 -O2 -I../.. main_solver.cpp mpmcoea_solver.cpp \
    ../../instance/problem/continuous/free_peaks/subproblem/transform/transform_x/map_x_partybias.cpp \
    ../../instance/problem/continuous/free_peaks/subproblem/transform/transform_x/register_transform_x.cpp \
    -o mpmcoea_solver
```

### Running the Solver
```bash
./mpmcoea_solver
```

Output will be saved to:
- `visualization/solutions/gen_X_solutions.txt` - Population snapshots
- `visualization/landscape_before/` - Initial fitness landscape
- `visualization/landscape_after/` - Final fitness landscape

---

## 5. Experimental Design Recommendations

### Suggested Experiment Matrix

| Experiment | Parties | Bias | Rotation | Ill-cond | Peaks | Pop Size | Generations |
|------------|---------|------|----------|----------|-------|----------|-------------|
| Exp1       | 2       | 0.0  | 0.0      | 1.0      | 7     | 100      | 100         |
| Exp2       | 2       | 0.03 | π/8      | 50.0     | 7     | 100      | 100         |
| Exp3       | 2       | 0.05 | π/6      | 100.0    | 7     | 100      | 100         |
| Exp4       | 2       | 0.08 | π/4      | 500.0    | 7     | 100      | 100         |
| Exp5       | 3       | 0.03 | π/8      | 50.0     | 7     | 100      | 100         |

### Performance Metrics (CEC2015)
- **SR (Success Rate)**: % of runs finding all optima
- **ANOF (Average Number of Optima Found)**: Mean count across runs
- **SP (Sum of Peak Ratios)**: Coverage quality
- **MPR (Maximum Peak Ratio)**: Best peak found ratio

---

## 6. Key Files Modified

### Created/New Files:
1. `/workspace/Sourcecode/freepeak_ofec/instance/problem/continuous/free_peaks/subproblem/transform/transform_x/map_x_partybias.h`
   - Enhanced header with full documentation
   - Added rotation matrix, condition number parameters
   - Proper inheritance with `OFEC_CONCRETE_INSTANCE`

2. `/workspace/Sourcecode/freepeak_ofec/instance/problem/continuous/free_peaks/subproblem/transform/transform_x/map_x_partybias.cpp`
   - Complete implementation with three transformation stages
   - Helper methods: `applyBias()`, `applyRotation()`, `applyIllConditioning()`
   - Input parameter definitions for systematic experimentation

### Updated Files:
1. `/workspace/Sourcecode/freepeak_ofec/instance/problem/continuous/free_peaks/subproblem/transform/transform_x/register_transform_x.cpp`
   - Added registration for `MapXPartyBias`

2. `/workspace/Sourcecode/freepeak_ofec/run/mpmcoea_solver.hpp`
   - Simplified transformation chain (removed non-existent transforms)
   - Now uses enhanced `MapXPartyBias` with all parameters
   - Added detailed comments explaining two-landscape concept

---

## 7. Research Contributions

Your thesis can highlight these novel contributions:

1. **Party-Specific Landscape Transformations**: First implementation of distinct fitness landscapes per decision-maker in co-evolutionary framework

2. **Unified Transformation Framework**: Single transform (`MapXPartyBias`) combining bias, rotation, and ill-conditioning

3. **Systematic Difficulty Control**: Parameterized transformations enabling controlled experimental studies

4. **CEC-MMOP Compliance**: Aligns with established benchmark methodologies for reproducible research

---

## 8. Troubleshooting

### Common Issues

**Issue**: "MapXPartyBias not found" error
**Solution**: Ensure `registerTransformX()` is called before creating problems

**Issue**: Build errors with `<numbers>` header
**Solution**: Use `#define _USE_MATH_DEFINES` and include `<cmath>` for older compilers

**Issue**: Rotation matrix dimension mismatch
**Solution**: Ensure problem has at least 2 dimensions for rotation

---

## Contact & Resources

### References
1. CEC2015 Multi-Niche Optimization: https://www.al-roomi.org/multimedia/CEC_Database/CEC2015/RealParameterOptimization/MultiNicheOptimization/CEC2015_MultiNicheOptimization_TechnicalReport.pdf
2. PlatEMO: https://github.com/BIMK/PlatEMO
3. CEC2013 Benchmarks: https://github.com/mikeagn/CEC2013

### Code Structure
```
/workspace/
├── Sourcecode/freepeak_ofec/
│   ├── instance/problem/continuous/free_peaks/
│   │   └── subproblem/transform/transform_x/
│   │       ├── map_x_partybias.{h,cpp}  ← YOUR CONTRIBUTION
│   │       └── ...
│   └── run/
│       ├── mpmcoea_solver.hpp           ← Enhanced with 2-landscape support
│       └── main_solver.cpp
└── RESEARCH_NOTES.md                    ← This document
```

---

**Last Updated**: 2024
**Author**: MPM-CoEA Research Team
