# HV NLoS Localization

## Overview
This repository contains MATLAB implementations for Hidden Vehicle (HV) localization using Non-Line-of-Sight (NLoS) signals. The project implements two different approaches to estimate the position and orientation of hidden vehicles using NLoS propagation paths through scatterers.

## Project Structure

### Grid Search Approach (`Grid search/`)
- **Main Script**: `SensingHV_montecarlo_d1known_filesave_final.m`
  - Monte Carlo simulation for HV localization
  - Grid search optimization for orientation parameters (ω, Q)
  - Performance evaluation across different HV positions
  - RMSE analysis and visualization

- **Support Functions**:
  - `matA_d1known.m`: Matrix A computation for 3D vehicle sensing equations
  - `matB_d1known.m`: Matrix B computation for TDoA-based localization
  - `xyzerr.m`: Error analysis and histogram visualization

- **Output**: `result/` folder containing simulation results

### Iterative Search Approach (`Iterative search/`)
- **Main Scripts**:
  - `SensingHV_Iterative_visual_compare.m`: Iterative algorithm with visualization
  - `SensingHV_Iterative_visual_compare_graphsplit.m`: Enhanced visualization with split graphs
  - `SensingHV_Iterative_visual_HVpos.m`: HV position estimation visualization

- **Support Functions**:
  - `calvp.m`: Calculate v parameters using pseudo-inverse method
  - `matA_caliter.m`, `matB_caliter.m`: Matrix computation for iterative algorithm
  - `matA_calvp.m`, `matB_calvp.m`: Matrix computation for v parameter calculation
  - `calnoisemodel.m`: Noise model calculation
  - `calnoisev.m`: Noise variance calculation
  - `triop.m`: Triangulation operations

## Algorithm Details

### Grid Search Method
1. **Setup**: Define sensing vehicle (SV), hidden vehicle (HV), and scatterer positions
2. **Parameter Space**: Exhaustive search over orientation angles ω (-179° to 180°) and Q (-9° to 10°)
3. **Optimization**: Minimize the objective function `w_est = sum(null(A.')'*B)`
4. **Position Estimation**: Calculate HV position using optimized parameters
5. **Monte Carlo**: Repeat with noise for statistical analysis

### Iterative Search Method
1. **Initialization**: Start with initial estimates for Q₀ and ω₀
2. **Iterative Update**:
   - Solve `A*X = B` where `X = [ΔQ; Δω]`
   - Update parameters: `Q(t+1) = Q(t) + ΔQ`, `ω(t+1) = ω(t) + Δω`
3. **Convergence**: Continue until changes are below tolerance (3e-3)
4. **Position Calculation**: Estimate HV position using converged parameters

## Key Features
- **3D Localization**: Full 3D position and orientation estimation
- **NLoS Propagation**: Utilizes scattered signals for localization
- **Noise Robustness**: Handles measurement noise in angles and time delays
- **Performance Analysis**: Comprehensive error analysis and visualization
- **Comparative Study**: Two different optimization approaches

## Parameters
- **P**: Number of NLoS paths (3-4 scatterers)
- **c**: Speed of light (3×10⁸ m/s)
- **Noise Levels**:
  - Angle noise: 1° standard deviation
  - Time delay noise: 1 ns standard deviation
  - Distance noise: 0.1 m standard deviation

## Usage
1. **Grid Search**: Run `SensingHV_montecarlo_d1known_filesave_final.m`
2. **Iterative Search**: Run `SensingHV_Iterative_visual_compare.m`

## Applications
- Autonomous vehicle sensing
- Through-wall radar systems
- Indoor positioning systems
- Surveillance and security applications

## Requirements
- MATLAB with Parallel Computing Toolbox (for `parfor` loops)
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Results
The algorithms provide:
- Position accuracy within meters for typical scenarios
- Orientation estimation accuracy within degrees
- Real-time capability for iterative approach
- Comprehensive performance maps for different scenarios
