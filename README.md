# PN Junction Simulation - Finite Difference Method

A MATLAB implementation for simulating current flow in a p-n diode using finite difference methods.

## 📖 Project Description

This project numerically solves the self-consistent system of equations governing the behavior of a p-n junction diode. It implements finite difference methods to calculate:

- Electric potential distribution V(x)
- Carrier concentrations n(x) and p(x)  
- Charge density ρ(x)
- Current density J(x)
- I-V characteristics of the diode

The simulation handles both equilibrium conditions and applied biases (forward and reverse), providing insights into semiconductor device physics beyond simplified analytical models.

## 🧮 Physical Models

The implementation solves the following key equations:

- **Poisson's Equation**: ΔV = -ρ(x)/εₛ
- **Carrier Statistics**: n(x) = n₀exp(qV/kₚT), p(x) = p₀exp(-qV/kₚT)
- **Charge Density**: ρ(x) = p(x) + n(x) + N_D - N_A
- **Drift-Diffusion Currents**: Includes both electron and hole contributions
- **Continuity Equation**: div(jₙ + jₚ) = 0 (stationary state)

## 🏗️ Code Structure

### Script Types

1. **Configuration Scripts** (`*_config.m`)
   - Define physical parameters and simulation settings
   - Easy adjustment without modifying core algorithms

2. **Function Scripts** (`Poisson.m`, `charge_bernoulli.m`, `Courants.m`, etc.)
   - Core numerical solvers and physics calculations
   - Modular and reusable components

3. **Execution Scripts** (`first_implementation.m`, `second_implementation.m`, `third_implementation.m`, `Caracteristique.m`)
   - Main simulation runners
   - Different implementations and analysis types

### Key Files

- `first_implementation.m` - Basic Poisson solver (diverges due to non-linearity)
- `second_implementation.m` - Newton-Raphson damped Poisson solver
- `third_implementation.m` - Full current calculation with Bernoulli functions
- `Caracteristique.m` - I-V characteristic generation
- `Poisson_NR.m` - Newton-Raphson Poisson discretization
- `charge_bernoulli.m` - Carrier concentration calculation with field effects
- `boundary_conditions.m` - Boundary condition handling

## 🚀 Getting Started

### Prerequisites

- MATLAB R2018b or later
- Basic understanding of semiconductor physics

### Usage

1. Clone the repository:
```bash
git clone https://github.com/your-username/pn-junction-simulation.git


