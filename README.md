# neon‑md‑simulation
Molecular dynamics simulation of neon atoms with Lennard‑Jones potential

# Neon Molecular Dynamics Simulation

Molecular dynamics simulation of 100 neon atoms in a cubic box with Lennard‑Jones interactions. This project implements the Verlet integration scheme with periodic boundary conditions and velocity scaling thermostat. The simulation computes the average energy at 300 K and the isochoric heat capacity $C_V$ over the temperature range 200‑400 K.

## Project Overview

This work simulates 100 neon atoms in a cubic box of side length $d = 1000 \times 10^{-10}\ \text{m}$. The interatomic interaction is modeled by the Lennard‑Jones (LJ) potential. The key tasks are:

1. Compute the average energy $\langle E \rangle$ of the system at $T = 300\ \text{K}$.
2. Compute the isochoric heat capacity $C_V$ over the temperature range 200‑400 K using the energy fluctuation formula.

## Physical Model
### Lennard‑Jones Potential

The interaction between neon atoms is described by the Lennard‑Jones potential:

$$
U(r_{ij}) = 4\varepsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6 \right]
$$

The total potential energy of the system is:

$$
U_{\text{total}} = \sum_{i=1}^{N} \sum_{j>i}^{N} U(r_{ij})
$$

The force on particle $i$ is given by:

$$
\mathbf{F}_i = -\nabla_i U_{\text{total}} = \sum_{j \neq i} 4\varepsilon \left[ 12 \frac{\sigma^{12}}{r_{ij}^{13}} - 6 \frac{\sigma^6}{r_{ij}^7} \right] \hat{\mathbf{r}}_{ij}
$$

### Parameters

| Parameter | Value | Description |
|:---|:---|:---|
| $N$ | 100 | Number of atoms |
| $m$ | $4 \times 1.667 \times 10^{-27}\ \text{kg}$ | Mass of neon atom |
| $\varepsilon$ | $1.48 \times 10^{-22}\ \text{J}$ | LJ potential well depth |
| $\sigma$ | $1.48 \times 10^{-10}\ \text{m}$ | LJ zero‑crossing distance |
| $d$ | $1000 \times 10^{-10}\ \text{m}$ | Simulation box side length |
| $r_{\text{cut}}$ | $3\sigma$ | Cutoff radius for interactions |
| $r_{\text{safe}}$ | $0.8\sigma$ | Safe distance to prevent numerical divergence |

### Equations of Motion

Newton's equations of motion are integrated using the Verlet algorithm:

$$
\mathbf{r}_i(t+\Delta t) = 2\mathbf{r}_i(t) - \mathbf{r}_i(t-\Delta t) + \mathbf{a}_i(t)\, \Delta t^2
$$

Velocities are computed from the position difference:

$$
\mathbf{v}_i(t) = \frac{\mathbf{r}_i(t+\Delta t)-\mathbf{r}_i(t-\Delta t)}{2\Delta t}
$$

### Boundary Conditions

Periodic boundary conditions are applied in all three dimensions to eliminate surface effects and simulate a bulk environment. When an atom crosses a boundary, it re‑enters from the opposite side.

### Cutoff Radius

A cutoff radius of $r_{\text{cut}} = 3\sigma$ is applied to reduce computational cost. A safe distance of $r_{\text{safe}} = 0.8\sigma$ prevents numerical divergence when atoms approach too closely.

## Simulation Workflow

### Initialization

1. **Parameter definition**: Set all physical parameters ($N$, mass, LJ parameters, box size, time step, total steps, etc.)
2. **Initial positions**: Place atoms on a $5 \times 5 \times 4$ uniform grid with small random perturbations
3. **Random seed**: Initialize the random number generator with a fixed seed for reproducibility
4. **Initial velocities**: Generate from Maxwell‑Boltzmann distribution using Box‑Muller transform
5. **Temperature scaling**: Apply velocity scaling to match target temperature ($300\ \text{K}$)
6. **History array initialization**: Compute $\mathbf{r}(t-\Delta t)$ for Verlet integration

### Main Loop

For each time step:

1. **Update positions** using Verlet algorithm
2. **Apply periodic boundary conditions**
3. **Compute forces** from LJ potential (force subroutine)
4. **Update velocities and accelerations**
5. **Temperature control**: Strong scaling during equilibration, gentle scaling during production
6. **Remove center-of-mass drift**
7. **Compute kinetic and potential energies**
8. **Collect thermodynamic data** (energy fluctuations for heat capacity)
9. **Output results**
   
### Heat Capacity Calculation

The isochoric heat capacity is computed from energy fluctuations:

$$
C_V = \frac{\langle E^2 \rangle - \langle E \rangle^2}{k_B T^2}
$$

The Welford algorithm is implemented to compute variance numerically, avoiding precision loss from subtracting large numbers.

## Results

### Average Energy at 300 K

At $T = 300\ \text{K}$, the system total energy stabilizes at approximately $2.305 \times 10^{-18}\ \text{J}$. The kinetic energy component is $6.210 \times 10^{-19}\ \text{J}$, agreeing with the theoretical value $6.213 \times 10^{-19}\ \text{J}$ (relative error: 0.05%). The potential energy component is relatively small, reflecting the weak interatomic interactions in this dilute gas system.

### Heat Capacity $C_V$ (200‑400 K)

| Temperature | $C_V$ (J/K) |
|:---|:---|
| 200 K | $9.533 \times 10^{-20}$ |
| 250 K | $7.791 \times 10^{-21}$ |
| 300 K | $5.414 \times 10^{-21}$ |
| 350 K | $1.254 \times 10^{-20}$ |
| 400 K | $4.890 \times 10^{-21}$ |

The theoretical value for a monoatomic ideal gas is $2.07 \times 10^{-21}\ \text{J/K}$.


### Sensitivity Analysis

- **Time step sensitivity**: Results show significant dependence on $\Delta t$. $\Delta t = 1.0 \times 10^{-16}\ \text{s}$ balances accuracy and computational cost.
- **Total step sensitivity**: Increasing steps from 100 000 to 500 000 reduces statistical fluctuations and improves convergence.

## Dependencies

- FORTRAN compiler (gfortran recommended)
- Python 3.x with numpy and matplotlib (for plotting)

## How to Run

1. Clone this repository:

```bash
git clone https://github.com/zhangzhao-1/neon-md-simulation.git
cd neon-md-simulation
```

2. Compile the FORTRAN source code:
   
```bash
gfortran -o md_simulation main.f90
```

3. Run the simulation:

```bash
./md_simulation
```

The program will output energy, temperature, and heat capacity data to the terminal and/or output files.

## Possible Future Extensions

Based on the current simulation results and observed limitations, several promising directions are available for further improvement and expansion of this molecular‑dynamics code:

- **System‑size and finite‑size scaling**: Study larger‑particle systems and perform finite‑size scaling analysis to extrapolate observables toward the thermodynamic limit and obtain more accurate heat‑capacity values.
- **Enhance algorithmic efficiency**: Adopt more efficient numerical algorithms to support larger‑system‑size and longer‑time simulations while maintaining numerical precision.
- **Extended physical observables**: Compute radial distribution function (RDF) and diffusion coefficients to characterize structural and transport properties of the system.
- **Parameter‑space exploration**: Refine physical parameters through literature review and parameter tuning to bring computed thermodynamic quantities (e.g., heat capacity) closer to theoretical reference values.
- **Generalization to other species**: Extend the simulation framework to other noble‑gas systems such as Ar, Kr, and Xe.
- **Parallelization**: Add MPI / OpenMP parallelization to enable large‑scale production‑grade simulations.

## References
1. Lennard‑Jones, J. E. (1924). On the Determination of Molecular Fields. *Proc. R. Soc. Lond. A*, 106(738), 463‑477.
2. Frenkel, D., & Smit, B. (2002). *Understanding Molecular Simulation: From Algorithms to Applications*. Academic Press.
3. Allen, M. P., & Tildesley, D. J. (2017). *Computer Simulation of Liquids*. Oxford University Press.
