# MacroMicroSimulator

## Macro simulator

This package provides a numerical solution to the viscous / inviscid quasilinear equation

$\displaystyle \frac{\partial \rho}{\partial t}(t, x) + \frac{\partial f(\rho)}{\partial x}(t, x) = \gamma \frac{\partial^2 \rho}{\partial x^2} + G(t, x, \rho(t, x))$

for $(t, x) \in [0, T] \times [0, L]$ and the dissipation coefficient $\gamma$ is positive or null. The distribution $G$ is defined as $G(t, x, \rho(t, x)) = \sum_i g_i(t, \rho(t, x)) \delta_{x_i}(x)$ where $g_i$ is an inflow / outflow at $x_i$. The Dirac at $u$ is represented by $\delta_u$.

The unknown is the density function $\rho: [0, T] \times [0, L] \to [0, 1]$. Initial condition and boundary conditions at $x = 0$ and $x = L$ must be provided. The flux function $f$ is defined such that
* $f$ is concave with a maximum at $\rho_c$;
* $f(\rho) = V_f\rho + o(\rho)$.

These conditions are coming from [[1]](http://proceedings.mlr.press/v144/barreau21a.html).

```julia
using MacroMicroSimulator

# Definition of the equation
flux = MacroMicroSimulator.Flux(ρ -> ρ * (1 - ρ), ρ_c=0.5f0, V_f=1.0f0)
g1 = (0.5f0, (t, ρ) -> (t >= 1.0f0) ? 0.05f0 * (1 + cos((t - 1) * 6)) : 0.0f0)
equation = MacroMicroSimulator.Equation(L=1.0f0, T=2.0f0, flux=flux, gs=[g1])

# Definition of the macro-simulator
simulator = MacroMicroSimulator.Simulator(equation, N_L=500)
initial_condition(simulator, x -> 0.8 * x)
top_boundary_condition(simulator, identity)
bottom_boundary_condition(simulator, t -> 0.9)
compute(simulator) # We compute the solution

MacroMicroSimulator.plot(simulator) |> display
```

For more information, please see the documentation.

## Micro simulator

Based on the density function, it is possible to generate trajectories of particles. Their speed in the flow is defined as $V(\rho) = f(\rho) / \rho$ where $\rho$ is the density at the particle position.

```julia
# Definition of the micro-simulator
probe_vehicles_initial_condition = [
    (t=0.0f0, x=0.15f0),
    (t=0.0f0, x=0.8f0),
    (t=0.45f0, x=0.0f0),
    (t=1.5f0, x=0.5f0)
]
probe_vehicles = MacroMicroSimulator.Sensors(probe_vehicles_initial_condition, simulator)
compute(probe_vehicles) # We compute positions of the sensors

MacroMicroSimulator.plot(probe_vehicles) |> display
```

The output of this example is given below.

![simulation image](figs/simulation.png)

For more information, please see the documentation.

# Aknowledgement

Please cite [[1]](http://proceedings.mlr.press/v144/barreau21a.html) if you use this script.

[1] [M. Barreau, J. Liu, K. H. Johansson. Learning-based State Reconstruction
for a Scalar Hyperbolic PDE under noisy Lagrangian Sensing, *Proceedings of the 3rd Conference on Learning for Dynamics and Control*, PMLR 144:34-46, 2021](http://proceedings.mlr.press/v144/barreau21a.html)
