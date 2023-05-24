using MacroMicroSimulator

# Definition of the equation
flux = MacroMicroSimulator.Flux(ρ -> ρ * (1 - ρ), ρ_c=0.5f0, V_f=1.0f0)
g1 = (0.5f0, (t, ρ) -> (t >= 1.0f0) ? 0.05f0 * (1 + cos((t - 1) * 6)) : 0.0f0)
equation = MacroMicroSimulator.Equation(L=1.0f0, T=2.0f0, flux=flux, gs=[g1])

# Definition of the macro-simulator
simulator = MacroMicroSimulator.Simulator(equation, N_L=500)
initial_condition(simulator, x -> 0.8 * x)
top_boundary_condition(simulator, identity)
bottom_boundary_condition(simulator, x -> 0.9)
compute(simulator) # We compute the solution

MacroMicroSimulator.plot(simulator) |> display

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