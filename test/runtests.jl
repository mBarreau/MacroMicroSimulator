using MacroMicroSimulator

# Definition of the equation
flux = MacroMicroSimulator.Flux(ρ -> ρ * (1 - ρ), 0.5, 1.0)
g1 = (0.5f0, (t, ρ) -> (t >= 1.0f0) ? 10.0f0 : 0.0f0)
equation = MacroMicroSimulator.Equation(L=1.0f0, T=2.0f0, flux=flux, γ=0, gs=[g1])

# Definition of the macro-simulator
simulator = MacroMicroSimulator.Simulator(equation, 900)
initial_condition(simulator, x -> 0.8 * x)
top_boundary_condition(simulator, identity)
bottom_boundary_condition(simulator, x -> 0.9)
compute(simulator) # We compute the solution

MacroMicroSimulator.plot(simulator) |> display

# Definition of the micro-simulator
probe_vehicles = MacroMicroSimulator.Sensors([0.1f0, 0.5f0, 0.8f0], simulator)
compute(probe_vehicles) # We compute positions of the sensors

MacroMicroSimulator.plot(probe_vehicles) |> display