using MacroMicroSimulator

flux = MacroMicroSimulator.Flux(ρ -> ρ * (1 - ρ), 0.5, 1.0)

simulator = MacroMicroSimulator.Simulator(1.0f0, 2.0f0, 1000, flux, γ=0)
initial_condition(simulator, x -> 0.8 * x)
top_boundary_condition(simulator, identity)
bottom_boundary_condition(simulator, x -> 0.9)

probe_vehicles = MacroMicroSimulator.Sensors([0.1f0, 0.5f0, 0.8f0], simulator)

compute(simulator)
compute(probe_vehicles)

MacroMicroSimulator.plot(probe_vehicles) |> display