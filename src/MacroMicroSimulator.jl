module MacroMicroSimulator

export Flux, Simulator, Sensors, initial_condition, top_boundary_condition, bottom_boundary_condition, compute, plot

include("simulator.jl")
include("sensors.jl")

end