using Plots, ProgressBars

include("simulator.jl")

struct Sensors
    positions::Vector{Vector{Float32}}
    times::Vector{Vector{Float32}}
    densities::Vector{Vector{Float32}}
    simulator::Simulator
end

function Sensors(initial_positions::Vector{Float32}, simulator::Simulator)
    positions = Vector{Float32}[]
    times = Vector{Float32}[]
    densities = Vector{Float32}[]
    for initial_position in initial_positions
        push!(positions, [initial_position])
        push!(times, [simulator.time.first])
        ρ = get_density(simulator, simulator.time.first, initial_position)
        push!(densities, [ρ])
        if ρ == undef
            @warn("Sensors are created before the initial condition of the simulator.")
        end
    end
    Sensors(positions, times, densities, simulator)
end

function get_speed(flux::Flux, ρ::Float32)
    (ρ > 0) ? flux.f(ρ) / ρ : flux.V_f
end

function compute(sensors::Sensors)
    Δ_T = sensors.simulator.Δ_T
    flux = sensors.simulator.flux
    time, _ = get_axes(sensors.simulator)
    sensors_data = zip(sensors.times, sensors.positions, sensors.densities)
    for (sensor_time, sensor_position, sensor_density) in ProgressBar(sensors_data)
        for t in time[1:(end-1)]
            x = sensor_position[end]
            if x > sensors.simulator.space.last
                break
            end
            ρ = get_density(simulator, t, x)
            v = get_speed(flux, ρ)
            push!(sensor_time, t + Δ_T)
            push!(sensor_position, x + Δ_T * v)
            push!(sensor_density, ρ)
        end
    end
end

function plot(sensors::Sensors)
    fig = plot(sensors.simulator)
    for (t, x) in zip(sensors.times, sensors.positions)
        Plots.plot!(t, x, legend=false, lc=:white)
    end
    fig
end

flux = Flux(ρ -> ρ * (1 - ρ), 0.5, 1.0)

simulator = Simulator(1.0f0, 2.0f0, 1000, flux, γ=0)
initial_condition(simulator, x -> 0.8 * x)
boundary_condition(simulator, identity, top_boundary_condition=true)
boundary_condition(simulator, x -> 0.9, top_boundary_condition=false)

probe_vehicles = Sensors([0.1f0, 0.5f0, 0.8f0], simulator)

compute(simulator)
compute(probe_vehicles)

plot(probe_vehicles) |> display