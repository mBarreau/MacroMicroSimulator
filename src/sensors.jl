using Plots, ProgressBars

struct Sensors
    positions::Vector{Vector{Float32}}
    times::Vector{Vector{Float32}}
    densities::Vector{Vector{Float32}}
    simulator::Simulator
end

const InitialPosition = NamedTuple{(:t, :x),Tuple{Float32,Float32}}

function Sensors(initial_positions::Vector{InitialPosition}, simulator::Simulator)
    positions = Vector{Float32}[]
    times = Vector{Float32}[]
    densities = Vector{Float32}[]
    (t_min, x_min) = (simulator.equation.time.first, simulator.equation.space.first)
    for initial_position in initial_positions
        push!(positions, [max(initial_position.x, x_min)])
        push!(times, [max(initial_position.t, t_min)])
        ρ = get_density(simulator, initial_position.t, initial_position.x)
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
    flux = sensors.simulator.equation.flux
    sensors_data = zip(sensors.times, sensors.positions, sensors.densities)
    iterator = ProgressBar(sensors_data)
    set_description(iterator, "MicroSimulator")
    for (sensor_time, sensor_position, sensor_density) in iterator
        t = sensor_time[end]
        x = sensor_position[end]
        while t + Δ_T < sensors.simulator.equation.time.last
            ρ = get_density(sensors.simulator, t, x)
            v = get_speed(flux, ρ)
            x = x + Δ_T * v
            x >= sensors.simulator.equation.space.last && break
            push!(sensor_position, x)
            t = t + Δ_T
            push!(sensor_time, t)
            ρ = get_density(sensors.simulator, t, x)
            push!(sensor_density, ρ)
        end
    end
end

function plot(sensors::Sensors)
    fig = plot(sensors.simulator)
    for (t, x) in zip(sensors.times, sensors.positions)
        Plots.plot!(t, x; legend=false, lc=:white)
    end
    fig
end
