using Intervals, Plots, ProgressBars

struct Flux
    f::Function
    ρ_c::Float32
    V_f::Float32
end

struct Simulator
    space::Interval{Float32,Closed,Closed}
    time::Interval{Float32,Closed,Closed}
    γ::Float32
    flux::Flux
    ρ::Matrix{Float32}
    Δ_T::Float32
    Δ_L::Float32
end

function Simulator(L::Float32, T::Float32, N_L::Int, flux::Flux; α::Float32=0.8f0, γ::Union{Float32,Int}=0)
    Simulator(Interval(0, L), Interval(0, T), N_L, flux; α=α, γ=γ)
end

function Simulator(space::Interval{Float32,Closed,Closed}, time::Interval{Float32,Closed,Closed}, N_L::Int, flux::Flux; α::Float32=0.8f0, γ::Union{Float32,Int}=0)
    Δ_L = (space.last - space.first) / N_L
    Δ_T = α * min(Δ_L / flux.V_f, Δ_L^2 / (2 * γ))
    N_T = ceil(Int32, (time.last - time.first) / Δ_T)
    ρ = Array{Float32}(undef, N_T, N_L)
    Simulator(space, time, γ, flux, ρ, Δ_T, Δ_L)
end

function change_density(simulator::Simulator, ρ; line::String)
    if ρ isa Function
        t, x = get_axes(simulator)
        var = (line == "all") ? x : t
        density = ρ.(var)
    else
        density = ρ
    end

    if line == "top"
        simulator.ρ[:, end] = density
    elseif line == "bottom"
        simulator.ρ[:, 1] = density
    else
        simulator.ρ[1, :] = density
    end

    clamp!(simulator.ρ, 0, 1)
end

function initial_condition(simulator::Simulator, ρ_0)
    change_density(simulator, ρ_0, line="all")
end

function top_boundary_condition(simulator::Simulator, ρ)
    change_density(simulator, ρ, line="top")
end

function bottom_boundary_condition(simulator::Simulator, ρ)
    change_density(simulator, ρ, line="bottom")
end

function g(flux::Flux, u, v)
    f = flux.f
    ρ_c = flux.ρ_c
    if u > v
        if v >= ρ_c
            out = f(v)
        elseif u <= ρ_c
            out = f(u)
        else
            out = f(ρ_c)
        end
    else
        out = min(f(u), f(v))
    end
    return out
end

function update(flux::Flux, ρ::Vector, Δ_T::Float32, Δ_L::Float32, γ::Union{Float32,Int})
    ρ_out = copy(ρ)
    for i in 2:(size(ρ)[1]-1)
        if γ == 0
            gpdemi = g(flux, ρ[i], ρ[i+1])
            gmdemi = g(flux, ρ[i-1], ρ[i])
            ρ_out[i] = ρ[i] - Δ_T * (gpdemi - gmdemi) / Δ_L
        else
            df = (flux.f(ρ[i+1]) - flux.f(ρ[i-1])) / (2 * Δ_L)
            ddρ = (ρ[i-1] - 2 * ρ[i] + ρ[i+1]) / (Δ_L^2)
            ρ_out[i] = ρ[i] + Δ_T * (-df + γ * ddρ)
        end
    end
    clamp!(ρ_out, 0, 1)
    return ρ_out[2:(end-1)]
end

function compute(simulator::Simulator)
    for t in ProgressBar(2:size(simulator.ρ)[1])
        simulator.ρ[t, 2:(end-1)] = update(simulator.flux, simulator.ρ[t-1, :], simulator.Δ_T, simulator.Δ_L, simulator.γ)
    end
    simulator.ρ
end

function get_density(simulator::Simulator, t::Float32, x::Float32)
    i = floor(Int, (t - simulator.time.first) / simulator.Δ_T) + 1
    j = floor(Int, (x - simulator.space.first) / simulator.Δ_L) + 1
    simulator.ρ[i, j]
end

function get_axes(simulator::Simulator)
    t = LinRange(simulator.time.first, simulator.time.last, size(simulator.ρ)[1])
    x = LinRange(simulator.space.first, simulator.space.last, size(simulator.ρ)[2])
    return t, x
end

function plot(simulator::Simulator)
    t, x = get_axes(simulator)
    figure = heatmap(t, x, transpose(simulator.ρ), colormap=:turbo)
    ylabel!("Space")
    xlabel!("Time")
    figure
end