using Intervals, Plots, ProgressBars

struct Flux
    f::Function
    ρ_c::Float32
    V_f::Float32
end

function Flux(f::Function; ρ_c::Float32, V_f::Float32)
    Flux(f, ρ_c, V_f)
end

const PointWiseFlow = Tuple{Float32,Function}

struct Equation
    space::Interval{Float32,Closed,Closed}
    time::Interval{Float32,Closed,Closed}
    flux::Flux
    γ::Union{Float32,Int}
    gs::Vector{PointWiseFlow}
end

function Equation(; L::Float32, T::Float32, flux::Flux, γ::Union{Float32,Int}=0, gs=PointWiseFlow[])
    Equation(Interval(0, L), Interval(0, T), flux, γ, gs)
end

struct Simulator
    equation::Equation
    ρ::Matrix{Float32}
    Δ_T::Float32
    Δ_L::Float32
end

function Simulator(equation::Equation, N_L::Int; α::Float32=0.8f0)
    Δ_L = (equation.space.last - equation.space.first) / N_L
    Δ_T = α * min(Δ_L / (2 * equation.flux.V_f), Δ_L^2 / (2 * equation.γ))
    N_T = ceil(Int32, (equation.time.last - equation.time.first) / Δ_T)
    ρ = Array{Float32}(undef, N_T, N_L)
    Simulator(equation, ρ, Δ_T, Δ_L)
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

function F(flux::Flux, u, v)
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

function G(gs::Vector{PointWiseFlow}, t, x::Interval{Float32,Open,Closed}, ρ)
    G = 0
    for (x_i, g_i) in gs
        if x_i in x
            G += g_i(t, ρ)
        end
    end
    G
end

function update(equation::Equation, ρ::Vector, Δ_T::Float32, Δ_L::Float32, t::Float32)
    flux = equation.flux
    ρ_out = copy(ρ)
    for i in 2:(size(ρ)[1]-1)
        pointwise_flow = G(equation.gs, t, Interval{Open,Closed}((i - 1) * Δ_L, i * Δ_L), ρ[i]) / Δ_L
        if equation.γ == 0
            gpdemi = F(equation.flux, ρ[i], ρ[i+1])
            gmdemi = F(equation.flux, ρ[i-1], ρ[i])
            ρ_out[i] = ρ[i] + Δ_T * (-(gpdemi - gmdemi) / Δ_L + pointwise_flow)
        else
            df = (flux.f(ρ[i+1]) - flux.f(ρ[i-1])) / (2 * Δ_L)
            ddρ = (ρ[i-1] - 2 * ρ[i] + ρ[i+1]) / (Δ_L^2)
            ρ_out[i] = ρ[i] + Δ_T * (-df + equation.γ * ddρ + pointwise_flow)
        end
    end
    clamp!(ρ_out, 0, 1)
    return ρ_out[2:(end-1)]
end

function compute(simulator::Simulator)
    if simulator.equation.γ == 0
        @info "MacroMicroSimulator: Finite volume method"
    else
        @info "MacroMicroSimulator: Finite difference method"
    end

    iterator = ProgressBar(2:size(simulator.ρ)[1])
    set_description(iterator, "MacroSimulator")
    for i in iterator
        simulator.ρ[i, 2:(end-1)] = update(simulator.equation, simulator.ρ[i-1, :], simulator.Δ_T, simulator.Δ_L, i * simulator.Δ_T)
    end
    simulator.ρ
end

function get_density(simulator::Simulator, t::Float32, x::Float32)
    time = simulator.equation.time
    space = simulator.equation.space
    i = floor(Int, (t - time.first) / simulator.Δ_T) + 1
    j = floor(Int, (x - space.first) / simulator.Δ_L) + 1
    simulator.ρ[i, j]
end

function get_axes(simulator::Simulator)
    time = simulator.equation.time
    space = simulator.equation.space
    t = LinRange(time.first, time.last, size(simulator.ρ)[1])
    x = LinRange(space.first, space.last, size(simulator.ρ)[2])
    return t, x
end

function plot(simulator::Simulator)
    t, x = get_axes(simulator)
    figure = heatmap(t, x, transpose(simulator.ρ), colormap=:turbo)
    ylabel!("Space")
    xlabel!("Time")
    figure
end
