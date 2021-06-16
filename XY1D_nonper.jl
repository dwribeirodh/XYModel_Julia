# Author: Daniel Ribeiro (ribei040@umn.edu)

using Distributions
using Plots
using ProgressBars
using FiniteDifferences
using LinearAlgebra
using SpecialFunctions
using DelimitedFiles
using PyCall
using QuadGK

"Lattice Methods"

function generate_lattice(L::Real, is_random::Bool)::Array
    """
    generates a 1d xy lattice
    """
    if is_random
        lattice = rand(-pi:pi, L)
    else
        lattice = pi * ones(L)
    end
    lattice = convert(Array{Float64,1}, lattice)
    return lattice
end

function find_nbrs_index(i::Int, L::Real)
    """
    given an index i in lattice of length L,
    returns the indices of the nearest neighbors of i
    """
    n1 = 0
    n2 = 0
    if i == 1
        n1 = 2
        n2 = -1 # not a real neighbor
    elseif i == L
        n1 = L - 1
        n2 = -1 # not a real neighbor
    else
        n1 = i - 1
        n2 = i + 1
    end
    return n1, n2
end

function find_nbrs(vec::Array, index::Int)
    """
    wrapper function for find_neighbors_index
    given a lattice and index i, returns the spins of nearest
    neighbors
    """
    L = length(vec)
    if index == 1 || index == L
        s1, s2 = find_nbrs_index(index, L)
        return vec[s1], -1
    else
        s1, s2 = find_nbrs_index(index, L)
        return vec[s1], vec[s2]
    end
end

"Numerical Thermo Methods"

function get_energy(vec)::Float64
    """
    computes the hamiltonian of the system given a
    configuration vec
    https://phas.ubc.ca/~berciu/TEACHING/PHYS502/PROJECTS/18BKT.pdf
    """
    energy = 0
    for (i, angle) in enumerate(vec)
        if i != 1
            s1 = find_nbrs(vec, i)[1]
            energy += cos(angle - s1)
        end
    end
    return -energy
end

function get_magnetization(vec)::Float64
    mag = 0
    for spin in vec
        mag += spin
    end
    return mag
end

function get_thermo_beta(T::Float64)::Float64
    return 1.0 / T
end

"Metropolis Methods"
function metropolis_step(vec::Array, energy::Float64, T::Float64)
    """
    performs one time step of metropolis algorithm
    """
    L = length(vec)
    β = get_thermo_beta(T)
    rand_spin = rand(1:L)
    nbrs = find_nbrs(vec, rand_spin)
    dθ = rand(-pi:pi)
    ΔE = 0
    for nn_angle in nbrs
        ΔE += cos(vec[rand_spin] + dθ - nn_angle) - cos(vec[rand_spin] - nn_angle)
    end
    ΔE = -ΔE
    y = exp(-β * ΔE)
    if rand() < y
        vec[rand_spin] = vec[rand_spin] + dθ
        energy = energy + ΔE
    end
    return vec, energy
end

function sweep_metropolis(T, epoch,
                        freq::Int64, L::Int64,
                        is_random::Bool, configs_path::String,
                        energy_path::String)
    """
    given temperature T, runs a simulation using the Metropolis algorithm
    Params
        epoch: number of spin flips
        freq: frequency with which to output data
        lattice_dim: lattice length
        is_random: generate random lattice or completely correlated lattice
    Returns
        E: numerical internal energy approximation
        Cv: numerical approximation for heat capacity
        M: numerical approximation for magnetization
    """
    β = get_thermo_beta(T)
    lattice = generate_lattice(L, is_random)
    energy = get_energy(lattice)
    time = 1:epoch
    cv = 0
    E = []
    M = []
    counter = 0
    for t in time
        lattice, energy = metropolis_step(lattice, energy, T)
        if t > 0.50 * epoch
            counter += 1
            if t % freq == 0
                mag = get_magnetization(lattice)
                push!(E, energy)
                push!(M, mag)
            end
            if counter == (epoch/10)
                save_configs(lattice, configs_path, t, T)
                counter = 0
            end
        end
    end
    cv = β^2 * var(E) / L
    E = E ./ L
    save_configs(E, energy_path, 0.0, T, is_energy=true)
    E = mean(E)
    M = M ./ L
    M = mean(M)
    return E, cv, M
end

function metropolis_wrapper(T, epoch,
                            freq::Int64, L::Int64,
                            is_random::Bool, configs_path::String,
                            energy_path::String)
    """
    generates thermodynamic data for 1D xy using Metropolis.
    """
    println("Running Metropolis simulation...")
    E = zeros(length(T))
    M = zeros(length(T))
    Cv = zeros(length(T))
    for (index, temp) in ProgressBar(enumerate(T))
        e, cv, mag = sweep_metropolis(temp, epoch,
                                    freq, L,
                                    is_random, configs_path,
                                    energy_path)
        E[index] = e
        M[index] = mag
        Cv[index] = cv
    end
    return E, Cv, M
end

"Exact Thermo Methods"
function get_exact_internal_energy(β::Float64)::Float64
    u = -besseli(1, β) / besseli(0, β)
end

function get_exact_cv(T::Float64)::Float64
    K = 1.0 / T
    μ = besseli(1, K) / besseli(0, K)
    cv = K^2 * (1 - μ / K - μ^2)
    return cv
end

function get_integrand(T::Float64)::Float64
    integrand = get_exact_cv(T) / T
    return integrand
end

function get_exact_entropy(T::Float64)::Float64
    s = quadgk(get_integrand, 0.01, T)[1]
end

function get_exact_properties(T, path::String)
    println("Calculating Exact Properties")
    cv = zeros(length(T))
    u = zeros(length(T))
    s = zeros(length(T))
    for (index, temp) in ProgressBar(enumerate(T))
        cv[index] = get_exact_cv(temp)
        u[index] = get_exact_internal_energy(1.0 / temp)
        s[index] = get_exact_entropy(temp)
    end
    return u, cv
end

"Plot Methods"
function plot_data(exact_data, metro_data, T_exact, T_sim, path, epoch, L)

    println("Plotting and saving figures to: " * path)

    e_plot = plot(
        T_exact,
        exact_data[1, :],
        title = "XY Energy",
        label = "exact",
        tick_direction = :out,
        legend = :bottomright,
        color = "black",
    )
    scatter!(T_sim, metro_data[1, :], label = "Metropolis")
    xlabel!("T")
    ylabel!("E/N")
    savefig(e_plot, path * "xy_energy_" * string(epoch) * "_" * string(L) * ".png")

    cv_plot = plot(
        T_exact,
        exact_data[2, :],
        title = "XY Cv",
        label = "exact",
        tick_direction = :out,
        legend = :best,
        color = "black",
    )
    scatter!(T_sim, metro_data[2, :], label = "Metropolis")
    xlabel!("T")
    ylabel!("Cv/N")
    savefig(cv_plot, path * "xy_cv_" * string(epoch) * "_" * string(L) * ".png")
end

"Save Config Methods"
function save_configs(vec::Array, path::String,
                    spins_flipped::Float64, T::Float64;
                    is_energy = false)
    a = "xy_config_" * string(T) * "_" * string(spins_flipped) * ".txt"
    b = "xy_energy_" * string(T) * "_" * ".txt"
    fname = ""
    is_energy ? fname=b : fname=a
    open(path * fname, "w") do io
        writedlm(io, vec)
    end
end

"Config File Method"
function get_params()
    params = Dict()
    open("config.txt") do f
        while !eof(f)
            line = readline(f)
            if '=' in line
                line = split(line, "=")
                key = line[1]
                val = line[2]
                val = convert_type(val)
                params[key] = val
            end
        end
    end
    T_sim = params["T_sim_initial_value"]:params["T_sim_step_size"]:params["T_sim_final_value"]
    T_exact = params["T_exact_initial_value"]:params["T_exact_step_size"]:params["T_exact_final_value"]
    epoch = params["epoch"]
    freq = params["freq"]
    L = params["L"]
    is_random = params["is_random"]
    plots_path = params["plots_path"]
    energy_path = params["energy_path"]
    configs_path = params["configs_path"]
    return (T_sim, T_exact,
            epoch, freq,
            L, is_random,
            plots_path, energy_path,
            configs_path)
end

function is_bool(name::SubString{String})::Bool
    if name == "true" || name == "false"
        return true
    else
        return false
    end
end

function is_int(name::SubString{String})::Bool
    if '.' in name || 'e' in name
        return false
    else
        return true
    end
end

function is_float(name::SubString{String})::Bool
    if '.' in name && name != "sweetsourcod.lempel_ziv"
        return true
    else
        return false
    end
end

function convert_type(name::SubString{String})
    if is_float(name)
        name = parse(Float64, name)
    elseif is_int(name)
        name = parse(Int64, name)
    elseif is_bool(name)
        name = parse(Bool, name)
    else
        name = convert(String, name)
    end
    return name
end

"Main Method"
function main()
    T_sim, T_exact, epoch, freq, L, is_random, plots_path, energy_path, configs_path = get_params()

    e, cv, m = metropolis_wrapper(T_sim, epoch, freq, L, is_random, configs_path, energy_path)
    metro_res = [e cv m]'


    u, cv_exact= get_exact_properties(T_exact, plots_path)
    exact_res = [u cv_exact]'

    plot_data(exact_res, metro_res, T_exact, T_sim, plots_path, epoch, L)

    println("---------- ### End of Program ### ----------")
end

cd("/Users/danielribeiro/XYModel_Julia")

main()
