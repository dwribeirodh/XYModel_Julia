# Author: Daniel Ribeiro (ribei040@umn.edu)

using Distributions
using Plots
using ProgressBars
using SpecialFunctions: besseli
using DelimitedFiles: readdlm, writedlm
using PyCall: pyimport
using QuadGK: quadgk
using Dates: today

function get_local_domain(k, n)
    """
    maps a global domain k to a local domain (i,j)
    """
    if k % n == 0
        i = k ÷ n
        j = n
    else
        i = (k ÷ n) + 1
        j = k % n
    end
    return i, j
end

function get_global_domain(i::Int, j::Int, L::Int)::Int
    """
    function to map index pair (i,j) to global index k.
    Used to map a matrix onto a vector
    """
    return j + (i - 1) * L
end

function get_top_nbr(L::Int, i::Int, j::Int)::Int
    """
    returns the top neighbor index k given an index pair (i,j)
    """
    return get_global_domain(ifelse(i == 1, L, i - 1), j, L)
end

function get_bottom_nbr(L::Int, i::Int, j::Int)::Int
    """
    returns the bottom neighbor index k given an index pair (i,j)
    """
    return get_global_domain(ifelse(i == L, 1, i + 1), j, L)
end

function get_left_nbr(L::Int, i::Int, j::Int)::Int
    """
    returns the left neighbor index k given an index pair (i,j)
    """
    return get_global_domain(i, ifelse(j == 1, L, j - 1), L)
end

function get_right_nbr(L::Int, i::Int, j::Int)::Int
    """
    returns the right neighbor index k given an index pair (i,j)
    """
    return get_global_domain(i, ifelse(j == L, 1, j + 1), L)
end

function get_nn_idx(L::Int, i::Int, j::Int)::Array
    """
    returns all 4 nearest neighbors of an index pair (i,j) mapped
    onto k global index
    """
    return [
        get_top_nbr(L::Int, i::Int, j::Int),
        get_right_nbr(L::Int, i::Int, j::Int),
        get_left_nbr(L::Int, i::Int, j::Int),
        get_bottom_nbr(L::Int, i::Int, j::Int),
    ]
end

function get_nn(lattice, L, i, j)
    """
    finds the spin values of the nearest neighbors
    of spin "angle"
    """
    nn_idx = get_nn_idx(L, i, j)
    nn = lattice[nn_idx]
    return nn
end

function generate_lattice(L::Real, d::Uniform{Float64}, is_random::Bool)::Array
    """
    generates 2D XY lattice
    """
    is_random ? lattice = rand(d, L^2) : lattice = pi * ones(L^2)
    return lattice
end

function reform_lattice(lattice, L::Int)
    """
    reconstruct square lattice from vector
    """
    sqrlat = zeros(Float64, (L, L))
    for idx = 1:L^2
        i, j = get_local_domain(idx, L)
        sqrlat[i, j] = lattice[idx]
    end
    return sqrlat
end

function get_energy(lattice, L)::Float64
    """
    computes the hamiltonian of the system given a
    configuration vec
    https://phas.ubc.ca/~berciu/TEACHING/PHYS502/PROJECTS/18BKT.pdf
    """
    energy = 0.0
    for (idx, angle) in enumerate(lattice)
        i, j = get_local_domain(idx, L)
        nn = get_nn(lattice, L, i, j)[1:2]
        energy += cos(angle - nn[1]) + cos(angle - nn[2])
    end
    return -energy
end

function metropolis_step(
    lattice,
    spin::Int,
    L::Int,
    energy::Float64,
    T::Float64,
    d::Uniform{Float64},
)
    """
    performs one time step of metropolis algorithm
    """
    β = 1.0 / T  # calculate thermo beta
    #rand_spin = rand(1:L^2) # select rand spin from vector
    randi, randj = get_local_domain(spin, L) # see sweep_metropolis for typewrite scheme implementation
    nnidx = get_nn_idx(L, randi, randj)
    nn = lattice[nnidx]
    dθ = pi * rand(d)
    ΔE = 0.0
    for nn_angle in nn
        ΔE += (
            cos(lattice[spin] + dθ - nn_angle) -
            cos(lattice[spin] - nn_angle)
        )
    end
    ΔE = -ΔE
    y = exp(-β * ΔE)
    flip = false
    if rand() < y
        lattice[spin] = lattice[spin] + dθ
        energy = energy + ΔE
        flip = true
    end
    return lattice, energy, flip
end

function sweep_metropolis(T, epoch, freq::Int64, L::Int64, d::Uniform{Float64})
    """
    given temperature T, runs a simulation using the Metropolis algorithm
    Params
        epoch: number of spin flips
        freq: frequency with which to output data
        lattice_dim: lattice length
    Returns
        E: numerical internal energy approximation
        Cv: numerical approximation for heat capacity
        M: numerical approximation for magnetization
    """
    β = 1.0 / T
    (
        T < 0.88 ? lattice = generate_lattice(L, d, false) :
        lattice = generate_lattice(L, d, true)
    )
    energy = get_energy(lattice, L)
    cv = 0.0
    E = []
    cv = []
    time = 0
    spin = 1
    while time < epoch
        lattice, energy, flip = metropolis_step(lattice, spin, L, energy, T, d)
        (spin == L^2) ? spin = 1 : spin += 1
        if flip
            time += 1
        end
        if (time > 0.5 * epoch) && (time % freq == 0)
            push!(E, energy)
        end
    end
    cv = (β^2 * var(E)) / L^2
    E = mean(E) / L^2
    return E, cv
end

function metropolis_simulation(
        T,
        epoch,
        freq::Int64,
        L::Int64,
        configs_path::String,
    )
    """
    generates thermodynamic data for 1D xy using Metropolis.
    """
    d = Uniform(-1.0, 1.0)
    println("Running Metropolis simulation...")
    E = zeros(length(T))
    Cv = zeros(length(T))
    for (idx, temp) in ProgressBar(enumerate(T))
        E[idx], Cv[idx] = sweep_metropolis(
            temp,
            epoch,
            freq,
            L,
            d
        )
    end
    return E, Cv
end

function sweep_metropolis_gif(
    T,
    epoch,
    freq::Int64,
    L::Int64,
    d::Uniform{Float64},
)
    β = 1.0 / T
    (
        T < 1.5 ? lattice = generate_lattice(L, d, false) :
        lattice = generate_lattice(L, d, true)
    )
    energy = get_energy(lattice, L)
    time = 0
    a = Animation()
    while time < epoch
        lattice, energy, flip = metropolis_step(lattice, L, energy, T, d)
        if flip
            time += 1
        end
        if time % freq == 0
            println("progresss = " * string((time / epoch) * 100) * "%")
            gif_vec = reform_lattice((lattice .% (2pi)), L)
            p = heatmap(gif_vec)
            frame(a, p)
        end
    end
    gif(a, "2dxy_metropolis.gif", fps = 100000)
end

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
    #is_random = params["is_random"]
    XY_path = params["XY_path"]
    return (T_sim, T_exact,
            epoch, freq,
            L, XY_path)
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

function main()
    T_sim, T_exact, epoch, freq, L, xy_repo_path = get_params()

    today_date = string(today())
    mkpath("Simulation_Results/"*today_date*"/configs/")
    mkpath("Simulation_Results/"*today_date*"/plots/")
    mkpath("Simulation_Results/"*today_date*"/energy/")
    configs_path = xy_repo_path*"/Simulation_Results/"*today_date*"/configs/"
    plots_path = xy_repo_path*"/Simulation_Results/"*today_date*"/plots/"

    e, cv = metropolis_simulation(T_sim, epoch, freq, L, configs_path)
    metro_res = [e cv]
    return metro_res, T_sim
    #TODO: implement exact results and plotting routines
end

data, T = main()
e = data[:, 1]
cv = data[:, 2]
scatter(T, e)
scatter(T, cv)
