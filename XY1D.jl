# Author: Daniel Ribeiro (ribei040@umn.edu)

using Distributions
using Colors
using Images
using Plots
using ProgressBars
using Elliptic
using HCubature
using FiniteDifferences
using Cuba
using DataStructures
using LinearAlgebra

function generate_lattice(L::Int, is_random::Bool)::Array
    """
    generates a 1D xy lattice
    """
    if is_random
        lattice = rand(-pi:pi, L)
    else
        lattice = pi * ones(L)
    end
    lattice = convert(Array{Float64, 1}, lattice)
    return lattice
end

function find_nbrs_index(i::Int, L::Int)
    """
    given an index i in lattice of length L,
    returns the indices of the nearest neighbors of i
    """
    n1 = 0
    n2 = 0
    if i == 1
        n1 = L
        n2 = 2
    elseif i == L
        n1 = L-1
        n2 = 1
    else
        n1 = i-1
        n2 = i+1
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
    s1, s2 = find_nbrs_index(index, L)
    s1 = vec[s1]
    s2 = vec[s2]
    return s1, s2
end

function get_energy(vec)::Float64
    """
    computes the hamiltonian of the system given a
    configuration vec
    https://phas.ubc.ca/~berciu/TEACHING/PHYS502/PROJECTS/18BKT.pdf
    """

    energy = 0
    for (i, angle) in enumerate(vec)
        s1, s2 = find_nbrs(vec, i)
        energy += cos(angle - s1) + cos(angle - s2)
    end
    return -energy/2
end

function get_magnetization(vec)::Float64
    mag = 0
    for spin in vec
        mag += spin
    end
    return mag
end

function get_thermo_beta(T::Float64)::Float64
    return 1/T
end

function metropolis_step(vec::Array, energy::Float64, T::Float64)
    """
    performs one time step of metropolis algorithm
    """
    # calculate length of lattice
    L = length(vec)
    #calculate β
    β = get_thermo_beta(T)
    # pick random spin in lattice
    rand_spin = rand(1:L)
    # find j neighbors of i
    nbrs = find_nbrs(vec, rand_spin)
    θj_1 = nbrs[1]
    θj_2 = nbrs[2]
    # generate random change in angle dθ
    dθ = rand(-pi:pi)
    θi_new = vec[rand_spin] + dθ
    θi_old = vec[rand_spin]
    #calculate ΔE
    ΔE = (cos(θi_new - θj_1) - cos(θi_old - θj_1) + cos(θi_new - θj_2) - cos(θi_old - θj_2))
    y = exp(-β*ΔE)
    if rand() < y
        vec[rand_spin] = θi_new
        energy = energy + ΔE
    end
    return vec, energy
end

function sweep_metropolis(T, epoch::Float64, freq::Int64, L::Int64, is_random::Bool)
    β = get_thermo_beta(T)
    lattice = generate_lattice(L, is_random)
    energy = get_energy(lattice)
    time = 1:epoch
    cv = 0
    E = []
    M = []
    for t in time
        lattice, energy = metropolis_step(lattice, energy, T)
        if t > 0.50*epoch && t % freq == 0
            mag = get_magnetization(lattice)
            push!(E, energy)
            push!(M, mag)
        end
    end
    cv = β^2 * var(E) / L
    E = mean(E) / L
    M = mean(M) / L
    return E, cv, M
end

function metropolis_wrapper(T, epoch::Float64, freq::Int64, L::Int64, is_random::Bool)
    println("Running Metropolis simulation...")
    E = zeros(length(T))
    M = zeros(length(T))
    Cv = zeros(length(T))
    for (index, temp) in ProgressBar(enumerate(T))
        energy, cv, mag = sweep_metropolis(temp, epoch, freq, L, is_random)
        E[index] = energy
        M[index] = mag
        Cv[index] = cv
    end
    return E, Cv, M
end

function test_deltaE(vec)
    energy_init = get_energy(vec)
    L = length(vec)
    rand_spin = rand(1:L)
    nbrs = find_nbrs(vec, rand_spin)
    dθ = rand(-pi:pi)
    ΔE = (cos(vec[rand_spin] + dθ - nbrs[1]) - cos(vec[rand_spin] - nbrs[1]) + cos(vec[rand_spin] + dθ - nbrs[2]) - cos(vec[rand_spin] - nbrs[2]))
    vec[rand_spin] = vec[rand_spin] + dθ
    test_energy = energy_init + ΔE
    energy_final = get_energy(vec)
    res = energy_final - test_energy
    return res
end

T = 0:0.2:5
epoch = 1e7
freq = 1000
L = 100
is_random = false

e, cv, m = metropolis_wrapper(T, epoch, freq, L, is_random)

plot(T, e)
plot(T, cv)
