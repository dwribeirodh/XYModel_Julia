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

function generate_lattice(L::Int, is_random::Bool)::Array
    """
    generates a 1D xy lattice
    """
    if is_random
        lattice = rand(-pi, pi, L)
    else
        lattice = pi * ones(L)
    end
    return lattice
end

function find_nbrs_index(index::Int, L::Int)::tuple
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

function find_nbrs(vec::Array, index::Int)::tuple
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

function get_energy(vec::Array)::Float64
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
