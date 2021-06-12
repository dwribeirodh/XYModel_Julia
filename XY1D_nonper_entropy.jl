# Author: Daniel Ribeiro (ribei040@umn.edu)

using Distributions
using Colors
using Images
using Plots
using ProgressBars
using Elliptic
using HCubature
using FiniteDifferences
using DataStructures
using LinearAlgebra
using SpecialFunctions
using DelimitedFiles
using PyCall


function read_config_file(fname::String, path::String)
    """
    function to read the text files. returns array with temperatures
    and array of configurations. make sure fname includes the path to
    the configurations
    """
    idx = split(fname, "_")
    temp = parse(Float64, idx[3])
    config = readdlm(path*fname, ',')
    config = map_vec(config)
    return config, temp
end

function dir_parser(path::String, sc)
    """
    parse directory of config data and return vector of entropy data
    and vector of temperature data
    """
    println("Parsing directory...")
    dir_list = readdir(path)
    configs = []
    T = []
    for fname in ProgressBar(dir_list)
        if fname != ".DS_Store"
            vec, temp = read_config_file(fname, path)
            vec = vec .+ 2pi
            vec = map_vec(vec)
            push!(configs, vec)
            if !in(temp, T)
                push!(T, temp)
            end
        end
    end
    nsamples = length(configs) / length(T)
    S = zeros(length(T))
    ctr = 0
    s = 0
    idx = 0
    println("Calculating entropy...")
    for lattice in ProgressBar(configs)
        s += get_entropy(lattice, sc)
        ctr += 1
        if ctr == nsamples
            idx += 1
            s = s/nsamples
            S[idx] = s
            counter = 0
            s = 0
        end
    end
    return S, T
end

function map2twopi(angle::Float64)::Float64
    """
    maps any angle back to [0, 2π]
    """
    angle = angle % 2pi
    return angle
end

function map_vec(vec::Array)::Array
    """
    wrapper function of map2twopi. This takes a vector
    and returns the vector with angles mapped back to [0, 2π]
    """
    for (idx,spin) in enumerate(vec)
        vec[idx] = map2twopi(spin)
    end
    return vec
end

function get_entropy(vec::Array, sc; niter = 100)
    """
    computes entropy of vec based on LZ77 compression
    """
    L = length(vec)
    cid_rand = get_cid_rand(L, niter, sc)
    cid = sc.lempel_ziv_complexity(vec, "lz77")[2]
    s = cid / cid_rand
    return s*log(2)
end

function get_cid_rand(L::Int, niter::Int, sc)::Float64
    """
    computes the CID of random binary sequence over niter
    iterations for statistical accuracy.
    """
    cid_rand = 0.0
    for i = 1:niter
        rand_seq = rand(-pi:pi, L)
        cid_rand += sc.lempel_ziv_complexity(rand_seq, "lz77")[2]
    end
    cid_rand = cid_rand / niter
    return cid_rand
end

function get_exact_entropy_(T::Float64)::Float64
    """
    calculates exact entropy of 1d xy model
    """
    s = -log(2pi*besseli(0, 1.0/T)) + 1.0/T*(besseli(1, 1.0/T) / besseli(0, 1.0/T))
    return s
end

function get_exact_entropy(T)::Array
    """
    wrapper function of get_exact_entropy_
    pass vector of temperatures, returns vector
    of exact normalize entropies
    """
    println("Calculating exact entropy...")
    S = zeros(length(T))
    for (idx, temp) in ProgressBar(enumerate(T))
        S[idx] = get_exact_entropy_(temp)
    end
    return S
end

function discretize_data()
    #TODO: implement
end

function plot_entropy(s_sim, s_exact, T_sim, T_exact, fpath)
    println("Plotting and daving figures to: " * fpath)
    s_plot = plot(T_exact, s_exact, title = "1D Ising Entropy", label = "exact",
                tick_direction = :out, legend = :bottomright, color = "black")
    scatter!(T_sim, s_sim, label = "LZ-77")
    xlabel!("T")
    ylabel!("S/N")
    savefig(s_plot, fpath*"1d_xy_entropy.png")
    println("---------- ### End of Program ### ----------")
end

path1 = "/Users/danielribeiro/XY_Results/06_11_21/configs/"
fpath = "/Users/danielribeiro/XY_Results/06_11_21/thermo_data/"

sc = pyimport("sweetsourcod.lempel_ziv")
T_exact = 0.01:0.01:5

S_sim, T_sim = dir_parser(path1, sc)
S_exact = get_exact_entropy(T_exact)

plot_entropy(S_sim, S_exact, T_sim, T_exact, fpath)
