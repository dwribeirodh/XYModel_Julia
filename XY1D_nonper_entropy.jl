# Author: Daniel Ribeiro (ribei040@umn.edu)

using Distributions
using Colors
#using Images
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
using StatsBase


function read_config_file(fname::String, path::String)
    """
    function to read the text files. returns array with temperatures
    and array of configurations. make sure fname includes the path to
    the configurations
    """
    idx = split(fname, "_")
    temp = parse(Float64, idx[3])
    config = readdlm(path * fname, ',')
    return config, temp
end

function dir_parser(path::String, sc)
    """
    parse directory of config data and return vector of entropy data
    and vector of temperature data
    """
    println("Parsing directory and calculating entropy...")
    dir_list = readdir(path)
    S = zeros(length(dir_list))
    T = []
    for (idx,fname) in ProgressBar(enumerate(dir_list))
        if fname != ".DS_Store"
            vec, temp = read_config_file(fname, path)
            vec = reshape(vec, length(vec))
            vec = discretize_data(vec)
            s = get_entropy(vec, sc)
            S[idx] = s
            if !in(temp, T)
                push!(T, temp)
            end
        end
    end
    println(string(length(S)))
    println(string(length(T)))
    nsample = length(S) / length(T)
    nsample = convert(Int64, nsample)
    S_avg = zeros(nsample)
    ctr = 0
    idx = 0
    s = 0.0
    for sample in S
        ctr += 1
        s +=  get_entropy(sample, sc)
        if ctr % nsample == 0
            idx += 1
            s = s / nsample
            S_avg[idx] = s
            s = 0.0
            ctr = 0
        end
    end
    return S_avg, T
end

function get_entropy(vec, sc; niter = 100)
    """
    computes entropy of vec based on LZ77 compression
    """
    L = length(vec)
    cid_rand = get_cid_rand(L, niter, sc)
    cid = sc.lempel_ziv_complexity(vec, "lz77")[2]
    s = cid / cid_rand
    return s * log(2)
end

function get_cid_rand(L::Int, niter::Int, sc)::Float64
    """
    computes the CID of random binary sequence over niter
    iterations for statistical accuracy.
    """
    cid_rand = 0.0
    for i = 1:niter
        rand_seq = rand([0 1], L)
        cid_rand += sc.lempel_ziv_complexity(rand_seq, "lz77")[2]
    end
    cid_rand = cid_rand / niter
    return cid_rand
end

function get_exact_free_energy(T::Float64)::Float64
    """
    computes the exact free energy of 1D XY model w/ zero external
    magnetic field.
    """
    f = -T * log(2pi*besseli(0, 1.0/T))
    return f
end

function get_exact_entropy_(T::Float64)::Float64
    """
    calculates exact entropy of 1d xy model
    """
    s = - central_fdm(10, 1)(get_exact_free_energy, T)
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

function discretize_data(vec::Array; binwidth = 0.025)::Vector
    """
    this method takes in a vec of data in radians and a desired
    binwidth, and returns the binned data
    """
    n = 255
    vec = vec .+ abs(minimum(vec))
    vec = vec .% (2pi)
    vec = vec ./ (2pi / n)
    vec = floor.(vec) .+ 1
    return vec
end

function plot_entropy(s_sim, s_exact, T_sim, T_exact, fpath)
    println("Plotting and daving figures to: " * fpath)
    s_plot = plot(
        T_exact,
        s_exact,
        title = "1D Ising Entropy",
        label = "exact",
        tick_direction = :out,
        legend = :bottomright,
        color = "black",
    )
    scatter!(T_sim, s_sim, label = "LZ-77")
    xlabel!("T")
    ylabel!("S/N")
    savefig(s_plot, fpath * "1d_xy_entropy.png")
    println("---------- ### End of Program ### ----------")
end

function vec_to_array(vec)::Array
    array = zeros(length(vec))
    for (idx, item) in vec
        array[idx] = item
    end
    return array
end
cd("/Users/danielribeiro/XYModel_Julia")
### Set params ###
path1 = "/Users/danielribeiro/XY_Results/06_14_21/configs/"
fpath = "/Users/danielribeiro/XY_Results/06_14_21/thermo/"
sc = pyimport("sweetsourcod.lempel_ziv")
T_exact = 0.01:0.01:5

### parse configs and calculate entropy ###
S_sim, T_sim = dir_parser(path1, sc)

### calculate exact entropy ###
S_exact = get_exact_entropy(T_exact)

### plpot data ###
plot_entropy(S_sim, S_exact, T_sim, T_exact, fpath)

plot(T_sim, S_sim)
