# Author: Daniel Ribeiro (ribei040@umn.edu)

using Distributions
using Colors
using Plots
using ProgressBars
using FiniteDifferences
using SpecialFunctions
using DelimitedFiles
using PyCall
using QuadGK

function read_config_file(fname::String, path::String)
    """
    function to read the text files. returns array with temperatures
    and array of configurations. make sure fname includes the path to
    the configurations
    """
    idx = split(fname, "_")
    config = readdlm(path * fname, ',')
    return config
end

function compress_directory(configs_path::String,
                    lz77_complexity_path::String,
                    sc,
                    n::Int)
    """
    parse directory of config data and return vector of entropy data
    and vector of temperature data
    """
    dir_list = readdir(configs_path)
    CID = zeros(length(dir_list))
    for (idx, fname) in enumerate(dir_list)
        if fname != ".DS_Store"
            vec = read_config_file(fname, configs_path)
            vec = discretize_data(vec, n = n)
            cid = get_cid(vec, sc)
            CID[idx] = cid
            save_cids(cid, lz77_complexity_path, idx)
            vec = nothing
            GC.gc()
        end
    end
    return CID
end

function get_cid(vec, sc)
    cid = sc.lempel_ziv_complexity(vec, "lz77")[2]
end

function get_entropy_(cid, sc,
                    n, L;
                    niter = 10)
    """
    computes entropy of vec based on LZ77 compression
    """
    cid_rand = get_cid_rand(L, niter, sc)
    s = cid / cid_rand
    return (s / log2(n)) * log(2pi)
end

function get_entropy(CID::Array, T, sc, n, L)
    S = zeros(length(CID))
    for (idx, cid) in enumerate(CID)
        S[idx] = get_entropy_(cid, sc, n, L)
    end
    nsample = length(S) / length(T)
    nsample = convert(Int64, nsample)
    ctr = 0
    idx = 0
    s = 0.0
    S_avg = zeros(length(T))
    for sample in S
        ctr += 1
        s += sample
        if ctr % nsample == 0
            idx += 1
            s = s / nsample
            S_avg[idx] = s
            s = 0.0
            ctr = 0
        end
    end
return S_avg
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

function get_exact_entropy_(T::Float64)::Float64
    s = quadgk(get_integrand, 0.01, T)[1]
end

function get_f(T::Float64)::Float64
    f = -T * log(2pi*besseli(0, 1.0/T))
end

function get_s(T::Float64)::Float64
    s = -central_fdm(3, 1)(get_f, T)
end

function get_exact_entropy(T)::Array
    """
    wrapper function of get_exact_entropy_
    pass range of temperatures, returns vector
    of exact normalize entropies
    """
    println("Calculating exact entropy...")
    S = zeros(length(T))
    for (idx, temp) in ProgressBar(enumerate(T))
        S[idx] = get_s(temp)
    end
    return S
end

function plot_entropy(s_sim, s_exact, T_sim, T_exact, n, plots_path)
    println("Plotting and saving figures to: " * plots_path)
    s_plot = plot(
        T_exact,
        s_exact,
        title = "1D XY Entropy",
        label = "exact",
        tick_direction = :out,
        legend = :bottomright,
        color = "black",
    )
    for (nval,s) in zip(n,s_sim)
        scatter!(T_sim, s, label = "n = "*string(nval))
    end
    xlabel!("T")
    ylabel!("S/N")
    savefig(s_plot, plots_path * "1d_xy_entropy.png")

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
    L = params["L"]
    plots_path = params["plots_path"]
    configs_path = params["configs_path"]
    lz77_complexity_path = params["lz77_complexity_path"]
    sweetsourcod = pyimport(params["sweetsourcod"])
    return (T_sim,
            T_exact,
            L,
            plots_path,
            configs_path,
            lz77_complexity_path,
            sweetsourcod)
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

function discretize_data(vec::Array; n=256)::Array
    """
    this method takes in a vec of data in radians and a desired
    binwidth, and returns the binned data
    """
    vec = vec .+ abs(minimum(vec))
    vec = vec .% (2pi)
    vec = vec ./ (2pi / n)
    vec = floor.(vec)
    vec = convert(Array{Int64, 2}, vec)
    return vec
end

function save_cids(cid::Float64, lz77_complexity_path::String,
                   idx::Int)
    fname = "xy_cid_"*"_"*string(idx)*".txt"
    open(lz77_complexity_path*fname, "w") do io
        writedlm(io, cid)
    end
end

function main()

    T_sim, T_exact, L, plots_path, configs_path, lz_complexity_path, sc = get_params()

    nbins = [2^i-1 for i = 1:8]
    S_sim = []
    println("compressing directory and computing entropy...")
    for n in ProgressBar(nbins)
        CID = compress_directory(configs_path,
                                    lz_complexity_path,
                                    sc,
                                    n)

        Sn = get_entropy(CID, T_sim, sc, n, L)
        push!(S_sim, Sn)
    end
    S_exact = get_exact_entropy(T_exact)

    plot_entropy(S_sim, S_exact, T_sim, T_exact, nbins, plots_path)

    println("---------- ### End of Program ### ----------")
end

cd("/Users/danielribeiro/XYModel_Julia")

main()


function get_f(T::Float64)::Float64
    f = -T * log(2pi*besseli(0, 1.0/T))
end

function get_s(T::Float64)::Float64
    s = -central_fdm(3, 1)(get_f, T)
end

println(get_s(100.0))
