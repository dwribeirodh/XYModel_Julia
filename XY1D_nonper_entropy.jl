# Author: Daniel Ribeiro (ribei040@umn.edu)

#using Distributions
using Colors
using Plots
using ProgressBars
using SpecialFunctions: besseli
using DelimitedFiles: readdlm, writedlm
using PyCall: pyimport
using QuadGK: quadgk
using Dates: today

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
                    sc,
                    L,
                    n::Int)
    """
    parse directory of config data and return vector of CID data
    """
    dir_list = readdir(configs_path)
    CID = zeros(length(dir_list))
    for (idx, fname) in enumerate(dir_list)
        if fname != ".DS_Store"
            vec = read_config_file(fname, configs_path)
            vec = discretize_data(vec, n = n)
            cid = get_cid(vec, sc)
            CID[idx] = cid
            #save_cids(cid, lz77_complexity_path, idx)
        end
    end
    return CID
end

function get_cid(vec, sc)::Float64
    """
    compute the CID of vec
    """
    #check_mem()
    cid = sc.lempel_ziv_complexity(vec, "lz77")[2]
end

function get_entropy_(cid, cid_rand,
                    sc, n, L)
    """
    computes entropy of vec based on LZ77 compression
    """
    s = cid / cid_rand
    return (s / log2(n)) * log(2pi) #this is the one we've been using

end

function get_entropy(CID::Array, CID_rand, T, sc, n, L)
    """
    wrapper function of get_entropy_
    computes vector of entropies given a vector CID
    """
    S = zeros(length(CID))
    for (idx, cid) in enumerate(CID)
        S[idx] = get_entropy_(cid, CID_rand, sc, n, L)
    end
    nsample = length(S) ÷ length(T)
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
    """
    computes exact heat capacity of 1D
    XY model given a temperature T
    https://en.wikipedia.org/wiki/Classical_XY_model#One_dimension
    """
    K = 1.0 / T
    μ = besseli(1, K) / besseli(0, K)
    cv = K^2 * (1 - μ / K - μ^2)
    return cv
end

function get_integrand(T::Float64)::Float64
    """
    returns the integrand of S = ∫(Cv/T)dT
    given a temperature T
    """
    integrand = get_exact_cv(T) / T
    return integrand
end

function get_exact_entropy_(T::Float64)::Float64
    """
    computes the entropy based on the formula
    S = ∫(Cv/T)dT. Numerical integration is employed
    from 0.01 to T
    """
    s = quadgk(get_integrand, 0.01, T)[1]
end

function get_exact_entropy(T)::Array
    """
    wrapper function of get_exact_entropy_
    pass range of temperatures, returns vector
    of exact normalized entropies
    """
    println("Calculating exact entropy...")
    S = zeros(length(T))
    for (idx, temp) in ProgressBar(enumerate(T))
        S[idx] = get_exact_entropy_(temp) / log2(exp(1))
    end
    return S
end

function plot_entropy(s_sim, s_exact, T_sim, T_exact, n, plots_path)
    """
    plots entropy data and saves figures to
    plots_path
    Plots are:
        1) 1d_xy_entropy.png --> entropy vs temperature for different
        number of bins n
        2) 1d_xy_entropyvsbins.png --> entropy vs n for different
        temperatures
    """
    println("Plotting and saving figures to: " * plots_path)
    s_plot = plot(
        T_exact,
        s_exact,
        title = "1D XY Entropy",
        label = "exact",
        tick_direction = :out,
        legend = :bottomright,
        color = "black",
        linewidth = 1.5
    )

    for (nval,s) in zip(n,s_sim)
        plot!(T_sim, s, label = "n = "*string(nval))
    end
    xlabel!("T")
    ylabel!("S/N")
    savefig(s_plot, plots_path * "1d_xy_entropy.png")

    i, j = length(s_sim), length(s_sim[1])
    S = zeros(i, j)
    for idx = 1:i
        S[idx, :] = s_sim[idx]
    end

    n_plot = plot(
        tick_direction = :out,
        legend = :bottomright,
        )
    for (idx,temp) in enumerate(collect(T_sim))
        if temp in [0.4, 0.8, 2.0, 2.8, 3.2, 4.8]
            plot!(n, S[:,idx], label = "T = "*string(temp))
        end
    end
    xlabel!("n")
    ylabel!("S/N")
    savefig(n_plot, plots_path * "1d_xy_entropyvsbins.png")
end

function get_params()
    """
    reads config.txt file and returns script parameters
    """
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
    xy_path = params["XY_path"]
    return (T_sim,
            T_exact,
            L,
            xy_path)
end

function is_bool(name::SubString{String})::Bool
    """
    helper function of get_params
    """
    if name == "true" || name == "false"
        return true
    else
        return false
    end
end

function is_int(name::SubString{String})::Bool
    """
    helper function of get_params
    """
    if '.' in name || 'e' in name
        return false
    else
        return true
    end
end

function is_float(name::SubString{String})::Bool
    """
    helper function of get_params
    """
    if '.' in name && name != "sweetsourcod.lempel_ziv"
        return true
    else
        return false
    end
end

function convert_type(name::SubString{String})
    """
    helper function of get_params
    """
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
    this method takes in a vec of data in radians and n bins
    and returns the discretized angle values
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
    """
    saves cid values
    """
    fname = "xy_cid_"*"_"*string(idx)*".txt"
    open(lz77_complexity_path*fname, "w") do io
        writedlm(io, cid)
    end
end

function check_mem()
    """
    performs full garbage collection
    ensures no segfaults take place
    """
    GC.gc(true)
end

function main()
    """
    main method of script
    """
    sc = pyimport("sweetsourcod.lempel_ziv")
    today_date = string(today())
    T_sim, T_exact, L, xy_path = get_params()
    #cd("Simulation_Results")
    configs_path = xy_path*"/Simulation_Results/"*today_date*"/"*"configs"*"/"
    plots_path = xy_path*"/Simulation_Results/"*today_date*"/"*"plots"*"/"


    nbins = [2*i for i = 4:2:15]
    S_sim = []
    cid_rand = get_cid_rand(L, 5, sc)
    println("compressing directory and computing entropy...")
    for (idx,n) in ProgressBar(enumerate(nbins))
        CID = compress_directory(configs_path,
                                sc,
                                L,
                                n)

        Sn = get_entropy(CID, cid_rand, T_sim, sc, n, L)
        push!(S_sim, Sn)
    end
    S_exact = get_exact_entropy(T_exact)
    plot_entropy(S_sim, S_exact, T_sim, T_exact, nbins, plots_path)
    println("---------- ### End of Program ### ----------")
end

main()
