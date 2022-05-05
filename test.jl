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

function generate_lattice(dim::Int, p)
    """
    Generates an Ising lattice based on a binomial distribution
    with probability p
    """

    d = Binomial(1, p)
    lattice = rand(d, dim, dim)
    for i = 1:dim
        for j = 1:dim
            if lattice[i, j] == 0
                lattice[i, j] = -1
            end
        end
    end
    return lattice
end


function vectorize_lattice(array::Array)::Array
    """
    vectorizes array row-wise
    """
    return vec(array')'
end

function reform_lattice(vector::Array)::Array
    """
    constructs sqaures lattice based on lattice vector
    """
    N = sqrt(length(vector))
    N = convert(Int64, N)
    lattice = zeros(N,N)
    for i = 1:N
        for j = 1:N
            k = get_global_domain(i, j, N)
            vector[k] = convert(Int64, vector[k])
            lattice[i,j] = vector[k]
        end
    end
    return convert(Array{Int64, 2}, lattice)
end

function get_local_domain(k, n)
    """
    maps a global domain k to a local domain (i,j)
    """
    if k%n == 0
        i = k ÷ n
        j = n
    else
        i = (k ÷ n) + 1
        j = k % n
    end
    return i,j
end

function get_global_domain(i::Int, j::Int, n::Int)::Int
    """
    function to map index pair (i,j) to global index k.
    Used to map a matrix onto a vector
    """
    return j + (i-1)*n
end

function get_top_nbr(m::Int, n::Int, i::Int, j::Int)::Int
    """
    returns the top neighbor index k given an index pair (i,j)
    """
    return get_global_domain(ifelse(i == 1, m, i-1), j, n)
end

function get_bottom_nbr(m::Int, n::Int, i::Int, j::Int)::Int
    """
    returns the bottom neighbor index k given an index pair (i,j)
    """
    return get_global_domain(ifelse(i == m, 1, i+1), j, n)
end

function get_left_nbr(m::Int, n::Int, i::Int, j::Int)::Int
    """
    returns the left neighbor index k given an index pair (i,j)
    """
    return get_global_domain(i, ifelse(j == 1, n, j-1), n)
end

function get_right_nbr(m::Int, n::Int, i::Int, j::Int)::Int
    """
    returns the right neighbor index k given an index pair (i,j)
    """
    return get_global_domain(i, ifelse(j == n, 1, j+1), n)
end

function get_nearest_nbrs(m::Int, n::Int, i::Int, j::Int)::Array
    """
    returns all 4 nearest neighbors of an index pair (i,j) mapped
    onto k global index
    """
    return [
    get_top_nbr(m::Int, n::Int, i::Int, j::Int),
    get_right_nbr(m::Int, n::Int, i::Int, j::Int),
    get_left_nbr(m::Int, n::Int, i::Int, j::Int),
    get_bottom_nbr(m::Int, n::Int, i::Int, j::Int)
    ]
end

function energy_ising(array::Array, vec::Array)::Float64
    """
    computes the energy configuration of Ising model following:
    H = -J*Σσσ
    """
    energy = 0
    m, n = size(array)


    for i = 1:m
        for j = 1:n
            k = get_global_domain(i, j, n)
            nbrs = get_nearest_nbrs(m, n, i, j) # nbrs is an index!!!!!
            energy_temp = 0
            for index in nbrs
                energy_temp += (vec[k] * vec[index])
            end
            energy += energy_temp
        end
    end
    return -energy/2
end

function get_magnetization(sigma_set, dim)
    """computes the normalized magnetization
    """
    M = 0
    for spin in sigma_set
        M += spin
    end
    return M / dim^2
end

function get_exact_energy(T)
    """
    computes the exact internal energy of system
    http://www.lps.ens.fr/~krzakala/ISINGMODEL.pdf
    """
    β = get_thermo_beta(T)
    m = 2 * sinh(2 * β) / cosh(2 * β)^2
    u =  - 1 / tanh(2 * β) * (1.0 + (2 * tanh(2 * β) ^ 2 - 1.0) * (2.0 / pi) * K(m ^ 2))
    return u
end

function get_exact_free_energy(T)
    """
    computes the exact free energy of the system
    http://www.lps.ens.fr/~krzakala/ISINGMODEL.pdf
    """
    β = get_thermo_beta(T)
    integral = hcubature(t -> log(cosh(2 * β) ^ 2 - sinh(2 * β) * (cos(t[1]) + cos(t[2]))),
        [0,0], [pi,pi])[1]
    f = (-1.0 / β) * (log(2) + integral / (2 * pi ^ 2))
    return f
end

function get_exact_magnetization(T)
    """
    computes exact magnetization of the system
    """
    β = get_thermo_beta(T)
    if T < 2.0 / log(1.0 + sqrt(2.0))
        m = (1.0 - sinh(2 * β) ^ (-4)) ^ (1.0 / 8)
    else
        m = 0
    end
    return m
end

function get_exact_entropy(T)
    """
    computes exact entropy of the system
    """
    s = (get_exact_energy(T) - get_exact_free_energy(T)) / T
    return s
end

function get_exact_cv(T)
    """
    approximates heat capacity with centered finite difference
    """
    cv = T * central_fdm(10, 1)(get_exact_entropy, T)
    return cv
end

function get_thermo_beta(T)
    """
    computes the thermodynamic beta given a temperature T
    """
    return 1 / T
end

function get_exact_properties(T)
    """
    wrapper function to compute exact thermodynamic quantities given
    a temperature range T
    """
    E = zeros(length(T))
    M = zeros(length(T))
    cv = zeros(length(T))
    s = zeros(length(T))
    println("Calculating exact results")
    for (i,temp) in ProgressBar(enumerate(T))
        s[i] = get_exact_entropy(temp)
        cv[i] = get_exact_cv(temp)
        E[i] = get_exact_energy(temp)
        M[i] = get_exact_magnetization(temp)
    end
    return E, cv, M
end

function markov_ising(sigma_set::Array, energy_init, β)
    """
    performs one time-step of the Metropolis algorithm
    Note that sigma_set is the vectorized Ising lattice
    and energy_init is the initial energy configuration of
    sigma_set
    """
    N = convert(Int, sqrt(length(sigma_set)))
    i,j = (rand(1:N), rand(1:N))
    k = get_global_domain(i, j, N)
    h = get_nearest_nbrs(N, N, i, j)
    h = sum(sigma_set[h])
    ΔE = 2 * h * sigma_set[k]
    Υ = exp(-β * ΔE)
    if rand() < Υ
        sigma_set[k] = -sigma_set[k]
        energy_init = energy_init + ΔE
    end

    return sigma_set, energy_init
end

function wolff_step_test(vec, n, T)
    i,j = (rand(1:n), rand(1:n))
    k = get_global_domain(i,j,n)
    state = vec[k]
    cluster_size = 0
    p = 1 - exp(-2*get_thermo_beta(T))
    cluster_queue = Queue{Int}()
    enqueue!(cluster_queue, k)
    vec[k] = -state
    while isempty(cluster_queue) == false
        idx = first(cluster_queue) #global domain
        i2, j2 = get_local_domain(idx, n)
        dequeue!(cluster_queue)
        #IDK what m_count is doing in original code
        cluster_size += 1
        nbrs = get_nearest_nbrs(n, n, i2, j2)
        for nnidx in nbrs
            if vec[nnidx] == state
                randn = rand()
                if randn < p
                    enqueue!(cluster_queue, nnidx)
                    vec[nnidx] = -state
                end
            end
        end
    end
    return vec, cluster_size
end


function wolff_cluster(vec, m, n, T)
    """
    performs one time steps of Wolff's clustering algorithm
    https://arxiv.org/pdf/cond-mat/0311623.pdf pg 7
    """
    β = get_thermo_beta(T)
    p = 1 - exp(-2*β)
    i,j = (rand(1:m), rand(1:n))
    k1 = get_global_domain(i,j,n)
    C = [k1]
    F_old = [k1]
    while isempty(F_old) == false
        F_new = []
        for item in F_old
            i2, j2 = get_local_domain(item, n)
            nbrs = get_nearest_nbrs(m,n,i2,j2)
            for neighbor in nbrs
             if vec[neighbor] == vec[item] && neighbor ∉ C
                if rand() < p
                    push!(F_new, neighbor)
                    push!(C, neighbor)
                end
             end
            end
        F_old = F_new
        end
        for spin in C
            vec[spin] *= -1
        end
    end
    energy = energy_ising(reform_lattice(vec), vec)
    return vec, energy
end

function sweep_metropolis(T::Float64, epoch::Float64, freq::Int64, lattice_dim::Int64, p::Float64)
    """
    runs simulation using Metropolis algorithm for one temperature value
    returns internal energy, heat capacity, and magnetization
    """

    lattice = generate_lattice(lattice_dim, p)
    lattice_set = vectorize_lattice(lattice)

    β = get_thermo_beta(T)

    energy = energy_ising(lattice, lattice_set)

    time = 1:epoch

    cv = 0
    E = []
    M = []

    for t in time
        lattice_set, energy = markov_ising(lattice_set, energy, β)
        if t > 0.45*epoch && t % freq == 0
            mag = get_magnetization(lattice_set, lattice_dim)
            push!(E, energy)
            push!(M, mag)
        end
    end
    cv = β^2*var(E)
    cv = cv / lattice_dim^2
    E = E ./ lattice_dim^2
    E = mean(E)
    M = mean(M)
    return E, cv, M
end

function metropolis_wrapper(T, epoch::Float64, freq::Int64, lattice_dim::Int64, p::Float64)
    """
    Generate thermodynamic data using Metropolis algorithm
    wrapper function of sweep_metropolis
    """
    E = zeros(length(T))
    M = zeros(length(T))
    Cv = zeros(length(T))
    println("Running Metropolis simulation...")
    for (index, temp) in ProgressBar(enumerate(T))
        energy, cv, mag = sweep_metropolis(temp, epoch, freq, lattice_dim, p)
        E[index] = energy
        Cv[index] = cv
        M[index] = mag
    end
    return E, Cv, M
end

function sweep_wolff(T::Float64, epoch::Float64, freq::Int64, lattice_dim::Int64, p::Float64)
    """
    runs simulation using Wolff algorithm for one temperature value
    returns internal energy, heat capacity, and magnetization
    """
    lattice = generate_lattice(lattice_dim, p)
    lattice = vectorize_lattice(lattice)
    energy = 0
    cv = 0
    spins_flipped = 0
    E = []
    M = []
    β = get_thermo_beta(T)

    counter = 0
    while spins_flipped < epoch
        lattice, cluster_size = wolff_step_test(lattice, lattice_dim, T)
        spins_flipped += cluster_size
        counter += 1
        if spins_flipped > 0.45*epoch && counter % freq == 0
            energy = energy_ising(reform_lattice(lattice), lattice)
            mag = get_magnetization(lattice, lattice_dim)
            push!(E, energy)
            push!(M, mag)
        end
    end

    cv = β^2*var(E)
    cv = cv / lattice_dim^2
    E = E ./ lattice_dim^2
    E = mean(E)
    M = mean(M)
    return E, cv, M
end

function wolff_wrapper(T, epoch::Float64, freq::Int64, lattice_dim::Int64, p::Float64)
    """
    Generate thermodynamic data using Wolff algorithm
    wrapper function of wolff_metropolis
    """
    E = zeros(length(T))
    M = zeros(length(T))
    Cv = zeros(length(T))
    println("Running Wolff simulation...")
    for (index, temp) in ProgressBar(enumerate(T))
        energy, cv, mag = sweep_wolff(temp, epoch, freq, lattice_dim, p)
        E[index] = energy
        Cv[index] = cv
        M[index] = mag
    end
    return E, Cv, M
end

function plot_results(metro_data, wolff_data, exact_data,
                    T_exact, T_sim, path::String, epoch::Float64, lattice_dim::Int64)
    """
    generates figures for data
    """
    println("Plotting and saving figures to: "*path)
    energy_plot = plot(T_exact, exact_data[1, :], title = "2D Ising Energy", label = "exact",
                    tick_direction = :out, legend = :bottomright, color = "black")
    scatter!(T_sim, metro_data[1, :], label = "Metropolis")
    scatter!(T_sim, wolff_data[1, :], label = "Wolff")
    xlabel!("T")
    ylabel!("E/N")
    savefig(energy_plot, path*"Ising_2D_energy_"*string(epoch)*"_"*string(lattice_dim)*".png")

    cv_plot = plot(T_exact, exact_data[2, :], title = "2D Ising Cv", label = "exact",
                    tick_direction = :out, legend = :best, color = "black")
    scatter!(T_sim, metro_data[2, :], label = "Metropolis")
    scatter!(T_sim, wolff_data[2, :], label = "Wolff")
    xlabel!("T")
    ylabel!("Cv/N")
    savefig(cv_plot, path*"Ising_2D_Cv_"*string(epoch)*"_"*string(lattice_dim)*".png")

    mag_plot = plot(T_exact, exact_data[3, :], title = "2D Ising Magnetization", label = "exact",
                    tick_direction = :out, legend = :best, color = "black")
    scatter!(T_sim, metro_data[3, :], label = "Metropolis")
    scatter!(T_sim, wolff_data[3, :], label = "Wolff")
    xlabel!("T")
    ylabel!("M/N")
    savefig(mag_plot, path*"Ising_2D_mag_"*string(epoch)*"_"*string(lattice_dim)*".png")
    println("---------- ### End of Program ### ----------")
end

### Set Simulation Params ###
path = "/Users/danielribeiro/Desktop/res/06_04_21/2D/"
T = 0.5:0.1:3.5
epoch = 1e6
freq = 100
lattice_dim = 10
p = 0.999

### Metropolis ###
energy_metro, cv_metro, mag_metro = metropolis_wrapper(T, epoch, freq, lattice_dim, p)
metro_data = [energy_metro cv_metro mag_metro]

### Wolff ###
# energy_w, cv_w, M_w = wolff_wrapper(T, epoch, freq, lattice_dim, p)
# wolff_data = [energy_w cv_w M_w]'

### Exact Solutions ###
# T_plot = 1:0.01:5
# exact_energy, exact_cv, exact_mag = get_exact_properties(T_plot)
# exact_data = [exact_energy exact_cv exact_mag]'


### Plot results ###
# plot_results(metro_data, wolff_data, exact_data, T_plot, T, path, epoch, lattice_dim)
