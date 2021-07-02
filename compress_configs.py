
import numpy as np
import matplotlib.pyplot as plt
from sweetsourcod.lempel_ziv import lempel_ziv_complexity
from os import listdir
from tqdm import tqdm
from scipy.special import iv as besseli
from scipy.integrate import quad

def read_file(fname, path):
    """
    reads a configuration txt file
    and opens it as a numpy array
    """
    config = np.loadtxt(path+fname, dtype=float)
    return config

def parse_directory(configs_path, n):
    """
    parse directory of config data and return vector of CID data
    """
    dir_list = listdir(configs_path)
    CID = np.zeros(len(dir_list), dtype = float)
    for (idx,fname) in tqdm(enumerate(dir_list)):
        if fname != ".DS_Store":
            vec = read_file(fname, configs_path)
            vec = discretize_data(vec, n = n)
            cid = get_cid(vec)
            CID[idx] = cid
    #del vec
    return CID

def discretize_data(vec, n = 255):
    """
    this method takes in a vec of data in radians and n bins
    and returns the discretized angle values
    """
    vec = vec + abs(min(vec))
    vec = vec % (2*np.pi)
    vec = vec / (2*np.pi/n)
    vec = np.int_(np.floor(vec))
    return vec

def get_cid(vec):
    """
    compute the CID of vec
    """
    cid = lempel_ziv_complexity(vec, "lz77")[1]
    return cid

def get_cid_rand(L, niter):
    """
    computes the CID of random binary sequence over niter
    iterations for statistical accuracy.
    """
    cid_rand = 0.0
    for i in range(1, niter+1, 1):
        rand_seq = np.random.randint(0, 2, size=L)
        cid_rand += lempel_ziv_complexity(rand_seq, "lz77")[1]
    cid_rand = cid_rand / niter
    return cid_rand

def get_entropy_(cid, cid_rand, n):
    """
    given a cid, cid_rand, and bin number n, computes the entropy
    of sequence
    """
    s = cid / cid_rand
    return (s/np.log2(n)) * np.log(2*np.pi)

def get_entropy(cid, cid_rand, T, n, L):
    """
    wrapper function of get_entropy_
    Takes in a vector of cids (cid)
    """
    S = np.zeros(len(cid), dtype = float)
    for (idx, cidval) in enumerate(cid):
        S[idx] = get_entropy_(cid, cid_rand, n)

    nsample = len(S) // len(T)
    ctr = 0
    idx = 0
    s = 0.0
    Savg = np.zeros(len(T), dtype = float)
    for sample in S:
        ctr += 1
        s += sample
        if ctr % nsample == 0:
            idx += 1
            s = s / nsample
            Savg[idx] = s
            s = 0.0
            ctr = 0

    return Savg

def get_exact_cv(T):
    """
    computes exact heat capacity of 1D
    XY model given a temperature T
    https://en.wikipedia.org/wiki/Classical_XY_model#One_dimension
    """
    K = 1.0 / T
    mu = besseli(1, K) / besseli(0, K)
    cv = K**2 * (1 - mu/K - mu**2)
    return cv

def get_integrand(T):
    """
    returns the integrand of S = ∫(Cv/T)dT
    given a temperature T
    """
    integrand = get_exact_cv(T) / T
    return integrand

def get_exact_entropy_(T):
    """
    computes the entropy based on the formula
    S = ∫(Cv/T)dT. Numerical integration is employed
    from 0.01 to T
    """
    s = quad(get_integrand, 0.01, T)[0]
    return s

def get_exact_entropy(T):
    """
    wrapper function of get_exact_entropy_
    pass range of temperatures, returns vector
    of exact normalized entropies
    """
    print("Calculating exact entropy...")
    S = zeros(len(T), dtype = float)
    for (idx, temp) in tqdm( enumerate(T) ):





path = "/panfs/roc/groups/7/mart5523/ribei040/XYModel_Julia/Simulation_Results/2021-07-01/configs/"

cid = parse_directory(path, 2)
