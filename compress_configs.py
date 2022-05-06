import numpy as np
# import matplotlib.pyplot as plt
from sweetsourcod.lempel_ziv import lempel_ziv_complexity
from os import listdir
from tqdm import tqdm

def read_file(fname, path):
    """
    reads a configuration txt file
    and opens it as a numpy array
    """
    config = np.loadtxt(path+fname, dtype='float64')
    return config

def parse_directory(configs_path, L, n, niter = 1000):
    """
    parse directory of config data and return vector of CID data
    """
    dir_list = listdir(configs_path)
    dir_list.sort()
    H = np.zeros(len(dir_list), dtype = float)
    cid_rand = get_cid_rand(L, niter)
    for (idx,fname) in tqdm(enumerate(dir_list)):
        if fname != ".DS_Store":
            lattice = read_file(fname, configs_path)
            lattice = discretize_data()
            cid = get_cid(lattice)
            H[idx] = get_entropy(cid, cid_rand)
    return H

def save_entropy(H, n, h_path):
    np.savetxt(h_path+"entropy_data"+"_n_"+str(n)+".txt", H, fmt='%f')

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

def get_entropy(cid, cid_rand):
    """
    given a cid, cid_rand, and bin number n, computes the entropy
    of sequence
    """
    s = cid / cid_rand
    return s

def discretize_data(vec, n):
    """
    this method takes in a vec of data in radians and n bins
    and returns the discretized angle values
    """
    vec = vec / (2*np.pi / n)
    vec = np.uint8(np.floor(vec))
    return vec

L = 512
configs_path = "/home/mart5523/ribei040/IsingModelJulia/Simulation_Results/2022-03-22/configs/"
h_path = "/home/mart5523/ribei040/IsingModelJulia/Simulation_Results/"
n = range(1, 255)
for n_bin in n:
    H = parse_directory(configs_path, n_bin, L**2)
    save_entropy(H, n_bin, h_path)
print("########## End of Program ##########")
