
import numpy as np
import matplotlib.pyplot as plt
from sweetsourcod.lempel_ziv import lempel_ziv_complexity
from os import listdir
from tqdm import tqdm

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
    for (idx,fname) in tqdm(enumerate(dir_list), ncols = 100):
        if fname != ".DS_Store":
            vec = read_file(fname, configs_path)
            vec = discretize_data(vec, n = n)
            cid = get_cid(vec)
            CID[idx] = cid
    del vec
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

# path = "/Users/danielribeiro/XY_Results/06_18_21/configs/"

# cid = parse_directory(path, 2)
