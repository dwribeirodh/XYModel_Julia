

#using Main.SCTester: __init__

#sc = __init__()
using PyCall

L = 1000000
nbins = [2*i for i = 1:2:15]
ctr = 0
sc = pyimport("sweetsourcod.lempel_ziv")
# loop through all n
for n in nbins
    # loop through all configurations
    for i = 1:65
        global ctr += 1
        println(ctr)
        # analogous to reading a configuration from scratch
        config = rand(1:255, L)
        # calculate entropy
        # 1.1300458785794484e6 --> cid of random sequence of same L
        entropy = sc.lempel_ziv_complexity(config, "lz77")[2] / 1.1300458785794484e6
    end
end
