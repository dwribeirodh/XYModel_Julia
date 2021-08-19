using Distributions
using ProgressBars

function get_local_domain(k, n)
    """
    maps a global domain k to a local domain (i,j)
    """
    if k % n == 0
        i = k ÷ n
        j = n
    else
        i = (k ÷ n) + 1
        j = k % n
    end
    return i, j
end

function find_idx(value, array)
    idxs = []
    for (idx,i) in enumerate(array)
        if i == value
            push!(idxs, idx)
        end
    end
    return idxs
end

function get_subarray(L, i, j, array)
    """
    returns the sub array of array
    indexed by [i:i+L-1, j:j+L-1]
    """
    return array[i:i+L-1, j:j+L-1]
end

function find_pattern(pattern, lattice, latlen, patlen)
    """
    looks for a square array pattern within
        square array. L is the side length of pattern
        and arrlen is the
    """
    is_pattern = false
    for i in 1:latlen-patlen
        for j in 1:latlen-patlen
            p = lattice[i:i+patlen-1, j:j+patlen-1]
            if p == pattern
                is_pattern = true
                break
            end
        end
    end
    return is_pattern
end

function find_pattern_length(arr1, arr2, pos, latlen)
    """
    select a global index pos from array 1 such that
    pos ∈ [1, latlen^2]. this function will look
    for the longest matching subarray present in arr2.
    """
    iidx,jidx = get_local_domain(pos, latlen)
    is_pattern = true
    patlen = 1
    while is_pattern
        pattern = get_subarray(patlen, iidx, jidx, arr1)
        is_pattern = find_pattern(pattern, arr2, latlen, patlen)
        patlen += 1
    end
    return patlen-1
end

d = Normal(10, 0.1)
a = floor.(rand(d, 100, 100))
b = floor.(rand(d, 100, 100))
find_pattern_length(a, b, 10, 100)

# a = floor.(10 .* rand(10, 10))
# b = floor.(10 .* rand(10, 10))

a = ones(10, 10)
b = zeros(10,10)

for i = 1:4
    for j = 1:4
        a[i,j] = 0
    end
end
find_pattern_length(a, b, 1, 10)
