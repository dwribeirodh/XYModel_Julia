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

function get_subarray(patlen, iidx, jidx, lat1, latlen)
    """
    returns the sub array of lat1
    indexed by [iidx:iidx+patlen-1, jidx:jidx+patlen-1].
    latlen is the side length of lat1.
    This method implements periodic BCs
    """
    store = []
    for i = iidx:(iidx+patlen-1)
        for j = jidx:(jidx+patlen-1)
            iidxmap, jidxmap = mod(i, latlen), mod(j, latlen)
            iidxmap = ifelse(iidxmap==0, latlen, iidxmap)
            jidxmap = ifelse(jidxmap==0, latlen, jidxmap)
            push!(store, lat1[iidxmap, jidxmap])
        end
    end
    return store
end

function find_pattern(pat1, patlen, lat2, latlen)
    """
    returns true if pat1 exists in lat2.
    Returns false otherwise.
    patlen is the side length of pattern and
    latlen is the side length of lat2
    """
    is_pattern = false
    for i = 1:latlen
    # for i = 1:latlen-patlen
        for j = 1:latlen
        # for j = 1:latlen-patlen
            pat2 = get_subarray(patlen, i, j, lat2, latlen)
            if pat2 == pat1
                is_pattern = true
                break
            end
        end
    end
    return is_pattern
end

function find_longest_match(iidx, jidx, lat1, lat2, latlen)
    """
    Returns the longest match length found in lat2
    Pass in an index pair (iidx, jidx) that exist in lat1.
    the algorithm then performs linear search on lat2 to find
    the longest possible match.
    """
    is_pattern = true
    patlen = 0
    while is_pattern
        patlen += 1
        pat1 = get_subarray(patlen, iidx, jidx, lat1, latlen)
        is_pattern = find_pattern(pat1, patlen, lat2, latlen)
    end
    return (patlen-1)^2
end

function get_entropy(lat1, lat2, latlen)
    """
    computes entropy
    """
    sum_len = 0
    nsamples = 0
    while sum_len <= latlen^2
        randi, randj = rand(1:latlen), rand(1:latlen)
        sum_len += find_longest_match(randi, randj, lat1, lat2, latlen)
        nsamples +=1
    end
    len_avg = sum_len / nsamples
    entropy = log2(latlen^2) / len_avg
end

latlen = 10
d = Bernoulli(0.5)
lat1 = rand(d, latlen, latlen)
lat2 = rand(d, latlen, latlen)
l = zeros(latlen^5)
for s = ProgressBar(1:latlen^5)
    randi, randj = get_local_domain(s, latlen)
    l[s] = find_longest_match(randi, randj, lat1, lat2, latlen)
end
entropy = get_entropy(lat1, lat2, latlen)
