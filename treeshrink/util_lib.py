from math import log

def minVar_bisect(L):
# cut a list L into two sets L1, L2 such that
# the var(L1) + var(L2) is minimized
# return (max(L1)+min(L2))/2 as the cut-off
    L.sort()
    LL = sorted([log(x+1e-5) for x in L])

    sum_left = 0
    sumsq_left = 0
    sum_right = 0
    sumsq_right = 0

    for x in LL:
        sum_left += x
        sumsq_left += x*x

    n = float(len(LL))

    minVar = sumsq_left/n - (sum_left/n)**2
    cutoff = L[-1]
    
    k = len(L)-1

    while k > 0:
        x = LL[k]    
        sum_left -= x
        sumsq_left -= x*x
        sum_right += x
        sumsq_right += x*x
     
        var_left = sumsq_left/float(k) - (sum_left/float(k))**2
        var_right = sumsq_right/(n-k) - (sum_right/(n-k))**2
        var = var_left+var_right
        if var < minVar:
            minVar = var
            cutoff = L[int(k-1)]  #(L[int(k-1)] + L[int(k)])/2
        k -= 1    

    print(sum(x > cutoff for x in L))

    return cutoff

