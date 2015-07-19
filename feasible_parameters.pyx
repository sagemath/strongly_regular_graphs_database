from sage.rings.arith import is_square

def eigenvalues(int v,int k,int l,int mu):
    r"""
    1.3.1 of [Distance-regular graphs]
    """
    b = (mu-l)
    c = (mu-k)
    D = b**2-4*c
    if not is_square(D):
        return [None,None]
    return [(-b+sqrt(D))/2.0,
            (-b-sqrt(D))/2.0]

cdef bint seems_feasible(int v, int k, int l, int mu):
    cdef int lambda_r, lambda_s,l1,l2,K2,D,F,e
    if (v<0 or k<=0 or l<0 or mu<0 or
        k>=v-1 or l>=k or mu>=k or
        v-2*k+mu-2 < 0 or # lambda of complement graph >=0
        v-2*k+l    < 0 or # mu of complement graph >=0
        mu*(v-k-1) != k*(k-l-1)):
        return False

    if (v-1)*(mu-l)-2*k == 0: # conference
        return True

    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return False

    # 1.3.1 of [Distance-regular graphs]
    if ((s+1)*(k-s)*k) % (mu*(s-r)):
        return False

    return True

    #lambda_r = ((s+1)*(k-s)*k) / (mu*(s-r))
    #lambda_s = v-1-lambda_r
    #
    ## Krein bound (Biggs p3)
    #e = r
    #l1,l2 = r,s
    #K2 = (k+l2)*(l1+1)**2-(l2+1)*(k+l2+2*l1*l2)
    #if K2<0: # symmetry hypothesis
    #    return False
    #if ((e <= 2 and mu > e**2+e+3*k) or
    #    (e  > 2 and mu > e**2+e+2*k)):
    #    print "Hey22"
    #    return False
    #return True

cdef int v,k,l,mu

vmax = 1301
feasible = set()
for v in range(vmax):
    for k in range(1,v-1):
        for l in range(k-1 ):

            mu = k*(k-l-1)/(v-k-1)
            if seems_feasible(v,k,l,mu):
                feasible.add((v,k,l,mu))

from sage.structure.sage_object import load
assert set(feasible) == set(load("brouwer.sobj"))
