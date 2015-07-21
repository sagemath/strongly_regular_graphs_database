vmax = 1301

def paley():
    ans = {}
    for q in range(vmax):
        if q%4 == 1 and is_prime_power(q):
            ans[q,(q-1)/2,(q-5)/4,(q-1)/4]= (lambda q : graphs.PaleyGraph(q),q)
    return ans

def orthogonal_arrays():
    # notations from
    # http://www.win.tue.nl/~aeb/graphs/OA.html
    ans = {}
    for v,k,l,mu in feasible:
        if not is_square(v):
            continue
        n = int(sqrt(v))
        if k % (n-1):
            continue
        m = k/(n-1)
        if (l  != (m-1)*(m-2)+n-2 or
            mu != m*(m-1)):
            continue
        if designs.orthogonal_arrays.exists(m,n):
            ans[v,k,l,mu] = (lambda m,n : graphs.OrthogonalArrayBlockGraph(m, n), m,n)
    return ans

def johnson():
    # Using notations of http://www.win.tue.nl/~aeb/graphs/Johnson.html
    #
    # J(n,m) has the parameters v = m(m – 1)/2, k = 2(m – 2), λ = m – 2, μ = 4.
    ans = {}
    for v,k,l,mu in sorted(feasible):
        m = l + 2
        if (mu != 4 or
            k != 2*(m-2) or
            v != m*(m-1)/2):
            continue
        ans[v,k,l,mu] = (lambda m: graphs.JohnsonGraph(m,2), m)
    return ans

def steiner():
    # Using notations from http://www.win.tue.nl/~aeb/graphs/S.html
    #
    # The block graph of a Steiner 2-design S(2,m,n) has parameters:
    # v = n(n-1)/m(m-1), k = m(n-m)/(m-1), λ = (m-1)^2 + (n-1)/(m–1)–2, μ = m^2.
    ans = {}
    for v,k,l,mu in sorted(feasible):
        if mu <= 1 or not is_square(mu):
            continue
        m = int(sqrt(mu))
        n = (k*(m-1))/m+m
        if (v == (n*(n-1))/(m*(m-1)) and
            k == m*(n-m)/(m-1) and
            l == (m-1)**2 + (n-1)/(m-1)-2 and
            designs.balanced_incomplete_block_design(n,m,existence=True)):
            ans[v,k,l,mu] = (lambda n,m: graphs.IntersectionGraph(map(frozenset,designs.balanced_incomplete_block_design(n,m))),n,m)
    return ans

def affine_polar():
    # Using notations from http://www.win.tue.nl/~aeb/graphs/S.html
    #
    # VO+(2e,q) has parameters: v = q^(2e), k = (q^(e−1) + 1)(q^e − 1), λ =
    # q(q^(e−2) + 1)(q^(e−1) − 1) + q − 2, μ = q^(e−1)(q^(e−1) + 1)
    #
    # VO−(2e,q) has parameters v = q^(2e), k = (q^(e−1) - 1)(q^e + 1), λ =
    # q(q^(e−2) - 1)(q^(e−1) + 1) + q − 2, μ = q^(e−1)(q^(e−1) - 1)
    ans = {}
    for v,k,l,mu in sorted(feasible):
        if (not is_square(v) or
            not is_prime_power(v)):
            continue
        prime,power = is_prime_power(v,get_data=True)
        if power%2:
            continue
        for e in divisors(power/2):
            q = prime**(power/(2*e))
            assert v == q**(2*e)
            if (k == (q**(e-1) + 1)*(q**e-1) and
                l == q*(q**(e-2) + 1)*(q**(e-1)-1)+q-2 and
                mu== q**(e-1)*(q**(e-1) + 1)):
                ans[v,k,l,mu] = (lambda d,q : graphs.AffineOrthogonalPolarGraph(d,q,sign='+'),2*e,q)
            if (k == (q**(e-1) - 1)*(q**e+1) and
                l == q*(q**(e-2)- 1)*(q**(e-1)+1)+q-2 and
                mu== q**(e-1)*(q**(e-1) - 1)):
                ans[v,k,l,mu] = (lambda d,q : graphs.AffineOrthogonalPolarGraph(d,q,sign='-'),2*e,q)
    return ans

def orthogonal_polar():
    # Using notations from http://www.win.tue.nl/~aeb/graphs/srghub.html
    #
    # VO+(2e,q) has parameters: v = q^(2e), k = (q^(e−1) + 1)(q^e − 1), λ =
    # q(q^(e−2) + 1)(q^(e−1) − 1) + q − 2, μ = q^(e−1)(q^(e−1) + 1)
    #
    # VO−(2e,q) has parameters v = q^(2e), k = (q^(e−1) - 1)(q^e + 1), λ =
    # q(q^(e−2) - 1)(q^(e−1) + 1) + q − 2, μ = q^(e−1)(q^(e−1) - 1)
    ans = {}
    for v,k,l,mu in sorted(feasible):
        if (not is_square(v) or
            not is_prime_power(v)):
            continue
        prime,power = is_prime_power(v,get_data=True)
        if power%2:
            continue


constructions = {}
for d in [paley(), orthogonal_arrays(), johnson(), steiner(), affine_polar()]:
    constructions.update(d)

for g in [graphs.SchlaefliGraph, graphs.HoffmanSingletonGraph, graphs.SimsGewirtzGraph, graphs.M22Graph, graphs.CameronGraph, graphs.McLaughlinGraph]:
    constructions[g().is_strongly_regular(parameters=True)] = [g]

for (v,k,l,mu),f in constructions.items():
    kk  = v-k-1
    ll  = v-2*k+mu-2
    mmu = v-2*k+l
    if (v,kk,ll,mmu) not in constructions:
        constructions[v,kk,ll,mmu] = (lambda x: f[0](*x).complement(),f[1:])

leftovers = {x:brouwer[x]['comments'] for x in set(brouwer).difference(constructions) if brouwer[x]['status']=='exists'}

for (v,k,l,mu),y in sorted(leftovers.items()):
    print "({:<4}, {:<4}, {:<4}, {:<4}): {}".format(v,k,l,mu,y)
