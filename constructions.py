def is_paley(v,k,l,mu):
    r"""
    Test if a Paley graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: t = is_paley(13,6,2,3); t
        sage: g = t[0](*t[1:]); g
        sage: g.is_strongly_regular(parameters=True)

        sage: t = is_paley(5,5,5,5); t
    """
    if (v%4 == 1 and is_prime_power(v) and
        k   == (v-1)/2 and
        l   == (v-5)/4 and
        mu  == (v-1)/4):
        return (lambda q : graphs.PaleyGraph(q),v)

def is_orthogonal_array_block_graph(v,k,l,mu):
    r"""
    Test if an Orthogonal Array graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: t = is_orthogonal_array_block_graph(64, 35, 18, 20); t
        sage: g = t[0](*t[1:]); g
        sage: g.is_strongly_regular(parameters=True)

        sage: t = is_orthogonal_array_block_graph(5,5,5,5); t
    """
    # notations from
    # http://www.win.tue.nl/~aeb/graphs/OA.html
    if not is_square(v):
        return
    n = int(sqrt(v))
    if k % (n-1):
        return
    m = k/(n-1)
    if (l  != (m-1)*(m-2)+n-2 or
        mu != m*(m-1)):
        return
    if designs.orthogonal_arrays.exists(m,n):
        return (lambda m,n : graphs.OrthogonalArrayBlockGraph(m, n), m,n)

def is_johnson(v,k,l,mu):
    r"""
    Test if a Johnson graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: t = is_johnson(10,6,3,4); t
        sage: g = t[0](*t[1:]); g
        sage: g.is_strongly_regular(parameters=True)

        sage: t = is_johnson(5,5,5,5); t
    """
    # Using notations of http://www.win.tue.nl/~aeb/graphs/Johnson.html
    #
    # J(n,m) has parameters v = m(m – 1)/2, k = 2(m – 2), λ = m – 2, μ = 4.
    m = l + 2
    if (mu == 4 and
        k  == 2*(m-2) and
        v  == m*(m-1)/2):
        return (lambda m: graphs.JohnsonGraph(m,2), m)

def is_steiner(v,k,l,mu):
    r"""
    Test if a Steiner graph is `(v,k,\lambda,\mu)`-strongly regular.

    A Steiner graph is the intersection graph of a Steiner set system. For more
    information, see http://www.win.tue.nl/~aeb/graphs/S.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: t = is_steiner(26,15,8,9); t
        sage: g = t[0](*t[1:]); g
        sage: g.is_strongly_regular(parameters=True)

        sage: t = is_steiner(5,5,5,5); t
    """
    # Using notations from http://www.win.tue.nl/~aeb/graphs/S.html
    #
    # The block graph of a Steiner 2-design S(2,m,n) has parameters:
    # v = n(n-1)/m(m-1), k = m(n-m)/(m-1), λ = (m-1)^2 + (n-1)/(m–1)–2, μ = m^2.
    if mu <= 1 or not is_square(mu):
        return
    m = int(sqrt(mu))
    n = (k*(m-1))/m+m
    if (v == (n*(n-1))/(m*(m-1)) and
        k == m*(n-m)/(m-1) and
        l == (m-1)**2 + (n-1)/(m-1)-2 and
        designs.balanced_incomplete_block_design(n,m,existence=True)):
        return (lambda n,m: graphs.IntersectionGraph(map(frozenset,designs.balanced_incomplete_block_design(n,m))),n,m)

def is_affine_polar(v,k,l,mu):
    r"""
    Test if an Affine Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/VO.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: t = is_affine_polar(81,32,13,12); t
        sage: g = t[0](*t[1:]); g
        sage: g.is_strongly_regular(parameters=True)

        sage: t = is_affine_polar(5,5,5,5); t
    """
    # Using notations from http://www.win.tue.nl/~aeb/graphs/VO.html
    #
    # VO+(2e,q) has parameters: v = q^(2e), k = (q^(e−1) + 1)(q^e − 1), λ =
    # q(q^(e−2) + 1)(q^(e−1) − 1) + q − 2, μ = q^(e−1)(q^(e−1) + 1)
    #
    # VO−(2e,q) has parameters v = q^(2e), k = (q^(e−1) - 1)(q^e + 1), λ =
    # q(q^(e−2) - 1)(q^(e−1) + 1) + q − 2, μ = q^(e−1)(q^(e−1) - 1)
    if (not is_square(v) or
        not is_prime_power(v)):
        return
    prime,power = is_prime_power(v,get_data=True)
    if power%2:
        return
    for e in divisors(power/2):
        q = prime**(power/(2*e))
        assert v == q**(2*e)
        if (k == (q**(e-1) + 1)*(q**e-1) and
            l == q*(q**(e-2) + 1)*(q**(e-1)-1)+q-2 and
            mu== q**(e-1)*(q**(e-1) + 1)):
            return (lambda d,q : graphs.AffineOrthogonalPolarGraph(d,q,sign='+'),2*e,q)
        if (k == (q**(e-1) - 1)*(q**e+1) and
            l == q*(q**(e-2)- 1)*(q**(e-1)+1)+q-2 and
            mu== q**(e-1)*(q**(e-1) - 1)):
            return (lambda d,q : graphs.AffineOrthogonalPolarGraph(d,q,sign='-'),2*e,q)

constructions = {
     ( 27,  16, 10,  8): [graphs.SchlaefliGraph],
     ( 50,   7,  0,  1): [graphs.HoffmanSingletonGraph],
     ( 56,  10,  0,  2): [graphs.SimsGewirtzGraph],
     ( 77,  16,  0,  4): [graphs.M22Graph],
     (231,  30,  9,  3): [graphs.CameronGraph],
     (275, 112, 30, 56): [graphs.McLaughlinGraph]}

test_functions = [is_paley, is_johnson,
                  is_orthogonal_array_block_graph,
                  is_steiner, is_affine_polar]

for params in feasible:
    for f in test_functions:
        ans = f(*params)
        if ans is not None:
            constructions[params] = ans
            break

# Complement SRGs
for (v,k,l,mu),f in constructions.items():
    kk  = v-k-1
    ll  = v-2*k+mu-2
    mmu = v-2*k+l
    if (v,kk,ll,mmu) not in constructions:
        constructions[v,kk,ll,mmu] = (lambda x: f[0](*x).complement(),f[1:])

leftovers = {x:brouwer[x]['comments'] for x in set(brouwer).difference(constructions) if brouwer[x]['status']=='exists'}

for (v,k,l,mu),y in sorted(leftovers.items()):
    print "({:<4}, {:<4}, {:<4}, {:<4}): {}".format(v,k,l,mu,y)
