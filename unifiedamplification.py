import scipy.stats as stats
from scipy.special import comb
import numpy as np
from datetime import datetime
import json
from time import process_time


def Delta(ep, alpha, r0, r1, n, p, tol=1e-12):
    rate = r0/(r0+r1)

    def fvlow(c):
        # vectorized computation for a list of c (required by scipy.stats.rv_discrete.expect)
        g01 = (1-alpha-alpha*p)*(n-(c+1))/(1-r0-r1)
        low01 = ((np.exp(ep)*p-1)*alpha*(c+1)/r1+(np.exp(ep)-1)*g01)/(alpha*(p/r0-1/r1+np.exp(ep)*(p/r1-1/r0)))
        g2 = (1-alpha-alpha*p)*(n-c)/(1-r0-r1)
        low2 = ((np.exp(ep)*p-1)*alpha*c/r1+(np.exp(ep)-1)*g2)/(alpha*(p/r0-1/r1+np.exp(ep)*(p/r1-1/r0)))

        v0 = np.array([1-stats.binom.cdf(int(np.ceil(low01[i]))-2, c[i], rate) for i in range(len(c))])
        v1 = np.array([1-stats.binom.cdf(int(np.ceil(low01[i]))-1, c[i], rate) for i in range(len(c))])
        v2 = np.array([1-stats.binom.cdf(int(np.ceil(low2[i]))-1, c[i], rate) for i in range(len(c))])

        # linear transformation (and weights are independent from c)
        p0 = p*alpha*v0+alpha*v1+(1-alpha-p*alpha)*v2
        p1 = alpha*v0+p*alpha*v1+(1-alpha-p*alpha)*v2
        return p0-np.exp(ep)*p1

    def fvhigh(c):
        g01 = (1-alpha-alpha*p)*(n-(c+1))/(1-r0-r1)
        high01 = ((np.exp(-ep)*p-1)*alpha*(c+1)/r1+(np.exp(-ep)-1)*g01)/(alpha*(p/r0-1/r1+np.exp(-ep)*(p/r1-1/r0)))
        g2 = (1-alpha-alpha*p)*(n-c)/(1-r0-r1)
        high2 = ((np.exp(-ep)*p-1)*alpha*c/r1+(np.exp(-ep)-1)*g2)/(alpha*(p/r0-1/r1+np.exp(-ep)*(p/r1-1/r0)))

        v0 = np.array([stats.binom.cdf(int(np.floor(high01[i]-1)), c[i], rate) for i in range(len(c))])
        v1 = np.array([stats.binom.cdf(int(np.floor(high01[i])), c[i], rate) for i in range(len(c))])
        v2 = np.array([stats.binom.cdf(int(np.floor(high2[i])), c[i], rate) for i in range(len(c))])

        p0 = p*alpha*v0+alpha*v1+(1-alpha-p*alpha)*v2
        p1 = alpha*v0+p*alpha*v1+(1-alpha-p*alpha)*v2
        return p1-np.exp(ep)*p0

    delta0 = stats.binom.expect(fvlow, args=(n-1, r0+r1), lb=0, ub=n-1, tolerance=tol, maxcount=n, chunksize=32)
    if r0 == r1:
        # by the symmetry of r0 and r1, delta0 always equals to delta1
        delta1 = delta0
    else:
        delta1 = stats.binom.expect(fvhigh, args=(n-1, r0+r1), lb=0, ub=n-1, tolerance=tol, maxcount=n, chunksize=32)
    #print("delta01", delta0, delta1)

    return max(delta0, delta1)


def amplificationUB(p, beta, q0, q1, n, delta, T):
    # privacy amplification upper bound
    alpha = beta/(p-1)
    r0 = alpha*p/q0
    r1 = alpha*p/q1
    return amplificationUBCore(p, alpha, r0, r1, n, delta, T)


def amplificationUBCore(p, alpha, r0, r1, n, delta, T):
    epL = 0
    epH = np.log(p) # one may also use the closed-form bound as a starting point

    # tolerance in scipy.stats.rv_discrete.expect
    tol = max(1e-19, min(1e-14, delta/n/10.0))
    for t in range(T):
        ep = (epL+epH)/2.0
        # add tol*n to delta, extremely conservative for rigidness, one may modify the scipy.stats.rv_discrete.expect to get tighter bounds
        if Delta(ep, alpha, r0, r1, n, p, tol)+tol*n   > delta:
            epL = ep
        else:
            epH = ep
    return epH


def amplificationLB(p, beta, q0, q1, n, delta, T):
    # privacy amplification lower bound
    alpha = beta/(p-1)
    r0 = alpha*p/q0
    r1 = alpha*p/q1
    return amplificationLBCore(p, alpha, r0, r1, n, delta, T)


def amplificationLBCore(p, alpha, r0, r1, n, delta, T):
    epL = 0
    epH = np.log(p)

    # tolerance in scipy.stats.rv_discrete.expect
    tol = max(1e-19, min(1e-14, delta/n/10.0))
    for t in range(T):
        ep = (epL+epH)/2.0
        # there is no tol*n, as scipy.stats.rv_discrete.expect always underestimate non-negative delta
        if Delta(ep, alpha, r0, r1, n, p, tol) > delta:
            epL = ep
        else:
            epH = ep
    return epL


def computeUBParameters(epsilon, mechanism, options):
    # variation-ratio parameters for upper bounds
    p = np.exp(epsilon)
    beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+1)
    q = np.exp(epsilon)
    if mechanism in ['laplace']:
        beta = 1-np.exp(-epsilon/2)
    elif mechanism in ['piecewise']:
        beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+np.exp(epsilon/2))
    elif mechanism in ['krr']:
        k = int(options)
        beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+k-1)
    elif mechanism in ['subset']:
        d, k = options[0], options[1]
        beta = (np.exp(epsilon)-1)*(comb(d-1, k-1)-comb(d-2, k-2))/(np.exp(epsilon)*comb(d-1, k-1)+comb(d-1, k))
    elif mechanism in ['localhash']:
        l = int(options)
        beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+l-1)
    elif mechanism in ['hardamard']:
        K, s = options[0], options[1]
        beta = (s*(np.exp(epsilon)-1)/2)/(s*np.exp(epsilon)+K-s)
    elif mechanism in ['hardamardB']:
        K, s = options[0], options[1]
        beta = (s*(np.exp(epsilon)-1))/(s*np.exp(epsilon)+K-s)
    elif mechanism in ['collision']:
        if hasattr(options, '__iter__'):
            s, l = options[0], options[1]
        else:
            s = options
            l = int(s*np.exp(epsilon)+2*s-1)
        beta = s*(np.exp(epsilon)-1)/(s*np.exp(epsilon)+l-s)
    elif mechanism in ['cheu']:
        m, epc, deltac = options[0], options[1], options[2]
        right = 33.0/5*(((np.exp(epc)+1)/(np.exp(epc)-1))**2)*np.log(4/deltac)
        assert m*m-4*m*right >= 0
        f = (m-np.sqrt(m*m-4*m*right))/(2*m) # the p in the original paper
        #print("mechanism", mechanism, m, 4*right, f)
        f = min(1/2-0.00000001, f)
        p = (1-f)**2/(f*f)
        beta = 1-2*f
        q = (1-f)/f
    elif mechanism in ['balls2bins']:
        d, s = options[0], options[1]
        p = 10**5  # simulate +inf
        beta = (p-1)/(p+1)
        q = d/s
    elif mechanism in ['krr_sep_best']:
        k = int(options)
        beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+k-1)
    elif mechanism in ['krr_sep_worst']:
        k = int(options)
        beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+2-1)
    elif mechanism in ['krr_para_basic']:
        k = int(options)
        beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+1)
    elif mechanism in ['krr_para_advanced']:
        k = int(options)
        beta = 0
        for h in range(int(np.log2(k))):
            beta = (np.exp(epsilon)-1)/(np.exp(epsilon)+2**(h+1)-1)/int(np.log2(k))
    return p, beta, q, q



if __name__ == '__main__':
    epsilons = np.array([0.1, 1.0, 3.0, 5.0]) # local epsilon for single-message protocols
    #epsilons = np.array([1.0, 2.0, 3.0, 4.0])
    #epsilons = np.arange(0.01, 1.02, 0.03) # global epsilon for multi-message protocols


    #ms = [None] # general LDP mechanisms
    ms = [None, "laplace", "piecewise", "krr", "subset", "localhash", "hardamard", "hardamardB", "collision"]
    #ms = [None, "collision"]
    #ms = ["krr_para_advanced", "krr_para_basic", "krr_sep_best", "krr_sep_worst"] # range queries with krr
    #ms = ["cheu", "balls2bins"] # multi-message protocols

    #bounds = ["variation_ratio", "variation_ratio_closed", "variation_ratio_tightclosed", "strong_clone", "strong_clone_closed", "clone", "clone_closed",
    #"erlingsson", "erlingsson_closed", "blanket_hoeffding", "blanket_bennett", "blanket_hoeffding_generic", "blanket_bennett_generic"]
    bounds = ["variation_ratio", "variation_ratio_closed", "variation_ratio_tightclosed", "strong_clone", "strong_clone_closed", "clone_closed"]
    #bounds = ["variation_ratio", "variation_ratio_tightclosed", "variation_ratio_closed"]

    # the dict for recording amplification settings/results
    results = {}
    results["epsilons"] = epsilons.tolist()
    results["ms"] = ms

    n = 1000000
    nusers = n
    delta = 0.01/n
    # number of binary search iterations
    T = 20

    d = 100 # 64, 2048, domain size
    s = 8 # sparsity parameter for set-valued data or key-value data
    results["n"] = n
    results["delta"] = delta
    results["d"] = d
    print("epsilons, n, ms, d, s", epsilons, n, ms, (d,s))
    print("bounds", bounds)

    #filename = "unified_n"+str(n)+"_d"+str(d)+"_s"+str(s)+".json"
    filename = None  # None means don't save results to disk

    for epsilon in epsilons:
        # for single-message
        options = [None, None, None, d, (d, int(np.ceil(d/(np.exp(epsilon)+0)))), int(np.exp(epsilon)+1), (d, d/2), (d, d/2, 8), (s, int(s*np.exp(epsilon)+2*s-1))]
        #options = [None, s]

        # for range query with parallel composition on krr
        #options = [d, d, d, d]

        # for multi-message protocol
        #options = [(n, epsilon, delta), (d, s)] # options for multi-message protocols

        for mi, m in enumerate(ms):
            if m in ["balls2bins"]:
                d, s = options[mi][0], options[mi][1]
                n = int(32*np.log(2/delta)*d/s/(epsilon*epsilon)) # number of covering messages
                print(m, n, d, s)
            elif m in ['krr_sep_best', 'krr_sep_worst']:
                d = int(options[mi])
                n = int(n/np.log2(d)) # separating into sub-populations for each hierarchy
            else:
                n = nusers

            if results.get(m) is None:
                results[m] = {}
                for bound in bounds:
                    results[m][bound] = []
            if "variation_ratio" in bounds:
                p, beta, q0, q1 = computeUBParameters(epsilon=epsilon, mechanism=ms[mi], options=options[mi])
                start = process_time()
                ub = amplificationUB(p, beta, q0, q1, n, delta, T)
                end = process_time()

                results[m]["variation_ratio"].append(ub)
                #print running time
                print("running time: m, epsilon, n, process time", m, epsilon, n, end-start)
                #results[m]["time_variation_ratio"].append(end-start)

            if "variation_ratio_closed" in bounds:
                p, beta, q0, q1 = computeUBParameters(epsilon=epsilon, mechanism=ms[mi], options=options[mi])
                q = max(q0, q1)
                r = beta*p/((p-1)*q)
                aub1 = epsilon
                if n > 8*np.log(4/delta)/r:
                    c = max(0, 4*(1-3*r)/(9*(1-2*r)))
                    aub1 = np.log(1+beta/((1-c)*(1+p)*beta/(p-1)+c)*(np.sqrt(32*np.log(4/delta)/(r*(n-1)))+4/(r*(n-1))))
                results[m]["variation_ratio_closed"].append(aub1)

            if "variation_ratio_tightclosed" in bounds:
                p, beta, q0, q1 = computeUBParameters(epsilon=epsilon, mechanism=ms[mi], options=options[mi])
                q = max(q0, q1)
                alpha = beta/(p-1)
                r = alpha*p/q
                aubx = epsilon
                Omega = (n-1)*(2*r)-np.sqrt(min(6*r,1/2)*(n-1)*np.log(4/delta))
                if Omega > 0 and (p+1)*alpha/2-(1-alpha-alpha*p)*r/(1-2*r) > 0.0 and Omega > (2*p*(beta+1+(beta-1)*p)*(n-1)+beta)/(q+p*(beta-1+(beta+1)*p)-p*q):
                    alpha = beta/(p-1)
                    aubx = np.log(1+beta*(2*np.sqrt(Omega*np.log(4/delta)/2)+1)/(alpha*Omega+beta*(Omega/2-np.sqrt(Omega*np.log(4/delta)/2))+(1-alpha-alpha*p)*(n-1-Omega)*r/(1-2*r)))
                results[m]["variation_ratio_tightclosed"].append(aubx)

            if "strong_clone" in bounds:
                p = np.exp(epsilon)
                alpha = 1.0/(np.exp(epsilon)+1)
                r = 1.0/(np.exp(epsilon)+1)
                ub_sc = amplificationUBCore(p, alpha, r, r, n, delta, T)
                results[m]["strong_clone"].append(ub_sc)

            if "strong_clone_closed" in bounds:
                alpha = 1.0/(np.exp(epsilon)+1)
                r = 1.0/(np.exp(epsilon)+1)
                aub2 = epsilon
                if n > 8*np.log(4/delta)/r:
                    aub2 = np.log(1+(np.exp(epsilon)-1)/(np.exp(epsilon)+1)*(np.sqrt(32*np.log(4/delta)/(r*(n-1)))+4/(r*(n-1))))
                results[m]["strong_clone_closed"].append(aub2)

            # clone reduction, Vitaly Feldman, Audra McMillan, and Kunal Talwar. Hiding among the clones: A simple and nearly optimal analysis of privacy amplification by shuffling. In 2021 IEEE 62nd Annual Symposium on Foundations of Computer Science (FOCS), pages 954–964. IEEE, 202
            # for comparison with numerical clone, download https://github.com/apple/ml-shuffling-amplification/
            """
            if "clone" in bounds:
                import computeamplification as CA
                r = 1.0/(2*np.exp(epsilon))
                ub_c = CA.numericalanalysis(n, 2*r, epsilon, delta, T, 10, True, coin=1/2, factor=1.0)[1]
                results[m]["clone"].append(ub_c)
            """

            if "clone_closed" in bounds:
                r = 1.0/(2*np.exp(epsilon))
                aub3 = epsilon
                if n > 8*np.log(4/delta)/r:
                    aub3 = np.log(1+(np.exp(epsilon)-1)/(np.exp(epsilon)+1)*(np.sqrt(32*np.log(4/delta)/(r*(n-1)))+4/(r*(n-1))))
                results[m]["clone_closed"].append(aub3)

            # privacy blanket, B. Balle, J. Bell, A. Gascon, and K. Nissim. The Privacy Blanket of the Shuffle Model, International Cryptology Conference (CRYPTO), 2019
            # for comparison with numerical privacy blanket, download https://github.com/BorjaBalle/amplification-by-shuffling
            """
            import shuffleddp
            from shuffleddp.amplification_bounds import *
            from shuffleddp.mechanisms import *
            if shuffleddp is not None and len({"erlingsson_closed", "blanket_hoeffding", "blanket_bennett", "blanket_hoeffding_generic", "blanket_bennett_generic"}.intersection(set(bounds))) > 0:
                blanketname = m
                if m in ["collision"]:
                    blanketname += str(options[mi])
                mks = [LDPMechanism(eps0=epsilon, name=blanketname), LDPMechanism(eps0=epsilon, name="generic")]

                # This bound does not use infomation about the underlying mechanism
                erlingsson = Erlingsson()
                all_bounds = [erlingsson]

                # The other two bounds can use a 'Generic' mechanism or a specific mechanism
                bound_types = [Hoeffding, BennettExact]
                for mk in mks:
                    for B in bound_types:
                        all_bounds.append(B(mk))
                epcbounds = []
                for bi, b in enumerate(all_bounds):
                    epc = b.get_eps(epsilon, n, delta)
                    epcbounds.append(epsilon if math.isnan(epc) else epc)

                results[m]["erlingsson_closed"].append(epcbounds[0])
                results[m]["blanket_hoeffding"].append(epcbounds[1])
                results[m]["blanket_bennett"].append(epcbounds[2])
                results[m]["blanket_hoeffding_generic"].append(epcbounds[3])
                results[m]["blanket_bennett_generic"].append(epcbounds[4])
            """

            print(epsilon, m, [results[m][bound][-1] for bound in bounds],  options[mi])

    if filename is not None:
        with open(datetime.now().isoformat().replace(':', '_')+'-'+filename, 'w') as outfile:
            json.dump(results, outfile)
