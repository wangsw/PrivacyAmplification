"""
Fourier composition analyses for generalized variation-ratio reduction in the shuffle model,
based on the code of https://github.com/DPBayes/numerical-shuffler-experiments
"""

import numpy as np
import scipy
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import time
import math
import itertools


#  Probabilities for the case a>0
def get_logP1(i, j, n, p, beta, q1, q2):

    a=j+1
    b=i-j

    alpha = beta / (p - 1)
    r1 = alpha * p / q1
    r2 = alpha * p / q2


    term0 = np.log((alpha*a/r1+p*alpha*b/r2+(1-alpha-alpha*p)*(n-a-b)/(1-r1-r2))/(a/r1))
    # Then term0 is multiplied by P(P_0=(a,b))
    term1 = scipy.special.gammaln(n) - scipy.special.gammaln(i + 1) - scipy.special.gammaln(n - i)
    term2 = scipy.special.gammaln(i + 1) - scipy.special.gammaln(j + 1) - scipy.special.gammaln(i - j + 1)

    return term0 + term1 + term2 + i*np.log(r1+r2)+(n-1-i)*np.log(1-r1-r2)+j*np.log(r1/(r1+r2))+(i-j)*np.log(r2/(r1+r2))


def get_logP2(i, j, n, p, beta, q1, q2):
    a = j + 1
    b = i - j

    alpha = beta/(p-1)
    r1 = alpha*p/q1
    r2 = alpha*p/q2

    term0 = np.log((p*alpha*a/r1+alpha*b/r2+(1-alpha-alpha*p)*(n-a-b)/(1-r1-r2))/(a/r1))
    # Then term0 is multiplied by P(P_0=(a,b))
    term1 = scipy.special.gammaln(n) - scipy.special.gammaln(i + 1) - scipy.special.gammaln(n - i)
    term2 = scipy.special.gammaln(i + 1) - scipy.special.gammaln(j + 1) - scipy.special.gammaln(i - j + 1)

    return term0 + term1 + term2 + i*np.log(r1+r2)+(n-1-i)*np.log(1-r1-r2)+j*np.log(r1/(r1+r2))+(i-j)*np.log(r2/(r1+r2))


def get_omega(n, k, p, beta, q1, q2, nx, L):

    P1 = []
    Lx = []

    alpha = beta / (p - 1)
    r1 = alpha * p / q1
    r2 = alpha * p / q2
    r = r1+r2


    tol_ = 60
    lower_i=int(max(0,np.floor((n-1)*(r-np.sqrt(tol_/(2*(n-1)))))))
    upper_i=int(min(n-1,np.ceil((n-1)*(r+np.sqrt(tol_/(2*(n-1)))))))

    #print("i", lower_i, upper_i)

    for i in range(lower_i, upper_i+1):

        lower_j=int(max(0,np.floor(i*(r1/(r1+r2)-np.sqrt(tol_/(2*(i+1)))))))
        upper_j=int(min(i,np.ceil(i*(r1/(r1+r2)+np.sqrt(tol_/(2*(i+1)))))))

        #for the cases b>0, a>0
        #print("j", (i, lower_j, upper_j), (alpha, r1, r2, r1+r2))

        for j in range(lower_j, upper_j+1):
            p_temp=get_logP2(i, j, n, p, beta, q1, q2)
            P1.append(p_temp)
            a=j+1
            b=i-j
            nom= p*alpha*a/r1+alpha*b/r2+(1-alpha-alpha*p)*(n-a-b)/(1-r1-r2)
            denom= alpha*a/r1+p*alpha*b/r2+(1-alpha-alpha*p)*(n-a-b)/(1-r1-r2)
            Lx.append(np.log(nom/denom))


    dx=2*L/nx
    grid_x=np.linspace(-L,L-dx,nx)
    omega_y=np.zeros(nx)
    len_lx=len(Lx)
    for i in range(0,len_lx):
        lx=Lx[i]
        if lx>-L and lx<L:
            # place the mass to the right end point of the interval -> upper bound for delta
            ii = min(nx-1, int(np.ceil((lx+L)/dx)))
            omega_y[ii]=omega_y[ii]+np.exp(P1[i])
        else:
            print('\tOUTSIDE OF [-L,L] - RANGE')

    # check the mass to be sure
    print('\tSum of Probabilities : ' + str(sum(omega_y)))



    # Compute the periodisation & truncation error term,
    # using the error analysis given in
    # Koskela, Antti, et al.
    # Tight differential privacy for discrete-valued mechanisms and
    # for the subsampled gaussian mechanism using fft.
    # International Conference on Artificial Intelligence and Statistics. PMLR, 2021.

    # the parameter lambda in the error bound can be chosen freely.
    # commonly the error term becomes negligible for reasonable values of L.
    # here lambda is just hand tuned to make the error term small.

    lambd=0.5*L
    lambda_sum_minus=0

    for i in range(len(Lx)):
        plf = -Lx[i]
        lambda_sum_minus+=np.exp(lambd*plf + P1[i])

    alpha_minus=np.log(lambda_sum_minus)

    print("\t", alpha_minus)

    #Then alpha plus
    lambda_sum_plus=0
    for i in range(len(Lx)):
        plf = -Lx[i]
        lambda_sum_plus+=np.exp(lambd*plf + P1[i])

    alpha_plus=np.log(lambda_sum_plus)

    #Evaluate the periodisation error bound using alpha_plus and alpha_minus
    T1=(np.exp(k*alpha_plus))
    T2=(np.exp(k*alpha_minus))
    error_term = (T1+T2)*(np.exp(-lambd*L)/(1-np.exp(-2*lambd*L)))
    print('\terror_term: ' + str(error_term))

    return omega_y, grid_x, error_term




def homogeneousComposition(n, p, beta, q1, q2, nx, L, epsilons, k):
    t1 = time.perf_counter()


    omega, grid_x, error_term = get_omega(n, 1, p, beta, q1, q2, nx, L)
    print('Time of forming PLD: ' + str(time.perf_counter()-t1))

    delta_table=[]
    deltas=[]

    for eps in epsilons:
        y = omega
        j_start = int(np.ceil((L+eps)/dx))

        delta_temp = np.sum((1-np.exp(eps-grid_x[j_start:]))*y[j_start:])
        deltas.append(delta_temp)
        print('eps: ' + str(eps) + ' delta : ' + str(delta_temp))

    delta_table.append(deltas)


    omega, grid_x, error_term = get_omega(n, 1, p, beta, q1, q2, nx, L)

    half = int(nx/2)

    fx=np.copy(omega)

    # Flip fx, i.e. fx <- D(fx), the matrix D = [0 I;I 0]
    temp = np.copy(fx[half:])
    fx[half:] = np.copy(fx[:half])
    fx[:half] = temp

    t1 = time.perf_counter()
    # Compute the DFT
    FF1 = np.fft.fft(fx)
    y = np.fft.ifft((FF1**k))
    print('Time of FFTs: ' + str(time.perf_counter()-t1))


    # Flip again, i.e. cfx <- D(cfx), D = [0 I;I 0]
    temp = np.copy(y[half:])
    y[half:] = y[:half]
    y[:half] = temp
    omega1=y

    deltas=[]

    for eps in epsilons:
        y = omega1
        j_start = int(np.ceil((L+eps)/dx))
        delta_temp = np.sum((1-np.exp(eps-grid_x[j_start:]))*y[j_start:])
        print('k:', k, 'eps: ' + str(eps) + ' delta : ' + str(delta_temp))
        deltas.append(delta_temp)

    #delta_table.append(deltas)



def heterougenousComposition(ns, ps, betas, q1s, q2s, nx, L, epsilons, startrepeat=None):
    # Heterougenous composition
    if startrepeat is None:
        startrepeat = len(ns)
    FF1s = []
    half = int(nx/2)
    for ni, n in enumerate(list(range(startrepeat))):
        omega, grid_x, error_term = get_omega(ns[ni], 1, ps[ni], betas[ni], q1s[ni], q2s[ni], nx, L)

        fx = np.copy(omega)

        # Flip fx, i.e. fx <- D(fx), the matrix D = [0 I;I 0]
        temp = np.copy(fx[half:])
        fx[half:] = np.copy(fx[:half])
        fx[:half] = temp

        # Compute the DFT
        FF1s.append(np.fft.fft(fx))

    for ni, n in enumerate(list(range(startrepeat, len(ns)))):
        FF1s.append(FF1s[ni%startrepeat])

    t1 = time.perf_counter()

    # save memory in np.prod(FF1s, axis=0), via a batch model
    batch = 32
    prods = []
    start = 0
    while start < len(FF1s):
        prods.append(np.prod(FF1s[start:min(start+batch, len(FF1s))], axis=0))
        start += batch
    y = np.fft.ifft(np.prod(prods, axis=0))
    print('Time of FFTs: ' + str(time.perf_counter()-t1))


    # Flip again, i.e. cfx <- D(cfx), D = [0 I;I 0]
    temp = np.copy(y[half:])
    y[half:] = y[:half]
    y[:half] = temp
    omega1=y

    deltas=[]

    for eps in epsilons:
        y = omega1
        j_start = int(np.ceil((L+eps)/dx))
        delta_temp = np.sum((1-np.exp(eps-grid_x[j_start:]))*y[j_start:])
        print('Heterogenous', len(ps), 'eps: ' + str(eps) + ' delta : ' + str(delta_temp))
        deltas.append(delta_temp)
    return deltas


# parameters 1
n=int(1E4)
p = np.exp(4.0)
beta = (p-1)/(p+1)
q1 = 1.001
q2 = p
n_eps = 30
epsilons = np.array(list(np.linspace(0.02, 0.5, 49))+list(np.linspace(0.52, 1.0, 25))+list(np.linspace(1.1, 4.0, n_eps)))


nx=int(1E7)
L=20
dx=2*L/nx
ks=[]

for k in ks:
    homogeneousComposition(n, p, beta, q1, q2, nx, L, epsilons, k)


#"""
print("-----multi-message-federated-learning-simulation-----")
def improved_tv_bound(epsilons):
    #import itertools
    d = len(epsilons)
    prod_term = np.prod([1 + np.exp(e) for e in epsilons])
    subsets = list(itertools.chain.from_iterable(itertools.combinations(range(d), r) for r in range(d+1)))

    sum_term = 0
    for V in subsets:
        sum_V = sum(epsilons[j] for j in V)
        sum_not_V = sum(epsilons[j] for j in range(d) if j not in V)
        term = max(np.exp(sum_V) - np.exp(sum_not_V), 0)
        sum_term += term

    total_variation = sum_term / prod_term
    return total_variation


def get_mms_inf_vector_params(n, s, m, v, T=20, reduction="advancedvaritionratio", reverse=False):
    # Multi-Message Shuffled Privacy in Federated Learning
    ns = []
    ps = []
    betas = []
    q1s = []
    q2s = []

    for t in range(T):
        denom = np.sum([4**(-mi/3) for mi in range(1, m)])+4**(-m/3+1/3)
        for mi in range(1, m+1):
            vmi = 4**(-mi/3)/denom*v if mi < m else 4**(-m/3+1/3)/denom*v
            pmi = (1-np.sqrt((vmi**2/s**2)/(vmi**2/s**2+4)))/2

            p = (1.0-pmi)/pmi
            beta = 1.0-2*pmi
            q1 = 1.0001
            q2 = p
            if s > 1:
                p = ((1.0-pmi)/pmi)**s
                beta = improved_tv_bound([np.log(p)/s]*s)
                q1 = p
                q2 = p
            if reduction in ["variationratio"]:
                p = ((1.0-pmi)/pmi)**s
                beta = improved_tv_bound([np.log(p)/s]*s)
                q1 = p
                q2 = p
            if reduction in ["strongerclone"]:
                p = ((1.0-pmi)/pmi)**s
                beta = (p-1)/(p+1)
                q1 = p
                q2 = p
            if reverse:
                temp = q1
                q1 = q2
                q2 = temp
            #print((n, s, m, v, t), (p, beta, q1, q2))
            ns.append(n)
            ps.append(p)
            betas.append(beta)
            q1s.append(q1)
            q2s.append(q2)
    return ns, ps, betas, q1s, q2s


n = 10000
s = 1
v = 4.0
m = int(np.ceil(v))
T = 1
epsilons = epsilons

print("settings", n, s, m, v, T)

print("-----stronger-clone-----")
ns, ps, betas, q1s, q2s = get_mms_inf_vector_params(n, s, m, v, T, "strongerclone")
heterougenousComposition(ns, ps, betas, q1s, q2s, nx, L, epsilons, startrepeat=len(ns)//T)
print("+++++stronger-clone+++++")

# print("-----variation-ratio-----")
# ns, ps, betas, q1s, q2s = get_mms_inf_vector_params(n, s, m, v, T, "variationratio")
# heterougenousComposition(ns, ps, betas, q1s, q2s, nx, L, epsilons, startrepeat=len(ns)//T)
# print("+++++variation-ratio+++++")

print("-----advanced-variation-ratio-----")
ns, ps, betas, q1s, q2s = get_mms_inf_vector_params(n, s, m, v, T, "advancedvariationratio")
heterougenousComposition(ns, ps, betas, q1s, q2s, nx, L, epsilons, startrepeat=len(ns)//T)

print("\t--reverse, due to insymmetry, take their maximum delta--")
ns, ps, betas, q1s, q2s = get_mms_inf_vector_params(n, s, m, v, T, "advancedvariationratio", reverse=True)
heterougenousComposition(ns, ps, betas, q1s, q2s, nx, L, epsilons, startrepeat=len(ns)//T)
print("\t++reverse, due to insymmetry, take their maximum delta++")

print("+++++advanced-variation-ratio+++++")

print("+++++multi-message-federated-learning-simulation+++++")
print("settings", n, s, m, v, T)

#"""

"""
print("-----multi-message-decision-trees-simulation-----")

def get_decision_trees_category_params(n, d, f, h, t, m, withdelta=False):
    ns = []
    ps = []
    betas = []
    q1s = []
    q2s = []
    deltaratios = []
    for ti in range(t):
        for hi in range(h):
            p = np.exp(20.0)
            beta = (p-1)/(p+1)
            q1 = f**(hi+1)
            q2 = q1

            nh = d-hi # number of histogram queries in the h-th level
            ns += [n*(m-1)+1]*nh
            ps += [p]*nh
            betas += [beta]*nh
            q1s += [q1]*nh
            q2s += [q2]*nh
            if withdelta is not False:
                deltaratios += [1.0]*nh #[np.exp(q1/2)/nh]*nh
    if not withdelta:
        return ns, ps, betas, q1s, q2s
    else:
        return ns, ps, betas, q1s, q2s, list(np.array(deltaratios)/np.sum(deltaratios))

n = 10000
d = 10
f = 2 # fan out
h = 5 # tree height
t = 1  # number of trees
m = 4 # number of messages per user

target_delta = 0.01/n


print("settings", n, d, f, h, t, m, target_delta)

print("-----advanced-composition-----")
from unifiedamplification import amplificationUB
ns, ps, betas, q1s, q2s, deltaratios = get_decision_trees_category_params(n, d, f, h, t, m, True)
nq = len(ns) # total number of queries

#print(ns, q1s)

epss = []
deltas = [target_delta/nq]*nq
for i in range(nq):
    oneround_epsilon = amplificationUB(ps[i], betas[i], q1s[i], q2s[i], ns[i], deltas[i], 15)
    epss.append(oneround_epsilon)
print("basic composition", (np.sum(epss), np.sum(deltas)), epss, deltas, q1s)



epss = []
deltas = list(np.array(deltaratios)*target_delta*0.5)
for i in range(nq):
    oneround_epsilon = amplificationUB(ps[i], betas[i], q1s[i], q2s[i], ns[i], deltas[i], 15)
    epss.append(oneround_epsilon)
print("advanced composition", (
                            np.sqrt(2*len(epss)*np.log(1/(target_delta-np.sum(deltas))))*np.min(epss)+len(epss)*np.min(epss)*(np.exp(np.min(epss))-1),
                            np.sqrt(2*len(epss)*np.log(1/(target_delta-np.sum(deltas))))*np.max(epss)+len(epss)*np.max(epss)*(np.exp(np.max(epss))-1),
                            target_delta), epss, deltas)
# np.min underestimates privacy, for just comparision, do not use in real world


# use heterogeneous composition
def compute_LHS(epsilon_js, epsilon_prime):
    m = len(epsilon_js)
    n = m  # Number of variables

    # Compute log_D = sum_{j=1}^{m} log(1 + e^{epsilon_j})
    log_D = sum(math.log1p(math.exp(epsilon_j)) for epsilon_j in epsilon_js)
    total_eps = sum(epsilon_js)

    sum_DV = 0.0

    # Iterate over all subsets V of [d]
    for bits in range(1 << n):
        S_V = sum(epsilon_js[j] for j in range(n) if (bits >> j) & 1)
        S_not_V = total_eps - S_V

        log_DV1 = S_V
        log_DV2 = epsilon_prime + S_not_V
        max_log = max(log_DV1, log_DV2)

        # Compute D_V = e^{S_V} - e^{epsilon' + S_not_V} safely
        exp_diff = math.exp(log_DV1 - max_log) - math.exp(log_DV2 - max_log)
        D_V = math.exp(max_log) * exp_diff if exp_diff > 0 else 0

        if D_V > 0:
            sum_DV += D_V

    LHS = math.exp(-log_D) * sum_DV
    return LHS

def compute_RHS(delta, delta_js):
    denom = math.prod(1 - delta_j for delta_j in delta_js)
    RHS = 1 - (1 - delta) / denom
    return RHS

def find_min_epsilon_prime(epsilon_js, delta_js, delta, tol=5e-2):
    m = len(epsilon_js)
    if m >= 20:
        raise ValueError("The number of epsilon_js must be less than 20 due to computational limitations.")

    lower = 0
    upper = sum(epsilon_js)*2  # Initial upper bound

    RHS = compute_RHS(delta, delta_js)

    while upper - lower > tol:
        mid = (lower + upper) / 2
        LHS = compute_LHS(epsilon_js, mid)

        if LHS <= RHS:
            upper = mid
        else:
            lower = mid
    return upper

# maxcomp = 10  # maximum number of composition, as it requires exponential computation costs
# start = 0
# currentep = 0.0
# currentdelta = 0.0
# while start < len(epss):
#     tempepss = [currentep]+epss[start:min(len(epss)+1, start+maxcomp)]
#     tempdeltas = [currentdelta]+deltas[start:min(len(epss)+1, start+maxcomp)]
#     tempdelta = (target_delta-currentdelta)/(np.ceil((len(epss)-start)/maxcomp))
#     currentep = find_min_epsilon_prime(tempepss, tempdeltas, currentdelta+tempdelta)
#     currentdelta += tempdelta
#     start += maxcomp
#     print((len(epss), start), currentep, currentdelta)
# print("heterougeneous advanced composition", (currentep, currentdelta), target_delta)


print("+++++advanced-composition+++++")


print("-----advanced-variation-ratio-----")
ns, ps, betas, q1s, q2s = get_decision_trees_category_params(n, d, f, h, t, m)
heterougenousComposition(ns, ps, betas, q1s, q2s, nx, L, epsilons, startrepeat=len(ns)//t)
print("+++++advanced-variation-ratio+++++")



print("+++++multi-message-decision-trees-simulation+++++")
print("settings", n, d, f, h, t, m, target_delta)
"""

#delta_table.append(deltas)
#pickle.dump(delta_table, open("./pickles/deltas_pld_gvrfa.p", "wb"))
