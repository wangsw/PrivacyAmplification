# Privacy Amplification via Shuffling

**Unified, Tight and Fast Privacy Amplification in the Shuffle Model of Differential Privacy**

&check; LDP randomizers, metric LDP randomizers, and some multi-message randomizers

&check; closed-form bounds and numerical bounds. The numerical bound takes 20 seconds for n=1,000,000 and 2 minutes for n=100,000,000. For better efficiency (but less tighter bound), use large tolerance (e.g., tol=\delta/n) or small number of iterations (e.g., T=10) 

&check; tight parallel composition

&check; tight sequential composition 

## Usage

- To compute tight privacy amplification upper bound, call:

    - **amplificationUB(p, beta, q0, q1, n, delta, T)**

    - for most local randomizers, q0=q1=q


- To compute tight privacy amplification lower bound, call:

    - **amplificationLB(p, beta, q0, q1, n, delta, T)**
    
    
    
   
- List of amplification parameters for LDP randomizers



| **randomizer**                                                           | **parameter $p$** | **parameter $\beta$**                                                                                                                                  | **parameter $q$** |
|------------------------------------------------------------------------------------|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| general mechanisms                                                                 | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+1}$                                                                                                                          | $e^{\epsilon}$              |
| Laplace mechanism for $[0,1]$ \cite{dwork2006calibrating}                          | $e^{\epsilon}$              | $1-e^{-\epsilon/2}$                                                                                                                                              | $e^{\epsilon}$              |
| Duchi et al. \cite{duchi2013local}  for $[-1,1]^d$                        | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+1}$                                                                                                                          | $e^{\epsilon}$              |
| Harmony mechanism for $[-1,1]^d$ \cite{nguyen2016collecting}                       | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+1}$                                                                                                                          | $e^{\epsilon}$              |
| Piecewise mechanism for $[-1,1]$ \cite{wang2019collecting}                         | $e^{\epsilon}$              | $1-e^{-\epsilon/2}$                                                                                                                                              | $e^{\epsilon}$              |
| randomized response (RR) on $\{0,1\}$ \cite{warner1965randomized}                  | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+1}$                                                                                                                          | $e^{\epsilon}$              |
| general RR on $d$ options \cite{kairouz2016discrete}                               | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+d-1}$                                                                                                                        | $e^{\epsilon}$              |
| binary RR (RAPPOR) on $d$ options \cite{duchi2013local,erlingsson2014rappor}       | $e^{\epsilon}$              | $\frac{e^{\epsilon/2}-1}{e^{\epsilon/2}+1}$                                                                                                                      | $e^{\epsilon}$              |
| $k$-subset mechanism on $d$ options \cite{wang2019local,ye2018optimal}             | $e^{\epsilon}$              | $\frac{(e^{\epsilon}-1)({d-1\choose k-1}-{d-2\choose k-2})}{e^{\epsilon}{d-1\choose k-1}+{d-1\choose k}}$                                                        | $e^{\epsilon}$              |
| local hash with length $l$ \cite{wang2017locally}                                  | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+l-1}$                                                                                                                        | $e^{\epsilon}$              |
| Hadamard response $(K,s)$ \cite{acharya2018hadamard}                               | $e^{\epsilon}$              | $\frac{s(e^{\epsilon}-1)/2}{s e^{\epsilon}+K-s}$                                                                                                                 | $e^{\epsilon}$              |
| Hadamard response $(K,s,B>1)$ \cite{acharya2018hadamard}                           | $e^{\epsilon}$              | $\frac{s(e^{\epsilon}-1)}{s e^{\epsilon}+K-s}$                                                                                                                   | $e^{\epsilon}$              |
| sampling RAPPOR on $s$ in $d$ options \cite{qin2016heavy}                          | $e^{\epsilon}$              | $\frac{s(e^{\epsilon/2}-1)}{d(e^{\epsilon/2}+1)}$                                                                                                                | $e^{\epsilon}$              |
| $k$-subset exponential on $s$ in $d$ options \cite{wang2018privset}                | $e^{\epsilon}$              | $\frac{(e^{\epsilon}-1)({d-s\choose k}-{d-2s\choose k})}{e^{\epsilon}({d\choose k}-{d-s\choose k})+{d-s\choose k}}$                                              | $e^{\epsilon}$              |
| PrivKV on $s$ in $d$ keys ($\epsilon_1+\epsilon_2=\epsilon$) \cite{Ye2019PrivKVKD} | $e^{\epsilon}$              | $\frac{2s\max\{\frac{e^{\epsilon_1}(e^{\epsilon_2}-1)}{e^{\epsilon_2}+1}, e^{\epsilon_1}-1+\frac{e^{\epsilon_2}-1}{2(e^{\epsilon_2}+1)}\}}{d(e^{\epsilon_1}+1)}$ | $e^{\epsilon}$              |
| PCKV-GRR on $s$ in $d$ keys \cite{gu2019pckv}                                      | $e^{\epsilon}$              | $\frac{s(e^{\epsilon}-1)}{s e^{\epsilon}+2d-s}$                                                                                                                  | $e^{\epsilon}$              |
| $(d,s)$-Collision with length $l$ \cite{wang2021hiding}                            | $e^{\epsilon}$              | $\frac{\min\{s,l-s\}(e^{\epsilon}-1)}{s e^{\epsilon}+l-s}$                                                                                                       | $e^{\epsilon}$              |


- List of amplification parameters for Local metric DP randomizers

| **randomizer**                                               | **parameter $p$** | **parameter $\beta$**                                                                                                   | **parameter $q$** |
|------------------------------------------------------------------------|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| general mechanisms                                                     | $e^{d_{01}}$                | $\frac{e^{d_{01}}-1}{e^{d_{01}}+1}$                                                                                               | $e^{d_{max}}$               |
| Laplace mechanism \cite{alvim2018local}, $\ell_1$-norm on $\mathbb{R}$ | $e^{d_{01}}$                | $1-e^{-d_{01}/2}$                                                                                                                 | $e^{d_{max}}$               |
| planar Laplace \cite{andres2013geo}, $\ell_2$-norm on $\mathbb{R}^2$   | $e^{d_{01}}$                | $2\int_{0}^{\frac{d_{01}}{2}}\int_{-\infty}^{+\infty}\frac{e^{-\sqrt{(x-\frac{d_{01}}{2})^2+y^2}}}{2\pi}\mathrm{d} y \mathrm{d}x$ | $e^{d_{max}}$               |
| $(B,m,F)$-WitchHat \cite{wang2023SMDP}, $\ell_1$-norm on $\mathbb{R}$  | $e^{d_{01}}$                | $\frac{2(e^{m}-e^{{d_{01}}/{F}}+{d_{01}}/{F}-m)}{F(e^{m}-1)+2B}$                                                                  | $e^{d_{max}}$               |



- List of amplification parameters for multi-message local randomizers

| **randomizer**                                                                                 | **parameter $p$** | **parameter $\beta$** | **parameter $q$** |
|----------------------------------------------------------------------------------------------------------|-----------------------------|---------------------------------|-----------------------------|
| Balcer et al. \cite{balcer2020separating} with coin prob. $p$ for binary summation             | $+\infty$          | $1$  |  $\max\{\frac{1}{p},\frac{1}{1-p}\}$ |
| Balcer et al. \cite{balcer2021connecting} with uniform coin for binary summation               | $+\infty$          | $1$  |  $2$ |
| Cheu et al. \cite{cheu2022differentially} with flip prob. $f\in [0,0.5]$ on $\{0,1\}^d$ vectors | $\frac{(1-f)^2}{f^2}$       | $1-2f$                          | $\frac{1-f}{f}$             |
| Balls-into-bins \cite{luo2022frequency} with $d$ bins (and $s$ special bins)                             | $+\infty$                   | $1$                             | $\frac{d}{s}$               |


