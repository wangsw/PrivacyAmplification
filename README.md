# Privacy Amplification via Shuffling

Unified, Tight, and Fast Privacy Amplification in the Shuffle Model of Differential Privacy

- To compute tight privacy amplification upper bound, call:

    - **amplificationUB(p, beta, q0, q1, n, delta, T, step)**

    - for most local randomizers, q0=q1=q, except 2-RR and exponential mechanism for metric privacy


- To compute tight privacy amplification lower bound, call:

    - **amplificationLB(p, beta, q0, q1, n, delta, T, step)**
    
    
    
   
- List of amplification parameters for LDP randomizers



| **randomizer**                                                           | **parameter $p$** | **parameter $\beta$**                                                                                                                                  | **parameter $q$** |
|------------------------------------------------------------------------------------|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| general mechanisms                                                                 | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+1}$                                                                                                                          | $e^{\epsilon}$              |
| Laplace mechanism for $[0,1]$ \cite{dwork2006calibrating}                          | $e^{\epsilon}$              | $1-e^{-\epsilon/2}$                                                                                                                                              | $e^{\epsilon}$              |
| Duchi \textit{et al.} \cite{duchi2013local}  for $[-1,1]^d$                        | $e^{\epsilon}$              | $\frac{e^{\epsilon}-1}{e^{\epsilon}+1}$                                                                                                                          | $e^{\epsilon}$              |
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
::: {#tab:parameterslmdp}
                             **randomizer**                             **parameter $p$**                                                         **parameter $\beta$**                                                         **parameter $q$**
  -------------------------------------------------------------------- ------------------- ----------------------------------------------------------------------------------------------------------------------------------- -------------------
                           general mechanisms                             $e^{d_{01}}$                                                     $\frac{e^{d_{01}}-1}{e^{d_{01}}+1}$                                                    $e^{d_{max}}$
   Laplace mechanism [@alvim2018local], $\ell_1$-norm on $\mathbb{R}$     $e^{d_{01}}$                                                              $1-e^{-d_{01}/2}$                                                             $e^{d_{max}}$
    planar Laplace [@andres2013geo], $\ell_2$-norm on $\mathbb{R}^2$      $e^{d_{01}}$      $2\int_{0}^{\frac{d_{01}}{2}}\int_{-\infty}^{+\infty}\frac{e^{-\sqrt{(x-\frac{d_{01}}{2})^2+y^2}}}{2\pi}\mathrm{d} y \mathrm{d}x$     $e^{d_{max}}$
   $(B,m,F)$-WitchHat [@wang2023SMDP], $\ell_1$-norm on $\mathbb{R}$      $e^{d_{01}}$                                      $\frac{2(e^{m}-e^{{d_{01}}/{F}}+{d_{01}}/{F}-m)}{F(e^{m}-1)+2B}$                                      $e^{d_{max}}$

  : A summary of amplification parameters of
  $\mech{S}(\mech{R}(x_1^0),..., \mech{R}(x_n))$ and
  $\mech{S}(\mech{R}(x_1^1),..., \mech{R}(x_n)))$ for local
  $d_\dom{X}$-DP randomizers, $d_{01}=d_\dom{X}(x_1^0,x_1^1)$ and
  $d_{max}=\max_{x\in \dom{X}} \max\{d_\dom{X}(x_1^0,x), d_\dom{X}(x_1^1,x)\}$.
:::


- List of amplification parameters for multi-message local randomizers
::: {#tab:parametersmulti}
                                             **randomizer**                                               **parameter $p$**     **parameter $\beta$**   **parameter $q$**
  ---------------------------------------------------------------------------------------------------- ----------------------- ----------------------- -------------------
   secret shares over group $\mathbb{Z}_{l}$ [@ghazi2019scalable; @ghazi2020fewer; @balle2020private]         $+\infty$                  $1$                   $l$
     Cheu *et al.* [@cheu2022differentially] with flip prob. $f\in [0,0.5]$ on $\{0,1\}^d$ vectors      $\frac{(1-f)^2}{f^2}$          $1-2f$            $\frac{1-f}{f}$
                Balls-into-bins [@luo2022frequency] with $d$ bins (and $s$ special bins)                      $+\infty$                  $1$              $\frac{d}{s}$

  : A summary of amplification parameters of some multi-message shuffle
  private protocols.
:::
