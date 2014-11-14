% title: omnia.md: Engineering a Full Python Stack for Biophysical Computation
% author: Kyle A. Beauchamp

---
title: Moore's Law



<center>
<img height=550 src=figures/moore_law.png />
</center>

<footer class="source"> 
http://en.wikipedia.org/wiki/Moore%27s_law
</footer>

---
title: Eroom's Law

<center>
<img height=450 src=figures/eroom.png />
</center>



<footer class="source"> 
Scannell, 2012
</footer>
---
title: How much do drugs cost?  $2,000,000

<center>
<img height=450 src=figures/cost_structure.jpg />
</center>

<footer class="source"> 
Paul, 2010
</footer>

---
title: Molecular Dynamics: Easy to sample, Hard to Normalize

We mostly understand the microscopic rules governing proteins and ligands:

$$ U(x) =  U_{bonds}(x) +  U_{angles}(x) +  U_{dihedrals}(x) +  U_{LJ}(x) + U_{electrostatic}(x)$$
$$U_{bonds} = \sum_{i} k_i (r - r_0)^2$$
$$U_{angles} = \sum_{i} k_i (\theta - \theta_0)^2$$
$$U_{dihedrals} = \sum_{i} k_i (1 + \cos(n_i \phi_i + \delta_i))$$
$$ U_{electrostatic} = \sum_{i > j} \frac{q_i q_j} {4 \pi \epsilon_0 r_{ij}} $$

---
title: Some notation

$$q_k(x) = \exp(-u_k(x))$$
$$Z_k = c_k = \int q_k(x) dx$$

$$P_k(x) = \frac{1}{c_k} \exp(-u_k(x)) = \frac{1}{c_k} q_k(x_n)$$

$$\langle A(x) \rangle_k = \frac{1}{c_k} \int A(x) q_k(x) dx$$

Free energy:

$$f_j - f_i = -\log \frac{c_j}{c_i}$$


---
title: Applications of Normalizing Constants

- Binding free energy / other chemistry applications
- Bayesian model comparison
- Monte Carlo integration

---
title: Alchemical Free Energy Calculation

Key Idea: Binding free energy is a log ratio of normalizing constants

<center>
<img height=325 src=figures/T4-Phenol.png />
<img height=325 src=figures/T4-Toluene.png />
</center>

<footer class="source"> 
alchemistry.org
</footer>


---
title: How to estimate normalizing constants of unnormalized distributions?
class: segue dark nobackground


---
title: Exponential Averaging (EXP)

Given unnormalized densities $q_1(x)$ and $q_2(x)$, notice the following "Zwanzig" identity:

$$\langle \frac{q_2(x)}{q_1(x)} \rangle_1 = \frac{1}{c_1} \int \frac{q_2(x)}{q_1(x)} q_1(x) dx = \frac{1}{c_1} \int q_2(x) dx = \frac{c_2}{c_1}$$

Given samples $x_n$ from $q_1(x)$, we can estimate the ratio $\frac{c_2}{c_1}$

$$ \frac{c_2}{c_1} = \langle \frac{q_2(x)}{q_1(x)} \rangle_1 \approx \frac{1}{n} \sum_n \frac{q_2(x_n)}{q_1(x_n)}$$ 

EXP has large bias and variance due to heavy weights in the tails.  

<footer class="source"> 
Zwanzig, 1954
</footer>

---
title: Bennet Acceptance Ratio (BAR)

EXP can be improved by using samples from $q_1(x)$ and $q_2(x)$.  Requires solving nonlinear equations:

$$0 = \sum_n^{N_1} [1 + \frac{N_1}{N_2} c \frac{q_1(x_n)}{q_2(x_n)}]^{-1} + \sum_n^{N_2} [1 + \frac{N_2}{N_1} \frac{q_2(x_n)}{q_1(x_n)} c^{-1}]^{-1}$$

BAR is the optimal estimator of its type--minimum variance and asymptotically unbiased.

Question: Is there an equivalent procedure using data from multiple states?


<footer class="source"> 
Bennett, 1975
</footer>

---
title: Multistate Bennett Acceptance Ratio (MBAR)


<center>
<img height=150 src=figures/mbar_title.png />

<img height=275 src=figures/mbar_figure.png />
</center>

- Asymptotically Unbiased
- Minimum variance among class of estimators


---
title: Deriving MBAR: Ingredients

- $\{x_n\}_{n=1}^{N}$ are independent samples drawn from states $s_n$
- $s_n$ known and fixed
- $K$ states: $s_n \in \{1, 2, ..., K\}$
- Unnormalized densities $q_k(x_n)$ are known for all $k$, $n$


---
title: Deriving MBAR: Approach

We want a model for computing arbitrary expectations via "quadrature":

$$\langle A(x) \rangle_k = \frac{1}{c_k} \int q_k(x)A(x) dx = \frac{1}{c_k} \sum_n \rho_n q_k(x_n) A(x_n)$$

- Requires introducing quadrature weights $\rho_n$ (masses)
- Nonparametric model, empirical measure, finite support $\{x_n\}$


---
title: Deriving MBAR: Properties

$$\langle A(x) \rangle_k = \frac{1}{c_k} \int q_k(x)A(x) dx = \frac{1}{c_k} \sum_n \rho_n q_k(x_n) A(x_n)$$

Notice that $1 = \langle 1 \rangle_k$ implies that

$$c_k = \sum_n \rho_n q_k(x_n)$$

Second, notice that the conditional probabilities pick up the quadrature weights due to the finite support:

$$p(x_n|s_n, \{c_k\}, \{\rho_n\}) = \langle \delta(x - x_n) \rangle_k = \frac{1}{c_{s_n}} \rho_n q_{s_n}(x_n)$$


---
title: Deriving MBAR: Likelihood

The likelihood of our dataset is given by a product:

$$\prod_n^N P(x_n|s_n, \{c_k\}, \{\rho_n\}) =  \prod_n \frac{q_{s_n}(x_n) \rho_n}{c_{s_n}}$$

Drop $q_k(x_n)$, as it is does not depend of the parameters $\rho_n$ or $c_k$:

$$\prod_n \frac{q_{s_n}(x_n) \rho_n}{c_{s_n}} \propto \prod_n \frac{\rho_n}{c_{s_n}} = \prod_n \rho_n \prod_n c_{s_n}^{-1}$$

<footer class="source"> 
Zan, 2000.  
Bartels, 2000.
Vardi, 1985.
Gelman, 1996.
</footer>


---
title: Deriving MBAR: Likelihood

$$\prod_n P(x_n|s_n, \{c_k\}, \{\rho_n\}) \propto \prod_n \rho_n \prod_n c_{s_n}^{-1}$$

Count and collect the normalizing constants:

$$\prod_n P(x_n|s_n, \{c_k\}, \{\rho_n\}) = \prod_n^N \rho_n \prod_k^K c_{k}^{-N_k}$$

Note that the state origin $s_n$ is in the likelihood ONLY via $N_k$!  Finally, take the log:

$$LL = \sum_n^N \log \rho_n - \sum_k^K N_k \log c_k$$

---
title: Deriving MBAR: MLE

Let's take the partial derivative of the log likelihood:

$$\frac{\partial LL}{\partial \rho_n} = \frac{1}{\rho_n} - \sum_k \frac{N_k}{c_k} \frac{\partial c_k}{\partial \rho_n}$$

From $c_i = \sum_n q_i(x_n) \rho_n$, we know that $\frac{\partial c_k}{\partial \rho_n} = q_i(x_n)$

$$\frac{\partial LL}{\partial \rho_n} = \frac{1}{\rho_n} - \sum_k \frac{N_k}{c_k} q_k(x_n) = 0$$

$$\rho_n = [\sum_k N_k c_k^{-1} q_k(x_n)]^{-1}$$

---
title: MLE solution to MBAR

Plugging the MLE $\rho_n$ into $c_k$, we recover the self-consistent equation from the paper:

$$c_i = \sum_n^N \frac{q_i(x_n)}{\sum_k^K N_k c_k^{-1} q_k(x_n)}$$

Plugging the MLE of $\rho_n$ into the LL, we can also formulate MBAR as an optimization problem with parameters $c_k$:

$$LL(c_k) = -\sum_k^K N_k \log c_k - \sum_n^N \log \sum_k^K N_k c_k^{-1} q_k(x_n)$$


---
title: Solving MBAR Quickly and Precisely

- Nonlinear optimization : maximize log-likelihood (fast, robust, imprecise)
- Nonlinear Equations : find roots of gradient (fast, fragile, precise)
- Self Consistent Iteration : self-consistent equation (slow, robust, precise)

Caveats:

- Precision loss--objective function reduces 10^6 numbers into a single number
- Line search failures in BFGS implementations
- Speed
- Scaling

Conclusion: for precision+speed, need to combine BFGS and NR-type.

---
title: Solving MBAR: Self Consistent Iteration

$$f_i^{n+1} = - \log \sum_n^N \frac{\exp(-u_i(x_n))}{\sum_k^K N_k \exp(f_k - u_k(x_n))}$$

Notice that this expression can be written as a sequence of two logsumexp operations

$$d_n = \log \sum_k^K \exp(f_k - u_k(x_n) + \log N_k)$$
$$f_i^{n+1} = -\log \sum_n^N \exp(-u_i(x_n) - d_n)$$

The MBAR objective function and gradient can also be written using logsumexp!

---
title: Fast logsumexp using NumExpr

$$ \log \sum_n b_n \exp a_n = c + \log \sum_n b_n \exp(a_n - c) $$

NumExpr is a python library for optimizing large-scale algebraic operations.

<pre class="prettyprint" data-lang="python">
import numexpr
import numpy as np

def logsumexp(a, axis):
    a_max = np.amax(a, axis=axis, keepdims=True)
    return a_max + np.log(numexpr.evaluate("exp(a - a_max)").sum(axis))
    # It's actually slightly more complicated than this, but you get the idea

def self_consistent_iteration(u_kn, N_k, f_k):
    log_denominator_n = logsumexp(f_k - u_kn.T, b=N_k, axis=1)
    return -1. * logsumexp(-log_denominator_n - u_kn, axis=1)

</pre>

---
title: Fast logsumexp using NumExpr

Under ideal conditions, outperforms Numpy and matches hand-written C!

<article>
<iframe  data-src="file:///home/kyleb/src/kyleabeauchamp/MBARJournalClub/notebook/Benchmark.html"></iframe>
</article>

---
title: Conclusions and Future Work

- MBAR is an estimator for combining samples from multiple distributions
- Estimate arbitrary expectations via summation / quadrature.
- Finite support on $\{x_n\}$
- MLE
- Bayesian MBAR--instead of MLE, sample the posterior.
- Gaussian Process MBAR--plugging in alternative kernels that are supported for all $x$, rather than just $\{x_n\}$?
- If you want domain users to understand a model, you *must* describe it without measure theory

<footer class="source"> 
Csanyi, 2014.  
Habeck, 2012
</footer>
