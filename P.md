# Introduction
#### Perform a Monte Carlo simulation on classical, 2-dimensional, spin-$1/2$ [Ising model](https://en.wikipedia.org/wiki/Ising_model). 
For classical square lattice Ising model, the energy configuration of $\{s_i\}$ is given by Hamiltonian, which is  $$H = -J\sum_{<i,j>} s_is_j+ B\sum_i s_i$$ where $<i,j>$ means the nearest neighbor in lattice site, and $s_i = \pm 1$. 
This model is a simplest model of a magnet. $J$ indicates the interaction between nearest neighbors, it is ferromagnetic for $J>0$, and antiferromagnetic for $J<0$. ($B$ stands for external magnetic field, which I will set as $0$.) Ising model has a critical point (second-transition point), and I'm going to invest the thermal properties near the point.
There exist some exact solution for 1D or 2D lattice, however, I'm going to examine this model of finite size using Monte Carlo simulation by C++. 
#### Self learning Monte Carlo method
For classical Ising model, there exist a global update (e.g. cluster update) that reduces auto-correlation time successfully. However, for other models, we can only apply local update, which is extremely slow near critical point. By SLMC, we can expect faster global update for any model given by Hamiltonian. 

# Background

### 1. Ising model



# Method

# Result & Discussion

# Reference
1.  M. E. J. Newman and G. T. Barkema, _Monte Carlo methods in statistical physics_. Oxford: Clarendon Press, 2001.
2. L. E. Reichl, _A modern course in statistical physics_. Weinheim: Wiley-VCH, 2017.
3.  C. J. Adkins, _An introduction to thermal physics_. Cambridge: Cambridge University Press, 2004.
4. D. P. Landau and K. Binder, _A guide to Monte Carlo simulations in statistical physics_. Cambridge, United Kingdom: Cambridge University Press, 2015.
5. K. Rummukainen, in _Monte Carlo simulation methods_, 2019.
6. F. Wood, “Markov Chain Monte Carlo (MCMC),” in _C19 MACHINE LEARNING - UNSUPERVISED_, Jan-2015.
7. D. Foreman-Mackey, “Autocorrelation time estimation,” _Dan Foreman-Mackey_, 17-Oct-2017. [Online]. Available: https://dfm.io/posts/autocorr/. [Accessed: 20-Aug-2020].
8. J. Liu, Y. Qi, Z. Y. Meng, and L. Fu, “Self-learning Monte Carlo method,” _Physical Review B_, vol. 95, no. 4, 2017.

> Written with [StackEdit](https://stackedit.io/).