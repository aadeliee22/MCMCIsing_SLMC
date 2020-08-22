# I. Introduction

#### 1. Perform a Monte Carlo simulation on  [Ising model](https://en.wikipedia.org/wiki/Ising_model). 
> Here, I'll discuss about square lattice spin-1/2 simple classical Ising model.

​	For this Ising model, the energy configuration of state $\{s_i\}$ is given by Hamiltonian,  
$$
H = -J\sum_{<i,j>} s_is_j+ B\sum_i s_i
$$
where $<i,j>$ means the nearest neighbor in lattice site, and $s_i = \pm 1$. 
	This model is a simplest model of a magnet. $J$ indicates the interaction between nearest neighbors, it is ferromagnetic for $J>0$, and anti-ferromagnetic for $J<0$. ($B$ stands for external magnetic field, which I will set as $0$.) Ising model has a critical point (second-transition point), and I'm going to invest the thermal properties near the point.

 <figure> 
     <img src="C:\Users\User\Desktop\mygit\others\graph40.png" alt="" style="zoom:110%;" align="center"/>
     <figcaption><bold>Lattice size = 40</bold></figcaption>
 </figure>

There exist some exact solution for 1 or 2-dimension lattice, however, I'm going to examine this model of finite size using Monte Carlo simulation by C++. 

#### 2.  Self learning Monte Carlo method
​	For classical Ising model, there exist a global update (e.g. cluster update) that reduces auto-correlation time successfully. However, for other models, we can only apply local update, which is extremely slow near critical point. By SLMC, we can expect faster global update for any model given by Hamiltonian. 

# II. Background

We'll discuss about some thermodynamic/statistical backgrounds as well as Monte Carlo methods. Also, we'll briefly talk about some solutions for Ising model.

## 1. Thermodynamics & Statistical physics
### 1.1. Partition function
 [Partition function](https://en.wikipedia.org/wiki/Partition_function_%28statistical_mechanics%29) contains important information for a given system. 
$$
Z = \sum_{\{s\}}\exp(-\mathcal{H}(s)/kT), \; \mathcal{H} \text{ : Hamiltonian}
$$
 By using $Z$, we can drive some important physical quantities. 

- Probability of state: $P_s = \exp(-\mathcal{H}(s)/kT)/Z$

However, it is usually hard to get exact function when the system size increases, as we cannot sum all existing states. It is hard for one to drive thermodynamic quantities from ambiguous partition function.
### 1.2. Free energy
 of a system is given by $$\mathcal{F} = -kT\ln Z$$ and other thermodynamic quantities can be obtained by simply calculation. As $\mathcal{F} = U - TS$, 
$$
\begin{aligned}

d\mathcal{F} &= dU - TdS -SdT 
\\ &= (TdS-PdV+\mu dN) -TdS-SdT
\\ &= -SdT - PdV+\mu dN

\end{aligned}
$$
Which implies 

$$S = -\left(\frac{\partial \mathcal{F}}{\partial T}\right)_{V, N}, P = -\left(\frac{\partial \mathcal{F}}{\partial V}\right)_{T, N}, \mu = -\left(\frac{\partial \mathcal{F}}{\partial N}\right)_{T, V}$$ 

Moreover, using these, we can obtain order parameters. Here, $\beta = kT$.
$$
U = -\frac{1}{Z}\frac{\partial Z}{\partial\beta} = -\frac{\partial}{\partial\beta}\ln Z = -T^2 \frac{\partial (\mathcal{F}/T)}{\partial T} = \langle E \rangle
\\ \begin{aligned} C &= T \frac{\partial S}{\partial T} = -\beta \frac{\partial S}{\partial \beta}= \frac{\partial U}{\partial T} = k\beta^2 \frac{\partial^2}{\partial \beta^2}\ln Z \\&=k\beta^2\left(\frac{1}{Z}\frac{\partial^2}{\partial\beta^2}Z - \left[ \frac{1}{Z}\frac{\partial}{\partial \beta}Z\right]^2\right) =k\beta^2(\langle E^2\rangle - \langle E\rangle^2) \end{aligned}
\\M = \frac{\partial \mathcal{F}}{\partial B} = \left\langle \sum_i s_i \right\rangle, \;\; \chi = \frac{\partial\langle M\rangle}{\partial B} = \beta (\langle M^2\rangle - \langle M \rangle ^2)
$$

where $U$ is internal energy, $C$ is [specific heat capacity](https://en.wikipedia.org/wiki/Specific_heat_capacity), $M$ is [magnetization](https://en.wikipedia.org/wiki/Magnetization), and $\chi$ is [magnetic susceptibility](https://en.wikipedia.org/wiki/Magnetic_susceptibility).

## 2. Some solution of classical Ising model

### 2.1. Exact solution for 1-dimension Lattice



### 2.2. Mean Field theory for d-dimension Lattice



### 2.3. Exact solution for 2-dimension Lattice



## 3. Markov Chain Monte Carlo (MCMC)



### 3.1. Local update: Metropolis-Hastings algorithm



### 3.2. Global update: Cluster update



# III. Method

> c++: parameters, steps, throwing steps

## 1. Ising model Analysis

### 1.1. Finite size scaling: critical exponent



### 1.2. Error analysis: Jack knife



### 1.3. Autocorrelation time



## 2. Self Learning update method

### 2.1. Consideration of NNN on formation of cluster

### 2.2. Change acceptance ratio of cluster flipping





# IV. Result & Discussion

## 1. Classical Ising model

### 1.1. Comparison of Local & Global update

#### 1.1.1. R square

#### 1.1.2. Integrated Autocorrelation time

### 1.2. Comparison of Lattice size

#### 1.2.1. Thermodynamic quantities

#### 1.2.2. Binder cumulant

#### 1.2.3. Finite size scaling



## 2. Self Learning Monte Carlo method

> SLMC on plaquette-Ising model: E vs T, E vs E_eff

### 2.1. Consideration of NNN on formation of cluster

### 2.2. Change acceptance ratio of cluster flipping

# V. Conclusion

> shortcoming: no consideration of autocorrelation time of SLMC


# VI. Reference
1.  M. E. J. Newman and G. T. Barkema, _Monte Carlo methods in statistical physics_. Oxford: Clarendon Press, 2001.
2. L. E. Reichl, _A modern course in statistical physics_. Weinheim: Wiley-VCH, 2017.
3.  C. J. Adkins, _An introduction to thermal physics_. Cambridge: Cambridge University Press, 2004.
4. D. P. Landau and K. Binder, _A guide to Monte Carlo simulations in statistical physics_. Cambridge, United Kingdom: Cambridge University Press, 2015.
5. K. Rummukainen, in _Monte Carlo simulation methods_, 2019.
6. F. Wood, “Markov Chain Monte Carlo (MCMC),” in _C19 MACHINE LEARNING - UNSUPERVISED_, Jan-2015.
7. D. Foreman-Mackey, “Autocorrelation time estimation,” _Dan Foreman-Mackey_, 17-Oct-2017. [Online]. Available: https://dfm.io/posts/autocorr/. [Accessed: 20-Aug-2020].
8. J. Liu, Y. Qi, Z. Y. Meng, and L. Fu, “Self-learning Monte Carlo method,” _Physical Review B_, vol. 95, no. 4, 2017.
