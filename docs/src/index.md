# Quasi-Monte Carlo Generators

Julia version of Dirk Nuyen's [Magic Point Shop](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/). 

QMCGenerators implements quasi-random (low discrepancy) sequence generators. Quasi-random points are carefully coordinated in a dependent manner to more evenly fill $[0,1]^s$ compared to pseudo-random (independent) points. 

Monte Carlo methods approximate the true mean

$$\mu = \mathbb{E}[f(U)] = \int_{[0,1]^s} f(u) \mathrm{d} u$$

for $f: [0,1]^s \to \mathbb{R}$ and $U \sim \mathcal{U}[0,1]^s$ by the sample mean

$$\hat{\mu} = \sum_{i=1}^n f(U_i)$$

where $U_1,U_2,\dots \sim \mathcal{U}[0,1]^s$. If $(U_i)_{i \geq 1}$ are chosen to be pseudo-random, then the sample average is a **Simple Monte Carlo** approximation with error $\mathcal{O}(1/\sqrt{n})$. If instead we choose $(U_i)_{i \geq 1}$ to be quasi-random then the sample average is a **Quasi-Monte Carlo** approximation with error approaching $\mathcal{O}(1/n)$. 

This package implements two flavors of quasi-random sequences: **Lattice rules** and **digital nets** in base 2. Independent randomizations may be applied to base sequences with random shifts for Lattices and random digital shifts for digital nets.  

## References

1. Kuo, F. Y., & Nuyens, D. (2016). [Application of quasi-Monte Carlo methods to elliptic PDEs with random diffusion coefficients: a survey of analysis and implementation](https://link.springer.com/article/10.1007/s10208-016-9329-5). Foundations of Computational Mathematics, 16, 1631-1696.
