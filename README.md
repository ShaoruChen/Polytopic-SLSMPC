# Polytopic-SLSMPC

This repository contains codes to implement the robust model predictive control (MPC) method, SLS MPC, proposed in 

[Robust Model Predictive Control with Polytopic Model Uncertainty through System Level Synthesis](https://arxiv.org/abs/2203.11375)\
Shaoru Chen, Victor M. Preciado, Manfred Morari, Nikolai Matni\
Under review of Automatica

and the robust MPC baselines therein. 

## Overview of SLS MPC

### Problem formulation
In robust MPC of an uncertain linear time-invariant system, we consider solving a finite-horizon robust optimal control problem with the current state $x(k)$ as follows

<img src="https://github.com/ShaoruChen/web-materials/blob/main/polytopic_SLS_MPC/mpc_formulation.png" width="450" height="200">

where the model uncertainty parameters $(\Delta_A, \Delta_B)$ lies in a polytopic set which means they belong to a convex hull of a finite number of vertices and the additive disturbances $w_t$ belongs to a polytope $\mathcal{W}$. We want to find a robust feedback controller $\pi_t(\cdot)$ such that the polytopic state and input constraints $\mathcal{X}$ and $\mathcal{U}$ are robustly satisfied. 

There is a naturally a tension in the parameterization of the feedback policy $\pi_t$ and the constraint tightening to guarantee robust constraint satisfaction of the above robust optimal control problem. A complex policy parameterization of $\pi_t$ makes it hard to tighten the constraints, while a simple policy parameterization is inherently conservative. We need to balance these two aspects in a good way. 

### Features of SLS MPC
Our proposed method, SLS MPC, has the following features:
1. It searches linear time-varying state feedback controller $\pi_t(x_{0:t}) = K^{t,t}x_0 + \cdots + K^{t,0}x_t$ where the feedback gains $K^{t,t-i}$ are optimized online. 
2. Instead of searching the controller parameter $K^{t,t-i}$ directly which often leads to nonconvex optimization, we use [System Level Synthesis](https://arxiv.org/abs/1904.01634) (SLS) to reparameterize the controller in the system response space which reveals useful structures about the robust optimal control problem. 
3. SLS MPC introduces a virtual additive disturabnce signal to over-approximate the uncertainty effects and simplifies constraint tightening. The reach set of the virtual additive disturbance signal is optimized online jointly with the feedback controller parameters through a set of linear constraints. 

### Effectiveness of SLS MPC
We compare SLS MPC with several baselines, including tube-based methods and methods that also optimize over LTV state feedback controllers, in solving the robust optimal control problem and compare their feasible domains with varying uncertainty parameters. 

<p float="left">
<img src="https://github.com/ShaoruChen/web-materials/blob/main/polytopic_SLS_MPC/coverage_comparison_eps_A.png" width="450" height="400">
<img src="https://github.com/ShaoruChen/web-materials/blob/main/polytopic_SLS_MPC/coverage_comparison_w.png" width="450" height="400">
</p>

In the above figures, coverage = (size of the feasible domain of each MPC method)/(size of the maximal robust control invariant set). The denominator is the theoretical upper bound on the feasible domain of any robust MPC method. $\epsilon_A$ denotes the level of model uncertainty and $\sigma_w$ denotes the magnitude of the additive disturbances $w_t$ considered. Note that SLS MPC always achieves more than 90% coverage even when the uncertainty parameters become large. More details about this example are included in the paper. 

