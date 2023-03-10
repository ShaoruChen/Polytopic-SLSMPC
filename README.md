# Polytopic-SLSMPC

This repository contains codes to implement the robust model predictive control (MPC) method, SLS MPC, proposed in 

[Robust Model Predictive Control with Polytopic Model Uncertainty through System Level Synthesis](https://arxiv.org/abs/2203.11375)\
Shaoru Chen, Victor M. Preciado, Manfred Morari, Nikolai Matni\
Under review of Automatica

and the robust MPC baselines therein. 

## Robust MPC
In robust MPC of an uncertain linear time-invariant system, we consider solving a finite-horizon robust optimal control problem with the current state $x(k)$ as follows

<img src="https://github.com/ShaoruChen/web-materials/blob/main/polytopic_SLS_MPC/mpc_formulation.png" width="450" height="200">

where the model uncertainty parameters $(\Delta_A, \Delta_B)$ lies in a polytopic set which means they belong to a convex hull of a finite number of vertices and the additive disturbances $w_t$ belongs to a polytope $\mathcal{W}$. We want to find a robust feedback controller $\pi_t(\cdot)$ such that the polytopic state and input constraints $\mathcal{X}$ and $\mathcal{U}$ are robustly satisfied. 
