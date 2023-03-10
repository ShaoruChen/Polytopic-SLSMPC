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
4. Using the virtual additive disturbances, we can easily derive an upper bound on the worst-case cost of the robust optimal control problem and minimize it.

### Effectiveness of SLS MPC
We compare SLS MPC with several baselines, including tube-based methods and methods that also optimize over LTV state feedback controllers, in solving the robust optimal control problem and compare their feasible domains with varying uncertainty parameters. 

<p float="left">
<img src="https://github.com/ShaoruChen/web-materials/blob/main/polytopic_SLS_MPC/coverage_comparison_eps_A.png" width="450" height="400">
<img src="https://github.com/ShaoruChen/web-materials/blob/main/polytopic_SLS_MPC/coverage_comparison_w.png" width="450" height="400">
</p>

In the above figures, coverage = (size of the feasible domain of each MPC method)/(size of the maximal robust control invariant set). The denominator is the theoretical upper bound on the feasible domain of any robust MPC method. $\epsilon_A$ denotes the level of model uncertainty and $\sigma_w$ denotes the magnitude of the additive disturbances $w_t$ considered. Note that SLS MPC always achieves more than 90% coverage even when the uncertainty parameters become large. More details about this example are included in the paper. 

## Robust MPC baselines
Different robust MPC methods are implemented in the mpc folder and are summarized below (naming follows from the SLS MPC paper).

### Tube-based methods

- (Tube-A)\
[Robust model predictive control using tubes](https://www.sciencedirect.com/science/article/abs/pii/S0005109803002838?casa_token=Af0HIAR3diAAAAAA:GBx9AXf6S43f7strRflfa-yPYxftN7A2oQKMz_tDXXn59TNMsvGPhLd7dCTFxC9PJ4MarINO8l0)\
Wilbur Langson, Ioannis Chryssochoos, S. V. Raković, and David Q. Mayne. \
Automatica 40, no. 1 (2004): 125-133.

- (Tube-B)\
[Robust MPC with recursive model update](https://www.sciencedirect.com/science/article/abs/pii/S0005109819300731?casa_token=wRxFJ9AGIHAAAAAA:D-Iw4Y9aObkUuWo1aWUeamVMnDwG354Y5A1kRjZv19hYFqqjGhsGNXkIkXK_Pjs7aq1yYM3iPjg)\
Matthias Lorenzen, Mark Cannon, and Frank Allgöwer.\
Automatica 103 (2019): 461-471.

- (Tube-C)\
[Linear robust adaptive model predictive control: Computational complexity and conservatism](https://ieeexplore.ieee.org/abstract/document/9028970) \
Johannes Köhler, Elisa Andina, Raffaele Soloperto, Matthias A. Müller, and Frank Allgöwer.\
In 2019 IEEE 58th Conference on Decision and Control (CDC), pp. 1383-1388. IEEE, 2019.

- (Tube-D)\
[Robust adaptive tube model predictive control](https://ieeexplore.ieee.org/abstract/document/8814456) \
Xiaonan Lu, and Mark Cannon. \
In 2019 American Control Conference (ACC), pp. 3695-3701. IEEE, 2019.

### Uncertainty over-approximation-based methods

- (Lumped-Disturbance-MPC)\
[A simple robust MPC for linear systems with parametric and additive uncertainty](https://ieeexplore.ieee.org/abstract/document/9482957) \
Monimoy Bujarbaruah, Ugo Rosolia, Yvonne R. Stürz, and Francesco Borrelli. \
In 2021 American Control Conference (ACC), pp. 2108-2113. IEEE, 2021.

- (Offline-Tightening-MPC)\
[Robust MPC for LPV systems via a novel optimization-based constraint tightening](https://www.sciencedirect.com/science/article/pii/S0005109822003156?casa_token=CiP-W89fO8UAAAAA:FyjgmtjKfeg5t5mxDe6LX2yyksPzWr7YjN142CzQV2TaD5g3Al4Ar2S4Gxk9we0A3JMxTw6_cHk) \
Monimoy Bujarbaruah, Ugo Rosolia, Yvonne R. Stürz, Xiaojing Zhang, and Francesco Borrelli. \
Automatica 143 (2022): 110459.

- (SLS-MPC) \
[Robust Model Predictive Control with Polytopic Model Uncertainty through System Level Synthesis](https://arxiv.org/abs/2203.11375)\
Shaoru Chen, Victor M. Preciado, Manfred Morari, Nikolai Matni\
Under review of Automatica

## Installation
Add the [mpc](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/tree/main/mpc) folder to MATLAB path and then you can run the [examples](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/tree/main/examples) in the paper. 

### Required toolboxes
[Yalmip](https://yalmip.github.io/) for formulating the control problems. [MOSEK](https://docs.mosek.com/9.3/toolbox/install-interface.html) is used as the default solver in the codes. 

[MPT3](https://www.mpt3.org/) for polyhedron operations. 

[MatlabProgressBar](https://www.mathworks.com/matlabcentral/fileexchange/57895-matlabprogressbar) for progress display (Not required if you remove the progress function in each for-loop, e.g. for i = progress(1:10) --> for i = 1:10).
