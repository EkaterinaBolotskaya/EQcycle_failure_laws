# EQcycle_failure_laws

## Background
1D spring-slider model is often used to simulate earthquake cycles. Failure law prescribed between the block and the rough surface plays an important role in the system behavior. This repository contains MATLAB scripts that numerically solve the equation of motion for 1D spring-slider system: 

$$
\ddot{D} = -\frac{1}{M}(K(D - V_{0} t) +\mu\sigma_n), \qquad \qquad \qquad \qquad (1) 
$$

where $M$ is the mass of the system, $K$ is the spring stiffness, $D$ is the slip, $V$ is the slip rate, $V_0$ is the load point velocity, $\sigma_n$ is the normal stress, and $\mu$ is the friction coefficient,

with several widely used and custom failure laws:

<ins>__Rate-and-state friction (RS) with aging law (RS):__</ins>

$$
\begin{cases}
      \mu = \mu_0 + a \ln \frac{V}{V_0} + b \ln \frac{V_0 \theta}{L}\\
      \dot \theta = 1 - \frac{V \theta}{L}\\
    \end{cases}
    ,  \qquad \qquad \qquad \qquad (2)
$$

where $\mu$ is the friction coefficient, $\mu_0$ is the reference friction coefficient, $V$ is the block slip rate, $\theta$ is the state variable, representing the sliding history, $V_0$ is the reference slip rate, $L$ is the characteristic length scale, and $a$ and $b$ are rate-and-state parameters.

<ins>__Slip-weakening friction (SW):__</ins>

$$
\mu =
    \begin{cases}
      \mu_s - (\mu_s-\mu_d)\frac{D}{D_c}               & D \le D_c\\
      \mu_d                                            & D > D_c\\
    \end{cases} 
    , \qquad \qquad \qquad \qquad (3)
$$

where $\mu$ is the friction coefficient, $\mu_s$ and $\mu_d$ are the static and dynamic friction coefficients, $D$ is the block slip, and $D_c$ is the slip-weakening distance.

<ins>__Exponential cohesive zone (ECZ):__</ins>

$$
\mu = \mu_d + (\mu_s-\mu_d)\frac{D+D_1}{D_2} e^{1-\frac{D+D_1}{D_2}} ,  \qquad \qquad \qquad \qquad (4)
$$

where $\mu$ is the friction coefficient, $\mu_s$ and $\mu_d$ are the static and dynamic friction coefficients, $D$ is the block slip, and $D_1$ and $D_2$ are the two slip distances which we call "slip-shift" and "slip-stretch" respectively.

<ins>__Parabolic cohesive zone (PCZ):__</ins>

$$
\mu =
    \begin{cases}
      \mu_s - (\mu_s-\mu_d) \left( \frac{D-D_1}{D_2} \right) ^2              & D \le D_1+D_2\\
      \mu_d                                            & D > D_1+D_2\\
    \end{cases} 
    ,  \qquad \qquad \qquad \qquad (5)
$$

where $\mu$ is the friction coefficient, $\mu_s$ and $\mu_d$ are the static and dynamic friction coefficients, $D$ is the block slip, and $D_1$ and $D_2$ are the the slip-strengthening and slip-weakening distances respectively.

<ins>__Double slip-weakening with initial strengthening (DSWIS) model:__</ins>

$$\mu =
    \begin{cases}
      \mu_i - (\mu_i-\mu_s)\frac{D}{D_s}               & D \le D_s\\
      \mu_s - (\mu_s-\mu_t)\frac{D-D_s}{D_{w1}}        & D_s < D \le D_s+D_{w1}\\
      \mu_t - (\mu_t-\mu_d)\frac{D-D_s-D_{w1}}{D_{w2}} & D_s+D_{w1} < D \le D_s+D_{w1}+D_{w2}\\
      \mu_d                                            & D > D_s+D_{w1}+D_{w2}\\
    \end{cases}
    , \qquad \qquad \qquad \qquad (6)
$$

where $\mu_i$, $\mu_s$, $\mu_t$, and $\mu_d$ are the initial, static, transitional, and dynamic friction coefficients, $D_{w1}$ and $D_{w2}$ are the slip-weakening distances, and $D_s$ is the slip-strengthening distance. The last segment of the failure law is horizontal with $\mu=\mu_d$.

<ins>__Reference model:__</ins>

This is the case of instantaneous transition from $\tau_s$ to $\tau_d$ over a negligibly small slip-weakening distance $(D_c=0)$.

$$
\mu =
    \begin{cases}
      \mu_s                & D = 0\\
      \mu_d                & D > 0\\
    \end{cases} 
    , \qquad \qquad \qquad \qquad (7)
$$

where $\mu$ is the friction coefficient, $\mu_s$ and $\mu_d$ are the static and dynamic friction coefficients, $D$ is the block slip.

The input parameters used in the simulations are chosen in such a way that the area under the weakening part of the shear stress vs. slip curve (energy curve) $E_w$ is approximately the same for all failure laws. The area under the strengthening part of the curve $E_s$ is approximately the same for most of the failure laws $(E_s = 0$ for SW by definition). Both $E_w$ and $E_s$ equal zero for the Reference model.

Different failure laws that produce very similar coseismic ruptures can have substantially different nucleation phases and pre-nucleation slip rate and slip evolution. This also results in differences in recurrence intervals, maximum coseismic slip rate, and cumulative slip amount per earthqake cycle.

## Repository contents
- MATLAB scripts:
  - *spring_slider_DSWIS*
  - *spring_slider_EXP* 
  - *spring_slider_PAR* 
  - *spring_slider_RS_aging_law*
  - *spring_slider_SW*
  - *Process_and_plot_data_for_several_FL*
- Data_FL repository
- README.md
- LICENSE

### *spring_slider_DSWIS*
Numerically solves non-dimensional equation (1) with a DSWIS failure law (6), plots the output for several earthquake cycles, saves the data in Data_FL folder (necessary to then plot comparative plots). Two data sets are produced by this script: Steep to gentle (DSWIS 1) and Gentle to steep (DSWIS 2).
### *spring_slider_EXP*
Numerically solves non-dimensional equation (1) with an ECZ failure law (4), plots the output for several earthquake cycles, saves the data in Data_FL folder (necessary to then plot comparative plots).
### *spring_slider_PAR*
Numerically solves non-dimensional equation (1) with a PCZ failure law (5), plots the output for several earthquake cycles, saves the data in Data_FL folder (necessary to then plot comparative plots).
### *spring_slider_RS_aging_law*
Numerically solves non-dimensional equation (1) with RS friction with aging state variable evolution law (2), plots the output for several earthquake cycles, saves the data in Data_FL folder (necessary to then plot comparative plots).
### *spring_slider_SW*
Numerically solves non-dimensional equation (1) with SW friction (3), plots the output for several earthquake cycles, saves the data in Data_FL folder (necessary to then plot comparative plots).
### *Process_and_plot_data_for_several_FL*
Processes the data from Data_FL repository for different failure laws. Plots slip and slip rate for the full cycle, energy curves, phase diagrams, spectra, slip and slip rate by phase etc. 

## Reference
The results obtained in these scripts are part of my Thesis: 
Bolotskaya, E., 2023. Effects of fault failure parameterization and bulk rheology on earthquake rupture (Doctoral dissertation, Massachusetts Institute of Technology).

and a publication in prep.
Please contact me if you use EQcycle_failure_laws for your research and would like a citation.

Release on Zenodo:

[![DOI](https://zenodo.org/badge/530423563.svg)](https://zenodo.org/badge/latestdoi/530423563)
