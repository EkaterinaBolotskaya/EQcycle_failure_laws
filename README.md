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

where $\mu$ is the friction coefficient, $\mu_s$ and $\mu_d$ are the static and dynamic friction coefficients, $D$ is the block slip, and $D_1$ and $D_2$ are the two slip distances which we call "slip-shift" and "slip-stretch" respectively.

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

The input parameters used in the simulations are chosen in such a way that the area under the weakening part of the shear stress vs. slip curve (energy curve) $E_w$ is approximately the same for all failure laws. The area under the strengthening part of the curve $E_s$ is approximately the same for most of the failure laws $(E_s = 0$ for SW by definition).

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
Analytically solves non-dimensional equation (1) with a generic linear friction segment, shows the three solution regimes: $K_k^f < K$ - harmonic oscillations, $K_k^f=K$ - cubic growth solution, and $K_k^f>K$ - exponential growth solution.
### *spring_slider_EXP*
Analytically solves the 1D dynamic spring slider equation with a double slip weakening with initial strengthening (DSWIS) failure law (2). Shows full analytic solutions for each segment with initial conditions from the previous segment.
### *spring_slider_PAR*
Mostly analytically (the equation to find the duration of different phases does not have analytical solutions, thus we solve for them numerically) solves equation (1) with DSWIS (2) and produces plots for a single set of failure law parameters: energy curves, phase diagrams, slip rate and slip vs. time for several earthquake cycles, slip and slip rate plot for different phases separately, spectra.
### *spring_slider_RS_aging_law*
Mostly analytically (same as above) solves equation (1) with DSWIS (2) and produces plots for several sets of failure law parameters (with the same axis scales) for comparison. Different sets of parameters (3 to 8 failure laws) are given as examples.
### *spring_slider_SW*
Estimates the lower bound on frequency of the oscillatoric solution of block slip with poly-linear friction for a range of fault lengths and slip-weakening distances, assuming the block goes through a single oscillation during the weakening process.
### *Process_and_plot_data_for_several_FL*
Analytically solves non-dimensional equation (1) with a generic linear friction segment, shows the three solution regimes: $K_k^f < K$ - harmonic oscillations, $K_k^f=K$ - cubic growth solution, and $K_k^f>K$ - exponential growth solution.

## Reference
Please refer the following article if you use EQcycle_polylinear for your research:

E. Bolotskaya and B.H. Hager; A 1D Spring‐Slider Model with a Simple Poly‐Linear Failure Law Produces Rich Variations in Slip Behavior. Bull. Seismol. Soc. Am. 2022; doi: https://doi.org/10.1785/0120220052

Release on Zenodo:

[![DOI](https://zenodo.org/badge/434003826.svg)](https://zenodo.org/badge/latestdoi/434003826)
