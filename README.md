# MATLAB Numerical Analysis Library

A high-performance, modular MATLAB toolkit designed for rigorous numerical computing, differential equation modeling, and scientific simulation. Built with a strong emphasis on vectorized performance, sparse matrix operations, and robust input sanitization, this library serves as a versatile backend for engineering analysis and mathematical modeling.

---

## Key Features & Modules

### 1. Nonlinear Equations & Optimization
* **Modified Newton's Method (`modified_newton.m`)**: Root-finding algorithm supporting custom multiplicity parameters for fast convergence on complex roots.
* **Levenberg-Marquardt Algorithm (`levenberg_marquardt.m`)**: Robust nonlinear least-squares solver featuring adaptive damping ($\lambda$) updates to transition smoothly between Gauss-Newton and gradient descent.
* **Fixed-Point Iteration (`fixed_point_iteration.m`)**: Vectorized fixed-point solver with infinity-norm stopping criteria.

### 2. Interpolation & Data Approximation
* **Lagrange Interpolation (`lagrange_interpolation.m`)**: Supports both analytical function evaluation over uniform nodes and discrete experimental data points, complete with automatic pointwise error analysis.
* **Lagrange Basis Polynomials (`lagrange_basis_polynomials.m`)**: Computes characteristic basis functions using anonymous functions and vectorized evaluations without requiring the Symbolic Math Toolbox.
* **Least-Squares Approximation (`least_squares_approximation.m`)**: Polynomial regression and curve fitting for noisy or discrete datasets.
* **Sum of Squared Errors (`sum_squared_error.m`)**: Utility function for residual evaluation.

### 3. Ordinary Differential Equations (ODEs) & Stability
* **Generic Runge-Kutta (`runge_kutta.m`)**: Solves arbitrary ODE systems using fully customizable Butcher tableaux.
* **Linear $\theta$-Method (`linear_theta_method.m`)**: Unifies Forward Euler ($\theta=0$), Backward Euler ($\theta=1$), and Crank-Nicolson ($\theta=0.5$) with implicit matrix factorization.
* **Stability Analysis (`check_stability.m`)**: Evaluates amplification factors $R(z)$ for various time-stepping schemes and plots theoretical stability regions in the complex plane.

### 4. Boundary Value Problems (BVPs) & Partial Differential Equations (PDEs)
* **1D BVP Solver (`solve_bvp_1d.m`)**: Solves general second-order linear BVPs with variable diffusion, advection, and reaction terms under mixed Dirichlet and Neumann boundary conditions (including ghost-node approximations).
* **1D Advection-Diffusion-Reaction (`variable_advection_diffusion_reaction_1d.m`)**: Implements adaptive Upwind differencing schemes to handle variable wind directions ($v(x)$) without spurious oscillations.
* **2D Dirichlet Poisson/Helmholtz Solver (`poisson_2d_dirichlet.m`)**: High-performance 2D solver utilizing Kronecker products (`kron`) and sparse matrices (`spdiags`) to handle spatially variable reaction terms $\sigma(x,y)$ on structured grids efficiently.

---

## Quick Start Example

Here is a quick example showing how to solve a 1D advection-diffusion-reaction boundary value problem using the library:

```matlab
% Define problem parameters
mu = 0.05;                        % Diffusion coefficient
eta_fun = @(x) 2.0 + 0*x;         % Constant advection velocity
sigma_fun = @(x) 1.0 + x.^2;      % Spatially variable reaction term
f_fun = @(x) sin(pi*x);           % Forcing source term
L = 1.0;                          % Domain length [0, 1]
h = 0.01;                         % Spatial step size
alpha = 0.0;                      % Left Dirichlet condition u(0)
beta = 1.0;                       % Right Dirichlet condition u(L)

% Solve the system
[U_com, x] = variable_advection_diffusion_reaction_1d(mu, eta_fun, sigma_fun, f_fun, L, h, alpha, beta);

% Plot results
figure('Color', 'w');
plot(x, U_com, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('u(x)');
title('1D Advection-Diffusion-Reaction Solution');
