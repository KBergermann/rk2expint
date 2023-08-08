% 
% Main function reproducing the numerical examples from [1]
% 1. Obtain required external software and network data
% 2. Choose poles and time step size
% 3. run_problem
% 
% Syntax:
% run_problem(pde, exptAb_routine, exp_rk_int, xi, xi_unique, dt, T_array,
% linear_system_solver, n, network, compute_error_to_ode15s_reference_solution,
% plot_solution)
% 
% Implemented options for inputs to run_problem:
% 
% pde:                   'AC_2D', 'GM_2D', 'AC_network', 'GM_network'
% exptAb_routine:        'phipm', 'kiops', 'rk2expint'
% exp_rk_int:            'SW2', 'ETD3RK', 'Krogstad4'
% xi, xi_unique:         1D-arrays of poles
% dt:                    time step size
% T_array:               1D-array with time points
% linear_system_solver:  'lu_Matlab', 'lu_pardiso', 'AGMG'
% network:               'minnesota', 'usroads_subset', 'ak', 'luxembourg_osm',
%                        'ny', 'roadNet_PA', 'brightkite'
% n:                     number of grid points per dimension (for AC_2D and GM_2D)
% compute_error_to_ode15s_reference_solution, plot_solution:
%                        boolean and self-explaining. For ode15s, the default
%                        tolerance is an absolute and relative tolerance of
%                        1e-06. This can be adapted in run_problem.
% 
% Reference:
% [1] K. Bergermann and M. Stoll, Adaptive rational Krylov methods for exponential Runge--Kutta integrators, arxiv preprint arXiv:2303.09482, (2023).
% [2] R.-U. Börner, O. G. Ernst, and S. Güttel, Three-dimensional transient electromagnetic modelling using rational Krylov methods, Geophysical Journal International, 202 (2015), pp. 2025–2043.
% [3] A. Carpenter, A. Ruttan, and R. Varga, Extended numerical computations on the “1/9” conjecture in rational approximation theory, in Rational Approximation and Interpolation, Springer, 1984, pp. 383–411.
% [4] M. Berljafa and S. Güttel, The RKFIT algorithm for nonlinear rational approximation, SIAM J. Sci. Comput., 39 (2017), pp. A2049–A2071.
% 

addpath('time_integrators')
addpath('system_solve_fhs')
addpath('subroutines')
addpath('kiops')
addpath('phipm')
addpath('rktoolbox')
addpath('AGMG_3.3.5-aca/Matlab')

%% Choose poles
%%% 1 cyclically repeated real pole, cf. [2]
% load('pole_files/one_repeated_real_pole.mat')
%%% 2 cyclically repeated real poles, cf. [2]
% load('pole_files/two_repeated_real_poles.mat')
%%% 4 cyclically repeated real poles, cf. [2]
% load('pole_files/four_repeated_real_poles.mat')
%%% Rational best approximation poles, cf. [3]
% load('pole_files/rational_best_approximation_poles.mat');
%%% Poles optimized by RKFIT [4] (compute_rkfit_poles.m in this directory)
load('pole_files/expint_poles.mat');

%% Set time array and run problem
T = 1; % final time
dt = 0.5; % time step size
T_array = 0:dt:T; % time grid

u = run_problem('AC_2D', 'rk2expint', 'SW2', xi, xi_unique, dt, T_array, 'AGMG', 100, 'minnesota', true, true);
