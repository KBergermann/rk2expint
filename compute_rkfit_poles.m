% Computes optimal poles for the function f(x)=e^(-x) using the RKFIT algorithm [1].
% In [1, Section 6], two sets of poles are used. The difference of their
% computation by this script are two lines:
% 
% - for the (35,30) poles on the interval [0, 1e06]:
%       ee = [0 , logspace(-6, 6, N-1)];
%       m = 30; k = 5;
% - for the (70,60) poles on the interval [0, 1e08]:
%       ee = [0 , logspace(-6, 8, N-1)];
%       m = 60; k = 10;
% 
% Code taken from http://guettel.com/rktoolbox/examples/html/example_expint.html
% 
% xi and xi_unique are saved in 'pole_files/expint_poles.mat' for further use in main.m
% 
% Input:
%       N:   number of sample points from the spectrum/approximation interval
%       ee:  vector of logspaced sample points in the spectrum/approximation interval
%       t:   vector of logspaced time steps
%       m:   denominator degree
%       m+k: numerator degree
% 
% Output:
%       xi:  Array of optimized poles (of size 1-by-m, possibly complex-valued)
% 
% Reference:
% [1] M. Berljafa and S. Güttel, The RKFIT algorithm for nonlinear rational approximation, SIAM J. Sci. Comput., 39 (2017), pp. A2049–A2071.
% 

addpath('rktoolbox')

N  = 500;
ee = [0 , logspace(-6, 6, N-1)];
A  = spdiags(ee(:), 0, N, N);
b  = ones(N, 1);
t  = logspace(-2, 0, 41);
for j = 1:length(t)
  F{j} = spdiags(exp(-t(j)*ee(:)), 0, N, N);
end

m = 30; k = 5;               % type (m+k, m)
xi = Inf(1,m);               % initial poles at infinity
param.k = k;                 % subdiagonal approximant  
param.maxit = 10;            % at most 10 RKFIT iterations
param.tol   = 0;             % exactly 10 iterations
param.real  = 1;             % data is real-valued
param.stable = 1;
[xi, ratfun, misfit, out] = rkfit(F, A, b, xi, param);
xi_unique = xi;
save('pole_files/expint_poles.mat', 'xi', 'xi_unique')
