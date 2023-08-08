function [w, m_ret, stats] = rk2expint(tau_out, A, u, tol, m_init, mmin, mmax, task1, xi, linear_system_solver, factorization, dt_2_flag)
% Rational Krylov Runge--Kutta exponential integrator (RK)^2EXPINT routine
% as introduced in [1] to adaptively compute linear combinations of
% \varphi-functions by means of approximating the action of the matrix
% exponential on a certain vector [1,2,3,4]. This code is an adaptation of
% the KIOPS routine (https://gitlab.com/stephane.gaudreault/kiops) [2],
% which, in turn, is based on the PHIPM and EXPMVP codes
% (http://www1.maths.leeds.ac.uk/~jitse/software.html) [3,4].
% 
% In essence, this routine replaces the incomplete orthogonalization
% procedure (IOP) of KIOPS by a rational Krylov procedure for which the
% rktoolbox implementation (http://guettel.com/rktoolbox/) [5] is used.
% This requires taking care of additional linear system solves as well as a
% novel a-posteriori error estimate for the rational Krylov case, 
% cf. [1, Subsection 4.3].
% 
% rk2expint is designed to be called by the same syntax as kiops, cf. the
% files in ../time_integrators
% 
% Input: 
%       tau_out:    time step size (h_i in [1, Theorem 3.2])
%       A:          discrete linear differential operator (sparse matrix)
%       u:          (block) vector (C in [1, Theorem 3.2])
%       tol:        desired accuracy of the rat. Krylov approximation to exp(A)c
%       m_init:     initial guess for the rat. Krylov subspace size
%       mmin:       minimum rat. Krylov subspace size
%       mmax:       maximum rat. Krylov subspace size
%       task1:      boolean, if true, allows for the approximation of a
%                   single \varphi-function at different points in time
%       xi:         1D-array of poles
%       linear_system_solver:
%                   string: 'lu_Matlab', 'lu_pardiso', or 'AGMG'
%       factorization:
%                   struct containing precomputed matrix decompositions of A
%       dt_2_flag:  boolean, if true, multiple matrix decompositions need
%                   to be accessed by direct linear system solvers
% 
% Output:
%       w:      desired linear combination of \varphi-functions
%       m_ret:  denominator degree used in rational Krylov subspace
%       stats:  some statistics on the procedure
% 
% References:
% [1] K. Bergermann and M. Stoll, Adaptive rational Krylov methods for exponential Runge--Kutta integrators, arxiv preprint arXiv:2303.09482, (2023).
% [2] Gaudreault, S., Rainwater, G. and Tokman, M., 2018. KIOPS: A fast adaptive Krylov subspace solver for exponential integrators. Journal of Computational Physics.
% [3] Niesen, J. and Wright, W.M., 2011. A Krylov subspace method for option pricing. SSRN 1799124
% [4] Niesen, J. and Wright, W.M., 2012. Algorithm 919: A Krylov subspace algorithm for evaluating the \varphi-functions appearing in exponential integrators. ACM Transactions on Mathematical Software (TOMS), 38(3), p.22
% [5] Berljafa, M., Elsworth, S., & Güttel, S. (2014). A rational Krylov toolbox for MATLAB.
% 


% n is the size of the original problem
% p is the highest indice of the phi functions
[n, ppo] = size(u);
p = ppo - 1;

if p == 0
   p = 1;
   % Add extra column of zeros
   u = [u, zeros(size(u))];
end

% Check inputs
if nargin < 12
   dt_2_flag = false;
   if nargin < 11
      factorization = struct;
      if nargin < 10
         linear_system_solver = 'lu_Matlab';
         if nargin < 9
            load('pole_files/expint_poles.mat');
            if nargin < 8
               task1 = true;
               if nargin < 7
                  mmax = 30;
                  if nargin < 6
                     mmin = 5;
                     if nargin < 5
                        m_init = mmin;
                        if nargin < 4
                           tol = 1.0e-7;
                           if nargin < 3
                              error('Not enough input arguments.');
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end
end

% We only allow m to vary between mmin and mmax
m = max(mmin, min(m_init, mmax));

% Preallocate matrix
V = zeros(n + p, mmax + 1);
H = zeros(mmax + 1, mmax + 1);
K = zeros(mmax + 1, mmax + 1);

step    = 0;
krystep = 0;
ireject = 0;
reject  = 0;
exps    = 0;
sgn     = -sign(tau_out(end));
tau_now = 0;
tau_end = abs(tau_out(end));
happy   = false;
j       = 0;

numSteps = size(tau_out, 2);

% Initial condition
w     = zeros(n, numSteps);
w_aug = zeros(p, 1);
w(:, 1) = u(:, 1);

% Normalization factors
normU = norm(u(:, 2:end),1);
if ppo > 1 && normU > 0
   ex = ceil(log2(normU));
   nu = 2^(-ex);
   mu = 2^(ex);
else
   nu = 1;
   mu = 1;
end

% Flip the rest of the u matrix
u_flip = nu*fliplr(u(:, 2:end));

% Compute and initial starting approximation for the step size
tau = tau_end;

% Setting the safety factors and tolerance requirements
if tau_end > 1
   gamma = 0.2;
   gamma_mmax = 0.1;
else
   gamma = 0.9;
   gamma_mmax = 0.6;
end
delta = 1.4;

% Used in the adaptive selection
oldm = NaN; oldtau = NaN; omega = NaN;
orderold = true; kestold = true;

l=1;
n_infty_poles = 0;

%% Extended matrix
Jp = spdiags([zeros(p,1), ones(p,1)], 0:1,p,p);
AA = -[A, u_flip; sparse(p,n), Jp];

%% rat_krylov multiply and solve
AB = struct;
AB.isreal = isreal(AA);
AB.A = AA;
AB.multiply = @(rho, eta, x) rho*(AA*x) - eta*x;

switch linear_system_solver
    case 'lu_Matlab'
        AB.solve = @(nu, mu, b) solve_lu_Matlab(Jp, u_flip, nu, mu, b, tau_out, factorization, dt_2_flag);
    case 'lu_pardiso'
        AB.solve = @(nu, mu, b) solve_lu_pardiso(A, Jp, u_flip, nu, mu, b, tau_out, factorization, dt_2_flag);
    case 'AGMG'
        AB.solve = @(nu, mu, b) solve_agmg_mult(A, u_flip, Jp, nu, mu, b);
end

while tau_now < tau_end

   if j==0
      % Update the last part of w
      for k=1:p-1
         i = p - k;
         w_aug(k) = (tau_now^i)/factorial(i) * mu;
      end
      w_aug(p) = mu;

      % Normalize initial vector (this norm is nonzero)
      beta = sqrt( w(:,l)' * w(:,l) + w_aug' * w_aug );

      % The first Krylov basis vector
      V(1:n, 1)     = (1/beta) * w(:,l);
      V(n+1:n+p, 1) = (1/beta) * w_aug;

      % Build rational Krylov subspace
      [V, K, H] = rat_krylov(AB, V(:,1), [xi(1:m), Inf]);
      
      % Keep track of the number of infinity poles, i.e., the difference
      % between numerator and denominator polynomial degree
      n_infty_poles = n_infty_poles + 1;
   else
      param.j = j;
      % Extend rational Krylov subspace
      [V, K, H] = rat_krylov(AB, V, K(1:j+1,1:j), H(1:j+1,1:j), [xi(oldm+1:m), Inf], param);
      
      % Keep track of the number of infinity poles, i.e., the difference
      % between numerator and denominator polynomial degree
      n_infty_poles = n_infty_poles + 1;
   end
   j = m + n_infty_poles;
   krystep = j;
   
   % Compute H_m K_m^(-1) for later use
   HK = zeros(j+1,j+1);
   HK(1:j,1:j) = H(1:j,1:j)*(K(1:j,1:j)\eye(j,j));

   % To obtain the phi_1 function which is needed for error estimate
   HK(1, j + 1) = 1;

   % Save h_j+1,j and remove it temporarily to compute the exponential of H
   nrm = H(j + 1, j);
   H(j + 1, j) = 0;

   % Compute the exponential of the augmented matrix
   F = expm_kiops(sgn * tau * HK);
   exps = exps + 1;

   % Restore the value of H_{m+1,m}
   H(j + 1, j) = nrm;

   if happy

      % Happy breakdown; wrap up
      omega   = 0;
      happy   = false;
      m_new   = m;
      tau_new = min(tau_end - (tau_now + tau), tau);

   else

      % Local a-posteriori error estimate, cf. [1, Subsection 4.3]
      w_ = F(1:j,end)/tau;
      u_ = K(1:j,1:j)\w_;
      err = abs(beta * tau * nrm * abs(u_(end)) * norm(V(:,j+1)));

      % Error for this step
      oldomega = omega;
      omega = tau_end * err / (tau * tol);
      
      % Estimate order
      if m == oldm && tau ~= oldtau && ireject >= 1
         order = max(1, log(omega/oldomega) / log(tau/oldtau));
         orderold = false;
      elseif orderold || ireject == 0
         orderold = true;
         order = j/4;
      else
         orderold = true;
      end
      % Estimate k
      if m ~= oldm && tau == oldtau && ireject >= 1
         kest = max(1.1, (omega/oldomega)^(1/(oldm-m)));
         kestold = false;
      elseif kestold || ireject == 0
         kestold = true;
         kest = 2;
      else
         kestold = true;
      end

      if omega > delta
         remaining_time = tau_end - tau_now;
      else
         remaining_time = tau_end - (tau_now + tau);
      end

      % Krylov adaptivity

      same_tau = min(remaining_time, tau);
      tau_opt  = tau * (gamma / omega)^(1 / order);
      tau_opt  = min(remaining_time, max(tau/5, min(5*tau, tau_opt)));

      m_opt = ceil(j + log(omega / gamma) / log(kest));
      m_opt = max(mmin, min(mmax, max(floor(3/4*m), min(m_opt, ceil(4/3*m)))));

      if j == mmax
         if omega > delta
            m_new = j;
            tau_new = tau * (gamma_mmax / omega)^(1 / order);
            tau_new = min(tau_end - tau_now, max(tau/5, tau_new));
         else
            tau_new = tau_opt;
            m_new = m;
         end
      else
         m_new = m_opt;
         tau_new = same_tau;
      end

   end

   % Check error against target
   if omega <= delta

      % Yep, got the required tolerance; update
      reject = reject + ireject;
      step = step + 1;

      % Udate for tau_out in the interval (tau_now, tau_now + tau)
      blownTs = 0;
      nextT = tau_now + tau;
      for k = l:numSteps
         if abs(tau_out(k)) < abs(nextT)
            blownTs = blownTs + 1;
         end
      end

      if blownTs ~= 0
         % Copy current w to w we continue with.
         w(:,l + blownTs) = w(:,l);

         for k = 0:blownTs - 1
            tauPhantom = tau_out(l+k) - tau_now;
            F2 = expm_kiops(sgn * tauPhantom * HK(1:j, 1:j));
            w(:, l+k) = beta * V(1:n, 1:j) * F2(1:j, 1);
         end

         % Advance l.
         l = l + blownTs;
      end

      % Using the standard scheme
      w(1:n, l) = beta * V(1:n, 1:j) * F(1:j, 1);

      % Update tau_out
      tau_now = tau_now + tau;

      j = 0;
      ireject = 0;

   else

      % Nope, try again
      ireject = ireject + 1;

      % Restore the original matrix
      H(1, j + 1) = 0;
   end

   oldtau = tau;
   tau    = tau_new;

   oldm = m;
   m    = m_new;

end

% Warn if we're trying to use more poles than are available and continue
% with infinity poles
if m>=mmax
   warning('Poles exhausted - continued with polynomial Krylov iterations.')
end

if tau_out(1)~=1 && task1
   if length(tau_out)==1
      w(:,l) = w(:,l)*(1/tau_out(l))^p;
   else
      phiHan = find(max(abs(u)));

      if isempty(phiHan)
         phiHan = length(u);
      else
         phiHan = phiHan-1;
      end

      for l = 1:numSteps
         w(:,l) = w(:,l)*(1/tau_out(l).^phiHan);
      end
   end
end

m_ret=m;

stats = [step, reject, krystep, exps];

end
