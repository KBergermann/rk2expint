function x = solve_lu_pardiso(A, Jp, C, nu, mu, b, tau_out, factorization, dt_2_flag)
% Solves the block linear system [1, Subsection 4.2] by backsubstitution
% and an LU or Cholesky decomposition precomputed by pardiso.
% 
% Input:
%       A:              discrete linear differential operator (sparse matrix)
%       Jp:             Jordan block to eigenvalue 0
%       C:              trajectory-dependent block vector, cf. [1, Theorem 3.2]
%       nu, mu:         Moebius transformation of pole xi
%       b:              right-hand-side vector
%       tau_out:        if ~=1, sub-time interval steps are performed
%                       (a decomposition of (xi*I-tau*A) with tau~=1 is required)
%       factorization:  factors L and U as well as amd permutation vector p
%                       of the discrete linear differential operator A
%       dt_2_flag:      boolean, if true, decomposition of (xi*I-(tau/2)*A)
%                       is required
% 
% Output:
%       x:      solution vector of the block linear system.
% 
% Reference:
% [1] K. Bergermann and M. Stoll, Adaptive rational Krylov methods for exponential Runge--Kutta integrators, arxiv preprint arXiv:2303.09482, (2023).
% 

    n = size(A,1);
    p = size(Jp,1);
    
    if nu==0 && mu==1
        % linear system solve with identity
        x = -b;
    else
        % solution of small bottom p-by-p system
        x2 = (-mu*eye(p)-Jp)\(b(n+1:end));
        % solution of [Equation (4.5), 1]. Depending on dt, the correct
        % decomposition must be chosen.
        if tau_out ~=1
            if issymmetric(A)
                [x1, ~] = pardisosolve(tril(-nu*A - mu*speye(size(A))), b(1:n) + C*x2, factorization.info{3*find(factorization.xi_unique==mu)-2}, true);
            else
                [x1, ~] = pardisosolve(-nu*A - mu*speye(size(A)), b(1:n) + C*x2, factorization.info{3*find(factorization.xi_unique==mu)-2}, true);
            end
        else
            if ~dt_2_flag
                if issymmetric(A)
                    [x1, ~] = pardisosolve(tril(-nu*A - mu*speye(size(A))), b(1:n) + C*x2, factorization.info{3*find(factorization.xi_unique==mu)-1}, true);
                else
                    [x1, ~] = pardisosolve(-nu*A - mu*speye(size(A)), b(1:n) + C*x2, factorization.info{3*find(factorization.xi_unique==mu)-1}, true);
                end
            else
                if issymmetric(A)
                    [x1, ~] = pardisosolve(tril(-nu*A - mu*speye(size(A))), b(1:n) + C*x2, factorization.info{3*find(factorization.xi_unique==mu)}, true);
                else
                    [x1, ~] = pardisosolve(-nu*A - mu*speye(size(A)), b(1:n) + C*x2, factorization.info{3*find(factorization.xi_unique==mu)}, true);
                end
            end
        end
        x = [x1; x2];
    end
end
