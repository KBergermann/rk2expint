function y = solve_agmg_mult(A, C, Jp, nu, mu, b)
% Solves the block linear system [1, Subsection 4.2] by backsubstitution
% and a flexible conjugate gradient preconditioned with aggregation-based
% multigrid (AGMG)
% 
% Input:
%       A:      discrete linear differential operator (sparse matrix)
%       C:      trajectory-dependent block vector, cf. [1, Theorem 3.2]
%       Jp:     Jordan block to eigenvalue 0
%       nu, mu: Moebius transformation of pole xi
%       b:      right-hand-side vector
% 
% Output:
%       y:      solution vector of the block linear system.
% 
% Reference:
% [1] K. Bergermann and M. Stoll, Adaptive rational Krylov methods for exponential Runge--Kutta integrators, arxiv preprint arXiv:2303.09482, (2023).
% 

    n = size(A,1);
    p = size(Jp,1);

    if nu==0 && mu==1
        % linear system solve with identity
        y = -b;
    else
        % get rid of potential complex-valued roundoff error
        if isreal(mu)
            b = real(b);
        end
        % solution of small bottom p-by-p system
        y2 = (-mu*eye(p)-Jp)\(b(n+1:end));
        % solution of [Equation (4.5), 1]
        agmg(-mu*speye(size(A))-nu*A,[],0,[],[],[],[],1);
        y1 = agmg(-mu*speye(size(A)) - nu*A, b(1:n) + C*y2, 0, 1e-07, 50, [], [], 2);
        y = [y1; y2];
    end
end