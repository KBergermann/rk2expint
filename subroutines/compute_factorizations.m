function factorization = compute_factorizations(Aop, dt, xi_unique, linear_system_solver, exp_rk_int)
% Precomputes LU or Cholesky decompositions of the matrices (xi*I-t*Aop)
% required in the direct solution of the block linear systems.
% 
% Input:
%       Aop:        discrete linear differential operator (sparse matrix),
%                   potentially including some multiplicative constants
%       dt:         time (sub-)step size in matrix specified above
%       xi_unique:  1D-array of unique poles (for repeated poles, the
%                   decomposition need only be computed once)
%       linear_system_solver:
%                   string, 'lu_Matlab' or 'lu_pardiso'
%       exp_rk_int: string, 'SW2', 'ETD3RK', or 'Krogstad4' (SW2 requires
%                   only one decomposition for t=1, ETD3RK requires two
%                   decompositions for t=1 and t=dt, and Krogstad4 requires
%                   three decompositions for t=1, t=dt, and t=dt/2).
%                   Indexing, {3*i-1} etc. is used to map the correct xi and t.
% 
% Output:
%       factorization: struct containing necessary ingredients of the
%                      required decompositions (depends on linear_system_solver)
% 

    % sets up factorization struct leaving space for the three
    % decompositions for t=1, t=dt, and t=dt/2
    factorization = struct;
    factorization.xi_unique = xi_unique;
    L = cell(3*length(xi_unique),1); U = cell(2*length(xi_unique),1); 
    p = cell(3*length(xi_unique),1);
    info = cell(3*length(xi_unique),1);

    bar=waitbar(0,'Precomputing matrix factorizations');
    tic
    
    for i=1:length(xi_unique)
        if strcmp(linear_system_solver,'lu_Matlab')
            % t=1 decomposition
            Aop_shifted = Aop - xi_unique(i)*speye(size(Aop));
            p{3*i-2} = amd(Aop_shifted);
            [L{3*i-2}, U{3*i-2}] = lu(Aop_shifted(p{3*i-2},p{3*i-2}));
            
            % ETD3RK and Krogstad4 contain a t=dt decomposition
            if strcmp(exp_rk_int,'ETD3RK') || strcmp(exp_rk_int,'Krogstad4')
                Aop_shifted = dt*Aop - xi_unique(i)*speye(size(Aop));
                p{3*i-1} = amd(Aop_shifted);
                [L{3*i-1}, U{3*i-1}] = lu(Aop_shifted(p{3*i-1},p{3*i-1}));
            else
                p{3*i-1} = 1; L{3*i-1} = 1; U{3*i-1} = 1;
            end
            
            % Krogstad4 additionally contains a t=dt/2 decomposition
            if strcmp(exp_rk_int,'Krogstad4')
                Aop_shifted = (dt/2)*Aop - xi_unique(i)*speye(size(Aop));
                p{3*i} = amd(Aop_shifted);
                [L{3*i}, U{3*i}] = lu(Aop_shifted(p{3*i},p{3*i}));
            else
                p{3*i} = 1; L{3*i} = 1; U{3*i} = 1;
            end

        elseif strcmp(linear_system_solver,'lu_pardiso')
            % t=1 decomposition
            Aop_shifted = Aop - xi_unique(i)*speye(size(Aop));
            % query complex-valuedness of poles and symmetry of Aop_shifted
            % to initialize the corresponding pardiso instance
            if ~isreal(xi_unique)
                if issymmetric(Aop_shifted)
                    info{3*i-2} = pardisoinit(6,0);
                    info{3*i-2} = pardisoreorder(tril(Aop_shifted),info{3*i-2},true);
                    info{3*i-2} = pardisofactor(tril(Aop_shifted),info{3*i-2},true);
                else
                    info{3*i-2} = pardisoinit(13,0);
                    info{3*i-2} = pardisoreorder(Aop_shifted,info{3*i-2},true);
                    info{3*i-2} = pardisofactor(Aop_shifted,info{3*i-2},true);
                end
            else
                if issymmetric(Aop_shifted)
                    info{3*i-2} = pardisoinit(2,0);
                    info{3*i-2} = pardisoreorder(tril(Aop_shifted),info{3*i-2},true);
                    info{3*i-2} = pardisofactor(tril(Aop_shifted),info{3*i-2},true);
                else
                    info{3*i-2} = pardisoinit(11,0);
                    info{3*i-2} = pardisoreorder(Aop_shifted,info{3*i-2},true);
                    info{3*i-2} = pardisofactor(Aop_shifted,info{3*i-2},true);
                end
            end
            
            % ETD3RK and Krogstad4 contain a t=dt decomposition
            if strcmp(exp_rk_int,'ETD3RK') || strcmp(exp_rk_int,'Krogstad4')
                Aop_shifted = dt*Aop - xi_unique(i)*speye(size(Aop));
                % query complex-valuedness of poles and symmetry of Aop_shifted
                % to initialize the corresponding pardiso instance
                if ~isreal(xi_unique)
                    if issymmetric(Aop_shifted)
                        info{3*i-1} = pardisoinit(6,0);
                        info{3*i-1} = pardisoreorder(tril(Aop_shifted),info{3*i-1},true);
                        info{3*i-1} = pardisofactor(tril(Aop_shifted),info{3*i-1},true);
                    else
                        info{3*i-1} = pardisoinit(13,0);
                        info{3*i-1} = pardisoreorder(Aop_shifted,info{3*i-1},true);
                        info{3*i-1} = pardisofactor(Aop_shifted,info{3*i-1},true);
                    end
                else
                    if issymmetric(Aop_shifted)
                        info{3*i-1} = pardisoinit(2,0);
                        info{3*i-1} = pardisoreorder(tril(Aop_shifted),info{3*i-1},true);
                        info{3*i-1} = pardisofactor(tril(Aop_shifted),info{3*i-1},true);
                    else
                        info{3*i-1} = pardisoinit(11,0);
                        info{3*i-1} = pardisoreorder(Aop_shifted,info{3*i-1},true);
                        info{3*i-1} = pardisofactor(Aop_shifted,info{3*i-1},true);
                    end
                end
            else
                % placeholder for SW2 to avoid uninitialized objects
                if ~isreal(xi_unique)
                    if issymmetric(Aop_shifted)
                        info{3*i-1} = pardisoinit(6,0);
                    else
                        info{3*i-1} = pardisoinit(13,0);
                    end
                else
                    if issymmetric(Aop_shifted)
                        info{3*i-1} = pardisoinit(2,0);
                    else
                        info{3*i-1} = pardisoinit(11,0);
                    end
                end
            end
            
            % Krogstad4 additionally contains a t=dt/2 decomposition
            if strcmp(exp_rk_int,'Krogstad4')
                Aop_shifted = (dt/2)*Aop - xi_unique(i)*speye(size(Aop));
                % query complex-valuedness of poles and symmetry of Aop_shifted
                % to initialize the corresponding pardiso instance
                if ~isreal(xi_unique)
                    if issymmetric(Aop_shifted)
                        info{3*i} = pardisoinit(6,0);
                        info{3*i} = pardisoreorder(tril(Aop_shifted),info{3*i},true);
                        info{3*i} = pardisofactor(tril(Aop_shifted),info{3*i},true);
                    else
                        info{3*i} = pardisoinit(13,0);
                        info{3*i} = pardisoreorder(Aop_shifted,info{3*i},true);
                        info{3*i} = pardisofactor(Aop_shifted,info{3*i},true);
                    end
                else
                    if issymmetric(Aop_shifted)
                        info{3*i} = pardisoinit(2,0);
                        info{3*i} = pardisoreorder(tril(Aop_shifted),info{3*i},true);
                        info{3*i} = pardisofactor(tril(Aop_shifted),info{3*i},true);
                    else
                        info{3*i} = pardisoinit(11,0);
                        info{3*i} = pardisoreorder(Aop_shifted,info{3*i},true);
                        info{3*i} = pardisofactor(Aop_shifted,info{3*i},true);
                    end
                end
            else
                % placeholderfor SW2 and ETD3RK to avoid uninitialized objects
                if ~isreal(xi_unique)
                    if issymmetric(Aop_shifted)
                        info{3*i} = pardisoinit(6,0);
                    else
                        info{3*i} = pardisoinit(13,0);
                    end
                else
                    if issymmetric(Aop_shifted)
                        info{3*i} = pardisoinit(2,0);
                    else
                        info{3*i} = pardisoinit(11,0);
                    end
                end
            end
        end
        waitbar(i/(length(xi_unique)-1),bar);
    end
    runtime = toc;
    fprintf('Computing matrix factorizations took %.4f seconds.\n', runtime)
    close(bar)
    
    % save computed decompositions to factorization struct
    if strcmp(linear_system_solver,'lu_Matlab')
        factorization.L = L;
        factorization.U = U;
        factorization.p = p;
    elseif strcmp(linear_system_solver,'lu_pardiso')
        factorization.info = info;
    end
end
