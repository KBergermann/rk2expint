function u = ETD3RK(g, T_array, u0, A, args)
% ETD3RK exponential integrator of stage 3 and stiff order 3 [1] that can
% be called similarly to Matlab's ODE solvers.
% 
% Input:
%       g:          function handle representing the non-linearity
%       T_array:    1D-array with time points
%       u0:         vector of initial conditions
%       A:          discrete linear differential operator (sparse matrix)
%       args:       additional arguments:
%           exptAb_routine:         string: 'phipm', 'kiops', or 'rk2expint'
%           xi:                     1D array of poles
%           linear_system_solver:   string: 'lu_Matlab', 'lu_pardiso', or 'AGMG'
%           factorization:          struct containing precomputed matrix decompositions of A
% 
% Output:
%       u:  trajectory of the computed solution.
% 
% Reference:
% [1] S. M. Cox and P. C. Matthews, Exponential time differencing for stiff systems, J. Comput. Phys., 176 (2002), pp. 430–455.
% 

    if nargin<5
        exptAb_routine = 'rk2expint';
        load('pole_files/expint_poles.mat');
        linear_system_solver = 'lu_Matlab';
    else
        exptAb_routine = args.exptAb_routine;
        xi = args.xi;
        linear_system_solver = args.linear_system_solver;
        factorization = args.factorization;
    end

    n = size(A,1);
    u = zeros(n,length(T_array)-1);
    u(:,1) = u0(:);

    % keep track of Krylov iteration numbers
    Krylov_iter = zeros(length(T_array)-1,3);
        
    t = T_array(1);
    bar=waitbar(0,'Time integration by ETD3RK');
    tic
    
    if strcmp(exptAb_routine,'phipm')
    
        for i=1:(length(T_array)-1)
            % current time step size
            h = T_array(i+1) - T_array(i);
            
            % preparations for block vectors
            Gn1 = g(T_array(i), u(:,i));
            Aun = A*u(:,i);
            
            % stage 1
            [w1, stats] = phipm(1, -(h/2)*A, [zeros(n,1), Gn1 - Aun], 1e-08);
            Krylov_iter(i,1) = stats(3);
            Un2 = u(:,i) + (h/2)*w1;
            Gn2 = g(T_array(i) + (1/2)*h, Un2);
            
            % stage 2
            [w2, stats] = phipm(1, -h*A, [zeros(n,1), -(Gn1-Aun) + 2*(Gn2-Aun)], 1e-08);
            Krylov_iter(i,2) = stats(3);
            Un3 = u(:,i) + h*w2;
            Gn3 = g(T_array(i) + h, Un3);
            
            % stage 3
            [w3, stats] = phipm(1, -h*A, [zeros(n,1), Gn1-Aun, -3*(Gn1-Aun) + 4*(Gn2-Aun) - (Gn3-Aun), 4*(Gn1-Aun) - 8*(Gn2-Aun) + 4*(Gn3-Aun)], 1e-08);
            Krylov_iter(i,3) = stats(3);
            u(:,i+1) = u(:,i) + h*w3;
            
            t = t+h;
            waitbar(i/(length(T_array)-1),bar);
        end
        runtime = toc;
        fprintf('\nphipm ETD3RK took %.2f seconds for %d time steps, i.e., %.4f seconds per time step.\n', runtime, (length(T_array)-1), runtime/(length(T_array)-1))
        close(bar)
        
        fprintf('\nAverage numbers of polynomial Krylov iterations per phipm call:\n %.0f  %.0f  %.0f\n\n', mean(Krylov_iter,1))

    elseif strcmp(exptAb_routine,'kiops')

        for i=1:(length(T_array)-1)
            % current time step size
            h = T_array(i+1) - T_array(i);
            
            % preparations for block vectors
            Gn1 = g(T_array(i), u(:,i));
            Aun = A*u(:,i);
            
            % stage 1
            [w1, ~, stats] = kiops(h/2, -A, [zeros(n,1), Gn1 - Aun], 1e-08);
            Krylov_iter(i,1) = stats(3);
            Un2 = u(:,i) + (h/2)*w1;
            Gn2 = g(T_array(i) + (1/2)*h, Un2);
            
            % stage 2
            [w2, ~, stats] = kiops(h, -A, [zeros(n,1), -(Gn1-Aun) + 2*(Gn2-Aun)], 1e-08);
            Krylov_iter(i,2) = stats(3);
            Un3 = u(:,i) + h*w2;
            Gn3 = g(T_array(i) + h, Un3);
            
            % stage 3
            [w3, ~, stats] = kiops(1, -h*A, [zeros(n,1), Gn1-Aun, -3*(Gn1-Aun) + 4*(Gn2-Aun) - (Gn3-Aun), 4*(Gn1-Aun) - 8*(Gn2-Aun) + 4*(Gn3-Aun)], 1e-08, 10, 10, 128, false);
            Krylov_iter(i,3) = stats(3);
            u(:,i+1) = u(:,i) + h*w3;
            
            t = t+h;
            waitbar(i/(length(T_array)-1),bar);
        end
       
        runtime = toc;
        fprintf('\nkiops ETD3RK took %.2f seconds for %d time steps, i.e., %.4f seconds per time step.\n', runtime, (length(T_array)-1), runtime/(length(T_array)-1))
        close(bar)
        
        fprintf('\nAverage numbers of polynomial Krylov iterations per kiops call:\n %.0f  %.0f  %.0f\n\n', mean(Krylov_iter,1))
     
    elseif strcmp(exptAb_routine,'rk2expint')

        for i=1:(length(T_array)-1)
            % current time step size
            h = T_array(i+1) - T_array(i);

            % preparations for block vectors
            Gn1 = g(T_array(i), u(:,i));
            Aun = A*u(:,i);

            % stage 1
            [w1, ~, stats] = rk2expint(h/2, -A, [zeros(n,1), Gn1 - Aun], 1e-08, 5, 5, length(xi), true, xi, linear_system_solver, factorization);
            Krylov_iter(i,1) = stats(3);
            Un2 = u(:,i) + (h/2)*w1;
            Gn2 = g(T_array(i) + (1/2)*h, Un2);

            % stage 2
            [w2, ~, stats] = rk2expint(h, -A, [zeros(n,1), -(Gn1-Aun) + 2*(Gn2-Aun)], 1e-08, 5, 5, length(xi), true, xi, linear_system_solver, factorization);
            Krylov_iter(i,2) = stats(3);
            Un3 = u(:,i) + h*w2;
            Gn3 = g(T_array(i) + h, Un3);

            % stage 3
            [w3, ~, stats] = rk2expint(1, -h*A, [zeros(n,1), Gn1-Aun, -3*(Gn1-Aun) + 4*(Gn2-Aun) - (Gn3-Aun), 4*(Gn1-Aun) - 8*(Gn2-Aun) + 4*(Gn3-Aun)], 1e-08, 5, 5, length(xi), false, xi, linear_system_solver, factorization);
            Krylov_iter(i,3) = stats(3);
            u(:,i+1) = u(:,i) + h*w3;

            t = t+h;
            waitbar(i/(length(T_array)-1),bar);
        end
        runtime = toc;
        fprintf('\nrk2expint ETD3RK took %.2f seconds for %d time steps, i.e., %.4f seconds per time step.\n', runtime, (length(T_array)-1), runtime/(length(T_array)-1))
        close(bar)

        fprintf('\nAverage numbers of rational Krylov iterations per rk2expint call:\n %.0f  %.0f  %.0f\n\n', mean(Krylov_iter,1))
    end
end
