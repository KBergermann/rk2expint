function u = Krogstad4(g, T_array, u0, A, args)
% Krogstad exponential integrator of stage 4 and stiff order 4 [1] that can
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
% [1] S. Krogstad, Generalized integrating factor methods for stiff PDEs, J. Comput. Phys., 203 (2005), pp. 72–88.
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
    Krylov_iter = zeros(length(T_array)-1,4);
    
    t = T_array(1);
    bar=waitbar(0,'Time integration by Krogstad4');
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
            [w2, stats] = phipm(1, -(h/2)*A, [zeros(n,1), (1/2)*(Gn1-Aun), -(Gn1-Aun) + (Gn2-Aun)], 1e-08);
            Krylov_iter(i,2) = stats(3);
            Un3 = u(:,i) + h*w2;
            Gn3 = g(T_array(i) + h/2, Un3);

            % stage 3
            [w3, stats] = phipm(1, -h*A, [zeros(n,1), (Gn1-Aun), -2*(Gn1-Aun) + 2*(Gn3-Aun)], 1e-08);
            Krylov_iter(i,3) = stats(3);
            Un4 = u(:,i) + h*w3;
            Gn4 = g(T_array(i) + h, Un4);

            % stage 4
            [w4, stats] = phipm(1, -h*A, [zeros(n,1), Gn1-Aun, -3*(Gn1-Aun) + 2*(Gn2-Aun) + 2*(Gn3-Aun) - (Gn4-Aun), 4*(Gn1-Aun) - 4*(Gn2-Aun) - 4*(Gn3-Aun) + 4*(Gn4-Aun)], 1e-08);
            Krylov_iter(i,4) = stats(3);
            u(:,i+1) = u(:,i) + h*w4;

            t = t+h;
            waitbar(i/(length(T_array)-1),bar);
        end
        runtime = toc;
        fprintf('\nphipm Krogstad4 took %.2f seconds for %d time steps, i.e., %.4f seconds per time step.\n', runtime, (length(T_array)-1), runtime/(length(T_array)-1))
        close(bar)

        fprintf('\nAverage numbers of polynomial Krylov iterations per phipm call:\n %.0f  %.0f  %.0f  %.0f\n\n', mean(Krylov_iter,1))
    
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
            [w2, ~, stats] = kiops(1, -(h/2)*A, [zeros(n,1), (1/2)*(Gn1-Aun), -(Gn1-Aun) + (Gn2-Aun)], 1e-08, 10, 10, 128, false);
            Krylov_iter(i,2) = stats(3);
            Un3 = u(:,i) + h*w2;
            Gn3 = g(T_array(i) + h/2, Un3);

            % stage 3
            [w3, ~, stats] = kiops(1, -h*A, [zeros(n,1), (Gn1-Aun), -2*(Gn1-Aun) + 2*(Gn3-Aun)], 1e-08, 10, 10, 128, false);
            Krylov_iter(i,3) = stats(3);
            Un4 = u(:,i) + h*w3;
            Gn4 = g(T_array(i) + h, Un4);

            % stage 4
            [w4, ~, stats] = kiops(1, -h*A, [zeros(n,1), Gn1-Aun, -3*(Gn1-Aun) + 2*(Gn2-Aun) + 2*(Gn3-Aun) - (Gn4-Aun), 4*(Gn1-Aun) - 4*(Gn2-Aun) - 4*(Gn3-Aun) + 4*(Gn4-Aun)], 1e-08, 10, 10, 128, false);
            Krylov_iter(i,4) = stats(3);
            u(:,i+1) = u(:,i) + h*w4;

            t = t+h;
            waitbar(i/(length(T_array)-1),bar);
        end
        runtime = toc;
        fprintf('\nkiops Krogstad4 took %.2f seconds for %d time steps, i.e., %.4f seconds per time step.\n', runtime, (length(T_array)-1), runtime/(length(T_array)-1))
        close(bar)

        fprintf('\nAverage numbers of polynomial Krylov iterations per kiops call:\n %.0f  %.0f  %.0f  %.0f\n\n', mean(Krylov_iter,1))
    
    elseif strcmp(exptAb_routine,'rk2expint')
        for i=1:(length(T_array)-1)
            % current time step size
            h = T_array(i+1) - T_array(i);

            % preparations for block vectors
            Gn1 = g(T_array(i), u(:,i));
            Aun = A*u(:,i);

            % stage 1
            [w1, ~, stats] = rk2expint(h/2, -A, real([zeros(n,1), Gn1 - Aun]), 1e-08, 5, 5, length(xi), true, xi, linear_system_solver, factorization);
            Krylov_iter(i,1) = stats(3);
            Un2 = u(:,i) + (h/2)*w1;
            Gn2 = g(T_array(i) + (1/2)*h, Un2);

            % stage 2
            [w2, ~, stats] = rk2expint(1, -(h/2)*A, real([zeros(n,1), (1/2)*(Gn1-Aun), -(Gn1-Aun) + (Gn2-Aun)]), 1e-08, 5, 5, length(xi), true, xi, linear_system_solver, factorization, true);
            Krylov_iter(i,2) = stats(3);
            Un3 = u(:,i) + h*w2;
            Gn3 = g(T_array(i) + h/2, Un3);

            % stage 3
            [w3, ~, stats] = rk2expint(1, -h*A, real([zeros(n,1), (Gn1-Aun), -2*(Gn1-Aun) + 2*(Gn3-Aun)]), 1e-08, 5, 5, length(xi), true, xi, linear_system_solver, factorization);
            Krylov_iter(i,3) = stats(3);
            Un4 = u(:,i) + h*w3;
            Gn4 = g(T_array(i) + h, Un4);

            % stage 4
            [w4, ~, stats] = rk2expint(1, -h*A, real([zeros(n,1), Gn1-Aun, -3*(Gn1-Aun) + 2*(Gn2-Aun) + 2*(Gn3-Aun) - (Gn4-Aun), 4*(Gn1-Aun) - 4*(Gn2-Aun) - 4*(Gn3-Aun) + 4*(Gn4-Aun)]), 1e-08, 5, 5, length(xi), false, xi, linear_system_solver, factorization);
            Krylov_iter(i,4) = stats(3);
            u(:,i+1) = u(:,i) + h*w4;

            t = t+h;
            waitbar(i/(length(T_array)-1),bar);
        end
        runtime = toc;
        fprintf('\nrk2expint Krogstad4 took %.2f seconds for %d time steps, i.e., %.4f seconds per time step.\n', runtime, (length(T_array)-1), runtime/(length(T_array)-1))
        close(bar)

        fprintf('\nAverage numbers of rational Krylov iterations per rk2expint call:\n %.0f  %.0f  %.0f  %.0f\n\n', mean(Krylov_iter,1))
    end
end
