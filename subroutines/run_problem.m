function u = run_problem(pde, exptAb_routine, exp_rk_int, xi, xi_unique, dt, T_array, linear_system_solver, n, network, compute_error_to_ode15s_reference_solution, plot_solution)
% Efficiently solves a stiff system of ODEs by exponential integrators,
% cf. [1,2,3,4,5].
% 
% Input:
% pde:                   string of prepared problems (choose from: 'AC_2D',
%                        'GM_2D', 'AC_network', 'GM_network')
% exptAb_routine:        string of implemented routines to adaptively
%                        approximate the action of matrix exponentials on
%                        certain vectors (choose from: 'phipm', 'kiops', 'rk2expint')
% exp_rk_int:            string of implemented exponential Runge--Kutta
%                        integrators (choose from: 'SW2', 'ETD3RK', 'Krogstad4')
% linear_system_solver:  string of strategy to solve linear systems arising
%                        in rational Krylov subspace methods (choose from: 
%                        'lu_Matlab', 'lu_pardiso', 'AGMG')
% network:               string of prepared example networks serving as
%                        inherently discrete domains on which to solve the
%                        ODE system (choose from: 'minnesota',
%                        'usroads_subset', 'ak', 'luxembourg_osm', 'ny',
%                        'roadNet_PA', 'brightkite')
% n:                     number of grid points per dimension (for
%                        pde=AC_2D or pde=GM_2D)
% compute_error_to_ode15s_reference_solution, plot_solution:
%                        boolean and self-explaining. For ode15s, the default
%                        tolerance is an absolute and relative tolerance of
%                        1e-06 that can be modified in this script by means
%                        of the variable ode15s_tol.
% 
% Output:
%       u:  trajectory of the computed solution to the system of ODEs
% 
% 
% References:
% [1] K. Bergermann and M. Stoll, Adaptive rational Krylov methods for exponential Runge--Kutta integrators, arxiv preprint arXiv:2303.09482, (2023).
% [2] Gaudreault, S., Rainwater, G. and Tokman, M., 2018. KIOPS: A fast adaptive Krylov subspace solver for exponential integrators. Journal of Computational Physics.
% [3] Niesen, J. and Wright, W.M., 2011. A Krylov subspace method for option pricing. SSRN 1799124
% [4] Niesen, J. and Wright, W.M., 2012. Algorithm 919: A Krylov subspace algorithm for evaluating the \varphi-functions appearing in exponential integrators. ACM Transactions on Mathematical Software (TOMS), 38(3), p.22
% [5] A. H. Al-Mohy and N. J. Higham, Computing the action of the matrix exponential, with an application to exponential integrators, SIAM J. Sci. Comput., 33 (2011), pp. 488–511.
% 

if nargin<12
    plot_solution = true;
    if nargin<11
        compute_error_to_ode15s_reference_solution = true;
        if nargin<10
            network = 'minnesota';
            warning('Network not recognized, using default value minnesota.')
            if nargin<9
                n=100;
                warning('Problem size n not specified, using default value 100.')
                if nargin<8
                    linear_system_solver = 'lu_Matlab';
                    warning('linear_system_solver not regconized, using default value lu_Matlab.')
                    if nargin<7
                        dt = 0.1;
                        warning('Time step size dt not specified, using default value 0.1.')
                        if nargin<6
                            T_array = 0:0.1:1;
                            warning('T_array not specified, using default value 0:0.1:1.')
                            if nargin<5
                                load('pole_files/expint_poles.mat');
                                warning('Unique poles xi_unique not specified, using default values from RKFIT poles.')
                                if nargin<4
                                    load('pole_files/expint_poles.mat');
                                    warning('Poles xi not specified, using default values from RKFIT poles.')
                                    if nargin<3
                                        exp_rk_int = 'ETD3RK';
                                        warning('exp_rk_int not recognized, using default value ETD3RK.')
                                        if nargin<2
                                            exptAb_routine = 'rk2expint';
                                            warning('exptAb_routine not recognized, using default value rk2expint.')
                                            if nargin<1
                                                error('Specify problem to be solved!')
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
    end
end

fprintf('#########################################\nSetting up problem...')

%% Set up problem
[Aop, u0, g, rhs, J, adj, usroads_coords] = setup_problem(pde, n, network);
n = size(Aop,1);

fprintf('done.\n\n')

% print basic information about problem
if strcmp(exptAb_routine,'rk2expint')
    if strcmp(pde,'AC_network')
        fprintf('Problem:                %s\nNetwork:                %s\nn =                     %d\nExponential integrator: %s\nexp(-tA)b routine:      %s\nLinear system solver:   %s\n\n', pde(1:end-8), network, n, exp_rk_int, exptAb_routine, linear_system_solver)
    else
        fprintf('Problem:                %s\nn =                     %d\nExponential integrator: %s\nexp(-tA)b routine:      %s\nLinear system solver:   %s\n\n', pde, n, exp_rk_int, exptAb_routine, linear_system_solver)
    end
else
    if strcmp(pde,'AC_network') || strcmp(pde,'GS_network') || strcmp(pde,'GM_network')
        fprintf('Problem:                %s\nNetwork:                %s\nn =                     %d\nExponential integrator: %s\nexp(-tA)b routine:      %s\n\n', pde(1:end-8), network, n, exp_rk_int, exptAb_routine)
    else
        fprintf('Problem:                %s\nn =                     %d\nExponential integrator: %s\nexp(-tA)b routine:      %s\n\n', pde, n, exp_rk_int, exptAb_routine)
    end
end

%% Precompute factorizations in case of direct solvers
if strcmp(exptAb_routine,'rk2expint') && (strcmp(linear_system_solver,'lu_Matlab') || strcmp(linear_system_solver,'lu_pardiso'))
    factorization = compute_factorizations(Aop, dt, xi_unique, linear_system_solver, exp_rk_int);
else
    factorization.xi_unique = xi_unique;
end

%% Solve problem
args = struct;
args.exptAb_routine = exptAb_routine;
args.xi = xi;
args.linear_system_solver = linear_system_solver;
args.factorization = factorization;

switch exp_rk_int
    case 'SW2'
        u = SW2(g, T_array, u0(:), Aop, args);
    case 'ETD3RK'
        u = ETD3RK(g, T_array, u0(:), Aop, args);
    case 'Krogstad4'
        u = Krogstad4(g, T_array, u0(:), Aop, args);
end

if compute_error_to_ode15s_reference_solution
    % Options for Matlab ODE solver
    ode15s_tol = 1e-06;
    options = odeset('Jacobian', J,'RelTol',ode15s_tol,'AbsTol',ode15s_tol,'Stats','off');

    fprintf('\n\node15s for reference solution using RelTol=AbsTol=%0.e...\n', ode15s_tol)
    
    tic
    [~, u_ode15] = ode15s(rhs, [0,T_array(end)], u0(:), options);
    runtime = toc;
    fprintf('\nRunning ode15s took %.2f seconds.\n', runtime)
    u_ode15 = u_ode15';

    fprintf('\nAbsolute infinity norm error between %s and ode15s:\n %e\n\n', exp_rk_int,  abs(max(u(:,end)-u_ode15(:,end))))
end


%% Plot solution
if plot_solution
    if strcmp(pde,'AC_2D')
        figure(1)
        set(gcf,'Position',[100 100 650 500])
        for i = 1:length(T_array)
            imagesc(reshape(real(u(:,i)),[round(sqrt(n)),round(sqrt(n))]))
%             colorbar
            sgtitle(['Time ',num2str(T_array(i))])
            pause(.0001)
        end
    elseif strcmp(pde,'GM_2D')
        figure(1)
        set(gcf,'Position',[100 100 1500 500])
        for i = 1:length(T_array)
            subplot(121)
            imagesc(reshape(real(u(1:round(n/2),i)),[round(sqrt(n/2)),round(sqrt(n/2))]))
%             colorbar
            subplot(122)
            imagesc(reshape(real(u(round(n/2)+1:end,i)),[sqrt(n/2),sqrt(n/2)]))
%             colorbar
            sgtitle(['Time ',num2str(T_array(i))])
            pause(.0001)
        end
    elseif strcmp(pde,'AC_network') && (strcmp(network,'brightkite') || strcmp(network,'minnesota') || strcmp(network,'ak') || strcmp(network,'luxembourg_osm') || strcmp(network,'ny') || strcmp(network,'roadNet_PA'))
        mSize = 2;
        
        G = graph(adj);
        set(gcf,'Position',[100 100 650 500])
        for i = 1:floor(length(T_array))
            G.Nodes.NodeColors = real(u(:,i));
            plot(G,'NodeCData',G.Nodes.NodeColors, 'MarkerSize', mSize, 'Marker', 'o', 'EdgeColor', [.7 .7 .7]);
            colorbar
            caxis([-1 1]);
            sgtitle(['Time ',num2str(T_array(i))])
            pause(.0001)
        end
    elseif strcmp(pde,'AC_network') && strcmp(network,'usroads_subset')
        G = graph(adj);
        set(gcf,'Position',[100 100 650 500])
        for i = 1:floor(length(T_array))
            G.Nodes.NodeColors = real(u(:,i));
            plot(G,'NodeCData',G.Nodes.NodeColors, 'MarkerSize', 8, 'Marker', 'o','XData',usroads_coords(:,1),'YData',usroads_coords(:,2));
%             colorbar
            sgtitle(['Time ',num2str(T_array(i))])
            pause(.0001)
        end
    elseif strcmp(pde,'GM_network') && (strcmp(network,'brightkite') || strcmp(network,'minnesota') || strcmp(network,'ak') || strcmp(network,'luxembourg_osm') || strcmp(network,'ny') || strcmp(network,'roadNet_PA'))
        mSize = 4;
        G = graph(adj);
        set(gcf,'Position',[100 100 1500 500])
        for i = 1:floor(length(T_array))
            subplot(121)
            G.Nodes.NodeColors = real(u(1:n/2,i));
            plot(G,'NodeCData',G.Nodes.NodeColors, 'MarkerSize', mSize, 'Marker', 'o', 'EdgeColor', [.7 .7 .7]);
            colorbar
            subplot(122)
            G.Nodes.NodeColors = real(u(n/2+1:end,i));
            plot(G,'NodeCData',G.Nodes.NodeColors, 'MarkerSize', mSize, 'Marker', 'o', 'EdgeColor', [.7 .7 .7]);
            colorbar
            sgtitle(['Time ',num2str(T_array(i))])
            pause(.0001)
        end
    elseif strcmp(pde,'GM_network') && strcmp(network,'usroads_subset')
        mSize = 4;
        G = graph(adj);
        set(gcf,'Position',[100 100 1500 500])
        for i = 1:floor(length(T_array))
            subplot(121)
            G.Nodes.NodeColors = real(u(1:n/2,i));
            plot(G,'NodeCData',G.Nodes.NodeColors, 'MarkerSize', mSize, 'Marker', 'o','XData',usroads_coords(:,1),'YData',usroads_coords(:,2));
            colorbar
            subplot(122)
            G.Nodes.NodeColors = real(u(n/2+1:end,i));
            plot(G,'NodeCData',G.Nodes.NodeColors, 'MarkerSize', mSize, 'Marker', 'o','XData',usroads_coords(:,1),'YData',usroads_coords(:,2));
            colorbar
            sgtitle(['Time ',num2str(T_array(i))])
            pause(.0001)
        end
    end
end
