function [Aop, u0, g, rhs, J, adj, usroads_coords] = setup_problem(pde, n, network)
% Sets up stiff systems of ODEs considered in [1, Section 6]. We consider
% the Allen--Cahn (AC) equation as well as the Gierer--Meinhardt (GM)
% equations on
% 1. finite difference discretizations of two-dimensional unit squares
% 2. graph Laplacians of inherently discrete network domains
% 
% Input:
%       pde:        string specifying the differential equation (AC or GM)
%                   and domain (2D or network)
%       n:          matrix size/number of degrees of freedom (total number
%                   of grid points or number of graph nodes)
%       network:    string indicating the network to be used as domain (if
%                   pde=AC_network or ode=GM_network)
% 
% Output:
%       Aop:            discrete linear differential operator (sparse matrix),
%                       potentially including some multiplicative constants
%       u0:             initial condition vector
%       g:              function handle of the non-linearity of the problem
%                       for use in the exponential integrator routines and ode15s
%       rhs:            function handle of the full right-hand side of the
%                       problem for use in ode15s
%       J:              function handle of the Jacobian of the problem for
%                       use in ode15s
%       adj:            (sparse) network adjacency matrix for plotting purposes
%       usroads_coords: n-by-2 coordinate array of the usroads_subset
%                       network for plotting purposes
% 
% Reference:
% [1] K. Bergermann and M. Stoll, Adaptive rational Krylov methods for exponential Runge--Kutta integrators, arxiv preprint arXiv:2303.09482, (2023).
% 
    
    switch pde
        case 'AC_2D'
            %%% Grid
            hx = 2/n; % mesh size
            n2 = n^2; % degrees of freedom

            %%% Initial conditions
            x = linspace(-1,1,n); y = linspace(-1,1,n); [X,Y] = meshgrid(x,y);
            u0 = 0.1*ones(n,n) + 0.1*cos(2*pi*X).*cos(2*pi*Y);

            %%% Interface parameter
            eps = .1;

            %%% Laplacian
            e = ones(n,1);
            Dx = spdiags([-e 2*e -e], -1:1, n, n); % 1D finite difference matrix
            I = speye(n);
            % Apply Neumann boundary in 1D
            Dx(1,1) = 1; Dx(1,2) = -1;
            Dx(end,end) = 1; Dx(end,end-1) = -1;
            A = 1/(hx^2)*(kron(I,Dx)+kron(Dx,I));

            %%% Discrete linear differential operator
            Aop = eps*A;

            %%% Function handles
            % Non-linearity
            g = @(t,u) -(u.^3 - u);
            % RHS
            rhs = @(t,u) -Aop*u+(u-(u.^3));
            % Jacobian
            J = @(t,u) -Aop+(speye(size(A))-spdiags(3.*u.^2,0,n2,n2));

        case 'GM_2D'
            %%% Grid
            L = 1; % domain length and width
            hx = L/n; % mesh size
            n2 = n^2; % degrees of freedom

            %%% initial conditions
            rng(0)
            a0x = ones(n,n)*.4 + rand(n,n)*0.2; a = reshape(a0x,[n2,1]);
            h0x = ones(n,n)*.2; h = reshape(h0x,[n2,1]);
            u0 = [a;h];

            %%% Parameters
            Da = .01; Dh = 1;
            rho = 1; mu = 1;
            rho_ = 1; nu = 1;

            %%% Laplacian with periodic boundary conditions
            e = ones(n,1);
            Dx = spdiags([-e 2*e -e], -1:1, n, n); % 1D finite difference matrix
            I = speye(n);
            A = 1/(hx^2)*(kron(I,Dx)+kron(Dx,I));
            % Apply periodic boundary conditions
            for i=1:(n-1)
                A(i*n,i*n+1) = -1/(hx^2);
                A(i*n+1,i*n) = -1/(hx^2);
            end
            A(n2-n+1:n2,1:n) = -1/(hx^2)*sparse(diag(ones(n,1))); A(n2,1) = -1/(hx^2);
            A(1:n,n2-n+1:n2) = -1/(hx^2)*sparse(diag(ones(n,1))); A(1,n2) = -1/(hx^2);

            %%% Discrete linear differential operator
            Aop = blkdiag(Da*A, Dh*A);

            %%% Function handles
            % Block non-linearities
            g = @(t,ah) [rho*ah(1:n2).^2./ah(n2+1:end) - mu*ah(1:n2); rho_*ah(1:n2).^2 - nu*ah(n2+1:end)];
            % Block RHS
            rhs = @(t,ah_vec) [rho*ah_vec(1:n2).^2./ah_vec(n2+1:end) - mu*ah_vec(1:n2) - Da*A*ah_vec(1:n2); rho_*ah_vec(1:n2).^2 - nu*ah_vec(n2+1:end) - Dh*A*ah_vec(n2+1:end)];
            % Block Jacobians
            J = @(t,ah_vec) [2*rho*spdiags(ah_vec(1:n2)./ah_vec(n2+1:end),0,n2,n2) - mu*speye(size(A)) - Da*A, -rho*spdiags(ah_vec(1:n2).^2./ah_vec(n2+1:end).^2,0,n2,n2); 2*rho_*spdiags(ah_vec(1:n2),0,n2,n2), - nu*speye(size(A)) - Dh*A];

        case 'AC_network'
            if strcmp(network,'brightkite')
                try
                    load network_data/loc-Brightkite.mat
                catch
                    error('Data set not found. See READme for instructions.')
                end
            elseif strcmp(network,'minnesota')
                load network_data/minnesota.mat
            elseif strcmp(network,'usroads_subset')
                load network_data/usroads.mat
                %%% select subset of nodes (longitudinal coords between -125 and -115)
                n_full = size(Problem.A,1);

                fileID = fopen('network_data/usroads_coord.txt','r');
                formatSpec = '%f';
                UScoords_ = fscanf(fileID,formatSpec);
                UScoords = reshape(UScoords_,[n_full,2]);
                ind = (UScoords(:,1)>-125) & (UScoords(:,1)<-115);

                Problem.A = Problem.A(ind,ind);
                
                usroads_coords(:,1) = UScoords(ind,1);
                usroads_coords(:,2) = UScoords(ind,2);

            elseif strcmp(network,'ak')
                load network_data/ak2010.mat
            elseif strcmp(network,'luxembourg_osm') 
                load network_data/luxembourg_osm.mat
            elseif strcmp(network,'ny')
                try
                    load network_data/ny2010.mat
                catch
                    error('Data set not found. See READme for instructions.')
                end
            elseif strcmp(network,'roadNet_PA')
                try
                    load network_data/roadNet-PA.mat
                catch
                    error('Data set not found. See READme for instructions.')
                end
            end

            %%% GCC
            G = graph(Problem.A);
            ind = max_connected_component(G);
            Problem.A = Problem.A(ind, ind);

            %%% Grid
            n = size(Problem.A,1);
    
            %%% Initial conditions
            rng(0)
            u0 = (rand(n,1)-0.5)*0.1;

            %%% Interface parameter & diffusion constant
            eps = .05;
            if strcmp(network,'minnesota')
                D_L = 5e03;
            elseif strcmp(network,'usroads_subset')
                D_L = 5e04;
            elseif strcmp(network,'ak')
                D_L = 1e-02;
            elseif strcmp(network,'luxembourg_osm')
                D_L = 5e06;
            elseif strcmp(network,'ny')
                D_L = 5;
            elseif strcmp(network,'roadNet_PA')
                D_L = 1e08;
            else
                D_L = 1e03;
            end

            %%% Graph Laplacian
            adj = Problem.A - spdiags(diag(Problem.A),0,n,n);
            D = spdiags(sum(adj,2),0,n,n);
            L = D - adj;

            %%% Discrete linear differential operator
            Aop = D_L*eps*L;
            
            %%% Function handles
            % Non-linearity
            g = @(t,u) -(u.^3 - u)/eps;
            % RHS
            rhs = @(t,u) -Aop*u+(u-(u.^3))/eps;
            % Jacobian
            J = @(t,u) -Aop+(speye(size(Aop))/eps-spdiags(3.*u.^2,0,n,n))/eps;

        case 'GM_network'
            if strcmp(network,'brightkite')
                try
                    load network_data/loc-Brightkite.mat
                catch
                    error('Data set not found. See READme for instructions.')
                end
            elseif strcmp(network,'minnesota')
                load network_data/minnesota.mat
            elseif strcmp(network,'usroads_subset')
                load network_data/usroads.mat
                %%% select subset of nodes (longitudinal coords between -125 and -115)
                n_full = size(Problem.A,1);

                fileID = fopen('network_data/usroads_coord.txt','r');
                formatSpec = '%f';
                UScoords_ = fscanf(fileID,formatSpec);
                UScoords = reshape(UScoords_,[n_full,2]);
                ind = (UScoords(:,1)>-125) & (UScoords(:,1)<-115);

                Problem.A = Problem.A(ind,ind);
                
                usroads_coords(:,1) = UScoords(ind,1);
                usroads_coords(:,2) = UScoords(ind,2);

            elseif strcmp(network,'ak')
                load network_data/ak2010.mat
            elseif strcmp(network,'luxembourg_osm') 
                load network_data/luxembourg_osm.mat
            elseif strcmp(network,'ny')
                try
                    load network_data/ny2010.mat
                catch
                    error('Data set not found. See READme for instructions.')
                end
            elseif strcmp(network,'roadNet_PA')
                try
                    load network_data/roadNet-PA.mat
                catch
                    error('Data set not found. See READme for instructions.')
                end
            end
            
            %%% GCC
            G = graph(Problem.A);
            ind = max_connected_component(G);
            Problem.A = Problem.A(ind, ind);
            n = size(Problem.A,1);

            %%% initial conditions
            rng(0)
            a0 = ones(n,1)*.4 + rand(n,1)*0.2;
            h0 = ones(n,1)*.2;
            u0 = [a0;h0];

            %%% Parameters
            Da = 10; Dh = 1000;
            rho = 8; mu = 8;
            rho_ = 8; nu = 8;

            %%% Graph Laplacian
            adj = Problem.A - spdiags(diag(Problem.A),0,n,n);
            D = spdiags(sum(adj,2),0,n,n);
            L = D - adj;

            %%% Discrete linear differential operator
            Aop = blkdiag(Da*L, Dh*L);

            %%% Function handles
            % Block non-linearities
            g = @(t,ah) [rho*ah(1:n).^2./ah(n+1:end) - mu*ah(1:n); rho_*ah(1:n).^2 - nu*ah(n+1:end)];
            % Block RHS
            rhs = @(t,ah_vec) [rho*ah_vec(1:n).^2./ah_vec(n+1:end) - mu*ah_vec(1:n) - Da*L*ah_vec(1:n); rho_*ah_vec(1:n).^2 - nu*ah_vec(n+1:end) - Dh*L*ah_vec(n+1:end)];
            % Block Jacobians
            J = @(t,ah_vec) [2*rho*spdiags(ah_vec(1:n)./ah_vec(n+1:end),0,n,n) - mu*speye(size(L)) - Da*L, -rho*spdiags(ah_vec(1:n).^2./ah_vec(n+1:end).^2,0,n,n); 2*rho_*spdiags(ah_vec(1:n),0,n,n), - nu*speye(size(L)) - Dh*L];

    end
    
    % assign arbitrary values for non-required variables to avoid errors
    if ~exist('adj','var')
        adj = 0;
    end
    if ~exist('usroads_coords','var')
        usroads_coords = 0;
    end
end