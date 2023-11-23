%% createMatrices. Creates differentiation matrices
% This script creates the matrices that represent the differential
% operators that are needed to numerically represent a PDE.

%% Obtain correct values
N = PDESetup.N; % Amount of grid points
L = PDESetup.L; % Length of the domain

h = L / (N-1); % The length between two grid points

%% Identy matrices
ex = sparse(1:N,1:N,1,N,N);

%% First derivative
Dx = sparse(1:N-1,2:N,1/2,N,N);
Dx = (Dx - Dx');
%Dx(1,N) = -1/2; Dx(N,1) = +1/2;% Periodic boundary conditions
Dx(1,2) = 0; Dx(N,N-1) = 0;


Dx = Dx/h;

%% Second derivative
D2x = sparse(1:N-1,2:N,1,N,N) - sparse(1:N,1:N,1,N,N);
D2x = (D2x + D2x');
%D2x(1,N) = 1; D2x(N,1) = 1; % Periodic boundary conditions
D2x(1,2) = 2; D2x(N,N-1) = 2; % Neumann Boundary Conditions



D2x = D2x/h^2;

%% Put everything in variable matrices
matrices.Dx = Dx;
matrices.D2x = D2x;