function [ F ] = setupPDEfunction( x, p, mat, PDESetup)
% setupPDEfunction Create the discretization of the PDE
% Uses the matrices that are created in createMatrices ('mat')
% Input: x: the variables, p: the parameters

%% obtain the differentiation matrices
D1x = mat.Dx;
D2x = mat.D2x;

N = length(D1x);
L = PDESetup.L;

%% Obtain variables and parameters

y = x(1:N);

mu = p(1);
D = p(2);

%% Example 1: Convection-noConvection
%k_T = 1;
%kappa_bar = 10;
%rho_r = -1;
%rho_A = 2;

k_T = 1;
kappa_bar = 10;
rho_r = -1;
rho_A = 2;

kappa = @(rho) kappa_bar/2 * (1 + tanh(rho-rho_r));

gg = @(x,mu) k_T * (rho_A +  mu * (1 + cos(x*pi/2)) );
ff = @(x,y,mu) - k_T * y - kappa(y) .* y + gg(x,mu);


%% Compute the spatially varying terms


XX = linspace(-L/2,L/2,N)';




%% The PDE
% INSTRUCTION
% Use point multiplication and powers, expect when multiplying with the
% differentiation operators/matrices.

f1 = D .* (D2x * y) + ff(XX,y,mu);


%% Construct function
F = [f1];

end

