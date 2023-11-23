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

%% Example: Atmospheric Boundary Layer


qq = @(x) 1 - 0.1 .* x;
ll = @(x) 1 + 0 .* x;
alpha = @(x) 3 + 0 .* x;
gg = @(x,y) exp(-2 * alpha(x) .* y);

cc = @(x,mu) 0.1 * x + mu;

ff = @(x,y,mu) qq(x) - ll(x) .* y - cc(x,mu) .* y .* gg(x,y);


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

