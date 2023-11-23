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

%gg = @(x) 0.25 * cos(x*pi) + 0.1*cos(10*x*pi);
%g = @(x) 0.25 * cos(2*x*pi);
g = @(x) cos(x*pi) .* cos(2*x*pi) .* sin(x*pi) ;
%g = @(x) cos(x*pi).*sin(3/2*x*pi);
%g = @(x) 0.5*cos(x*pi);
%g = @(x) 0.25*cos(pi*x);
%ff = @(x,y,mu) 1 .* y - y.^3 + mu + g(x).*y;
%g = @(x) 0.25 * cos(pi*x);
%g = @(x) x;
%g = @(x) 0.5*cos(pi*x).*sin(3/2*pi*x);
g = @(x) 0.5*cos(10*pi*x);

%g = @(x) 0.5*cos(pi*x);
%h = @(x) 0.5*cos(2*pi*x);
%ff = @(x,y,mu) 1 .* y - y.^3 + mu + g(x).*y + h(x);
%ff = @(x,y,mu) g(x) .* y - y.^3 + mu;
ff = @(x,y,mu) y - y.^3 + mu + g(x);

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

