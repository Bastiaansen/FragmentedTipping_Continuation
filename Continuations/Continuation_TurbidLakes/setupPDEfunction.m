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

%% Example: Turbid Lakes

r = 10;
hv = 0.1;
p = 4;

%cc = @(x) 1 + 0.1*cos(x*pi);
cc = @(x) 1;
pp = @(x) 3 + 0.1 * cos(x*pi);
vv = @(x,y)  1/hv .* 1./(1 + y.^pp(x)) ;
mus = @(x) mu * (1 - 0.1*cos(x*pi));
rr = @(x) r;

%kk = @(x) 1.5 + 0.1 * cos(x*pi);
%vv = @(x,y) 1/hv .* (1 - tanh(kk(x).*(y-1)))/2;


ff = @(x,y,mu) rr(x) .* y .* mus(x) .* 1 ./ (1 + vv(x,y)) - cc(x) .* y.^2;



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

