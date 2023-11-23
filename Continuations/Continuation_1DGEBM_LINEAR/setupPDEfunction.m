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

T = x(1:N);

mu = p(1);
D = p(2);
Q0 = p(3);
Q1 = p(4);


%% Fixed parameters
alpha_1 = 0.7;
alpha_2 = 0.289;
T_1 = 260;
T_2 = 293;
K = 0.1;
A = -339.647;
B = 2.218;

%% Compute the spatially varying terms


XX = linspace(-L/2,L/2,N)';
THETA = asin(XX);

%Q = Q0 * (1 + Q1 * cos(2*THETA));
Q = Q0 * (1 - 0.241 * (3*XX.^2-1));
alpha = alpha_1 + (alpha_2-alpha_1) * (1 + tanh(K * (T - (T_1+T_2)/2)))/2;

zD = (1 - XX.^2);
zA = -2*XX;



%% The PDE
% INSTRUCTION
% Use point multiplication and powers, expect when multiplying with the
% differentiation operators/matrices.

f1 = D * zD .* (D2x * T) + D * zA .* (D1x * T) + Q .* (1 - alpha) - (A+B*T) + mu;


%% Construct function
F = [f1];

end

