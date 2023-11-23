function [y,p] = setupPDEinitial(PDESetup)
%setupPDEinitial Creates initial condition for the continatuion problem
%   currently should be done manually

%% Obtain the standard & Usually needed stuff
N = PDESetup.N;
L = PDESetup.L;

h = L / (N-1);

%% Here you should put in your initial conditions
% INSTRUCTION: y should be N * #var's and p should be #par's in length
% Make sure the parameters are put in in the right order (as they are
% obtained in the function 'setupPDEfunction.m'



mu = 0.9;
mu = 0.05;
D = 0.01;

%mu = 0;
%D = 0.01;

p = [mu; D];


%% Initial state

XX = linspace(-L/2,L/2,N)';
%y = ones(length(XX),1) * 7;
y = ones(length(XX),1) * 1;


end

