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



mu = -10;
%D = 0.01;
%D = 100;
D = 0.30;
Q0 = 341.3;
Q1 = 1/4;

p = [mu; D; Q0; Q1];


%% For a perfect homoclinic pulse
XX = linspace(-L/2,L/2,N)';

    IC_type = 'cold';

    %u = 283;
    u1 = 240;
    u2 = 300;
    %if x < 0.75 && x > 0.25
    
    if strcmp(IC_type, 'warm')
        u = u2 * ones(length(XX),1);
    elseif strcmp(IC_type, 'cold')
        u = u1 * ones(length(XX),1);
    elseif strcmp(IC_type, '2-front')
        u = u2 .* (XX > -0.5) .* (XX < 0.5) + ...
            u1 .* (XX <=-0.5) + u1 .* (XX >= 0.5);
    elseif strcmp(IC_type, '1-front')
        u = u2 .* (XX > -0.5) + ...
            u1 .* (XX <= -0.5);
    elseif strcmp(IC_type, '1-front2')
        u = u2 .* (XX < 0.5) + ...
            u1 .* (XX >= 0.5);
    end


y = [u];


end

