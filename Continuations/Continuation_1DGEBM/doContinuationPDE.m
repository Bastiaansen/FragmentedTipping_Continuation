%% doContinuation
% This script runs a parameter continuation

%% OPTIONAL: TURN OFF WARNINGS
warning('off','all')
warning

clear all
close all
%% Set-up

PDESetup.N = 251; % amount of discretization points (ideally odd)
PDESetup.L = 2; % The length of the domain under consideration

% 1 First: create the necessary differentiation matrices
createMatrices();


% 2 Second: create the function f that we need
f = @(x,p) setupPDEfunction(x,p,matrices, PDESetup);

% 3 Third: Put in a good enough approximation of a starting point
[x0,p] = setupPDEinitial(PDESetup);


%% Options

options.continuationParameter = 1; % which par is the bifurcation parameter
options.parameterMinimum = -20;
options.parameterMaximum = 200;
options.parameterIncrement = -0.01;
options.parameterMaxSteps = 10^5; % Maximum number of par. cont. steps

options.trackParameterEachNSteps = 200; % Save how many steps? (als outputted)

options.trackStability = 0; % Keep track of stability of fixed points? Time-intensive for larger problem as eigenvalues need to be computed
options.trackSaddleNode = 0; % Keep track of saddle-node bif's?
options.trackHopf = 0; % Keep track of Hopf bifurcations?
options.trackBranchPoints = 0; % Keeps track of branch points

options.searchBranches = 0; % Search branches as well?

options.amountOfPlotPoints = 100; % How many plot points per branch (to prevent long figure creation and lack of visual clues for stability)

options.threshold = 10e-5; % How precise do we need our root |f(x)| < threshold
options.delta = 0.001; % Used for numerical derivatives
options.maxSteps = 10^2; % How many steps are allowed in Newton's method?

% Which continatuion method to use. Implemented are 'natural',
% 'pseudo-arclength', ...
options.continuationMethod = 'pseudo-arclength';

%% Initialization
% First:  find a good root with f(x,p) = 0;


% First write a function g that only has one input for internal use
g = @(x) f(x(1:length(x0)),x(length(x0)+1:end));

[xInit, ~] = NewtonScheme( g, x0, p, options );

%% Parameter continuation

[ys, ps, ubstables, bifs, bps] = parameterContinuation(f, xInit, p, options,0);

%% Search other branches as well

if(options.searchBranches)
    
    
    
    yBranches1 = {};
    pBranches1 = {};
    unstBranches1 = {};
    bifsBranches1 = {};
    bpsBranches1 = {};
    
    yBranches2 = {};
    pBranches2 = {};
    unstBranches2 = {};
    bifsBranches2 = {};
    bpsBranches2 = {};
    
    for k=1:length(bps)
        display(['Now Searching from branch point at ' num2str(bps{k}.values(length(x0)+1)) ' > direction 1'])
        [ tauNew ] = branchSwitch( bps{k}.F, bps{k}.values, bps{k}.tau, bps{k}.J, bps{k}.D0, options);
        xInitNew = bps{k}.values(1:length(x0));
        pBC = bps{k}.values(length(x0)+1);
        pRest = bps{k}.values(length(x0)+2:end);
        m = options.continuationParameter;
        pNew = [pRest(1:m-1);pBC;pRest(m:end)];
        
        [yBranches1{k}, pBranches1{k}, unstBranches1{k}, bifsBranches1{k}, bpsBranches1{k}] = parameterContinuation(f, xInitNew, pNew, options, tauNew);
        %options.parameterIncrement = - options.parameterIncrement; % so we effectively go other dir
        display(['Now Searching from branch point at ' num2str(bps{k}.values(length(x0)+1)) ' > direction 2'])
        [yBranches2{k}, pBranches2{k}, unstBranches2{k}, bifsBranches2{k}, bpsBranches2{k}] = parameterContinuation(f, xInitNew, pNew, options, -tauNew);
        %options.parameterIncrement = - options.parameterIncrement;
    end
end

%% Plotting

plotDiagramPDE();

%% Back to original variables

changeBackToPDEVars();