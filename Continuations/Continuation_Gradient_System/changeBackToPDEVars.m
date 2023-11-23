%% changeBackToPDEVars changes the ys back to the original vars
% changes back to the original u-kind of variables mostly used in PDEs

N = PDESetup.N;
L = PDESetup.L;

us = {};

index = 0;
for i=1:N:length(x0)
	index = index + 1;
	us{index} = ys(i:i+N-1,:);   
end

xs = linspace(-L/2,L/2,N);

%% Plot configuration at the end

colourJet = jet(length(us));

figure()
hold on

for i=1:length(us)
    plot(xs,us{i}(:,end), 'color', colourJet(i,:))
end

xlabel('$x$', 'Interpreter', 'latex')
ylabel('$T(x)$', 'Interpreter', 'latex')