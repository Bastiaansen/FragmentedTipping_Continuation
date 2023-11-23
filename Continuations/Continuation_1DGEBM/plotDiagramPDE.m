%% Plotting

N = PDESetup.N;
L = PDESetup.L;
h = L/(N-1);

%% L1 norm computation
% Compute the (approximate) L^1 norm of the solution
index = 0;
clear zs
for i=1:N:length(x0)
   index = index + 1;
   %zs(index,:) = ( sum( ((ys(i:i+N-2,:)+ys(i+1:i+N-1,:))/2).^2, 1 )*h ).^(1/2);
   %zs(index,:) = ( sum( ((ys(i:i+N-2,:)+ys(i+1:i+N-1,:))/2), 1 )*h );
   zs(index,:) = mean( ys(i:i+N-1,:), 1);
end

%% Plotting
colorJet = jet(size(zs,1));

% Only plot so many points
plotPoints = options.amountOfPlotPoints;
M = ceil(length(ps)/plotPoints);

if(options.trackStability)
    figure()
    for j=1:size(zs,1)
        for i=1:M:length(ps)-M
           if(ubstables(i) > 0)
              plot([ps(i), ps(i+M)],[zs(j,i),zs(j,i+M)], ':', 'color', colorJet(j,:))
              hold on
           else
              plot([ps(i), ps(i+M)],[zs(j,i),zs(j,i+M)], '-', 'color', colorJet(j,:))
              hold on
           end
        end
    end    
else
    plot(ps(1:M:end),zs(:,1:M:end), 'color', colorJet(1,:))
    hold on
end

xlabel(['$p_' num2str(options.continuationParameter) '$'], 'Interpreter', 'latex')
ylabel('$\|x\|_2^2$', 'Interpreter', 'latex')



%% Plotting of bifurcation points

for j=1:length(bifs)
    if(strcmp(bifs{j}.type, 'SN'))
       ybuf = bifs{j}.values(1:length(x0));
       for i=1:N:length(x0)
            zbuf = ( sum( ((ybuf(i:i+N-2)+ybuf(i+1:i+N-1))/2).^2, 1 )*h ).^(1/2);
            plot(bifs{j}.values(length(x0)+1),zbuf, 'r*', 'MarkerSize', 12, 'LineWidth', 3)
       end
    end
    if(strcmp(bifs{j}.type, 'HF'))
       ybuf = bifs{j}.values(1:length(x0));
       for i=1:N:length(x0)
            zbuf = ( sum( ((ybuf(i:i+N-2)+ybuf(i+1:i+N-1))/2).^2, 1 )*h ).^(1/2);
            plot(bifs{j}.values(length(x0)+1),zbuf, 'bo', 'MarkerSize', 12, 'LineWidth', 3)
       end
    end
end