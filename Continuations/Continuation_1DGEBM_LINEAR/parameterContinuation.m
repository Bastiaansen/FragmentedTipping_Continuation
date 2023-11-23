function [ys, ps, unstables, bifs, bps] = parameterContinuation(f, xInit, p, options, tauStart)
%parameterContinuation Performs continuation in the required parameter
%   options.continuationParameter: which is the continuation parameter
%   options.continuationMethod: which continuation method to use.
%   Implemented are: (1) natural, (2) pseudo-arclength, ...

m = options.continuationParameter;
method = options.continuationMethod;
maxSteps = options.parameterMaxSteps;

trackEachNSteps = options.trackParameterEachNSteps;


parMin = options.parameterMinimum;
parMax = options.parameterMaximum;
increment = options.parameterIncrement;

trackStability = options.trackStability;
trackSaddleNode = options.trackSaddleNode;
trackHopf = options.trackHopf;
trackBranchPoints = options.trackBranchPoints;

delta = options.delta;


%% Set-up stuff

unstables = 0; % keeps track of the amount of unstable eigenvalues
bifs = {};
bps = {};
k = 1;
j = 1;
index = 1;

%% PARAMETER CONTINUATION
% Make a new function g, with only x dependence, but such that x(n+1) is
% the bifurcation parameter. Only used internally
n = length(xInit);
g = @(x) f(x(1:n), [x(n+2:n+m);x(n+1);x(n+m+1:end)]);



if(strcmp(method, 'natural'))
    x = xInit;
    p = [p(m); p(1:m-1); p(m+1:end)];
    iter = 1;
    pBC = p(1);
    gV = true;
    ps(index) = p(1);
    ys(:,index) = x;
    while( (pBC > parMin) && (pBC < parMax) && gV && iter < maxSteps)
        xold = x;
        [ x, gV, Dg] = NewtonScheme( g, x, p, options );
        pBC = p(1);
        if(gV)
            if(mod(iter,trackEachNSteps)==0)
                index = index + 1;
                ps(index) = pBC;
                ys(:,index) = x;
            end
            
            if(trackStability)
               unstables(index) = sum(real(eig(Dg))>0);
            end
            if(trackSaddleNode)
                if(~exist('DOLD'))
                    DOLD = sign(det(Dg));
                end
                DNEW = sign(det(Dg));
                if(DNEW*DOLD ~= 1) % then a saddle-node bif has occured 
                    display(['SN - Bifurcation parameter value is ' num2str(pBC)])
                    bifs{j}.type = 'SN';
                    bifs{j}.values = [x;p];
                    j = j+1;
                    
                    index = index + 1;
                    ps(index) = pBC;
                    ys(:,index) = x;
                end
                DOLD = DNEW;
            end
            if(trackHopf)
                if(~exist('HFOLD'))
                    HFOLD = sign(det(bialternateI(Dg)));
                end
                HFNEW = sign(det(bialternateI(Dg)));
                if(HFNEW*HFOLD ~= 1) % then a Hopf bif has occured
                   display(['HF - Bifurcation parameter value is ' num2str(pBC)])
                   bifs{j}.type = 'HF';
                   bifs{j}.values = [x;p];
                   j = j + 1;
                   
                   index = index + 1;
                   ps(index) = pBC;
                   ys(:,index) = x;
                end
                HFOLD = HFNEW;
            end
            if(trackBranchPoints)
               if(~exist('BPOLD'))
                   BPOLD = sign(det([Dg,CompDer( g, [x;p], options.delta*[zeros(length(x),1);1;zeros(length(p)-1,1)]);[(x-xold)',increment]]));
               end
               BPNEW = sign(det([Dg,CompDer( g, [x;p], options.delta*[zeros(length(x),1);1;zeros(length(p)-1,1)]);[(x-xold)',increment]]));
               if(BPNEW * BPOLD ~= 1)
                  display(['BP - Bifurcation parameter value is ' num2str(pBC)])
                  bps{k}.type = 'BP';
                  bps{k}.values = [x;p];
                  k = k + 1;
                  
                  index = index + 1;
                  ps(index) = pBC;
                  ys(:,index) = x;
               end  
               BPOLD = BPNEW;
            end
        end
        p(1) = p(1) + increment;
        
        iter = iter + 1;
        if(mod(iter,trackEachNSteps)==0)
           display(['     Bifurcation parameter value is ' num2str(pBC)]) 
        end
    end
    
elseif(strcmp(method, 'pseudo-arclength'))
    x = [xInit; p(m)];
    p = [p(1:m-1); p(m+1:end)];
    
    iter = 1;
    pBC = x(end);
    gV = true;
    
    ps(index) = x(end);
    ys(:,index) = x;
    % Initial direction computation
    if(sum(abs(tauStart))~=0)
        tau = tauStart;
    else
        tau = [(g([x;p]+delta*[zeros(n,1);1;zeros(length(p),1)])-g([x;p]))/delta;1];
        tau = tau / norm(tau);
    end
    tauPrev = tau;
    gArcLength = @(y) tau' * (y(1:n+1)-x) + increment;
    xold = x;
    x = x - tau * increment;
    
    while( (pBC > parMin) && (pBC < parMax) && gV && iter < maxSteps)

        
        g2 = @(y) [g(y);gArcLength(y)];
        
        [ x, gV, Dg] = NewtonScheme( g2, x, p, options );
                
        pBC = x(end);
        if(gV)
            if(mod(iter,trackEachNSteps)==0)
                index = index + 1;
                ps(index) = pBC;
                ys(:,index) = x;
                if(trackStability)
                    unstables(index) = sum(real(eig(Dg(1:end-1,1:end-1)))>0);
                end
            end
            if(trackSaddleNode)
                if(~exist('DOLD'))
                    DOLD = sign(det(Dg(1:end-1,1:end-1)));
                end
                DNEW = sign(det(Dg(1:end-1,1:end-1)));
                if(DNEW*DOLD ~= 1) % then a saddle-node bif has occured
                    display(['SN - Bifurcation parameter value is ' num2str(pBC)])
                    bifs{j}.type = 'SN';
                    bifs{j}.values = [x; p];
                    j = j+1;
                    
                    index = index + 1;
                    ps(index) = pBC;
                    ys(:,index) = x;
                end
                DOLD = DNEW;
            end
            if(trackHopf)
                if(~exist('HFOLD'))
                    HFOLD = sign(det(bialternateI(Dg(1:end-1,1:end-1))));
                end
                HFNEW = sign(det(bialternateI(Dg(1:end-1,1:end-1))));
                if(HFNEW*HFOLD ~= 1) % then a Hopf bif has occured
                   display(['HF - Bifurcation parameter value is ' num2str(pBC)])
                   bifs{j}.type = 'HF';
                   bifs{j}.values = [x;p];
                   j = j + 1;
                   index = index + 1;
                   ps(index) = pBC;
                   ys(:,index) = x;
                end
                HFOLD = HFNEW;
            end
            if(trackBranchPoints)
               if(~exist('BPOLD'))
                   BPOLD = sign(det([Dg(1:end-1,:);(x-xold)']));
               end
               BPNEW = sign(det([Dg(1:end-1,:);(x-xold)']));
               if(BPNEW * BPOLD ~= 1)
                  display(['BP - Bifurcation parameter value is ' num2str(pBC)])
                  bps{k}.type = 'BP';
                  bps{k}.values = [x;p];
                  tauBPS = (tau+tauPrev)/2;
                  tauBPS = tauBPS / norm(tauBPS);
                  bps{k}.tau = tauBPS;
                  bps{k}.F = g;
                  bps{k}.D0 = [Dg(1:end-1,:);tauBPS'];
                  bps{k}.J = Dg(1:end-1,:);
                  k = k + 1;
                  
                  index = index + 1;
                  ps(index) = pBC;
                  ys(:,index) = x;
               end
               BPOLD = BPNEW;
            end
        end
        
        tau = Dg \ [zeros(length(x)-1,1);1];
        tau = tau / norm(tau);
        tauPrev = tau;
        xold = x;
        x = x - tau * increment;
        
        
        gArcLength = @(y) tau' * (y(1:n+1)-x) + increment;
        
        iter = iter + 1;
        if(mod(iter,trackEachNSteps)==0)
           display(['     Bifurcation parameter value is ' num2str(pBC)]) 
        end
    end
else
   
    error('Parameter continuation method not implemented')
    
end

end

