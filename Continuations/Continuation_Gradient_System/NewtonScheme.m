function [ xstar, gV, Df] = NewtonScheme( f, x0, p, options )
%NewtonScheme Computes the numerical root xstar of function f

threshold = options.threshold;
delta = options.delta;
maxSteps = options.maxSteps;

x = x0;
ftest = sum(abs(f([x;p])));
steps = 0;
gV = true;
[x,Df] = NewtonStep(f,x,p,delta);

while( ftest > threshold)
   
    [x,Df] = NewtonStep(f,x,p,delta);
    ftest = sum(abs(f([x;p])));
    steps = steps + 1;
    
    if(steps > maxSteps)
       display(['current pBC value: ' num2str(p(1))])
       display('Too many steps needed - No good value found; Continuation step stopped')
       gV = false;
       xstar = x0;
       return
    end
    
end
    
xstar = x;


end

