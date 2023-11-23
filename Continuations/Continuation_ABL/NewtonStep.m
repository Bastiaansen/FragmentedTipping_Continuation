function [ xNew, Df] = NewtonStep( f, xOld, pars, delta )
%NewtonStep Performs one Newton step
%   xNew = xOld - inv(Df(xOld)) * f(xOld)

Df = zeros(length(xOld));

for i=1:length(xOld)
   del = delta * [zeros(i-1,1);1;zeros(length(xOld)-i,1);zeros(length(pars),1)];
   Df(:,i) = CompDer(f,[xOld;pars],del); 
end

xNew = xOld - Df \ f([xOld; pars]);


end

