function [ fder ] = CompDer( f, x, delta )
%CompDer Numerically computes derivative of f at x with f(x+del)-f(x)/|del|
% f: R^m \rightarrow R^n, x,delta \in R^m

if(length(x) ~= length(delta))
    error('x and delta should be the same length');
end

fder = ( f(x+delta)-f(x) ) / norm(delta);

end

