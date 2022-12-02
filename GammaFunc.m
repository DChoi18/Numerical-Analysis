function [G,nval] = GammaFunc(x,a,b,tol)
% Function to approximate the gamma function with composite trapezoidal
% rule and a truncation of the infinite interval of integration to a finite
% interval
%
% Inputs:
%        x = value to evaluate gamma function at
%        a = lower bound of integration (will be zero based on definition)
%        b = upper bound of integration
%        tol = desired tolerance to evaluate the integral
% Outputs:
%        G = Value of the Gamma function at x
%        nval = number of nodes/function evaluations used for integration
% Author: Derrick Choi

n = b/10; %number of intervals
step = 1;

fun = @(y) CompTrap(@(t) t.^(y-1).*exp(-t),[a b],n);

G = fun(x);

count = 1;
while abs((G-gamma(x))/gamma(x))>=tol
    fprintf('Iteration = %d\n',count)
    if mod(n,10) == 0
        fprintf('n = %d, Relative Error = %0.6f\n',n,abs((G-gamma(x))/gamma(x)))
    end
    n = n+step;
    fun = @(y) CompTrap(@(t) t.^(y-1).*exp(-t),[a b],n);
    G = fun(x);
    count = count +1;
end
nval = n;

end