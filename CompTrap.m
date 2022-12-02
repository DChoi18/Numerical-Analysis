function I = CompTrap(f,interval,n)
% function to approximate an integral of function f over an interval with
% n+1 points composite trapezoidal Rule
%
% Inputs: 
%        f = function handle for function to integrate
%        interval = [a b] bounds of integration
%        n = number of points
% Ouputs: 
%        I = approximation of the integral
% Author: Derrick Choi

% Determine beginning and end points
a = interval(1);
b = interval(2);

% point spacing
h = (b-a)/n;

% generate all the points used for quadrature
x = a:h:b;
% function value at these points
y = f(x);
% Use trapezoidal rule  = h/2*( f(a)+f(b)+2*sum(f(interior points)) )
I = 1/2*[1,1]*[y(1:end-1);y(2:end)]*h*ones(length(x)-1,1);
end