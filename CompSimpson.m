function I = CompSimpson(f,interval,n)
% function to approximate an integral of function f over an interval with
% n+1 points composite Simpson's Rule
%
% Inputs: 
%        f = function handle for function to integrate
%        interval = [a b] bounds of integration
%        n = number of points (even number)
% Ouputs: 
%        I = approximation of the integral
% Author: Derrick Choi

a = interval(1);
b = interval(2);
h = (b-a)/n;

x = a:h:b;
y = f(x);

I = h/3*( y(1) + 2*sum(y(3:2:end-1))+4*sum(y(2:2:end))+y(end) );
end