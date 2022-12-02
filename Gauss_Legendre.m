function I = Gauss_Legendre(f,interval,n)
% function to approximate an integral over an interval with
% gaussian-legendre quadrature using n nodes
% 
% Inputs: 
%        f = function handle for the function to integrate
%        interval = [a b] bounds of integration
%        n = number of nodes
% Outputs:
%        I = approximation of the integral
% Interval of integration
a = interval(1);
b = interval(2);
% Determine weights and nodes for Gaussian Quadrature for given n value
[xi,wi] = lgwt(n,a,b);
% Approximate the integral = sum( f(xi)*wi )
I = f(xi)'*wi;
end