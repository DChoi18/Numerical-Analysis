function [xout,yout] = OneDim_SpectralCollocation_DirichletBVP(p,q,f,N,a,b,alpha,beta)
% function to solve the following BVP with Chebyshev spectral collocation
% y'' + p(x)y' + q(x)y = f(x)
% y(a) = alpha, y(b) = beta
%
% Inputs:
%        q = function handle for variable coefficient to y
%        f = function handle for forcing term
%        N = Number of nodes - 1;
%        a = left endpoint
%        b = right endpoint
%        alpha = BC at x = a
%        beta = BC at x = b
% Output:
%        xout = chebyshev collocation points
%        yout = approximation of solution
% Author: Derrick Choi
% Note: cheb function taken from spectral methods in matlab by Trefethen


% get differentiation matrix
[D,x] = cheb(N);

% map to interval of interest
x = (b-a)/2*x+(b+a)/2;
J = diag((b-a)/2*ones(length(x),1)); % jacobian

% First Derivative in interval of interest
D1 = D*inv(J);
% Second derivative matrix
D2 = D*D*inv(J)*inv(J);

% right hand side
rhs = f(x(2:N));

% approximation of left hand side
A = [1 zeros(1,N);
    D2(2:N,1)+D1(2:N,1)*p(x(2)) D2(2:N,2:N)+D1(2:N,2:N)*diag(p(x(2:N)))+diag(q(x(2:N))) D2(2:N,N+1)+D1(2:N,N+1)*p(x(N+1));
    zeros(1,N) 1];
rhs = [beta;rhs;alpha];

% solve linear system
yout = A\rhs;
xout = x;
end