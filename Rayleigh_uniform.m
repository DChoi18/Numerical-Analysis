function [x,w] = Rayleigh_uniform(N)
% this applies Rayleigh-Ritz with piecewise linear functions
% defined on an equispaced grid to a boundary value problem
% 
% Inputs: 
%        N = number of Nodes
% Outputs:
%        x = points at which solution to the bvp is approximated
%        w = values of the approximation at the corresponding points
% Note: Code modified from example code posted in class

x = linspace(0,1,N+2); % define an equispaced grid
h = x(2)-x(1);

% evaluate the coefficients at each of the nodes.
[p,q,f] = pqffun(x);

% Element integrals
Q1 = h/12   *(q(2:N)     + q(3:N+1));
Q2 = h/12   *(3*q(2:N+1) + q(1:N));
Q3 = h/12   *(3*q(2:N+1) + q(3:N+2));
Q4 = 1/(2*h)*(p(2:N+2)   + p(1:N+1));
Q5 = h/6    *(2*f(2:N+1) + f(1:N));
Q6 = h/6    *(2*f(2:N+1) + f(3:N+2));


% Assemble the matrix and rhs
A = diag(Q4(1:N)+Q4(2:N+1)+Q2+Q3) + ...
    diag(-Q4(2:N)+Q1,1) + diag(-Q4(2:N)+Q1,-1);

b = Q5'+Q6';


w = A\b; % solve the linear system

w = [0;w;0]; % pad to include boundary data
x = x'; % output a column vector
end


function[p,q,f] = pqffun(x)

% evaluate the functions p(x), q(x), f(x) in the case of
% solving the ode -(x^2*y')'+4y= sin(pi*x); y(0) = y(1) =0

% Problem 1 p,q,f functions
% p = -1*ones(1,length(x));
% q = -4*ones(1,length(x));
% f = 4*x;

% Problem 2 p,q,f functions
p = exp(-x);
q = exp(-x);
f = ((x-1)-(x+1).*exp(-(x-1)));

end


function L = interpolation_matrix(x,pointx)

n  = size(pointx,2);
L2 = ones(n,size(x,2));
for i=1:n
    for j=1:n
        if (i~=j)
            L2(i,:)=L2(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
        end
    end
end
L = L2.';

end
