function [x,y] = FD_DirichletBVP(p,q,r,xspan,BC,h)
% function uses a second order Finite Difference scheme to solve the 
% Dirichlet boundary value problem of the form
% 
% y" = p(x)y'+q(x)y+r(x) for x in (a,b)
% y(a) = alpha
% y(b) = beta
%
% Inputs: 
%        p,q,r = function handles to define the differential equation
%        xspan = endpoints of domain
%        BC = dirichlet boundary condition values
%        h = node/mesh spacing
% Output:
%        x = points at which the solution to the BVP is approximated
%        y = value of the approximation of the solution at corresponding x
% Author: Derrick Choi

%% Build Left-Hand Side

% Extract BCs
alpha = BC(1);
beta = BC(2);
% output nodes
x = xspan(1):h:xspan(2);
N = length(x);

% contribution from second derivative
A = 1/(h^2)*(-2*diag(ones(N-2,1))+diag(ones(N-3,1),1)+diag(ones(N-3,1),-1));

% Contribution from p(x)
A = A+1/h*diag(p(x(2:end-1)))*(diag(0.5*ones(N-3,1),1)+diag(-0.5*ones(N-3,1),-1));
A = A';
% Contribution from q(x)
A = A-diag(q(x(2:end-1)));

%% Right Hand Side

f = transpose(r(x(2:end-1)));

% add BC contributions
f(1) = f(1)-(alpha/(h^2)+alpha/(2*h)*p(x(2)));
f(end) = f(end)-(beta/(h^2)-beta/(2*h)*p(x(end-1)));

%% Solve linear system

y = A\f;

%% Add BC to solution vector

y = [alpha;y;beta];
x = x';
end