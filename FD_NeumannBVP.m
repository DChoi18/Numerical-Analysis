function [x,y] = FD_NeumannBVP(p,q,r,xspan,BC,h,method)
% function uses a second order Finite Difference scheme to solve the
% 1D - Neumann boundary value problem of the form
%
% y" = p(x)y'+q(x)y+r(x) for x in (a,b)
% y'(a) = alpha
% y'(b) = beta
%
% Inputs:
%        p,q,r = function handles to define the differential equation
%        xspan = endpoints of domain
%        BC = Neumann boundary condition values
%        h = node/mesh spacing
%        method = finite difference scheme for approximating derivative
%                 Boundary Condition
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

switch method
    case 'Forward'
        % contribution from second derivative
        A = 1/(h^2)*(-2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1));

        % Contribution from p(x)
        A = A+1/h*diag(p(x))*(diag(-0.5*ones(N-1,1),1)+diag(0.5*ones(N-1,1),-1));

        % Contribution from q(x)
        A = A-diag(q(x));
        
        % Right Hand Side

        f = transpose(r(x(2:end-1)));

        % add BC contributions
        A(1,1:2) = [-1/h 1/h];
        A(end,end-1:end) = [-1/h 1/h];

        % Solve
        f = [alpha;f;beta];
        y = A\f;
        x = x';
    case 'Centered'
        % contribution from second derivative
        A = 1/(h^2)*(-2*diag(ones(N+2,1))+diag(ones(N+1,1),1)+diag(ones(N+1,1),-1));

        % Contribution from p(x)
        A = A+1/h*diag(p(x))*(diag(-0.5*ones(N+1,1),1)+diag(0.5*ones(N+1,1),-1));

        % Contribution from q(x)
        A(2:N+1,2:N+1) = A(2:N+1,2:N+1)-diag(q(x));
        
        % RHS
        b = transpose(r(x));
        
        % Add BC contributions
        A(1,1:3) = [-1/(2*h) 0 1/(2*h)];
        A(end,end-2:end) = [-1/(2*h) 0 1/(2*h)];
        
        % Solve
        b = [alpha;b;beta];
        y = A\b;
        % Output
        y = y(2:end-1);
        x = x';
end