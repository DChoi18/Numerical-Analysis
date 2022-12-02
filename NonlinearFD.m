function [x,y,iter] = NonlinearFD(y0,Fk,M,xspan,BC,N,tol)
% function to solve a nonlinear dirichlet boundary value problem with
% Finite Differences and Newton's method
%
% Inputs:
%        y0 = initial guess for the Newton method iteration
%        Fk = function handle for the nonlinear equations from discretizing
%             the differential equation
%        M = function handle for Jacobian (without the boundary conditions
%            implemented)
%        N = number of nodes
%        xspan = [a b] (endpoints of the domain)
%        BC = [alpha beta] (values of the boundary conditions)
%        tol = tolerance for Newton's method for root finding

% Nodes and BCs
x = linspace(xspan(1),xspan(2),N+2);
alpha = BC(1);
beta = BC(2);

% First iteration of Newton's method
J = M(x,y0);
J(1,:) = [1 zeros(1,N+1)];
J(end,:) = [zeros(1,N+1) 1];
Delta_y = inv(J)*[alpha;Fk(x,y0);beta];
count = 1; % iteration counter

% update this guess
yn = y0+Delta_y;

% Newton's Method 
while norm(Delta_y)>tol 
    fprintf('Iteration = %d, Error in Newton''s Method = %f \n',count,norm(Delta_y))
    % Get Jacobian Matrix
    J = M(x,yn);
    J(1,:) = [1 zeros(1,N+1)];
    J(end,:) = [zeros(1,N+1) 1];
    
    % Solve linear system for the change in approximation of solution
    Delta_y = inv(J)*-[0;1*Fk(x,yn);0];
    
    % update guess
    yn = yn+Delta_y;
    
    % update iteration counter
    count = count+1;
end

% outputs
x = x';
y = [alpha;yn(2:end-1);beta];
iter = count;
end