function [r,info,x_i,err] = NewtonSystem(x0,F,J,tol,Nmax)
% Function to compute root of a system of equations with initial guess x and
% accuracy given by tol with Newton's Method
%
% Inputs:
%         F = functions to find root
%         J = Jacobian matrix
%         x = initial guess (vector input)
%         tol = desired accuracy (based on relative error)
%         Nmax = max number of iterations (in cases that algorithm diverges)
%         norm_method = type of vector norm
% Outputs:
%         r = approximation of the fixed point
%         info = error message (0 - successful, 1-unsuccessful)
%         x_i = approximation of root at ith iteration
% Author: Derrick Choi

x_i(:,1) = x0;
xn_p1 = x0-(J(x0))\F(x0);
iterations = 1;
while norm(xn_p1-x0)/norm(xn_p1)>tol
    tic
    %check if iterations exceeded max amount
    if iterations > Nmax
        info = 1;
        iterations = Nmax;
        r = xn_p1;
        return
    end
    x_i(:,iterations+1) = xn_p1;
    x0 = xn_p1; %update guess
    xn_p1 = xn_p1-J(xn_p1)\F(xn_p1); %update next point
    iterations = iterations+1; %update number of iterations
    toc
end
x_i(:,iterations+1) = xn_p1;
fprintf('Iterations = %d\n',iterations)
r = xn_p1;
info = 0;
%compute absolute error at each iterations
err = zeros(1,length(x_i));
for i = 1:length(x_i)
    err(i) = norm(x_i(:,i)-r);
end
end