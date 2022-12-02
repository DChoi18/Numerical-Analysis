function [r,info,errors,iterates] = NewtonMethod(f,fprime,x0,tol,Nmax)
% Function to compute root of a function with initial guess x0 and
% accuracy given by tol with Newton's Method
% 
% Inputs:
%         f = function to find root
%         fprime = first derivative of function
%         x = initial guess
%         tol = desired accuracy (based on relative error)
%         Nmax = max number of iterations (in cases that algorithm diverges)
% Outputs: 
%         r = approximation of the root
%         info = error message (0 - successful, 1-unsuccessful)
%         errors = relative error associated with each iteration
%         iterates = approximation of root at each iteration
% Author: Derrick Choi

iterates(1) = x0;
xn_p1 = x0-f(x0)/fprime(x0); %next point
iterations = 1;
while abs((xn_p1-x0)/x0)>tol
    %check if iterations exceeded max amount
    if iterations > Nmax
        info = 1;
        break
    end
    iterates(iterations+1) = xn_p1;
    errors(iterations) = abs((xn_p1-x0)/x0);
    x0 = xn_p1; %update guess
    xn_p1 = xn_p1-f(xn_p1)/fprime(xn_p1); %update next point
    iterations = iterations+1; %update number of iterations
end

iterates(iterations) = xn_p1;
errors(iterations) = abs((xn_p1-x0)/x0);
fprintf('Iterations = %d\n',iterations)
r = xn_p1;
info = 0;
end