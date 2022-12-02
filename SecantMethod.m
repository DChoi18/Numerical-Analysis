function [r,info,errors,iterates] = SecantMethod(f,x0,x1,tol,Nmax)
% Function to compute root of a function with initial guesses x0,x1
% accuracy given by tol with Secant Method
% 
% Inputs:
%         f = function to find root
%         x0 = initial guess
%         x1 = next guess
%         tol = desired accuracy (based on relative error)
%         Nmax = max number of iterations (in cases that algorithm diverges)
% Outputs: 
%         r = approximation of the root
%         info = error message (0 - successful, 1-unsuccessful)
%         errors = relative error associated with each iteration
%         iterates = approximation of root at each iteration
% Author: Derrick Choi

iterates(1) = x0;
iterates(2) = x1;

xn_p1 = x1-f(x1)*(x1-x0)/(f(x1)-f(x0)); %next point
iterations = 1;
while abs((xn_p1-x1)/x1)>tol
    %check if iterations exceeded max amount
    if iterations > Nmax
        info = 1;
        break
    end
    iterates(iterations+2) = xn_p1;
    errors(iterations) = abs((xn_p1-x1)/x1);
       x0 = x1;%update guess
    x1 = xn_p1; %update guess
    xn_p1 = x1-f(x1)*(x1-x0)/(f(x1)-f(x0)); %update next point
    iterations = iterations+1; %update number of iterations
end

iterates(iterations+2) = xn_p1;
errors(iterations) = abs((xn_p1-x1)/x1);
fprintf('Iterations = %d\n',iterations)
r = xn_p1;
info = 0;
end