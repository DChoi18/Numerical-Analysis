function [r,info,errors] = FixedPtSystem(f,g,x0,tol,Nmax,norm_method)
% Function to compute root of a system of equations with initial guess x and
% accuracy given by tol with Newton's Method
%
% Inputs:
%         f,g = functions to find root
%         jacobian = Jacobian matrix
%         x = initial guess (vector input)
%         tol = desired accuracy (based on relative error)
%         Nmax = max number of iterations (in cases that algorithm diverges)
% Outputs:
%         r = approximation of the fixed point
%         info = error message (0 - successful, 1-unsuccessful)
% Author: Derrick Choi

xn_p1 = jacobian(x0(1),x0(2))*[f(x0(1),x0(2));g(x0(1),x0(2))];
iterations = 1;
while norm(abs((xn_p1-x0)./x0),norm_method)>tol
    errors(iterations) = norm(abs((xn_p1-x0)./x0),norm_method);
    %check if iterations exceeded max amount
    if iterations > Nmax
        info = 1;
        iterations = Nmax;
        break
    end
    x0 = xn_p1; %update guess
    xn_p1 = xn_p1-inv(jacobian(xn_p1(1),xn_p1(2)))*[f(xn_p1(1),xn_p1(2));g(xn_p1(1),xn_p1(2))]; %update next point
    iterations = iterations+1; %update number of iterations
end
fprintf('Iterations = %d\n',iterations)
r = xn_p1;
info = 0;
end