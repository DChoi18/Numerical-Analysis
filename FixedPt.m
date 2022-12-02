function [x,er] = FixedPt(fun,x0,tol,Nmax)
% Function to compute fixed point of a function with initial guess x and
% accuracy given by tol
% 
% Inputs:
%         fun = function to find fixed point
%         x = initial guess
%         tol = desired accuracy (based on relative error)
%         Nmax = max number of iterations (in cases that algorithm diverges)
% Outputs: 
%         x = approximation of the fixed point
%         er = error message (0 - successful, 1-unsuccessful)
% Author: Derrick Choi

iterations = 1;
xn_p1 = fun(x0);
while abs((xn_p1-x0)/xn_p1)>tol %relative error stopping criteria
    if iterations > Nmax %stop if past the Nmax iterations
        er = 1;
        break
    end
    x0 = xn_p1; %update the guess
    xn_p1 = fun(xn_p1);
    iterations = iterations+1; %increase iteration counter
end
fprintf('Iterations = %d\n',iterations-1)
x = xn_p1;
er = 0;
end