function [r,msg,x_i,err] = QuasiNewton(x,f,J,tol,Nmax)
% function to approximate root of system of nonlinear equations by a
% Quasi-Newton Method 
% 
% Inputs: 
%         x = vector of initial guesses
%         f = vector of functions to find root of
%         J = inverse of Jacobian of f at the initial guess
%         tol = tolerance of the approximation
%         Nmax = max number of iterations to perform if iteration is divergent
% Outputs:
%         r = vector of roots to f
%         msg = error message of code running (0 - successful, 1 -
%               unsuccessful)
%         x_i = value of the root at the ith iteration
% Author: Derrick Choi

%Store known iterates
x_i(:,1) = x;
invJ = inv(J(x));
%find next iterate
xn_p1 = x-invJ*f(x);
%set iteration counter
iterations = 1;
while norm(xn_p1-x)/norm(xn_p1) > tol %start iterations
    tic
    %stop if more than max allowed iterations
    if iterations > Nmax 
        r = x;
        iterations = Nmax;
        fprintf('Iterations = %d\n, Nmax reached',iterations);
        msg = 1;
        return
    end
    %store iterations
    x_i(:,iterations+1) = xn_p1;
    %New updates
    x = xn_p1;
    xn_p1 = xn_p1-invJ*f(xn_p1);
    iterations = iterations+1;
    toc
end
x_i(:,iterations+1) = xn_p1;
fprintf('Iterations = %d\n', iterations)
r = xn_p1;
msg = 0;
%compute absolute error at each iterations
err = zeros(1,length(x_i));
for i = 1:length(x_i)
    err(i) = norm(x_i(:,i)-r);
end
end
