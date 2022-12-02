function [r,msg,x_i,err] = BroydenMethod(x,F,J,tol,Nmax)
% function to approximate root of system of nonlinear equations by a
% Quasi-Newton Method (code based on pseudo-code of the textbook)
% 
% Inputs: 
%         x = vector of initial guesses
%         f = vector of functions to find root of
%         J = Jacobian of f at the initial guess
%         tol = tolerance of the approximation
%         Nmax = max number of iterations to perform if iteration is divergent
% Outputs:
%         r = vector of roots to f
%         msg = error message of code running (0 - successful, 1 -
%               unsuccessful)
%         x_i = value of the root at the ith iteration
% Author: Derrick Choi


v = F(x);
A = inv(J(x));
s = -A*v;
%Store known iterates
x_i(:,1) = x;
%next point
xn_p1 = x+s;
%set iteration counter
iterations = 1;
while norm(xn_p1-x)/norm(xn_p1) > tol %start iterations
    tic
    %stop if more than max allowed iterations
    if iterations > Nmax 
        r = xn_p1;
        iterations = Nmax;
        fprintf('Iterations = %d\n, Nmax reached',iterations);
        msg = 1;
        return
    end
    %store iterations
    x_i(:,iterations+1) = xn_p1;
    
    %Update inverse of jacobian
    w = v;
    v = F(xn_p1);
    y = v-w;
    z = -A*y;
    p = -s'*z;
    
    u = s'*A;
    u = u';
    
    A = A+1/p*(s+z)*u';
    s = -A*v;
    % update points
    
    x = xn_p1;
    xn_p1 = xn_p1+s;
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