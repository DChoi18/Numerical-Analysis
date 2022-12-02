function [r,info] = bisect(f,a,b,tol)
% Function to find root of a function to a certain accuracy in a given
% interval using bisection method
%
% Inputs: 
%        fun = function to find root of
%        a =  start of the interval
%        b =  end of interval           
%        tol = accuracy desired for using the bisection method
% Outputs:
%        r = approximate value of the root
%        info = flag (0 - successful, 1 - unsuccessful)
% Author: Derrick Choi

%% Check Input
if f(a)*f(b)>0
    error('No sign change in the given interval')
end
%% Start Bisection Method
c = (a+b)/2; %midpoint

%check if this midpoint is the root or with in the specified tolerance
if f(c) == 0 || abs(a-c)/abs(c)<tol
    r = c;
    info = 0;
    return
end

%bisect interval until root is within specified tolerance 
iterations = 1;
while abs(a-c)/abs(c) > tol
    %choose the interval with the root
    if f(a)*f(c) < 0
        b = c;
    else
        a = c;
    end
    %find new midpoint
    c = (b+a)/2;
    %update iteration counter
    iterations = iterations+1;
end
%% Return Output
fprintf('Number of iterations = %d\n', iterations)
r = c;
info = 0;
end