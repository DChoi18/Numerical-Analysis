function [point,error] = Move2LevelCurve(funs,x0,tol,Nmax)
% function to implement iteration scheme in Problem 4 of Homework 4
% Inputs: 
%         funs = functions and the first derivatives
%         x0 = initial guess
%         tol = tolerance on relative error
%         Nmax = max number of iterations
% Outputs:
%         point = point on level curve
%         error = absolute error from each iteration
% Author: Derrick Choi

%look at what the next point would be
values = funs(x0); 
xn_p1 = zeros(size(x0,1),1);
denom = 0;
for i = 2:length(values)
    denom = denom+values(i)^2; %sum of first derivatives squared
end

for j = 1:length(xn_p1)
    xn_p1(j) = x0(j)-values(1)/denom*values(j+1);
end

%check to see if this is less than tolerance. If not,start iteration scheme
%until satisfied
iterations = 1;
while norm(xn_p1-x0)/norm(x0)>tol
   if iterations > Nmax
       iterations = Nmax;
       break
   end
   error(iterations) = norm(xn_p1-x0)/norm(x0);
   %update guesses
   x0 = xn_p1;
   values = funs(x0);
   xn_p1 = zeros(size(x0,1),1);
   denom = 0;
   for i = 2:length(values)
       denom = denom+values(i)^2;
   end
   
   for j = 1:length(xn_p1)
       xn_p1(j) = x0(j)-values(1)/denom*values(j+1);
   end
   %update iteration counter
   iterations = iterations+1;
end

error(iterations) = norm(xn_p1-x0)/norm(x0);
%assign output
fprintf('Iterations = %d\n',iterations)
point = xn_p1;
end