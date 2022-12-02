function d2fdx2 = CenteredDifference(f,x0,h,order)
% function to approximate the second derivative of the function f at a point x0 using a
% step-size h and a centered difference formula with a certain order method
%
% Inputs: 
%        f = function handle
%        x0 = point to evaluate 
%        h = step size
%        order = order of formula to use for approximation
% Outputs:
%         d2fdx2 = approximation of second derivative
% Author: Derrick Choi

switch order
    case 2
        d2fdx2 = (f(x0+h)-2*f(x0)+f(x0-h))./(h.^2);
    case 4
        d2fdx2 = (-f(x0+2*h)+16*f(x0+h)-30*f(x0)+16*f(x0-h)-f(x0-2*h))./(12*h.^2);
end

end