function [tout,yout] = TaylorSecondOrder(f,DfDt,y0,h,tspan)
% function to approximate solution to an IVP with second order Taylor's
% method with step size h and initial condition y0 over a span of time
%
% Inputs:
%       f = derivative (function handle)
%       DfDT = total derivative of f (function handle)
%       y0 = initial condition
%       h = step-size
%       tspan = [a b] time interval for the IVP
% Outputs:
%       tout = times at which solution is approximated
%       yout = approximation of the solution at corresponding times
% Author: Derrick Choi

% Build time vector
tout = tspan(1):h:tspan(2);

% Preallocate
yout = zeros(1,length(tout));

% Initialize
yout(1) = y0;

% Second-order Taylor Method
for i = 2:length(yout)
    yout(i) = yout(i-1)+h*f(tout(i-1),yout(i-1))+h^2/2*DfDt(tout(i-1),yout(i-1));
end

end