function [x_interp,y_interp] = function_eval(f,interval)
% function to evaluate an interpolating function on a fine grid defined by
% the given interval
%
% Inputs: 
%         f = interpolated function function handle
%         interval = interval to evaluate the interpolated function
% Outputs:
%         x_interp = grid
%         y_interp = values of interpolated function on the grid
% Author: Derrick Choi

% Define Grid
x_interp = linspace(interval(1),interval(2),1001);
% Preallocate
y_interp = zeros(1,length(x_interp));
% evaluate
for i = 1:length(x_interp)
    y_interp(i) = f(x_interp(i));
end

end