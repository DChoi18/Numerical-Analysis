function [t,y] = EulerMethod(dy,y0,h,tspan,mode)
% function to numerically approximate solution to first order initial value
% problem with explicit Euler's method given a step size h and the interval of integration
%
% Inputs:
%        dy = function handle for the derivative f(t,y)
%        y0 = initial condition
%        h = step-size
%        tspan = time interval to solve the differential equation
%        mode = 1 or 2, 1 approximates for the entire interval, 2
%               approximates for a particular point (assumes that the point
%               coincides with the end of the tspan variable)
% Outputs:
%        t = times at approximated solution
%        y = value of the approximate solution at corresponding time
% Author: Derrick Choi

% Determine mode
switch mode
    case 1
        % Time vector
        t = tspan(1):h:tspan(2);

        % Initialize approximate solution vector
        y = zeros(length(y0),length(t));
        y(:,1) = y0;

        % Step through the times
        for i = 2:length(t)
            yprime = dy(t(i-1),y(:,i-1));
            y(:,i) = y(:,i-1)+h*yprime;
        end

    case 2
        % Initial time
        t0 = tspan(1);
        while t0<tspan(2)
            % make the last step such that it lands on the end time of the
            % interval
            if t0+h>tspan(2)
                h = tspan(2)-t0;
            end
            % Evaluate derivative
            yprime = dy(t0,y0);
            
            % Update
            y0 = y0+h*yprime;
            t0 = t0+h;
        end
        % output
        t = t0;
        y = y0;
end

end





