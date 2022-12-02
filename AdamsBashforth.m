function [tout,yout] = AdamsBashforth(f,h,w0,tspan,m)
% function to approximate solution to an IVP with an m-step Adams-Bashforth
% method with step size h over a interval tspan. Initial steps of required
% for the method are found with fixed time step fourth-order Runge Kutta method
%
% Inputs:
%          f = derivative function handle
%          h = step size
%          wi = initial condition
%          tspan = time interval 
%          m = number of steps
% Outputs:
%          tout = time at which solution is approximated
%          yout = approximation of solution at corresponding times
% Author: Derrick Choi

% Output time vector
tout = (tspan(1):h:tspan(2));
% Output solution vector
yout = zeros(length(w0),length(tout));
switch m
    case 2
        % RK4 to get initial steps
        [~,wi] = RK4(f,h,w0,[tout(1) tout(m)],1);
        
        % Assign initial steps to output
        yout(:,[1:2]) = [w0 wi(2:end)];
        
        % Initial function evaluatons
        currentstep = f(tout(2),yout(:,2));
        step_m1 = f(tout(1),yout(:,1));
        
        % 2-step Adams-Bashforth method
        for i = 3:size(yout,2)
            yout(:,i) = yout(i-1)+h/2*[currentstep step_m1]*[3;-1];
            % update values
            step_m1 = currentstep;
            currentstep = f(tout(i),yout(:,i));

        end
    case 4
        % RK4 to get initial steps
        [~,wi] = RK4(f,h,w0,[tout(1) tout(m)],1);

        % Assign initial steps to output
        yout(:,[1:4]) = [w0 wi(2:end)];
        
        % Initial function evaluations
        currentstep = f(tout(4),yout(:,4));
        step_m1 = f(tout(3),yout(:,3));
        step_m2 = f(tout(2),yout(:,2));
        step_m3 = f(tout(1),yout(:,1));
        
        % 4-step Adams-Bashforth method
        for i = 5:size(yout,2)
            yout(:,i) = yout(:,i-1)+h/24*[currentstep step_m1 step_m2 step_m3]*[55;-59;37;-9];
            % update
            step_m3 = step_m2;
            step_m2 = step_m1;
            step_m1 = currentstep;
            currentstep = f(tout(i),yout(:,i));
        end

end

end