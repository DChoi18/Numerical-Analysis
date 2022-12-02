function [t,xout] = RK4(dxidt,deltaT,IC,tspan,mode)
% function to solve higher/systems of differential equations with 4th-order
% Runge-Kutta method. (higher order ODEs need to first be converted to a
% system of first order ODEs)
%
% Inputs:
%        dxidt = function handle that outputs a vector (each component is a
%                differential equation)
%        deltaT = time step size
%        IC = vector of initial conditions
%        tspan = time span to solve the set of DEs for
%        mode = 1 or 2, 1 approximates for the whole interval, 2
%               approximates at a point (assumed to be at the end of the time
%               interval)
% Output:
%        t = time at which solution is approximated
%        xout = value of the solution at each corresponding time output
% Author: Derrick Choi

switch mode
    case 1
        % Time vector
        t = tspan(1):deltaT:tspan(2);
        % Preallocate
        xout = zeros(length(IC),length(t));
        % Initialize
        xout(:,1) = IC;

        %Runge-Kutta
        for i = 2:length(t)
            % K-values
            k1 = dxidt(t(i-1),xout(:,i-1));
            k2 = dxidt(t(i-1)+0.5*deltaT,xout(:,i-1)+0.5*deltaT*k1);
            k3 = dxidt(t(i-1)+0.5*deltaT,xout(:,i-1)+0.5*deltaT*k2);
            k4 = dxidt(t(i-1)+deltaT,xout(:,i-1)+deltaT*k3);
            % update
            xout(:,i) = xout(:,i-1)+deltaT/6*(k1+2*k2+2*k3+k4);
        end
    case 2
        t = tspan(1);
        xout = IC;
        while t<tspan(2)
            % make the last step such that it lands on the end time of the
            % interval
            if t+deltaT>tspan(2)
                deltaT = tspan(2)-t;
            end

            % K-values
            k1 = dxidt(t,xout);
            k2 = dxidt(t+0.5*deltaT,xout+0.5*deltaT*k1);
            k3 = dxidt(t+0.5*deltaT,xout+0.5*deltaT*k2);
            k4 = dxidt(t+deltaT,xout+deltaT*k3);
            % update
            xout = xout+deltaT/6*(k1+2*k2+2*k3+k4);            
            t = t+deltaT;
        end
end

end
















