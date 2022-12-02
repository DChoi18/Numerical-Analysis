function [tout,yout] = ShootingMethod(dy,xint,beta,iBC,ICshoot1,ICshoot2,ystep,guessID,tol)
% function that solves a nonlinear boundary value problem using the 
% shooting method with a RK4 numerical integration scheme and a secant 
% method nonlinear equation iterator
% 
% Inputs: 
%        dy = function handle representing equivalent first order system of
%             ODEs
%        xint     = [a, b] interval over which boundary value problem is to
%                   be solved
%        beta     = boundary condition at x (or eta) = b
%        iBC      = index of the boundary condition to satisfy at b
%        ICshoot1 = initial guess of initial conditions for RK4 integration
%        ICshoot2 = second initial guess of initial conditions for RK4
%                   integration
%        ystep    = step size for RK4 integration (also the same as the
%                   spacing of nodes at which solution is approximated)
%        guessID  = index of the initial condition variable that is guessed
%        tol      = solver tolerance in satisfying the boundary condition at
%                   x = b
% Outputs:
%        tout = points at which solution is approximated 
%        yout = approximate solution at corresponding tout values
%
% Author: Derrick Choi

% Fourth order Runge-Kutta (RK4) integrator set up
a = xint(1) ;
b = xint(2) ;
tspan = [a b];

% first two iterations since secant method needs two initial guesses
[~,yout0] = RK4(dy,ystep,ICshoot1,tspan,1);
[tout,yout] = RK4(dy,ystep,ICshoot2,tspan,1);
 
niter = 2; % iteration counter

% update the guess if needed
while abs((beta-yout(iBC,end))/beta)>tol

    % print to status command line 
    fprintf('Iteration %d, Rel. Error %0.3e\n',niter+1,(beta-yout(iBC,end))/beta)
    
    % update guess with secant method
    temp = ICshoot2;
    ICshoot2(guessID) = ICshoot2(guessID)-((yout(iBC,end)-beta)*(ICshoot2(guessID)-ICshoot1(guessID)))/(yout(iBC,end)-yout0(iBC,end));
    ICshoot1 = temp;
    
    % update other variables (i.e yout0 now becomes yout)
    yout0 = yout;

    % Rerun time integrator to get new yout
    [tout,yout] = RK4(dy,ystep,ICshoot2,tspan,1);

    % update counter for stats purposes
    niter = niter+1;
end

% note once while loop is exited, the output is the yout and tout at niter+1
% (i.e. one more than what is shown in command line)
end
