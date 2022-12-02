function [app,ier] = adaptive_quad(f,a,b,quadrature,tol,Nmax)
% function to approximate an integral from a to b with adaptive quadrature
% to a desired tolerance with N levels of refinement
%
% Inputs: 
%        f = function handle for function to integrate
%        a = starting value for integration
%        b = end value value for integration
%        quadrature = type of quadrature scheme to use for the adaptive quadrature
%        tol = desired tolerance
%        Nmax = max levels of refinement
% Outputs:
%        app = approximation of the integral
%        ier = error message
%
% Code modified from example code provided in class 

% number of nodes - 1 in the composite quadrature
n = 4;
% Determine the type of quadrature scheme to use
switch quadrature
    case 'Composite Trapezoid'
        q = @CompTrap;
    case 'Composite Simpson'
        q = @CompSimpson;
    case 'Gaussian'
        q = @Gauss_Legendre;
end
% Intlist  = interval list to be processed
% Intsave = integral approximations for the corresponding intervals
% Levels = number of level of the interval

% initialize everything - these are for processing
Intlist(:,1) = [a,b];
Intsave(1) = q(f,[a b],n);

% create arrays for storing adaptive mesh
Intkeep = [];
Intvalkeep = [];

app = 0; % our approximation using the adaptive scheme

i = 1;  % counter for interval that need to be processed.
Levels(1) = 0;

while i>0
    %     split interval in half
    aloc = Intlist(1,1); bloc = Intlist(2,1); %reads endpoints
    xmid = (aloc+bloc)/2; % computes midpoint
    Q1 = q(f,[aloc,xmid],n);
    Q2 = q(f,[xmid,bloc],n);
    if (abs(Intsave(1)-(Q1+Q2))<tol)
        app = app + (Q1+Q2);
        %         Store the adaptive mesh
        Intkeep = [Intkeep; [aloc xmid; xmid bloc]];
        Intvalkeep = [Intvalkeep; Q1; Q2];
        %         Delete information
        Intlist(:,1) = [];
        Intsave(1) = [];
        Levels(1) = [];
        i = i -1; % reduce interval count by 1
    else
        Nlevel = Levels(1) + 1;
        if Nlevel >Nmax
            disp('error: did not converge')
            ier = 1;
            return
        end
        %         add the sub intervals to the list to be processed.
        Intlist= [Intlist [aloc;xmid]];
        Intsave = [Intsave; Q1];
        Levels = [Levels; Levels(1)+1];
        Intlist= [Intlist [xmid;bloc]];
        Intsave = [Intsave; Q2];
        Levels = [Levels; Levels(1)+1];
        i = i + 1;

        %         Remove the information for the interval we just checked.
        Intlist(:,1) = [];
        Intsave(1) = [];
        Levels(1) = [];
    end
end



% plot adaptive mesh
xvals = Intkeep(:,1);
xvals = unique(xvals);
xvals = [xvals; b];

fx = f(xvals);
xx = linspace(a,b,1000);
figure
plot(xx,f(xx),'-',xvals,fx,'ro')


ier = 0;
disp('Sucess: adaptive quad')
fprintf('number of intervals is  %d .\n',size(Intkeep,1));

end