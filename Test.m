% number of intervals
n = 8;
%reaction coeffecient
sigma = 1;
%interval size
h = 1/n;
%make second derivative Dirichlet matrix
DIFF = 2*diag(sparse(ones(n-1,1))) ...
-1*diag(sparse(ones(n-2,1)),-1) ...
-1*diag(sparse(ones(n-2,1)), 1);
%change upper-left entry for left Neumann b.c.
DIFF(1,1) = 1;
%change lower-right entry for right Neumann b.c.
DIFF(n-1,n-1) = 1;
%build the dicrete system
A=(1/h^2)*DIFF + sigma*diag(sparse(ones(n-1,1)));
%display matrix
full(A)
%%
clc; clear; close all
npoints = 10;
xleft = 0.0;
xright = 10.0;
deltax = (xright-xleft)/(npoints-1);
x = zeros(npoints,1);
xhalf = zeros(npoints-1,1);
A = zeros(npoints);
b = zeros(npoints,1);
for k = 1:npoints
  x(k) = xleft + (k-1)*deltax;
end
for k = 1:npoints-1
  xhalf(k) = (x(k)+x(k+1))/2.0;
end
A(1,1) = -3/(xhalf(1)^3-x(1)^3)*xhalf(1)^2/(x(2)-x(1));
A(1,2) = +3/(xhalf(1)^3-x(1)^3)*xhalf(1)^2/(x(2)-x(1));
for k = 2:npoints-1
  A(k,k-1) = +3/(xhalf(k)^3-xhalf(k-1)^3)*xhalf(k-1)^2/(x(k)-x(k-1));
  A(k,k)   = -3/(xhalf(k)^3-xhalf(k-1)^3)*(xhalf(k)^2/(x(k+1)-x(k))+xhalf(k-1)^2/(x(k)-x(k-1)));
  A(k,k+1) = +3/(xhalf(k)^3-xhalf(k-1)^3)*xhalf(k)^2/(x(k+1)-x(k));
end
A(npoints,npoints) = 1.0;
for k = 1:npoints-1
  b(k) = -4.0*pi*exp(-x(k)^2);
end
b(npoints) = 0.0;
y = A\b;
plot(x,y)