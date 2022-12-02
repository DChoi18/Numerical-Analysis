
function derive_test1
% derive_test1.m - Chebyshev differentation of a smooth function
xx = [-1:.01:1]; uu = exp(xx).*sin(5*xx); clf
for N = [10 20]
[D,x] = cheb(N);
u = exp(x).*sin(5*x);
subplot('position',[ .15 .66-.4*(N==20) .31 .28])
plot(x,u,'.','markersize',14), grid on
line(xx,uu,'linewidth',.8)
title([ 'u(x), N=' int2str(N)])
error = D*u - exp(x).*(sin(5*x)+5*cos(5*x));
subplot('position', [.55 .66-.4*(N==20) .31 .28])
plot(x,error,'.','markersize',14), grid on
line(x,error,'linewidth',.8)
title( [' error in u''(x), N=' int2str(N)])
pause
keyboard
end

% test second derivative

figure(2)
for N = [10 20]
[D,x] = cheb(N);
u = exp(x).*sin(5*x);
upp = exp(x).*(5*cos(5*x)-25*sin(5*x))+exp(x).*(sin(5*x)+5*cos(5*x));
subplot('position',[ .15 .66-.4*(N==20) .31 .28])
plot(x,u,'.','markersize',14), grid on
line(xx,uu,'linewidth',.8)
title([ 'u(x), N=' int2str(N)])
error = D^2*u - upp;
subplot('position', [.55 .66-.4*(N==20) .31 .28])
plot(x,error,'.','markersize',14), grid on
line(x,error,'linewidth',.8)
title( [' error in u''(x), N=' int2str(N)])
pause
end


% from Spectral methods in Matlab
return

% CHEB
% compute D = differentiation matrix, x = Chebyshev grid
function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D - diag(sum(D'));
% off-diagonal entries
% diagonal entries

return