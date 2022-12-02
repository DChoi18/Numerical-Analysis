function ODE_Dirichlet

% uxx = f(x)
%u(-1) = u(1) = 0
% u(-1) = alpha
% u(1) = beta

alpha = 0 ;
beta = 0;

f = @(x) exp(4*x);
uex = @(x) (exp(4*x)-x*sinh(4)-cosh(4))/16;

NN = 5*2.^[0:3];
%NN = 1000;
err = zeros(1,length(NN));

for j = 1:length(NN)
    
N = NN(j);

[D,x] = cheb(N);
D2 = D*D;

rhs1 = f(x(2:N));

A = [1 zeros(1,N);D2(2:N,:); zeros(1,N) 1];

rhs =[beta; rhs1; alpha];

uapp = A\rhs;

norm(uex(x)-uapp)/norm(uex(x))

err(j) = norm(uex(x)-uapp);

end


semilogy(NN,err,'x-','LineWidth',3)

keyboard
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
