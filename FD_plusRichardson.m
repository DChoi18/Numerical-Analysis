function FD_plusRichardson

% y'' = p(x)y'+q(x)y+r(x)
% y(a) = alpha, y(b) = beta


a = 0; b = 1;
alpha = 0; beta = 2;


%Data for Richardson is made with h = 0.1, 0.05, 0.025
h= 0.1;

N = (b-a)/h;

x = [a:h:b];


p = @(x) 0*x;
q = @(x) 4*ones(length(x),1);
r = @(x) -4*x;
Matypp = 1/h^2*(diag(2*ones(N-1,1))-diag(ones(N-2,1),-1)-diag(ones(N-2,1),1));
Matyp = 1/(2*h)*(-diag(ones(N-2,1),-1)+diag(ones(N-2,1),1));

A = Matypp+diag(p(x(2:end-1)))*Matyp+ diag(q(x(2:end-1)));
    
rhs = -r(x(2:end-1));
rhs(1) = rhs(1) + (1/h^2-1/(2*h)*-p(x(2)))*alpha;
rhs(end) = rhs(end) +(1/h^2+1/(2*h)*-p(x(end-1)))*beta;

sol = A\rhs';

yapp =[alpha; sol; beta];

c2 = 1/70*(8-12*sin(log(2))-4*cos(log(2)));
c1 = 11/10-c2;

y = @(x) c1*x+c2./(x.^2)-3/10*sin(log(x))-1/10*cos(log(x));

yex = y(x);

figure(1)
plot(x,yapp,x,yex)
figure(2)
plot(x,abs(yapp-yex'),'o-')
abs(yapp-yex')
keyboard



%  for creating the Richardson extrapolation examples
% picking x = 1.5;
n = 3;
Nlevel = 3;
Data = zeros(n,Nlevel);
Data(1,1) = 1.481120262112219;
Data(2,1) = 1.481149589575161;
Data(3,1) = 1.481156957695558;

Data = Richardson(Data,n,Nlevel);

keyboard

return



function Data = Richardson(Data,n,Nlevel)

for k = 2:Nlevel
    for j = k:n
        Data(j,k) = (2^(k-1)*Data(j,k-1)-Data(j-1,k-1))/(2^(k-1)-1);
    end
end

% For fun try to make this faster.  Note the inner loop can 
% be avoided.

return
