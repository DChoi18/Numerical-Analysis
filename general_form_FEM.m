function [xh,uapp] = general_form_FEM(k,q,f,N,a,b,alpha,beta)
% this builds a first order approximation to the following 
% BVP
% -d/dx(k(x)du/dx) + q(x) u(x) = f(x) for x \in (a,b)
% u(a) = alpha u(b) = beta
% code adapted from the one in class to handle nonhomogeneous BCs

% N = number of elements = number of nodes -1;
h = (b-a)/(N);
xh = a+[0:N]*h;


A = make_system(k,q,xh,N);
b = make_rhs(f,xh,N);
b = modify_rhs(k,q,b,xh,alpha,beta);

uapp = A\b;


uapp = [alpha; uapp; beta];
xh = xh';
end

function b = make_rhs(f,xh,N)

b = zeros(N-1,1);

for k = 1:N-1
    hm = xh(k+1)-xh(k);
    hp = xh(k+2)-xh(k+1);
    b(k) = 1/hm* integral(@(x) f(x).*(x-xh(k)),xh(k),xh(k+1)) +...
           1/hp* integral(@(x) f(x).*(xh(k+2)-x),xh(k+1),xh(k+2));
end

end

function A = make_system(k,q,xh,N)

A = sparse(N-1);

% make the mass matrix
for j = 1:N-1
    hm = xh(j+1)-xh(j);
    hp = xh(j+2)-xh(j+1);
    A(j,j) = hm^(-2)*integral(k,xh(j),xh(j+1))+hp^(-2)*integral(k,xh(j+1),xh(j+2));
    if (j>1)
      A(j,j-1) = -hm^(-2)*integral(k,xh(j),xh(j+1));
    end
    if (j<N-1)
        A(j,j+1) = -hp^(-2)*integral(k,xh(j+1),xh(j+2));
    end
end

% now make stiffness matrix
for j = 1:N-1
    hm = xh(j+1)-xh(j);
    hp = xh(j+2)-xh(j+1);
    A(j,j) = A(j,j) + hm^(-2)*integral(@(x) q(x).*(x-xh(j)).^2,xh(j),xh(j+1)) +...
                 hp^(-2)*integral(@(x) q(x).*(xh(j+2)-x).^2,xh(j+1),xh(j+2));
    if (j>1)
      A(j,j-1) = A(j,j-1) + hm^(-2)*integral(@(x) q(x).*(x-xh(j)).*(xh(j+1)-x),xh(j),xh(j+1));
    end
    if (j<N-1)
        A(j,j+1) = A(j,j+1) + hp^(-2)*integral(@(x) q(x).*(x-xh(j+1)).*(xh(j+2)-x),xh(j+1),xh(j+2));
    end
end

end

function b = modify_rhs(k,q,v,xh,alpha,beta)
% modify first and last entries for the nonhomogenous Dirichlet BCs
    b = v;
    h = xh(2)-xh(1);
    b(1) = v(1) - alpha*((-1/h^2)*integral(k,xh(1),xh(2))+integral(@(x) q(x).*(xh(2)-x)/(xh(2)-xh(1)).*(x-xh(1))./(xh(2)-xh(1)),xh(1),xh(2)));
    b(end) = v(end)-beta*((-1/h^2)*integral(k,xh(end-1),xh(end))+integral(@(x) q(x).*(xh(end)-x)/(xh(end)-xh(end-1)).*(x-xh(end-1))./(xh(end)-xh(end-1)),xh(end-1),xh(end)));
end