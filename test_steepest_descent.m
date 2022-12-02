function [x,ier,x_i] = test_steepest_descent(x0,tol,Nmax)


F = @(x,y,z) [x+cos(x*y*z)-1;
         (1-x)^(1/4) + y+ 0.05*z^2-0.15*z-1;
         -(x)^2-0.1*(y)^2+0.01*y+z-1];
     
J = @(x,y,z) [1-y*z*sin(x*y*z) -x*z*sin(x*y*z) -x*y*sin(x*y*z);
          -1/4*(1-x)^(-3/4) 1 0.1*z-0.15;
          -2*x -0.2*y+0.01 1];
      




[x,iter,ier,x_i] = steepest_descent(F,J,x0,tol,Nmax); 

fprintf('Iterations = %d\n',iter)



return

function [x,iter, ier,x_i ] = steepest_descent(F,J,x0,tol,Nmax)

%tol = desired accuracy
% Nmax = max number of iterations

x_i(:,1) = x0;
for k = 1:Nmax
    g1 = evalg(F,x0);
    z = eval_gradient(x0,F,J);
    z0 = norm(z);
    if z0 == 0
        ier = 1;
        %zero gradient = might be a minimum
        x = x0;
        iter = k;
        return
    end
    z = z/z0; % make unit vector
    alpha1 = 0;
    alpha3 = 1;
    g3 = evalg(F,x0-alpha3*z);
    
    while (g3>= g1)
        alpha3 = alpha3/2;
        g3 = evalg(F,x0-alpha3*z);
        if alpha3<tol/2
            ier = 1; 
            disp('not likely getting a better approximation')
            x = x0;
            iter = k;
            return
        end
    end
    
    alpha2 = alpha3/2;
    g2 = evalg(F,x0-alpha2*z);
    
%     define coefficients for the quadratic
    h1 = (g2-g1)/alpha2;
    h2 = (g3-g2)/(alpha3-alpha2);
    h3 = (h2-h1)/alpha3;
    
%   the minimimum of the quadratric is 
    alpha0 = 0.5*(alpha2-h1/h3); 
    g0 = evalg(F,x0-alpha0*z);
    
%   This is the part that you need to do.
%   find alpha in (alpha0, alpha3) so that g = g(x0-alpha*z) = min{g0,g3}
    ga= min(g0,g3);
    Nsecmax = 100;
    [alpha, ier_sec] = secant(F,alpha0, alpha3,x0,ga,z,tol,Nsecmax);
    
    if ier_sec ==1
        disp('secant did not ')
        x = x0;
        iter = k;
        return
    end


%   Update 
    x0 = x0-alpha*z;
    g = evalg(F,x0);
    x_i(:,k+1) = x0;
    
    if abs(g-g1)<tol 
        ier = 0;
        x = x0;
        iter = k;
        return
    end
        
end

ier = 2;
x = x0;
iter = k;


return

function [alpha, ier] = secant(F,alpha0,alpha3,x0,ga,z,tol,Nmax)

for iter = 1:Nmax
    f3 = evalg(F,x0-alpha3*z)-ga;
    f0 = evalg(F,x0-alpha0*z)-ga;
    alphan = alpha3-f3*(alpha3-alpha0)/(f3-f0);
    
    fn = evalg(F,x0-alphan*z)-ga;
    if abs(fn)==0 || abs(alphan-alpha3)/abs(alphan) <tol
        alpha = alphan;
        ier = 0;
        return
    end
    alpha0 = alpha3;
    alpha3 = alphan;
    
end

% method failed
alpha = alphan;
ier = 1;

return



function g = evalg(F,x0)


n = length(x0);

Ftmp = F(x0(1),x0(2),x0(3));
 
g = 0;
for j = 1:n
    g = g+Ftmp(j)^2;
end

return


function gradval = eval_gradient(x0,F,J)

Ftmp = F(x0(1),x0(2),x0(3));
Jtmp = J(x0(1),x0(2),x0(3));

gradval = 2*Jtmp'*Ftmp;
 

return