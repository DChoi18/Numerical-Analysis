function p = HermiteInterp(xi,yi,yi_prime)
% function to interpolate a function using Hermite Polynomial
% Interpolation. Method based on pseudocode from textbook
%
% Inputs: 
%        yi = function value at the nodes
%        yi_prime = derivative value at the nodes
%        xi = interpolation nodes
% Outputs:
%        p = polynomial approximation
% Author: Derrick Choi

%set up lower triangular matrix
n = length(xi);
Q = zeros(2*n, 2*n);

for i = 1:n
   z(2*i-1) = xi(i);
   z(2*i) = xi(i);
   Q(2*i-1,1) = yi(i);
   Q(2*i,1) = yi(i);
   Q(2*i,2) = yi_prime(i);
   if i~=1
       Q(2*i-1,2) = (Q(2*i-1,1)-Q(2*i-2,1))/(z(2*i-1)-z(2*i-2));
   end
end

for i = 3:2*n
    for j = 3:i
        Q(i,j) = (Q(i,j-1)-Q(i-1,j-1))/(z(i)-z(i-j+1));
    end
end

%extract coefficients to hermite polynomial
c = diag(Q);
%form the polynomial
p = @(x) 0;
count = 1;
termcount = 1;
f{1} = @(x) (x-xi(count));
for i = 2:length(c)
    p = @(x) p(x) + c(i)*f{i-1}(x);
    if mod(termcount,2) ==0
        count = count+1;
    end
    termcount = termcount+1;
    f{i} = @(x) f{i-1}(x).*(x-xi(count));
end
%add first term
p = @(x) p(x)+c(1);
end
