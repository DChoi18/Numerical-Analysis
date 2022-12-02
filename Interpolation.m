function c = Interpolation(xi,yi)
% function to interpolate data (xi,yi) with a polynomial using Vandermonde
% matrix
%
% Inputs: 
%         xi = interpolation nodes
%         yi = data at the nodes
% Outputs:
%         c = vector of coefficients to the monomial basis for polynomials
% Author: Derrick Choi

% preallocate
V = zeros(length(xi),length(xi));
% Build Matrix
V(:,1) = 1;
for n = 2:length(xi)
   V(:,n) = xi.^(n-1); 
end
% Solve linear system
c = V\yi;
end