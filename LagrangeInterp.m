function p = LagrangeInterp(nodes,data)
% function to interpolate a function using Lagrange Polynomial
% Interpolation (barycentric form)
%
% Inputs: 
%        data = function value at the nodes
%        nodes = interpolation nodes
% Outputs:
%        p = polynomial approximation
% Author: Derrick Choi

% Product term
phi = @(x) prod(x-nodes);
% weights
w = zeros(1,length(nodes));
for i = 1:length(w)
    w(i) = 1/(prod(nodes(i)-nodes(nodes~=nodes(i))));
end

p = @(x) phi(x).*sum(w./(x-nodes).*data);

end