function A = make_dirichlet(N,a,b)

% function adapted form HW_prob2 in HW 5 solutions to make sparse matrix
% representation of A
h = (b-a)/N;
N = double(int64(N)); %ensure that N is an integer


D = sparse(1:N+1,1:N+1,-2/h^2*ones(N+1,1),N+1,N+1);
E = sparse(2:N+1,1:N,1/h^2*ones(N,1),N+1,N+1);
Ah = E'+D+E;

A = sparse(N+1,N+1);
A(1,1) = 1;
A(2:end-1, :) = Ah(2:end-1,:);
A(end,end) = 1;
return