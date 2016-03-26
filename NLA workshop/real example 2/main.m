clear
clc
BC = 2;     % 2 for NBC

global psi phi psi2 w_t w
N = 1000;   % Time grid number
M = 1000;  % 	Space grid number
[x,w]= legslb(M+1);
[t,w_t]= legslb(N+1);

[ phi,Lambda ] = genbasis( BC , M);    % phi is the space basis in Vm
[ B, psi ,Dpsi, psi2 ] = gen_time_matrix( N );

%% Extend to Time-space case:
[xx,tt] = meshgrid(x,t);
u0 = cos(pi.*xx);
gamma = 0.1;
Exact = (2+tt).*cos(pi*xx);
global NonlinearTerm DiffOperator    % This is f'(u_n)

RHS= cos(pi*xx)+gamma*pi^2.*(tt+2).*cos(pi*xx)+xx.*(tt+2).*cos(pi*xx)-gamma*pi^2.*u0-xx.*u0;
F = psi2*diag(w_t)*RHS*diag(w)*phi';
F=F(:);
NonlinearTerm = xx;    
%% Build Precondition Matrix

Constant = 1000;  % This is the constant for preconditioning. Make sure it is bigger than cond(A)

Ututa = zeros(N,M-1);
Ututa = Ututa(:);

DiffOperator = kron(2.*Lambda,eye(N,N)) + kron(eye(M-1,M-1),gamma.*B);
PrecondMatrix = DiffOperator + Constant.*kron(Lambda,B);

%%%% afun*Ututa = F %%% Here afun is a function of Ututa, because the
%%%% problem is nonlinear, we cannot seperate it into Ax=b, bicgstab allows
%%%% us to solve Ututa as long as we know how the operator acts on Ututa.
Ututa = bicgstab(@afun,F,1e-5,1000,PrecondMatrix);


Ututa = reshape(Ututa,N,M-1);
un = psi'*Ututa*phi + u0;
norm(un-Exact)
mesh(xx,tt,un)
