% leinit:  Compute  nodes x (Legendre points) and weights w
%          for Legendre-Gauss-Lobatto quadrature
%          and compute s_{jk}=L_k(x_j)

  function [x,w,s] = leinit(N)

for i=1:N-2
    beta(i)=.5*sqrt(i*(i+2)/((i+.5)*(i+1.5)));
end

  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  x(2:N) = diag(D);

% adding the two end points
  x(1)=-1.;
  x(N+1)=1.;

  [x,i] = sort(x);

%%%%%% compute s_{jk}=L_k(x_j)

for i=1:N+1

    s(i,1)=1;   s(i,2)=x(i);

    for k=1:N-1
        s(i,k+2)=((2*k+1)*x(i)*s(i,k+1)-k*s(i,k))/(k+1);
    end
end


%%%% compute the weight

  w(1)=2/(N*(N+1));
  w(N+1)=w(1);
for j=2:N
  w(j)=w(1)/legen(N,0,x(j))^2;
end
