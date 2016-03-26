function [ z ] = dyrbf( xt,x )
%BRBF generates the matrix of boundary collocation points which satisfies
%the Neumann Condition.
global p


m=length(xt);
n=length(x);

z=zeros(m,n);

for j=1:n
    z(:,j)=sign(xt(:,2)-x(j,2))*(2*p-1).*(abs(xt(:,1)-x(j,1))).^(2*p-1).*(xt(:,2)-x(j,2)).^(2*p-2);
end


end

