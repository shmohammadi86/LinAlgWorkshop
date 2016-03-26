function [ evM ] = evaluateMatrix( epoints,ctrs )
%EVALUATE generates the evaluate matrix in the pde for final answer

global p

m=length(epoints);
n=length(ctrs);
evM=zeros(m,n);

for j=1:n
    evM(:,j)=(abs(epoints(:,1)-ctrs(j,1)).*abs(epoints(:,2)-ctrs(j,2))).^(2*p-1);
end
