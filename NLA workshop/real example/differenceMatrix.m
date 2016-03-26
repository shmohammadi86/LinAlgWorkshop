function [ DM ] = differenceMatrix( intdata,ctrs )
%DIFFERENCEMATRIX get the matrix of RBF of intdata and ctrs
global p
m=length(intdata);
n=length(ctrs);   %number of centers    column
DM=zeros(m,n);
for j=1:n
    DM(:,j)=(abs(intdata(:,1)-ctrs(j,1)).*abs(intdata(:,2)-ctrs(j,2))).^(2*p-1);
end

end

