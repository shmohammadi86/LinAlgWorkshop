function [vec,value]=powernom(start,A,n)
%Power method for computing eigenvalues
%start - starting vector
%A- matrix to compute eigenvalue 
%n-num of iteration desired
x=start;
for i=1:n
y=A*x;%
n=norm(y);%
x=y/n;%
%pause
end
vec=x;
value=n;