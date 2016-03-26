function [ extra ] = extraterm( xt,x,d )
%Calculate the extra term when RL definition is transformed into Caputo
%form
% x and xj is the data point when we calculate Dx, r is the radial of
% another variable, say y.
global beta
global p
m=length(xt);  %number of collocation points  row
n=length(x);   %number of centers    column 

z=zeros(m,n);
if d==1
    for xx=1:m
        for j=1:n
            if xt(xx,1)==0
                z(xx,j)=0;
            end

            if xt(xx,1)~=0
                z(xx,j)=x(j,1)^(2*p-1)*(abs(xt(xx,2)-x(j,2)))^(2*p-1)/gamma(1-beta)*xt(xx,1)^(-beta)+...
                    (1-2*p)*(x(j,1))^(2*p-2)*(abs(xt(xx,2)-x(j,2)))^(2*p-1)/gamma(2-beta)*xt(xx,1)^(1-beta);
            end
        end
    end
   extra=z; 
end


if d==2
    for xx=1:m
        for j=1:n
            if xt(xx,2)==0
                z(xx,j)=0;
            end

            if xt(xx,2)~=0
                z(xx,j)=x(j,2)^(2*p-1)*(abs(xt(xx,1)-x(j,1)))^(2*p-1)/gamma(1-beta)*xt(xx,2)^(-beta)+...
                   (1-2*p)*(x(j,2))^(2*p-2)*(abs(xt(xx,1)-x(j,1)))^(2*p-1)/gamma(2-beta)*xt(xx,2)^(1-beta); 
            end
        end
    end
   extra=z; 
end

