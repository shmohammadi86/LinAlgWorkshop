function z=approximate_functiony(xt,x)  

global beta
global p
A=zeros(1,2*p-2);
for k=0:(2*p-3)
    A(k+1)=nchoosek(2*p-3,k);
end

m=length(xt);  %number of collocation points  row
n=length(x);   %number of centers    column 



z=zeros(m,n);

for j=1:n
        tmp=0;
        for k=0:(2*p-3)
            tmp=tmp+ A(k+1)/(k-beta+2) * (-1)^(2*p-3-k).* (xt(:,2)-x(j,2)).^(2*p-3-k).*( ...
                                    xt(:,2).^(k-beta+2) -2*(max(0,xt(:,2)-x(j,2))).^(k-beta+2) );
        end
        tmp=(2*p-1)*(2*p-2)/ gamma(2-beta).*tmp.*((abs(xt(:,1)-x(j,1))).^(2*p-1));
        z(:,j)=tmp;
end