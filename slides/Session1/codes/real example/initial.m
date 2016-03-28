function [ inival ] = initial( intdata,ctrs,idmat )
%INITIAL generate the initial value of lambda (lambda_0) under the initial condition===rhs
global bctype

bdydata=setdiff(ctrs,intdata,'rows');
n=length(bdydata);


if strcmp(bctype,'dbc')
    DM=idmat;
%     rhs=sin(pi*intdata(:,1)).*sin(pi*intdata(:,2));
%     rhs=1000*(intdata(:,1).^4-2*intdata(:,1).^3+intdata(:,1).^2).*(intdata(:,2).^4-2*intdata(:,2).^3+intdata(:,2).^2);
    rhs=(intdata(:,1).^2.*(pi-intdata(:,1))).*(intdata(:,2).^2.*(pi-intdata(:,2)));
BM=differenceMatrix(bdydata,ctrs);
    rhs=[rhs;zeros(n,1)];
    CM=[DM;BM];
    inival=CM\rhs;
    
elseif strcmp(bctype,'delta')
    DM=differenceMatrix(ctrs,ctrs);
    rhs=10^(-25)*exp(-1./(ctrs(:,1).^2-1/4)).*exp(-1./(ctrs(:,2).^2-1/4));
    inival=DM\rhs;
    
elseif strcmp(bctype,'nbc')
    N=100;
    DM=differenceMatrix(ctrs,ctrs);
    rhs=cos(pi*ctrs(:,1)).*cos(pi*ctrs(:,2));
%     rhs=1000*(ctrs(:,1).^4-2*ctrs(:,1).^3+ctrs(:,1).^2).*(ctrs(:,2).^4-2*ctrs(:,2).^3+ctrs(:,2).^2);
    inival=DM\rhs;
    
end

end

