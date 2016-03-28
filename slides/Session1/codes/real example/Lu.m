function [ Lu ] = Lu( xt )
%LU is the right hand side of the system  which is the function
%F=(x^(-1/2)*y+y^(-1/2)*x)/gamma(2-beta)  beta=1.5
global beta

% Lu=2.*xt(:,1).^(2-beta).*xt(:,2).^2/gamma(3-beta)+2.*xt(:,2).^(2-beta).*xt(:,1).^2/gamma(3-beta)...
%     +25*xt(:,1).^2.*xt(:,2).^2;

Lu=2.*xt(:,1).^(2-beta).*xt(:,2).^2/gamma(3-beta)...
    +2.*xt(:,2).^(2-beta).*xt(:,1).^2/gamma(3-beta);

end

