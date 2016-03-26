function [ v ] = initialSol( M,N,x,w,phi,Lambda,gamma)
%INITIALSOL Use a 2nd order implicit time scheme to generate a guess for Newton
%method

v = zeros(N+1,M+1);
u0 = cos(pi.*x);
Vn = u0;   % u_n
Vn_1 = u0;  % u_{n-1}
v(1,:) = u0;
s = 2; % Time-step parameter, s>0
dt = 1/N;
k = 2;

while k< N+2
    rhs = 4*(1+dt*s+dt).*Vn - (1+2*dt*s+2*dt).*Vn_1 - 2*dt.*(2.*Vn.^3-Vn_1.^3);
    for i = 0 : (M-2)
        F(i+1,1) = (rhs'.*phi(i+1,:) )*w;
        Vtuta(i+1,1) = F(i+1,1)/(2*dt*gamma+(3+2*s*dt)*Lambda(i+1,i+1));
    end
    V = Vtuta' * phi;
    v(k,:) = V';
    
    Vn_1 = Vn;
    Vn = V';
    k = k+1;
end
end


% function [ v ] = initialSol( M,N,x,w,phi,Lambda,gamma ,epsilon)
% %INITIALSOL Use a 1st order time scheme to generate a guess for Newton
% %method
% 
% v = zeros(N+1,M+1);
% u0 = cos(pi.*x);
% v(1,:) = zeros(1,M+1);
% 
% vn = v(1,:)';
% dt = 0.01;
% k = 2;
% 
% while k< N+2
%     rhs = vn-dt*gamma/epsilon^2.*((vn+u0).^3-(vn+u0))-dt*gamma*pi^2.*cos(pi.*x);
%     for i = 0 : (M-2)
%         F(i+1,1) = (rhs'.*phi(i+1,:) )*w;
%         Vtuta(i+1,1) = F(i+1,1)/( Lambda(i+1,i+1) + dt*gamma);
%     end
%     V = Vtuta' * phi;
%     vn = V';
%     v(k,:) = vn;
%     k = k+1;
% end
% end

