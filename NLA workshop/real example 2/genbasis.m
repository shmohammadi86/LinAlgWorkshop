function [ phi,Lambda ] = genbasis( BC , N )
%GENBASIS generate the eigenfunction as basis function
% Input : BC: 1 for Homogeneous Dirichlet,2 for Homogeneous Neumann
%         N : Number of grid points
%         


    [x,w]= legslb(N+1);
    [DLpoly,Lpoly]=lepolym(N,x);
        % Need n-1 polynomials as basis, k=0,1,..N-2
    %% Different BCs
    if BC == 1
        legbasis = Lpoly(1:N-1,:) - Lpoly(3:end,:);
        %% Scale legbasis to make S an identity matrix
        for k = 0 : N-2
            legbasis(k+1,:) = legbasis(k+1,:)/(sqrt(4*k+6));  % Scale factor
        end
    end

    if BC == 2
        legbasis(1,:) = 0.5*ones(1,N+1);
        for k = 1 : N-2
            bk(k+1,1) = (k+1)*k/((k+2)*(k+3));  
            legbasis(k+1,:) = (Lpoly(k+1,:) - bk(k+1,:).*Lpoly(k+3,:))/(sqrt((4*k+6)*bk(k+1)));
        end
    end

    %% Calculate the Mass Matrix M
    for i = 0 : (N-2)
        for j = 0 : (N-2)
            M(i+1,j+1) = (legbasis(j+1,:).*legbasis(i+1,:))*w;
        end
    end
    M(abs(M)<1e-15) = 0;

    %% Calculate eig and construct eigenfunctions phi
    [E,Lambda] = eig(M);
    for k = 0 : N-2 
        phi(k+1,:) = E(:,k+1)'*legbasis;   % Eigenfuntion
    end

end

