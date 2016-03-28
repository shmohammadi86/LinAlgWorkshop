function [ B, psi,Dpsi , psi2 ] = gen_time_matrix( N )
%GEN_TIMEBASIS Generate the time basis function using Legendre polynomials
%with DBC in interval [-1 1]
% psi is in S_N ; psi2 is in S_N^*
% B is the mass matrix for time direction

    [x,w]= legslb(N+1);
    [DLpoly,Lpoly]=lepolym(N,x);
        % Need n-1 polynomials as basis, k=0,1,..N-2

    psi = Lpoly(1:N,:) + Lpoly(2:end,:);
    Dpsi = DLpoly(1:N,:) + DLpoly(2:end,:);
    psi2 = Lpoly(1:N,:) - Lpoly(2:end,:);

    for i = 0 : (N-1)
        for j = 0 : (N-1)
            B(i+1,j+1) = (psi(j+1,:).*psi2(i+1,:))*w;
        end
    end
    B(abs(B)<1e-15) = 0;


end

