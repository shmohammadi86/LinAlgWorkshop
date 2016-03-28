function [ DLpoly2 ] = lepolydiff2( DLpoly )
%LEPOLYDIFF2 Generate the second order derivative of Legendre polynomials
[M,N] = size(DLpoly);   % N is the polynomial order   M is the number of Gauss points

DLpoly2 = zeros(N,M);

for k = 2 : N-1;
    DLpoly2(k+1,:) = (2*k-1).*DLpoly(k,:) + DLpoly2(k-1,:) ;
end



end

