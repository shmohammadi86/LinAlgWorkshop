
function [evec,eval,Q,T] = lanczosqr(A,k,v,p)
%
%   Input:  A -- an n by n matrix  (A = A' assumed)
%           k -- a positive integer (k << n assumed)
%           v -- an n  vector (v .ne. 0 assumed)
%
%   Output: Q -- an n by k orthogonal matrix
%           T -- a k by k symmetric tridiagonal matrix
%           z -- an n vector
%
% 
%           with   AQ = QT + ze_k'

    n = length(v);
    T = zeros(k);
    Q = zeros(n,k);
    v1 = v/norm(v);
    z = A*v1;
    alpha = v1'*z;
    z = z - v1*alpha;
    Q(:,1) = v1; 
    T(1,1) = alpha;
    for j = 2:k,
        beta = norm(z); 
        v0 = v1; v1 = z/beta;        
        z = A*v1 - v0*beta;
        alpha = v1'*z;
        z = z - v1*alpha;
        T(j,j-1) = beta; T(j-1,j) = beta; T(j,j) = alpha;
        Q(:,j)   = v1;
    end 
 %QR algorithm to compute eigenpairs       
T
 %size(T)
 TQ=T;
evec=1;
for k=1:p
    [q,r]=qr(TQ);
    evec=evec*q;
    TQ=r*q;
end
eval=diag(TQ);