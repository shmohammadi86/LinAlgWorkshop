function y = afun (x)
    global psi psi2 w_t NonlinearTerm phi w
    global DiffOperator
    [N,~] = size(psi);
    [M,~] = size(phi);
    M = M +1 ;
    temp1 = DiffOperator*x;
    x = reshape(x,N,M-1);
    temp2 = psi'* x * phi;
    temp3 = psi2*diag(w_t)*(NonlinearTerm.* temp2)*diag(w)*phi'; 
    y = temp1 + temp3(:);
end