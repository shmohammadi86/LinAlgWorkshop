% The function s=ledisctran(n,x,w,p,iflag) performs the discrete 2D-Legendre transforms 
% between the physical space (i.e., physical values) and frequence space
% (Legendre expansion coefficients) at the Legendre-Gauss-Lobatto points
% Input:
%  n, x,w--- number of LGL points in x, where (x,w) can be computed by
%          [x,w]=legslb(n). Note: x,w are column vectors  
%           x_t x is the Gauss point in time and space respectively
%  iflag==0--- forward transform  
%    p--- (input) physical values at collocation points
%    s--- (output) expansion coefficients 
%  iflag not= 0--- backward transform  
%    p--- (input) expansion coefficients 
%    s--- (output) physical values at collocation points 
%
%  See Page 101 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepolym()   
% N   % Time grid number
% M   % 	Space grid number
%  Written by Duo Cao

function s=ledisctran2d(M,N,x,w,t,w_t,un,iflag)
% compute the Legndre polynomials up to order M-1. Note: T(i,j)=L_i(x_j)
 T_s=lepolym(M,x);  
 T_t=lepolym(N,t);
 
 if iflag==0, 
     temp=T_t*diag(w_t)*un*diag(w)*T_s';
     s= temp.*([[0:N-1]'+0.5;(N)/2]*[[0:M-1]'+0.5;(M)/2]'); % see (3.193)
     return;
 end
 
 s = T_t'*un*T_s; 
 return
 
 
 