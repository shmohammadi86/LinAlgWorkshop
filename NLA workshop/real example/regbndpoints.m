function bndpoints=regbndpoints(n) % defines "regular" boundary points
% using approx. n-1 equidistant angles
% in [0,2*pi] and calling bndfct.m.

 h=2*pi/n;        % to make the angles lie at equidistant values
 phi=0:h:(2*pi-h+0.00000001) ;% careful: no duplication of last=first point
 [nx,nphi]=size(phi);
 z=zeros(nphi,2);
 [r,phi]=bndfct(phi);
 [x,y]=pol2cart(phi,r);
 bndpoints=[x' y'];
