function [ rms_err ] = errev( pf )
%ERREV Summary of this function goes here
%   Detailed explanation goes here
global neval
global RBFdomain

grid = linspace(-1,1,neval); 
[xe,ye] = meshgrid(grid);
r=domainfct(xe,ye);
xe=(xe+1)/2;
ye=(ye+1)/2;

exact=xe.^2.*ye.^2;
error=exact-pf;
error(find(r>0))=0;

if RBFdomain==7
    a=length(xe)^2;
    for i=1:a
    if xe(i)>0.5 && ye(i)>0.5
        error(i)=0;
    end
    end
end
figure(1)
pause(1)
mesh(xe,ye,error)
rms_err = norm(pf-exact)/neval;
fprintf('RMS error:     %e\n', rms_err);
end

