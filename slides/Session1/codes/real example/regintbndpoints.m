function intpoints=regintbndpoints(h) % generate regular 
% interior, boundary, and outside data points
% with approximate fill distance h
n=ceil(sqrt(8)/h);   % adjust stepsizes for regular grid on [-1,1]^2
dh=2./n;
% dh=h;
[X,Y]=meshgrid(-4:dh:4+0.000001,-4:dh:4+0.000001);
[nx,ny]=size(X);
xr=reshape(X,[(nx)*(nx) 1]);
yr=reshape(Y,[(nx)*(nx) 1]);
zr=domainfct(xr,yr);
regpoints=[xr yr];     % this ends the generation of regular points

count=0;               % we count interior points first
for i=1:nx.*nx
    if zr(i,1)< 1000.*dh  % do take points on the boundary and outside 
        count=count+1;
    end
end
intpoints=zeros(count,2); % we now can allocate and store

count=1;                
for i=1:nx.*nx
    if zr(i,1)< 1000.*dh
        intpoints(count,:)=regpoints(i,:);  % and now we store.
        count=count+1;
    end
end