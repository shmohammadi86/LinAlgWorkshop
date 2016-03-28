function points=getpointND(numpoints, dim, nrand)
% generates roughly numpoints points in dim dimensions
% random, if nrand>0, else regular in [-1,1]^dim
if nargin==1,
    dim=2;
    nrand=1;
end
if nargin==2,
    nrand=0;
end
if nrand~=0 % random case
    rand('state',nrand);
    points=rand(numpoints,dim); 
    return
end

% regular case
np_oneD=ceil(numpoints.^(1/dim)); % treat all dimensions alike
numpoints=np_oneD^dim;
h=1/(np_oneD-1);
if dim==1
    points=(0:h:1)';
    return
end
if dim==2
    [X,Y]=meshgrid(0:h:1,0:h:1);                 
    points=[X(:) Y(:)];
    return
end


