% Script that performs Kansa collocation for 2D Laplace equation
% with Dirichlet BCs
% Test function:  f=x^2*y^2

clear all
global p
global beta
global dim
beta=1.8;
dim=2;
p=2;
point_type=0;
  a=6;
  B=zeros(a,1);
for xxx=1:a
    B(xxx)=2^(xxx+5)+1;
end
errmat=zeros(a,1);
H=zeros(a,1);
%   Number and type of collocation points
for NNN=6:a
    N=B(NNN);
    gridtype = 'h';
    dsites=getpointND(N,dim,point_type);
    neval = 50;
    intdata = pi*dsites;
    
    %% construct boundary points
    sn = ceil(sqrt(N)); bdylin = linspace(0,1,sn)';
    bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
    bdydata = pi*[bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1); ...
        flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
    h=1/(sn-1);
    
    %%
    ctrs = intdata;
    
    index0=find(intdata(:,1)==0 | intdata(:,2)==0 | intdata(:,1)==pi |intdata(:,2)==pi ) ;
    int0=intdata(index0,:);
    intdata = setdiff( intdata, int0, 'rows');
    % Create neval-by-neval equally spaced evaluation locations
    % in the unit square
    grid = pi*linspace(0,1,neval);
    [xe,ye] = meshgrid(grid);
    epoints = [xe(:) ye(:)];
    
    %   % Compute evaluation matrix
    EM= evaluateMatrix(epoints,ctrs);
    
    
    % Compute blocks for collocation matrix
    %%Interior points%%
    DM_intdata1 = approximate_functionx(intdata,ctrs);
    Ex_x=extraterm(intdata,ctrs,1);
    DM_intdata2 = approximate_functiony(intdata,ctrs);
    Ex_y=extraterm(intdata,ctrs,2);
    DM=DM_intdata1+DM_intdata2+Ex_x+Ex_y;
    %%boundary points%%
    BM=brbf(bdydata,ctrs);
    %   CM = [DM; BCM1; BCM2; BCM3; BCM4];
    CM=[DM;BM];
    %% Create right-hand side    %%for u=x^2*y^2
    rhs=[Lu(intdata);bdydata(:,1).^2.*bdydata(:,2).^2];
    %% RBF solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CM*lambda = rhs   =>  Ax=b %%%%%%%%%%%%%%%%%% 
    lambda=CM\rhs;
    
    Pf = EM * lambda;
    pf=reshape(Pf,neval,neval);
    exact=xe.^2.*ye.^2;
    
    error=exact-pf;
    figure(1)
    mesh(xe,ye,error); pause(1)
    % Compute maximum error on evaluation grid
    rms_err = norm(pf-exact)/neval;
    fprintf('RMS error:     %e\n', rms_err)
    fview = [-30,30];  % viewing angles for plot
    
    % error plot
    errmat(NNN,1)=rms_err;
    H(NNN,1)=h;
    % H(NNN,1)=1/T;
end
