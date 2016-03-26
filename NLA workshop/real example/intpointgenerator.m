function [ intdata ] = intpointgenerator( N )
%POINTGENERATOR Summary of this function goes here
%   Detailed explanation goes here
global RBFdomain

h=2*pi/N;


%% interior point
% intdata=(regintpoints(h)+1)/2;
intdata2=(regintbndpoints(h)+1)/2;
r=domainfct(intdata2(:,1),intdata2(:,2));
intdata=(intdata2([find(r<0)],:)+1)/2;

%% plot
% plot(intdata(:,1),intdata(:,2),'rx');
plot(intdata(:,1),intdata(:,2),'gx');
hold on
figure(1)
axis equal;

end

