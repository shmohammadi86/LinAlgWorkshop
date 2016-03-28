clear all

alpha=0.8;
beta=1.8;
line=linspace(0,pi,100);
[x,y]=meshgrid(line);
t=0;

tempx=0;
for n=4000:-1:1   %10^-14
    tempx=tempx+mlf(alpha,1,-n^beta*t^alpha,10)*sin(n*x)*(8*(-1)^(n+1)-4)/n^3;
end

tempy=0;
for n=4000:-1:1   %10^-14
    tempy=tempy+mlf(alpha,1,-n^beta*t^alpha,10)*sin(n*y)*(8*(-1)^(n+1)-4)/n^3;
end
temp=tempx.*tempy;
mesh(x,y,temp)