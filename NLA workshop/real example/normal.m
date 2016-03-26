function [nx,ny]=normal(x,y,RBFdomain)   % This defines domains by specifying the
% boundary radius r in polar coordinates.
% No other routine is needed for domain definition.
% Make sure that the domain is in [-1,1]
% and that it is star-shaped wrt. the origin.
switch RBFdomain
    case 1   % circle
        nx = 2*x;
        ny = 2*y;
    case 2   % the square [-1,1]^2
        nx = zeros(size(x));
        ny = zeros(size(y));
        idx = find (abs(x-1)<10*eps);  nx(idx) = 1; 
        idx = find (abs(x)<10*eps);  nx(idx) = -1;         
        idx = find (abs(y-1)<10*eps);  ny(idx) = 1; 
        idx = find (abs(y)<10*eps);  ny(idx) = -1;        
    case 3   % cardioid
        [t,r] = cart2pol(x,y);
        nx = 0.4*cos(t)+0.8*cos(t).^2-0.4;
        ny = 0.4*sin(t).*(2*cos(t)+1);        
    case 4   % five- pointed star
        [t,r] = cart2pol(x,y);
        nx = -1.50*sin(5.0*t).*sin(t)+(.45+.3*cos(5.0*t)).*cos(t);
        ny = 1.50*sin(5.0*t).*cos(t)+(.45+.3*cos(5.0*t)).*sin(t);
    case 5   % smooth peanut
        [t,r] = cart2pol(x,y);
        nx =   -0.60* sin(2.0 *t) .*sin(t) + (0.45 + 0.3 *cos(2.0 *t)).* cos(t);
        ny = .60*sin(2.0*t).*cos(t)+(.45+.3*cos(2.0*t)).*sin(t)
    case 6  % peanut with corners
        [t,r] = cart2pol(x,y);
        nx = -.30*sign(cos(1.0*t)).*sin(1.0*t).*sin(t)+(.15+.3*abs(cos(1.0*t))).*cos(t);
        ny = 0.30*sign(cos(1.0 *t)).* sin(1.0 *t).* cos(t) + (0.15 + 0.3 *abs(cos(1.0 *t))).* sin(t)
end

% Normalize and return
c = sqrt( nx.^2 + ny.^2 );
if c~=0
    nx = nx./c;
    ny = ny./c;
end