
function [uevct,vevct] = getvel_evect(evect, lev, len, m, n)


u1_1 = transpose(reshape(evect(1:(m+1)*(n)),m+1,n));
v1_1 = transpose(reshape(evect((m+1)*(n)+1:(m+1)*(n)+(n+1)*(m)),m,n+1));

% ---- Compute grid related parameters ------------------------------------

% If we are not considering the smallest grid level, we need to compute the
% grid spacing and x and y offsets for the current grid level. 
%
% Note:  the cells in each grid level are twice as large in both directions
% as the ones of the previous grid level



offsetx = 1;
offsety = 1.5;

fac = 2^(lev-1);                                        
delta = len ./ m *fac;                                  % Grid spacing in both directions for current grid
offx = 2^(lev-1) * len/2 - len/2 + offsetx;             % Offset in x direction for current grid
offy = 2^(lev-1) * (n*len/m)/2 - (n*len/m)/2 + offsety; % Offset in y direction for current grid


% Need to divide by cell face length 'delta' to convert to velocity

    u = u1_1 ;
    v = v1_1 ;




% ---- Interpolate all variables to the same grid -------------------------

% Note: All variables are interpolated to the cell vertices

% Create the grid that will be used for all variables
[xn, yn] = meshgrid(delta:delta:(m-1)*delta, delta:delta:(n-1)*delta);
xn = xn - offx;
yn = yn - offy;

% Grid for x-velocities (vertical cell faces)
[xu, yu] = meshgrid(0:delta:m*delta, delta/2:delta:(n-0.5)*delta); 
xu = xu - offx;
yu = yu - offy;

% Grid for y-velocities (horizontal cell faces)
[xv, yv] = meshgrid(delta/2:delta:(m-0.5)*delta, 0:delta:n*delta);    
xv = xv - offx;
yv = yv - offy;


% Interpolate all variables accordingly to xn, yn
uevct = interp2(xu,yu,u,xn,yn);
vevct = interp2(xv,yv,v,xn,yn);



% % ---- Rotate the grid and the bodies to the lab-frame --------------------
% 
% xnp = rox + (xn-rox)*cos(theta) - (yn-roy)*sin(theta);
% ynp = roy + (xn-rox)*sin(theta) + (yn-roy)*cos(theta); 
% xbp = rox + (xb-rox)*cos(theta) - (yb-roy)*sin(theta);
% ybp = roy + (xb-rox)*sin(theta) + (yb-roy)*cos(theta);
% xn = xnp; yn = ynp; xb = xbp; yb = ybp;
% 
% unp = un*cos(theta) - vn*sin(theta);
% vnp = un*sin(theta) + vn*cos(theta); 
% u0pnp = u0pn*cos(theta) - v0pn*sin(theta);
% v0pnp = u0pn*sin(theta) + v0pn*cos(theta); 
% u0rnp = u0rn*cos(theta) - v0rn*sin(theta);
% v0rnp = u0rn*sin(theta) + v0rn*cos(theta);
% un = unp; vn = vnp;
% u0pn = u0pnp; v0pn = v0pnp;
% u0rn = u0rnp; v0rn = v0rnp;

end
