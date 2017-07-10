function dist = build_cylinder(n)

%dist is the grid spacing between the cylinder surface


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
% Define the IB to be a cylinder with radius 1
% centered at 0
ds = 2*pi/n;
spt = 0:ds:(n-1)*ds;
xfun = @(z) 0.5*cos(z);
yfun = @(z) 0.5*sin(z);
xhat = xfun(spt);
yhat = yfun(spt);

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

%Write the points to a data file:

fileID = fopen('body.001.inp','w');
fprintf(fileID,'%-6d \n',length(xhat));
fprintf(fileID,'%-1s \n','F');
for j = 1 : n
    if j == 2
        dist = sqrt( (xhat(j) - xhat(j-1))^2 + (yhat(j) - yhat(j-1))^2 )
    end  
fprintf(fileID,'%-20.16f %-20.16f\n',xhat(j), yhat(j));
end
fclose(fileID);

plot(xhat,yhat)

