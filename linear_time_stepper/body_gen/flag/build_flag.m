function [ds, xhat, yhat] = build_flag(ds)



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%Build a plate of length 1
xtemp = 0 : ds : 1;
ytemp = zeros(size(xtemp));

n = length(xtemp);

xhat = zeros(size(xtemp)); yhat = xhat;
% th = 5.9*pi/180; %rotation angle ymax = 0.1

% th = 18*pi/180; %ymax = 0.3
% th = 0.7*pi/180; %ymax = 0.01
% th = 0.07*pi/180; %ymax = 0.001
th = 0;
for j = 1 : length(xtemp)
    
   R = [cos(th) -sin(th); sin(th) cos(th)];
   vect = [xtemp(j); ytemp(j)];
   vrot = R*vect;
   xhat(j) = vrot(1); yhat(j) = vrot(2);
   
end
length(xhat)

plot(xhat,yhat,'b')
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

% xhat2 = zeros(size(xtemp)); yhat2 = xhat2;
% 
% for j = length(xtemp) : -1 : 1
%    
%     xhat2(j) = (1 - (length(xtemp) - j)/length(xtemp)) * cos(0.1*pi);
%     yhat2(j) = (1 - (length(xtemp) - j)/length(xtemp)) * sin(0.1*pi);
%     
% end
% 
% hold on
% plot(xhat2,yhat2,'r--')

% xhat = xhat2; yhat = yhat2;

%Write the points to a data file:

fileID = fopen('body.001.inp','w');
fprintf(fileID,'%-6d \n',n);
fprintf(fileID,'%-1s \n','T');
for j = 1 : n
fprintf(fileID,'%-20.16f %-20.16f\n',xhat(j), yhat(j));
end
fclose(fileID);

% plot(xhat,yhat)

