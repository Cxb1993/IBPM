clear all, close all, clc
 

%%
dispfile='/disp.dat';                      % Filename for Force Data
disps=importdata(strcat(pwd,dispfile));   % Import all data from file

dt=1e-3;                      
% Size of time step

vels = zeros(size(disps));
vels(:,1) = disps(:,1);

vels(2:end,2) = disps(2:end,2) - disps(1:end - 1,2);
vels(:,2) = vels(:,2) / dt;

vels(2:end,3) = disps(2:end,3) - disps(1:end - 1,3);
vels(:,3) = vels(:,3)/ dt; % Plot tip disp


figure(2), hold on
% subplot(2,1,1)

plot(disps(1:end,1)*dt,disps(1:end,3),'-');  % x-force in the body-fixed frame
plot(vels(1:end,1)*dt,vels(1:end,3),'-');  % x-force in the body-fixed frame

grid off;
set(gca,'ticklabelinterpreter','latex','fontsize',14)

dy = disps(40000:end,3);
uy = vels(40000:end,3);

figure(3), hold on
plot(disps(:,3), vels(:,3))
set(gca,'ticklabelinterpreter','latex','fontsize',14)


figure(4), hold on
poinc_sect = dy(uy >= -1e-2 & uy <= 1e-2);

plot(zeros(size(poinc_sect)), poinc_sect, 'o' )
set(gca,'ticklabelinterpreter','latex','fontsize',14)
