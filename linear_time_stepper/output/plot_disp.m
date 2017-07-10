clear all, close all, clc


%%
dispfile='/disp.dat';                      % Filename for Force Data
disps=importdata(strcat(pwd,dispfile));   % Import all data from file

dt=1e-3;                      
% Size of time step                         

% Cut out the initial transients from the plot
t0=50;                                      % Start time for output
i0=1;                                       % Find corresponding index
while disps(i0,1)*dt<t0
    i0=i0+1;
end
 
% Plot tip disp
figure(2), hold on
% subplot(2,1,1)

plot(disps(1:end,1)*dt,disps(1:end,3),'-');  % x-force in the body-fixed frame

grid on


grid off
% ylabel('$\Delta y_{tip}$','interpreter','latex','fontsize',16)

xlabel('$t$','interpreter','latex','fontsize',16)

set(gca,'ticklabelinterpreter','latex','fontsize',14)

% subplot(2,1,2)
    
    
    
forcefile='/force.dat';                      % Filename for Force Data
forces=importdata(strcat(pwd,forcefile));   % Import all data from file

dt=1e-3;                                 % Size of time step

% Cut out the initial transients from the plot
t0=50;                                      % Start time for output
i0=1;                                       % Find corresponding index
while forces(i0,1)*dt<t0
    i0=i0+1;
end

plot(forces(10:end,1)*dt,forces(10:end,2),'-');  % y-force in the body-fixed frame

plot(forces(10:end,1)*dt,forces(10:end,3),'-');  % y-force in the body-fixed frame

grid off
% ylabel('$C_L$','interpreter','latex','fontsize',16)

xlabel('$t$','interpreter','latex','fontsize',16)



set(gca,'ticklabelinterpreter','latex','fontsize',14)


figure(10) , subplot(2,1,1)
[w,e] = ezfft(disps(i0:end,1)*dt,disps(i0:end,3));
% 
loglog(w/(2*pi),e./max(abs(e)) ), hold on
set(gca,'ticklabelinterpreter','latex','fontsize',14)


figure(10) , subplot(2,1,2)
[w,e] = ezfft(forces(i0:end,1)*dt,forces(i0:end,3));
% 
loglog(w/(2*pi),e./max(abs(e)) ), hold on

[w,e] = ezfft(forces(i0:end,1)*dt,forces(i0:end,2));
% 
loglog(w/(2*pi),e./max(abs(e)) ), hold on
set(gca,'ticklabelinterpreter','latex','fontsize',14)

% axis([0 20 0 1])
% ax = gca;
% ax.XTick = 0 : 4 : 20;

