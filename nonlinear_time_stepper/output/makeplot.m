% Plots the Vorticity field and velocity fields around the body
% Plots the x and y forces on the body
clear all, close all, clc
 
%%
% -------- General Parameters ---------------------------------------------
dt=1e-3;


% Number of the time step of the file: this identifies which file should be read
istart = 10000; iend = 12000; inv = 500;
 

% Contour maximum values and number of contour levels
% Vorticity
cmax_w = 5;
clev_w = 18; %use even #
clevs = linspace( -cmax_w, cmax_w, clev_w );
cmap = build_cmap( length(clevs) );

mg = 5; %# of grid levels used for simulation

% Range for plots
range = [-2 10 -4 4];

% Plot the pressure field?
pressure=0;
 



for it = istart : inv : iend

    figure(1), clf
    
    for j = mg : -1 : 1
        
        % Get data from specified file (see getdata for details)
        [xb,yb,codeb,xn,yn,un,vn,u0pn,v0pn,u0rn,v0rn,wn,sn,pn] ...
            = getdata(pwd,it,j,pressure);
        
        colormap(cmap);
        
        wn( wn > cmax_w ) = cmax_w;
        wn( wn < -cmax_w ) = -cmax_w;
%         wn( abs(wn) < ( min(abs(clevs)) - min(abs(clevs))/2 ) ) = 0;
        
        
        figure(1)
        contourf(xn,yn,wn,'edgecolor','none'); shading flat;
        hold on
    end
    
    % Plot body 
    fill(xb(find(codeb==1)),yb(find(codeb==1)),'k');
   
    axis equal
    axis(range);
    set(gca,'ticklabelinterpreter','latex',...
        'fontsize',14)
    box on
    
    pause
    
end


% -------- Plot Forces ----------------------------------------------------
%%
forcefile='/force.dat';                      % Filename for Force Data
forces=importdata(strcat(pwd,forcefile));   % Import all data from file

dt=3e-3;                                 % Size of time step                         

% Cut out the initial transients from the plot
t0=0;                                      % Start time for output
i0=1;                                       % Find corresponding index
while forces(i0,1)*dt<t0
    i0=i0+1;
end
    
figure(3), hold on

% Plot X force
subplot(2,1,1), hold on
plot(forces(i0:end,1)*dt,forces(i0:end,2),'-');  % x-force in the body-fixed frame
ylim([0 5])
grid on
title('$C_D$','interpreter','latex')

% Plot Y Force
subplot(2,1,2), hold on
plot(forces(i0:end,1)*dt,forces(i0:end,3),'-');  % y-force in the body-fixed frame
grid on
ylim([-0.3 0.3])
title('$C_L$','interpreter','latex')

xlabel('t/T','fontsize',16)





