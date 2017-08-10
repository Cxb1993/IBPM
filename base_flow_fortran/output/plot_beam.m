clear all,  clc, close all

for j = 40000 : 1000 : 66750

    t = j
    
    %identify force file
    no_smooth_str = strcat('total_force_',num2str(t,'%07i'),'.dat');
    
    %load it
    load(no_smooth_str)
    
    no_smooth_var = strcat('total_force_',num2str(t,'%07i'));
    no_smooth_true = eval(no_smooth_var);
    
    xb= no_smooth_true(:,1);
    yb = no_smooth_true(:,2);
    
    
    
    figure(3), hold on
    plot(xb, yb, 'k-','linewidth',2)
    
    axis([-0.5 1.5 -1 1])
    
    set(gca,'Xtick',[],'Ytick',[],'linewidth',2)
    box on
    
    
end
