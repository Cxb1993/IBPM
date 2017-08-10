clear all, close all, clc


for j = 40000 : 250 : 66750

    t = j
    
    %identify force file
    no_smooth_str = strcat('total_force_',num2str(t,'%07i'),'.dat');
    
    %load it
    load(no_smooth_str)
    
    no_smooth_var = strcat('total_force_',num2str(t,'%07i'));
    no_smooth_true = eval(no_smooth_var);
    
    xhat = no_smooth_true(1,1);
    yhat = no_smooth_true(1,2);
    
    ux = no_smooth_true(1,3); 
    uy = no_smooth_true(1,4); 
    
    figure(1), hold on
    plot(yhat,uy,'bx')
   
%     pause
end