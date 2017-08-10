clear all, close all, clc
%%
for j = 1 : 3
    
    f = load(['force_lin_',num2str(j),'.dat']);
    ind = f(end,1) * (j-1);
    
    figure(1), hold on
    subplot(2,1,1), hold on
    plot((f(:,1) +ind) * 0.002 + 0.198, f(:,2) )
    
    subplot(2,1,2), hold on
    plot( (f(:,1)+ind) * 0.002 + 0.198, f(:,3))
    
    
    d = load(['disp_',num2str(j),'.dat']);
    
    
    figure(10), hold on
    plot( (d(:,1)+ind)*0.002, d(:,3) )
    
end
