clear all, close all, clc



cmap = zeros(64,3);
cmap(31:34,:) = [1 1 1;1 1 1;1 1 1;1 1 1];
cmap(1:30,1) = linspace(0,1,30);
cmap(1:30,2) = linspace(0,1,30);
cmap(1:30,3) = 1;
cmap(35:end,1) = 1;
cmap(35:end,2) = linspace(1,0,30);
cmap(35:end,3) = linspace(1,0,30);


save('cmap.mat','cmap');