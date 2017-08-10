clear all, close all, clc


writerObj = VideoWriter('beamvid.avi');
writerObj.FrameRate=5;
open(writerObj);
for K = 40 : 150
  filename = sprintf('flag_%d.png', K);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);