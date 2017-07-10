clear all, close all, clc

% Re = 200;
% Reh = 2;

L = 3; %length of first subdomain in x dir


circum = pi; %Circumference of the circle

M = 100

h = L/M


n = get_cyl_points(circum,h)

ds = build_cylinder(n)