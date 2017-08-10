clear all, close all, clc

%Test to make sure the linear solver from this directory is giving me
%physical results...

system('rm output/ib_lin*.var');
system('rm output/*.dat');
[~,~] = system('rm output/*.chd');

%Compute the linear and adjoint modes of the linearized and adjoint N-S
%equations
params.Re = 200;

params.m = 400; %# of x points
params.n = 200; %# of y points
params.len = 2;

params.nq = (params.m+1)*params.n + (params.n+1)*params.m; %Number of points for velocity 

params.nb = 100 ; %Number of body points
params.nf = 2*params.nb; %Number of points in force vector

params.mg = 5;

params.offx = 0.2;
params.offy = 0.5;

params.R_rho = 10;
params.R_E = 1e5;
params.R_sh = 3.3e-6;
params.R_th = 0.01;

% n = nq*mg + nomga * mg + nf; %Total # of points in eigenvector
n = params.nq*params.mg +  params.nf + 3*(3*params.nb); %Total # of points in eigenvector
 
params.dt = 0.001;
params.time_inc = 0.01; %Sample flow field at some fraction of a convective time unit.

m = params.m; %# of x points
n = params.len ;

nq = params.nq ; %Number of points for velocity 

nb = params.nb ; %Number of body points
nf = params.nf; %Number of points in force vector

mg = params.mg ;

dt = params.dt ;
time_inc = params.time_inc ;
len = params.len;

params.standard = 'F';
params.pinned = 'F';


params.lin = 'T';
params.adjoint = 'F';

[~,~] = system('rm output/q_init.var');

var0trial = zeros( nq, mg );

fileID = fopen('output/q_init.var','w','b');
for j = 1 : nq
    fwrite(fileID,var0trial(j,:),'real*8','b');
end
fclose(fileID);


uib0 = zeros(3*nb,1);
uib0(2:3:end) = linspace(0.01,0, nb);

uib0

[~,~] = system('rm output/uib_init.var');

fileID = fopen('output/uib_init.var','w','b');
fwrite(fileID,uib0,'real*8','b');
fclose(fileID);


udib0 = zeros(3*nb,1);


[~,~] = system('rm output/udib_init.var');

fileID = fopen('output/udib_init.var','w','b');
fwrite(fileID,udib0,'real*8','b');
fclose(fileID);

uddib0 = zeros(3*nb,1);


[~,~] = system('rm output/uddib_init.var');

fileID = fopen('output/uddib_init.var','w','b');
fwrite(fileID,uddib0,'real*8','b');
fclose(fileID);



fb0 = zeros(2*nb,1);

fileID = fopen('output/fb_init.var','w','b');
fwrite(fileID,fb0,'real*8','b');
fclose(fileID);




%Build input file for code...
    %We are flicking the code with a bodyforce from time 0.02 to 0.04. Run
    %until 0.06 to feel effect of body force.
timestp = round( 10000 ); 
isv = 100 ;
build_inp( timestp, isv, params );


% system('rm src/user.f90');
% system('cp user_pert.f90 src/');
% system('mv src/user_pert.f90 src/user.f90');


system('make clean');
system('make ib');
system('bin/ib');

%%
run output/makeplot.m



