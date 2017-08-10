clear all, close all, clc
format long

%%
tic

%Compute the linear and adjoint modes of the linearized and adjoint N-S
%equations
params.Re = 200;

params.m = 200; %# of x points
params.n = 64; %# of y points
params.len = 2;

params.indbgx = 1 : ( params.m+1)*params.n; 
params.indbgy = 1 : ( params.n+1)*params.m;

params.nq = (params.m+1)*params.n + (params.n+1)*params.m; %Number of points for velocity 

params.nb = 50 ; %Number of body points
params.nf = 2*params.nb; %Number of points in force vector

params.mg = 5;

params.offx = 0.2;
params.offy = 0.32;

params.R_rho = 30;
params.R_E = 1e3;
params.R_sh = 1e-7;
params.R_th = 0.01;



params.lin = 'T';
params.adjoint = 'F';
params.pinned = 'T';
params.standard = 'T'; %standard flag or inverted configuration?

params = get_inds( params );

% n = nq*mg + nomga * mg + nf; %Total # of points in eigenvector
params.ntot = params.nvel_tot +  params.nf + 3*(3*params.nb); %Total # of points in eigenvector
 
params.dt = 0.001;
params.time_inc = 0.01; %Sample flow field at some fraction of a convective time unit.

%Want the k largest modes
k = 1;


%Generate a compatible vector by using randn to create a random vector
%and advancing 1 time step in the linear code to get the compatible force and
%velocity...

vars0 = get_varinit( params );


cp = 1

for j = 1 : 3

    %Compute the eigenvectors and eigenvalues...
    vars0 =  lin_advance( vars0, params ) ;
    [~,~] = system(['mv output/force_lin.dat force_lin_',num2str(j),'.dat']);
    [~,~] = system(['mv output/disp_lin.dat disp_',num2str(j),'.dat']);
end


toc
