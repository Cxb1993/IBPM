clear all, close all, clc

addpath('build_mats')

%compute a base flow satisfying the steady equations of motion using a
%Newton-Raphson method.

%NOTE: requires an initial guess obtained through a previous base
%flow computation either in Matlab or fortran. If obtained in Matlab the
%file must be called "base.mat", and if from fortran "base.var"

file_start_mat = 'base.mat';
file_start_fort = 'base.var';


if exist( file_start_mat, 'file') == 2 
    load( 'base.mat' )
    
elseif exist(file_start_fort,'file') == 2
    
    [parms, soln] = read_fort_base( file_start_fort );
    
else
    error( ['Requires an initial guess for base flow! ',...
        'Provide base flow either as base.var (fortran output) or', ...
        ' base.mat (matlab output)'])
end
    

%If desired, overwrite parameters from file
parms.Re = 200; %Reynolds #
parms.U_body = -1; %velocity of body in lab frame
parms.tol = 1e-8; %tolerance for terminating N-R

%Parameters for flag (only needed if parms.deform == 'T')
parms.deform = 'T'; %Deforming body?
parms.R_rho = 50;
parms.R_E = 1e5;
parms.R_sh = 4.2e-6;
parms.R_th = 0.01;
parms.inverted = 'T'; %inverted configuration?
parms.clamped = 'T'; %clamped or pinned?


%run main function
main_fun( parms , soln )

