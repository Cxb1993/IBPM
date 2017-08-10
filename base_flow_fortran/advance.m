function [g, soln_vars] = advance( vars0, params )

%Call the IBPM fortran code to advance so that it returns the solution

% [~,~] = system('rm output/*.chd');

dt = params.dt ;
time_inc = params.time_inc ;

%Write current variables to binary files for fortran code to read...
write_2_var( vars0, params)

%Build input file for code that only runs 1 time step:
it_stop = round( time_inc/ dt ); %stop the run at the end of the time increment
                                 %where we want to take a snapshot
it_sv = it_stop;  %Save a variable at that time

build_inp( it_stop, it_sv, params );

%Now run it

[~,~] = system('rm output/ib*.var');
[~,~] = system('rm output/*.dat');

% [~,~] = system('make clean');
% [~,~] = system('make ib');
% [~,~] = system('bin/ib');
system('bin/ib')

[q, f, uib, udib] = getvars( it_stop );


qsm = truncate( q, params );


soln_vars = [qsm; uib; udib; f];


g = vars0 - soln_vars;





