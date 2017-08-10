function [] = forward( it_stop, it_sv, params  )

%Build input file for code that only runs 1 time step:

build_inp( it_stop, it_sv, params );

%Now run it

[~,~] = system('rm output/ib_lin*.var');
[~,~] = system('rm output/*.dat');

[~,~] = system('rm src/user.f90');
[~,~] = system('cp user.f90 src/');


[~,~] = system('make clean');
[~,~] =system('make ib');
[~,~] =system('bin/ib');