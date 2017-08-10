function var0 = get_varinit( params )

if ~exist('eig_guess.mat')


    %Build initial velocity and run for a time step
    m = params.m; %# of x points
    n = params.len ;

    nq = params.nq ; %Number of points for velocity 

    nb = params.nb ; %Number of body points
    nf = params.nf; %Number of points in force vector

    mg = params.mg ;

    dt = params.dt ;
    time_inc = params.time_inc ;
    len = params.len;



    [~,~] = system('rm output/q_init.var');

    var0trial = zeros( nq, mg );
    var0trial(:, 1 : mg-2) = randn( nq, mg-2);
    
    var0trial = var0trial / norm(var0trial) ; 

    fileID = fopen('output/q_init.var','w','b');
    for j = 1 : nq
        fwrite(fileID,var0trial(j,:),'real*8','b');
    end
    fclose(fileID);


    uib0 = zeros(3*nb,1);
%     uib0(2:3:end) = linspace(0,0.01, nb);


    [~,~] = system('rm output/uib_init.var');

    fileID = fopen('output/uib_init.var','w','b');
    fwrite(fileID,uib0,'real*8','b');
    fclose(fileID);


    udib0 = zeros(3*nb,1);


    [~,~] = system('rm output/udib_init.var');

    fileID = fopen('output/udib_init.var','w','b');
    fwrite(fileID,udib0,'real*8','b');
    fclose(fileID);


%     fb0 = zeros(2*nb,1);
% 
%     fileID = fopen('output/fb_init.var','w','b');
%     fwrite(fileID,fb0,'real*8','b');
%     fclose(fileID);

    uddib0 = zeros(3*nb,1);


    [~,~] = system('rm output/uddib_init.var');

    fileID = fopen('output/uddib_init.var','w','b');
    fwrite(fileID,uddib0,'real*8','b');
    fclose(fileID);


    %Build input file for code...
    timestp = round( 2000 ); 
    isv = timestp ;
    build_inp( timestp, isv, params );

    %Now run it
    [~,~] = system('rm output/ib_lin*.var');
    [~,~] = system('rm output/total_force_lin*.dat');
    [~,~] = system('rm output/total_disp_lin*.dat');
    [~,~] = system('rm output/cfl_lin.dat');
    [~,~] = system('rm output/slip_lin.dat');
    [~,~] = system('rm output/force_lin.dat');
    [~,~] = system('rm output/*.chd');



    [~,~] = system('rm src/user.f90');
    [~,~] = system('cp user.f90 src/');
    % [~,~] = system('cp user_pert.f90 src/');
    % [~,~] = system('mv src/user_pert.f90 src/user.f90');


%     system('make clean')
%     system('make ib')
    system('bin/ib')


    %Extract data after 1 time step and return this as IC for eigs:
    [q, f, uib, udib, uddib] = getvars( timestp );


    qsm = truncate( q, params );

%     var0 = [qsm ; uib; udib; f; uddib];

     var0 = [qsm ; uib; udib; uddib];


    var0 = var0/ max(max(abs((var0))));
    
    
    system(['mv output/ib_lin',num2str(timestp,'%07i'),'.var output/ib_lin0000000.var'])
    
    run output/makeplot.m
    
    pause

    save('eig_guess.mat', 'var0')
    
else
    load('eig_guess.mat')
    
    var0 = var0 / norm( var0 );
    
end
