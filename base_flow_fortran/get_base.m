clear all, close all, clc 
format long

%%
tic

%Don't forget to run 
%       make clean
%       make ib
%in terminal to compile fortran code.

if ~exist( 'base_new.mat' )

    %Load guess
    if ~exist( 'ib0000000guess.var'  )
        error('Requires initial guess from fortran output entitled ib0000000.var')
    else
        params = read_parms( 'ib0000000guess.var' )
    end


    %parameters for structure
    params.R_rho = 50;
    params.R_E = 1e5;
    params.R_sh = 4.1e-6;
    params.R_th = 0.01;
    params.standard = 'F';
    params.pinned = 'F';

    %build indices to convert fortran variables to matlab variables
    params.indbgx = 1 : ( params.m+1)*params.n; 
    params.indbgy = 1 : ( params.n+1)*params.m;
    params = get_inds( params );

    %move initial guess to output to restart simulation from this guess...
    [~,~] = system('cp ib0000000guess.var output/');
    [~,~] = system('mv output/ib0000000guess.var output/ib0000000.var');

    %convert fortran variables to matlab variables
    [q0, f0, uib0, udib0] = getvars( 0 );
    qold = truncate( q0, params ); %Truncate velocity to remove 
                                   %overlapping regions and scale
                                   %to get physical velocity



    varold = [qold; uib0; udib0; f0];



else

    %If we already have initial guess in matlab format then start with that
    load( 'base_new.mat')
    
    %Overload structure variables if desired
    params.R_sh = 4.1e-6;
    
	varold = var;

end


%Get old q and f:
params.time_inc = 1000 * params.dt; %How long in time to run...
gold = advance( varold, params ); %Compute g using initial guess 
                                                                        
%Parameters for NR iteration                                        
tol = 1e-4;
params.err = max(abs(gold))/max(abs(varold))
rstrt = 20; gmtol = 1e-4; maxit_gm = 50; %Parameters for gmres

while params.err > tol
    
    %Rescale tolerance to be appropriate relative to vel:
    gmtolrel = gmtol * norm( varold ) / norm( gold ) 
    Dg_times = @(var) Jac_times( var, varold, gold, params );
    
    [Dginvxg, flg, relres, itr, resv] = gmres( Dg_times, gold, rstrt, ...
        gmtolrel, maxit_gm);
    
    
    relres
    
    varnew = varold - Dginvxg; %Update guess for equilibrium point
    
    gnew = advance( varnew, params ); 

    params.err = max(abs( gnew ) ) / max(abs( varold ) ); %Check error

    err_show = params.err

    %Prepare for next iteration:
    varold = varnew;
    gold = gnew;

    var = varold;

    save('base_new.mat','var','params')
    
end

var = varold;

save('base_new.mat','var','params')


%save output from fortran run as base that the matlab base flow solver can
%use...
ind = params.time_inc/params.dt;
str = ['ib' num2str(ind,'%7.7i') '.var'];
[~,~] = system(['cp output/', str, ' ./']);
[~,~] = system(['mv ', str, ' base.var']);

toc
