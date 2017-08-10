clear all, close all, clc


addpath('build_mats')

parms.Re = 200; %Reynolds #
parms.deform = 'T'; %Deforming body?
parms.R_rho = 50;
parms.R_E = 1e5;
parms.R_sh = 0.35 * 1e-5;
parms.R_th = 0.01;
parms.inverted = 'T'; %inverted configuration?
parms.clamped = 'T'; %clamped or pinned?
parms.U_body = -1; %velocity of body in lab frame


%---Build data matrix 

    %snapshots for POD modes
    tvect = 100000 : 50 : 164750;
    parms.tvect = tvect;

    for j = 1 : length( tvect )

        j
        ind = tvect( j )

        
        filename = ['output/ib' num2str(ind,'%7.7i') '.var'];

        [parms,soln] = read_fort_base( filename, parms );
        
        if j == 1
           
            display( 'getting matrices ...' )
            mats = get_mats_preproc( parms, soln );

        end

        %Get vel flux and circulation from stfn
        soln.q = mats.C * soln.s; %get vel flux

        vel = mats.M_vel * soln.q;
        
        xb = soln.xb - parms.xb0;

        nvel = length( vel );
        nb = length( soln.xb ) / 2;

        n = nvel + 4 * nb;

        if j == 1
           
            Xbg = zeros( n, length( tvect ) );
            
        end


        bgvect = [ vel; xb ; soln.vb ];

        Xbg(:, j ) = bgvect;

    end


%---

save('datamat.mat','Xbg','parms','-v7.3')