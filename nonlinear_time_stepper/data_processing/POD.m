clear all, close all, clc

    addpath('build_mats')
    
%--get data matrix X and subtract mean
    load('datamat.mat')

    X = Xbg;
    parms.Kb = 0.35;
    parms.M_rho = 0.5;
    parms.Ke = 1000;

    %Subtract mean:
    Xtil = X;

    mean_X = mean( X, 2 );
    for j = 1 : length(parms.tvect)
        Xtil(:,j) = Xtil(:,j) - mean_X;
    end
%--

%---get weighting matrix for norm

    [W, Wsq] = build_W( parms );

%---      

%--get singular vectors V

    smallmat = (Xtil' * Wsq * Xtil );
    [V,Sig_sq] = eig( smallmat ); 
    
    %sort sing vals in order of decreasing magnitude
    [Sigsq, ind_srt] = sort( diag( Sig_sq ), 'descend' );
        
    V = V(:, ind_srt );
    
    
    max(max(abs( V * diag(Sigsq )  * V' - smallmat ) ) )
    
    
    Sig = real(sqrt( Sigsq ));
        
    Siginv = Sig;
    Siginv( Siginv ~= 0 ) = 1./ Siginv( Siginv ~= 0 );
        
    Siginv = diag(Siginv);
    
    
    %%
    close all, clc
    figure(4323422)
    semilogy( real( Sig ) / max(abs(real(Sig))), 'k.','markersize', 30)
    axis([0 50 1e-4 1])
    set(gca,'Ytick',[1e-4 1e-2 1], 'Xtick', [0 25 50])
    ylabel('$\sigma$','interpreter','latex','fontsize',28)
    xlabel('index','interpreter','latex','fontsize',28)
    set(gca,'ticklabelinterpreter','latex','fontsize',28)
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'Color', [1 1 1]); % Sets figure background
    set(gca, 'Color', [1 1 1]); % Sets axes background
    set(gcf, 'PaperUnits', 'centimeters' )
    set(gcf, 'PaperSize', [30 10.5] )
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 30 10.5] )
    set( gcf, 'PaperPosition', [0 0 15 10.5] )


    dirsv = '~/Google Drive/Data processing/paper/Figures/';
    print('-dpdf',[dirsv,'sig_invLC.pdf'],'-r300')
    
    %%
    U = W * Xtil * V * Siginv;
    
%     max(max(abs( U' * U - speye( length( U(1,:)), length(U(1,:) ) ) ) ) )
    
    Uhat = Xtil * V * Siginv;
    
%     max(max(abs( Uhat' * Wsq * Uhat - speye( length( U(1,:)), length(U(1,:) ) ) ) ) )
    
    
%     figure(432343232)
%     plot(V(:,1)), hold on
%     plot(V(:,2))
%     plot(V(:,3))
%     plot(V(:,4))
%--

%%
close all, clc
%# of x-vel (flux) points
nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
%# of y-vel (flux) points
nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
%Total # of vel (flux) points
nq = nu + nv;
nb = parms.nb;

[mats.M_vel, mats.M_vel_inv] = get_M_vel( parms );
mats.M_vort = get_M_vort( parms );
[mats.C, mats.R] = get_curls( parms );


x0 = X(:,1);
b = Uhat' * Wsq * x0;

U_use = zeros(size(Uhat));

%---postprocessing: plot POD modes
    for j =  1 : 4 %1 : 1 : r
        
        U_use(:,j) = Uhat(:,j) * b(j);
        
        Uvect = real(U_use(:, j ));
        

        qv = mats.M_vel_inv * Uvect( 1 : nq );
        
        %Convert vel to vorticity
        gamm = mats.R * qv;
        wnv = mats.M_vort * gamm;
        
        xbv = Uvect( nq + (1 : nb ) );
        ybv = Uvect( nq + nb + (1: nb ) );
        
        figsavename = 'POD_invLC';
        title = ['$\hat{\textbf{u}}_',num2str(j),'$'];

        plot_mode_combo( wnv, xbv, ybv, parms, j, figsavename, title )
        
%         plot_mode_flag( xbv, ybv, parms, j, 'POD' )
        
%         plot_mode( wnv, xbv, ybv, parms, j, 'POD' )
%         
%         plot_mode_flag( xbv, ybv, parms, j, 'POD' )

    end


%%

close all, clc

addpath('build_mats')

%# of x-vel (flux) points
nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
%# of y-vel (flux) points
nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
%Total # of vel (flux) points
nq = nu + nv;
nb = parms.nb;

[mats.M_vel, mats.M_vel_inv] = get_M_vel( parms );
mats.M_vort = get_M_vort( parms );
[mats.C, mats.R] = get_curls( parms );

itvect = [100000, 101300, 101600, 102800] / 50 - 2000 + 1;

mod_nu = 1 : 6;
for j = 1 : length( itvect )
    
    it = itvect( j );
    
    Uvect = mean_X + Uhat(:, mod_nu) * ( Uhat( :, mod_nu)' ...
        * ( Wsq * Xtil(:, it ) ) );
    
    
    qv = mats.M_vel_inv * Uvect( 1 : nq );
        
    %Convert vel to vorticity
    gamm = mats.R * qv;
    wnv = mats.M_vort * gamm;

    xbv = Uvect( nq + (1 : nb ) );
    ybv = Uvect( nq + nb + (1: nb ) );

    plot_mode_reconstruct( wnv, xbv, ybv, parms, j, 'POD_reconstruct_invLC' )

%     plot_mode_flag( xbv, ybv, parms, 1, 'POD_reconstruct' )
    
%     pause
    
end




