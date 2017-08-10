clear all, close all, clc

load('datamat.mat')
dt = parms.dt;

%---get DMD modes from data matrix

    %Rescale X for covariance matrix:
    X = Xbg(:, 1 : end - 1);
    X2 = Xbg(:, 2:end );
    
    [U, S, V] = svd( X, 'econ' );

    %Compute DMD (Phi are eigenvectors ) 
    r = 15; % truncate at r ;

    U = U(:, 1 : r);
    S = S(1 : r, 1 : r);
    V = V(:, 1 : r);
    
    Atilde = U' * X2 * V * inv(S);

    [W, eigs] = eig(Atilde);
    Phi = X2 * V * inv(S) * W;
    
    %sort eigenvalues in order of decreasing magnitude
    [eig_srt, ind_srt] = sort( diag( eigs ), 'descend' );
    Phi_srt = Phi(:, ind_srt );
          
%---      


 
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
b = Phi_srt \ x0;

Phi_srt_use = zeros(size(Phi_srt));

%---postprocessing: plot POD modes
    count = 1;
    figcount = 1;
    for j = [ 4, 2 ]%,1, 3, 5, 6, 7 ]%1 : 1 : r
        
        Phi_srt_use(:,j) = Phi_srt(:,j) * b(j);
        
        for jj = 1 : 2
        
            if jj == 1
                Uvect = imag(Phi_srt_use(:, j ));
                title = ['Imag({\boldmath${\phi}$}$_',num2str(figcount),')$'];
            else
                Uvect = real(Phi_srt_use(:,j));
                title = ['Real({\boldmath${\phi}$}$_',num2str(figcount),')$'];

            end
                    

        qv = mats.M_vel_inv * Uvect( 1 : nq );
        
        %Convert vel to vorticity
        gamm = mats.R * qv;
        wnv = mats.M_vort * gamm;
        
        xbv = Uvect( nq + (1 : nb ) );
        ybv = Uvect( nq + nb + (1: nb ) );
        
        figsavename = 'DMD_invLC';
        plot_mode_combo( wnv, xbv, ybv, parms, count, figsavename, title )

        count = count + 1

%         plot_mode( wnv, xbv, ybv, parms, j )
%         
%         plot_mode_flag( xbv, ybv, parms, j )

        end
        
        figcount = figcount + 1;


    end
    
%     figure(123456)
%     theta = (0 : 1 : 100) * 2 *pi/100;
%     plot(cos(theta), sin(theta), 'k--');
%     hold on, grid on
%     scatter(real(eig_srt), imag(eig_srt), 'ok')
%         
%%
    omega =  log( eig_srt) / 0.05;
    imom = imag(omega ) / (2*pi);
    reom = real(omega);
    
    om_sh = reom + 1i * imom
    
    figure(4323422)
    plot( reom, imom ,'k.','markersize', 30), hold on
    vlin = linspace(-5,5,20);
    plot( zeros(size(vlin)), vlin ,'k--','linewidth', 0.4)
    plot( vlin, zeros(size(vlin)) ,'k--','linewidth', 0.4)
    axis([-0.5 .5 -1.25 1.25])
    set(gca,'Xtick',[-.45 0 0.45], 'Ytick', [-1.25 0 1.25])
    xlabel('Real($\omega$)','interpreter','latex','fontsize',28)
    ylabel('Imag($\omega$)','interpreter','latex','fontsize',28)
    set(gca,'ticklabelinterpreter','latex','fontsize',28)
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'Color', [1 1 1]); % Sets figure background
    set(gca, 'Color', [1 1 1]); % Sets axes background
    set(gcf, 'PaperUnits', 'centimeters' )
    set(gcf, 'PaperSize', [10.5 10.5] )
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 10.5 10.5] )
    set( gcf, 'PaperPosition', [0 0 10.5 10.5] )


    dirsv = '~/Google Drive/Data processing/paper/Figures/';
    print('-dpdf',[dirsv,'eigs_invLC.pdf'],'-r300')

        
%---
%%
close all, clc

% omega =  log( eig_srt) / dt;
%     imom = imag(omega ) / (2*pi);
%     reom = real(omega);
%     
%     om_sh = reom + 1i * imom
%     
    
    %# of x-vel (flux) points
nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
%# of y-vel (flux) points
nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
%Total # of vel (flux) points
nq = nu + nv;
%# of vort (circ) points
ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
nb = parms.nb ;

[mats.M_vel, mats.M_vel_inv] = get_M_vel( parms );
mats.M_vort = get_M_vort( parms );
[mats.C, mats.R] = get_curls( parms );


%---time simulation of system

    time = [0, 1.3, 1.6, 2.8];
    x0 = X(:,1);
    b = Phi_srt \ x0;
    ind = [1 : 6];%, 2 : 11];
    
    for j = 1 : length( time )
       
        soln = (Phi_srt(:,ind)) * ( exp( omega(ind) * time(j) ) .* b(ind) );
        
        Uvect = real(soln);
        
        qv = mats.M_vel_inv * Uvect( 1 : nq );
        
        %Convert vel to vorticity
        gamm = mats.R * qv;
        wnv = mats.M_vort * gamm;
        
        xbv = Uvect( nq + (1 : nb ) );
        ybv = Uvect( nq + nb + (1: nb ) );

        
%         figure(1), clf
%         plot_mode( wnv, xbv, ybv, parms, 1 )
%         plot_mode_flag( xbv, ybv, parms, 1 )
%         
%         pause%(0.05) 
        

        plot_mode_reconstruct( wnv, xbv, ybv, parms, j, 'DMD_reconstruct_invLC' )


        
    end



%---
          
          