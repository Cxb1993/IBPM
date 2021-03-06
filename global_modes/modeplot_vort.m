clear all, close all, clc


%%
addpath('../base_flow/build_mats/')

% Contour maximum values and number of contour levels
% Vorticity
% cmax_w = 3;
% clev_w = 20;
% clevs = linspace( -cmax_w, cmax_w, clev_w );

% Range for plots
range = [-0.5 5 -0.5 0.5];

load('cmap.mat')

load('../rigid_base/base.mat')

load('modes.mat')

D

k = 1;

%--Various variables 
    %# of x-vel (flux) points
    nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
    %# of y-vel (flux) points
    nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
    %Total # of vel (flux) points
    nq = nu + nv;
    %# of vort (circ) points
    ngam = get_vort_ind( parms.m-1, parms.n-1, parms.mg, parms );
    %# of body points
    nb = parms.nb;
    %# of surface stress points
    nf = nb * 2;
%--

%--specify matrices for ease in ensuing code

    if exist('RC.mat', 'file') ~= 2
        mats = get_RC( parms );
        save( 'RC.mat', 'mats')
    else
        load('RC.mat')
    end

    C = mats.C; R = mats.R; 

%--

for jj = 1 : 2
    
    if jj == 1
        Vv = real( V( 1 : ngam, k) );
    else
        Vv = imag( V( 1 : ngam, k ) );
    end


    
    gamma = R * C * Vv;
    
    Xv = zeros( parms.n-1, parms.m-1, parms.mg );
    Yv = Xv;
    Omega = Xv;

    %Get x-y points and vorticity
    for lev = 1 : parms.mg

        fac = 2^(lev-1);

        % Grid spacing in both directions for current grid
        delta = parms.len ./ parms.m *fac;

        % Offset in x direction for current grid
        offx = 2^(lev-1) * parms.len/2 - parms.len/2 + parms.offx ;

        % Offset in y direction for current grid
        offy = 2^(lev-1) * (parms.n*parms.len/parms.m)/2 - ...
            (parms.n*parms.len/parms.m)/2 + parms.offy ;


        %--get grid points


            xv = delta : delta : (parms.m-1) * delta;
            yv = delta : delta : (parms.n-1) * delta;

            xv = xv - offx;
            yv = yv - offy;

            [Xv(:,:,lev), Yv(:,:,lev)] = meshgrid( xv, yv );
 
        %--
    
        %--get vorticity

        if lev == 1

            ind_s = 1;
            ind_e = get_vort_ind( parms.m-1, parms.n-1, lev, parms );
            omega = gamma(ind_s : ind_e ) / delta^2;

%             omega(omega > cmax_w ) = cmax_w;
%             omega(omega < -cmax_w ) = -cmax_w;

            Omega(:,:,lev) = transpose( reshape( omega, parms.m-1, parms.n-1 ) );                                  

            %store omega for coarse-grid interpolation
            omega_s = omega;
        else
            
            omega = zeros( (parms.m-1) * (parms.n-1), 1 );
            
            %--need some fancy indexing for what follows:
            
                %get indices of current grid (indices for overlapping
                %region are zero)
                indvortx = repmat(1 : parms.m-1,[1,parms.n-1]);
                indvorty = repelem(1 : parms.n-1, parms.m-1);
                vort_ind = get_vort_ind(indvortx,indvorty,lev,parms);
            
                %get all indices on gridlevel 1
                indvortx = repmat(1 : parms.m-1,[1,parms.n-1]);
                indvorty = repelem(1 : parms.n-1, parms.m-1);
                vort_ind_f = get_vort_ind(indvortx,indvorty,1,parms);
                
                %Get indices on gridlevel 1 just for overlapping part
                indvortx = repmat(2 : 2 : parms.m-2,[1,parms.n/2-1]);
                indvorty = repelem(2 : 2 : parms.n-2, parms.m/2-1);
                vort_ind_s = get_vort_ind(indvortx,indvorty, 1, parms);
                
                %indexing on the vector gamma starts at 1...
                n_sub = get_vort_ind( parms.m-1, parms.n-1, lev-1, parms);
                
                %non-overlapping indices:
                ind_no_over = ( (vort_ind - n_sub) > 0  );
                ind_no_over = vort_ind_f( ind_no_over ~= 0 );
                
                %overlapping indices:
                ind_over = ( (vort_ind - n_sub) < 0  ); 
                ind_over = vort_ind_f( ind_over ~= 0 );
            
            %--
            
            
            
            %Get part that doesn't overlap with finer grid
            omega( ind_no_over ) = gamma( vort_ind( vort_ind ~= 0 ) )/ delta^2;
            
            
%             %for overlapping part...
            omega( ind_over ) = omega_s( vort_ind_s );
%                         
            Omega(:,:,lev) = transpose( reshape( omega, parms.m-1, parms.n-1 ) );
            
            %save omega for next coarse grid
            omega_s = omega;

        end

        %--



    end
    
    %plot them
    
    figure(1), subplot(2,1,jj)
    
    cmax_w = max(max(max(abs( Omega(:,:,1) ))));
    clevs = linspace( -0.2, 0.2, 20 ) ;
    Omega( Omega > max(clevs) ) = max(clevs);
    Omega( Omega < min(clevs) ) = min( clevs );

    for j = parms.mg : -1 : 1
        
        contourf(Xv(:,:,j), Yv(:,:,j), Omega(:,:,j), clevs, ...
            'edgecolor','none'); shading flat;
        
        colormap( cmap )
        axis equal
        axis(range)
        set(gca, 'TickLabelInterpreter','latex','fontsize',14)
        box on
        if jj == 1
            str = 'Re($\omega_p$)';
        else
            str = 'Im($\omega_p$)';
        end
        
        title(str, 'interpreter','latex','fontsize',18)
        hold on        
    end
    
    %plot body
    plot(soln.xb( 1 : parms.nb ), soln.xb( 1+parms.nb : 2*parms.nb ),'k', ...
        'linewidth', 2)

end
