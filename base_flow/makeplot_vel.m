clear all, close all, clc


%%
addpath('./build_mats/')

% Contour maximum values and number of contour levels
% Vorticity
cmax_w = 2;
clev_w = 20;
clevs = linspace( -cmax_w, cmax_w, clev_w );

% Range for plots
range = [-4 10 -5 5];

load('cmap.mat')


    load(['base.mat'])

    q = soln.q + parms.q0;
    q = q( 1 : get_velx_ind( parms.m-1, parms.n, parms.mg, parms ) );
    
    Xv = zeros( parms.n, parms.m-1, parms.mg );
    Yv = Xv;
    U = Xv;

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
            yv = delta/2 : delta : (parms.n-1/2) * delta;

            xv = xv - offx;
            yv = yv - offy;

            [Xv(:,:,lev), Yv(:,:,lev)] = meshgrid( xv, yv );
 
        %--
    
        %--get velocity

        if lev == 1

            ind_s = 1;
            ind_e = get_velx_ind( parms.m-1, parms.n, lev, parms );
            velu = q( ind_s : ind_e ) / delta;

            velu(velu > cmax_w ) = cmax_w;
%             velu(velu < -cmax_w ) = -cmax_w;

            U(:,:,lev) = transpose( reshape( velu, parms.m-1, parms.n ) );                                  

            %store omega for coarse-grid interpolation
            velu_s = velu;
        else
            
            velu = zeros( (parms.m-1) * (parms.n), 1 );
            
            %--need some fancy indexing for what follows:
            
                %get indices of current grid (indices for overlapping
                %region are zero)
                indvelx = repmat(1 : parms.m-1,[1,parms.n]);
                indvely = repelem(1 : parms.n, parms.m-1);
                vel_ind = get_velx_ind(indvelx,indvely,lev,parms);
            
                %get all indices on gridlevel 1
                indvelx = repmat(1 : parms.m-1,[1,parms.n]);
                indvely = repelem(1 : parms.n, parms.m-1);
                vel_ind_f = get_velx_ind(indvelx,indvely,1,parms);
                
                %Get indices on gridlevel 1 just for overlapping part
                indvelx = repmat(2 : 2 : parms.m-2,[1,parms.n/2]);
                indvely = repelem(2 : 2 : parms.n, parms.m/2-1);
                vel_ind_s1 = get_velx_ind(indvelx,indvely, 1, parms);
                
                indvelx = repmat(2 : 2 : parms.m-2,[1,parms.n/2]);
                indvely = repelem(1 : 2 : parms.n-1, parms.m/2-1);
                vel_ind_s2 = get_velx_ind(indvelx,indvely, 1, parms);
                
                %indexing on the vector gamma starts at 1...
                n_sub = get_velx_ind( parms.m-1, parms.n, lev-1, parms);
                
                %non-overlapping indices:
                ind_no_over = ( (vel_ind - n_sub) > 0  );
                ind_no_over = vel_ind_f( ind_no_over ~= 0 );
                
                %overlapping indices:
                ind_over = ( (vel_ind - n_sub) < 0  ); 
                ind_over = vel_ind_f( ind_over ~= 0 );
            
            %--
            
            
            
            %Get part that doesn't overlap with finer grid
            velu( ind_no_over ) = q( vel_ind( vel_ind ~= 0 ) )/ delta;
            
            
%             %for overlapping part...
            velu( ind_over ) = 1/2* ( velu_s( vel_ind_s1 ) + ...
                velu_s( vel_ind_s2 ) );
%                         
            U(:,:,lev) = transpose( reshape( velu, parms.m-1, parms.n ) );
            
            %save omega for next coarse grid
            velu_s = velu;

        end

        %--



    end
    
    %plot them
    
    figure(1), clf

    for j = parms.mg : -1 : 1
        
        figure(1), hold on
        contourf(Xv(:,:,j), Yv(:,:,j), U(:,:,j), clevs, ...
            'edgecolor','none'); shading flat;

        colormap( cmap )
        axis equal
%         axis(range)
                
    end
    
    %plot body
    fill(soln.xb( 1 : parms.nb ), soln.xb( 1+parms.nb : 2*parms.nb ),'k'  )


   
