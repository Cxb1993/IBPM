clear all, close all, clc

%---Build data matrix

    %snapshots for POD modes
    tvect = 35000 : 250 : 66750;

    %params for getdata...
    dir = './';
    lev = 3;
    pressure = 0;

    for j = 1 : length( tvect )

        ind = tvect( j );


        [xb,yb,codeb,xn,yn,un,vn,u0pn,v0pn,u0rn,v0rn,wn,sn,pn] = ...
            getdata(dir,ind,lev,pressure);

        n1 = length( wn(:,1) );
        n2 = length( wn(1,:) );
        nb = length( xb );
        
        nw = n1 * n2;
        n = nw + 2 * nb ;
        
        if j == 1
           
            Xbg = zeros( n, length( tvect ) );
            
        end

        wnvect = reshape( wn, [nw, 1] );

        bgvect = [ wnvect; xb; yb];

        Xbg(:, j ) = bgvect;

    end


%---

%---get DMD modes from data matrix

    %Rescale X for covariance matrix:
    X = Xbg(:, 1 : end - 1);
    X2 = Xbg(:, 2:end );
    
    [U, S, V] = svd( X, 'econ' );

    %Compute DMD (Phi are eigenvectors ) 
    r = 6; % truncate at r ;

    U = U(:, 1 : r);
    S = S(1 : r, 1 : r);
    V = V(:, 1 : r);
    
    Atilde = U' * X2 * V * inv(S);

    [W, eigs] = eig(Atilde);
    Phi = X2 * V * inv(S) * W;
          
%---          
          
%---postprocessing: plot POD modes
    load('cmap.mat')
    for j = 1 : 1 : r
        
        Uvect = real(Phi(:, j ));
        
        wnv = Uvect( 1 : nw) ;
        
        xb = linspace(0, 1, nb)'; yb = zeros(nb, 1);
        xbv = Uvect( nw + (1 : nb ) ) + xb;
        ybv = Uvect( nw + nb + (1 : nb) ) + yb;
        
        wnmat = reshape( wnv, [n1, n2] );
        
        figure(j)%, subplot(4,1,count)
        cmax_w = 0.02;
        clev_w = 20;
        clev = -cmax_w : 2 * cmax_w / clev_w : cmax_w;
        
        wnmat( wnmat > cmax_w ) = cmax_w;
        wnmat( wnmat < -cmax_w ) = -cmax_w;
        
        contourf(xn,yn,wnmat,clev,...
            'edgecolor','none'); shading flat;
        
        colormap(cmap)
        axis equal
        
        hold on
        
        %     contour(xn,yn,wn, cinc_w : cinc_w : cmax_w)
        caxis([-cmax_w cmax_w])
        
        %     title(['$t = ', num2str((it)*dt),'$'],'interpreter','latex')
        set(gca,'ticklabelinterpreter','latex')
        
        % Plot body number 1
        plot(xbv(find(codeb==1)),ybv(find(codeb==1)),'k','linewidth',2);
%         plot(xb(find(codeb==1)),yb(find(codeb==1)),'k','linewidth',2);

        box on
        
        
%         
%         figure(200),
%         plot(xbv, ybv)
%         
        
        xlabel('$x$','interpreter','latex','fontsize',16)
        ylabel('$y$','interpreter','latex','fontsize',16)
          
        
%         pause

        
        
        
    end
    
    figure(123456)
    theta = (0 : 1 : 100)*2 *pi/100;
    plot(cos(theta), sin(theta), 'k--');
    hold on, grid on
    scatter(real(diag(eigs)), imag(diag(eigs)), 'ok')
        

    omega = diag( log( eigs) / 0.25)

%---


          
          