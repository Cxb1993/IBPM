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
           
            X = zeros( n, length( tvect ) );
            
        end

        wnvect = reshape( wn, [nw, 1] );

        bgvect = [ wnvect; xb; yb];

        X(:, j ) = bgvect;

    end


%---

%---get POD modes from data matrix

    %Rescale X for covariance matrix:
    X = 1/sqrt( length(tvect) - 1 ) * X;
    XT = X'; %Store X'
    
    %Subtract mean:
    mean_X = mean( X, 2);
    for j = 1 : length(tvect)
        X(:,j) = X(:,j) ;
    end
        
    %Find nonzero singular value and right singular vector
    [V, Sig_sq] = eig( X' * X ); 
    
    Sig = sqrt( Sig_sq );
    
%     figure(1000)
    loglog( diag( Sig ), 'x' )
    
    VTXT = V' * XT;
    
    %Get POD modes U:
    U = zeros(n, length( tvect ) );
    for j = 1 : length( tvect )
        U(:, j) = (VTXT(j,:) ./ Sig(j,j) ).' ; 
    end
    
%---

%---postprocessing: plot POD modes
    load('cmap.mat')
    for j = length(tvect) : -1 : 1
        
        Uvect = U(:, j );
        
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
        
        box on
%         axis([-2 4.5 -1 1])
        
        
        
        figure(200),
        plot(xbv, ybv)
        
        pause
        
        xlabel('$x$','interpreter','latex','fontsize',16)
        ylabel('$y$','interpreter','latex','fontsize',16)
                
    end


%---


