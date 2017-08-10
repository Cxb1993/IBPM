function [] = plot_mode( wnv, xbv, ybv, parms, fignum, figname )


nb = parms.nb;
load('cmap.mat')


[Xv, Yv, Omega] = get_mats_vort( parms, wnv );


cmax_w = min( max(max(max(Omega(:,:,1)))), abs(min(min(min(Omega(:,:,1)))) ) )/50;

figure(fignum)


ax1 = axes;

for jj = parms.mg : -1 : 1



    om_u = Omega(:,:,jj);   

    clev_w = 20;
    om_u( om_u > cmax_w ) = cmax_w;
    om_u( om_u < -cmax_w ) = -cmax_w;


    clev = -cmax_w : 2 * cmax_w / clev_w : cmax_w;

    contourf(ax1,Xv(:,:,jj),Yv(:,:,jj),om_u,clev,...
        'edgecolor','none'); shading flat;

    colormap(ax1, cmap)
    axis equal
    set(gca,'ticklabelinterpreter','latex','fontsize',16)
    xlabel('$x$','interpreter','latex','fontsize',16)
    ylabel('$y$','interpreter','latex','fontsize',16)

    box on
    hold on


    plot(xbv + parms.xb0( 1 : nb ),ybv+ parms.xb0( nb + (1 : nb ) ),...
        'k','linewidth',2)


    axis(ax1,[-0.5 4 -1 1])
    
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'Color', [1 1 1]); % Sets figure background
        set(gca, 'Color', [1 1 1]); % Sets axes background
        set(gcf, 'PaperUnits', 'centimeters' )
        set(gcf, 'PaperSize', [15 12] )
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [0 0 15 12] )
        set(gcf, 'PaperPosition', [0 0 15 12] )
        
        
%         axis equal
%         axis(range);
        set(gca,'ticklabelinterpreter','latex',...
            'fontsize',24)
        
        
        ylabel('$y$','interpreter','latex','fontsize',24)
        
%         set(gca,'Xtick',[0 1 2 3 4 5], 'Ytick',[-2 -1 0 1 2])
        xlabel('$x$','interpreter','latex','fontsize',24)

%         dirsv = '~/Google Drive/figs/';

        dirsv = '~/Desktop/';



end

print('-dpng',[dirsv,figname,num2str(fignum),'.png'],'-r200')

