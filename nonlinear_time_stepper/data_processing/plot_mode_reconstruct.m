function [] = plot_mode_reconstruct( wnv, xbv, ybv, parms, fignum, figname )


nb = parms.nb;
% load('cmap.mat')


[Xv, Yv, Omega] = get_mats_vort( parms, wnv );


cmax_w = 5;
clev_w = 20;
clev = -cmax_w : 2 * cmax_w / clev_w : cmax_w;

cmap = build_cmap( length(clev) );

figure(fignum)


ax1 = axes;

for jj = parms.mg : -1 : 1



    om_u = Omega(:,:,jj);   

    om_u( om_u > cmax_w ) = cmax_w;
    om_u( om_u < -cmax_w ) = -cmax_w;

    contourf(ax1,Xv(:,:,jj),Yv(:,:,jj),om_u,clev,...
        'edgecolor','none'); shading flat;

    colormap(ax1, cmap)
    axis equal
    set(gca,'ticklabelinterpreter','latex','fontsize',16)
    xlabel('$x$','interpreter','latex','fontsize',16)
    ylabel('$y$','interpreter','latex','fontsize',16)

    box on
    hold on

    axis(ax1,[-0.5 5 -3 3])
    

end


 plot(xbv + parms.xb0( 1 : nb ),ybv+ parms.xb0( nb + (1 : nb ) ),...
        'k','linewidth',2)


set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'Color', [1 1 1]); % Sets figure background
set(gca, 'Color', [1 1 1]); % Sets axes background
set(gcf, 'PaperUnits', 'centimeters' )
set(gcf, 'PaperSize', [15 15] )
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15] )
set( gcf, 'PaperPosition', [0 0 15 15] )


set(gca,'ticklabelinterpreter','latex',...
    'fontsize',26)


ylabel('$y$','interpreter','latex','fontsize',30)

set(gca,'Ytick',[-3 -1.5 0 1.5 3])
set(gca,'Xtick',[0 1 2 3 4 5])
xlabel('$x$','interpreter','latex','fontsize',30)


dirsv = '~/Google Drive/Data processing/paper/Figures/';
print('-dpdf',[dirsv,figname,num2str(fignum),'.pdf'],'-r300')

