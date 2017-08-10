function [] = plot_mode_combo( wnv, xbv, ybv, parms, fignum, ...
    figsavename, figtitle )


nb = parms.nb;
% load('cmap.mat')


[Xv, Yv, Omega] = get_mats_vort( parms, wnv );


cmax_w = min( max(max(max(Omega(:,:,1)))), abs(min(min(min(Omega(:,:,1)))) ) )/10
clev_w = 20;
clev = -cmax_w : 2 * cmax_w / clev_w : cmax_w;

cmap = build_cmap( length(clev) );

figure(fignum)

for jj = parms.mg : -1 : 1



    om_u = Omega(:,:,jj);   

    
    om_u( om_u > cmax_w ) = cmax_w;
    om_u( om_u < -cmax_w ) = -cmax_w;

    contourf(Xv(:,:,jj),Yv(:,:,jj),om_u,clev,...
        'edgecolor','none'); shading flat;

    colormap(cmap)
    axis equal
    set(gca,'ticklabelinterpreter','latex','fontsize',16)
    xlabel('$x$','interpreter','latex','fontsize',16)
    ylabel('$y$','interpreter','latex','fontsize',16)

    box on
    hold on


    plot(xbv + parms.xb0( 1 : nb ),ybv+ parms.xb0( nb + (1 : nb ) ),...
        'k','linewidth',3), hold on


    axis([-0.5 5 -3 3])
   

end


set(gca,'ticklabelinterpreter','latex','fontsize',30,'Xtick',[0 1 2 3 4 5],...
    'Ytick',[-3 -1.5 0 1.5 3])

ylabel('$y$','interpreter','latex','fontsize',32)
xlabel('$x$','interpreter','latex','fontsize',32)

title(figtitle,'interpreter','latex','fontsize',30)



axes('Position',[.25 .745 0.2 0.14]);
box on
nb = parms.nb;
plot(xbv + parms.xb0( 1 : nb ),ybv ...
    + parms.xb0( nb + (1 : nb ) ), 'k-','linewidth',2)
axis([0 1 -1, 1])
% set(gca, 'ticklabelinterpreter','latex','fontsize', 14, 'Ytick', [-0.04 0 0.04])
set(gca, 'ticklabelinterpreter','latex','Xtick', [], 'Ytick', [])


        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'Color', [1 1 1]); % Sets figure background
        set(gca, 'Color', [1 1 1]); % Sets axes background
        set(gcf, 'PaperUnits', 'centimeters' )
        set(gcf, 'PaperSize', [15 15] )
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [0 0 15 15] )
        set(gcf, 'PaperPosition', [0 0 15 15] )

dirsv = '~/Google Drive/Data processing/paper/Figures/';
% dirsv = '~/Desktop/';

print('-dpdf',[dirsv,figsavename,num2str(fignum),'.pdf'],'-r350')

