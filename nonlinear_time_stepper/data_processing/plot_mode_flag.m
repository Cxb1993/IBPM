function [] = plot_mode_flag( xbv, ybv, parms, fignum, figname )


        load('cmap.mat')

        figure(fignum+ 1000), hold on

        nb = parms.nb;
        plot(xbv + parms.xb0( 1 : nb ),ybv ...
            + parms.xb0( nb + (1 : nb ) ), 'k-','linewidth',2)
        
        axis([0 1 -0.1 0.1])
        set(gca,'ticklabelinterpreter','latex','fontsize',16)
        xlabel('$x$','interpreter','latex','fontsize',16)
        ylabel('$y$','interpreter','latex','fontsize',16)
        set(gca,'Xtick',[0 0.5 1], 'Ytick',[-0.1 0 0.1])    
            
        box on

        figure(fignum+1000)
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'Color', [1 1 1]); % Sets figure background
        set(gca, 'Color', [1 1 1]); % Sets axes background
        set(gcf, 'PaperUnits', 'centimeters' )
        set(gcf, 'PaperSize', [15 10.5] )
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [0 0 15 10.5] )
        set( gcf, 'PaperPosition', [0 0 15 10.5] )
        
          
        dirsv = '~/Google Drive/figs/';
        print('-dpng',[dirsv,figname,'_flag',num2str(fignum+1000),'.png'],'-r200')
        
        
%         pause
