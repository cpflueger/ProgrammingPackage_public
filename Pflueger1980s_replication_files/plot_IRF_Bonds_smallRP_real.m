%Plot: 2x3 panel of impulse responses to Short-Term MP shock.
% Subplot 1: MP Rate
% Subplot 2: Output Gap
% Subplot 3: Equity Return (PD)
% Subplot 4: (Nominal) Bond Yield

function plot_IRF_Bonds_smallRP_real(asset1, asset2, figNameFile)
%%
        % Plot figure
        h = figure;
        pause(0.00001);
%         frame_h = get(handle(gcf),'JavaFrame');
%         set(frame_h,'Maximized',1);
        set(gcf,'WindowState','Maximized');
        set(gcf,'color','w');
        
        x_axis = [0:16];
             
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
        
        subplot(3,3,1);
        p = plot(x_axis,-asset1.Irf1.PD(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        p = plot(x_axis,-asset2.Irf1.PD(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf1.y10nom(1:numel(x_axis)),'-k');
        axis([0,16, -20,10])
        hold off
        ylabel('Stock Div. Yield','fontweight','normal','fontsize',20)
        title('Demand Shock','fontsize',24);
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
        
        
        subplot(3,3,4);
        p = plot(x_axis,asset1.Irf1.y10nom(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        axis([0,16, -1.7,1.5])
        p = plot(x_axis,asset2.Irf1.y10nom(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf1.y10nom(1:numel(x_axis)),'-k');
        hold off
        ylabel('Nom. Yield','fontweight','normal','fontsize',20)
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
        

        subplot(3,3,7);
        p = plot(x_axis,asset1.Irf1.y10real(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        axis([0,16, -0.3,0.3])
        p = plot(x_axis,asset2.Irf1.y10real(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf1.y10nom(1:numel(x_axis)),'-k');
        hold off
        ylabel('Real Yield','fontweight','normal','fontsize',20)
        set(gca, 'FontSize', 18); % to increase ticks font size on axes


        subplot(3,3,2);
        p = plot(x_axis,-asset1.Irf2.PD(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        p = plot(x_axis,-asset2.Irf2.PD(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf2.y10nom(1:numel(x_axis)),'-k');
        axis([0,16, -20,10])
        hold off
        title('Supply Shock','fontsize',24);
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
        
        
        subplot(3,3,5);
        p = plot(x_axis,asset1.Irf2.y10nom(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        axis([0,16, -1.7,1.5])
        p = plot(x_axis,asset2.Irf2.y10nom(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf2.y10nom(1:numel(x_axis)),'-k');
        hold off
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
 

        subplot(3,3,8);
        p = plot(x_axis,asset1.Irf2.y10real(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        %axis([0,16, -0.1,1])
        axis([0,16, -0.3,0.3])
        p = plot(x_axis,asset2.Irf2.y10real(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf2.y10nom(1:numel(x_axis)),'-k');
        hold off
        set(gca, 'FontSize', 18); % to increase ticks font size on axes


        subplot(3,3,3);
        p = plot(x_axis,-asset1.Irf3.PD(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        p = plot(x_axis,-asset2.Irf3.PD(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf2.y10nom(1:numel(x_axis)),'-k');
        axis([0,16, -20,10])
        hold off
        title('MP Shock','fontsize',24);
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
        
        subplot(3,3,6);
        p = plot(x_axis,asset1.Irf3.y10nom(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        axis([0,16, -1.7,1.5])
        p = plot(x_axis,asset2.Irf3.y10nom(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf3.y10nom(1:numel(x_axis)),'-k');
        hold off
        set(gca, 'FontSize', 18); % to increase ticks font size on axes
       

        subplot(3,3,9);
        p = plot(x_axis,asset1.Irf3.y10real(1:numel(x_axis)), '-', 'linewidth', 3);
        set(p,'Color', [0,0,0])
        hold on
        axis([0,16, -0.3,0.3])
        p = plot(x_axis,asset2.Irf3.y10real(1:numel(x_axis)), '--', 'linewidth', 3);
        plot(x_axis,0*asset2.Irf3.y10nom(1:numel(x_axis)),'-k');
        hold off
        legend('1980s Calibration','2000s Calibration','fontsize',14, 'Location','NorthEast')
        set(gca, 'FontSize', 18); % to increase ticks font size on axes

        saveas(gcf,char(figNameFile),'png')
        
end
