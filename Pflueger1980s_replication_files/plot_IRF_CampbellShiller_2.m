%Plot: 2x3 panel of impulse responses to Short-Term MP shock.
% Subplot 1: MP Rate
% Subplot 2: Output Gap
% Subplot 3: Equity Return (PD)
% Subplot 4: (Nominal) Bond Yield

function plot_IRF_CampbellShiller_2(asset1, asset2, figNameFile)
%%
        % Plot figure
        h = figure;
        pause(0.00001);
%         frame_h = get(handle(gcf),'JavaFrame');
%         set(frame_h,'Maximized',1);
        set(gcf,'WindowState','Maximized');
        set(gcf,'color','w');
        x_axis = [0:16];
             
       minaxis3=min([asset1.Irf3.shatsim(1:numel(x_axis))', asset2.Irf3.shatsim(1:numel(x_axis))']);
       maxaxis3=max([asset1.Irf3.shatsim(1:numel(x_axis))', asset2.Irf3.shatsim(1:numel(x_axis))']);
       
       minaxisynom1=min([asset1.Irf3.y10nom(1:numel(x_axis)), asset1.Irf3.y10nom_rn(1:numel(x_axis)), asset1.Irf3.y10nom_rp(1:numel(x_axis))]);
       minaxisynom2=min([asset2.Irf3.y10nom(1:numel(x_axis)), asset2.Irf3.y10nom_rn(1:numel(x_axis)), asset2.Irf3.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom1=max([asset1.Irf3.y10nom(1:numel(x_axis)), asset1.Irf3.y10nom_rn(1:numel(x_axis)), asset1.Irf3.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom2=max([asset2.Irf3.y10nom(1:numel(x_axis)), asset2.Irf3.y10nom_rn(1:numel(x_axis)), asset2.Irf3.y10nom_rp(1:numel(x_axis))]);
       yaxisy10nom3=[min(minaxisynom1,minaxisynom2),max(maxaxisynom1, maxaxisynom2)];
       
       minaxisynom1=min([asset1.Irf3.y10nom(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis)), asset1.Irf3.y10nom_rn(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis)), asset1.Irf3.y10nom_rp(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis))]);
       minaxisynom2=min([asset2.Irf3.y10nom(1:numel(x_axis))-asset2.Irf3.i(1:numel(x_axis)), asset2.Irf3.y10nom_rn(1:numel(x_axis))-asset2.Irf3.i(1:numel(x_axis)), asset2.Irf3.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom1=max([asset1.Irf3.y10nom(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis)), asset1.Irf3.y10nom_rn(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis)), asset1.Irf3.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom2=max([asset2.Irf3.y10nom(1:numel(x_axis))-asset2.Irf3.i(1:numel(x_axis)), asset2.Irf3.y10nom_rn(1:numel(x_axis))-asset2.Irf3.i(1:numel(x_axis)), asset2.Irf3.y10nom_rp(1:numel(x_axis))]);
        
       yaxisspread3=[min(minaxisynom1,minaxisynom2),max(maxaxisynom1, maxaxisynom2)];
       
       minaxis2=min([asset1.Irf2.shatsim(1:numel(x_axis))', asset2.Irf2.shatsim(1:numel(x_axis))']);
       maxaxis2=max([asset1.Irf2.shatsim(1:numel(x_axis))', asset2.Irf2.shatsim(1:numel(x_axis))']);
       
       minaxisynom1=min([asset1.Irf2.y10nom(1:numel(x_axis)), asset1.Irf2.y10nom_rn(1:numel(x_axis)), asset1.Irf2.y10nom_rp(1:numel(x_axis))]);
       minaxisynom2=min([asset2.Irf2.y10nom(1:numel(x_axis)), asset2.Irf2.y10nom_rn(1:numel(x_axis)), asset2.Irf2.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom1=max([asset1.Irf2.y10nom(1:numel(x_axis)), asset1.Irf2.y10nom_rn(1:numel(x_axis)), asset1.Irf2.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom2=max([asset2.Irf2.y10nom(1:numel(x_axis)), asset2.Irf2.y10nom_rn(1:numel(x_axis)), asset2.Irf2.y10nom_rp(1:numel(x_axis))]);
       yaxisy10nom2=[min(minaxisynom1,minaxisynom2),max(maxaxisynom1, maxaxisynom2)];
       
       minaxisynom1=min([asset1.Irf2.y10nom(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis)), asset1.Irf2.y10nom_rn(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis)), asset1.Irf2.y10nom_rp(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis))]);
       minaxisynom2=min([asset2.Irf2.y10nom(1:numel(x_axis))-asset2.Irf2.i(1:numel(x_axis)), asset2.Irf2.y10nom_rn(1:numel(x_axis))-asset2.Irf2.i(1:numel(x_axis)), asset2.Irf2.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom1=max([asset1.Irf2.y10nom(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis)), asset1.Irf2.y10nom_rn(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis)), asset1.Irf2.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom2=max([asset2.Irf2.y10nom(1:numel(x_axis))-asset2.Irf2.i(1:numel(x_axis)), asset2.Irf2.y10nom_rn(1:numel(x_axis))-asset2.Irf2.i(1:numel(x_axis)), asset2.Irf2.y10nom_rp(1:numel(x_axis))]);
        
       yaxisspread2=[min(minaxisynom1,minaxisynom2),max(maxaxisynom1, maxaxisynom2)];
       
       minaxis1=min([asset1.Irf1.shatsim(1:numel(x_axis))', asset2.Irf1.shatsim(1:numel(x_axis))']);
       maxaxis1=max([asset1.Irf1.shatsim(1:numel(x_axis))', asset2.Irf1.shatsim(1:numel(x_axis))']);
       
       minaxisynom1=min([asset1.Irf1.y10nom(1:numel(x_axis)), asset1.Irf1.y10nom_rn(1:numel(x_axis)), asset1.Irf1.y10nom_rp(1:numel(x_axis))]);
       minaxisynom2=min([asset2.Irf1.y10nom(1:numel(x_axis)), asset2.Irf1.y10nom_rn(1:numel(x_axis)), asset2.Irf1.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom1=max([asset1.Irf1.y10nom(1:numel(x_axis)), asset1.Irf1.y10nom_rn(1:numel(x_axis)), asset1.Irf1.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom2=max([asset2.Irf1.y10nom(1:numel(x_axis)), asset2.Irf1.y10nom_rn(1:numel(x_axis)), asset2.Irf1.y10nom_rp(1:numel(x_axis))]);
       yaxisy10nom1=[min(minaxisynom1,minaxisynom2),max(maxaxisynom1, maxaxisynom2)];
       
       minaxisynom1=min([asset1.Irf1.y10nom(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis)), asset1.Irf1.y10nom_rn(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis)), asset1.Irf1.y10nom_rp(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis))]);
       minaxisynom2=min([asset2.Irf1.y10nom(1:numel(x_axis))-asset2.Irf1.i(1:numel(x_axis)), asset2.Irf1.y10nom_rn(1:numel(x_axis))-asset2.Irf1.i(1:numel(x_axis)), asset2.Irf1.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom1=max([asset1.Irf1.y10nom(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis)), asset1.Irf1.y10nom_rn(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis)), asset1.Irf1.y10nom_rp(1:numel(x_axis))]);
       maxaxisynom2=max([asset2.Irf1.y10nom(1:numel(x_axis))-asset2.Irf1.i(1:numel(x_axis)), asset2.Irf1.y10nom_rn(1:numel(x_axis))-asset2.Irf1.i(1:numel(x_axis)), asset2.Irf1.y10nom_rp(1:numel(x_axis))]);
        
       yaxisspread1=[min(minaxisynom1,minaxisynom2),max(maxaxisynom1, maxaxisynom2)];
       
       
        subplot(2,3,1);
        p = plot(x_axis,asset1.Irf1.y10nom(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis)), 'linewidth',2);
        
        set(p,'Color', [0,0,0]);
        hold on 
        p = plot(x_axis,asset1.Irf1.y10nom_rn(1:numel(x_axis))-asset1.Irf1.i(1:numel(x_axis)), '--', 'linewidth', 2);
        set(p,'Color','red');
        p = plot(x_axis,asset1.Irf1.y10nom_rp(1:numel(x_axis)), ':', 'linewidth', 2);
        set(p,'Color','blue');
        plot(x_axis,0*asset1.Irf1.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([yaxisspread1(1),yaxisspread1(2)])
        ylim([-2,0.75])
        title('Demand Shock')
        hold off   
        
        
        subplot(2,3,2);
        p = plot(x_axis,asset1.Irf2.y10nom(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis)), 'linewidth',2);
        
        set(p,'Color', [0,0,0]);
        hold on 
        p = plot(x_axis,asset1.Irf2.y10nom_rn(1:numel(x_axis))-asset1.Irf2.i(1:numel(x_axis)), '--', 'linewidth', 2);
        set(p,'Color','red');
        p = plot(x_axis,asset1.Irf2.y10nom_rp(1:numel(x_axis)), ':', 'linewidth', 2);
        set(p,'Color','blue');
        plot(x_axis,0*asset1.Irf2.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([yaxisspread2(1),yaxisspread2(2)])
        ylim([-2,0.75])
        title({'Panel A: 1979.Q4-2001.Q1 Calibration','Supply Shock'}) 
        %title('Calibration: 1979.Q4-2001.Q1')
        hold off   
        
        subplot(2,3,3);
        p = plot(x_axis,asset1.Irf3.y10nom(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis)), 'linewidth',2);
        
        set(p,'Color', [0,0,0]);
        hold on 
        p = plot(x_axis,asset1.Irf3.y10nom_rn(1:numel(x_axis))-asset1.Irf3.i(1:numel(x_axis)), '--', 'linewidth', 2);
        set(p,'Color','red');
        p = plot(x_axis,asset1.Irf3.y10nom_rp(1:numel(x_axis)), ':', 'linewidth', 2);
        set(p,'Color','blue');
        plot(x_axis,0*asset1.Irf1.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([yaxisspread3(1),yaxisspread3(2)])
        ylim([-2,0.75])
        title('MP Shock')

        hold off   
        
        subplot(2,3,4);
        p = plot(x_axis,asset2.Irf1.y10nom(1:numel(x_axis))-asset2.Irf1.i(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on 
        p = plot(x_axis,asset2.Irf1.y10nom_rn(1:numel(x_axis))-asset2.Irf1.i(1:numel(x_axis)), '--', 'linewidth', 2);
        set(p,'Color','red');
        p = plot(x_axis,asset2.Irf1.y10nom_rp(1:numel(x_axis)), ':', 'linewidth', 2);
        set(p,'Color', 'blue');
        plot(x_axis,0*asset2.Irf1.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([yaxisspread1(1),yaxisspread1(2)])
        ylim([-2,0.75])
        title('Demand Shock')
        hold off
        
       
        subplot(2,3,5);
        p = plot(x_axis,asset2.Irf2.y10nom(1:numel(x_axis))-asset2.Irf2.i(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on 
        p = plot(x_axis,asset2.Irf2.y10nom_rn(1:numel(x_axis))-asset2.Irf2.i(1:numel(x_axis)), '--', 'linewidth', 2);
        set(p,'Color','red');
        p = plot(x_axis,asset2.Irf2.y10nom_rp(1:numel(x_axis)), ':', 'linewidth', 2);
        set(p,'Color', 'blue');
        plot(x_axis,0*asset2.Irf2.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([yaxisspread2(1),yaxisspread2(2)])
        ylim([-2,0.75])
        %ylabel('Yield Spread','fontweight','normal','fontsize',20);
        title({'Panel B: 2001.Q2-2019.Q4 Calibration','Supply Shock'}) 
        hold off 
        
        
       
        subplot(2,3,6);
        p = plot(x_axis,asset2.Irf3.y10nom(1:numel(x_axis))-asset2.Irf3.i(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on 
        p = plot(x_axis,asset2.Irf3.y10nom_rn(1:numel(x_axis))-asset2.Irf3.i(1:numel(x_axis)), '--', 'linewidth', 2);
        set(p,'Color','red');
        p = plot(x_axis,asset2.Irf3.y10nom_rp(1:numel(x_axis)), ':', 'linewidth', 2);
        set(p,'Color', 'blue');
        plot(x_axis,0*asset2.Irf3.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([yaxisspread3(1),yaxisspread3(2)])
        ylim([-2,0.75])
        legend({'Overall Yield Spread','Risk Neutral', 'Risk Premium'}, 'Location', 'southeast', 'FontSize', 14);

        title('MP Shock')
        hold off 
        
        saveas(gcf,char(figNameFile),'png')
        
end
