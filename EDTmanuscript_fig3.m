%this script can be run only after manuscript fig 1 code has been run since
%it sets up the colormaps, etc. 

%this script makes the DOS figures--figure 3. 


%first, plot the traditional method

for ii = 1 
    figure()
    hold on;

    for k = 12:15 %index over temperatures. Note that we only go to the offeset point in specified in the AS script
        
        plot(abs(Nt_s_detune{1,ii}(:,k))/1e18,Ew_detune{1,ii}(:,k) , 'LineWidth', 3,'Color', colSet_green(k,:));

    end

    box on;

    s2 = gca;
    pbaspect(s2, [2 1 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [0 10];
    s2.XTick = [0 2 4 6 8 10];
    s2.XTickLabel =[];
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
%     s2.XLabel.String = 'N_{t} (\times 10^{18} cm^{-3}eV^{-1})';

    s2.YLabel.String = 'E_{\omega} (eV)';
%     s2.YScale = 'log';
    s2.YLim = [.26 .38];
    s2.YTick = [.26 .3  .34 .38];
%     s2.YTickLabel = {'0.26' '0.3' '0.34' '0.38'};
    
    hold off
    
%     
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_old_100kHZ_',num2str(Va(1,ii)*1000),'mV_detune_nolabel')];
    savefig(gcf, fullfile(fig3Dir,frameName),'compact');
    
end
%%

%next plot the new method

for ii = 1
    figure()
    hold on;
    for k = 12:15
        
        plot(abs(dos_vals_detune{k,ii}(:,1))/1e18,Ew_detune{1,ii}(:,k), 'LineWidth', 3,'Color', colSet_green(k,:));
    
    end
    
    box on;

    s2 = gca;
    pbaspect(s2, [2 1 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [0 10];
    s2.XTick = [0 2 4 6 8 10];
    % s2.XTickLabel =;
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = 'N_{t} (\times 10^{18} cm^{-3}eV^{-1})';

    s2.YLabel.String = 'E_{\omega} (eV)';
%     s2.YScale = 'log';
    s2.YLim = [.26 .38];
    s2.YTick = [.26 .3  .34 .38];
    s2.YTickLabel = {'0.26' '0.3' '0.34' '0.38'};
    
    hold off
    
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_new_100kHZ_',num2str(Va(1,ii)*1000),'mV_detune')];
    savefig(gcf, fullfile(fig3Dir,frameName),'compact');
end

%%

%make composite DOS fig

for ii = 1
    figure()
    hold on;
    for k = 12:15
        
        errorbar(ew_q,old_agg_detune(1,:),old_agg_detune(2,:),'Color',colSetGreen(15,:),'Marker','o','LineWidth', 2,'MarkerFaceColor',colSetGreen(15,:),'MarkerEdgeColor',colSet_green(1,:),'MarkerSize',20);
        errorbar(ew_q,new_agg_detune(1,:),new_agg_detune(2,:),'Color',colSetGreen(15,:),'Marker','o','LineWidth', 2,'MarkerFaceColor',colSetGreen(15,:),'MarkerEdgeColor','k','MarkerSize',20);
        
    end
    axis square;
    box on;

    s2 = gca;
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [.3 .4];
    s2.XTick = [0.3 0.35 0.4];
    s2.XTickLabel = {'0.3' '0.35' '0.4'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = 'E_{\omega} (eV)';

    s2.YLabel.String = 'N_{t} (cm^{-3}eV^{-1})';
    s2.YScale = 'log';
    s2.YLim = [1e16 1e19];
    s2.YTick = [1e16 1e17 1e18 1e19];
    s2.YTickLabel = {'10^{16}' '10^{17}' '10^{18}' '10^{19}'};
        
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_aggregate_tot_100kHz',num2str(Va(1,ii)*1000),'mV_detune')];
    savefig(gcf, fullfile(fig3Dir,frameName),'compact');
end
