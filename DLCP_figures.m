%create nice looking DLCP figures

figuresdir2A = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\09_2017_m1\analysis\figures\DLCP\new\nice';

%Ndl-x figure

figure;
hold on
for k = 12:15
    for iii = 4:5 
      

        h2 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
            set(h2,{'markers'},{20})

    end
end

s2 = gca;
box on
pbaspect(s2, [2 1 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [20 100];
s2.XTick = [20 40 60 80 100];
% s2.XTickLabel =;
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = '<x> (nm)';

s2.YLabel.String = 'N_{DL} (cm^{-3})';
s2.YScale = 'log';
s2.YLim = [1e16 1e20];
s2.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};

size = s2.OuterPosition;

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlX_highT_highF')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

% close(gcf)


%%

figure; 
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for iii = 1:3 
        hold on
        
        h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
            set(h1,{'markers'},{20})
        
%         s1 = get(gca);
    end
end

s1 = gca;
box on
pbaspect(s1, [2 1 1]);
s1.LineWidth = 2;
s1.FontSize = 44;
s1.XLim = [20 100];
s1.XTick = [20 40 60 80 100];
s1.FontName = 'Helvetica';
s1.XTickLabel =[];
s1.TickLength = [.02 .02];

s1.YLabel.String = 'N_{DL} (cm^{-3})';
s1.YScale = 'log';
s1.YLim = [1e16 1e20];
s1.YTick = [1e16 1e18 1e20];
s1.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlX_lowT_lowF')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

% close(gcf)


%%

%Ndl-E figure

figure;
hold on
for k = 12:15
    for iii = 1:5 
  
        errorbar(E_DLCP_detune{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err_detune{1,k}(1,iii),E_DLCP_err_detune{1,k}(1,iii),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',3,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgecolor','k');
    
    end
end

s2 = gca;
box on
pbaspect(s2, [2 1 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [0.1 0.5];
s2.XTick = [0.2 0.3 0.4];
% s2.XTickLabel =;
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'E_{\omega} (eV)';

s2.YLabel.String = 'N_{DL} (cm^{-3})';
s2.YScale = 'log';
s2.YLim = [1e16 1e20];
s2.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};

size = s2.OuterPosition;

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlE_highT_allF')];
% % print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2A,frameName),'compact');
% 
% % close(gcf)


%%

figure;
hold on
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for iii = 1:5
  
        errorbar(E_DLCP_detune{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err_detune{1,k}(1,iii),E_DLCP_err_detune{1,k}(1,iii),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',3,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgecolor','k');
    
    end
end

s1 = gca;
box on
pbaspect(s1, [2 1 1]);
s1.LineWidth = 2;
s1.FontSize = 44;
s1.XLim = [0.1 0.5];
s1.XTick = [0.2 0.3 0.4];
s1.XTickLabel =[];
s1.FontName = 'Helvetica';
s1.TickLength = [.02 .02];
% s1.XLabel.String = 'E_{\omega} (eV)';

s1.YLabel.String = 'N_{DL} (cm^{-3})';
s1.YScale = 'log';
s1.YLim = [1e16 1e20];
s1.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};

% size = s2.OuterPosition;

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlE_allT_allF')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

% close(gcf)

%%

%Ndl-T

figure;
hold on

for k = 12:15
    for iii = 1:5
        
    h1 = errorbar(T(k,1),DLN_tot{k,iii}(1,1),DLN_tot{k,ii}(1,2),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
    
    end
end

s2 = gca;
box on
pbaspect(s2, [2 1 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [180 300];
s2.XTick = [180 220 260 300];
% s2.XTickLabel =;
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'E_{\omega} (eV)';

s2.YLabel.String = 'N_{DL} (cm^{-3})';
s2.YScale = 'log';
s2.YLim = [1e16 1e20];
s2.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};

size = s2.OuterPosition;

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlT_highT_allF')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

% close(gcf)

%%


figure;
hold on

for k = 1+offset_ref(1,1):offset_ref(1,2)
    for iii = 1:3
        
    h1 = errorbar(T(k,1),DLN_tot{k,iii}(1,1),DLN_tot{k,ii}(1,2),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
    
    end
end

s1 = gca;
box on
pbaspect(s1, [2 1 1]);
s1.LineWidth = 2;
s1.FontSize = 44;
s1.XLim = [180 300];
s1.XTick = [180 220 260 300];
s1.XTickLabel = [];
s1.FontName = 'Helvetica';
s1.TickLength = [.02 .02];
% s1.XLabel.String = 'E_{\omega} (eV)';

s1.YLabel.String = 'N_{DL} (cm^{-3})';
s1.YScale = 'log';
s1.YLim = [1e16 1e20];
s1.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};

% size = s2.OuterPosition;

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlT_allT_lowF')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

% close(gcf)

