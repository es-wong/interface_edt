%create nice looking DLCP figures

figuresdir2A = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\09_2017_m1\analysis\figures\DLCP\new\nice';

%Ndl-x figure

f1 = figure;


%add first plot
s1 = subplot_tight(2,1,1,.04);
% s1 = subplot(2,1,1);  
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for iii = 1:3 
        hold on
        
        h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
            set(h1,{'markers'},{20})

    end
end

box on
% pbaspect(s1, [2 1 1]);
s1.LineWidth = 1.5;
s1.FontSize = 30;
s1.XLim = [20 100];
s1.XTick = [20 40 60 80 100];
s1.FontName = 'Helvetica';
s1.XTickLabel =[];
s1.TickLength = [.02 .02];

s1.YLabel.String = 'N_{DL} (cm^{-3})';
s1.YScale = 'log';
s1.YLim = [1e16 1e20];
s1.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};


% get(gca,'OuterPosition')
% 
% v = get(gca,'Position');
% set(gca,'Position',[v(1) v(2)*4 v(3) v(4)])

%add second plot in 2x1 grid
s2 = subplot_tight(2,1,2,.04);
% s2 = subplot_tight(2,1,2); 
    for k = 12:15
        for iii = 4:5 
            hold on
        
            h2 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
                set(h2,{'markers'},{20})

        end
    end
    
box on
% pbaspect(s2, [2 1 1]);
s2.LineWidth = 1.5;
s2.FontSize = 30;
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

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% s1.OuterPosition(3) = ;
% s2.OuterPosition(3) = 

% frameName = [strcat('test')];

% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

% close(gcf)