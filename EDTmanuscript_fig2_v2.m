%this script can be run only after manuscript fig 1 code has been run since
%it sets up the colormaps, etc. 

%this script makes the DLCP figures--figure 2. 

%green-->blue colormap
colSetGB = zeros(length(1+offset_ref(1,1):offset_ref(1,2)),3);

%colSetGB(1:6,:) = colSetGreen(1+offset_ref(1,1):10,:);
%colSetGB(7:end,:) = colSetBlue(11:15,:);

colSetGB(1:6,:) = colSetGreen(1:2:11,:);
colSetGB(7:end,:) = colSetBlue(11:15,:);

%blue-->green colormap
colSetBG = zeros(length(1+offset_ref(1,1):offset_ref(1,2)),3);

%colSetBG(1:6,:) = colSetBlue(1+offset_ref(1,1):10,:);
%colSetBG(7:end,:) = colSetGreen(11:15,:);

colSetBG(1:6,:) = colSetBlue(1:2:11,:);
colSetBG(7:end,:) = colSetGreen(11:15,:);

markzV2 = {'o' 'd'};

%%
%Ndl-x figure--high T/highF

figure()
hold on
for k = 12:15
    
    for iii = 1:3
        
        h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSetGreen(k,:),'Marker', markzV2{1,1}, 'LineWidth',1,'MarkerFaceColor',colSet_blue(k,:),'MarkerEdgeColor','k');    
        set(h1,{'markers'},{20})
        
    end
    
    for iii = 4:5 
      
        h2 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSetGreen(k,:),'Marker', markzV2{1,2}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
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

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlX_composite')];
savefig(gcf, fullfile(fig2Dir,frameName),'compact');

%%

%Ndl-E figure composite

figure;
hold on
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for iii = 1:5 
        if iii < 4
            h1 = errorbar(E_DLCP_detune{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err_detune{1,k}(1,iii),E_DLCP_err_detune{1,k}(1,iii),'Color',colSetGB(k-4,:),'Marker',markzV2{1,1},'LineWidth',3,'MarkerFaceColor',colSetGB(k-4,:), 'MarkerSize',20,'MarkerEdgecolor','k');
        elseif iii > 3
            h1 = errorbar(E_DLCP_detune{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err_detune{1,k}(1,iii),E_DLCP_err_detune{1,k}(1,iii),'Color',colSetBG(k-4,:),'Marker',markzV2{1,2},'LineWidth',3,'MarkerFaceColor',colSetBG(k-4,:), 'MarkerSize',20,'MarkerEdgecolor','k');
        end
    end
end

s1 = gca;
box on
pbaspect(s1, [2 1 1]);
s1.LineWidth = 2;
s1.FontSize = 44;
s1.XLim = [0.1 0.5];
s1.XTick = [0.1 0.2 0.3 0.4 0.5];
% s1.XTickLabel =[];
s1.FontName = 'Helvetica';
s1.TickLength = [.02 .02];
s1.XLabel.String = 'E_{\omega} (eV)';

s1.YLabel.String = 'N_{DL} (cm^{-3})';
s1.YScale = 'log';
s1.YLim = [1e16 1e20];
s1.YTickLabel = {'10^{16}' '10^{18}' '10^{20}'};


set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlE_composite')];
savefig(gcf, fullfile(fig2Dir,frameName),'compact');


%%

%Ndl-T line plots

figure;
hold on

for k = 10:14
    for iii = 1:4:5
        if iii == 1
            h1 = errorbar(T(k,1),DLN_tot{k,iii}(1,1),DLN_tot{k,ii}(1,2)/4,DLN_tot{k,ii}(1,2)/4,'Color',colSetBlue(k,:),'Marker',markzV2{1,1},'LineWidth',1,'MarkerFaceColor',colSetBlue(k,:),'MarkerSize',20,'MarkerEdgeColor','k');
        elseif iii == 5
            h1 = errorbar(T(k,1),DLN_tot{k,iii}(1,1),DLN_tot{k,ii}(1,2)/4,DLN_tot{k,ii}(1,2)/4,'Color',colSetGreen(k,:),'Marker',markzV2{1,2},'LineWidth',1,'MarkerFaceColor',colSetGreen(k,:),'MarkerSize',20,'MarkerEdgeColor','k');
        end
    end
end

s2 = gca;
box on
pbaspect(s2, [1 2 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [220 300];
s2.XTick = [230 265 300];
% s2.XTickLabel =;
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'T (K)';

s2.YLabel.String = 'N_{DL} (cm^{-3})';
s2.YScale = 'log';
s2.YLim = [1e17 1e19];
s2.YTickLabel = {'10^{17}' '10^{18}' '10^{19}'};

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NdlT_composite')];
savefig(gcf, fullfile(fig2Dir,frameName),'compact');



