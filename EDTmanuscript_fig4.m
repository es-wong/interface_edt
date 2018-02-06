%this script can be run only after manuscript fig 1 code has been run since
%it sets up the colormaps, etc. 

%this script makes the vacuum figures--figure 4.

var1 = 'AS';
var2 = 'AS_norm';

var3 = 'time_prime';
var4 = 'Ew';
var5 = 'EwErr';

var6 = 'time';
var7 = 'DLN_tot';

var8 = 'i_max';
var9 = 'f_max';

var10 = 'distAv';
var11 = 'NdlAv';
var12 = 'NdlAvErr';
var13 = 'distAvErr';


cd(vacDir)
T1 = load('AS_DLCP_vac_2',var1);
T2 = load('AS_DLCP_vac_2', var2);
T3 = load('AS_DLCP_vac_2',var3);
T4 = load('AS_DLCP_vac_2',var4);
T5 = load('AS_DLCP_vac_2',var5);
T6 = load('AS_DLCP_vac_2',var6);
T7 = load('AS_DLCP_vac_2',var7);
T8 = load('AS_DLCP_vac_2',var8);
T9 = load('AS_DLCP_vac_2',var9);
T10 = load('AS_DLCP_vac_2',var10);
T11 = load('AS_DLCP_vac_2',var11);
T12 = load('AS_DLCP_vac_2',var12);
T13 = load('AS_DLCP_vac_2',var13);

vacAS = T1.(var1);
vacAS_norm = T2.(var2);

time_prime = T3.(var3);
vacEw = T4.(var4);
vacEwErr = T5.(var5);

time = T6.(var6);
vacDLNtot = T7.(var7);

i_max = T8.(var8);
vacFmax = T9.(var9);

distAv = T10.(var10);
NdlAv = T11.(var11);
NdlAvErr = T12.(var12);
distAvErr = T13.(var13);

vacT = find(T(:,1) == 300);

%%

%cap-freq

for ii = 1:b_max %index across V_{bias}
    figure(); %Plot the capacitance data grouped according to bias level to look for trends across this parameter
    hold on;
    for i = 1:i_max %plot the remaining temperature curves and edit the plot accordingly
        plot(vacAS{i,ii}(:,1),vacAS_norm{1,ii}(:,i),'LineWidth',2,'Color',colSetGreen(vacT,:));  
    end
    axis square;
    box on;

    s2 = gca;
    box on
    pbaspect(s2, [1 1 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [6e2 1e6];
    s2.XTick = [1e3 1e4 1e5 1e6];
    s2.XTickLabel = {'10^{3}' '10^{4}' '10^{5}' '10^{6}'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '\omega (rad/s)';
    s2.XScale = 'log';

    s2.YLabel.String = 'C (nF cm^{-2})';
%     s2.YScale = 'log';
    s2.YLim = [0 650];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};


    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('C_vac_',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig4Dir,frameName),'compact');
end

%%

%w0-energy
%     
% for ii = 1:3
    figure() %create the Arrhenius plot
    for ii = 1
    hold on
        
    g1 = errorbar(time_prime(1,:),vacEw(:,ii),(vacEwErr(:,ii)./2),'Marker',markzlowF{1,ii},'color',colSetGreen(vacT,:),'MarkerFacecolor',colSetGreen(vacT,:),'MarkerEdgeColor','k'); 
    set(g1,{'markers'},{20},{'Linewidth'},{1});
    
    end
    
    box on;
    
    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [0 1e3];
%     s2.XTick = [0 5e2 1e3];
    s2.XTickLabel = {'0' '500' '1000'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = 't (min)';
%     s2.XScale = 'log';

    s2.YLabel.String = 'E_{\omega_{0}} (eV)';
    s2.YScale = 'log';
    s2.YLim = [.35 .45];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};

    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('energy_time_totSI')];
    savefig(gcf, fullfile(fig4Dir,frameName),'compact');

%%

figure()
hold on;
for i = 1:i_max
    
    h2 = errorbar(time(1,i),vacDLNtot{i,4}(1,1),vacDLNtot{i,4}(1,2),'Color',colSetGreen(vacT,:),'Marker','d','LineWidth',1,'MarkerFaceColor',colSetGreen(vacT,:),'MarkerSize',20,'MarkerEdgeColor','k');
    
end

box on

 box on;
    
    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [-50 1050];
%     s2.XTick = [0 5e2 1e3];
    s2.XTickLabel = {'0' '500' '1000'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = 't (min)';
%     s2.XScale = 'log';

    s2.YLabel.String = 'N_{DL} (cm^{-3})';
    s2.YScale = 'log';
    s2.YLim = [1e17 3e18];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};

    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL-time_log_',num2str(freq(4,1)),'Hz')];
    savefig(gcf, fullfile(fig4Dir,frameName),'compact');
% 
%%
%Ndl-X at 7kHz (for SI)


figure()
hold on
for k = 1:i_max %index over temperatures
    for iii = vacFmax %index over frequencies

        h1 = errorbar(1e3*distAv(k,iii),NdlAv(k,iii),NdlAvErr(k,iii),NdlAvErr(k,iii),distAvErr(k,iii),distAvErr(k,iii),'Marker', 'o','LineWidth',2,'color',colSetGreen(vacT,:),'MarkerFaceColor',colSetGreen(vacT,:),'MarkerEdgeColor','k');
        set(h1,{'markers'},{20});

    end
end

box on

s2 = gca;
box on
pbaspect(s2, [1 2 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [40 60];
% s2.XTick = [0 5e2 1e3];
% s2.XTickLabel = {'0' '500' '1000'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = '<x> (nm)';
% s2.XScale = 'log';

s2.YLabel.String = 'N_{DL} (cm^{-3})';
s2.YScale = 'log';
s2.YLim = [1e17 2e18];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};
% 

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NDL-X_log_7kHz')];
savefig(gcf, fullfile(fig4Dir,frameName),'compact');



