%this script can be run only after manuscript fig 1 code has been run since
%it sets up the colormaps, etc. 

%this script makes the JV figures--figure 5.

%make solar sim curve
figure
hold on

h1 = plot(dark(:,2),(dark(:,4)),'Color',colSetGreen(vacT,:),'MarkerFaceColor',colSet_green(vacT,:),'LineWidth',2,'Marker','o','MarkerEdgeColor','k','MarkerSize',20);    
h2 = plot(light2(:,2),(light2(:,4)),'Color',colSetGreen(3,:),'MarkerFaceColor',colSet_green(3,:),'LineWidth',2,'Marker','o','MarkerEdgeColor','k','MarkerSize',20);

p1 = plot(light2(30,2),(light2(30,4)),'Color','r','MarkerFaceColor','r','LineWidth',2,'Marker','o','MarkerEdgeColor','k','MarkerSize',20);

% grid on
x1 = line([-.3 1.3], [0 0], 'LineWidth', 3, 'Color','k','LineStyle','--');
% y1 = line([0.6, 0.6], [-50, 150], 'LineWidth', 3, 'Color','r','LineStyle','--');
y1 = line([0.7, 0.7], [-50, 150], 'LineWidth', 3, 'Color','r','LineStyle','--');



s2 = gca;
box on
pbaspect(s2, [1 1 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [-.3 1.3];
s2.XTick = [-.3 0 .3 .6 .9 1.2];
s2.XTickLabel ={'-0.3' '0' '0.3' '0.6' '0.9' '1.2'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'V (volt)';

s2.YLabel.String = 'J (mA cm^{-2})';
% s2.YScale = 'log';
s2.YLim = [-50 150];
s2.YTick = [-50 0 50 100 150];
s2.YTickLabel = {'-50' '0' '50' '100' '150'};
 
hold off
% 
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('JV_solarSim2_gridOff_Eit')];
savefig(gcf, fullfile(fig5Dir,frameName),'compact');


%%
%JV-T

figure
hold on
for k = 2:t_max+1

    h2 = plot(JVD(:,1),abs(JVD(:,k)),'Color',colSetGreen(k-1,:),'LineWidth',2,'Marker','o','MarkerSize',20,'MarkerFaceColor',colSet_green(k-1,:),'MarkerEdgecolor','k');


end
hold off
box on

s2 = gca;
box on
pbaspect(s2, [1 1 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [-.3 .8];
s2.XTick = [-.2 0 .2 .4 .6];
s2.XTickLabel ={'-0.2' '0' '0.2' '0.4' '0.6'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'V (volt)';

s2.YLabel.String = 'J (mA cm^{-2})';
s2.YScale = 'log';
s2.YLim = [1e-3 1e3];
s2.YTick = [1e-3 1e-1 1e1 1e3];
s2.YTickLabel = {'10^{-3}' '10^{-1}' '10^{1}' '10^{3}'};
  
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('JV_temp_label')];
savefig(gcf, fullfile(fig5Dir,frameName),'compact');

%%
%make r0-T fig

figure()
hold on
for k = 12:15

    h1 = plot(T(k,1),1e3*r_ext(1,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');

end
hold off
box on

s2 = gca;
box on
pbaspect(s2, [1 2 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [250 300];
s2.XTick = [250 275 300];
s2.XTickLabel ={'250' '275' '300'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'T (K)';

s2.YLabel.String = 'R_{0} (\Omega cm^{-2})';
% s2.YScale = 'log';
% s2.YLim = [0.2 0.6];
% s2.YTick = [1e-3 1e-1 1e1 1e3];
% s2.YTickLabel = {'10^{-3}' '10^{-1}' '10^{1}' '10^{3}'};
  
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('r0-T_temp')];
savefig(gcf, fullfile(fig5Dir,frameName),'compact');

%%

var1 = 'jv';
var2 = 'colSet_green';

var3 = 'i_sel';
var4 = 'time_sel';
var5 = 'r_ext';
var6 = 'colSet_green_sel';

cd(vacDir)
T1 = load('AS_DLCP_vac_2',var1);
T2 = load('AS_DLCP_vac_2', var2);
T3 = load('AS_DLCP_vac_2',var3);
T4 = load('AS_DLCP_vac_2',var4);
T5 = load('AS_DLCP_vac_2',var5);
T6 = load('AS_DLCP_vac_2',var6);


jv = T1.(var1);
colSetGreen2 = T2.(var2);

i_sel = T3.(var3);
time_sel = T4.(var4);
vacRext = T5.(var5);
colSetGreen2A = T6.(var6);

%%

figure
hold on
h1 = plot(jv(:,1),abs(jv(:,2)),'Color',colSetGreen2(1,:),'Marker','o','MarkerFaceColor',colSetGreen2(1,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',20);
for k = 10:10:i_max
    h2 = plot(jv(:,1),abs(jv(:,k)),'Color',colSetGreen2(k,:),'Marker','o','MarkerFaceColor',colSetGreen2(k,:),'LineWidth',2,'MarkerEdgeColor','k','MarkerSize',20);
end
axis square
box on

s2 = gca;
box on
pbaspect(s2, [1 1 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [-.3 .8];
s2.XTick = [-.2 0 .2 .4 .6];
s2.XTickLabel ={'-0.2' '0' '0.2' '0.4' '0.6'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'V (volt)';

s2.YLabel.String = 'J (mA cm^{-2})';
s2.YScale = 'log';
s2.YLim = [1e-3 1e3];
s2.YTick = [1e-3 1e-1 1e1 1e3];
s2.YTickLabel = {'10^{-3}' '10^{-1}' '10^{1}' '10^{3}'};
  
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('JV_vacSel')];
savefig(gcf, fullfile(fig5Dir,frameName),'compact');


%%

figure()
hold on
for k = 1:i_sel-1

    h1 = plot(time_sel(1,k),1e3*vacRext(1,k),'Color',colSetGreen2A(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSetGreen2A(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');

end

axis square
box on

s2 = gca;
box on
pbaspect(s2, [1 2 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
s2.XLim = [-50 1050];
s2.XTick = [0 500 1000];
s2.XTickLabel ={'0' '500' '1000'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'time (min)';

s2.YLabel.String = 'R_{0} (\Omega cm^{-2})';
s2.YScale = 'log';
s2.YLim = [1e1 1e4];
% s2.YTick = [1e-3 1e-1 1e1 1e3];
% s2.YTickLabel = {'10^{-3}' '10^{-1}' '10^{1}' '10^{3}'};
  
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('r0-time_vacSel')];
savefig(gcf, fullfile(fig5Dir,frameName),'compact');

%%

%make j0-time figures

var1 = 'j0';

cd(vacDir)
T1 = load('AS_DLCP_vac_2',var1);

j0vac = T1.(var1);

figure;
hold on
    
for k = 1:i_sel-1
    
    h1 = plot(time_sel(1,k),j0vac(1,k),'Color',colSet_green_sel(i_sel-1,:),'Marker','o','LineWidth',1,'MarkerFaceColor',colSet_green_sel(i_sel-1,:), 'MarkerSize',20, 'MarkerEdgeColor','k');

end

axis square
box on

x = xlabel ('Time (min)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
% xlim([0 .5])
% set(gca,'xtick',[240 270 300 330])

y = ylabel ('J_{0} (mAcm^{-2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
% ylim([1e-3 1e2])
set(gca,'yscale','log');
% set(gca,'ytick',[1e-3 1e-1 1e1])

pbaspect([1 1.6 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('j0-time_vacSel')];
savefig(gcf, fullfile(fig5Dir,frameName),'compact');

%%
%make figure for j0 vs. T. fit to get phi_b

figure;
hold on
    
for k = 12:14
    
    h1 = plot(T(k,1),log(j0(1,k)),'Color',colSet_green(k,:),'Marker','o','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');

end

axis square
box on

x = xlabel ('T (K)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
% xlim([0 .5])
% set(gca,'xtick',[240 270 300 330])

y = ylabel ('J_{0} (mAcm^{-2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
% ylim([1e-3 1e2])
% set(gca,'yscale','log');
% set(gca,'ytick',[1e-3 1e-1 1e1])

pbaspect([1 1.6 1]);
set(gca,'linewidth',1.5);
hold off

% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('j0_time')];
% print(gcf, '-dpng', strcat(figuresdir3,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir3,frameName),'compact');


% print(gcf, '-dpng', strcat(figuresdir3A,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir3A,frameName),'compact');

% close(gcf)






