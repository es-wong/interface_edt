%this script aggregates and makes nice figs for the PbS-EDT manuscript,
%tentatively titled, 'Chemical Origin and Defect Dynamcis of Bulk and
%Interface Traps in PbS QD Solar Cells'

%set up directories

fig1Dir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\manuscript_figs\figures\AS_conductance';
fig2Dir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\manuscript_figs\figures\DLCP';
fig3Dir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\manuscript_figs\figures\DOS';
fig4Dir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_int. erface\analysis_work\EDT_schottky\manuscript_figs\figures\vacuum';
fig5Dir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\manuscript_figs\figures\JV';

thermalDir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\09_2017_m1\analysis\code';
vacDir = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\10-2017_m3\analysis\code';

%set up colormaps and symbols. blue is for interface states, green is for deep states. recall t_max = 17 for thermal data, i_max = 50 for
%vacuum data. i_sel = 6. f_max(thermal) = 5. f_max(vac) = 4. 

colSetBlue = makeColorMap([207.188 207.18 214.889]/255, [53 55 160]/255, 17); 
colSetGreen = makeColorMap([217.694 207.1875 207.1875]/255, [1 66 2]/255,17);

markzlowF = {'o', 'd', 's'};
markzhighF = {'*', '<'};

%%

%Figure 1--AS and conductance data

cd(thermalDir)
load('AS_JV_DLCP_2.mat')
%%
peakIdx = zeros(length(ASpeaks(1+offset_ref(1,1):offset_ref(1,2))),1);
w0freq = zeros(length(ASpeaks(1+offset_ref(1,1):offset_ref(1,2))),1);

%find indices of w0 in original frequency array to highlight frequency cutoff in TAS fig
for k = 1+offset_ref(1,1):offset_ref(1,2)
    
    freqHL = AS{k,1}(:,1) - ASpeaks(k,1);
    [minfreq,peakIdx(k,1)] = min(abs(freqHL));
    w0freq(k,1) = AS{k,1}(peakIdx(k,1),1);
    
    
    
end
 

%cap-freq fig with w0 highlights
for ii = 1:b_max 
    figure; 
    hold on;
    for k = 1+offset_ref(1,1):offset_ref(1,2)
        plot(AS{k,ii}(:,1),AS_norm{1,ii}(:,k),'LineWidth',3,'Color',colSetGreen(k,:));
        
        plot(AS{k,ii}(peakIdx(k,1),1),AS_norm{1,ii}(peakIdx(k,1),k),'Marker','o','MarkerFaceColor',colSetGreen(k,:),'MarkerEdgeColor','k','MarkerSize',20);
        
        
%         for m = 1:length(AS{k,ii}(:,1))
%             if AS{k,ii}(m,1) == w0freq(k,1)
%                 plot
%             else
%             end
%         end
    end
    
    axis square;
    box on;

    s2 = gca;
    box on
    pbaspect(s2, [1.5 1 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [2e2 5e6];
    s2.XTick = [1e3 1e4 1e5 1e6];
    s2.XTickLabel = {'10^{3}' '10^{4}' '10^{5}' '10^{6}'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '\omega (rad/s)';
    s2.XScale = 'log';

    s2.YLabel.String = 'C (nF cm^{-2})';
%     s2.YScale = 'log';
    s2.YLim = [0 800];
    s2.YTickLabel = {'0' '200' '400' '600' '800'};


    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('C_w0_',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact');
    
end

%%
%capacitane derivative

for ii = 1:b_max
        figure();%Plot the derivative data according to bias level to look for trends accross this parameter
        hold on;
        for k = 1+offset_ref(1,1):offset_ref(1,2)
            plot(AS{k,ii}(:,1),ASderiv_norm{ii,k}(:,2),'LineWidth',3, 'Color', colSetGreen(k,:)); %Plot the remaining temperatures and edit the plot accordingly
        end
        axis square
        box on;

        s2 = gca;
        box on
        pbaspect(s2, [1 2 1]);
        s2.LineWidth = 2;
        s2.FontSize = 44;
        s2.XLim = [3e3 3e6];
        s2.XTick = [1e4 1e5 1e6];
        s2.XTickLabel = {'10^{4}' '10^{5}' '10^{6}'};
        s2.FontName = 'Helvetica';
        s2.TickLength = [.02 .02];
        s2.XLabel.String = '\omega (rad/s)';
        s2.XScale = 'log';

        s2.YLabel.String = '\omega^{dC}/_{d \omega} (nF s cm^{-1})';
    %     s2.YScale = 'log';
        s2.YLim = [0 350];
        s2.YTick = [0 175 350];
        s2.YTickLabel = {'0' '175' '350'};


        set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
        frameName = [strcat('smoothed_derivaive_sub',num2str(B.(BNames{ii,1})),'mV')];
        savefig(gcf, fullfile(fig1Dir,frameName),'compact'); 
       
end

%%
%capacitance arrhenius

for ii = 1:b_max
    figure() 
    hold on
    g1 = errorbar(1e3*AS_inv{1,ii+b_max}(:,1),AS_inv{1,ii+b_max}(:,2),gauss_err(1+offset_ref(ii,1):offset_ref(ii,2),ii),'o','MarkerFaceColor',colSetGreen(12,:),'MarkerEdgeColor','k'); %create the scatter plot of ln(1/f) vs. 1/T

        set(g1,{'markers'},{20},{'Linewidth'},{1});

    plot(X_vals(1,:),Yhat(:,ii),'--r','LineWidth', 3); %plot the best fit line over it

    
    axis square
    box on

    s2 = gca;
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [3 7];
    s2.XTick = [3 5 7];
    s2.XTickLabel = {'3' '5' '7'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '1000/T^{-1} (K^{-1})';
%     s2.XScale = 'log';

    s2.YLabel.String = 'ln({\omega_0}T^{-2}) (rad s^{-1}K^{-2})';
%     s2.YScale = 'log';
    s2.YLim = [-3 3];
    s2.YTick = [-3 -1.5 0 1.5 3];
    s2.YTickLabel = {'-3' '-1.5' '0' '1.5' '3'};
    hold off

    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Arrhenius_sub',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact'); 
    
end

%%
%capacitance de-tuning

fermi = fittype('A./(1+exp(E.*x_fermi))','independent',{'x_fermi'},'coefficients',{'A','E'},'dependent',{'f_fermi'});


x_fermiT = 1e3./(T(1+offset_ref(1,1):12,1));
f_fermiT = dc_dom(1:8,1);

%approximate the fermi fit by first fitting a simple exponential
f6 = fit(x_fermiT,f_fermiT,'exp1');

figure
plot(f6,(1e3./(T(1+offset_ref(1,1):offset_ref(1,2),1))),dc_dom(:,1));

detune_expT = coeffvalues(f6); %parameter estimates
detune_err_expT = confint(f6,.6827); %one sigma confidence interval on the parameter estimates for the exponential approximation
detune_energy_errorT = detune_err_expT(1,2) - detune_expT(1,2); %calculate the error in the energy estimate

x_fermi_valsT = linspace(1e3./(T(offset_ref(1,1),1)),1e3./(T(offset_ref(1,2),1)));
f_fermi_vals_expT = detune_expT(1).*exp(x_fermi_valsT.*detune_expT(2));

for ii = 1 
    figure() 
    
    g1 = plot(1e3.*(1./(T(1+offset_ref(1,1):offset_ref(1,2),1))),dc_dom(:,ii),'o','color',colSetGreen(12,:),'MarkerEdgeColor','k','MarkerFaceColor',colSetGreen(12,:)); 
    set(g1,{'markers'},{20},{'Linewidth'},{1});
    
    hold on
    
    plot(x_fermi_valsT(1,:),f_fermi_vals_expT(1,:),'--r','LineWidth', 3); 

    box on

    s2 = gca;
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [3 6];
    s2.XTick = [3 4.5 6];
    s2.XTickLabel = {'3' '4.5' '6'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '1000/T (K^{-1})';
%     s2.XScale = 'log';

    s2.YLabel.String = '\omega ^{dC}/_{d\omega} (F)';
%     s2.YScale = 'log';
%     s2.YLim = [-3 3];
%     s2.YTick = [-3 -1.5 0 1.5 3];
%     s2.YTickLabel = {'-3' '-1.5' '0' '1.5' '3'};
    hold off

    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('amp_fermiFit_scale_',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact'); 
   
end

%%
%conductance figures
%find the index of the G-omega peaks


Gw0idx = zeros(length(Gpeaks_w0),1);


for k = 12:15
   Gw0tmp = AS{1,1}(:,1) - Gpeaks_w0(k,1);
   [val, Gw0idx(k,1)] = min(abs(Gw0tmp));
   
end
%%
%G/omega

for ii = 1 
    figure();
    hold on;
    for k = 8:15
%      for k = 
        plot(AS{k,ii}(:,1),(1e9*AS{k,ii}(:,8)),'LineWidth',3,'Color',colSetGreen(k,:)); 
        for j = 12:15
            p1 = plot(AS{j,ii}(Gw0idx(j,1),1),(1e9*AS{j,ii}(Gw0idx(j,1),8)),'Color',colSetGreen(j,:),'MarkerFaceColor',colSetGreen(j,:),'LineWidth',2,'Marker','d','MarkerEdgeColor','k','MarkerSize',20);
        end
    end
    axis square
    box on

    
    axis square;
    box on;

    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [1e5 1e6];
    s2.XTick = [1e5 1e6];
    s2.XTickLabel = {'10^{5}' '10^{6}'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '\omega (rad/s)';
    s2.XScale = 'log';

    s2.YLabel.String = '^{G}/_{\omega} (nSs^{-1})';
%     s2.YScale = 'log';
%     s2.YLim = [0 20];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};


    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('G-omegaSub1_',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact');
    
end

%%

%G-omega arrhenius plot

for ii = 1 
    figure
    
    g1 = errorbar(X_sigma(:,ii),Y_sigma(:,ii),sigma_0_err(12:15,ii),'d','MarkerFaceColor',colSetGreen(k,:));
    set(g1,{'markers'},{20},{'Linewidth'},{1},{'MarkerEdgeColor'},{'k'});
    
    hold on;
    plot(X_vals_sigma(1,:),Yhat_sigma(:,ii),'--r','LineWidth', 3); %plot the best fit line over it

    axis square;
    box on;
    
    
    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [3.3 4.0];
    s2.XTick = [3.3 3.65 4.0];
    s2.XTickLabel = {'3.3' '3.65' '4.0'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '1000/T (K^{-1})';
%     s2.XScale = 'log';

    s2.YLabel.String = 'ln(\omega_{0}T^{-2}) (rad s^{-1}K^{-2})';
%     s2.YScale = 'log';
%     s2.YLim = [0 20];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};


    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('G-omega_Arrhenius_',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact');
 
end

%%

%Gp-T

for ii = 1 
    figure();
    hold on;
    for k = 12:15
        
    h1 = plot(T(k,1),1e3*Gp_w0(1,k),'Marker','d','MarkerFaceColor',colSetGreen(k,:),'Color',colSetGreen(k,:),'MarkerEdgeColor','k'); 
        set(h1,{'MarkerSize'},{20});
        
    end
    axis square;
    box on;

    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [250 300];
    s2.XTick = [250 275 300];
    s2.XTickLabel = {'250' '275' '300'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = 'T (K)';
%     s2.XScale = 'log';

    s2.YLabel.String = 'G_{p} (mS)';
%     s2.YScale = 'log';
    s2.YLim = [2 8];
%     s2.YTickLabel = {'0' '200' '400' '600' '800'};


    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Gp_omega0-T_',num2str(B.(BNames{ii,1})),'mV')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact');
    
end

%%

%G/omega--interface (low frequency)

%I didn't fit these peaks to Gaussians (they actually look rather
%Lorentzian...). I read them off the plot. So the arrays below contain both
%the temperatures and the peak positions that I obtained manually


omIs = [540.4 603.2 929.9 1093 1288 1596]; %array for interface state resonant frequency
tIS = 220:10:270; %array for temperatures at which we meaure interface state responses


lowFidx = zeros(length(omIs),1);


for j = 1:length(omIs)
   lowFtmp = AS{1,1}(:,1) - omIs(1,j);
   [val, lowFidx(j,1)] = min(abs(lowFtmp));
   
end

for ii = 1 
    figure();
    hold on;
    for k = 1+offset_ref(1,1):offset_ref(1,2)
%      for k = 
        plot(AS{k,ii}(:,1),(1e9*AS{k,ii}(:,8)),'LineWidth',3,'Color',colSetBlue(k,:)); 
        for j = 1:length(omIs)
            p1 = plot(AS{j+7,ii}(lowFidx(j,1),1),(1e9*AS{j+7,ii}(lowFidx(j,1),8)),'Color',colSetBlue(j+7,:),'MarkerFaceColor',colSetBlue(j+7,:),'LineWidth',2,'Marker','o','MarkerEdgeColor','k','MarkerSize',20);
        end
    end
    axis square
    box on

    
    axis square;
    box on;

    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [4e2 1e4];
    s2.XTick = [4e2 1e4];
    s2.XTickLabel = {'10^{2}' '10^{3}' '10^{4}'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '\omega (rad/s)';
    s2.XScale = 'log';

    s2.YLabel.String = '^{G}/_{\omega} (nSs^{-1})';
%     s2.YScale = 'log';
    s2.YLim = [300 620];
    s2.YTick = [300 400 500 600];
    s2.YTickLabel = {'300' '400' '500' '600'};


    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('omegaIs')];
    savefig(gcf, fullfile(fig1Dir,frameName),'compact');
    
end


nSS = zeros(1,length(wTest)); %array to store densities of interface states

for j = 1:length(wTest)
    nSS(1,j) = (2.5/q) * AS{j+7,ii}(lowFidx(j,1),8);
end
%%
%nSS plot for SI

figure();
hold on;
for j = 1:length(omIs)
    p1 = plot(tIS(1,j),nSS(1,j),'Color',colSetBlue(j+7,:),'MarkerFaceColor',colSetBlue(j+7,:),'LineWidth',2,'Marker','o','MarkerEdgeColor','k','MarkerSize',20);
end
    
box on

axis square;
box on;

s2 = gca;
box on
pbaspect(s2, [1 2 1]);
s2.LineWidth = 2;
s2.FontSize = 44;
% s2.XLim = [4e2 1e4];
% s2.XTick = [4e2 1e4];
% s2.XTickLabel = {'10^{2}' '10^{3}' '10^{4}'};
s2.FontName = 'Helvetica';
s2.TickLength = [.02 .02];
s2.XLabel.String = 'T (K)';
% s2.XScale = 'log';

s2.YLabel.String = 'N_{SS}(cm^{-3}eV^{-1})';
% s2.YScale = 'log';
% s2.YLim = [300 620];
% s2.YTick = [300 400 500 600];
% s2.YTickLabel = {'300' '400' '500' '600'};

% 
% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('nss_SI')];
% savefig(gcf, fullfile(fig1Dir,frameName),'compact');





