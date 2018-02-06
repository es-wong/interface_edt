%This script re-analyzes the DLCP data taking into account the energy de-tuning.

%plot Ndl vs <x>

for k = 1:t_max 
    figure()
    for ii = 1:v_max2 %index over voltages
        for iii = 1:5 %index over frequencies
            hold on
            h1 = plot(pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_greenSub(iii,:),'Marker', '.', 'LineWidth',3);
%             h2 = plot(IRpos_DLN{i,iii}(:,2),IRpos_DLN{i,iii}(:,3),'Color',colSet_green(iii,:),'Marker', 'd', 'LineWidth',3); %plot the filtered DLN vs. profile distance
                set(h1,{'markers'},{50});
%                 set(h2,{'markers'},{10})
        end
    end

    axis square;

    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
%     xlim([.1 .3]);
    
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
%     ylim([1e16 1e18]); 
    
    hold off 
    box on
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL_x_log_allF_',num2str(T(k,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     close(gcf)
end

%%
%This section of the code determines the energy of the applied AC frequency
%used in the DLCP measurement

E_DLCP_detune = cell(1,t_max); 
E_DLCP_err_detune = cell(1,t_max);

%calculate the energy and coresponding error
for k = 1:t_max
    for iii = 1:f_max
%         E_DLCP_detune{1,k}(1,iii) = abs(detune_exp(1,2)) + kB*T(k,1)*log(2*nu_0{1,1}(k,1)/(2*pi*freq(iii,1))); 
%         E_DLCP_err_detune{1,k}(1,iii) = detune_energy_error/detune_exp(1,2)*E_DLCP_detune{1,k}(1,iii); 
    
        %E_DLCP_detune{1,k}(1,iii) = abs(E_a_sigma(1,1)) + kB*T(k,1)*log(2*nu_0{1,1}(k,1)/(2*pi*freq(iii,1))); 
        %E_DLCP_err_detune{1,k}(1,iii) = sqrt( (E_a_sigma(1,2)/E_a_sigma(1,1))^2 + ((nu_0{1,1}(k,2)/nu_0{1,1}(k,1))^2))*E_DLCP_detune{1,k}(1,iii); 
  
        E_DLCP_detune{1,k}(1,iii) = abs(kB*T(k,1)*log(2*nu0tot(k,1)/(2*pi*freq(iii,1)))); 
        E_DLCP_err_detune{1,k}(1,iii) = nu0tot(k,2)/nu0tot(k,1) * E_DLCP_detune{1,k}(1,iii);
        
    end
end

%re-scale the N_{dl} vs. <x> plot to determine V_{dl} vs E_{\omega}
for k = 1+offset_ref(1,1):offset_ref(1,2)
    figure
    hold on
    for iii = 1:5
        errorbar(E_DLCP_detune{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err_detune{1,k}(1,iii),E_DLCP_err_detune{1,k}(1,iii),'Color',colSet_greenSub(iii,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_blue(t_max,:), 'MarkerSize',40);  
    end
    hold off
    axis square;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.15 .45]);
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    hold off
     
    set(gca,'linewidth',1.5);
    box on
%     
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL_E_log_allF_',num2str(T(k,1)),'K_detune')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end
%%
%make a total Ndl-E plot
markz = {'o','d','s', 'p', '>'};

figure
hold on
for k = 12:15
%     figure
    hold on
    for iii = 1:5
        errorbar(E_DLCP_detune{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err_detune{1,k}(1,iii),E_DLCP_err_detune{1,k}(1,iii),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',3,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgecolor','k');
    end
    hold off
    axis square;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.15 .45]);
    set(gca,'xtick',[.1 .2 .3 .4]);
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
%     hold off
    
    pbaspect([1 .6 1]);
    set(gca,'linewidth',1.5); 
    box on
end
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NDL_E_log_highT_detune')];
print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');
close(gcf)
hold off


%%
%re-plot the relevant data for the manuscript figures

DLCP_norm = cell(t_max,v_max2);
C_line_norm = cell(t_max,v_max2);

for k = 1:t_max 
    for ii = 1:v_max2 
        for iii = 1:f_max
            
            DLCP_norm{k,ii}{1,iii}(:,1) = DLCP{k,ii}{1,iii}(:,1);
            DLCP_norm{k,ii}{1,iii}(:,2) = (DLCP{k,ii}{1,iii}(:,2).*1e9)./A_cm;
            
            C_line_norm{k,ii}{1,iii}(:,1) = C_line{k,ii}{1,iii}(:,1);
            C_line_norm{k,ii}{1,iii}(:,2) = (C_line{k,ii}{1,iii}(:,2).*1e9)./A_cm;
        end
    end
end


for k = 1+offset_ref(1,1):offset_ref(1,2) 
    figure()
    for ii = 1:v_max2 
        hold on
        subplot(round(v_max2/2),2,ii); 
        for iii = 1:2
            hold on
            scatter(DLCP_norm{k,ii}{1,iii}(:,1),DLCP_norm{k,ii}{1,iii}(:,2), 200, colSet_greenSub(iii,:), 'filled'); %for each subplot, plot the test frequencies
            plot(C_line_norm{k,ii}{1,iii}(:,1),C_line_norm{k,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            t = title(strcat('V_{max} = -', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'FontName','Calibri')
            hold off
        end
        axis square;
        box on;
        x = xlabel ('V_{RMS} (V)');
        set(x,'FontName', 'Calibri');
        set(gca,'FontSize',16);
        y = ylabel ('C (nF/cm^{2})');
        set(y,'FontName','Calirbi');
        set(gca,'FontSize',16);
        hold off
    end
end

%%
%plot a sample figure consistent across all samples so that we can show an
%example of the experimental output in the text

T_DLCP_detune = find(T(:,1) == input('What temperature should be used as a sample? '));

% dc_DLCP = find(dc(1,:) == input('What is the magnitude of the maximum voltage that should be used as a sample? '));

%%

figure()
for k = T_DLCP_detune
    for ii = dc_DLCP
        hold on
        for iii = 1:5
            hold on
            
            h1 = plot(DLCP_norm{k,ii}{1,iii}(:,1),DLCP_norm{k,ii}{1,iii}(:,2),'Marker',markz{1,iii},'MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k');
                set(h1,{'MarkerSize'},{20});
                set(h1,{'LineWidth'},{1});
            
                plot(C_line_norm{k,ii}{1,iii}(:,1),C_line_norm{k,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            
        end
        axis square;
        box on;
        x = xlabel ('V_{RMS} (V)');
        set(x,'FontName', 'Calibri');
        set(gca,'FontSize',46);
        xlim([0 .3])
        set(gca,'xtick',[0 0.15 0.3])
         
        y = ylabel ('C (nF/cm^{2})');
        set(y,'FontName','Calirbi');
        set(gca,'FontSize',46);
        ylim([200 700]);
        set(gca,'ytick',[200 450 700])
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5); 
    hold off;
    
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('cap_amp_fit_samp_',num2str(B.(BNames{ii,1})),'mV_detuneAll')];
% %     print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% %     savefig(gcf, fullfile(figuresdir2A,frameName),'compact');
%     
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     
%     close(gcf)
    end
end
%%
%Make the sample N_{DL} vs. <x> plots

for k = T_DLCP_detune 
    figure()
    for ii = 1:v_max2 
        for iii = 1:5
            hold on
            
            h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii},'MarkerEdgeColor','k','MarkerFacecolor',colSet_green(k,:));
            set(h1,{'MarkerSize'},{20});
            set(h1,{'LineWidth'},{1});
      
        end
    end
    axis square;
     box on

    x = xlabel ('<x> (nm)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
%     xlim([32 36]);
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    ylim([1e17 2e18]);
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5); 
    hold off;
 

%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL_x_samp_',num2str(T(k,1)),'K_detune')];
% %     print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% %     savefig(gcf, fullfile(figuresdir2A,frameName),'compact');
%     
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     
%     close(gcf)
end

%%

%plot a distance profile of the carriers at different temperatures to
%determine where the response is localized within the film

% for iii = 1:
figure()
% for k = 12:15
for k = 1+offset_ref(1,1):offset_ref(1,2)
%     for iii = 1:f_max
        for iii = 1:3 %index over voltages
        hold on
        
        h2 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
            set(h2,{'markers'},{20})

        end
end
    axis square
    box on

    x = xlabel ('<x> (nm)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([20 100]);
%     set(gca,'xtick',[20  100])
        
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    ylim([1e16 1e20]); 
    set(gca,'ytick',[1e16 1e18 1e20])

    pbaspect([1 .6 1]);
    set(gca,'linewidth',1.5);
    hold off
% 
%     end

% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('N_DL-x_lowT_lowF_detune')];
% % print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% % savefig(gcf, fullfile(figuresdir2A,frameName),'compact');
% 
% print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2,frameName),'compact');
% 
% close(gcf)

%%
%Make carrier density vs. temperature plots

figure();
hold on;
for k = 12:15
    for iii = 1:5
    h1 = errorbar(T(k,1),DLN_tot{k,iii}(1,1),DLN_tot{k,ii}(1,2),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
%     h2 = errorbar(T(k,1),DLN_tot{k,2}(1,1),DLN_tot{k,2}(1,2),'Color',colSet_green(k,:),'Marker','d','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgeColor','k');
%     h3 = errorbar(T(k,1),DLN_tot{k,3}(1,1),DLN_tot{k,3}(1,2),'Color',colSet_green(k,:),'Marker','s','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgeColor','k');
    end
    axis square
    box on
    
    x = xlabel ('T (K)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([250 300])
    set(gca,'xtick',[250 275 300])
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    ylim([1e16 1e19])
    set(gca,'ytick',[1e15 1e17 1e18 1e19]);
    set(gca,'yscale','log');
    
end


pbaspect([1 1.6 1]);
set(gca,'linewidth',1.5);
hold off

% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('N_DL-T_highT_allF_3')];
% % % print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% % % savefig(gcf, fullfile(figuresdir2A,frameName),'compact'); 
% 
% print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2,frameName),'compact');
% 
% close(gcf)


%%
