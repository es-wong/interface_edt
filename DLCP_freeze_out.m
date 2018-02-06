%re-calculate the DLCP data at the closest frequency to omega_{d} using the hand-calculated equivalent dielectric

epsilon_equiv = 11.5; %input calculated equivalent epsilon

N_dl_equiv = cell(t_max,v_max); %holds the data for the carrier density


% for k = 1:offset_ref(1,2) %index over temperatures
for k = 1:t_max
    for ii = 1:v_max2 %index over max-voltages
        for iii =  omega_t %index over frequencies
            N_dl_equiv{k,ii}(1,omega_idx(1,k)) = (-C_fits{k,ii}{1,omega_idx(1,k)}(1,3)^3/(2*q*epsilon_equiv*eps_0*(A^2)*C_fits{k,ii}{1,omega_idx(1,k)}(1,2))*1e-6); %re-calculate calculate ionized gap density in cm^{-3}
        end
    end
end

pos_Ndl_equiv = cell(t_max,v_max); %create a cell to hold the positive values of N_dl

% for k = 1:offset_ref(1,2) %index over temperatures
for k = 1:t_max
    for ii = 1:v_max2 %index over maximum voltages
        testN = N_dl_equiv{k,ii}; %extract the N_dl values into a test array
        testN(testN <0) = 0; %assign zeros to the negative values in the test array
        pos_Ndl_equiv{k,ii} = testN; %re-import the test array into the positive cell
    end
end

DLN_equiv = cell(t_max,1); %initialize a cell that will hold the data to plot drive level density vs. profile distance
pos_DLN_equiv = cell(t_max,1); %initialize a cell that will hold the cleaned drive level density vs. profile distance 
DLN_test_equiv = zeros(v_max2,3); %initialize a test array on which  to run the cleaning functions below


% for k = 1:offset_ref(1,2) %index over temperatures
for k = 1:t_max
    for ii = 1:v_max2 %index over V_{max}
%         for iii = omega_t %index over frequencies
            DLN_equiv{k,omega_idx(1,k)}(:,1) = -dc'; %write the V_{max} to the array for reference
            DLN_equiv{k,omega_idx(1,k)}(ii,2) = ((epsilon_equiv*eps_0*A)/C_fits{k,ii}{1,omega_idx(1,k)}(1,3))*1e6; %calculate the drive level <x> in units of microns 
            DLN_equiv{k,omega_idx(1,k)}(ii,3) = N_dl_equiv{k,ii}(1,omega_idx(1,k)); %import the N_dl data into the array
%         end
    end
end


% for k =1:offset_ref(1,2) %index over temperature
for k =1:t_max
    for ii = 1:v_max2 %index over voltage
%         for iii = omega_t %index over frequency
            DLN_test_equiv(:,1:3) = DLN_equiv{k,omega_idx(1,k)}(:,1:3); %import the unfiltered DLN sub arrays into a test array 
            DLN_pos_equiv = DLN_test_equiv(DLN_test_equiv(:,3) >= 0,:); %import only the data that contains positive carrier values in column 3 into a new test array DLN pos
            pos_DLN_equiv{k,omega_idx(1,k)} = DLN_pos_equiv; %re-import the test array into the new cell, filtering the data
%         end
    end
end

% for k =1:offset_ref(1,2) %index over temperature
for k = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over voltages
        for iii = 1:2%index over frequencies
            hold on
            h1 = plot(pos_DLN_equiv{k,omega_idx(1,k)}(:,2),pos_DLN_equiv{k,omega_idx(1,k)}(:,3),'Color',colSet_green(omega_idx(1,k),:),'Marker', '.', 'LineWidth',3);
            h2 = plot(pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(iii,:),'Marker', '.', 'LineWidth',3); %plot the filtered DLN vs. profile distance
                set(h1,{'markers'},{50});
                set(h2,{'markers'},{50})
                if iii >= free_idx(k,1)
                    break
                end
        end
    end
    axis square;
%     t = title(strcat('$N_{DL}$ vs. Distance,', num2str(T(i,1)), 'K'), 'FontSize',40); %Title the graph
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.1 .35]);
%     ylim([1e16 1e18]); 
    box on
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL_x_re-calculated_',num2str(T(k,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end

%%
colSet_bluered = makeColorMap([53 55 160]/255, [1 0 0], t_max);

figure()
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for ii = 1:v_max2 %index over voltages
        hold on
        h2 = plot(pos_DLN{k,free_idx(k,1)}(:,2),pos_DLN{k,free_idx(k,1)}(:,3),'Color',colSet_bluered(k,:),'Marker', '.', 'LineWidth',3); %plot the filtered DLN vs. profile distance
%         h2 = plot(pos_DLN{k,omega_idx(1,k)}(:,2),pos_DLN{k,omega_idx(1,k)}(:,3),'Color',colSet_bluered(k,:),'Marker', '.', 'LineWidth',3);    
    set(h2,{'markers'},{50})
    end
    axis square;
%     t = title(strcat('$N_{DL}$ vs. Distance,', num2str(T(i,1)), 'K'), 'FontSize',40); %Title the graph
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.09 .16]);
%     ylim([1e16 1e19]); 
    box on
end
% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('N_free-T_spatial_TAS-range')];
% print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2,frameName),'compact');
% close(gcf)


%%

figure()
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for ii = 1:v_max2 %index over voltages
        hold on
%         h2 = plot(pos_DLN{k,free_idx(k,1)}(:,2),pos_DLN{k,free_idx(k,1)}(:,3),'Color',colSet_bluered(k,:),'Marker', '.', 'LineWidth',3); %plot the filtered DLN vs. profile distance
        h2 = plot(pos_DLN{k,omega_idx(1,k)}(:,2),pos_DLN{k,omega_idx(1,k)}(:,3),'Color',colSet_bluered(k,:),'Marker', 'd', 'LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:));    
    set(h2,{'markers'},{20})
    end
    axis square;
%     t = title(strcat('$N_{DL}$ vs. Distance,', num2str(T(i,1)), 'K'), 'FontSize',40); %Title the graph
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
    xlim([.09 .17]);
    ylim([1e16 1e19]); 
    box on
end
% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('N_omega-T_spatial_TAS-range')];
% print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2,frameName),'compact');
% close(gcf)

%%

figure();
hold on;
for k = 1 + offset_ref(1,1):offset_ref(1,2)
%     h1 = errorbar(T(k,1),E_DLCP{1,k}(1,free_idx(k,1)),E_DLCP_err{1,k}(1,free_idx(k,1)),'Color',colSet_bluered(k,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',40);
%     h2 = errorbar(T(k,1),E_DLCP{1,k}(1,omega_idx(1,k)),E_DLCP_err{1,k}(1,omega_idx(1,k)),'Color',colSet_bluered(k,:),'Marker','d','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',20);
    h1 = plot(T(k,1),E_DLCP{1,k}(1,free_idx(k,1)),'Color',colSet_bluered(k,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',40);
%     h2 = plot(T(k,1),E_DLCP{1,k}(1,omega_idx(1,k)),'Color',colSet_bluered(k,:),'Marker','d','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',20);
    axis square;
    x = xlabel ('T (K)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    y = ylabel ('E_{\omega} (eV)');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
%     set(gca,'yscale','log');
end
hold off
box on
% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('E_omega-T_free-omega')];
% print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2,frameName),'compact');
% close(gcf)

%%


figure();
hold on;
for k = 1 + offset_ref(1,1):offset_ref(1,2)
    h1 = errorbar(T(k,1),DLN_tot{k,free_idx(k,1)}(1,1),DLN_tot{k,free_idx(k,1)}(1,2),'Color',colSet_bluered(k,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',40);
%     h2 = errorbar(T(k,1),DLN_tot{k,omega_idx(1,k)}(1,1),DLN_tot{k,omega_idx(1,k)}(1,2),'Color',colSet_bluered(k,:),'Marker','d','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',20);
    axis square;
    x = xlabel ('T (K)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
end
hold off
box on
% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('N-T_free-omega')];
% print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2,frameName),'compact');
% close(gcf)





%%
%based on analysis on what energy ranges are swept out (modifications of
%the code of the second 2 sections above), we hand-make arrays to track the
%deep carrier energies and densities

%bounds of deep temp analysis
deep_temp1 = 240;
deep_temp2 = 290;

dt_idx1 = find(T(:,1) == deep_temp1);
dt_idx2 = find(T(:,1) == deep_temp2);

deep_trap = zeros(5,1+(dt_idx2-dt_idx1));

for k = 1:length(deep_trap)
    deep_trap(1,k) = T(k+9,1);
end

%%

for k = 1:length(deep_trap)
    deep_trap(3,k) = E_DLCP{1,k+9}(1,deep_trap(2,k));
    deep_trap(4,k) = DLN_tot{k+9,deep_trap(2,k)}(1,1);
    deep_trap(5,k) = DLN_tot{k+9,deep_trap(2,k)}(1,2);
end

figure()
hold on
for k = dt_idx1:dt_idx2
%     h1 = errorbar(T(k,1),E_DLCP{1,k}(1,free_idx(k,1)),E_DLCP_err{1,k}(1,free_idx(k,1)),'Color',colSet_bluered(k,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',40);
%     h2 = errorbar(T(k,1),E_DLCP{1,k}(1,omega_idx(1,k)),E_DLCP_err{1,k}(1,omega_idx(1,k)),'Color',colSet_bluered(k,:),'Marker','d','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',20);
    h1 = plot(deep_trap(1,k-9),deep_trap(3,k-9),'Color',colSet_bluered(k,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',40);
    h2 = plot(T(k,1),E_DLCP{1,k}(1,free_idx(k,1)),'Color',colSet_bluered(k,:),'Marker','d','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',20);
    axis square;
    x = xlabel ('T (K)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    y = ylabel ('E_{\omega} (eV)');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    xlim([230 300]) 
    set(gca,'xtick',[230 265 300]);
%     set(gca,'yscale','log');
end
hold off
box on
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('Ef_Edt_energies-temp')];
print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');
close(gcf)


figure();
hold on;
for k = dt_idx1:dt_idx2
    h1 = errorbar(deep_trap(1,k-9),deep_trap(4,k-9),deep_trap(5,k-9),'Color',colSet_bluered(k,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',40);
    h2 = errorbar(T(k,1),DLN_tot{k,free_idx(k,1)}(1,1),DLN_tot{k,free_idx(k,1)}(1,2),'Color',colSet_bluered(k,:),'Marker','d','LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:), 'MarkerSize',20);
    axis square;
    x = xlabel ('T (K)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    xlim([230 300]) 
    set(gca,'xtick',[230 265 300]);
    set(gca,'yscale','log');
end
hold off
box on
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('NEf_Ntd-temp')];
print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');
close(gcf)

%%


figure()
for k = dt_idx1:dt_idx2
    for ii = 1:v_max2 %index over voltages
        hold on
%         h2 = plot(pos_DLN{k,free_idx(k,1)}(:,2),pos_DLN{k,free_idx(k,1)}(:,3),'Color',colSet_bluered(k,:),'Marker', '.', 'LineWidth',3); %plot the filtered DLN vs. profile distance
        h2 = plot(pos_DLN{k,free_idx(k,1)}(:,2),pos_DLN{k,free_idx(k,1)}(:,3),'Color',colSet_bluered(k,:),'Marker', 'd', 'LineWidth', 3, 'MarkerFaceColor',colSet_bluered(k,:));    
    set(h2,{'markers'},{20})
    end
    axis square;
%     t = title(strcat('$N_{DL}$ vs. Distance,', num2str(T(i,1)), 'K'), 'FontSize',40); %Title the graph
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.09 .16]);
%     ylim([1e16 1e19]); 
    box on
end
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('N_free-T_spatial')];
print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');
close(gcf)


%%

figure()
for k = dt_idx1:dt_idx2
    for ii = 1:v_max2 %index over voltages
        hold on
        h2 = plot(pos_DLN{k,deep_trap(2,k-9)}(:,2),pos_DLN{k,deep_trap(2,k-9)}(:,3),'Color',colSet_bluered(k,:),'Marker', '.', 'LineWidth',3); %plot the filtered DLN vs. profile distance
%         h2 = plot(pos_DLN{k,deep_trap(2,k-9)}(:,2),pos_DLN{k,omega_idx(1,k)}(:,3),'Color',colSet_bluered(k,:),'Marker', 'd', 'LineWidth',3,'MarkerFaceColor',colSet_bluered(k,:));    
    set(h2,{'markers'},{50})
    end
    axis square;
%     t = title(strcat('$N_{DL}$ vs. Distance,', num2str(T(i,1)), 'K'), 'FontSize',40); %Title the graph
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.09 .17]);
%     ylim([1e16 1e19]); 
    box on
end
set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('N_deep_trap-T_spatial')];
print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');
close(gcf)

