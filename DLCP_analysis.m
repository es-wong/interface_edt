%This script analyzes the DLCP data and outputs the carrier number. The
%depletion width can then be determined from the carrier number and is
%output.

dc1 = input('What is the magnitude of the first DC bias in mV? ');
dc2 = input('What is the magnitude of the second DC bias in mV? ');
dc3 = input('What is the magnitude of the third DC bias in mV? ');
dc4 = input('What is the magnitude of the fourth DC bias in mV? ');
dc5 = input('What is the magnitude of the fifth DC bias in mV? ');
dc6 = input('What is the magnitude of the sixth DC bias in mV? ');

dc = [dc1 dc2 dc3 dc4 dc5 dc6];

v_max2 = size(dc,2);

field1 = 'index_1'; %create a structure to map the maximum voltage variables. Here, B is the strcucture containing the values for voltage.
value1 = dc(1,1);
field2 = 'index_2';
value2 = dc(1,2);
field3 = 'index_3';
value3 = dc(1,3);
field4 = 'index_4';
value4 = dc(1,4);
field5 = 'index_5';
value5 = dc(1,5);
field6 = 'index_6';
value6 = -300;
C = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6); %create the structure
CNames = fieldnames(C); %create a cell array with the field names so that we can call on them inside the plot for-loop
%%
fq1 = input('What is the first, lowest applied frequency? ');
fq2 = input('What is the second applied frequency? ');
fq3 = input('What is the third applied frequency? ');
fq4 = input('What is the fourth applied frequency? ');
% fq5 = input('What is the fifth applied frequency? ');
% fq6 = input('What is the sixth applied frequency? ');


freq = [fq1 fq2 fq3 fq4]'; %create an array for the frequencies used in the measurement.

f_max = size(freq,1);

field1 = 'index_1'; %create a structure to map the maximum voltage variables. Here, B is the strcucture containing the values for voltage.
value1 = freq(1,1);
field2 = 'index_2';
value2 = freq(2,1);
field3 = 'index_3';
value3 = freq(3,1);
field4 = 'index_4';
value4 = freq(4,1);
% field5 = 'index_5';
% value5 = freq(5,1);
% field6 = 'index_6';
% value6 = freq(6,1);
D = struct(field1,value1,field2,value2,field3,value3,field4,value4); %create the structure
DNames = fieldnames(D); %create a cell array with the field names so that we can call on them inside the plot for-loop

%%

%This section of the code renames the files into a manageable string so
%that we can import the data properly.

for j = t_start:t_step:t_stop %master index over temperature
    for m = 1:v_max2 %index over maximum bias v_max2
        for h = 1:f_max %index over frequencies
            oldname1 = strcat('dev',num2str(dev_number),'_T' , num2str(j),'K_F', num2str(D.(DNames{h,1})),'HZ_Vmax-', num2str(C.(CNames{m,1})), 'mV_DLCP.txt');    % Variable 'oldname' is a string that identifies the file
            newname1 = strcat('T_',num2str(j),'K_', num2str(h), '_', num2str(m),'.txt');   % Variable 'newname' generates a shorter name. 
            movefile(oldname1,newname1);  %Rewrite the longer filename with the newer, shorter name defined above.
        end
    end     
end

%do the same for the voltage check data

for j = t_start:t_step:t_stop %master index over temperature
    for m = dc_start:dc_step:dc_stop %index over frequencies
        for h = 1:f_max %index over maximum bias V_max
            oldname3 = strcat('dev',num2str(dev_number),'_T' , num2str(j),'K_F', num2str(D.(DNames{h,1})),'HZ_Vmax-', num2str(m), 'mV_DLCP_VCheck.txt');    % Variable 'oldname' is a string that identifies the file
            newname3 = strcat('T_',num2str(j),'K_', num2str(h), '_', num2str(m/100),'check.txt');   % Variable 'newname' generates a shorter name. 
            movefile(oldname3,newname3);  %Rewrite the longer filename with the newer, shorter name defined above.
        end
    end
end

%%

%This section imports the data into the workspace. First we build the outer
%layer of the cell with vertical axis temperature and horizontal axis
%frequeny. We then nest the maximum voltage data within than cell
%structure, so that we end up with a 2-tier cell holding all the raw data. 

DLCP = cell(t_max,v_max2);
check = cell(t_max,v_max2);

for j = 1:t_max 
    for m = 1:v_max2
        for h = 1:f_max 
            
            filename1 = strcat('T_',num2str(t_start+t_step*(j-1)),'K_', num2str(h), '_', num2str(m),'.txt'); %create for each 0 bias temperature point
            delimiter1 = ',';
            G = tdfread(filename1,delimiter1); %Import dataset into the active workspace; assign columns in each array 
            DLCP{j,m}{1,h}(:,1) = G.Amplitude;
            DLCP{j,m}{1,h}(:,2) = G.Capacitance;
            DLCP{j,m}{1,h}(:,3) = G.Dissipation;
      
            
        end
    end
end
%%
colSet_greenSub = makeColorMap([150 218 160]/255, [1 66 2]/254, f_max); %make a green colormap spanning the frquency space

T = zeros(t_max,1);
for k = 1:t_max
    T(k,1) = t_start+(k-1)*t_step;
end
%%
%Define constants that will be used later in our looped calculation. Note
%that the semiconductor radius is in m, as are the device area (m^2) and
%volume(m^3). Energies are in eV.
 
epsilon = 25.1;
A = 4e-6;  
eps_0 = 8.8541878e-12;
q = 1.60217657e-19;
kB = 8.6173324*(10^-5); 
A_cm = 4e-2; %area of the device in square cm
eps_0_cm =  8.8541878e-14; %vacuum permittivity in cm 


figuresdir2 = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\09_2017_m1\analysis\figures\DLCP\new';
% figuresdir2A = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\TBAI\06_2016_m1_d0\analysis\figures\DLCP\carrier_analysis';

%%

%Plot the data to check it before analysis 

for k = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max
            hold on
            scatter(DLCP{k,ii}{1,iii}(:,1),DLCP{k,ii}{1,iii}(:,2), 200, colSet_greenSub(iii,:), 'filled'); %for each subplot, plot the test frequencies
            t = title(strcat('V_{max} =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'FontName','Calibri')
            hold off
        end
        axis square;
        x = xlabel ('V_{RMS} (Volts)');
        set(x,'FontName', 'Calibri');
        set(gca,'FontSize',16);
        y = ylabel ('C (Farads)');
        set(y,'FontName','Calibri');
        set(gca,'FontSize',16);
        hold off
    end
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Cap_Amp_',num2str(T(k,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end

%%

%This section of the code fits the data with a decaying quadratic and
%outputs the fit parameters into a cell for the fit parameters C_fit. Each
%of the parameters is stored in 3 x 1 array in the cell. 

C_fits = cell(t_max,v_max2); 
C_line = cell(t_max,v_max2);

for k = 1:t_max
    for ii = 1:v_max2
        for iii = 1:f_max
            C_fits{k,ii}{1,iii}(:,:) = polyfit(DLCP{k,ii}{1,iii}(:,1),DLCP{k,ii}{1,iii}(:,2),2); %create the fit parameters and store them in a 1 x 3 array
            C_line{k,ii}{1,iii}(:,1) = (.01:.01:.3)'; %create the voltage domain for the fit line
            C_line{k,ii}{1,iii}(:,2) = C_fits{k,ii}{1,iii}(1,3) + (C_fits{k,ii}{1,iii}(1,2).* C_line{k,ii}{1,iii}(:,1)) + (C_fits{k,ii}{1,iii}(1,1).* C_line{k,ii}{1,iii}(:,1).^2); %calculate the best-ft function based on the calculated parameters
        end
    end
end

%replot the data with the best fit function

for k = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max
            hold on
            scatter(DLCP{k,ii}{1,iii}(:,1),DLCP{k,ii}{1,iii}(:,2), 200, colSet_greenSub(iii,:), 'filled'); %for each subplot, plot the test frequencies
            plot(C_line{k,ii}{1,iii}(:,1),C_line{k,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            t = title(strcat('V_{max} =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'FontName','Calibri')
            hold off
        end
        axis square;
        x = xlabel ('V_{RMS} (Volts)');
        set(x,'FontName', 'Calibri');
        set(gca,'FontSize',16);
        y = ylabel ('C (Farads)');
        set(y,'FontName','Calibri');
        set(gca,'FontSize',16);
        hold off
    end
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Cap_Amp_fit_',num2str(T(k,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end

%%

%This section of the code calculates N_dl and stores it as a number in the
%cell. 

N_dl = cell(t_max,v_max); %holds the data for the free carrier density


for k = 1:t_max 
    for ii = 1:v_max2 
        for iii =  1:f_max   
            N_dl{k,ii}(1,iii) = (-C_fits{k,ii}{1,iii}(1,3)^3/(2*q*epsilon*eps_0*(A^2)*C_fits{k,ii}{1,iii}(1,2))*1e-6); %calculate ionized gap density in cm^{-3}
        end
    end
end

%This section of the code plots the drive level density N_dl as a function
%of profile distance <x> calculated from the fit coefficients obtained
%above.

DLN = cell(t_max,f_max); 
pos_DLN = cell(t_max,f_max); 
DLN_test = zeros(v_max2,3);

for k = 1:t_max 
    for ii = 1:v_max2 
        for iii = 1:f_max 
            
            DLN{k,iii}(:,1) = -dc'; %write the V_{max} to the array for reference
            DLN{k,iii}(ii,2) = ((epsilon*eps_0*A)/C_fits{k,ii}{1,iii}(1,3))*1e6; %calculate the drive level <x> in units of microns 
            DLN{k,iii}(ii,3) = N_dl{k,ii}(1,iii); %import the N_dl data into the array
            
            DLN_test(:,1:3) = DLN{k,iii}(:,1:3); %import the unfiltered DLN sub arrays into a test array 
            DLN_pos = DLN_test(DLN_test(:,3) >= 0,:); %import only the data that contains positive carrier values in column 3 into a new test array DLN pos
            pos_DLN{k,iii} = DLN_pos; %re-import the test array into the new cell, filtering the data
            
        end
    end
end

%This section of the code averages over all the bias points to calculate the total carrier density 

DLN_tot = cell(t_max,f_max); %create a cell to store the drive-level density averaged over DC-bias

for k = 1:t_max
    for iii = 1:f_max
        
        DLN_tot{k,iii}(1,1) = mean(pos_DLN{k,iii}(:,3)); 
        DLN_tot{k,iii}(1,2) = std(pos_DLN{k,iii}(:,3))/numel(pos_DLN{k,iii}(:,3)); 
         
    end    
end

%plot Ndl vs <x>

for k = 1:t_max 
    figure()
    for ii = 1:v_max2 %index over voltages
        for iii = 1:3 %index over frequencies
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
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL_x_log_lowF_',num2str(T(k,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end

%%
%This section of the code determines the energy of the applied AC frequency
%used in the DLCP measurement

E_DLCP = cell(1,t_max); 
E_DLCP_err = cell(1,t_max);

%calculate the energy and coresponding error
for k = 1:t_max
    for iii = 1:f_max
        E_DLCP{1,k}(1,iii) = kB*T(k,1)*log(2*nu_0{1,1}(k,1)/(2*pi*freq(iii,1))); 
        E_DLCP_err{1,k}(1,iii) = E_DLCP{1,k}(1,iii)*(nu_0{1,1}(k,2)/nu_0{1,1}(k,1)); 
    end
end

%re-scale the N_{dl} vs. <x> plot to determine V_{dl} vs E_{\omega}
for k = 1+offset_ref(1,1):offset_ref(1,2)
    figure
    hold on
    for iii = 1:3
        errorbar(E_DLCP{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err{1,k}(1,iii),E_DLCP_err{1,k}(1,iii),'Color',colSet_greenSub(iii,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_blue(t_max,:), 'MarkerSize',40);  
    end
    hold off
    axis square;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.1 .35]);
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    hold off
     
    set(gca,'linewidth',1.5);
    box on
    
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL_E_log_both_',num2str(T(k,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     close(gcf)
end
%%
%make a total Ndl-E plot
markz = {'o','d','s', 'p', '>'};

figure
hold on
for k = 1+offset_ref(1,1):offset_ref(1,2)
%     figure
    hold on
    for iii = 1:3
        errorbar(E_DLCP{1,k}(1,iii),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),DLN_tot{k,iii}(1,2),E_DLCP_err{1,k}(1,iii),E_DLCP_err{1,k}(1,iii),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',3,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgecolor','k');
    end
    hold off
    axis square;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.1 .4]);
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
frameName = [strcat('NDL_E_log_total')];
print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');
close(gcf)
hold off


%%

%This section of the code selects the proper temperature and frequency
%ranges to continue the analysis

omega_dlcp = freq*2*pi;
freq_comp = zeros(t_max,length(omega_dlcp));
Ttmp = zeros(t_max,length(omega_dlcp));
t_cut = zeros(2,length(omega_dlcp));

for k = 1:t_max
    for j = 1:length(omega_dlcp)
        freq_comp(k,j) = ASpeaks(k,1) - omega_dlcp(j,1);
    end
end

%since we want omega_dlcp < omega_0, we need to select positive values
omega_mask = freq_comp>0;

%determine the temperature range over which DLCP data is valid
for j = 1:length(omega_dlcp)-1
    Ttmp(:,j) = omega_mask(:,j).*T(:,1);
    t_cut(1,j) = find(Ttmp(:,j)>0,1,'first');
    t_cut(2,j) = t_max;
end

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

T_DLCP = find(T(:,1) == input('What temperature should be used as a sample? '));

dc_DLCP = find(dc(1,:) == input('What is the magnitude of the maximum voltage that should be used as a sample? '));

%%

figure()
for k = T_DLCP
    for ii = dc_DLCP
        hold on
        for iii = 1:3
            hold on
            
            if iii == 1
                h1 = plot(DLCP_norm{k,ii}{1,iii}(:,1),DLCP_norm{k,ii}{1,iii}(:,2),'Marker','o','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k');
                set(h1,{'MarkerSize'},{20});
                set(h1,{'LineWidth'},{1});
            
                plot(C_line_norm{k,ii}{1,iii}(:,1),C_line_norm{k,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
                
            elseif iii == 2
                h1 = plot(DLCP_norm{k,ii}{1,iii}(:,1),DLCP_norm{k,ii}{1,iii}(:,2),'Marker','d','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k');
                set(h1,{'MarkerSize'},{20});
                set(h1,{'LineWidth'},{2});
            
                plot(C_line_norm{k,ii}{1,iii}(:,1),C_line_norm{k,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            
            
            elseif iii == 3
                h1 = plot(DLCP_norm{k,ii}{1,iii}(:,1),DLCP_norm{k,ii}{1,iii}(:,2),'Marker','s','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k');
                set(h1,{'MarkerSize'},{20});
                set(h1,{'LineWidth'},{2});
            
                plot(C_line_norm{k,ii}{1,iii}(:,1),C_line_norm{k,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            
            end
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
%         ylim([0 220]);
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5); 
    hold off;
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('cap_amp_fit_samp_',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2A,frameName),'compact');
    
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    
    close(gcf)
    end
end
%%
%Make the sample N_{DL} vs. <x> plots

for k = T_DLCP 
    figure()
    for ii = 1:v_max2 
        for iii = 1:3
            hold on
            
            if iii == 1
            h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', 'o','MarkerEdgeColor','k','MarkerFacecolor',colSet_green(k,:));
            set(h1,{'MarkerSize'},{20});
            set(h1,{'LineWidth'},{1});
            
            elseif iii == 2
            h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', 'd','MarkerEdgeColor','k','MarkerFacecolor',colSet_green(k,:));
            set(h1,{'MarkerSize'},{20});
            set(h1,{'LineWidth'},{3});
            
            elseif iii == 3
            h1 = plot(1e3*pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', 's','MarkerEdgeColor','k','MarkerFacecolor',colSet_green(k,:));
            set(h1,{'MarkerSize'},{20});
            set(h1,{'LineWidth'},{3});
            
            end       
        end
    end
    axis square;
     box on

    x = xlabel ('<x> (nm)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([32 36]);
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
%     ylim([1e16 1e18]);
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5); 
    hold off;
 

    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL_x_samp_',num2str(T(k,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2A,frameName),'compact');
    
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    
    close(gcf)
end

%%

%plot a distance profile of the carriers at different temperatures to
%determine where the response is localized within the film

for iii = 1
figure()
% for k = t_cut(1,1):t_max
for k = 1+offset_ref(1,1):offset_ref(1,2)
% for k = 1:4:5
    for ii = 1:v_max2 %index over voltages
        hold on
        

        h2 = plot(pos_DLN{k,iii}(:,2),pos_DLN{k,iii}(:,3),'Color',colSet_green(k,:),'Marker', markz{1,iii}, 'LineWidth',1,'MarkerFaceColor',colSet_green(k,:),'MarkerEdgeColor','k');    
            set(h2,{'markers'},{20})

    end
    axis square
    box on

    x = xlabel ('<x> (\mum)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.03 .09]);
    set(gca,'xtick',[.03 .05 .07 .09])
        
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    ylim([1e17 1e19]); 
%     set(gca,'xtick',[1e16 1e17 3e17])

    pbaspect([1 .6 1]);
    set(gca,'linewidth',1.5);
    hold off
% 
end

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('N_DL-x_offset_',num2str(freq(iii,1)/1e3),'_kHz')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2A,frameName),'compact');

print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');

close(gcf)
end
%%
%Make carrier density vs. temperature plots

figure();
hold on;
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for iii = 4:5
    
    h1 = errorbar(T(k,1),DLN_tot{k,iii}(1,1),DLN_tot{k,iii}(1,2),'Color',colSet_green(k,:),'Marker',markz{1,iii},'LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
%     h2 = errorbar(T(k,1),DLN_tot{k,2}(1,1),DLN_tot{k,2}(1,2),'Color',colSet_green(k,:),'Marker','d','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgeColor','k');
%     h3 = errorbar(T(k,1),DLN_tot{k,3}(1,1),DLN_tot{k,3}(1,2),'Color',colSet_green(k,:),'Marker','s','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20,'MarkerEdgeColor','k');
    
    end
    axis square
    box on
    
    x = xlabel ('T (K)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([180 300])
    set(gca,'xtick',[180 220 260 300])
    
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    ylim([1e17 2e19])
    set(gca,'ytick',[1e17 1e18 1e19]);
    set(gca,'yscale','log');
    
end


pbaspect([1.6 1 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('N_DL-T_off1_highF')];
% print(gcf, '-dpng', strcat(figuresdir2A,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir2A,frameName),'compact'); 

print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir2,frameName),'compact');

close(gcf)


%%
