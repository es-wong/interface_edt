%This script analyzes the DLCP data and outputs the carrier number. The
%depletion width can then be determined from the carrier number and is
%output.

dc_start = input('What is the magnitude of the start DC bias in mV? ');
dc_step = input('What is the magnitude of the step DC bias in mV? ');
dc_stop = input('What is the magnitude of the final DC bias in mV? ');

v_max2 = ((dc_stop-dc_start) + dc_step)/dc_step;

fq1 = input('What is the first, lowest applied frequency? ');
fq2 = input('What is the second applied frequency? ');
fq3 = input('What is the third applied frequency? ');
fq4 = input('What is the fourth applied frequency? ');
fq5 = input('What is the fifth applied frequency? ');
fq6 = input('What is the sixth applied frequency? ');


freq = [fq1 fq2 fq3 fq4 fq5 fq6]'; %create an array for the frequencies used in the measurement.

f_max = size(freq,1);

field1 = 'index_1'; %create a structure to map the maximum voltage variables. Here, B is the strcucture containing the values for voltage.
value1 = freq(1,1);
field2 = 'index_2';
value2 = freq(2,1);
field3 = 'index_3';
value3 = freq(3,1);
field4 = 'index_4';
value4 = freq(4,1);
field5 = 'index_5';
value5 = freq(5,1);
field6 = 'index_6';
value6 = freq(6,1);
D = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6); %create the structure
DNames = fieldnames(D); %create a cell array with the field names so that we can call on them inside the plot for-loop



%%

%This section of the code renames the files into a manageable string so
%that we can import the data properly.


for j = t_start:t_step:t_stop %master index over temperature
    for m = dc_start:dc_step:dc_stop %index over maximum bias v_max2
        for k = 1:f_max; %index over frequencies
            oldname1 = strcat('dev',num2str(dev_number),'_T' , num2str(j),'K_F', num2str(D.(DNames{k,1})),'HZ_Vmax-', num2str(m), 'mV_DLCP.txt');    % Variable 'oldname' is a string that identifies the file
            newname1 = strcat('T_',num2str(j),'K_', num2str(k), '_', num2str(m/dc_step),'.txt');   % Variable 'newname' generates a shorter name. 
            movefile(oldname1,newname1);  %Rewrite the longer filename with the newer, shorter name defined above.
        end;
    end;      
end;
%%
for j = t_start:t_step:t_stop %master index over temperature
    for m = dc_start:dc_step:dc_stop %index over maximum bias v_max2
        for k = 1:f_max; %index over frequencies
            oldname2 = strcat('dev',num2str(dev_number),'_IR__T' , num2str(j),'K_F', num2str(D.(DNames{k,1})),'HZ_Vmax-', num2str(m), 'mV_DLCP.txt');    % Variable 'oldname' is a string that identifies the file
            newname2 = strcat('IR_T_',num2str(j),'K_', num2str(k), '_', num2str(m/dc_step),'.txt');   % Variable 'newname' generates a shorter name. 
            movefile(oldname2,newname2);  %Rewrite the longer filename with the newer, shorter name defined above.
        end
    end;
end;

%%
%do the same for the voltage check data

for j = t_start:t_step:t_stop %master index over temperature
    for m = dc_start:dc_step:dc_stop %index over frequencies
        for k = 1:f_max; %index over maximum bias V_max
            oldname3 = strcat('dev',num2str(dev_number),'_T' , num2str(j),'K_F', num2str(D.(DNames{k,1})),'HZ_Vmax-', num2str(m), 'mV_DLCP_VCheck.txt');    % Variable 'oldname' is a string that identifies the file
            newname3 = strcat('T_',num2str(j),'K_', num2str(k), '_', num2str(m/100),'check.txt');   % Variable 'newname' generates a shorter name. 
            movefile(oldname3,newname3);  %Rewrite the longer filename with the newer, shorter name defined above.
        end
    end;
end;
%%
for j = t_start:t_step:t_stop %master index over temperature
    for m = dc_start:dc_step:dc_stop %index over frequencies
        for k = 1:f_max; %index over maximum bias V_max
            oldname4 = strcat('dev',num2str(dev_number),'_IR__T' , num2str(j),'K_F', num2str(D.(DNames{k,1})),'HZ_Vmax-', num2str(m), 'mV_DLCP_VCheck.txt');    % Variable 'oldname' is a string that identifies the file
            newname4 = strcat('IR_T_',num2str(j),'K_', num2str(k), '_', num2str(m/100),'check.txt');   % Variable 'newname' generates a shorter name. 
            movefile(oldname4,newname4);  %Rewrite the longer filename with the newer, shorter name defined above.
        end
    end;
end;

%%

%This section imports the data into the workspace. First we build the outer
%layer of the cell with vertical axis temperature and horizontal axis
%frequeny. We then nest the maximum voltage data within than cell
%structure, so that we end up with a 2-tier cell holding all the raw data. 

DLCP = cell(t_max,v_max2);
% IRDLCP = cell(t_max,v_max2);
check = cell(t_max,v_max2);

for j = 1:t_max; %index across temperature
    for m = 1:v_max2; %index across maximum voltage
        for k = 1:f_max %index across frequency
            filename1 = strcat('T_',num2str(t_start+t_step*(j-1)),'K_', num2str(k), '_', num2str(m),'.txt'); %create for each 0 bias temperature point
            delimiter1 = ',';
            C = tdfread(filename1,delimiter1); %Import dataset into the active workspace; assign columns in each array 
            DLCP{j,m}{1,k}(:,1) = C.Amplitude;
            DLCP{j,m}{1,k}(:,2) = C.Capacitance;
            DLCP{j,m}{1,k}(:,3) = C.Dissipation;
            
%             
%             filename2 = strcat('IR_T_',num2str(t_start+t_step*(j-1)),'K_', num2str(k), '_', num2str(m),'.txt'); %create for each 0 bias temperature point
%             delimiter2 = ',';
%             E = tdfread(filename2,delimiter2); %Import dataset into the active workspace; assign columns in each array 
%             IRDLCP{j,m}{1,k}(:,1) = E.Amplitude;
%             IRDLCP{j,m}{1,k}(:,2) = E.Capacitance;
%             IRDLCP{j,m}{1,k}(:,3) = E.Dissipation;
            
        end;
    end;
end;

% for j = 1:t_max; %index across voltage
%     for m = 1:v_max2; %index across maximum voltage
%         for k = 1:f_max; %index over frequency
%             filename = strcat('T_',num2str(180+10*(j-1)),'K_', num2str(k), '_', num2str(m),'check.txt'); %create for each 0 bias temperature point
%             delimiter = ',';
%             D = tdfread(filename,delimiter); %Import dataset into the active workspace; assign columns in each array 
%             check{j,m}{1,k}(:,1) = D.Amplitude;
%             check{j,m}{1,k}(:,2) = D.DC_Bias;
%         end;
%     end;
% end;

colSet_green = makeColorMap([150 218 160]/255, [1 66 2]/254, f_max); %make a green colormap spanning the frquency space


field1 = 'index_1'; %create a structure to map the maximum voltage variables. Here, B is the strcucture containing the values for voltage.
value1 = -50;
field2 = 'index_2';
value2 = -100;
field3 = 'index_3';
value3 = -150;
field4 = 'index_4';
value4 = -200;
field5 = 'index_5';
value5 = -250;
field6 = 'index_6';
value6 = -300;
C = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6); %create the structure
CNames = fieldnames(C); %create a cell array with the field names so that we can call on them inside the plot for-loop


T = zeros(t_max,1);
for i = 1:t_max;
    T(i,1) = t_start+(i-1)*t_step;
end;

A = 4e-6;  %Define constants that will be used later in our looped calculation. Note that the semicondcutor radius is in m, as are the device area (m^2) and volume (m^3). Energy parameters are in eV.  
epsilon = 11.5;
eps_0 = 8.8541878e-12;
q = 1.60217657e-19;
kB = 8.6173324*(10^-5); 
A_cm = 4e-2; %area of the device in square cm
eps_0_cm =  8.8541878e-14; %vacuum permittivity in cm 


figuresdir2 = 'C:\Users\Eric\Desktop\Admittance\ZnO-PbS\light_soak_study\IR_soak\sanity_checks\EDT_schottky_junction\try2\analysis\DLCP';


%%
%Plot the data to check it before analysis 

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max;
            hold on
            scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet_green(iii,:), 'filled'); %for each subplot, plot the test frequencies
            t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'Interpreter','Latex')
%             if ii == v_max2;
%                 subplot(round(v_max2/2), 2,[ii,ii+1]); %divide the figure into a grid 
%                 scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet1a(iii,:), 'filled'); %for each subplot, plot the test frequencies
%             end;
            hold off
        end;
        axis square;
        x = xlabel ('AC Amplitude, $V_{AC}$ (Volts)');
        set(x,'Interpreter', 'Latex');
        set(gca,'FontSize',16);
        y = ylabel ('Capacitance, C (Farads)');
        set(y,'Interpreter','Latex');
        set(gca,'FontSize',16);
%         suptitle(strcat('Capacitance vs AC Amplitude,_',num2str(T(i,1)),'K'));
        hold off
    end;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('Cap_Amp_',num2str(T(i,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)
end;

%%

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max;
            hold on
            scatter(IRDLCP{i,ii}{1,iii}(:,1),IRDLCP{i,ii}{1,iii}(:,2), 200, colSet_green(iii,:), 'd','filled'); %for each subplot, plot the test frequencies
            t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'Interpreter','Latex')
%             if ii == v_max2;
%                 subplot(round(v_max2/2), 2,[ii,ii+1]); %divide the figure into a grid 
%                 scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet1a(iii,:), 'filled'); %for each subplot, plot the test frequencies
%             end;
            hold off
        end;
        axis square;
        x = xlabel ('AC Amplitude, $V_{AC}$ (Volts)');
        set(x,'Interpreter', 'Latex');
        set(gca,'FontSize',16);
        y = ylabel ('Capacitance, C (Farads)');
        set(y,'Interpreter','Latex');
        set(gca,'FontSize',16);
%         suptitle(strcat('Capacitance vs AC Amplitude,_',num2str(T(i,1)),'K'));
        hold off
    end;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Cap_Amp_IR_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end

%%
%This section of the code fits the data with a decaying quadratic and
%outputs the fit parameters into a cell for the fit parameters C_fit. Each
%of the parameters is stored in 3 x 1 array in the cell. 

C_fits = cell(t_max,v_max2); %holds data for fit parameters
C_line = cell(t_max,v_max2); %holds data for the line calculated from the fits


for i = 1:t_max;
    for ii = 1:v_max2;
        for iii = 1:f_max;
            C_fits{i,ii}{1,iii}(:,:) = polyfit(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2),2); %create the fit parameters and store them in a 1 x 3 array
            C_line{i,ii}{1,iii}(:,1) = (.01:.01:.3)'; %create the voltage domain for the fit line
            C_line{i,ii}{1,iii}(:,2) = C_fits{i,ii}{1,iii}(1,3) + (C_fits{i,ii}{1,iii}(1,2).* C_line{i,ii}{1,iii}(:,1)) + (C_fits{i,ii}{1,iii}(1,1).* C_line{i,ii}{1,iii}(:,1).^2); %calculate the best-ft function based on the calculated parameters
        end;
    end;  
end;
%%
IRC_fits = cell(t_max,v_max2); %holds data for fit parameters
IRC_line = cell(t_max,v_max2); %holds data for the line calculated from the fits


l_start1 = input('At which AC amplitude should the fit to the illuminated capacitance curve start? ');

l_start = find(IRDLCP{1,1}{1,1}(:,1) == l_start1);


for i = 1:t_max;
    for ii = 1:v_max2;
        for iii = 1:f_max;
            IRC_fits{i,ii}{1,iii}(:,:) = polyfit(IRDLCP{i,ii}{1,iii}(l_start:end,1),IRDLCP{i,ii}{1,iii}(l_start:end,2),2); %create the fit parameters and store them in a 1 x 3 array
            IRC_line{i,ii}{1,iii}(:,1) = (l_start1:.01:.3)'; %create the voltage domain for the fit line
            IRC_line{i,ii}{1,iii}(:,2) = IRC_fits{i,ii}{1,iii}(1,3) + (IRC_fits{i,ii}{1,iii}(1,2).* IRC_line{i,ii}{1,iii}(:,1)) + (IRC_fits{i,ii}{1,iii}(1,1).* IRC_line{i,ii}{1,iii}(:,1).^2); %calculate the best-ft function based on the calculated parameters
        end;
    end;  
end;

%%
%replot the data with the best fit function

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max;
            hold on
            scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet_green(iii,:), 'filled'); %for each subplot, plot the test frequencies
            plot(C_line{i,ii}{1,iii}(:,1),C_line{i,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'Interpreter','Latex')
%             if ii == v_max2;
%                  subplot(round(v_max2/2),2,[ii,ii+1]);
%                  scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet1(iii,:), 'filled'); %for each subplot, plot the test frequencies
%                  plot(C_line{i,ii}{1,iii}(:,1),C_line{i,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
%                  t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
%                  set(t,'Interpreter','Latex');
%             end;
            hold off
        end;
        axis square;
        x = xlabel ('AC Amplitude, $V_{AC}$ (Volts)');
        set(x,'Interpreter', 'Latex');
        set(gca,'FontSize',16);
        y = ylabel ('Capacitance, C (Farads)');
        set(y,'Interpreter','Latex');
        set(gca,'FontSize',16);
        hold off
    end;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Cap_Amp_fit_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end;
%%

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max;
            hold on
            scatter(IRDLCP{i,ii}{1,iii}(:,1),IRDLCP{i,ii}{1,iii}(:,2), 200, colSet_green(iii,:),'d', 'filled'); %for each subplot, plot the test frequencies
            plot(IRC_line{i,ii}{1,iii}(:,1),IRC_line{i,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
            t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'Interpreter','Latex')
%             if ii == v_max2;
%                  subplot(round(v_max2/2),2,[ii,ii+1]);
%                  scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet1(iii,:), 'filled'); %for each subplot, plot the test frequencies
%                  plot(C_line{i,ii}{1,iii}(:,1),C_line{i,ii}{1,iii}(:,2),'r', 'LineWidth', 3);
%                  t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
%                  set(t,'Interpreter','Latex');
%             end;
            hold off
        end;
        axis square;
        x = xlabel ('AC Amplitude, $V_{AC}$ (Volts)');
        set(x,'Interpreter', 'Latex');
        set(gca,'FontSize',16);
        y = ylabel ('Capacitance, C (Farads)');
        set(y,'Interpreter','Latex');
        set(gca,'FontSize',16);
        hold off
    end;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Cap_Amp_fit_IR_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end;

%%

%This section of the code calculates N_dl and stores it as a number in the
%cell. 

N_dl = cell(t_max,v_max); %holds the data for the free carrier density
% IRN_dl = cell(t_max,v_max);

for i = 1:t_max; %index over temperatures
    for ii = 1:v_max2; %index over max-voltages
        for iii =  1:f_max; %index over frequencies
            N_dl{i,ii}(1,iii) = (-C_fits{i,ii}{1,iii}(1,3)^3/(2*q*epsilon*eps_0*(A^2)*C_fits{i,ii}{1,iii}(1,2))*1e-6); %calculate ionized gap density in cm^{-3}
            
%             IRN_dl{i,ii}(1,iii) = (-IRC_fits{i,ii}{1,iii}(1,3)^3/(2*q*epsilon*eps_0*(A^2)*IRC_fits{i,ii}{1,iii}(1,2))*1e-6); %calculate ionized gap density in cm^{-3}
        end;
    end;
end;

pos_Ndl = cell(t_max,v_max); %create a cell to hold the positive values of N_dl
IRpos_Ndl = cell(t_max,v_max); %create a cell to hold the positive values of N_dl

for i = 1:t_max %index over temperatures
    for ii = 1:v_max2; %index over maximum voltages
        testN = N_dl{i,ii}; %extract the N_dl values into a test array
        testN(testN <0) = 0; %assign zeros to the negative values in the test array
        pos_Ndl{i,ii} = testN; %re-import the test array into the positive cell
%         
%         IRtestN = IRN_dl{i,ii}; %extract the N_dl values into a test array
%         IRtestN(testN <0) = 0; %assign zeros to the negative values in the test array
%         IRpos_Ndl{i,ii} = IRtestN; %re-import the test array into the positive cell
    end;
end;

%%
%This section of the code plots the drive level density N_dl as a function
%of profile distance <x> calculated from the fit coefficients obtained
%above.

DLN = cell(t_max,f_max); %initialize a cell that will hold the data to plot drive level density vs. profile distance
pos_DLN = cell(t_max,f_max); %initialize a cell that will hold the cleaned drive level density vs. profile distance 
DLN_test = zeros(v_max2,3); %initialize a test array on which  to run the cleaning functions below


% IRDLN = cell(t_max,f_max); 
% IRpos_DLN = cell(t_max,f_max); 
% IRDLN_test = zeros(v_max2,3);


for i = 1:t_max; %index over temperatures
    for ii = 1:v_max2 %index over V_{max}
        for iii = 1:f_max; %index over frequencies
            DLN{i,iii}(:,1) = (-dc_start:-dc_step:-dc_stop)'; %write the V_{max} to the array for reference
            DLN{i,iii}(ii,2) = ((epsilon*eps_0*A)/C_fits{i,ii}{1,iii}(1,3))*1e6; %calculate the drive level <x> in units of microns 
            DLN{i,iii}(ii,3) = N_dl{i,ii}(1,iii); %import the N_dl data into the array
            
%             IRDLN{i,iii}(:,1) = (-dc_start:-dc_step:-dc_step)'; %write the V_{max} to the array for reference
%             IRDLN{i,iii}(ii,2) = ((epsilon*eps_0*A)/IRC_fits{i,ii}{1,iii}(1,3))*1e6; %calculate the drive level <x> in units of microns 
%             IRDLN{i,iii}(ii,3) = IRN_dl{i,ii}(1,iii); %import the N_dl data into the array
        end;
    end;
end;


for i =1:t_max %index over temperature
    for ii = 1:v_max2 %index over voltage
        for iii = 1:f_max %index over frequency
            DLN_test(:,1:3) = DLN{i,iii}(:,1:3); %import the unfiltered DLN sub arrays into a test array 
            DLN_pos = DLN_test(DLN_test(:,3) >= 0,:); %import only the data that contains positive carrier values in column 3 into a new test array DLN pos
            pos_DLN{i,iii} = DLN_pos; %re-import the test array into the new cell, filtering the data
            
%             IRDLN_test(:,1:3) = IRDLN{i,iii}(:,1:3); %import the unfiltered DLN sub arrays into a test array 
%             IRDLN_pos = IRDLN_test(IRDLN_test(:,3) >= 0,:); %import only the data that contains positive carrier values in column 3 into a new test array DLN pos
%             IRpos_DLN{i,iii} = IRDLN_pos; %re-import the test array into the new cell, filtering the data
        end;
    end;
end;

%%

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2; %index over voltages
        for iii = 1:f_max-1 %index over frequencies
            hold on
            h1 = plot(pos_DLN{i,iii}(:,2),pos_DLN{i,iii}(:,3),'Color',colSet_green(iii,:),'Marker', '.', 'LineWidth',3);
%             h2 = plot(IRpos_DLN{i,iii}(:,2),IRpos_DLN{i,iii}(:,3),'Color',colSet_green(iii,:),'Marker', 'd', 'LineWidth',3); %plot the filtered DLN vs. profile distance
                set(h1,{'markers'},{50});
%                 set(h2,{'markers'},{10})
        end;
    end;

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
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL_x_log_both_',num2str(T(i,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     close(gcf)
end;

%%

%This section of the code calculates the measurment energy E_{omega}
%(E_DLCP in the code) of the measurement and plots the defect density N_dl 
%against E_{omega} so that the energy of the defect can be deduced.

DLN_tot = cell(t_max,f_max); %create a cell to store the drive-level density averaged over DC-bias
N_diff = zeros(t_max,2); %create a cell to store the difference in trapped and free carriers

IRDLN_tot = cell(t_max,f_max); %do the same for the IR data
IRN_diff = zeros(t_max,2);

free_index = input('Which frequency index is the free carrier index?: ');

for i = 1:t_max
    for iii = 1:f_max;
        
        DLN_tot{i,iii}(1,1) = mean(pos_DLN{i,iii}(:,3)); %average the drive level densities over DC-bias
        DLN_tot{i,iii}(1,2) = std(pos_DLN{i,iii}(:,3)); %calculate the 1-sigma variation in the drive level density
        
%         IRDLN_tot{i,iii}(1,1) = mean(IRpos_DLN{i,iii}(:,3)); %do the same for the IR data
%         IRDLN_tot{i,iii}(1,2) = std(IRpos_DLN{i,iii}(:,3));
        
    end;
    
    N_diff(i,1) = DLN_tot{i,1}(1,1) - DLN_tot{i,free_index}(1,1); %calculate the difference in the free and trapped carriers
    N_diff(i,2) = sqrt(DLN_tot{i,1}(1,2)^2 + DLN_tot{i,free_index}(1,2)^2); %propogate the 1-sgima error through for the difference
    
     
%     IRN_diff(i,1) = IRDLN_tot{i,1}(1,1) - IRDLN_tot{i,free_index}(1,1); %do the same for the IR data
%     IRN_diff(i,2) = sqrt(IRDLN_tot{i,1}(1,2)^2 + IRDLN_tot{i,free_index}(1,2)^2);
    
end;

for iii = 1:f_max;
figure()
hold on;
    for i = 1:t_max; 
        errorbar(T(i,1),DLN_tot{i,iii}(1,1),DLN_tot{i,iii}(1,2),'Color',colSet_green(iii,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_blue(t_max,:), 'MarkerSize',40);
        
%         errorbar(T(i,1),IRDLN_tot{i,iii}(1,1),IRDLN_tot{i,iii}(1,2),'Color',colSet_green(iii,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);
%         scatter(T(i,1),DLN_tot{i,iii}(1,1),250,colSet3(iii,:),'filled');
    end;
    axis square;
    t = title(strcat('$N_{DL}$ vs. Temperature'), 'FontSize',40); %Title the graph
    set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('$<T>$ ($K$)');
    set(x,'Interpreter', 'Latex');
    set(gca,'FontSize',40);
    y = ylabel ('$N_{DL}$, ($cm^{-3}$)');
    set(y,'Interpreter','Latex');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.1 .35]);
    xlim([195 315]); 
    box on
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL-T_log_both_',num2str(freq(iii,1)./1e3),'kHz')];
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     close(gcf)
end;
%%
E_DLCP = cell(1,t_max); %create an array to store the measurement energy for DLCP
E_DLCP_err = cell(1,t_max); %create an array to store the uncertaity in the measurement energy calculation

for i = 1:t_max;
    for iii = 1:f_max;
        E_DLCP{1,i}(1,iii) = kB*T(i,1)*log(nu_0(1,1)/(2*pi*freq(iii,1))); %calculate the measurement energy for DLCP
%         E_DLCP{1,i}(2,iii) = kB*T(i,1)*log(IR_nu_0(1,1)/(2*pi*freq(iii,1))); %do the same for the IR data
    end;
end;

for i = 1:t_max;
    for iii = 1:f_max;
        E_DLCP_err{1,i}(1,iii) = kB*T(i,1)*log(1/(2*pi*freq(iii,1)))*(nu_0(2,1)/nu_0(1,1)); %calculate the measurement energy for DLCP
%         E_DLCP_err{1,i}(2,iii) = kB*T(i,1)*log(1/(2*pi*freq(iii,1)))*(IR_nu_0(2,1)/IR_nu_0(1,1)); %do the same for the IR data
    end;
end;



%%
%plot the measurement energy vs. the total carrier number
for i = 1:t_max
    figure
    hold on
    for iii = 1:f_max;
        errorbar(E_DLCP{1,i}(1,iii),DLN_tot{i,iii}(1,1),DLN_tot{i,iii}(1,2),'Color',colSet_blue(t_max,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_blue(t_max,:), 'MarkerSize',40);
%         errorbar(E_DLCP{1,i}(2,iii),IRDLN_tot{i,iii}(1,1),IRDLN_tot{i,iii}(1,2),'Color',colSet_red(t_max,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);
    end;
    hold off
    axis square;
%     t = title(strcat('$N_{DL}$ vs. Energy'), 'FontSize',40); %Title the graph
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('Energy, E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('N_{DL} (cm^{-3})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.1 .35]);
%     xlim([255 335]); 
    box on
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL_E_log_both_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end;

%%
tau = zeros(1,f_max); %create an array to store the equivalent carrier lifetime for each of the applied frequencies 

for iii = 1:f_max;
    tau(1,iii) = 1/freq(iii,1); %calculate the equivalent carrier lifetime for each of the applied frequencies  
end;

%plot the number of carriers responding vs equivalent carrier lifetime for each of the applied frequencies

for i = 1:t_max
    figure
    hold on
    for iii = 1:f_max;
        errorbar(tau(1,iii),DLN_tot{i,iii}(1,1),DLN_tot{i,iii}(1,2),'Color',colSet_blue(t_max,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_blue(t_max,:), 'MarkerSize',40);
%         errorbar(tau(1,iii),IRDLN_tot{i,iii}(1,1),IRDLN_tot{i,iii}(1,2),'Color',colSet_red(t_max,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);
    end;
    hold off
    axis square;
    t = title(strcat('$N_{DL}$ vs. Lifetime'), 'FontSize',40); %Title the graph
    set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('Lifetime, $\tau$ ($s$)');
    set(x,'Interpreter', 'Latex');
    set(gca,'FontSize',40);
    y = ylabel ('$N_{DL}$ ($cm^{-3}$)');
    set(y,'Interpreter','Latex');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    hold off
%     xlim([.1 .35]);
%     xlim([255 335]); 
    box on
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('NDL_tau_log_both_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end;

%%
%This section of the code plots the dissipation of both the dark and
%IR-illuminated curves to ensure that it is below 5 and the capacitance
%data points are (roughly) valid

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max;
            hold on
            scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,3), 200, colSet_green(iii,:), 'filled'); %for each subplot, plot the test frequencies
            t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'Interpreter','Latex')
%             if ii == v_max2;
%                 subplot(round(v_max2/2), 2,[ii,ii+1]); %divide the figure into a grid 
%                 scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet1a(iii,:), 'filled'); %for each subplot, plot the test frequencies
%             end;
            hold off
        end;
        axis square;
        x = xlabel ('AC Amplitude, $V_{AC}$ (Volts)');
        set(x,'Interpreter', 'Latex');
        set(gca,'FontSize',16);
        y = ylabel ('D ($\frac{G}{\omega C}$)');
        set(y,'Interpreter','Latex');
        set(gca,'FontSize',16);
%         suptitle(strcat('Capacitance vs AC Amplitude,_',num2str(T(i,1)),'K'));
        hold off
    end;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('D_amp_',num2str(T(i,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir2,frameName),'compact');
%     close(gcf)
end;

%%

for i = 1:t_max %index over temperatures
    figure()
    for ii = 1:v_max2 %index over maximum voltages
        hold on
        subplot(round(v_max2/2),2,ii); %divide the figure into a grid 
        for iii = 1:f_max-1;
            hold on
            scatter(IRDLCP{i,ii}{1,iii}(:,1),IRDLCP{i,ii}{1,iii}(:,3), 200, colSet_green(iii,:), 'd','filled'); %for each subplot, plot the test frequencies
            t = title(strcat('$V_{max}$ =', num2str(C.(CNames{ii,1})), 'mV'), 'FontSize', 24);
            set(t,'Interpreter','Latex')
%             if ii == v_max2;
%                 subplot(round(v_max2/2), 2,[ii,ii+1]); %divide the figure into a grid 
%                 scatter(DLCP{i,ii}{1,iii}(:,1),DLCP{i,ii}{1,iii}(:,2), 200, colSet1a(iii,:), 'filled'); %for each subplot, plot the test frequencies
%             end;
            hold off
        end;
        axis square;
        x = xlabel ('AC Amplitude, $V_{AC}$ (Volts)');
        set(x,'Interpreter', 'Latex');
        set(gca,'FontSize',16);
        y = ylabel ('D ($\frac{G}{\omega C}$)');
        set(y,'Interpreter','Latex');
        set(gca,'FontSize',16);
%         suptitle(strcat('Capacitance vs AC Amplitude,_',num2str(T(i,1)),'K'));
        hold off
    end;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('D_amp_IR_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end
%%
diss_amp = cell(t_max,f_max); %create a cell to store the disspation values 
% IR_diss_amp = cell(t_max,f_max);%create a cell to store the illuminated dissipation values

diss_avg = cell(t_max,f_max); %create a cell that averages the dissipation over voltages
% IR_diss_avg = cell(t_max,f_max); %do the same for the illuminated values

%%

for i = 1:t_max
    for ii = 1:v_max2;
        for iii = 1:f_max;
            
           diss_amp{i,iii}(ii,1) = mean(DLCP{i,ii}{1,iii}(:,3)); %average the dissipation values over amplitude
           diss_amp{i,iii}(ii,2) = std(DLCP{i,ii}{1,iii}(:,3)); %calculate the 1-sigma error in dissipation over amplitude
%            IR_diss_amp{i,iii}(ii,1) =  mean(IRDLCP{i,ii}{1,iii}(:,3)); %do the same for the illuminated data
%            IR_diss_amp{i,iii}(ii,2) =  std(IRDLCP{i,ii}{1,iii}(:,3));
           
           diss_avg{i,iii}(1,1) = mean(diss_amp{i,iii}(:,1)); %average the disipation values for each frequency
           diss_avg{i,iii}(1,2) = std(diss_amp{i,iii}(:,1)); %calculate the 1-sigma error in disipation for each frequency
%            IR_diss_avg{i,iii}(1,1) = mean(IR_diss_amp{i,iii}(:,1)); %do the same for the illuminated data
%            IR_diss_avg{i,iii}(1,2) = std(IR_diss_amp{i,iii}(:,1));
           
        end;
    end;
end;


%%
%plot the averaged dissipation value at each energy
for i = 1:t_max
    figure
    hold on
    for iii = 1:f_max;
        errorbar(E_DLCP{1,i}(1,iii),diss_avg{i,iii}(1,1),diss_avg{i,iii}(1,2),'Color',colSet_blue(t_max,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet_blue(t_max,:), 'MarkerSize',40);
%         errorbar(E_DLCP{1,i}(2,iii),IR_diss_avg{i,iii}(1,1),IR_diss_avg{i,iii}(1,2),'Color',colSet_red(t_max,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);

    end;
    hold off
    axis square;
    t = title(strcat('Dissipation vs. Energy'), 'FontSize',40); %Title the graph
    set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('Energy, $E_{\omega}$ ($eV$)');
    set(x,'Interpreter', 'Latex');
    set(gca,'FontSize',40);
    y = ylabel ('Dissipation, $D$, ($\frac{S}{HzK}$)');
    set(y,'Interpreter','Latex');
    set(gca,'FontSize',40);
    set(gca,'yscale','log');
    hold off
%     xlim([.1 .35]);
%     xlim([255 335]); 
    box on
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('D-E_log_both_',num2str(T(i,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir2,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir2,frameName),'compact');
    close(gcf)
end;

%%

DLN_min = zeros(2,t_max); %array to store the minimum of the dark data
% IR_DLN_min = zeros(2,t_max); %array to store the minimum of the light data

DLN_search = zeros(t_max,f_max); %searchable array for the dark N_{DL} values
% IR_DLN_search = zeros(t_max,f_max); %searchable array for the light N_{DL} values

DLN_search_err = zeros(t_max,f_max); %searchable array for N_{DL} error
% IR_DLN_search_err = zeros(t_max,f_max); %searchable array for the light N_{DL} error

for i = 1:t_max;
    for iii = 1:f_max;
        DLN_search(i,iii) = DLN_tot{i,iii}(1,1); %extract the dark N_{DL} values 
%         IR_DLN_search(i,iii) = IRDLN_tot{i,iii}(1,1); %extract the light N_{DL} values 
        
        
        DLN_search_err(i,iii) = DLN_tot{i,iii}(1,2); %do the same for the error values
%         IR_DLN_search_err(i,iii) = IRDLN_tot{i,iii}(1,2); 
    end;
end;

min_idx1 = zeros(1,t_max); %array to store the index of the frequency at which the minimum occurs in the dark N_{DL} values
% min_idx2 = zeros(1,t_max); %array to store the index of the frequency at which the minimum occurs in the light N_{DL} values

for i = 1:t_max;
    
    [dmin1,min_idx1(1,i)] = min(DLN_search(i,:)); %extract the frequency index of the minimum and the value for the dark data 
%     [dmin2,min_idx2(1,i)] = min(IR_DLN_search(i,:)); %extract the frequency index of the minimum and the value for the light data 
    
    DLN_min(1,i) = dmin1; %store the minimum value
    DLN_min(2,i) = DLN_tot{i,min_idx1(1,i)}(1,2); %extract the standard deviation of the associated minimum N_{DL} value

%     IR_DLN_min(1,i) = dmin2; %do the same for the light data
%     IR_DLN_min(2,1) = IRDLN_tot{i,min_idx2(1,i)}(1,2);
    
end;

N_free = DLN_min';
% IRN_free = IR_DLN_min';

%%

