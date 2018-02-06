%%This program analyzes the output of the modified admittance program. Note
%%that we are now measuring C_{p} and G_{p} directly. 

t_start = input('What is the start measurement temperature? ');
t_stop = input('What is the final measurement temperature? ');
t_step = input('What is the step in between measurement temperatures? ');

v_start = input('What is the start measurement bias? ');
v_stop = input('What is the final measurement bias? ');
v_step = input('What is the step in between measurement biases? ');


dev_number = input('Which device number on the chip was measured? ');


%%
%This section of the code renames the files into a manageable string so
%that we can import the data properly.


for j = t_start:t_step:t_stop %master index over temperature
    for k = 1:(v_stop/v_step)+1; %index over maximum bias V_max
            oldname1 = strcat('dev',num2str(dev_number),'_T' , num2str(j),'K_Vbias', num2str(abs(v_step)-(k*abs(v_step))),'mV_All_Modified.txt'); %identify the file
            newname1 = strcat('T_',num2str(j),'K_', num2str(k),'.txt'); %generate a shorter name. 
            movefile(oldname1,newname1); %rewrite the longer filename with a shorter name defined above.
    end;
%     for k = 4:5;
%         oldname2 = strcat('dev2_T' , num2str(j),'K_Vbias', num2str(k*100),'mV_All_Modified.txt'); %do the same for the positive bias point
%         newname2 = strcat('T_',num2str(j),'K_', num2str(k),'.txt' ); 
%         movefile(oldname2,newname2);  
%     end;
end;



for j = t_start:t_step:t_stop %master index over temperature
    for k = 1:(v_stop/v_step)+1; %index over maximum bias V_max
            oldname1 = strcat('dev',num2str(dev_number),'_IR__T' , num2str(j),'K_Vbias', num2str(abs(v_step)-(k*abs(v_step))),'mV_All_Modified.txt'); %identify the file
            newname1 = strcat('IR_T_',num2str(j),'K_', num2str(k),'.txt'); %generate a shorter name. 
            movefile(oldname1,newname1); %rewrite the longer filename with a shorter name defined above.
    end;
%     for k = 4:5;
%         oldname2 = strcat('dev2_IR__T' , num2str(j),'K_Vbias', num2str(k*100),'mV_All_Modified.txt'); %do the same for the positive bias point
%         newname2 = strcat('IR_T_',num2str(j),'K_', num2str(k),'.txt' ); 
%         movefile(oldname2,newname2);  
%     end;
end;

%%

%This section stores all the data from the experiment into 1 cell array
%and configures the algorithm with proper constants, library structures,
%and paths to save figures

AS = cell(((t_stop-t_start)/t_step)+1,(v_stop/v_step)+1);
IRAS = cell(((t_stop-t_start)/t_step)+1,(v_stop/v_step)+1);

t_max = size(AS,1);
v_max = size(AS,2);

for j = 1:t_max; %index over temperature
    for k = 1:v_max; %index over V_{bias}
            
            filename1 = strcat('T_',num2str(t_start+(j-1)*t_step),'K_', num2str(k),'.txt'); %create for each 0 bias temperature point
            delimiter1 = '\t';
            A = tdfread(filename1,delimiter1); %Import dataset into the active workspace; assign columns in each array 
            AS{j,k}(:,1) = A.Frequency_rads;
            AS{j,k}(:,2) = A.Cp_Farad;
            AS{j,k}(:,3) = A.Bp_S;
            AS{j,k}(:,4) = A.Dissipation;
            AS{j,k}(:,5) = A.Gp_S;
            
%             filename2 = strcat('IR_T_',num2str(t_start+(j-1)*t_step),'K_', num2str(k),'.txt'); %create for each 0 bias temperature point
%             delimiter2 = '\t';
%             A = tdfread(filename2,delimiter2); %Import dataset into the active workspace; assign columns in each array 
%             IRAS{j,k}(:,1) = A.Frequency_rads;
%             IRAS{j,k}(:,2) = A.Cp_Farad;
%             IRAS{j,k}(:,3) = A.Bp_S;
%             IRAS{j,k}(:,4) = A.Dissipation;
%             IRAS{j,k}(:,5) = A.Gp_S;
            
    end;
end;


colSet_red = makeColorMap([217.694 207.1875 207.1875]/255, [191 0 0]/255, t_max); %make a red colormap for soaked data 
colSet_blue = makeColorMap([207.188 207.18 214.889]/255, [53 55 160]/255, t_max); %make a blue colormap for the reference data (if applicable)
colSet_purple = makeColorMap([210 201 209]/255, [87 42 88]/255, t_max); %make a purple colormap for the low intensity light

field1 = 'index_1'; %create a structure to map the voltage parameter to the corresponding data set.
value1 = 0;
field2 = 'index_2';
value2 = -100;
field3 = 'index_3';
value3 = -200;
% field4 = 'index_4';
% value4 = 400;
% field5 = 'index_5';
% value5 = 500;
B = struct(field1,value1,field2,value2,field3,value3); %create the structure
BNames = fieldnames(B); %create a cell array with the field names so that we can call on them inside the plot for-loop

T = zeros(t_max,1); %store the experiment temperature as an array
for k = 1:t_max;
    T(k,1) = t_start+(k-1)*t_step;
end;

A = 4e-6;  %Define constants that will be used later in our looped calculation. Note that the semicondcutor radius is in m, as are the device area (m^2) and volume (m^3). Energy parameters are in eV.  
epsilon = 20;
eps_0 = 8.8541878e-12;
q = 1.60217657e-19;
kB = 8.6173324*(10^-5);  


figuresdir1 = 'C:\Users\Eric\Desktop\Admittance\ZnO-PbS\light_soak_study\IR_soak\sanity_checks\real_cell\light_soak_temp\Analysis\TAS';

%%

for k = 1:t_max;
    for ii = 1:v_max;
        AS{k,ii}(~any(isnan(AS{k,ii}),2),:);
    end;
end;

%%
% This section calculates the derivatives of the data by using the diff
% command and stores it in the original data array in a new column in each
% cell

for k = 1:t_max; %index over temperatures
    for ii = 1:v_max; %index over voltages
        
            AS{k,ii}(:,6) = smooth(AS{k,ii}(:,2),9,'sgolay');%smooth the data with a Savitzky-Golay Filter
            AS{k,ii+v_max}(:,1) = diff(AS{k,ii}(:,1)); %calculate difference in the frequency data across all temperature points and biases
            AS{k,ii+v_max}(:,2) = diff(AS{k,ii}(:,6)); %calculate the difference in the capacitance data across all temperature points and biases
            AS{k,ii+v_max}(:,3) = AS{k,ii+v_max}(:,2)./AS{k,ii+v_max}(:,1); %calculate the numnerical derivative
            AS{k,ii+v_max}(:,4) = -AS{k,ii}(2:end,1).*AS{k,ii+v_max}(:,3);  %normalize the differential capacitance by frequency. Note that the indices in the frequency variable have shifted. 
            
%             IRAS{k,ii}(:,6) = smooth(IRAS{k,ii}(:,2),9,'sgolay');%smooth the data with a Savitzky-Golay Filter
%             IRAS{k,ii+v_max}(:,1) = diff(IRAS{k,ii}(:,1)); %calculate difference in the frequency data across all temperature points and biases
%             IRAS{k,ii+v_max}(:,2) = diff(IRAS{k,ii}(:,6)); %calculate the difference in the capacitance data across all temperature points and biases
%             IRAS{k,ii+v_max}(:,3) = IRAS{k,ii+v_max}(:,2)./IRAS{k,ii+v_max}(:,1); %calculate the numnerical derivative
%             IRAS{k,ii+v_max}(:,4) = -IRAS{k,ii}(2:end,1).*IRAS{k,ii+v_max}(:,3);  %normalize the differential capacitance by frequency. Note that the indices in the frequency variable have shifted. 
                     
    end;
end;


%%
%This section of the code plots the raw output of the experiment. The first
%block plots capacitance, the second plots conductance, and the final block
%is for dissipation.

%plot the capacitance as a function of frequency and save it to the proper
%path

for ii = 1:v_max; %index across V_{bias}
    figure(); %Plot the capacitance data grouped according to bias level to look for trends across this parameter
    hold on;
    for k = 1:t_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(AS{k,ii}(:,1),AS{k,ii}(:,6),'LineWidth',3,'Color',colSet_blue(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('C vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel('Capacitance, C (Farads)', 'FontSize', 40);
    set(y,'FontName','Calibri');
    xlim([3e2 10e6]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([.5e-8 1.5e-8]);    
    hold off;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('C_ref_',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)
end;
%%
for ii = 1:v_max; %index across V_{bias}
    figure(); %Plot the capacitance data grouped according to bias level to look for trends across this parameter
    hold on;
    for k = 1:t_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(IRAS{k,ii}(:,1),IRAS{k,ii}(:,6),'LineWidth',3,'Color',colSet_red(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('C vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel('Capacitance, C (Farads)', 'FontSize', 40);
    set(y,'FontName','Calibri');
    xlim([3e2 10e6]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([.5e-8 1.5e-8]);    
    hold off;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('C_IR_',num2str(B.(BNames{ii,1})),'mV')];
    print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end;

%%
%next, plot the conductance as a function of frequency and store it to the
%proper path

for ii = 1:v_max; %index across V_{bias}
    figure(); %Plot the conductance data grouped according to bias level to look for trends across this parameter
    hold on;
    for k = 1:t_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(AS{k,ii}(:,1),AS{k,ii}(:,5),'LineWidth',3,'Color',colSet_blue(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('G vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel('Conductance, G (Siemens)', 'FontSize', 40);
    set(y,'FontName','Calibri');
    xlim([3.0e2 10e6]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([0 .01]);    
    hold off;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('G_ref_',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)
end;

%next, plot the dissipation as a function of frequency and store it to the
%proper path

for ii = 1:v_max; %index across V_{bias}
    figure(); %Plot the dissipation data grouped according to bias level to look for trends across this parameter
    hold on;
    for k = 1:t_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(AS{k,ii}(:,1),AS{k,ii}(:,4),'LineWidth',3,'Color',colSet_blue(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('D vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 4); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',48);
    y = ylabel('Dissipation, D', 'FontSize', 48);
    set(y,'FontName','Calibri');
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
    xlim([3.0e2 10e6]);
%     ylim([0 4]);    
%     hold off;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('D_ref_',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)
end;

%%
for ii = 1:v_max; %index across V_{bias}
    figure(); %Plot the conductance data grouped according to bias level to look for trends across this parameter
    hold on;
    for k = 1:t_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(IRAS{k,ii}(:,1),IRAS{k,ii}(:,5),'LineWidth',3,'Color',colSet_red(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('G vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel('Conductance, G (Siemens)', 'FontSize', 40);
    set(y,'FontName','Calibri');
    xlim([3.0e2 10e6]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([0 .01]);    
    hold off;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('G_IR_',num2str(B.(BNames{ii,1})),'mV')];
    print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end;

%next, plot the dissipation as a function of frequency and store it to the
%proper path

for ii = 1:v_max; %index across V_{bias}
    figure(); %Plot the dissipation data grouped according to bias level to look for trends across this parameter
    hold on;
    for k = 1:t_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(IRAS{k,ii}(:,1),IRAS{k,ii}(:,4),'LineWidth',3,'Color',colSet_red(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('D vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel('Dissipation, D', 'FontSize', 40);
    set(y,'FontName','Calibri');
    set(gca,'xscale','log');
    xlim([3.0e2 10e6]);
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([0 4]);    
    hold off;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('D_IR_',num2str(B.(BNames{ii,1})),'mV')];
    print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end;

%%       
%This section plots the smoothed differential capacitance as a function of frequency
%so that we can extract the activation energy E_t.

%This script differentiates the capacitace-frequency data from the TAS
%experiment. It then fits the peaks in the derivative to a Gaussian fit,
%and calculates the error in the fit. The derived resonant frequency is
%then fit to an Arrhenius plot from which the trap energy and attempt to
%escape frequency are determined


diffx = cell(v_max,t_max); %create a cell to hold the smoothed and differentiated x-values
diffy = cell(v_max,t_max); %create a cell to hold the smoothed and differentiated y-values
ASderiv = cell(v_max,t_max); %create a cell to hold the functional derivative
ASderiv_norm = cell(v_max,t_max); %create a cell to store the values of the normalized derivative

IRdiffx = cell(v_max,t_max); %do the same for the illuminated values
IRdiffy = cell(v_max,t_max);
IRASderiv = cell(v_max,t_max);
IRASderiv_norm = cell(v_max,t_max);


for k = 1:t_max; %index over temperatures
    for ii = 1:v_max; %index over voltages
            diffx{ii,k}(:,1) = filter(smooth_diff(10),1,AS{k,ii}(:,1)); %differentiate and smooth the frequency
            diffy{ii,k}(:,1) = filter(smooth_diff(10),1,AS{k,ii}(:,6)); %differentiate and smooth the capacitance
            ASderiv{ii,k}(:,1) = diffy{ii,k}(:,1)./diffx{ii,k}(:,1); %calculate the functional derivative
            ASderiv_norm{ii,k}(:,1) = -AS{k,ii}(1:end,1).*ASderiv{ii,k}(:,1);  %normalize the differential capacitance by frequency. Note that the indices in the frequency variable have shifted. 
    
%             IRdiffx{ii,k}(:,1) = filter(smooth_diff(10),1,IRAS{k,ii}(:,1)); %do the same for the illuminated values
%             IRdiffy{ii,k}(:,1) = filter(smooth_diff(10),1,IRAS{k,ii}(:,6)); 
%             IRASderiv{ii,k}(:,1) = IRdiffy{ii,k}(:,1)./IRdiffx{ii,k}(:,1); 
%             IRASderiv_norm{ii,k}(:,1) = -IRAS{k,ii}(1:end,1).*IRASderiv{ii,k}(:,1);   
    end;
end;

%%


for ii = 1:v_max
        figure();%Plot the derivative data according to bias level to look for trends accross this parameter
        hold on;
        for k = 1:t_max;
            plot(AS{k,ii}(1:end,1),ASderiv_norm{ii,k}(:,1),'LineWidth',3, 'Color', colSet_blue(k,:)); %Plot the remaining temperatures and edit the plot accordingly
        end;
        axis square
        box on;
%         t = title (strcat('$\omega\frac{dC}{d\omega}$ vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',20);
%         set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%         set(t, 'units', 'normalized')
        x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
        set(x,'FontName','Latex');
        set(gca,'FontSize',40);
        y = ylabel('\omega^{dC}/_{d \omega} (Farads)','FontSize', 40);
        set(y,'FontName','Calibri');
        set(gca,'xscale','log');
        set(gca,'XMinorTick','on');
        xlim([2e3 10e6]);
        set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%         xlim([2092 4.827e4]);
%         ylim([0 5e-9]); 
        hold off;
%         set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%         frameName = [strcat('smoothed_derivative_',num2str(B.(BNames{ii,1})),'mV')];
%         print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%         savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%         close(gcf)
end;

%%

for ii = 1:v_max
        figure(); %Plot the derivative data according to bias level to look for trends accross this parameter
        hold on;
        for k = 1:t_max;
            plot(IRAS{k,ii}(1:end,1),IRASderiv_norm{ii,k}(:,1),'LineWidth',3, 'Color', colSet_red(k,:)); %Plot the remaining temperatures and edit the plot accordingly
        end;
        axis square
        box on;
%         t = title (strcat('$\omega\frac{dC}{d\omega}$ vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',20);
%         set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%         set(t, 'units', 'normalized')
        x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
        set(x,'FontName','Calibri');
        set(gca,'FontSize',40);
        y = ylabel('\omega^{dC}/_{d\omega} (Farads)','FontSize', 40);
        set(y,'FontName','Calibri');
        set(gca,'xscale','log');
        set(gca,'XMinorTick','on');
        xlim([2e3 10e6]);
        set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%         xlim([2092 4.827e4]);
%         ylim([0 5e-9]); 
        hold off;
        set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
        frameName = [strcat('smoothed_derivative_IR_',num2str(B.(BNames{ii,1})),'mV')];
        print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
        savefig(gcf, fullfile(figuresdir1,frameName),'compact');
        close(gcf)
end;



%%

diff_fit = cell(v_max,t_max);

fit_round = zeros(size(AS{3,1},1),3); %create an array to round the frequency values so that we can read the relevant frequency range directly from the plot
fit_round(:,1) = round(AS{3,1}(:,1),-0);
fit_round(:,2) = round(AS{3,1}(:,1),-1);
fit_round(:,3) = round(AS{3,1}(:,1),-2);
fit_round(:,4) = round(AS{3,1}(:,1),-3);

fitidx = zeros(2,2); %array to store the bounds of the frequency range of interest
fitidx(1,1) = find(fit_round(:,1) == 3405); %input the first boundary condition
fitidx(1,2) = find(fit_round(:,4) == 2.13e6); %input the final boundary condition

fitidx(2,1) = find(fit_round(:,1) == 7672); %do the same for the illuminated data
fitidx(2,2) = find(fit_round(:,4) == 1.539e6);

ASpeaks = zeros(t_max,v_max); %create arrays to store the peak frequencies as fit below
IRASpeaks = zeros(t_max,v_max);

sigma = zeros(t_max,v_max);
IRsigma = zeros(t_max,v_max);

gauss_err = zeros(t_max,v_max); %create arrays to store the error associated with each gaussian fit 
IRgauss_err = zeros(t_max,v_max); 

%%
pts = input('How many points on either side of the max should be fit as a Gaussian for the second fit? ');

for k = 1:t_max;
    for ii = 1:v_max;
    
    f1 = fit(log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),ASderiv_norm{ii,k}(fitidx(1,1):fitidx(1,2),1),'gauss1'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,log(AS{i,ii}(fitidx(1,1):fitidx(1,2),1)),ASderiv_norm{ii,i}(fitidx(1,1):fitidx(1,2),1)); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);

    val = exp(p1(2)); %extract the peak fit value
    tmp = (AS{k,ii}(fitidx(1,1):fitidx(1,2),1)-val); %subtract the peak value from the frequency range in the experiment
    [idx idx] = min(abs(tmp)); %find the minimum, picking out the peak in the experimental data
    g_idx = round((fitidx(1,1) + idx)-1,0); %find the index of the peak in the global data range
    closest = AS{k,ii}(g_idx,1); %extract the frequency of the peak

    X_sub = AS{k,ii}((g_idx-pts):(g_idx+pts),1); %extract a sub-range of x-data
    Y_sub = ASderiv_norm{ii,k}((g_idx-pts):(g_idx+pts),1); %extract a sub-range of y-data
    Y_max = max(Y_sub); %extract the max of the y-range
    Y_sub_prime = Y_sub./Y_max; %normalize the y-data against the max
    
    options = fitoptions('gauss1','StartPoint', [ASderiv_norm{ii,k}(g_idx,1), log(closest), log(AS{k,ii}((g_idx+pts),1))-log(AS{k,ii}((g_idx-pts),1))]); %using the first fit, input fit parameters for the second fit
    [f1a,test] = fit(log(X_sub),Y_sub_prime,'gauss1',options); %perform the second fit
    p1a = coeffvalues(f1a); %extract the coefficients of the fit
    e1a = confint(f1a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    ASpeaks(k,ii) = exp(p1a(2)); %obtain the centroid of the gaussian
    sigma(k,ii) =  (exp(e1a(2,2))-exp(e1a(1,2)))/2; %calculate the std deviation associated with the centroid fit in a 'x +/- y' format, with y as the std deviation
    gauss_err(k,ii) = sigma(k,ii)/ASpeaks(k,ii); %normalize the std deviation by the mean
    
    figure()
    hold on;
    plot(f1a,log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),(ASderiv_norm{ii,k}(fitidx(1,1):fitidx(1,2),1)./Y_max)); 
    axis square
    t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
%     plot(f1a,log(AS{i,ii}(fitidx(1,1):fitidx(1,2),1)),(ASderiv_norm{ii,i}(fitidx(1,1):fitidx(1,2),1)./Y_max),'residuals')
    
    
    %do the same for the IR-illuminated data
    
    f2 = fit(log(IRAS{k,ii}(fitidx(2,1):fitidx(2,2),1)),IRASderiv_norm{ii,k}(fitidx(2,1):fitidx(2,2),1),'gauss1'); %do the same for the illuminated data 
    p2 = coeffvalues(f2); 

%     figure()
%     hold on;
%     plot(f2,log(IRAS{i,ii}(fitidx(2,1):fitidx(2,2),1)),IRASderiv_norm{ii,i}(fitidx(2,1):fitidx(2,2),1));
%     t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV IR,1'), 'FontSize',30);
    
    IR_val = exp(p2(2)); %extract the peak fit value
    IR_tmp = (IRAS{k,ii}(fitidx(2,1):fitidx(2,2),1)-IR_val); %subtract the peak value from the frequency range in the experiment
    [idx idx] = min(abs(IR_tmp)); %find the minimum, picking out the peak in the experimental data
    IR_g_idx = round((fitidx(2,1) + idx)-1,0); %find the index of the peak in the global data range
    IR_closest = IRAS{k,ii}(IR_g_idx,1); %extract the frequency of the peak

    IR_X_sub = IRAS{k,ii}((IR_g_idx-pts):(IR_g_idx+pts),1); %extract a sub-range of x-data
    IR_Y_sub = IRASderiv_norm{ii,k}((IR_g_idx-pts):(IR_g_idx+pts),1); %extract a sub-range of y-data
    IR_Y_max = max(IR_Y_sub); %extract the max of the y-range
    IR_Y_sub_prime = IR_Y_sub./IR_Y_max; %normalize the y-data against the max
    
    IR_options = fitoptions('gauss1','StartPoint', [IRASderiv_norm{ii,k}(g_idx,1), log(IR_closest), log(IRAS{k,ii}((IR_g_idx+pts),1))-log(IRAS{k,ii}((IR_g_idx-pts),1))]); %using the first fit, input fit parameters for the second fit
    f2a = fit(log(IR_X_sub),IR_Y_sub_prime,'gauss1',IR_options); %perform the second fit
    p2a = coeffvalues(f2a); %extract the coefficients of the fit
    e2a = confint(f2a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    IRASpeaks(k,ii) = exp(p2a(2)); %obtain the centroid of the gaussian
    IRsigma(k,ii) =  (exp(e2a(2,2))-exp(e2a(1,2)))/2; %calculate the std deviation associated with the centroid fit in a 'x +/- y' format, with y as the std deviation
    IRgauss_err(k,ii) = IRsigma(k,ii)/IRASpeaks(k,ii); %normalize the std deviation by the mean
    
    figure()
    hold on;
    plot(f2a,log(IRAS{k,ii}(fitidx(2,1):fitidx(2,2),1)),(IRASderiv_norm{ii,k}(fitidx(2,1):fitidx(2,2),1)./IR_Y_max)); 
    axis square
    t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV IR,2'), 'FontSize',30);
       
    end;
end;

%%
offset_ref = zeros(v_max,2); %create an array to store the offset values for each voltage range. The first column represents the offset in the temperature curves to start looking for peaks.
offset_IR = zeros(v_max,2);

for ii = 1:v_max; %index across voltages
    offset_ref(ii,1) = 0; %specify the offset in curves at which to start analysis 
    offset_ref(ii,2) = t_max; %specify curve at which to stop analysis
    
    offset_IR(ii,1) = 0;
    offset_IR(ii,2) = t_max;
end;

AS_inv = cell(1,2*v_max); %create a cell to hold the 1/(\omega) and ln(1/T)
AS_fits = cell(1,v_max); %create a cell to hold the fit parameters
AS_line = cell(1,v_max); %create a cell to hold the calculated best fit line
E_t = zeros(2,v_max); %create an array to store the trap energy (slope)
nu_0 = zeros(2,v_max); %create an array to store the reduced attempt escape frequency 
nu_00 = zeros(t_max,v_max); %create an array to store the temperature-dependent attempt to escape frequency



IRAS_inv = cell(1,2*v_max); %create a cell to hold the 1/(\omega) and ln(1/T)
IRAS_fits = cell(1,v_max); %create a cell to hold the fit parameters
IRAS_line = cell(1,v_max); %create a cell to hold the calculated best fit line
IR_E_t = zeros(2,v_max); %create an array to store the trap energy (slope)
IR_nu_0 = zeros(2,v_max); %create an array to store the reduced attempt escape frequency 
IR_nu_00 = zeros(t_max,v_max); %create an array to store the temperature-dependent attempt to escape frequency


%%

for k = 1:t_max; %index over temperatures
    for ii = 1:v_max; %index over voltages
        AS_inv{1,ii}(k+offset_ref(ii,1),1) = 1./T(k+offset_ref(ii,1),1); %calculate inverse temp
        AS_inv{1,ii}(k+offset_ref(ii,1),2) = reallog(1./ASpeaks(k+offset_ref(ii,1),ii)); %calculate ln(1/(\omega))
    end;
    if k >= (offset_ref(ii,2)-offset_ref(ii,1));
        break;
    end;
end;


for k = 1:t_max; %index over temperatures
    for ii = 1:v_max; %index over voltages
        IRAS_inv{1,ii}(k+offset_IR(ii,1),1) = 1./T(k+offset_IR(ii,1),1); %calculate inverse temp
        IRAS_inv{1,ii}(k+offset_IR(ii,1),2) = reallog(1./IRASpeaks(k+offset_IR(ii,1),ii)); %calculate ln(1/(\omega))
    end;
    if k >= (offset_IR(ii,2)-offset_IR(ii,1));
        break;
    end;
end;

%%

for ii = 1:v_max; %index over voltages
        temp1 = AS_inv{1,ii}(:,1); %create a mask that removes zeros from the cell array   
        temp2 = AS_inv{1,ii}(:,2);
        temp1 = temp1(temp1>0);
        temp2 = temp2(abs(temp2)>0);
        AS_inv{1,ii+v_max}(:,1) = temp1; %create a mask that removes zeros from ln(1/f)
        AS_inv{1,ii+v_max}(:,2) = temp2; %apply the mask to 1/T 
        
        temp3 = IRAS_inv{1,ii}(:,1); %create a mask that removes zeros from the cell array   
        temp4 = IRAS_inv{1,ii}(:,2);
        temp3 = temp3(temp3>0);
        temp4 = temp4(abs(temp4)>0);
        IRAS_inv{1,ii+v_max}(:,1) = temp3; %create a mask that removes zeros from ln(1/f)
        IRAS_inv{1,ii+v_max}(:,2) = temp4; %apply the mask to 1/T 
     
end;

%%

X = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store the x-values for the Arrenhius plot
Y = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store the y-values for the Arrhenius plot
Yhat = zeros(100,v_max); %create an array to store the predicted y-values for the plot
Yhat_resid = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store the predicted y-values for the residual calculation
beta = zeros(2,v_max); %crete an array to store the fit parameters
errbeta = cell(1,v_max); %create an array to store the covariance matrix associated with each fit
W = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store the inverse variance for weights
Wprime = zeros(offset_ref(1,2)-offset_ref(1,1),1); %create a temporary array to store the weights for the particular fit being performed
fiterr = zeros(2,v_max); %create an array to calculate the error for each paramter
resid = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store the residuals of each fit
unit_test1 = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store sum of the square of the residulas divided by the number of points
unit_test2 = zeros(offset_ref(1,2)-offset_ref(1,1),v_max); %create an array to store the sum of the variances for each data point

X_vals = linspace(1/200, 1/350, 100); %create an array of x-values for the prediction line

for ii = 1:v_max; %index over voltages
    
    X(:,ii) = AS_inv{1,ii+v_max}(:,1); %extract the inverse temperatures
    Y(:,ii) = AS_inv{1,ii+v_max}(:,2); %extract ln(1/freq)
    Xtilde = [ones(length(X(:,ii)),1),X(:,ii)]; %create an array of ones and the x-values  for the fit specified by bias 'ii' 
    W(:,ii) = 1; %give the experimental data points equal weight, assuming that the error in the trend is greater than the error in each individual measurement
%     W(:,ii) = (gauss_err(1+offset_ref(1,1):offset_ref(1,2),ii).^2); %create the weighting array. Note that we are interested in 1/sigma^2 for each data point
    Wprime(:,1) = W(:,ii); %extract the variances for the fit being run
    
    beta(:,ii) =  (Xtilde.' * diag(Wprime)^(-1) * Xtilde)^(-1) * (Xtilde.' * diag(Wprime)^(-1) * Y(:,ii)); %calculate the regression coefficients
    Yhat(:,ii) = X_vals(1,:).*beta(2,ii) + beta(1,ii); %generate the array of predicted values for the plot

    Yhat_resid(:,ii) = X(:,ii).*beta(2,ii) + beta(1,ii); %generate the array of predicted values for the residuals calculation
    
    resid(:,ii) = Y(:,ii) - Yhat_resid(:,ii); %calculate the residuals
    ei_sq = mean(resid(:,ii).^2);
    
    errbeta{1,ii} = ei_sq.*(Xtilde.' *diag(Wprime)^(-1)* Xtilde)^(-1); %calculate the variance in the regression coefficients
%     errbeta{1,ii} = (Xtilde.' *diag(Wprime)^(-1)* Xtilde)^(-1); %calculate the variance in the regression coefficients
    fiterr(1,ii) = sqrt(errbeta{1,ii}(1,1)/(offset_ref(1,2)-offset_ref(1,1))); %normalize the error in the constant properly
    fiterr(2,ii) = sqrt(errbeta{1,ii}(2,2)/(offset_ref(1,2)-offset_ref(1,1))); %normalize the error in the slope properly
     
    E_t(1,ii) = kB*beta(2,ii); %calculate the trap energy
    E_t(2,ii) = kB*fiterr(2,ii); %calculate the 1-sigma error in the energy
    
    nu_0(1,ii) = exp(abs(beta(1,ii)))/2; %calculate the attempt to escape frequency
    nu_0(2,ii) = nu_0(1,ii)*fiterr(1,ii)/2; %calculate the 1-sigma error in the attempt to escape frequency
    
  
    unit_test1(:,ii) = (sum(resid(:,ii).^2)/length(Y)); %calculate the sum of the square of the residuals
    unit_test2(:,ii) = (sum(W(:,ii))/length(Y)); %calculate the sum of the variances of the data points
    
    %note that the unit tests should be of the same order of magnitude if W
    %is specified correctly
    
end;  

%%

%do the same for the IR data

IR_X = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store the x-values for the Arrenhius plot
IR_Y = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store the y-values for the Arrhenius plot
IR_Yhat = zeros(100,v_max);
IR_Yhat_resid = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store the predicted y-values for the fit
IR_beta = zeros(2,v_max); %crete an array to store the fit parameters
IR_errbeta = cell(1,v_max); %create an array to store the covariance matrix associated with each fit
IR_W = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store the inverse variance for weights
IR_Wprime = zeros(offset_IR(1,2)-offset_IR(1,1),1); %create a temporary array to store the weights for the particular fit being performed
IR_fiterr = zeros(2,v_max); %create an array to calculate the error for each paramter
IR_resid = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store the residuals of each fit
IR_unit_test1 = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store sum of the square of the residulas divided by the number of points
IR_unit_test2 = zeros(offset_IR(1,2)-offset_IR(1,1),v_max); %create an array to store the sum of the variances for each data point

IR_X_vals = linspace(1/200, 1/350, 100); %create an array of x-values for the prediction line

for ii = 1:v_max; %index over voltages
    
    IR_X(:,ii) = IRAS_inv{1,ii+v_max}(:,1); %extract the inverse temperatures
    IR_Y(:,ii) = IRAS_inv{1,ii+v_max}(:,2); %extract ln(1/freq)
    IR_Xtilde = [ones(length(IR_X(:,ii)),1),IR_X(:,ii)]; %create an array of ones and the x-values  for the fit specified by bias 'ii' 
    IR_W(:,ii) = 1; %give the experimental data points equal weight, assuming that the error in the trend is greater than the error in each individual measurement
%     IR_W(:,ii) = (IRgauss_err(1+offset_IR(1,1):offset_IR(1,2),ii).^2); %create the weighting array. Note that we are interested in 1/sigma^2 for each data point
    IR_Wprime(:,1) = IR_W(:,ii); %extract the variances for the fit being run
    
    IR_beta(:,ii) =  (IR_Xtilde.' * diag(IR_Wprime)^(-1) * IR_Xtilde)^(-1) * (IR_Xtilde.' * diag(IR_Wprime)^(-1) * IR_Y(:,ii)); %calculate the regression coefficients
    IR_Yhat(:,ii) = IR_X_vals(1,:).*IR_beta(2,ii) + IR_beta(1,ii); %generate the array of predicted values
    
    IR_Yhat_resid(:,ii) = IR_X(:,ii).*beta(2,ii) + beta(1,ii); %calculate the residuals
    
    IR_resid(:,ii) = IR_Y(:,ii) - IR_Yhat_resid(:,ii); %calculate the residuals
    IR_ei_sq = mean(resid(:,ii).^2); %find the mean of the square of the residuals
    
    IR_errbeta{1,ii} = IR_ei_sq.*(IR_Xtilde.' *diag(IR_Wprime)^(-1)* IR_Xtilde)^(-1); %calculate the variance in the regression coefficients
%     IR_errbeta{1,ii} = (IR_Xtilde.' *diag(IR_Wprime)^(-1)* IR_Xtilde)^(-1); %calculate the variance in the regression coefficients
    IR_fiterr(1,ii) = sqrt(IR_errbeta{1,ii}(1,1)/(offset_IR(1,2)-offset_IR(1,1))); %normalize the error in the constant properly
    IR_fiterr(2,ii) = sqrt(IR_errbeta{1,ii}(2,2)/(offset_IR(1,2)-offset_IR(1,1))); %normalize the error in the slope properly
   
    IR_E_t(1,ii) = kB*IR_beta(2,ii); %calculate the trap energy
    IR_E_t(2,ii) = kB*IR_fiterr(2,ii); %calculate the 1-sigma error in the energy.00
    
    IR_nu_0(1,ii) = exp(abs(IR_beta(1,ii)))/2; %calculate the attempt to escape frequency
    IR_nu_0(2,ii) = IR_nu_0(1,ii)*fiterr(1,ii)/2; %calculate the 1-sigma error in the attempt to escape frequency
    
    IR_unit_test1(:,ii) = (sum(IR_resid(:,ii).^2)/length(IR_Y)); %calculate the sum of the square of the residuals
    IR_unit_test2(:,ii) = (sum(IR_W(:,ii))/length(IR_Y)); %calculate the sum of the variances of the data points
    
    %note that the unit tests should be of the same order of magnitude if W
    %is specified correctly
    
end;  


%%


for ii = 1:v_max; %indices for the voltage range
    figure() %create the Arrhenius plot
    g1 = errorbar(AS_inv{1,ii+v_max}(:,1),AS_inv{1,ii+v_max}(:,2),gauss_err(1+offset_ref(ii,1):offset_ref(ii,2),ii),'.','color',colSet_blue(t_max,:)); %create the scatter plot of ln(1/f) vs. 1/T
    hold on;
    g2 = errorbar(IRAS_inv{1,ii+v_max}(:,1),IRAS_inv{1,ii+v_max}(:,2),IRgauss_err(1+offset_IR(ii,1):offset_IR(ii,2),ii),'.','color',colSet_red(t_max,:));
        set(g1,{'markers'},{50},{'Linewidth'},{4});
        set(g2,{'markers'},{50},{'Linewidth'},{4});
    plot(X_vals(1,:),Yhat(:,ii),'--k','LineWidth', 3); %plot the best fit line over it
    plot(IR_X_vals(1,:),IR_Yhat(:,ii),':k','LineWidth', 3)
    
    
    axis square;
    box on;
    t = title (strcat('Arrhenius Plot,', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize', 40);
    set(t,'Interpreter', 'Latex')
    x = xlabel ('$T^{-1} (K^{-1})$', 'FontSize', 40);
    set(x,'Interpreter','Latex');
    set(gca,'FontSize',40);
    y = ylabel ('ln($\frac{1}{\omega_0}$)', 'FontSize',40);
    set(y,'Interpreter', 'Latex');
    hold off
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Arrhenius_',num2str(B.(BNames{ii,1})),'mV')];
    print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir1,frameName),'compact');
    close(gcf)
end;

