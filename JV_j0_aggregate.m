%This script analyzes the dark J-V characteristics and attempts to relate
%the saturation current of the device to the trap state detected in the TAS
%and DLCP measurements

figuresdir7 = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\09_2017_m1\analysis\figures\JV\new';

% d_max = 3;

A_cm = 4e-2;

% colSet_red = makeColorMap([217.694 207.1875 207.1875]/255, [191 0 0]/255, t_max); %make a red colormap for soaked data 
% colSet_blue = makeColorMap([207.188 207.18 214.889]/255, [53 55 160]/255, t_max); %make a blue colormap for the reference data (if applicable)
% colSet_purple = makeColorMap([210 201 209]/255, [87 42 88]/255, t_max); %make a purple colormap for the low intensity light
% colSet_black = makeColorMap([226 226 226]/255, [0 0 0]/255, t_max); %make a black colormap


for j = t_start:t_step:t_stop %master index over temperature
    
    oldname1 = strcat('dev',num2str(dev_number),'_IVDark_' , num2str(j),'K.txt'); %identify the file
    newname1 = strcat('T_',num2str(j),'K_D.txt'); %generate a shorter name. 
    movefile(oldname1,newname1); %rewrite the longer filename with a shorter name defined above.
    
    oldname2 = strcat('dev',num2str(dev_number),'_IVLight_' , num2str(j),'K.txt'); %identify the file
    newname2 = strcat('T_',num2str(j),'K_L.txt'); %generate a shorter name. 
    movefile(oldname2,newname2); %rewrite the longer filename with a shorter name defined above.
    
end

%%
IVD = zeros(48,t_max+1);
IVL = zeros(48,t_max+1);
jvm = size(IVD,1);

for j = 1:t_max %index over temperature
    
            filename1 = strcat('T_',num2str(t_start+t_step*(j-1)),'K_D.txt'); 
            delimiter = ',';
            D = tdfread(filename1,delimiter); %import the dark data into the active workspace; assign columns in each array 
            IVD(:,1) = D.VDS;
            IVD(:,j+1) = D.ID;
            
            filename2 = strcat('T_',num2str(t_start+t_step*(j-1)),'K_L.txt'); 
            delimiter = ',';
            L = tdfread(filename2,delimiter); %import the illuminated into the active workspace; assign columns in each array 
            IVL(:,1) = L.VDS;
            IVL(:,j+1) = L.ID;          
end

JVD = zeros(size(IVD)); 
JVL = zeros(size(IVL));

for k = 1:t_max
    JVD(:,1) = IVD(:,1);
    JVD(:,k+1) = (IVD(:,k+1).*1e3)./A_cm;
    
    JVL(:,1) = IVL(:,1);
    JVL(:,k+1) = (IVL(:,k+1).*1e3)./A_cm;
end


jsc = zeros(t_max,1); %create an array to store the short-circuit current values at each temperature

for k = 1:t_max
    xq = -.01:.001:.01; %define a test interval of voltage values to interpolate between
    vq = interp1(JVD(:,1),JVD(:,k+1),xq); %for each curve, interpolate the curve to find the current at short circuit conditions 
    jsc_idx = find(xq == 0); %extract the index of the interpolated 0V point
    jsc(k,1) = vq(1,jsc_idx);%extract the jsc value
end


%%

%This section stores all the data from the experiment into 1 cell array
%and configures the script with proper constants, library structures,
%and paths to save figures


figure
hold on
for k = 2:t_max+1

    h2 = plot(JVD(:,1),abs(JVD(:,k)),'Color',colSet_green(k-1,:),'LineWidth',2,'Marker','o','MarkerSize',20,'MarkerFaceColor',colSet_green(k-1,:),'MarkerEdgecolor','k');
%     h3 = plot(JVL(:,1),abs(JVL(:,k)),'Color',colSet_blue(k-1,:),'LineWidth',2,'Marker','o','MarkerSize',20,'MarkerFaceColor',colSet_blue(k-1,:),'MarkerEdgecolor','k');


end
axis square
box on

x = xlabel ('V (volt)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
xlim([-.3 .8])
% set(gca,'xtick',[-.2 -.1 -.0 .1 .2])

y = ylabel ('J (^{mA}/_{cm^2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',46);
ylim([1e-3 1e2])
set(gca,'yscale','log');
set(gca,'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2])

pbaspect([1.4 1 1.4]);
set(gca,'linewidth',1.5);
hold off
% 
% set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
% frameName = [strcat('JV_log')];
% print(gcf, '-dpng', strcat(figuresdir7,'\',frameName),'-r0');
% savefig(gcf, fullfile(figuresdir7,frameName),'compact');

% close(gcf)

%%


rev_idx = find(JVD(:,1) == -.25);
fwd_idx = find(JVD(:,1) == .11);
zero_idx = find(JVD(:,1) == -.01);


dj = filter(smooth_diff(10),1,JVD(:,2:end));
dv = filter(smooth_diff(10),1,JVD(:,1));

djdv = dj./dv;

shunt_g = zeros(1,size(djdv,2));

%%
figure();
hold on;
for k = 4:t_max

    h2 = plot(JVD(rev_idx:zero_idx,1),djdv(rev_idx:zero_idx,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');

%     shunt_g(1,k) = (djdv(rev_idx,k));
    shunt_g(1,k) = mean((djdv(rev_idx:zero_idx,k)));

end
axis square
box on

x = xlabel ('V (volt)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
xlim([-.26 0])
% set(gca,'xtick',[240 270 300 330])

y = ylabel ('dJ/dV (^{mS}/_{cm^2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',46);
% ylim([1e-3 1e2])
set(gca,'yscale','log');
% set(gca,'ytick',[1e-3 1e-1 1e1])

pbaspect([1 1 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('shunt_g')];
print(gcf, '-dpng', strcat(figuresdir7,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir7,frameName),'compact');

close(gcf)



%%

fwd_idx2 = find(JVD(:,1) == .15);
%%
dvdj = dv./dj;

jgv = zeros(jvm,t_max);

for k = 1:t_max
    jgv(:,k) = JVD(:,k+1)-(shunt_g(1,k).*JVD(:,1));
end

%%
figure();
hold on;
for k = 12:15
    h2 = plot((1./jgv(zero_idx+1:fwd_idx2,k)),dvdj(zero_idx+1:fwd_idx2,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
end
axis square
box on

x = xlabel ('J^{-1} (mA^{-1}cm^{2})');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
% xlim([0 .04])
% set(gca,'xtick',[240 270 300 330])

y = ylabel ('dV/dJ ({\Omega}{cm^2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
% ylim([1e-3 1e2])
% set(gca,'yscale','log');
% set(gca,'ytick',[1e-3 1e-1 1e1])

pbaspect([1 1 1]);
set(gca,'linewidth',1.5);
hold off

%%

rfit = zeros(2,t_max);

% rfit(1,12) = .5418;
% rfit(2,12) = .1184;
% 
% rfit(1,13) = .4156;
% rfit(2,13) = .1071;
% 
% rfit(1,14) = .4558;
% rfit(2,14) = .09583;
% 
% rfit(1,15) = .3284;
% rfit(2,15) = .008617;

% rfit(1,12) = .1184;
% rfit(2,12) = .04552;
% 
% rfit(1,13) = .1186;
% rfit(2,13) = .04126;
% 
% rfit(1,14) = .1181;
% rfit(2,14) = .03733;
% 
% rfit(1,15) = .1331;
% rfit(2,15) = .03412;

rfit(1,12) = 4.966;
rfit(2,12) = 1.175;

rfit(1,13) = 6.263;
rfit(2,13) = 1.364;

rfit(1,14) = 6.391;
rfit(2,14) = 1.496;

rfit(1,15) = 5.238;
rfit(2,15) = 1.493;


% rfit(1,1) = 4741;
% rfit(2,1) = 52.51;
% 
% rfit(1,2) = 4696;
% rfit(2,2) = 53.57;
% 
% rfit(1,3) = 4608;
% rfit(2,3) = 65.12;

rfit_idx1 = zeros(1,t_max);
rfit_idx2 = zeros(1,t_max);

for k = 12:15
    tmp1 = 1./jgv(:,k) - rfit(1,k);
    [min1,rfit_idx1(1,k)] = min(abs(tmp1));

    tmp2 = 1./jgv(:,k) - rfit(2,k);
    [min2,rfit_idx2(1,k)] = min(abs(tmp2));
end

rfit_line = cell(1,t_max);
rfits = zeros(t_max,2);
rfit_err = cell(1,t_max);
r_ext = zeros(1,t_max);

for k = 12:15
    [rfits(k,:), rfit_err{1,k}] = polyfit((1./jgv(rfit_idx1(1,k):rfit_idx2(1,k),k)),dvdj(rfit_idx1(1,k):rfit_idx2(1,k),k),1);
    rfit_line{1,k}(:,1) = linspace((1./jgv(rfit_idx1(1,k),k)),(1./jgv(rfit_idx2(1,k),k)),100);
    rfit_line{1,k}(:,2) = rfits(k,2) + rfit_line{1,k}(:,1).*rfits(k,1);
    r_ext(1,k) = rfits(k,2);
end
%%
figure();
hold on

for k = 12:15

    h1 = plot((1./jgv(zero_idx+1:fwd_idx2,k)),dvdj(zero_idx+1:fwd_idx2,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
    h2 = plot(rfit_line{1,k}(:,1),rfit_line{1,k}(:,2),'Color','r', 'LineWidth',3);

end

axis square
box on

x = xlabel ('(J-GV)^{-1} (mA^{-1}cm^{2})');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
xlim([0 8])
% set(gca,'xtick',[0 2e3 4e3 6e3])

y = ylabel ('dV/dJ ({\Omega}{cm^2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
ylim([.02 .06])
% set(gca,'yscale','log');
% set(gca,'ytick',[0 0.025 .05])

pbaspect([1 1 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('r_extract')];
print(gcf, '-dpng', strcat(figuresdir7,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir7,frameName),'compact');

close(gcf)

%%
%make plot of r_0 vs. temp
figure()
hold on
for k = 12:15

    h1 = plot(T(k,1),1e3*r_ext(1,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
%     h2 = plot(rfit_line{1,k}(:,1),rfit_line{1,k}(:,2),'Color','r', 'LineWidth',3);

end

axis square
box on

x = xlabel ('T (K)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
xlim([250 300])
set(gca,'xtick',[250 275 300])

y = ylabel ('R_{0} (\Omega cm^{-2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
% ylim([.02 .06])
% set(gca,'yscale','log');
% set(gca,'ytick',[0 0.025 .05])

pbaspect([1 1.6 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('r0-T')];
print(gcf, '-dpng', strcat(figuresdir7,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir7,frameName),'compact');

close(gcf)





%%

vrj = zeros(length(JVD),t_max);

for k = 12:15
    vrj(:,k) = JVD(:,1) - r_ext(1,k).*JVD(:,k+1);  
end


figure();
hold on;
for k = 12:15
%     h1 = plot(vrj(zero_idx:end,k),jgv(zero_idx:end,k),'Color',colSet_red(k,:),'Marker',markz{1,k},'LineWidth',1,'MarkerFaceColor',colSet_red(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
    
%     h1 = plot(vrj(zero_idx:fwd_idx2,k),jgv(zero_idx:fwd_idx2,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');
    
    h1 = plot(vrj(zero_idx:44,k),jgv(zero_idx:44,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',1,'MarkerFaceColor',colSet_green(k,:), 'MarkerSize',20, 'MarkerEdgeColor','k');



end

axis square
box on

x = xlabel ('V-RJ (V)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
% xlim([-.01 .05])
% set(gca,'xtick',[240 270 300 330])

y = ylabel ('J-GV (mAcm^{-2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
% ylim([1e-3 1e2])
set(gca,'yscale','log');
% set(gca,'ytick',[1e-3 1e-1 1e1])

pbaspect([1 1 1]);
set(gca,'linewidth',1.5);
hold off


%%

jfit = zeros(2,t_max);

jfit(1,12) = .1544;
jfit(2,12) = .1072;

jfit(1,13) = .1572;
jfit(2,13) = .08778;

jfit(1,14) = .1405;
jfit(2,14) = .08266;

jfit(1,15) = .1701;
jfit(2,15) = .08642;

% jfit(1,1) = .1564;
% jfit(2,1) = .04832;
% 
% jfit(1,2) = .1618;
% jfit(2,2) = .04859;
% 
% jfit(1,3) = .1619;
% jfit(2,3) = .04856;


jfit_idx1 = zeros(1,t_max);
jfit_idx2 = zeros(1,t_max);


for k = 12:15
    tmp1 = vrj(:,k) - jfit(1,k);
    [min1,jfit_idx1(1,k)] = min(abs(tmp1));

    tmp2 = vrj(:,k) - jfit(2,k);
    [min2,jfit_idx2(1,k)] = min(abs(tmp2));
end

jfit_line = cell(1,2);
jfits = zeros(t_max,2);
jfit_err = cell(1,t_max);
j0 = zeros(1,t_max);

for k = 12:15
    
    [jfits(k,:), jfit_err{1,k}] = polyfit(vrj(jfit_idx2(1,k):jfit_idx1(1,k),k), log(jgv(jfit_idx2(1,k):jfit_idx1(1,k),k)), 1);
%     jfit_line{1,k}(:,1) = linspace( vrj(jfit_idx2(1,k),k), vrj(jfit_idx1(1,k),k), 100);
    jfit_line{1,k}(:,1) = linspace(0, vrj(jfit_idx1(1,k),k), 100);
    jfit_line{1,k}(:,2) = exp(jfits(k,2) + jfit_line{1,k}(:,1).*jfits(k,1));
    j0(1,k) = exp(jfits(k,2));
    
end    


figure;
hold on;
for k = 12:15
    
    h1 = plot(vrj(zero_idx:fwd_idx2,k),jgv(zero_idx:fwd_idx2,k),'Color',colSet_green(k,:),'Marker','o','LineWidth',2,'MarkerFaceColor',colSet_green(k,:),'MarkerSize',20, 'MarkerEdgeColor','k');
    h2 = plot(jfit_line{1,k}(:,1),jfit_line{1,k}(:,2),'Color','r', 'LineWidth',3);
    
end

axis square
box on

x = xlabel ('V-RJ (V)');
set(x,'FontName', 'Calibri');
set(gca,'FontSize',46);
% xlim([0 .25])
% set(gca,'xtick',[240 270 300 330])

y = ylabel ('J-GV (mAcm^{-2})');
set(y,'FontName','Calibri');
set(gca,'FontSize',40);
% ylim([1e-6 1e-3])
set(gca,'yscale','log');
% set(gca,'ytick',[1e-3 1e-1 1e1])

pbaspect([1 1 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('j0_extract')];
print(gcf, '-dpng', strcat(figuresdir7,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir7,frameName),'compact');

%%

d_rec = tau_calc*(j_0*DLN_tot{13,freq_cutoff}(1,1)/(q*n_i^2))^2;

d_calc = tau_rec*(j_0*DLN_tot{13,freq_cutoff}(1,1)/(q*n_i^2))^2;

d_est = (mu(13,1)*kB*300)/q;

