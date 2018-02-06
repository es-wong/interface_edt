%This script analyzes the output of the freeze out script and extracts
%frequcnies of both conductance peaks at low temperature when both are well
%defined. It then finds the conductance values at those frequencies, and
%calculates the resulting resistances. Using those resistances, we can then
%determine the conductivity of the PbS region. 

figuresdir6 = 'C:\Users\Eric\Desktop\Lab_Work\admittance\ligand_and_interface\analysis_work\EDT_schottky\09_2017_m1\analysis\figures\conductance';

for k = 1:t_max
    for ii = 1:b_max
        AS{k,ii}(:,7) = AS{k,ii}(:,5)./AS{k,ii}(:,1); %calculate G/omega 
        AS{k,ii}(:,8) = smooth(AS{k,ii}(:,7),9,'sgolay'); %smooth G/omega
    end
end


%plot the data
%%
for ii = 1 
    figure();
    hold on;
    for k = 8:15
%      for k = 
        plot(AS{k,ii}(:,1),(1e9*AS{k,ii}(:,8)),'LineWidth',3,'Color',colSet_green(k,:));   
     end
    axis square;
    box on;

    x = xlabel ('\omega (rad/s)', 'FontSize', 46); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
    xlim([1e5 1e6]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
    
    y = ylabel('^{G}/_{\omega} (nSs^{-1})', 'FontSize', 46);
    set(y,'FontName','Calibri');
    ylim([0 20]);    
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('G-omega__sub1',num2str(0),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir6,frameName),'compact');
%     close(gcf)
end

%%
%Fit the peaks assuming 2 Gaussian peaks
% 
freq_G1 = input('What is the lower bound frequency for which to start the Gaussian fits to the derivative spectra? ');
freq_G2 = input('What is the upper bound frequency for which to start the Gaussian fits to the derivative spectra? ');

freq_test = zeros(size(AS{1,1},1),3);
freq_test = AS{1,1}(:,1);
freq_test(:,2) = abs(freq_test(:,1) - freq_G1);
freq_test(:,3) = abs(freq_test(:,1) - freq_G2);


fitidx = zeros(1,2); %array to store the bounds of the frequency range of interest
[min1,fitidx(1,1)] = min(freq_test(:,2));
[min2,fitidx(1,2)] = min(freq_test(:,3));

G_gauss2fit1 = cell(t_max,b_max);
G_gauss2fit2 = cell(t_max,b_max);


G_Y_max = zeros(t_max,b_max);
G_Y_norm = cell(t_max,b_max);

G_w0_idx = zeros(t_max,v_max);
G_wd_idx = zeros(t_max,v_max);

Gpeaks_w0 = zeros(t_max,b_max);
Gpeaks_wd = zeros(t_max,b_max);

Gpeaks_err2 = cell(t_max,b_max);
Gpeaks_w0_err = zeros(t_max,b_max);
Gpeaks_wd_err = zeros(t_max,b_max);

g1 = zeros(t_max,b_max);
g12 = zeros(t_max,b_max);

%%
%perform the fit assuming 2 Gaussian peaks

pts = input('How many points on either side of the max should be fit as a Gaussian for the second fit? ');


for k = 1+offset_ref(1,1):offset_ref(1,2)
    for ii = 1
    
    f1 = fit(log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),AS{k,ii}(fitidx(1,1):fitidx(1,2),8),'gauss2'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),ASderiv_norm{ii,k}(fitidx(1,1):fitidx(1,2),1)); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);
    
    G_gauss2fit1{k,ii} = p1;
    
 
    G_Y_max(k,ii) = max(AS{k,ii}(fitidx(1,1):fitidx(1,2),8)); %extract the max of the y-range
    G_Y_norm{k,ii}(:,1) = AS{k,ii}(fitidx(1,1):fitidx(1,2),8)./G_Y_max(k,ii); %normalize the y-data against the max
    
    G_gauss2fit_scale = G_gauss2fit1;
    G_gauss2fit_scale{k,1}(1,1) = G_gauss2fit_scale{k,1}(1,1)/G_Y_max(k,ii);
    G_gauss2fit_scale{k,1}(1,4) = G_gauss2fit_scale{k,1}(1,4)/G_Y_max(k,ii);
    
    
    options = fitoptions('gauss2','StartPoint', G_gauss2fit_scale{k,ii}); %using the first fit, input fit parameters for the second fit
    [f1a,test] = fit(log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),G_Y_norm{k,ii},'gauss2',options); %perform the second fit
    G_gauss2fit2{k,ii} = coeffvalues(f1a); %extract the coefficients of the fit
    e1a = confint(f1a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    Gpeaks_err2{k,ii} = e1a;
   
    
 
    %determine which set of coefficients correpond to the AS peak (lower
    %frequency peak)
    [min_G2,G_w0_idx(k,ii)] = min([1e7 G_gauss2fit2{k,ii}(1,2) 1e7 1e7 G_gauss2fit2{k,ii}(1,5) 1e7]); 
    
    Gpeaks_w0(k,ii) = exp(G_gauss2fit1{k,ii}(1,G_w0_idx(k,ii)));%extract the value of the peak frequency
    tmp_0 = (AS{k,ii}(fitidx(1,1):fitidx(1,2),1)-Gpeaks_w0(k,ii)); %subtract the peak value from the frequency range in the experiment
    [val idx_0] = min(abs(tmp_0)); %find the minimum, picking out the peak in the experimental data
    w0_idx = round((fitidx(1,1) + idx_0)-1,0); %find the index of the peak in the global data range
    g1(k,ii) = AS{k,ii}(w0_idx,5); %extract the conductance of the peak
    
    
    %set the other set of coefficients to the conductance peak
    if G_w0_idx(k,ii) == 2
        G_wd_idx(k,ii) = 5;
    elseif G_w0_idx(k,ii) == 5
        G_wd_idx(k,ii) = 2;
    end
   
    %perform the same conductance extraction for the freeze-out peak
    Gpeaks_wd(k,ii) = exp(G_gauss2fit1{k,ii}(1,G_wd_idx(k,ii)));
    tmp_d = (AS{k,ii}(fitidx(1,1):fitidx(1,2),1)-Gpeaks_wd(k,ii));  
    [val idx_d] = min(abs(tmp_d)); 
    wd_idx = round((fitidx(1,1) + idx_d)-1,0); 
    g12(k,ii) = AS{k,ii}(wd_idx,5);  
   
    
    Gpeaks_w0_err(k,ii) = (exp(Gpeaks_err2{k,ii}(2,G_w0_idx(k,ii))) - exp(Gpeaks_err2{k,ii}(1,G_w0_idx(k,ii))))/2;
    Gpeaks_wd_err(k,ii) = (exp(Gpeaks_err2{k,ii}(2,G_wd_idx(k,ii))) - exp(Gpeaks_err2{k,ii}(1,G_wd_idx(k,ii))))/2;
    

    
    figure()
    hold on;
    plot(f1a,log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),(AS{k,ii}(fitidx(1,1):fitidx(1,2),8)./G_Y_max(k,ii))); 
    axis square
    t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
    
    end
end


%%
%perform the fit assuming the existence of only 1 gaussian peak

G_gauss1fit1 = cell(t_max,b_max);
G_gauss1fit2 = cell(t_max,b_max);

G_Y_max1 = zeros(t_max,b_max);
G_Y_norm1 = cell(t_max,b_max);

G_w0_idx1 = zeros(t_max,v_max);

Gpeaks_w01 = zeros(t_max,b_max);

Gpeaks_err1 = cell(t_max,b_max);
Gpeaks_w01_err = zeros(t_max,b_max);

cond_peak = zeros(t_max,b_max);

peak_guess = zeros(t_max,b_max);

%%
pts = input('How many points on either side of the max should be fit as a Gaussian for the second fit? ');

for k = 1+offset_ref(1,1):offset_ref(1,2)
    for ii = 1
    
    f1 = fit(log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),AS{k,ii}(fitidx(1,1):fitidx(1,2),8),'gauss1'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),ASderiv_norm{ii,k}(fitidx(1,1):fitidx(1,2),1)); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);
    
    G_gauss1fit1{k,ii} = p1;
    peak_tmp = exp(G_gauss1fit1{k,ii}(1,2));
    tmp = (AS{k,ii}(fitidx(1,1):fitidx(1,2),1)-peak_tmp); 
    [idx idx] = min(abs(tmp));
    guess_idx = round((fitidx(1,1) + idx)-1,0); 
    peak_guess(k,ii) = AS{k,ii}(guess_idx,1); 
    
    
    gX_sub = AS{k,ii}((guess_idx-pts):(guess_idx+pts),1); %extract a sub-range of x-data
    gY_sub = AS{k,ii}((guess_idx-pts):(guess_idx+pts),8); %extract a sub-range of y-data
    gY_max = max(gY_sub); %extract the max of the y-range
    gY_sub_prime = gY_sub./gY_max; %normalize the y-data against the max
    
    options = fitoptions('gauss1','StartPoint', G_gauss1fit_scale{k,ii}); %using the first fit, input fit parameters for the second fit
    [f1a,test] = fit(log(gX_sub),gY_sub_prime,'gauss1',options); %perform the second fit
    G_gauss1fit2{k,ii} = coeffvalues(f1a); %extract the coefficients of the fit
    e1a = confint(f1a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    Gpeaks_err1{k,ii} = e1a;
     
    Gpeaks_w01(k,ii) = exp(G_gauss1fit2{k,ii}(1,2));%extract the value of the peak frequency
    tmp_0 = (AS{k,ii}(fitidx(1,1):fitidx(1,2),1)-Gpeaks_w01(k,ii)); %subtract the peak value from the frequency range in the experiment
    [val idx_0] = min(abs(tmp_0)); %find the minimum, picking out the peak in the experimental data
    w0_idx1 = round((fitidx(1,1) + idx_0)-1,0); %find the index of the peak in the global data range
    cond_peak(k,ii) = AS{k,ii}(w0_idx1,5); %extract the conductance of the peak
    
    Gpeaks_w01_err(k,ii) = exp(Gpeaks_err1{k,ii}(2,2) - Gpeaks_err1{k,ii}(1,2))/2;
    
    figure()
    hold on
    plot(f1a,log(AS{k,ii}(fitidx(1,1):fitidx(1,2),1)),(AS{k,ii}(fitidx(1,1):fitidx(1,2),8)./G_Y_max1(k,ii))); 
    axis square
    t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
    
    end
end

%%
for ii = 1 
    figure();
    hold on;
%     for k = 1 + offset_ref(1,1):offset_ref(1,2)
      for k = 12:15
%         h1 = plot(T(k,1),Gpeaks_w01(k,ii),'Marker','o','MarkerFaceColor',colSet_blue(k,:),'Color',colSet_blue(k,:),'MarkerEdgeColor','k'); 
%             set(h1,{'MarkerSize'},{20});
        h1 = plot(T(k,1),Gpeaks_w0(k,ii),'Marker','o','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k'); 
            set(h1,{'MarkerSize'},{20});
%         h2 = plot(T(k,1),Gpeaks_wd(k,ii),'Marker','d','MarkerFaceColor',colSet_blue(k,:),'Color',colSet_blue(k,:),'MarkerEdgeColor','k'); 
%             set(h2,{'MarkerSize'},{20});
       end
    axis square;
    box on;

    x = xlabel ('T (K)', 'FontSize', 46); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'xscale','log');
    xlim([250 300]);
    set(gca,'xtick',[250 275 300]);
     
    y = ylabel('\omega{_0} (rad/s)', 'FontSize', 46);
    set(y,'FontName','Calibri');
    set(gca,'yscale','log');
    ylim([2e5 5e5]);   
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off;
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('conudcance_omega_0_limited_T_',num2str(0),'mV')];
    print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir6,frameName),'compact');
    close(gcf)
end

%%
%find the indices of the conductance peaks to extract the non-frequency
%normalized conductance (Gp). this will be a good comparison
%for our calculation. 

sigIdx = zeros(1,t_max);
sigIdx_test = zeros(length(freq_test),t_max);
Gp_w0 = zeros(1,t_max); %array to store Gp

for k = 12:15
    sigIdx_test(:,k) = abs(freq_test(:,1) - Gpeaks_w0(k,1));
    [min1,sigIdx(1,k)] = min(sigIdx_test(:,k));
    Gp_w0(1,k) = AS{k,1}(sigIdx(1,k),5);
end

for ii = 1 
    figure();
    hold on;
    for k = 12:15
        
        h1 = plot(T(k,1),Gp_w0(1,k),'Marker','o','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k'); 
            set(h1,{'MarkerSize'},{20});
    end
    axis square;
    box on;

    x = xlabel ('T (K)', 'FontSize', 46); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'xscale','log');
    xlim([250 300]);
    set(gca,'xtick',[250 275 300]);
     
    y = ylabel('G_{p} (S)', 'FontSize', 46);
    set(y,'FontName','Calibri');
%     set(gca,'yscale','log');
    ylim([2e-3 8e-3]);   
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off;
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('Gp_w0_',num2str(0),'mV')];
    print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir6,frameName),'compact');
    close(gcf)
end




%%

w0_energy = zeros(t_max,b_max);
% wd_energy = zeros(t_max,b_max);

for k = 12:15
    w0_energy(k,ii) = kB*T(k,1)*log((2*nu_00(1,1)*T(k,1)^2)/(Gpeaks_w0(k,ii)));
%     wd_energy(k,ii) = kB*T(k,1)*log((2*nu_0(1,1))/(Gpeaks_wd(k,ii)));
end


for ii = 1 
    figure();
    hold on;
    for k = 12:15
%       for k = 1+ offset_ref(1,1):5
        h1 = plot(T(k,1),w0_energy(k,ii),'Marker','.','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:)); 
            set(h1,{'MarkerSize'},{50});
%         h2 = plot(T(k,1),wd_energy(k,ii),'Marker','d','MarkerFaceColor',colSet_blue(k,:),'Color',colSet_blue(k,:)); 
%             set(h2,{'MarkerSize'},{10});
    end
    axis square;
    box on;
%     t = title(strcat('C vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('T (K)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
    y = ylabel('E_{\omega}', 'FontSize', 46);
    set(y,'FontName','Calibri');
%     set(gca,'yscale','log');
%     xlim([1.8e2 1e6]);
%     set(gca,'xscale','log');
%     set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([.5e-8 1.5e-8]);    
    hold off;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('omega_0_energy_',num2str(0),'mV')];
    print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir6,frameName),'compact');
    close(gcf)
end

%%

%This section of the code calculates the conductivity based on the dielectric relaxation frequency and fits the
%conductivity to an Arrhenius model

sigma_0 = zeros(t_max,b_max);
sigma_0_err = zeros(t_max,b_max);
sigma_0_ln = zeros(t_max,b_max);

E_a_sigma = zeros(2,b_max);
sigma_00 = zeros(2,b_max);


for k = 12:15
    for ii = 1
        sigma_0(k,ii) = epsilon*eps_0*Gpeaks_w0(k,ii);
%         sigma_0_err(k,ii) = sigma_0(k,ii)*(Gpeaks_w0_err(k,ii)/Gpeaks_w0(k,ii));
        sigma_0_err(k,ii) = .1*sigma_0(k,ii);
        sigma_0_ln(k,ii) = log(sigma_0(k,ii)/(T(k,1)^2));
    end
end

%create the arrhenius plot
% X_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store the x-values for the Arrenhius plot
% Y_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store the y-values for the Arrhenius plot
% Yhat_sigma = zeros(100,b_max); %create an array to store the predicted y-values for the plot
% Yhat_resid_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store the predicted y-values for the residual calculation
% beta_sigma = zeros(2,b_max); %crete an array to store the fit parameters
% errbeta_sigma = cell(1,b_max); %create an array to store the covariance matrix associated with each fit
% W_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store the inverse variance for weights
% Wprime_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),1); %create a temporary array to store the weights for the particular fit being performed
% fiterr_sigma = zeros(2,b_max); %create an array to calculate the error for each paramter
% resid_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store the residuals of each fit
% unit_test1_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store sum of the square of the residulas divided by the number of points
% unit_test2_sigma = zeros(offset_ref(1,2)-offset_ref(1,1),b_max); %create an array to store the sum of the variances for each data point

X_sigma = zeros(4,b_max); 
Y_sigma = zeros(4,b_max); 
Yhat_sigma = zeros(100,b_max); 
Yhat_resid_sigma = zeros(4,b_max); 
beta_sigma = zeros(2,b_max); 
errbeta_sigma = cell(1,b_max); 
W_sigma = zeros(4,b_max); 
Wprime_sigma = zeros(4,1); 
fiterr_sigma = zeros(2,b_max); 
resid_sigma = zeros(4,b_max); 
unit_test1_sigma = zeros(4,b_max); 
unit_test2_sigma = zeros(4,b_max); 



X_vals_sigma = linspace(1e3*(1/250), 1e3*(1/300), 100); %create an array of x-values for the prediction line

for ii = 1 
    
%     X_sigma(:,ii) = 1./(T(1+offset_ref(1,1):offset_ref(1,2),1)); %extract the inverse temperatures
    X_sigma(:,ii) = 1e3*(1./(T(12:15,1)));
    
%     Y_sigma(:,ii) = sigma_0_ln(1+offset_ref(1,1):offset_ref(1,2),ii); %extract ln(sigma_0)
%     Y_sigma(:,ii) = sigma_0_ln(12:15,ii);
%     Y_sigma(:,ii) = Gp_w0(1,12:15);
    Y_sigma(:,ii) = log(Gpeaks_w0(12:15,ii)./(T(12:15,1).^2));
    
    Xtilde_sigma = [ones(length(X_sigma(:,ii)),1),X_sigma(:,ii)]; %create an array of ones and the x-values  for the fit specified by bias 'ii' 
    W_sigma(:,ii) = 1; %give the experimental data points equal weight, assuming that the error in the trend is greater than the error in each individual measurement
%     W(:,ii) = (gauss_err(1+offset_ref(1,1):offset_ref(1,2),ii).^2); %create the weighting array. Note that we are interested in 1/sigma^2 for each data point
    Wprime_sigma(:,1) = W_sigma(:,ii); %extract the variances for the fit being run
    
    beta_sigma(:,ii) =  (Xtilde_sigma.' * diag(Wprime_sigma)^(-1) * Xtilde_sigma)^(-1) * (Xtilde_sigma.' * diag(Wprime_sigma)^(-1) * Y_sigma(:,ii)); %calculate the regression coefficients
    Yhat_sigma(:,ii) = X_vals_sigma(1,:).*beta_sigma(2,ii) + beta_sigma(1,ii); %generate the array of predicted values for the plot

    Yhat_resid_sigma(:,ii) = X_sigma(:,ii).*beta_sigma(2,ii) + beta_sigma(1,ii); %generate the array of predicted values for the residuals calculation
    
    resid_sigma(:,ii) = Y_sigma(:,ii) - Yhat_resid_sigma(:,ii); %calculate the residuals
    ei_sq_sigma = mean(resid_sigma(:,ii).^2);
    
    errbeta_sigma{1,ii} = ei_sq_sigma.*(Xtilde_sigma.' *diag(Wprime_sigma)^(-1)* Xtilde_sigma)^(-1); %calculate the variance in the regression coefficients
%     errbeta{1,ii} = (Xtilde.' *diag(Wprime)^(-1)* Xtilde)^(-1); %calculate the variance in the regression coefficients
    fiterr_sigma(1,ii) = sqrt(errbeta_sigma{1,ii}(1,1)/(offset_ref(1,2)-offset_ref(1,1))); %normalize the error in the constant properly
    fiterr_sigma(2,ii) = sqrt(errbeta_sigma{1,ii}(2,2)/(offset_ref(1,2)-offset_ref(1,1))); %normalize the error in the slope properly
     
    E_a_sigma(1,ii) = -1e3*kB*beta_sigma(2,ii); %calculate the trap energy
    E_a_sigma(2,ii) = -1e3*kB*fiterr_sigma(2,ii);%calculate the 1-sigma error in the energy
    %%
    sigma_00(1,ii) = exp(beta_sigma(1,ii))/2; %calculate the attempt to escape frequency in rad/s
    sigma_00(2,ii) = sigma_00(1,ii)*fiterr_sigma(1,ii); %calculate the 1-sigma error in the attempt to escape frequency in rad/s
    %%
%     for k = 1:t_max;
%         nu_0{k,ii}(1,1) = T(k,1)^2*nu_00(1,ii); %calculate the
%         temperature-dependent attempt to escape frequency in rad/(sK^2)
%         nu_0{k,ii}(1,2) = T(k,1)^2*nu_00(2,ii); %calculat the 1-sigma
%         error in the attempt to escape frequency in rad/(sK^2)
%     end;
  
    unit_test1_sigma(:,ii) = (sum(resid_sigma(:,ii).^2)/length(Y_sigma)); %calculate the sum of the square of the residuals
    unit_test2_sigma(:,ii) = (sum(W_sigma(:,ii))/length(Y_sigma)); %calculate the sum of the variances of the data points
    
    %note that the unit tests should be of the same order of magnitude if W
    %is specified correctly
    
end
%%
%create the Arrhenius plot
for ii = 1 
    figure
    
    g1 = errorbar(X_sigma(:,ii),Y_sigma(:,ii),sigma_0_err(12:15,ii),'o','MarkerFaceColor',colSet_green(k,:));
%     g1 = plot(X_sigma(:,ii),Y_sigma(:,ii),'o','MarkerFaceColor',colSet_green(k,:)); %create the scatter plot of ln(1/f) vs. 1/T
%     g1 = errorbar(X_sigma(:,ii),Y_sigma(:,ii),sigma_0_err(1+offset_ref(1,1):offset_ref(1,2),ii),'o','MarkerFaceColor',colSet_blue(k,:)); %create the scatter plot of ln(1/f) vs. 1/T
    set(g1,{'markers'},{20},{'Linewidth'},{1},{'MarkerEdgeColor'},{'k'});
    
    hold on;
    plot(X_vals_sigma(1,:),Yhat_sigma(:,ii),'--r','LineWidth', 3); %plot the best fit line over it

    axis square;
    box on;
    x = xlabel ('1000/T (K^{-1})', 'FontSize', 46);
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
    xlim([3.3 4.0]);
    set(gca,'xtick',[3.3 3.65 4]);
    
    y = ylabel ('ln(\omega_{0}T^{-2})) (rad s^{-1}K^{-2})', 'FontSize',46);
    set(y,'FontName', 'Calibri')
%     ylim([-21.5 -20.5]);
%     set(gca,'ytick',[-21.5 -21 -20.5 -20]);
    
    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off

%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('Gw0_arrhenius_',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir6,frameName),'compact');
%     close(gcf)
end
%%
for ii = 1 
    figure() 
    
       h1 = errorbar (T(12:15,1),sigma_0(12:15,ii),sigma_0_err(12:15,ii),'Marker','o','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k');
%     h1 = plot(T(12:15,1),sigma_0(12:15,ii),'Marker','o','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k'); 
%     h1 = errorbar (T(1+offset_ref(1,1):offset_ref(1,2),1),sigma_0(1+offset_ref(1,1):offset_ref(1,2),ii),sigma_0_err(1+offset_ref(1,1):offset_ref(1,2),ii),'Marker','o','MarkerFaceColor',colSet_blue(k,:),'Color',colSet_blue(k,:),'MarkerEdgeColor','k'); 
            set(h1,{'MarkerSize'},{20});
            
    axis square;
    box on;
    x = xlabel ('T (K)', 'FontSize', 46);
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
    xlim([250 300]);
    set(gca,'xtick',[250 275 300]);
    
    y = ylabel ('\sigma_R (S/m)', 'FontSize',46);
    set(y,'FontName', 'Calibri')
    set(gca,'yscale','log')
    ylim([4e-5 1.2e-4]);

    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('sigma_',num2str(B.(BNames{ii,1})),'mV')];
    print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir6,frameName),'compact');
    close(gcf)
end
%%

sigma_project = (300^2)*exp((1/300)*beta_sigma(2,1) + beta_sigma(1,1))
% sigma_project_err = sigma_project.* sqrt((E_a_sigma(2,1)/E_a_sigma(1,1))^2 + (sigma_00(2,1)/sigma_00(1,1))^2)

%%
%calculate the composite attempt frequency by considering nu_00 and
%sigma_00 together

nu00tot = zeros(2,3);
nu0tot = zeros(t_max,2);

for ii = 1:b_max
   
    nu00tot(1,ii) = (exp(abs(beta(1,ii)+beta_sigma(1,1)))/2);
    nu00tot(2,ii) = nu00tot(1,ii) *sqrt(fiterr_sigma(1,1)^2 + fiterr(1,ii)^2);
    
    for k = 1:t_max 
        nu0tot(k,1) = nu00tot(1,1)*T(k,1)^2;
        nu0tot(k,2) = nu00tot(2,1)/nu00tot(1,1) * nu0tot(k,1);
    end
    
end
%%
%need to come up with an energy scale, so let's play around with this a
%bit 

sigma0 = zeros(t_max,2);

for k = 1:t_max 
    sigma0(k,1) = sigma_00(1,1)*T(k,1)^2;
    sigma0(k,2) = (sigma_00(2,1)/sigma_00(1,1)) * sigma0(k,1);
end

test1 = kB*T(12:15,1).*log(2*nu0tot(12:15,1)./(Gpeaks_w0(12:15,1)))

test2 = kB*T(12:15,1).*log(2*sigma0(12:15,1)./(Gpeaks_w0(12:15,1)))


