%repeat the DOS_final calculation, but add the de-tuning energy in, just to
%rescale the energy. I don't know if this is kosher, but it's worth a shot.

%Calculates the measurement energy E_{\omega}. This requires a conversion
%form eV to joules

Ew_detune = cell(1,b_max);
Ew_err_detune = cell(1,b_max);

Ew_j_detune = cell(1,b_max);
Ew_j_err_detune = cell(1,b_max);

Efn = zeros(t_max,1);

kBJ = 1.3806485e-23; %Boltzmann constant in joules
eV_j = 1.60218e-19; %conversion factor for eV to joules
%%
for ii = 1
    for k = 1:t_max 
        %Ew_detune{1,ii}(:,k) = .08 + kB*T(k,1)*log((2*nu_0{1,ii}(k,1))./AS{k,ii}(:,1));
        %Ew_err_detune{1,ii}(:,k) = nu_0{1,ii}(k,2)/nu_0{1,ii}(k,1).*Ew_detune{1,ii}(:,k);
        
        Ew_detune{1,ii}(:,k) =  kB*T(k,1)*log((2*nu0tot(k,1))./AS{k,ii}(:,1));
        Ew_err_detune{1,ii}(:,k) = (nu0tot(k,2)/nu0tot(k,1)).*Ew_detune{1,ii}(:,k);
        
        Efn(k,1) = Ef(k,1) - E_v;  


    end
end
%%
for ii =1:b_max    
    Ew_j_detune{1,ii} = Ew_detune{1,ii}.*eV_j;
    Ew_j_err_detune{1,ii} = Ew_err_detune{1,ii}.*eV_j; 
    
    Efn_j = Efn.*eV_j; 
end

%%
%calculate the trap density N_{t} and plot 
%against the measurement energy E_{\omega} calculated above. note
%the units

Nt_detune = cell(2,b_max); 
Nt_s_detune = cell(1,b_max); 
Nt_int_detune = cell(1,b_max); 


%calculate the trap density in m^{-3} (cell row 1) and cm^{-3} cell row 2
for ii = 1
    for k = 1:t_max
        
        Nt_detune{1,ii}(:,k) = ((V_bi - Va(1,ii))^2)./(A*W.*(q*V_bi-(Efn_j(k,1) - Ew_j_detune{1,ii}(:,k)))).*ASderiv_norm{ii,k}(:,1)./(kB*T(k,1)); %calculate N_t,the density of trap states
        Nt_detune{2,ii} = Nt_detune{1,ii}*1e-6; 
        Nt_s_detune{1,ii}(:,k) = smooth(Nt_detune{2,ii}(:,k),9,'sgolay'); 
        Nt_int_detune{1,ii}(:,k) = trapz(Ew_detune{1,ii}(:,k),Nt_s_detune{1,ii}(:,k));
        
    end
end

%%
%now calculate the DOS based on the new method using the DLCP data. first
%set up the data for the gaussian fits to the TAS derivative

%first re-plot the derivate spectra to visualize the frequency bounds

for ii = 1
        figure();%Plot the derivative data according to bias level to look for trends accross this parameter
        hold on;
        for k = 1+offset_ref(1,1):offset_ref(2,2)
%         for k = offset_ref(2,2)-1
            plot(AS{k,ii}(:,1),ASderiv_norm{ii,k}(:,1),'LineWidth',3, 'Color', colSet_blue(k,:)); %Plot the remaining temperatures and edit the plot accordingly
        end
        axis square
        box on;
%         t = title (strcat('$\omega\frac{dC}{d\omega}$ vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',20);
%         set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%         set(t, 'units', 'normalized
        x = xlabel ('Frequency, \omega (rad/s)', 'FontSize', 40); %continue formatting the axes
        set(x,'FontName','Latex');
        set(gca,'FontSize',40);
        y = ylabel('\omega^{dC}/_{d \omega} (Farads)','FontSize', 40);
        set(y,'FontName','Calibri');
        set(gca,'xscale','log');
        set(gca,'XMinorTick','on');
        xlim([4e2 3e6]);
        set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%         xlim([2092 4.827e4]);
%         ylim([0 2.7e-9]); 
        hold off;
end

%%
% freq_dos1 = input('What is the lower frequency bound  for which to start the Gaussian fits to the derivative spectra? ');
freq_dos2 = input('What is the upper frequency bound for which to start the Gaussian fits to the derivative spectra? ');

freq_test = zeros(size(AS{1,1},1),3);
freq_test = AS{1,1}(:,1);
freq_test(:,2) = abs(freq_test(:,1) - freq_dos1);
freq_test(:,3) = abs(freq_test(:,1) - freq_dos2);

%array to store the bounds of the frequency range of interest
fitidx2 = zeros(1,2); 
[min1,fitidx2(1,1)] = min(freq_test(:,2));
[min2,fitidx2(1,2)] = min(freq_test(:,3));
%%
%array to store fit data
dos_fits_detune = cell(1,t_max);
dos_peaks_detune = zeros(t_max,b_max);
dos_amp_detune = zeros(t_max,b_max);
fit_width_detune = zeros(t_max,2);
width_scale_detune = zeros(t_max,2);
freq_width_detune = zeros(t_max,2);
energy_width_detune = zeros(t_max,2);
energy_stdv_detune = zeros(t_max,1);

%array for ratio of height to width of Gaussian fit
xi_detune = zeros(t_max,b_max);

%array to store peak trap density in cm^{-3}ev^{-1}
N_peak_detune = zeros(t_max,b_max);

%array to store scaling factor between derivative curve and N_peak
dos_scaling_detune = zeros(t_max,b_max);

%array to store normalized derivative values
dos_vals_detune = cell(t_max,b_max);
%%
pts_dos = input('How many points on either side of the max should be fit as a Gaussian for the second fit? ');

%%
%perform the fits

% for k = offset_ref(1,2)-1
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for ii = 1
    
    f1 = fit(log(AS{k,ii}(fitidx2(1,1):fitidx2(1,2),1)),ASderiv_norm{ii,k}(fitidx2(1,1):fitidx2(1,2),1),'gauss1'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,log(AS{k,ii}(fitidx2(1,1):fitidx2(1,2),1)),ASderiv_norm{ii,k}(fitidx2(1,1):fitidx2(1,2),1)); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);

    val = exp(p1(2)); %extract the peak fit value
    tmp = (AS{k,ii}(fitidx2(1,1):fitidx2(1,2),1)-val); %subtract the peak value from the frequency range in the experiment
    [idx idx] = min(abs(tmp)); %find the minimum, picking out the peak in the experimental data
    g_idx = round((fitidx2(1,1) + idx)-1,0); %find the index of the peak in the global data range
    closest = AS{k,ii}(g_idx,1); %extract the frequency of the peak

    dos_range = log(AS{k,ii}(fitidx2(1,1):fitidx2(1,2),1));
    range_max = max(dos_range);
    
    %extract ranges and max of both x and y data
    X_sub = log(AS{k,ii}((g_idx-pts_dos):(g_idx+pts_dos),1));
%     X_max = max(X_sub);
    Y_sub = ASderiv_norm{ii,k}((g_idx-pts_dos):(g_idx+pts_dos),1); 
    Y_max = max(Y_sub);
    
    %normalize data
    Y_sub_prime = Y_sub./Y_max; 
    X_sub_prime = X_sub./range_max;
    
    options = fitoptions('gauss1','StartPoint', [ASderiv_norm{ii,k}(g_idx,1), log(closest), log(AS{k,ii}((g_idx+pts_dos),1))-log(AS{k,ii}((g_idx-pts_dos),1))]); %using the first fit, input fit parameters for the second fit
    [f1a,test] = fit(X_sub_prime,Y_sub_prime,'gauss1',options); 
%     [f1a,test] = fit(log(X_sub_prime),Y_sub_prime,'gauss1');
    dos_fits_detune{k,ii} = coeffvalues(f1a); %extract the coefficients of the fit
    e1a = confint(f1a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    dos_peaks_detune(k,ii) = exp(dos_fits_detune{k,ii}(2)); %obtain the centroid of the gaussian
    
    dos_tmp = (AS{k,ii}(fitidx2(1,1):fitidx2(1,2),1)-val);
    [dos_idx dos_idx] = min(abs(dos_tmp));
    dos_g_idx = round((fitidx2(1,1) + dos_idx)-1,0);
    dos_amp_detune(k,ii) = ASderiv_norm{ii,k}(dos_g_idx,1);
    
    %calculate the width of the gaussian fit in normalized units, frequency, and energy   
    fit_width_detune(k,1) =  dos_fits_detune{k,ii}(2) + (0.5*(dos_fits_detune{k,ii}(1,3)));
    fit_width_detune(k,2) =  dos_fits_detune{k,ii}(2) - (0.5*(dos_fits_detune{k,ii}(1,3)));
    
    width_scale_detune(k,1) = fit_width_detune(k,1)*range_max;
    width_scale_detune(k,2) = fit_width_detune(k,2)*range_max;
    
    freq_width_detune(k,1) = exp(width_scale_detune(k,1));
    freq_width_detune(k,2) = exp(width_scale_detune(k,2));
    
    energy_width_detune(k,1) = kB*T(k,1)*log((2*nu0tot(k,1))/freq_width_detune(k,1));
    energy_width_detune(k,2) = kB*T(k,1)*log(2*nu0tot(k,1)/freq_width_detune(k,2));
    
    energy_stdv_detune(k,1) = abs(energy_width_detune(k,1) - energy_width_detune(k,2));
        
    figure()
    hold on;
%     plot(f1a,log(AS{k,ii}(fitidx2(1,1):fitidx2(1,2),1)),(ASderiv_norm{ii,k}(fitidx2(1,1):fitidx2(1,2),1)./Y_max)); 
    plot(f1a,dos_range./range_max,(ASderiv_norm{ii,k}(fitidx2(1,1):fitidx2(1,2),1)./Y_max)); 
    axis square
    t = title(strcat('T =', num2str(T(k,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
    
    
    N_peak_detune(k,ii) = DLN_tot{k,5}(1,1)/(energy_stdv_detune(k,1)*sqrt(pi));
    dos_scaling_detune(k,ii) = N_peak_detune(k,ii)/dos_amp_detune(k,ii);
    
    dos_vals_detune{k,ii}(:,1) = ASderiv_norm{ii,k}(:,1).*dos_scaling_detune(k,ii);
    
    
    end
end

%%
%plot the output of the 2 different DOS calculations

%first, plot the traditional method

for ii = 1 
    figure()
    hold on;

    for k = 12:15 %index over temperatures. Note that we only go to the offeset point in specified in the AS script
        
        plot(abs(Nt_s_detune{1,ii}(:,k))/1e18,Ew_detune{1,ii}(:,k) , 'LineWidth', 3,'Color', colSet_green(k,:));

    end

    box on;

    y = ylabel ('E_{\omega} (eV)');
    set(y,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    ylim([.26 .38]);
    set(gca,'ytick',[.26 .3  .34 .38]);
    
    x = xlabel('N_{t} (\times 10^{18} cm^{-3}eV^{-1})');
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
%     set(gca,'xscale','log');
    xlim([0 10]);
    set(gca,'xtick',[0 2 4 6 8 10]);


    pbaspect([1 .5 1]);
    set(gca,'linewidth',1.5);
    hold off;
    hold off
    box on
    
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_old_100kHZ_',num2str(Va(1,ii)*1000),'mV_detune')];
    print(gcf, '-dpng', strcat(figuresdir4,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir4,frameName),'compact');
      
    close(gcf)
end

%%

%next plot the new method

for ii = 1
    figure()
    hold on;
    for k = 12:15
        
        plot(abs(dos_vals_detune{k,ii}(:,1))/1e18,Ew_detune{1,ii}(:,k), 'LineWidth', 3,'Color', colSet_green(k,:));
    
    end
    axis square;
    box on;

    y = ylabel ('E_{\omega} (eV)');
    set(y,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    ylim([.26 .38]);
    set(gca,'ytick',[.26 .3  .34 .38]);
    
    x = xlabel('N_{t} (\times 10^{18} cm^{-3}eV^{-1})');
    set(x,'FontName','Calibri');
    set(gca,'FontSize',46);
%     set(gca,'xscale','log');
    xlim([0 10]);
    set(gca,'xtick',[0 2 4 6 8 10]);


    pbaspect([1 .5 1]);
    set(gca,'linewidth',1.5);
    hold off;
    hold off
    box on
    
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_new_100Khz_',num2str(Va(1,ii)*1000),'mV_detune')];
    print(gcf, '-dpng', strcat(figuresdir4,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir4,frameName),'compact');
    
    close(gcf)
end

%%
%create composite DOS plots with error bars

%input global endpoints within Ew that will define a window on which to perform
%the analysis
ew1_detune = .3;
ew2_detune= .4;

%create an array for interpolation across this domain so that we have the
%same number of points in each window for each temperature
ew_q = ew1_detune:1e-3:ew2_detune;

ew_idx = zeros(t_max,2);

old_spline_detune = zeros(t_max,length(ew_q));
new_spline_detune = zeros(t_max,length(ew_q));

old_agg_detune = zeros(2,length(ew_q));
new_agg_detune = zeros(2,length(ew_q));

%find the indices of the endpoints for each curve in the range of analysis and use a cubicspline to (essentially) fit the function over this domain 
for k = 12:15
   
   ew_tmp1 = abs(Ew_detune{1,1}(:,k) - ew1_detune);
   ew_tmp2 = abs(Ew_detune{1,1}(:,k) - ew2_detune);
 
   [min1,ew_idx(k,1)] = min(ew_tmp1);
   [min2,ew_idx(k,2)] = min(ew_tmp2);
    
   old_spline_detune(k,:) = interp1(Ew_detune{1,1}(ew_idx(k,2):ew_idx(k,1),k),Nt_s_detune{1,ii}(ew_idx(k,2):ew_idx(k,1),k),ew_q,'spline');
   new_spline_detune(k,:) = interp1(Ew_detune{1,1}(ew_idx(k,2):ew_idx(k,1),k),dos_vals_detune{k,ii}(ew_idx(k,2):ew_idx(k,1),1),ew_q,'spline');
   
end

%calculate the composite curve
for int_idx = 1:length(ew_q)

  old_agg_detune(1,int_idx) = mean(abs(old_spline_detune(1+offset_ref(1,1):offset_ref(1,2),int_idx))); 
  old_agg_detune(2,int_idx) = std(abs(old_spline_detune(1+offset_ref(1,1):offset_ref(1,2),int_idx)))/length(1+offset_ref(1,1):offset_ref(1,2));

  new_agg_detune(1,int_idx) = mean(abs(new_spline_detune(1+offset_ref(1,1):offset_ref(1,2),int_idx))); 
  new_agg_detune(2,int_idx) = std(abs(new_spline_detune(1+offset_ref(1,1):offset_ref(1,2),int_idx)))/length(1+offset_ref(1,1):offset_ref(1,2));

end

%%

%re-plot the output of the 2 DOS calculations with the average to see just
%how shitty a job we did


for ii = 1 
    figure()
    hold on;

    for k = 12:15
        
%         plot(Ew{1,ii}(:,k),abs(Nt_s{1,ii}(:,k)), 'LineWidth', 3,'Color', colSet_blue(k,:));

        errorbar(ew_q,old_agg_detune(1,:),old_agg_detune(2,:),'Color',colSet_green(15,:),'Marker','o','LineWidth', 2,'MarkerFaceColor',colSet_green(15,:),'MarkerEdgeColor',colSet_green(1,:),'MarkerSize',20);
       
    end

    box on;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.3 .4]);
    set(gca,'xtick',[.3 .35 .4]);
    
    y = ylabel('N_{t} (cm^{-3}eV^{-1})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    ylim([5e16 5e18])
    set(gca,'ytick',[1e16 1e17 1e18 1e19]);


    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off;
    hold off
    box on
    
    
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_aggregate_old_100kHz',num2str(Va(1,ii)*1000),'mV_detune')];
    print(gcf, '-dpng', strcat(figuresdir4,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir4,frameName),'compact');
      
    close(gcf)
end

%next plot the new method

for ii = 1
    figure()
    hold on;
    for k = 12:15
        
%         plot(Ew{1,ii}(:,k),abs(dos_vals{k,ii}(:,1)), 'LineWidth', 3,'Color', colSet_blue(k,:));
        
        errorbar(ew_q,new_agg_detune(1,:),new_agg_detune(2,:),'Color',colSet_green(15,:),'Marker','o','LineWidth', 2,'MarkerFaceColor',colSet_green(15,:),'MarkerEdgeColor','k','MarkerSize',20);
        
    end
    axis square;
    box on;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.3 .4]);
    set(gca,'xtick',[.3 .35 .4]);
    
    y = ylabel('N_{t} (cm^{-3}eV^{-1})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    ylim([5e16 5e18])
    set(gca,'ytick',[1e16 1e17 1e18 1e19]);


    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off;
    hold off
    box on
        
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_aggregate_new_100kHz',num2str(Va(1,ii)*1000),'mV_detune')];
    print(gcf, '-dpng', strcat(figuresdir4,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir4,frameName),'compact');
    
    close(gcf)
end


for ii = 1
    figure()
    hold on;
    for k = 12:15
        
%         plot(Ew{1,ii}(:,k),abs(dos_vals{k,ii}(:,1)), 'LineWidth', 3,'Color', colSet_blue(k,:));
        errorbar(ew_q,old_agg_detune(1,:),old_agg_detune(2,:),'Color',colSet_green(15,:),'Marker','o','LineWidth', 2,'MarkerFaceColor',colSet_green(15,:),'MarkerEdgeColor',colSet_green(1,:),'MarkerSize',20);
        errorbar(ew_q,new_agg_detune(1,:),new_agg_detune(2,:),'Color',colSet_green(15,:),'Marker','o','LineWidth', 2,'MarkerFaceColor',colSet_green(15,:),'MarkerEdgeColor','k','MarkerSize',20);
        
    end
    axis square;
    box on;

    x = xlabel ('E_{\omega} (eV)');
    set(x,'FontName', 'Calibri');
    set(gca,'FontSize',46);
    xlim([.3 .4]);
    set(gca,'xtick',[.3 .35 .4]);
    
    y = ylabel('N_{t} (cm^{-3}eV^{-1})');
    set(y,'FontName','Calibri');
    set(gca,'FontSize',46);
    set(gca,'yscale','log');
    ylim([5e16 5e18])
    set(gca,'ytick',[1e16 1e17 1e18 1e19]);


    pbaspect([1 1.6 1]);
    set(gca,'linewidth',1.5);
    hold off;
    hold off
    box on
        
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('DOS_aggregate_tot_100kHz',num2str(Va(1,ii)*1000),'mV_detune')];
    print(gcf, '-dpng', strcat(figuresdir4,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir4,frameName),'compact');
    
    close(gcf)
end




%%
%V_bi estimation

ew_align = cell(1,length(1+offset_ref(1,1):offset_ref(1,2)));
dos_align = cell(1,length(1+offset_ref(1,1):offset_ref(1,2)));

%extract the proper DOS and Ew data, noting that we need to mirror 
for k = 12:15
    
    ew_align{1,k} = flipud(Ew_detune{1,1}(ew_idx(k,2):ew_idx(k,1),k));

    dos_align{1,k} = flipud(dos_vals_detune{k,1}(ew_idx(k,2):ew_idx(k,1),1));
    
end

%parameters and structures needed to store the V_bi calculation
depCm = 9.0e-6;

vbiEst = cell(1,length(2+offset_ref(1,1):offset_ref(1,2)));
vbi_tot = zeros(2,length(2+offset_ref(1,1):offset_ref(1,2)));

%back-calculate what V_{bi} is based on equation 22 in the Walter paper
for k = 12:15
    
    vbiEst{1,k} = 1/q.*(kBJ*T(k,1)*q*depCm).*dos_align{1,k}(:,1).*(ASderiv_norm{1,k}(ew_idx(k,2):ew_idx(k,1),2).*1e-9).^-1;
    
    vbi_tot(1,k) = mean(vbiEst{1,k});
    vbi_tot(2,k) = std(vbiEst{1,k})/length(12:15);
    
end

 disp([mean(vbi_tot(1,12:15)) std(vbi_tot(1,12:15))/length(12:15)])
 