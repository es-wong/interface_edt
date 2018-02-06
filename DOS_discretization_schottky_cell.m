T_test = find(T(:,1) == input('Which temperature in Kelvin should be considered when determining the bounds of integration for the DOS? '));

%%

dos_fit = zeros(2,2); %array to store the boundaries of the energy range over which we will fit Gaussians to each DOS curve

dos_fit(1,1) = .2409; %upper energetic bound
dos_fit(1,2) = .1672; %lower energetic bound

% dos_fit(2,1) = .3473;
% dos_fit(2,2) = .1654;
dos_idx = zeros(t_max,2); %create an array to store the indices

%find the indices of the DOS boundary energies
for i = 1:t_max
    Ew_tmp1 = abs(Ew{1,1}(:,i) - dos_fit(1,1)); %find the index of the upper energetic bound (should be lower index) 
    Ew_tmp2 = abs(Ew{1,1}(:,i) - dos_fit(1,2)); %find the index of the lower energetic bound (should be higher index)

    [d_idx1,d_idx1] = min(Ew_tmp1);
    [d_idx2,d_idx2] = min(Ew_tmp2);

    d_closest1 = Ew{1,1}(d_idx1,i); %find the closest energy to the upper energetic bound
    d_closest2 = Ew{1,1}(d_idx2,i); %find the closest energy to the lower energetic bound
 
    dos_idx(i,1) = d_idx1; %extract the index of the upper energetic bound
    dos_idx(i,2) = d_idx2; %extract the index of the lower energetic bound
    
end;

dos_max = zeros(t_max,v_max); %array to store the maximum of the DOS
dos_peaks = zeros(t_max,v_max); %array to store the peak value of the DOS
dos_sigma = zeros(t_max,v_max); %array to store the standard deviation of the DOS 
g_idx = zeros(t_max,v_max); %array to store the first guess of the max

hm = zeros(t_max,v_max); %arry to store the half maximum of the DOS
fwhm = cell(t_max,v_max); %array to store the FWHM

dos_int = zeros(t_max,v_max); %array to store the integrated DOS at each temperature
dos_dlcp = cell(t_max,v_max); %array to store the DLCP check


dos_pts = input('How many points on either side of the DOS max should be fit as a Gaussian for the second fit? ');

E_step = input('What is the discretization energy in eV you would like to use to match the DOS and DLCP curves? ');

%%
%This part of the script performs a 2 tier Gaussian fit to the DOS data. 

for i = 5:t_max;
    for ii = 1;
    
    f1 = fit(Ew{1,ii}(dos_idx(i,1):dos_idx(i,2),i),(Nt_s{1,ii}(dos_idx(i,1):dos_idx(i,2),i)),'gauss1'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,Ew{1,ii}(dos_idx(1,1):dos_idx(1,2),i),abs(Nt_s{1,ii}(dos_idx(1,1):dos_idx(1,2),i))); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);

    dos_val = p1(2); %extract the peak fit value
    tmp = (Ew{1,ii}(dos_idx(i,1):dos_idx(i,2),i)-dos_val); %subtract the peak value from the frequency range in the experiment
    [idx idx] = min(abs(tmp)); %find the minimum, picking out the peak in the experimental data
%     g_idx(i,ii) = round((dos_idx(i,1) + idx)-1,5); %find the index of the peak in the global data range
    closest3 = Ew{1,ii}(g_idx(i,ii),i); %extract the energy of the peak

    dos_subx = Ew{1,ii}((g_idx(i,ii)-dos_pts):(g_idx(i,ii)+dos_pts),i); %extract a sub-range of x-data
    dos_suby = Nt_s{1,ii}((g_idx(i,ii)-dos_pts):(g_idx(i,ii)+dos_pts),i); %extract a sub-range of y-data
    data_max = max(dos_suby); %extract the max of the y-range
    dos_suby_prime = dos_suby./data_max; %normalize the y-data against the max

    dos_options = fitoptions('gauss1','StartPoint', [Nt_s{1,ii}(g_idx(i,ii),i), closest3, abs(Ew{1,ii}((g_idx(i,ii)+dos_pts),i)-Ew{1,ii}((g_idx(i,ii)-dos_pts),i))]); %using the first fit, input fit parameters for the second fit
    f1a = fit(dos_subx, dos_suby, 'gauss1', dos_options); %perform the second fit
    p1a = coeffvalues(f1a); %extract the coefficients of the fit
    e1a = confint(f1a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    dos_max(i,ii) = p1a(1);
    dos_peaks(i,ii) = p1a(2); %obtain the centroid of the gaussian
    dos_sigma(i,ii) =  (e1a(2,2)-e1a(1,2))/2; %calculate the std deviation associated with the centroid fit in a 'x +/- y' format, with y as the std deviation
%     dos_gauss_err(i,ii) = dos_sigma(i,ii)/dos_peaks(i,ii); %normalize the std deviation by the mean
   
    figure()
    hold on;
    plot(f1a,Ew{1,ii}(dos_idx(i,1):dos_idx(i,2),i),Nt_s{1,ii}(dos_idx(i,1):dos_idx(i,2),i)); 
    axis square
    t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
%     plot(f1a,log(AS{i,ii}(fitidx(1,1):fitidx(1,2),1)),(ASderiv_norm{ii,i}(fitidx(1,1):fitidx(1,2),1)./Y_max),'residuals')

  
    hm(i,ii) = dos_max(i,ii)/2; %find the half max of the data

    fwhm_tmp1 = abs(hm(i,ii) - Nt_s{1,ii}(1:g_idx(i,ii),i)); %find the value of the half max in the high energy side of the max
    fwhm_tmp2 = abs(hm(i,ii) - Nt_s{1,ii}(g_idx(i,ii):end,i)); %do the same for the low energy side of the data


    [f_idx1,f_idx1] = min(fwhm_tmp1); 
    [f_idx2,f_idx2] = min(fwhm_tmp2);
    
    disp([f_idx1, f_idx2]),
    
    fwhm{i,ii}(1,1) = f_idx1; %extract the index of the half max for the high energy side
    fwhm{i,ii}(1,2) = f_idx2+g_idx(i,ii); %extract the half max for the low energy side, noting that we are extracting the index in the global energy range

    fwhm{i,ii}(1,3) = Ew{1,ii}(f_idx1,i); %find the energy at the value of the half-max for the high energy side
    if i < 15;      
        fwhm{i,ii}(1,4) = Ew{1,ii}(f_idx2+g_idx(i,ii),i); %do the same for the low-energy side
    elseif i > 14; 
         fwhm{i,ii}(1,4) = Ew{1,ii}(end,i); %do the same for the low-energy side
    end;
    

    dos_int(i,ii) = quad(f1a,fwhm{i,ii}(1,4),fwhm{i,ii}(1,3),1); %integrate the DOS across the FWHM
    
    dos_dlcp{i,ii}(:,1) = E_DLCP{1,i}(1,:); %find the energies interrogated in DLCP
    dos_dlcp{i,ii}(:,2) = dos_dlcp{i,ii}(:,1) - E_step; %create a lower bound for the DLCP energy
    dos_dlcp{i,ii}(:,3) = dos_dlcp{i,ii}(:,1) + E_step; %create an upper bound for the DLCP energy
    
    for iii = 1:size(E_DLCP{1,i},2)
        dos_dlcp{i,ii}(iii,4) = quad(f1a,dos_dlcp{i,ii}(iii,2),dos_dlcp{i,ii}(iii,3),1); %integrate the DOS curve around the DLCP energies
    end;
    
        
    end;
end;

%%

%This section of the code does the same process for the illuminated curves


dos_fit = zeros(2,2); %array to store the boundaries of the energy range over which we will fit Gaussians to each DOS curve

dos_fit(1,1) = .3473; %upper energetic bound
dos_fit(1,2) = .1654; %lower energetic bound

% dos_fit(2,1) = .3473;
% dos_fit(2,2) = .1654;

%find the indices of the DOS boundary energies

Ew_tmp1 = abs(Ew{1,1}(:,T_test) - dos_fit(1,1)); %find the index of the upper energetic bound (should be lower index) 
Ew_tmp2 = abs(Ew{1,1}(:,T_test) - dos_fit(1,2)); %find the index of the lower energetic bound (should be higher index)

[d_idx1,d_idx1] = min(Ew_tmp1);
[d_idx2,d_idx2] = min(Ew_tmp2);

d_closest1 = Ew{1,1}(d_idx1,T_test); %find the closest energy to the upper energetic bound
d_closest2 = Ew{1,1}(d_idx2,T_test); %find the closest energy to the lower energetic bound

dos_idx = zeros(2,2); %create an array to store the indices
 
dos_idx(1,1) = d_idx1; %extract the index of the upper energetic bound
dos_idx(1,2) = d_idx2; %extract the index of the lower energetic bound


%%
%this part of the script does the same as above for the IR-illuminated data
%
IR_dos_fit = zeros(2,2); %array to store the boundaries of the energy range over which we will fit Gaussians to each DOS curve

IR_dos_fit(1,1) = .3473; %upper energetic bound
IR_dos_fit(1,2) = .1654; %lower energetic bound

% dos_fit(2,1) = .3473;
% dos_fit(2,2) = .1654;

%find the indices of the DOS boundary energies

Ew_tmp1 = abs(IREw{1,1}(:,T_test) - IR_dos_fit(1,1)); %find the index of the upper energetic bound (should be lower index) 
Ew_tmp2 = abs(IREw{1,1}(:,T_test) - IR_dos_fit(1,2)); %find the index of the lower energetic bound (should be higher index)

[IR_d_idx1,IR_d_idx1] = min(Ew_tmp1);
[IR_d_idx2,IR_d_idx2] = min(Ew_tmp2);

IR_d_closest1 = Ew{1,1}(IR_d_idx1,T_test); %find the closest energy to the upper energetic bound
IR_d_closest2 = Ew{1,1}(IR_d_idx2,T_test); %find the closest energy to the lower energetic bound

IR_dos_idx = zeros(2,2); %create an array to store the indices
 
IR_dos_idx(1,1) = IR_d_idx1; %extract the index of the upper energetic bound
IR_dos_idx(1,2) = IR_d_idx2; %extract the index of the lower energetic bound

IR_dos_max = zeros(t_max,v_max); %array to store the maximum of the DOS
IR_dos_peaks = zeros(t_max,v_max); %array to store the peak value of the DOS
IR_dos_sigma = zeros(t_max,v_max); %array to store the standard deviation of the DOS 
IR_g_idx = zeros(t_max,v_max); %array to store the first guess of the max

IR_hm = zeros(t_max,v_max); %arry to store the half maximum of the DOS
IR_fwhm = cell(t_max,v_max); %array to store the FWHM

IR_dos_int = zeros(t_max,v_max); %array to store the integrated DOS at each temperature
IR_dos_dlcp = cell(t_max,v_max); %array to store the DLCP check


IR_dos_pts = input('How many points on either side of the DOS max should be fit as a Gaussian for the second fit? ');

E_step = input('What is the discretization energy in eV you would like to use to match the DOS and DLCP curves? ');

%%
%This part of the script performs a 2 tier Gaussian fit to the DOS data. 

for i = 1:t_max;
    for ii = 1;
    
    f1 = fit(IREw{1,ii}(IR_dos_idx(1,1):IR_dos_idx(1,2),i),(IRNt_s{1,ii}(IR_dos_idx(1,1):IR_dos_idx(1,2),i)),'gauss1'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,Ew{1,ii}(dos_idx(1,1):dos_idx(1,2),i),abs(Nt_s{1,ii}(dos_idx(1,1):dos_idx(1,2),i))); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);

    IR_dos_val = p1(2); %extract the peak fit value
    tmp = (IREw{1,ii}(IR_dos_idx(1,1):IR_dos_idx(1,2),i)-IR_dos_val); %subtract the peak value from the frequency range in the experiment
    [idx idx] = min(abs(tmp)); %find the minimum, picking out the peak in the experimental data
    IR_g_idx(i,ii) = round((IR_dos_idx(1,1) + idx)-1,5); %find the index of the peak in the global data range
    closest4 = IREw{1,ii}(IR_g_idx(i,ii),i); %extract the energy of the peak

    IR_dos_subx = IREw{1,ii}((IR_g_idx(i,ii)-IR_dos_pts):(IR_g_idx(i,ii)+IR_dos_pts),i); %extract a sub-range of x-data
    IRdos_suby = IRNt_s{1,ii}((IR_g_idx(i,ii)-IR_dos_pts):(IR_g_idx(i,ii)+IR_dos_pts),i); %extract a sub-range of y-data
    data_max = max(IRdos_suby); %extract the max of the y-range
    dos_suby_prime = IRdos_suby./data_max; %normalize the y-data against the max

    dos_options = fitoptions('gauss1','StartPoint', [IRNt_s{1,ii}(IR_g_idx(i,ii),i), closest4, abs(IREw{1,ii}((IR_g_idx(i,ii)+IR_dos_pts),i)-IREw{1,ii}((IR_g_idx(i,ii)-IR_dos_pts),i))]); %using the first fit, input fit parameters for the second fit
    f1a = fit(IR_dos_subx, IRdos_suby, 'gauss1', dos_options); %perform the second fit
    p1a = coeffvalues(f1a); %extract the coefficients of the fit
    e1a = confint(f1a,.6827); %generate the confidence intervals assocaited with the fit (1 sigma = 68.27%)
    IR_dos_max(i,ii) = p1a(1);
    IR_dos_peaks(i,ii) = p1a(2); %obtain the centroid of the gaussian
    IR_dos_sigma(i,ii) =  (e1a(2,2)-e1a(1,2))/2; %calculate the std deviation associated with the centroid fit in a 'x +/- y' format, with y as the std deviation
%     dos_gauss_err(i,ii) = dos_sigma(i,ii)/dos_peaks(i,ii); %normalize the std deviation by the mean
   
    figure()
    hold on;
    plot(f1a,IREw{1,ii}(IR_dos_idx(1,1):IR_dos_idx(1,2),i),IRNt_s{1,ii}(IR_dos_idx(1,1):IR_dos_idx(1,2),i)); 
    axis square
    t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
%     plot(f1a,log(AS{i,ii}(fitidx(1,1):fitidx(1,2),1)),(ASderiv_norm{ii,i}(fitidx(1,1):fitidx(1,2),1)./Y_max),'residuals')
    
    IR_hm(i,ii) = IR_dos_max(i,ii)/2; %find the half max of the data

    fwhm_tmp1 = abs(IR_hm(i,ii) - IRNt_s{1,ii}(1:IR_g_idx(i,ii),i)); %find the value of the half max in the high energy side of the max
    fwhm_tmp2 = abs(IR_hm(i,ii) - IRNt_s{1,ii}(IR_g_idx(i,ii):end,i)); %do the same for the low energy side of the data
        
    [IR_f_idx1,IR_f_idx1] = min(fwhm_tmp1); 
    [IR_f_idx2,IR_f_idx2] = min(fwhm_tmp2);
        
    IR_fwhm{i,ii}(1,1) = IR_f_idx1; %extract the index of the half max for the high energy side
    IR_fwhm{i,ii}(1,2) = IR_f_idx2+IR_g_idx(i,ii); %extract the half max for the low energy side, noting that we are extracting the index in the global energy range

    IR_fwhm{i,ii}(1,3) = IREw{1,ii}(IR_f_idx1,i); %find the energy at the value of the half-max for the high energy side
    IR_fwhm{i,ii}(1,4) = IREw{1,ii}(IR_f_idx2+IR_g_idx(i,ii),i); %do the same for the low-energy side

    IR_dos_int(i,ii) = quad(f1a,IR_fwhm{i,ii}(1,4),IR_fwhm{i,ii}(1,3),1); %integrate the DOS across the FWHM
    
    IR_dos_dlcp{i,ii}(:,1) = E_DLCP{1,i}(2,:); %find the energies interrogated in DLCP
    IR_dos_dlcp{i,ii}(:,2) = IR_dos_dlcp{i,ii}(:,1) - E_step; %create a lower bound for the DLCP energy
    IR_dos_dlcp{i,ii}(:,3) = IR_dos_dlcp{i,ii}(:,1) + E_step; %create an upper bound for the DLCP energy
    
    for iii = 1:size(E_DLCP{1,i},2)
        IR_dos_dlcp{i,ii}(iii,4) = quad(f1a,IR_dos_dlcp{i,ii}(iii,2),IR_dos_dlcp{i,ii}(iii,3),1); %integrate the DOS curve around the DLCP energies
    end;
    
        
    end;
end;

