% This script combines the DLCP carrier data and the AS energy data to
% generate the density of states plots for each voltage in the admittance
% measurement. 

V_bi = input('What is the built-in voltage of the solar cell? ');
% N_free = input('What is the free carrier density in m^{-3} in the dark as determined by DLCP? ');
% IRN_free = input('What is the free carrier density in m^{-3} under illumination as determined by DLCP? ');

%%
% This section of the code calculates the depletion width W for each of the
% applied voltages using the free carrier density determined in the
% DLCP measurements. Similarly, it calculates the Fermi Level E_{F} from
% the intrinsic levels E_{i} and n_{i}.


%input constants needed for the calculation

A = 4e-6;  %Define constants that will be used later in our looped calculation. Note that the semicondcutor radius is in m, as are the device area (m^2) and volume (m^3). Energy parameters are in eV.  
epsilon = 10;
eps_0 = 8.8541878e-12;
q = 1.60217657e-19;
r_nc = 1.55e-9;
v_nc = (4/3)*pi*(r_nc)^3;
chi = .6;  %packing fraction of the film
v_film = A*3e-7;
N_particles = (v_film/v_nc)*chi;   
N_v = (2/v_nc)*chi;  %density of states for a single nanocrystal
E_v = -5.3;
M_0 = 9.10938356e-31; %electron mass in kg
M_e = .3*M_0;
M_h = .34*M_0;
E_c = -3.8;
kB = 8.6173324*(10^-5); %Boltzmann constant in eV

figuresdir4 = 'C:\Users\Eric\Desktop\Admittance\ZnO-PbS\light_soak_study\IR_soak\try18\600mW_3\DOS';

W = zeros(t_max,v_max); %initialize an array to hold the depletion with W as a function of applied bias
IRW = zeros(t_max,v_max); %initialize an array to hold the depletion with W as a function of applied bias
W_error = zeros(t_max,v_max); %initialize an array to hold the error in the depletion with W as a function of applied bias
IRW_error = zeros(t_max,v_max); %initialize an array to hold the error in the depletion with W as a function of applied bias
Ei = zeros(t_max,1); %initialize an array to store the intrinsic energy E_i as a function of temperature
ni = zeros(t_max,1); %intialize an array to store the intrinsic carrier numner n_i as afunction of temperature
Ef = zeros(t_max,1); %initialize an array to store the Fermi level as a function of temperature
IREf = zeros(t_max,1);
Va = zeros(1,v_max); %initialize an array to store the applied bias levels

for i = 1: t_max;
    for ii = 1:v_max; %index over voltages
        Va(1,ii) = (B.(BNames{ii,1})/1000); %extract the applied voltages in the AS measurement in volts
%         W(i,ii) = sqrt(2*epsilon*eps_0*(V_bi - Va(1,ii))./(q*N_free(i,1)*1e6)); %calculate the depletion width W in meters
%         IRW(i,ii) = sqrt(2*epsilon*eps_0*(V_bi - Va(1,ii))./(q*IRN_free(i,1)*1e6)); %calculate the depletion width W in meters


        W(i,ii) = 180e-9;
        IRW(i,ii) = 180e-9;
        W_error(i,ii) = (1/2)*(N_free(i,2)/N_free(i,1))*W(i,ii);
        IRW_error(i,ii) = (1/2)*(IRN_free(i,2)/IRN_free(i,1))*IRW(i,ii);
    end;
end;


for i = 1:t_max; %index over tempertuares
    Ei(i,1) = (E_c + E_v)/2 + .75*kB*(T(i,1))*reallog(M_h/M_e); %calculate E_{i} for each temperature 
    ni(i,1) = N_v*exp((E_v-Ei(i,1))/(kB*T(i,1))); %calculate n_{i} for each temperature 
    Ef(i,1) = -kB*T(i,1)*reallog((N_free(i,1))/ni(i,1))+Ei(i,1); %calculate E_{F} for each temperature 
    IREf(i,1) = -kB*T(i,1)*reallog((IRN_free(i,1))/ni(i,1))+Ei(i,1);
end;



%%

%This section of the code calculates the measurement energy E_{\omega}.
%Note that the calculation requires a conversion from eV into Joules and
%then back. 

Ew = cell(1,v_max);
Ew_j = cell(1,v_max);
IREw = cell(1,v_max);
IREw_j = cell(1,v_max);
% Ew = zeros(size(AS{1,1},1),t_max); %initialize an array to sore the E_{\omega} values as a function of temperature and frequency
Efn = zeros(t_max,1);
IREfn = zeros(t_max,1);
kBJ = 1.3806485e-23; %Boltzmann constant in joules
eV_j = 1.60218e-19; %conversion factor for eV to joules

% for iii = 1:size(AS{1,1},1);
    for ii = 1:v_max %index over voltages
        for i = 1:t_max; %index over temperatures
            Ew{1,ii}(:,i) = kB*T(i,1)*log((nu_0(1,ii))./AS{i,ii}(:,1)); %calculate E_{\omega} in eV using xi_{0} (this includes temperature dependence)
            Efn(i,1) = Ef(i,1) - E_v; %calculate E_{fn(inf)} 
            
            IREw{1,ii}(:,i) = kB*T(i,1)*log((IR_nu_0(1,ii))./IRAS{i,ii}(:,1)); %calculate E_{\omega} in eV using xi_{0} (this includes temperature dependence)
            IREfn(i,1) = IREf(i,1) - E_v; %calculate E_{fn(inf)} 
       end;
    end;
% end;

for ii =1:v_max;    
    Ew_j{1,ii}(:,:) = Ew{1,ii}(:,:).*eV_j; %convert E_{\omega} into joules and store in a separate array
    Efn_j = Efn.*eV_j; %convert E_{fn(inf)} into joules and store as a separate array
    
    IREw_j{1,ii}(:,:) = IREw{1,ii}(:,:).*eV_j; %convert E_{\omega} into joules and store in a separate array
    IREfn_j = IREfn.*eV_j; %convert E_{fn(inf)} into joules and store as a separate array
end;

%%

%This section of the code calculates the trap density N_{t} and plots it
%against the measurement energy E_{\omega} calculated above. Carefully note
%the units in the calculation below.

        
Nt = cell(2,v_max); %intialize an array to store the trap density N_{t}
Nt_s = cell(1,v_max); %initialize an array to store the smoothed trap data
Nt_int = cell(1,v_max); %initialize an array to store the integrated trap density

IRNt = cell(2,v_max); %intialize an array to store the trap density N_{t}
IRNt_s = cell(1,v_max); %initialize an array to store the smoothed trap data
IRNt_int = cell(1,v_max); %initialize an array to store the integrated trap density

for ii = 1:v_max;
    for i = 1:t_max;
        Nt{1,ii}(:,i) = ((V_bi - Va(1,ii))^2)./(A*W(1,ii).*(q*V_bi-(Efn_j(i,1) - Ew_j{1,ii}(:,i)))).*ASderiv_norm{ii,i}(:,1)./(kB*T(i,1)); %calculate N_t,the density of trap states
        
        IRNt{1,ii}(:,i) = ((V_bi - Va(1,ii))^2)./(A*IRW(1,ii).*(q*V_bi-(IREfn_j(i,1) - IREw_j{1,ii}(:,i)))).*IRASderiv_norm{ii,i}(:,1)./(kB*T(i,1));
    
    end;
end;


for i = 1:t_max; %index over temperatures
    for ii = 1:v_max; %index over bias points
        Nt{2,ii}(:,:) = Nt{1,ii}(:,:).*1e-6; %put the trap density in cm^{-3}
        Nt_s{1,ii}(:,i) = smooth(Nt{2,ii}(:,i),9,'sgolay'); %smooth the data for plotting by using Savitzky-Golay smoothing
        Nt_int{1,ii}(:,i) = trapz(Ew{1,ii}(:,i),Nt_s{1,ii}(:,i)); %calculate the total density of traps at each temperature by integrating the area of each curve
        
        IRNt{2,ii}(:,:) = IRNt{1,ii}(:,:).*1e-6; %put the trap density in cm^{-3}
        IRNt_s{1,ii}(:,i) = smooth(IRNt{2,ii}(:,i),9,'sgolay'); %smooth the data for plotting by using Savitzky-Golay smoothing
        IRNt_int{1,ii}(:,i) = trapz(IREw{1,ii}(:,i),IRNt_s{1,ii}(:,i)); %calculate the total density of traps at each temperature by integrating the area of each curve
        
    end;
end;

%%
for ii = 1:v_max %index over bias points
    figure()
    hold on;
    for i = 1:t_max; %index over temperatures. Note that we only go to the offeset point in specified in the AS script
        plot(Ew{1,ii}(:,i),abs(Nt_s{1,ii}(:,i)), 'LineWidth', 3,'Color', colSet2(i,:));
        plot(IREw{1,ii}(:,i),abs(IRNt_s{1,ii}(:,i)), 'LineWidth', 3,'Color', colSet1(i,:));
    end;
    axis square;
    box on;
%     t = title(strcat('Density of Trap States vs. Measurement Energy, ', num2str(Va(1,ii)*1000), 'mV'));
%     set(t,'Interpreter','Latex') %Set the font to LaTex font
    x = xlabel ('$E_{\omega}$ (eV)');
    set(x,'Interpreter', 'Latex');
    set(gca,'FontSize',50);
    y = ylabel('$N_{t}$, ($cm^{-3}eV^{-1}$)');
    set(y,'Interpreter','Latex');
    set(gca,'FontSize',50);
    set(gca,'yscale','log');
    xlim([.1 .35]);
    ylim([1e14 2e20]);
    set(gca,'ytick',[1e15 1e17 1e19 1e20]);
    hold off
    box on
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('DOS_both_all_',num2str(Va(1,ii)*1000),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir4,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir4,frameName),'compact');
%     close(gcf)
end;

  
%%

%This section of the code integrates the density of states of both dark and
%illuminated IR by performing a Gaussian fit and then integrating over the
%FWHM of the fit

T_test = find(T(:,1) == input('Which temperature in Kelvin should be considered when determining the bounds of integration for the DOS? '));


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

for i = 1:t_max;
    for ii = 1;
    
    f1 = fit(Ew{1,ii}(dos_idx(1,1):dos_idx(1,2),i),(Nt_s{1,ii}(dos_idx(1,1):dos_idx(1,2),i)),'gauss1'); %fit the normalized derivative data with a gaussian model
    p1 = coeffvalues(f1); %extract the coefficients of each value
    
%     figure()
%     hold on
%     plot(f1,Ew{1,ii}(dos_idx(1,1):dos_idx(1,2),i),abs(Nt_s{1,ii}(dos_idx(1,1):dos_idx(1,2),i))); %plot the figure for inspection
%     axis square
%     t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,1'), 'FontSize',30);

    dos_val = p1(2); %extract the peak fit value
    tmp = (Ew{1,ii}(dos_idx(1,1):dos_idx(1,2),i)-dos_val); %subtract the peak value from the frequency range in the experiment
    [idx idx] = min(abs(tmp)); %find the minimum, picking out the peak in the experimental data
    g_idx(i,ii) = round((dos_idx(1,1) + idx)-1,5); %find the index of the peak in the global data range
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
    plot(f1a,Ew{1,ii}(dos_idx(1,1):dos_idx(1,2),i),Nt_s{1,ii}(dos_idx(1,1):dos_idx(1,2),i)); 
    axis square
    t = title(strcat('T =', num2str(T(i,1)), 'K,', num2str(B.(BNames{ii,1})), 'mV,2'), 'FontSize',30);
%     plot(f1a,log(AS{i,ii}(fitidx(1,1):fitidx(1,2),1)),(ASderiv_norm{ii,i}(fitidx(1,1):fitidx(1,2),1)./Y_max),'residuals')
    
    hm(i,ii) = dos_max(i,ii)/2; %find the half max of the data

    fwhm_tmp1 = abs(hm(i,ii) - Nt_s{1,ii}(1:g_idx(i,ii),i)); %find the value of the half max in the high energy side of the max
    fwhm_tmp2 = abs(hm(i,ii) - Nt_s{1,ii}(g_idx(i,ii):end,i)); %do the same for the low energy side of the data
        
    [f_idx1,f_idx1] = min(fwhm_tmp1); 
    [f_idx2,f_idx2] = min(fwhm_tmp2);
        
    fwhm{i,ii}(1,1) = f_idx1; %extract the index of the half max for the high energy side
    fwhm{i,ii}(1,2) = f_idx2+g_idx(i,ii); %extract the half max for the low energy side, noting that we are extracting the index in the global energy range

    fwhm{i,ii}(1,3) = Ew{1,ii}(f_idx1,i); %find the energy at the value of the half-max for the high energy side
    fwhm{i,ii}(1,4) = Ew{1,ii}(f_idx2+g_idx(i,ii),i); %do the same for the low-energy side

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

IR_dos_fit = zeros(2,2); %array to store the boundaries of the energy range over which we will fit Gaussians to each DOS curve

IR_dos_fit(1,1) = .3113; %upper energetic bound
IR_dos_fit(1,2) = .1379; %lower energetic bound

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

%%



figure()
hold on;

plot(T(:,1),N_free(:,1),'Color',colSet2(t_max,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet2(t_max,:), 'MarkerSize',40);        
plot(T(:,1),dos_int(:,1),'Color',colSet2(t_max-5,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet2(t_max-5,:), 'MarkerSize',40);        


plot(T(:,1),IRN_free(:,1),'Color',colSet1(t_max,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);
plot(T(:,1),IR_dos_int(:,1),'Color',colSet1(t_max-5,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);


axis square;
t = title(strcat('$N_{DL}$ vs. Temperature'), 'FontSize',40); %Title the graph
set(t,'Interpreter','Latex') %Set the font to LaTex font
x = xlabel ('Temperature, $<T>$ ($K$)');
set(x,'Interpreter', 'Latex');
set(gca,'FontSize',40);
y = ylabel ('$N_{DL}$ ($cm^{-3}$)');
set(y,'Interpreter','Latex');
set(gca,'FontSize',40);
set(gca,'yscale','log');
hold off
%     xlim([.1 .35]);
xlim([235 335]); 
box on
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL-T_log_both_',num2str(T(i,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)


%%


figure()
hold on;

% plot(T(:,1),N_free(:,1),'Color',colSet2(t_max,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet2(t_max,:), 'MarkerSize',40);        
plot(T(:,1),W(:,1),'Color',colSet2(t_max-5,:),'Marker','.','LineWidth',3,'MarkerFaceColor',colSet2(t_max-5,:), 'MarkerSize',40);        


% plot(T(:,1),IRN_free(:,1),'Color',colSet1(t_max,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);
plot(T(:,1),IRW(:,1),'Color',colSet1(t_max-5,:),'Marker','d','LineWidth',3, 'MarkerFaceColor','w', 'MarkerSize',10);


axis square;
t = title(strcat('$N_{DL}$ vs. Temperature'), 'FontSize',40); %Title the graph
set(t,'Interpreter','Latex') %Set the font to LaTex font
x = xlabel ('Temperature, $<T>$ ($K$)');
set(x,'Interpreter', 'Latex');
set(gca,'FontSize',40);
y = ylabel ('$N_{DL}$ ($cm^{-3}$)');
set(y,'Interpreter','Latex');
set(gca,'FontSize',40);
set(gca,'yscale','log');
hold off
%     xlim([.1 .35]);
xlim([235 335]); 
box on
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('NDL-T_log_both_',num2str(T(i,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf

