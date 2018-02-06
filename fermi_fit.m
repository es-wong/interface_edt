
dc_dom = zeros(offset_ref(1,2) - offset_ref(1,1),b_max);
%%
for k = 1+offset_ref(1,1):offset_ref(1,2)
    for ii = 1:b_max
        
        tmp_om = AS{k,ii}(:,1) - ASpeaks(k,ii);
        [min_om, idx_om] = min(abs(tmp_om));
        dc_dom(k-offset_ref(1,1),ii) = ASderiv_norm{ii,k}(idx_om,2); 
    
    end
end

%%
    
for ii = 1 %indices for the voltage range
    figure() %create the Arrhenius plot
    g1 = plot(1./T(1+offset_ref(1,1):offset_ref(1,2),1),dc_dom(:,ii),'.','color',colSet_green(15,:)); %create the scatter plot of ln(1/f) vs. 1/T
    hold on;
%     g2 = errorbar(IRAS_inv{1,ii+v_max}(:,1),IRAS_inv{1,ii+v_max}(:,2),IRgauss_err(1+offset_IR(ii,1):offset_IR(ii,2),ii),'.','color',colSet_red(t_max,:));
        set(g1,{'markers'},{50},{'Linewidth'},{4});
%         set(g2,{'markers'},{50},{'Linewidth'},{4});
%     plot(X_vals(1,:),Yhat(:,ii),'--k','LineWidth', 3); %plot the best fit line over it
%     plot(IR_X_vals(1,:),IR_Yhat(:,ii),':k','LineWidth', 3)
    
    
    axis square;
    box on;

    x = xlabel ('T^{-1} (K^{-1})', 'FontSize', 40);
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel ('\omega ^{dC}/_{d\omega} (F)', 'FontSize',40);
    set(y,'FontName', 'Calibri')
%     xlim([X_vals(1,end) X_vals(1,1)]);
    hold off
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('derivative_amp_',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)
end

%%
%define the fermi model for fitting the data
fermi = fittype('A./(1+exp(E.*x_fermi))','independent',{'x_fermi'},'coefficients',{'A','E'},'dependent',{'f_fermi'});


x_fermi = 1./(kB*T(1+offset_ref(1,1):12,1));
f_fermi = dc_dom(1:8,1);

%approximate the fermi fit by first fitting a simple exponential
f6 = fit(x_fermi,f_fermi,'exp1');

figure
plot(f6,(1./(kB.*T(1+offset_ref(1,1):offset_ref(1,2),1))),dc_dom(:,1));

detune_exp = coeffvalues(f6); %parameter estimates
detune_err_exp = confint(f6,.6827); %one sigma confidence interval on the parameter estimates for the exponential approximation
detune_energy_error = detune_err_exp(1,2) - detune_exp(1,2); %calculate the error in the energy estimate

%%

%run the fermi fit
f7 = fit(x_fermi,f_fermi,fermi,'problem',{},'StartPoint',[19*dc_dom(end,1) abs(detune_exp(2))]);

figure
plot(f7,(1./(kB.*T(1+offset_ref(1,1):offset_ref(1,2),1))),dc_dom(:,1));

detune_fermi = coeffvalues(f7);
detune_err_fermi = confint(f7,.6827); %one sigma confidence interval on the parameter estimates for the fermi fit


%%

x_fermi_vals = linspace(1./(kB*T(offset_ref(1,1),1)),1./(kB*T(offset_ref(1,2),1)));
f_fermi_vals_exp = detune_exp(1).*exp(x_fermi_vals.*detune_exp(2));
% f_fermi_vals_fermi = detune_fermi(1)./(1 + exp(detune_fermi(2).*x_fermi_vals));
%%
for ii = 1 
    figure() 
    g1 = plot(1./(kB*T(1+offset_ref(1,1):offset_ref(1,2),1)),dc_dom(:,ii),'.','color',colSet_green(t_max,:)); %create the scatter plot of ln(1/f) vs. 1/T
    hold on;
%     g2 = errorbar(IRAS_inv{1,ii+v_max}(:,1),IRAS_inv{1,ii+v_max}(:,2),IRgauss_err(1+offset_IR(ii,1):offset_IR(ii,2),ii),'.','color',colSet_red(t_max,:));
        set(g1,{'markers'},{50},{'Linewidth'},{4});
%         set(g2,{'markers'},{50},{'Linewidth'},{4});
    plot(x_fermi_vals(1,:),f_fermi_vals_exp(1,:),'--r','LineWidth', 3); %plot the best fit line over it
%     plot(x_fermi_vals(1,:),f_fermi_vals(1,:),':k','LineWidth', 3)
    
    
    axis square;
    box on;
%     t = title (strcat('Arrhenius Plot,', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize', 40);
%     set(t,'Interpreter', 'Latex')
    x = xlabel ('(k_{B}T)^{-1} (eV^{-1})', 'FontSize', 46);
    set(x,'FontName','Calibri');
    xlim([42 60]);
    set(gca,'xtick',[42 51 60])
    
    set(gca,'FontSize',46);
    y = ylabel ('\omega ^{dC}/_{d\omega} (F)', 'FontSize',46);
    set(y,'FontName', 'Calibri')
    
    pbaspect([1 1.6 1])
    set(gca,'linewidth',1.5);
    hold off
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('amp_fermi_fit',num2str(B.(BNames{ii,1})),'mV')];
%     print(gcf, '-dpng', strcat(figuresdir1,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir1,frameName),'compact');
%     close(gcf)
end
