%G/omega

wTest = [540.4 603.2 929.9 1093 1288 1596];
tTest = 220:10:270;


lowFidx = zeros(length(wTest),1);


for j = 1:length(wTest)
   lowFtmp = AS{1,1}(:,1) - wTest(1,j);
   [val, lowFidx(j,1)] = min(abs(lowFtmp));
   
end

for ii = 1 
    figure();
    hold on;
    for k = 1+offset_ref(1,1):offset_ref(1,2)
%      for k = 
        plot(AS{k,ii}(:,1),(1e9*AS{k,ii}(:,8)),'LineWidth',3,'Color',colSetBlue(k,:)); 
        for j = 1:length(wTest)
            p1 = plot(AS{j+7,ii}(lowFidx(j,1),1),(1e9*AS{j+7,ii}(lowFidx(j,1),8)),'Color',colSetBlue(j+7,:),'MarkerFaceColor',colSetBlue(j+7,:),'LineWidth',2,'Marker','o','MarkerEdgeColor','k','MarkerSize',20);
        end
    end
    axis square
    box on

    
    axis square;
    box on;

    s2 = gca;
    box on
    pbaspect(s2, [1 2 1]);
    s2.LineWidth = 2;
    s2.FontSize = 44;
    s2.XLim = [4e2 1e4];
    s2.XTick = [4e2 1e4];
    s2.XTickLabel = {'10^{2}' '10^{3}' '10^{4}'};
    s2.FontName = 'Helvetica';
    s2.TickLength = [.02 .02];
    s2.XLabel.String = '\omega (rad/s)';
    s2.XScale = 'log';

    s2.YLabel.String = '^{G}/_{\omega} (nSs^{-1})';
%     s2.YScale = 'log';
    s2.YLim = [300 620];
    s2.YTickLabel = {'300' '400' '500' '600'};


%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('G-omegaSub1_',num2str(B.(BNames{ii,1})),'mV')];
%     savefig(gcf, fullfile(fig1Dir,frameName),'compact');
    
end


%%
% wTest = [540.4 603.2 929.9 1093 1288 1596 1879 2212];
% tTest = 220:10:290;

% wTest = [929.9 1093 1288 1596];
% tTest = 240:10:270;

wTestLn = log((wTest.*3)./(tTest.^(2)));
% blabla = log(ASpeaks(5:end,1)./(T(5:end,1).^2));
tTestAr = 1e3./tTest;

fTest = polyfit(tTestAr,wTestLn,1);
eTest = -fTest(1).*kB*1e3;

nu_00test(1,1) = (exp(abs(fTest(2)))/2); %calculate the attempt to escape frequency in rad/s
% nu_00test(2,ii) = (nu_00test(1,ii)*fiterr(1,ii)/2); %calculate the 1-sigma error in the attempt to escape frequency in rad/s
    
for k = 1:length(tTest)
    nu_0test(k,1) = tTest(1,k)^2*nu_00test(1,ii); %calculate the temperature-dependent attempt to escape frequency in rad/s
%     nu_0test{1,ii}(k,2) = tTest(k,1)^2*nu_00test(2,ii); %calculat the 1-sigma error in the attempt to escape frequency in rad/s
end

XvalTest = linspace(1e3*(1/200), 1e3*(1/300), 100); %create an array of x-values for the prediction line
YhatTest = fTest(1).* XvalTest + fTest(2);

figure
scatter(tTestAr,wTestLn,'Filled')
hold on
% plot(polyval(fTest,tTestAr),'color','k')
plot(XvalTest,YhatTest);

disp([eTest])


%%
ii = 1;

Xtest = zeros(1,length(tTestAr));
Ytest = zeros(1,length(wTestLn));


Xtest = tTestAr; %extract the inverse temperatures
Ytest = wTestLn; %extract ln(1/freq)
XtildeTest = [ones(length(Xtest(:,ii)),1),Xtest(:,ii)]; %create an array of ones and the x-values  for the fit specified by bias 'ii' 
Wtest(:,ii) = 1; %give the experimental data points equal weight, assuming that the error in the trend is greater than the error in each individual measurement
%     W(:,ii) = (gauss_err(1+offset_ref(1,1):offset_ref(1,2),ii).^2); %create the weighting array. Note that we are interested in 1/sigma^2 for each data point
WprimeTest(:,1) = Wtest(:,ii); %extract the variances for the fit being run

betaTest(:,ii) =  (XtildeTest.' * diag(WprimeTest)^(-1) * XtildeTest)^(-1) * (XtildeTest.' * diag(WprimeTest)^(-1) * Ytest(:,ii)); %calculate the regression coefficients
YhatTest(:,ii) = X_vals(1,:).*betaTest(2,ii) + betaTest(1,ii); %generate the array of predicted values for the plot

Yhat_residTest(:,ii) = Xtest(:,ii).*betaTest(2,ii) + betaTest(1,ii); %generate the array of predicted values for the residuals calculation

residTest(:,ii) = Ytest(:,ii) - Yhat_residTest(:,ii); %calculate the residuals
ei_sqTest = mean(residTest(:,ii).^2);

errbetaTest{1,ii} = ei_sqTest.*(XtildeTest.' *diag(WprimeTest)^(-1)* XtildeTest)^(-1); %calculate the variance in the regression coefficients
%     errbeta{1,ii} = (Xtilde.' *diag(Wprime)^(-1)* Xtilde)^(-1); %calculate the variance in the regression coefficients
fiterrTest(1,ii) = sqrt(errbetaTest{1,ii}(1,1)/(offset_ref(1,2)-offset_ref(1,1))); %normalize the error in the constant properly
fiterrTest(2,ii) = sqrt(errbetaTest{1,ii}(2,2)/(offset_ref(1,2)-offset_ref(1,1))); %normalize the error in the slope properly

E_tTest(1,ii) = -kB*betaTest(2,ii)*1e3; %calculate the trap energy
E_tTest(2,ii) = -kB*fiterrTest(2,ii)*1e3;%calculate the 1-sigma error in the energy

%%


%estimate what the ratio \omega_{0}' and \omega_{0} are based on detuning energy
%and nu_00 of original. see if this matches up with anything we have
%measured (see Lab Notebook V for details)

ii = 1;

omR = zeros(1,length(1+ offset_ref(1,1):offset_ref(1,2)));
omR0 = zeros(1,length(1+ offset_ref(1,1):offset_ref(1,2)));

for k = 1+ offset_ref(1,1):offset_ref(1,2)
    omR(1,k) = (nu_0{1,ii}(k,1))*exp(-detune_exp(2)/(kB*T(k,1))); 
    omR0(1,k) = (nu_00{1,ii}(k,1))*exp(-detune_exp(2)/(kB*T(k,1))); 

end
    

%%
yCalcT = zeros(2,1);
beta2 = zeros(2,1);
xCalcT = zeros(2,1);




yCalcT = log([2.195e5 2.87700e5]./[12889 1596])';

beta2 = log(nu_0{1,ii}(12:13,1));   

% xCalcT = 1e3./kB*T(12:13,1);

xCalcT = kB*T(12:13,1);

Etf = xCalcT.*(yCalcT - beta2);

%%

nSS = zeros(1,length(wTest));

for j = 1:length(wTest)
    nSS(1,j) = (2.5/q) * AS{j+7,ii}(lowFidx(j,1),8);
end


%%

%i need to formalize this and add it back into the original script so that
%I can make nice figures!


