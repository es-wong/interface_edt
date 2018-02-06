%calculate the mobilities and perform a fit to grab mu_{0}
mu = zeros(t_max,2);

%calculate geometric ratio for conductance in cm
geoFact = 9e-6/4e-2;

%calculate mobilities in cm^2/(Vs)

for k = 12:15
    %mu(k,1) = (sigma_0(k,1)*geoFact)/(q*DLN_tot{k,1}(1,1));
    mu(k,1) = (Gp_w0(1,k)*geoFact)/(q*DLN_tot{k,5}(1,1));
    mu(k,2) = mu(k,1) * sqrt((DLN_tot{k,1}(1,2)/DLN_tot{k,1}(1,1))^2 + (sigma_0_err(k,1)/sigma_0(k,1))^2);
end

%%

figure() 
hold on
for k = 12:15

g1 = errorbar((T(k,1)),mu(k,1),mu(k,2),'Marker','o','MarkerFaceColor',colSet_green(k,:),'Color',colSet_green(k,:),'MarkerEdgeColor','k'); 
        set(g1,{'MarkerSize'},{20},'LineWidth',2);
end

axis square;
box on;
x = xlabel ('T (K)', 'FontSize', 46);
set(x,'FontName','Calibri');
set(gca,'FontSize',46);
xlim([250 300])
set(gca,'xtick',[250 275 300]);

y = ylabel ('\mu (cm^{2}V^{-1}s^{-1})', 'FontSize',46);
set(y,'FontName', 'Calibri')
set(gca,'yscale','log')
ylim([1e-5 8e-5])
% set(gca,'ytick',[1e-6 1e-5]);

pbaspect([1 1.6 1]);
set(gca,'linewidth',1.5);
hold off

set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
frameName = [strcat('mu_Gp_total_100kHz')];
print(gcf, '-dpng', strcat(figuresdir6,'\',frameName),'-r0');
savefig(gcf, fullfile(figuresdir6,frameName),'compact');
close(gcf)
