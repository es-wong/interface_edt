%This script performs the freeze out analysis outlined by JW Lee et al. It fits the capacitnace curves to obtain the resistivity, which can then be combined with the DLCP data to determine the mobility  

%plot the capacitances measured at the reverse bias points as a fucntion of
%temperature to determine visually if freeze-out is occuring 

for k = 1:t_max; %index across V_{bias}
    figure(); %Plot the capacitance data grouped according to bias level to look for trends across this parameter
    hold on;
    for ii = 1:v_max; %plot the remaining temperature curves and edit the plot accordingly
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
    xlim([6e2 12e6]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([.5e-8 1.5e-8]);    
    hold off;
%     set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
%     frameName = [strcat('C_FO_',num2str(T(k,1)),'K')];
%     print(gcf, '-dpng', strcat(figuresdir5,'\',frameName),'-r0');
%     savefig(gcf, fullfile(figuresdir5,frameName),'compact');
%     close(gcf)
end;

%%
%plot G/omega to determine a similar phenomena


for k = 1:t_max; %index across V_{bias}
    figure(); %Plot the capacitance data grouped according to bias level to look for trends across this parameter
    hold on;
    for ii = 1:v_max; %plot the remaining temperature curves and edit the plot accordingly
        plot(AS{k,ii}(:,1),(AS{k,ii}(:,5)./AS{k,ii}(:,1)),'LineWidth',3,'Color',colSet_blue(k,:));   
    end;
    axis square;
    box on;
%     t = title(strcat('C vs. $\omega$,', ' ', num2str(B.(BNames{ii,1})), 'mV'), 'FontSize',40); %Title the graph
%     set(t, 'Interpreter', 'Latex', 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
    x = xlabel ('\omega (rad/s)', 'FontSize', 40); %continue formatting the axes
    set(x,'FontName','Calibri');
    set(gca,'FontSize',40);
    y = ylabel('^{G}/_{\omega}', 'FontSize', 40);
    set(y,'FontName','Calibri');
    xlim([1e4 1e5]);
    set(gca,'xscale','log');
    set(gca,'xtick',[1e3 1e4 1e5 1e6]);
%     ylim([.5e-8 1.5e-8]);    
    hold off;
    set(gcf, 'color','white', 'Position',[1 -80 1600 900], 'PaperPosition', [.25 .25 10 8], 'inverthardcopy','off')
    frameName = [strcat('G-omega_',num2str(T(k,1)),'K')];
    print(gcf, '-dpng', strcat(figuresdir5,'\',frameName),'-r0');
    savefig(gcf, fullfile(figuresdir5,frameName),'compact');
    close(gcf)
end;
