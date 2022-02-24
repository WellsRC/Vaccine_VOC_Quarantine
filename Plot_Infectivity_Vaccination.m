clear;
clc;
close all;

t=linspace(0,20,1000001);
ts=[5.7];
td=ts+20;
Vac=0;
R0=3;

for VOC=0:3
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    subplot('Position',[0.084868421052632,0.147859922178988,0.895131578947368,0.832140077821012]);
    for ii=1:2
       plot(t,R0.*ViralShedding_Symptomatic(t,td,ts,ii-1,VOC),'LineWidth',2); hold on
    end

    set(gca,'LineWidth',2,'tickdir','out','Fontsize',24,'XTick',0:20);
    box off;
    xlabel('Days post-infection','Fontsize',28);
    ylabel('Infectivity','Fontsize',28);
    legend({'Unvaccinated','Vaccinated'},'Fontsize',24);
    xlim([0 15]);
    print(gcf,['Vaccination_Infectivity_VOC=' num2str(VOC) '.png'],'-dpng','-r600');
end