clear;
clc;
close all;

t=linspace(0,20,1000001);
ts=[5.7];
td=ts+20;
Vac=1;
R0=3;

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
subplot('Position',[0.084868421052632,0.147859922178988,0.895131578947368,0.832140077821012]);
for VOC=1:3
   plot(t,TestSensitivity(t,ts,[],Vac,VOC),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',24,'XTick',0:20);
box off;
xlabel('Days post-infection','Fontsize',28);
ylabel('RT-PCR diagnostic sensitivity','Fontsize',28);
legend({'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'},'Fontsize',24);

print(gcf,'Sensitivity_VOC.png','-dpng','-r600');