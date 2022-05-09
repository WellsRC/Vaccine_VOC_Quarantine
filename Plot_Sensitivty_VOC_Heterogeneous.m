clear;
clc;
close all;

t=linspace(0,20,1000001);
ts=[6.3 5.0 4.3 3.2];
td=ts+20;
Vac=0;
VOC=0;
R0=[2.79 4.19 5.08 6.93];

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
subplot('Position',[0.084868421052632,0.147859922178988,0.895131578947368,0.832140077821012]);
for ii=1:4
   plot(t,TestSensitivity(t,ts(ii),[],Vac,ii-1),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',24,'XTick',0:20);
box off;
xlabel('Days post-infection','Fontsize',28);
ylabel('RT-PCR diagnostic sensitivity','Fontsize',28);
legend({'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'},'Fontsize',24);

print(gcf,'Sensitivity_VOC_Heterogeneous.png','-dpng','-r600');