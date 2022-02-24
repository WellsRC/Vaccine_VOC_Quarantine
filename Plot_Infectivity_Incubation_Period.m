clear;
clc;
close all;

t=linspace(0,20,1000001);
ts=[3.2 4.4 5.7 8.3];
td=ts+20;
Vac=0;
VOC=0;
R0=3;

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
subplot('Position',[0.084868421052632,0.147859922178988,0.895131578947368,0.832140077821012]);
for ii=1:4
   plot(t,R0.*ViralShedding_Symptomatic(t,td(ii),ts(ii),Vac,VOC),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',24,'XTick',0:20);
box off;
xlabel('Days post-infection','Fontsize',28);
ylabel('Infectivity','Fontsize',28);
legend({[num2str(ts(1)) '-day incubation period'],[num2str(ts(2)) '-day incubation period'],[num2str(ts(3)) '-day incubation period'],[num2str(ts(4)) '-day incubation period']},'Fontsize',24);

print(gcf,'Incubation_period_Infectivity.png','-dpng','-r600');