clear;
clc;
close all;


ts=5.723;
td=ts+20;
Vac=0;
VOC=0;

opts=optimset('TolX',10^(-16));
TPeak=fminbnd(@(x)-ViralShedding_Symptomatic(x,td,ts,0,0),4.8,5.1,opts);
tP=linspace(0,TPeak,1000001);
tA=linspace(TPeak,50,1000001);
load('MLE-Estimate-RTPCR.mat','beta');

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
subplot('Position',[0.084868421052632,0.147859922178988,0.895131578947368,0.832140077821012]);

   plot(Relative_Infection_PCR(tP,ts,Vac,VOC),PCRSens(tP,beta),'k',Relative_Infection_PCR(tA,ts,Vac,VOC),PCRSens(tA,beta),'b','LineWidth',2); hold on

set(gca,'LineWidth',2,'tickdir','out','Fontsize',24);
box off;
xlabel('Relative infectivity','Fontsize',28);
ylabel('RT-PCR diagnostic sensitivity','Fontsize',28);
legend({['Before peak infectivity (f_P)'],['After peak infectivity (f_A)']},'Fontsize',24,'Location','NorthWest');

print(gcf,'Relative_Infectivity_vs_Sensitivty.png','-dpng','-r600');