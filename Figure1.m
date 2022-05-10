clear;
clc;
close all;

addpath('Alternative_Test_Results');
figure('units','normalized','outerposition',[0 0.05 0.7 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infectivity (Unvaccinated vs Vaccinated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts=[6.3];
td=ts+20;
VOC=0;

subplot('Position',[0.105,0.7527,0.35,0.24]);

load('Peak_Infection.mat','mmv','tsv');
mm=pchip(tsv,mmv,ts);
% Prior to the peak

t_PRE=linspace(0,mm,1001);
m_nVAC=(40-20.7)/(-3.5-0);
Ct_PRE=m_nVAC.*(t_PRE-mm)+20.7;

% After the peak
t_POST=linspace(mm,20,1001);
m_nVAC=(20.7-40)/(0-7.5);
Ct_POST=m_nVAC.*(t_POST-mm)+20.7;
       
plot(Ct_PRE,ViralShedding_Symptomatic(t_PRE,inf,ts,0,0)./ViralShedding_Symptomatic(mm,inf,ts,0,0),'-','color',hex2rgb('CF3721'),'LineWidth',2); hold on
plot(Ct_POST,ViralShedding_Symptomatic(t_POST,inf,ts,0,0)./ViralShedding_Symptomatic(mm,inf,ts,0,0),'-.','color',hex2rgb('CF3721'),'LineWidth',2); hold on

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'Xminortick','off','Ytick',[0:0.2:1],'Yminortick','off','xdir','reverse');
box off;
xlabel('Ct value','Fontsize',20,'Position',[40,-0.1250,-1]);
ylabel('Relative infectivity','Fontsize',20);
legend({'Pre-peak','Post-peak'},'Fontsize',14,'Location','NorthWest');
legend boxoff;
xlim([20.7 55]);
% ylim([0 0.5])
text(65.081,0.983122362869198,'A','Fontsize',28,'FontWeight','bold');


subplot('Position',[0.56121669777,0.7527,0.35,0.24]);
t=linspace(0,20,1000001);
R0=2.79;
for ii=1:2
   plot(t,R0.*ViralShedding_Symptomatic(t,td,ts,ii-1,VOC),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:2:20,'Xminortick','off','Ytick',[0:0.1:0.5],'Yminortick','off');
box off;
xlabel('Days post-infection','Fontsize',20,'Position',[8,-0.065970462185924,-1]);
ylabel('Infectivity','Fontsize',20);
legend({'Unvaccinated','Vaccinated'},'Fontsize',14);
legend boxoff;
xlim([0 16]);
ylim([0 0.5])
text(-4.635992550338415,0.488,'B','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity (Unvaccinated vs Vaccinated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot('Position',[0.105,0.446605876393113,0.35,0.24]);

tsB=5.723;
load('Peak_Infection.mat','mmv','tsv');
mm=pchip(tsv,mmv,tsB);
% Prior to the peak
t_PRE=linspace(0,mm,1001);
S_PRE = TestSensitivity(t_PRE,tsB,[],0,0);

% After the peak
t_POST=linspace(mm,50,1001);
S_POST = TestSensitivity(t_POST,tsB,[],0,0);
       
plot(ViralShedding_Symptomatic(t_PRE,inf,tsB,0,0)./ViralShedding_Symptomatic(mm,inf,tsB,0,0),S_PRE,'-','color',hex2rgb('#486824'),'LineWidth',2); hold on
plot(ViralShedding_Symptomatic(t_POST,inf,tsB,0,0)./ViralShedding_Symptomatic(mm,inf,tsB,0,0),S_POST,'-.','color',hex2rgb('#486824'),'LineWidth',2); hold on

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'Xminortick','off','Ytick',[0:0.2:1],'Yminortick','off','Xtick',[0:0.1:1]);
box off;
xlabel({'Relative infectivity'},'Fontsize',20);
ylabel({'RT-PCR diagnostic','sensitivity'},'Fontsize',20);
legend({'Pre-peak','Post-peak'},'Fontsize',14,'Location','NorthWest');
legend boxoff;
xlim([0 1]);
% ylim([0 0.5])
text(-0.30107526881720,0.983122362869198,'C','Fontsize',28,'FontWeight','bold');


subplot('Position',[0.56121669777,0.446605876393113,0.35,0.24]);
load('Abbot PanBio_LR_Parameters.mat','beta')
for ii=1:2
   plot(t,TestSensitivity(t,ts,beta,ii-1,VOC),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:2:20,'Xminortick','off','Ytick',[0:0.1:1],'Yminortick','off');
box off;
xlabel('Days post-infection','Fontsize',20,'Position',[10,-0.131940924371849,-1]);
ylabel({'Diagnostic','sensitivity'},'Fontsize',20);
legend({'Unvaccinated','Vaccinated'},'Fontsize',14);
legend boxoff;
xlim([0 20]);
ylim([0 1])
text(-4.6978069368,0.976,'D','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% PQT (VAC vs Non_VAC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pA,~]=VariantParameters(0);
pA=pA(1:2);
R_Immunity=zeros(15,2);
for VAC=0:1
        load(['TestingonExit_Abbot PanBio_Vaccinated=' num2str(VAC) '.mat'],'VOCv','qv','IDSLA','IDSLS');
    R_Immunity(:,VAC+1)=(pA(VAC+1).*IDSLA(VOCv==0)+(1-pA(VAC+1)).*IDSLS(VOCv==0));
end


PR=Probability_Onward(R_Immunity,'Negative Binomial');



subplot('Position',[0.105,0.082,0.35,0.283687943262412]);

plot([0:14],sqrt(PR),'LineWidth',2); 
xlim([0 14]);
box off;
ytick=[0 0.002 0.01 0.05 0.1 0.15 0.2 0.25 0.35];
set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:14,'Xminortick','off','Ytick',sqrt(ytick),'yticklabel',num2str(ytick'),'Yminortick','off');
xlabel({'Duration of quarantine (days)'},'Fontsize',20)
ylabel({'Probability of','post-quarantine transmission'},'Fontsize',18);

legend({'Unvaccinated','Vaccinated'},'Fontsize',14);
legend boxoff;
text(-4.056578039996936,0.5830,'E','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Infectivity (VOC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=linspace(0,20,1000001);
ts=[6.3 5.0 4.3 3.2];
td=ts+20;
Vac=0;
R0=[2.79 4.19 5.08 6.93];



subplot('Position',[0.56121669777,0.082,0.35,0.283687943262412]);
for ii=1:4
   pp(ii)=plot(t,R0(ii).*ViralShedding_Symptomatic(t,td(ii),ts(ii),Vac,ii-1),'LineWidth',2); hold on
end

pp(1).LineStyle=':';
pp(2).LineStyle='--';
pp(3).LineStyle='-.';

pp(1).Color=hex2rgb('#1b9E77');
pp(2).Color=hex2rgb('#d95f02');
pp(3).Color=hex2rgb('#7570b3');
pp(4).Color=hex2rgb('#e7298a');

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:1:10,'Xminortick','off','Ytick',[0:0.25:2.25],'Yminortick','off');
box off;
xlabel('Days post-infection','Fontsize',20,'Position',[5,-0.34,-1]);
ylabel('Infectivity','Fontsize',20);
legend({'Original','Alpha','Delta','Omicron'},'Fontsize',14);
legend boxoff;
xlim([0 10]);
ylim([0 2.25])
text(-2.3347050,2.227,'F','Fontsize',28,'FontWeight','bold');
print(gcf,['Figure1_a.png'],'-dpng','-r300');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quarantine (Background immunity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('units','normalized','outerposition',[0 0.05 0.7 1]);

subplot('Position',[0.105,0.082,0.35,0.283687943262412]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity three doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
[pA,epsI]=VariantParameters(3);

R_Immunity=zeros(15,3);
for V=1:3
    if(V==1)
        load(['TestingonExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
    R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
end


RIm=zeros(15,12);

for jj=1:11
    vacu=(jj-1)./10;
    V_Im_Pop=(1-vacu).*epsI(1)+vacu.*epsI(2);
    
    
    w1=0;
    w2=(1-vacu).*(1-epsI(1));
    w3=vacu.*(1-epsI(2));
    w=w1+w2+w3;
    w1=w1./w;
    w2=w2./w;
    w3=w3./w;
    RIm(:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
end

RIm(:,12)=R_Immunity(:,1);

PR=Probability_Onward(RIm,'Negative Binomial');

pp=plot([0:14],sqrt(PR(:,[1 4 8 11 12])'),'LineWidth',2); 

pp(1).Color=hex2rgb('67001f');
pp(2).Color=hex2rgb('ce1256');
pp(3).Color=hex2rgb('df65b0');
pp(4).Color=hex2rgb('d4b9da');
pp(5).Color=[0.75 0.75 0.75];
pp(5).LineStyle='-.';

xlim([0 14]);
ylim([0 sqrt(0.4)]);
box off;
ytick=[0 0.002 0.01 0.05 0.1 0.15 0.2 0.3 0.4];
set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:14,'Xminortick','off','Ytick',sqrt(ytick),'yticklabel',num2str(ytick'),'Yminortick','off');
xlabel({'Duration of quarantine (days)'},'Fontsize',20)
ylabel({'Probability of','post-quarantine transmission'},'Fontsize',18);
text(10.2624703087886,0.618504947750682,'Booster uptake','Horizontalalignment','center','Fontsize',16);
legend({'0%','30%','70%','100%','No vaccination'},'Fontsize',14);
legend boxoff;
text(-4.056578039996936,0.637210474,'G','Fontsize',28,'FontWeight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Serial testing (Background immunity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('Position',[0.56121669777,0.082,0.35,0.283687943262412]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
[pA,epsI]=VariantParameters(3);

R_Immunity=zeros(15,3);
VOC=3;
for V=1:3
    if(V==1)
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(VOC) '_VAC= 0.mat'],'RTotA','RTotS','R0','td','ts');
        VAC=0;
    else
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(VOC) '_VAC= 1.mat'],'RTotA','RTotS','R0','td','ts');
        VAC=1;
    end
    R_Immunity(1:14,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
    
    RSONS=R0(VOC+1).*integral(@(t)ViralShedding_Symptomatic(t,td,ts,VAC,VOC),0,ts);
    RAONS=R0(VOC+1).*integral(@(t)ViralShedding_Asymptomatic(t,td,ts,VAC,VOC),0,inf);
    R_Immunity(15,V)=pA(V).*RAONS+(1-pA(V)).*RSONS;
end

RIm=zeros(15,1001);

for jj=1:1001
    vacu=(jj-1)./1000;
    V_Im_Pop=(1-vacu).*epsI(1)+vacu.*epsI(2);
    
    
    w1=0;
    w2=(1-vacu).*(1-epsI(1));
    w3=vacu.*(1-epsI(2));
    w=w1+w2+w3;
    w1=w1./w;
    w2=w2./w;
    w3=w3./w;
    RIm(:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
end

contourf([0:0.1:100],[1:14],RIm(1:14,:),[0.05:0.1:2.05],'LineStyle','none')
hold on
contourf([0:0.1:100],[14.35 15.35],[RIm(15,:);RIm(15,:)],[0.05:0.1:2.05],'LineStyle','none')
caxis([0.05 2.05]);

TFRE1=zeros(1001,1);
for ii=1:1001
    if(max(RIm(1:14,ii)>=1))
        TFRE1(ii)=pchip(RIm(1:14,ii),[1:14],1);
    else
        TFRE1(ii)=NaN;
    end
end
plot([0:0.1:100],TFRE1,'k','LineWidth',2);

if(min(RIm(15,:)<=1))
    plot(pchip(RIm(15,:),[0:0.1:100],1).*ones(101,1),linspace(14.35,15.35,101),'k','LineWidth',2);
end

load('Heatmap_Fig1F.mat','z');
colormap(z);

ylabel({'Frequency of','testing (days^{-1})'},'Fontsize',20,'Position',[-9.430107065426462,8.175007033348267,1.000000000000014])
xlabel({'Booster uptake'},'Fontsize',20,'Position',[50.00004768371581,-1.469395744745933,1.000000000000014]);

box off;
XTL=cell(15,1);
for ii=1:14
   XTL{ii}=num2str(ii); 
end
XTL{15}='No test';

set(gca,'LineWidth',2,'tickdir','out','Fontsize',17,'XTick',[0:10:100],'Xminortick','off','Ytick',[1:14 14.85],'YTickLabel',XTL);

h=colorbar;
h.Position=[0.91566265060241,0.083013171225937,0.01355421686747,0.279702127659575];
h.Label.String='Effective reproduction number';
h.Label.Rotation=270;
h.Label.Position=[5.2637033992343,1.050000953674316,0];
h.Ticks=[0.05 0.25:0.25:1.75 2.05];
xtickformat('percentage');
xtickangle(90);

text(-23.48606129,15,'H','Fontsize',28,'FontWeight','bold');
rmpath('Alternative_Test_Results');
print(gcf,['Figure1_b.png'],'-dpng','-r300');