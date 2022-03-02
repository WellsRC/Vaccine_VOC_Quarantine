clear;
clc;
close all;

addpath('Alternative_Test_Results');
figure('units','normalized','outerposition',[0 0.05 0.6 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infectivity (Unvaccinated vs Vaccinated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=linspace(0,20,1000001);
ts=[6.3];
td=ts+20;
VOC=0;
R0=2.79;

subplot('Position',[0.125,0.7527,0.370588732394367,0.24]);
for ii=1:2
   plot(t,R0.*ViralShedding_Symptomatic(t,td,ts,ii-1,VOC),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:2:20,'Xminortick','on','Ytick',[0:0.1:0.5],'Yminortick','on');
box off;
xlabel('Days post-infection','Fontsize',20,'Position',[8,-0.065970462185924,-1]);
ylabel('Infectivity','Fontsize',20);
legend({'Unvaccinated','Vaccinated'},'Fontsize',14);
legend boxoff;
xlim([0 16]);
ylim([0 0.5])
text(-5.324164593349167,0.488,'A','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity (Unvaccinated vs Vaccinated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('Position',[0.610915492957746,0.7527,0.370588732394367,0.24]);
load('Abbot PanBio_LR_Parameters.mat','beta')
for ii=1:2
   plot(t,TestSensitivity(t,ts,beta,ii-1,VOC),'LineWidth',2); hold on
end

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:2:20,'Xminortick','on','Ytick',[0:0.1:1],'Yminortick','on');
box off;
xlabel('Days post-infection','Fontsize',20,'Position',[10,-0.131940924371849,-1]);
ylabel({'Diagnostic','sensitivity'},'Fontsize',20);
legend({'Unvaccinated','Vaccinated'},'Fontsize',14);
legend boxoff;
xlim([0 20]);
ylim([0 1])
text(-5.945118764845,0.976,'B','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% PQT (VAC vs Non_VAC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pA=[0.351 0.9611];
R_Immunity=zeros(15,2);
for VAC=0:1
        load(['TestingonExit_Abbot PanBio_Vaccinated=' num2str(VAC) '.mat'],'VOCv','qv','IDSLA','IDSLS');
    R_Immunity(:,VAC+1)=(pA(VAC+1).*IDSLA(VOCv==0)+(1-pA(VAC+1)).*IDSLS(VOCv==0));
end


PR=Probability_Onward(R_Immunity,'Negative Binomial');

subplot('Position',[0.125,0.446605876393113,0.370588732394367,0.24]);

plot([0:14],sqrt(PR),'LineWidth',2); 
xlim([0 14]);
box off;
ytick=[0 0.002 0.01 0.05 0.1 0.15 0.2 0.25 0.35];
set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:14,'Xminortick','off','Ytick',sqrt(ytick),'yticklabel',num2str(ytick'),'Yminortick','off');
xlabel({'Duration of quarantine (days)'},'Fontsize',20)
ylabel({'Probability of','post-quarantine transmission'},'Fontsize',18);

legend({'Unvaccinated','Vaccinated'},'Fontsize',14);
legend boxoff;
text(-4.688836104513065,0.5830,'C','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Infectivity (VOC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=linspace(0,20,1000001);
ts=[6.3 5.0 4.3 3.0];
td=ts+20;
Vac=0;
R0=[2.79 4.19 5.08 6.93];

subplot('Position',[0.610915492957746,0.446605876393113,0.370588732394367,0.24]);
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

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:1:10,'Xminortick','on','Ytick',[0:0.25:2.25],'Yminortick','on');
box off;
xlabel('Days post-infection','Fontsize',20,'Position',[5,-0.34,-1]);
ylabel('Infectivity','Fontsize',20);
legend({'Original','Alpha','Delta','Omicron'},'Fontsize',14);
legend boxoff;
xlim([0 10]);
ylim([0 2.25])
text(-2.958361,2.227,'D','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quarantine (Background immunity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('Position',[0.125,0.082,0.370588732394367,0.283687943262412]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity three doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

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
epsI=[0.38 0.82];
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
text(-4.688836104513065,0.637210474,'E','Fontsize',28,'FontWeight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Serial testing (Background immunity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('Position',[0.610915492957746,0.082,0.370588732394367,0.283687943262412]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

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

RIm=zeros(15,11);

epsI=[0.38 0.82];
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


RImNV=R_Immunity(:,1);

RImT=RIm(:,[1 4 8 11]);
VV=[0 30 70 100];
for ii=1:4
    if((ii==1)||(ii==4))
        plot(linspace(1,6.75,101),pchip([1:14],RImT(1:14,ii),linspace(1,6.75,101)),'k','LineWidth',2); hold on;
        
        plot(linspace(8.25,14,101),pchip([1:14],RImT(1:14,ii),linspace(8.25,15,101)),'k','LineWidth',2);
    else
        plot(linspace(1,6.75,101),pchip([1:14],RImT(1:14,ii),linspace(1,6.75,101)),'k-.','LineWidth',1.5);
        
        plot(linspace(8.25,14,101),pchip([1:14],RImT(1:14,ii),linspace(8.25,15,101)),'k-.','LineWidth',1.5);
    end
    text(mean([6.75 8.25]),pchip([1:14],RImT(1:14,ii),mean([6.75 8.25])),[num2str(VV(ii)) '%'],'Fontsize',16,'HorizontalAlignment','center');
end

plot(linspace(1,6.75,101),pchip([1:14],RImNV(1:14),linspace(1,6.75,101)),'color',[0.75 0.75 0.75],'LineWidth',2); hold on;
plot(linspace(8.25,14,101),pchip([1:14],RImNV(1:14),linspace(8.25,15,101)),'color',[0.75 0.75 0.75],'LineWidth',2);


    text(mean([6.75 8.25]),pchip([1:14],RImNV(1:14),mean([6.75 8.25])),[{'No','vaccination'}],'Fontsize',16,'HorizontalAlignment','center','color',[0.75 0.75 0.75]);



plot([1 15.5],[1 1],'-.','color',[0.75 0.75 0.75],'LineWidth',2);

MM=min(pchip([1:14],max(RImT(1:14,:),[],2),linspace(1,14,1001)),1);
MM2=flip(pchip(1:14,min(RImT(1:14,:),[],2),linspace(1,14,1001)));
patch([linspace(1,14,1001) flip(linspace(1,14,1001))],[MM MM2],'g','facealpha',0.1,'LineStyle','none');


MM=max(pchip([1:14],max(RImT(1:14,:),[],2),linspace(1,14,1001)),1);
MM2=flip(max(pchip([1:14],min(RImT(1:14,:),[],2),linspace(1,14,1001)),1));

patch([linspace(1,14,1001) flip(linspace(1,14,1001))],[MM MM2],'r','facealpha',0.1,'LineStyle','none');

for ii=1:4
    if((ii==1)||(ii==4))
        plot(linspace(14.25,15.5,101),RImT(15,ii).*ones(101,1),'k','LineWidth',2); hold on;
    else
        plot(linspace(14.25,15.5,101),RImT(15,ii).*ones(101,1),'k-.','LineWidth',1.5);
    end
end

plot(linspace(14.25,15.5,101),RImNV(15).*ones(101,1),'color',[0.75 0.75 0.75],'LineWidth',1.5);

patch([14.25 15.5 15.5 14.25],[max(min(RImT(15,:),[],2),1) max(min(RImT(15,:),[],2),1) max(RImT(15,:),[],2) max(RImT(15,:),[],2)],'r','facealpha',0.1,'LineStyle','none');

patch([14.25 15.5 15.5 14.25],[min(max(RImT(15,:),[],2),1) min(max(RImT(15,:),[],2),1) min(RImT(15,:),[],2) min(RImT(15,:),[],2)],'g','facealpha',0.1,'LineStyle','none');

plot(14.25.*ones(101,1), linspace(0,4,101),'-.','color',[0.75 0.75 0.75],'LineWidth',2);
% text(14.4,2,'No test','Horizontalalignment','center','verticalalignment','middle','rotation',90,'Fontsize',16,'color',[0.75 0.75 0.75]);
% plot(14.5.*ones(101,1), linspace(2.5,4,101),'-.','color',[0.75 0.75 0.75],'LineWidth',2);
box off;
XTL=cell(15,1);
for ii=1:14
   XTL{ii}=num2str(ii); 
end
XTL{15}='';

text(14.896674584323046,-0.671,'No test','Horizontalalignment','center','verticalalignment','middle','rotation',90,'Fontsize',16);

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',1:1:15,'XTickLabel',XTL,'Xminortick','off','Ytick',[0:0.5:4],'Yminortick','on');
ylim([0 4]);
xlim([1 15.5])
xlabel({'Frequency of testing (days^{-1})'},'Fontsize',20,'Position',[7.5,-0.4813,-1])
ylabel({'Effective','reproduction number'},'Fontsize',20);

text(-2.8409,4.03,'F','Fontsize',28,'FontWeight','bold');
rmpath('Alternative_Test_Results');
print(gcf,['Figure1.png'],'-dpng','-r300');