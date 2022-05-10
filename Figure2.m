clear;
clc;
close all;

% addpath('Alternative_Test_Results');
figure('units','normalized','outerposition',[0 0.05 0.6 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Quarantine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Background immunity two doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subplot('Position',[0.101232394366199,0.5947,0.39,0.39]);

RIm=zeros(4,15,101);
for VOC=0:3
    [pA,epsI]=VariantParameters(VOC);
    R_Immunity=zeros(15,3);
    for V=1:3
        if(V==1)
            load(['TestingonExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
        else
            load(['TestingonExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
        end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==VOC)+(1-pA(V)).*IDSLS(VOCv==VOC));
    end


    for jj=1:101
        vacu=(jj-1)./100;
        V_Im_Pop=vacu.*epsI(1);


        w1=1-vacu;
        w2=vacu.*(1-epsI(1));
        w3=0;
        w=w1+w2+w3;
        w1=w1./w;
        w2=w2./w;
        w3=w3./w;
        RIm(VOC+1,:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
    end
end

PR=Probability_Onward(RIm,'Negative Binomial');

PR_Threshold=squeeze(PR(1,qv(VOCv==0)==7,1));

QD=zeros(5,101);
qq=unique(qv);
for VOC=0:3
    for ii=1:101
        qd=qq(PR(VOC+1,:,ii)<=PR_Threshold);
        if(~isempty(qd))
           QD(VOC+1,ii)=min(qd); 
        else
           QD(VOC+1,ii)=15; 
        end
    end
end
QD(5,:)=7.*ones(1,101);
dy=0.1;
QD(2,:)=QD(2,:)+dy;
QD(3,:)=QD(3,:)+2.*dy;
QD(4,:)=QD(4,:)+3.*dy;
ppt=plot([0:1:100],QD(5,:),':','color',[0.75 0.75 0.75],'LineWidth',2); hold on;
pp=plot([0:1:100],QD(1:4,:),'LineWidth',2); 

pp(1).LineStyle=':';
pp(2).LineStyle='--';
pp(3).LineStyle='-.';

pp(1).Color=hex2rgb('#1b9E77');
pp(2).Color=hex2rgb('#d95f02');
pp(3).Color=hex2rgb('#7570b3');
pp(4).Color=hex2rgb('#e7298a');


set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:10:100,'Xminortick','off','Ytick',[0:1:14],'Yminortick','off');
box off;
xlabel('Two-dose uptake','Fontsize',20,'Position',[50 -2 -1]);
ylabel({'Minimum quarantine','duration (days)'},'Fontsize',20);
legend([pp;ppt],{'Original','Alpha','Delta','Omicron','Original (No vaccination)'},'Fontsize',14);
legend boxoff;
xlim([0 100]);
ylim([0 14])
xtickformat('percentage');
xtickangle(90);
text(-22.25418716671259,13.818,'A','Fontsize',28,'FontWeight','bold');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Background immunity three doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.598591549295775,0.5947,0.39,0.39]);


RIm=zeros(4,15,101);
for VOC=0:3
    [pA,epsI]=VariantParameters(VOC);
    R_Immunity=zeros(15,3);
    for V=1:3
        if(V==1)
            load(['TestingonExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
        else
            load(['TestingonExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
        end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==VOC)+(1-pA(V)).*IDSLS(VOCv==VOC));
    end


    for jj=1:101
        vacu=(jj-1)./100;
        V_Im_Pop=(1-vacu).*epsI(1)+vacu.*epsI(2);


        w1=0;
        w2=(1-vacu).*(1-epsI(1));
        w3=vacu.*(1-epsI(2));
        w=w1+w2+w3;
        w1=w1./w;
        w2=w2./w;
        w3=w3./w;
        RIm(VOC+1,:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
    end
end

PR=Probability_Onward(RIm,'Negative Binomial');

QD=zeros(5,101);
qq=unique(qv);
for VOC=0:3
    for ii=1:101
        qd=qq(PR(VOC+1,:,ii)<=PR_Threshold);
        if(~isempty(qd))
           QD(VOC+1,ii)=min(qd); 
        else
           QD(VOC+1,ii)=15; 
        end
    end
end

QD(5,:)=7.*ones(1,101);

dy=0.1;
QD(2,:)=QD(2,:)+dy;
QD(3,:)=QD(3,:)+2.*dy;
QD(4,:)=QD(4,:)+3.*dy;

ppt=plot([0:1:100],QD(5,:),':','color',[0.75 0.75 0.75],'LineWidth',2); hold on;
pp=plot([0:1:100],QD(1:4,:),'LineWidth',2); 
pp(1).LineStyle=':';
pp(2).LineStyle='--';
pp(3).LineStyle='-.';

pp(1).Color=hex2rgb('#1b9E77');
pp(2).Color=hex2rgb('#d95f02');
pp(3).Color=hex2rgb('#7570b3');
pp(4).Color=hex2rgb('#e7298a');

set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:10:100,'Xminortick','off','Ytick',[0:1:14],'Yminortick','off');
box off;
xlabel('Booster uptake','Fontsize',20,'Position',[50 -2 -1]);
ylabel({'Minimum quarantine','duration (days)'},'Fontsize',20);

legend([pp;ppt],{'Original','Alpha','Delta','Omicron','Original (No vaccination)'},'Fontsize',14);
legend boxoff;
xlim([0 100]);
ylim([0 14])
xtickformat('percentage');
xtickangle(90);
text(-22.25418716671259,13.818,'B','Fontsize',28,'FontWeight','bold');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% Serial testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Background immunity two doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.101232394366199,0.103657548125633,0.39,0.39]);

RIm=zeros(4,15,101);
for VOC=0:3
    
    R_Immunity=zeros(15,3);
    [pA,epsI]=VariantParameters(VOC);
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

    for jj=1:101
        vacu=(jj-1)./100;
        V_Im_Pop=vacu.*epsI(1);


        w1=(1-vacu);
        w2=vacu.*(1-epsI(1));
        w3=0;
        w=w1+w2+w3;
        w1=w1./w;
        w2=w2./w;
        w3=w3./w;
        RIm(VOC+1,:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
    end
end

Freq=[1:15];
TF=zeros(4,101);
for VOC=0:3
    for ii=1:101
        qd=Freq(RIm(VOC+1,:,ii)<1);
        if(~isempty(qd))
           TF(VOC+1,ii)=max(qd); 
        else
           TF(VOC+1,ii)=0; 
        end
    end
end


dy=0.1;
TF(2,:)=TF(2,:)-dy;
TF(3,:)=TF(3,:)-2.*dy;
TF(4,:)=TF(4,:)-3.*dy;

% patch([0 0 100 100],[-3.*dy 0.5 0.5 -3.*dy],'k','facealpha',0.1,'LineStyle','none'); hold on;

for ii=1:4
    X=[0:1:100];
    Y=TF(ii,:);
    if(~isempty(Y(Y<=14)))
        pp(ii)=plot(X(Y<=14),Y(Y<=14),'LineWidth',2); hold on;
        pp2(ii)=pp(ii);
        if(~isempty(Y(Y>14)))
            Xt=X(Y<=14);
            Yt=Y(Y<=14);
            ss(ii)=scatter(Xt(end),Yt(end),40);
        
            pp2(ii)=plot(X(Y>14),Y(Y>14),'LineWidth',2); hold on;
            Xt=X(Y>14);
            Yt=Y(Y>14);
            ss2(ii)=scatter(Xt(1),Yt(1),40,'filled'); hold on;
        end
    else        
        pp2(ii)=plot(X(Y>14),Y(Y>14),'LineWidth',2); hold on;
        pp(ii)=pp2(ii);
    end
end


pp(1).LineStyle=':';
pp(2).LineStyle='--';
pp(3).LineStyle='-.';

pp(1).Color=hex2rgb('#1b9E77');
pp(2).Color=hex2rgb('#d95f02');
pp(3).Color=hex2rgb('#7570b3');
pp(4).Color=hex2rgb('#e7298a');

pp2(1).LineStyle=':';
pp2(2).LineStyle='--';
pp2(3).LineStyle='-.';

pp2(1).Color=hex2rgb('#1b9E77');
pp2(2).Color=hex2rgb('#d95f02');
pp2(3).Color=hex2rgb('#7570b3');
pp2(4).Color=hex2rgb('#e7298a');

ss2(1).MarkerEdgeColor=hex2rgb('#1b9E77');
ss2(2).MarkerEdgeColor=hex2rgb('#d95f02');
% ss2(3).MarkerEdgeColor=hex2rgb('#7570b3');
% ss2(4).MarkerEdgeColor=hex2rgb('#e7298a');

ss2(1).MarkerFaceColor=hex2rgb('#1b9E77');
ss2(2).MarkerFaceColor=hex2rgb('#d95f02');
% ss2(3).MarkerFaceColor=hex2rgb('#7570b3');
% ss2(4).MarkerFaceColor=hex2rgb('#e7298a');

ss(1).MarkerEdgeColor=hex2rgb('#1b9E77');
ss(2).MarkerEdgeColor=hex2rgb('#d95f02');
% ss(3).MarkerEdgeColor=hex2rgb('#7570b3');
% ss(4).MarkerEdgeColor=hex2rgb('#e7298a');


set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:10:100,'Xminortick','off','Ytick',[0:1:15],'Yminortick','off','YTickLabel',{'','1','2','3','4','5','6','7','8','9','10','11','12','13','14','No test'},'YDir','reverse');
box off;
xlabel('Two-dose uptake','Fontsize',20,'Position',[50 17.7 -1]);
ylabel({'Minimum frequency','of testing  (days^{-1})'},'Fontsize',20,'Position',[-9.231752673157025,7.5,-0.1]);
legend boxoff;
legend(pp,{'Original','Alpha','Delta','Omicron'},'Fontsize',14,'location','SouthWest');

xlim([0 100]);
ylim([-3.*dy 15])
xtickformat('percentage');
xtickangle(90);

text(-25.8659,-0.218,'C','Fontsize',28,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Background immunity two doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.598591549295775,0.103657548125633,0.39,0.39]);

RIm=zeros(4,15,101);
for VOC=0:3
    
    
    R_Immunity=zeros(15,3);
    [pA,epsI]=VariantParameters(VOC);
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

    for jj=1:101
        vacu=(jj-1)./100;
        V_Im_Pop=(1-vacu).*epsI(1)+vacu.*epsI(2);


        w1=0;
        w2=(1-vacu).*(1-epsI(1));
        w3=vacu.*(1-epsI(2));
        w=w1+w2+w3;
        w1=w1./w;
        w2=w2./w;
        w3=w3./w;
        RIm(VOC+1,:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
    end
end

Freq=[1:15];
TF=zeros(4,101);
for VOC=0:3
    for ii=1:101
        qd=Freq(RIm(VOC+1,:,ii)<1);
        if(~isempty(qd))
           TF(VOC+1,ii)=max(qd); 
        else
           TF(VOC+1,ii)=0; 
        end
    end
end


dy=0.1;
TF(2,:)=TF(2,:)-dy;
TF(3,:)=TF(3,:)-2.*dy;
TF(4,:)=TF(4,:)-3.*dy;

% patch([0 0 100 100],[-3.*dy 0.5 0.5 -3.*dy],'k','facealpha',0.1,'LineStyle','none'); hold on;
for ii=1:4
    X=[0:1:100];
    Y=TF(ii,:);
    if(~isempty(Y(Y<=14)))
        pp(ii)=plot(X(Y<=14),Y(Y<=14),'LineWidth',2); hold on;
        
        if(~isempty(Y(Y>14)))
            Xt=X(Y<=14);
            Yt=Y(Y<=14);
            ss(ii)=scatter(Xt(end),Yt(end),40);
        
            pp2(ii)=plot(X(Y>14),Y(Y>14),'LineWidth',2); hold on;
            Xt=X(Y>14);
            Yt=Y(Y>14);
            ss2(ii)=scatter(Xt(1),Yt(1),40,'filled'); hold on;
        end
    else        
        pp2(ii)=plot(X(Y>14),Y(Y>14),'LineWidth',2); hold on;
        pp(ii)=pp2(ii);
    end
end


pp(1).LineStyle=':';
pp(2).LineStyle='--';
pp(3).LineStyle='-.';

pp(1).Color=hex2rgb('#1b9E77');
pp(2).Color=hex2rgb('#d95f02');
pp(3).Color=hex2rgb('#7570b3');
pp(4).Color=hex2rgb('#e7298a');

pp2(1).LineStyle=':';
pp2(2).LineStyle='--';
pp2(3).LineStyle='-.';

pp2(1).Color=hex2rgb('#1b9E77');
pp2(2).Color=hex2rgb('#d95f02');
pp2(3).Color=hex2rgb('#7570b3');
pp2(4).Color=hex2rgb('#e7298a');

ss2(1).MarkerEdgeColor=hex2rgb('#1b9E77');
ss2(2).MarkerEdgeColor=hex2rgb('#d95f02');
ss2(3).MarkerEdgeColor=hex2rgb('#7570b3');
% ss2(4).MarkerEdgeColor=hex2rgb('#e7298a');

ss2(1).MarkerFaceColor=hex2rgb('#1b9E77');
ss2(2).MarkerFaceColor=hex2rgb('#d95f02');
ss2(3).MarkerFaceColor=hex2rgb('#7570b3');
% ss2(4).MarkerFaceColor=hex2rgb('#e7298a');

ss(1).MarkerEdgeColor=hex2rgb('#1b9E77');
ss(2).MarkerEdgeColor=hex2rgb('#d95f02');
ss(3).MarkerEdgeColor=hex2rgb('#7570b3');
% ss(4).MarkerEdgeColor=hex2rgb('#e7298a');


set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',0:10:100,'Xminortick','off','Ytick',[0:1:15],'Yminortick','off','YTickLabel',{'','1','2','3','4','5','6','7','8','9','10','11','12','13','14','No test'},'YDir','reverse');
box off;
xlabel('Booster uptake','Fontsize',20,'Position',[50 17.7 -1]);
ylabel({'Minimum frequency','of testing  (days^{-1})'},'Fontsize',20,'Position',[-9.231752673157025,7.5,-0.1]);
legend boxoff;
legend(pp,{'Original','Alpha','Delta','Omicron'},'Fontsize',14,'Position',[0.607981220657277,0.338230330794659,0.107394364334538,0.102836876533918]);

xlim([0 100]);
ylim([-3.*dy 15])
xtickformat('percentage');
xtickangle(90);
text(-25.865925315696778,-0.181636363636366,'D','Fontsize',28,'FontWeight','bold');

% rmpath('Alternative_Test_Results');
print(gcf,['Figure2.png'],'-dpng','-r300');