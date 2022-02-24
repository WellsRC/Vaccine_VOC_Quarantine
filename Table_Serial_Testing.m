clear;
clc;
R0=[2.79 4.19 5.08 6.93];
pAv=[0.351 0.351 0.2204 0.275];
R_No_Vac=zeros(14,4);
for VOC=0:3
    
    load(['VOC=' num2str(VOC) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=0.mat'],'RTotA','RTotS');
    pA=pAv(VOC+1);
    R_No_Vac(:,VOC+1)=pA.*RTotA+(1-pA).*RTotS;
end
Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(R_No_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];
writetable(T_ST,'Serial_Testing.xlsx','Sheet','RT-PCR_24h_Non_Vac');

pAv=[0.9611 0.9546 0.9064 0.275];
R_Vac=zeros(14,4);
for VOC=0:3
    load(['VOC=' num2str(VOC) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=1.mat'],'RTotA','RTotS');
    pA=pAv(VOC+1);
    R_Vac(:,VOC+1)=pA.*RTotA+(1-pA).*RTotS;
end
Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(R_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];
writetable(T_ST,'Serial_Testing.xlsx','Sheet','RT-PCR_24h_Vac');


pAv=[0.9789 0.9754 0.9493 0.7629];
R_Third_Vac=zeros(14,4);
for VOC=0:3
    load(['VOC=' num2str(VOC) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=1.mat'],'RTotA','RTotS');
    pA=pAv(VOC+1);
    R_Third_Vac(:,VOC+1)=pA.*RTotA+(1-pA).*RTotS;
end

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(R_Third_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];
writetable(T_ST,'Serial_Testing.xlsx','Sheet','RT-PCR_24h_Third_Vac');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Delta  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.2204 0.9064 0.9493];

R_Immunity=zeros(14,3);
for V=1:3
    if(V==1)
        load(['VOC=' num2str(2) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=0.mat'],'RTotA','RTotS');
    else
        load(['VOC=' num2str(2) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=1.mat'],'RTotA','RTotS');
    end
    R_Immunity(:,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(14,11);
epsI=[0.76 0.94];
for jj=1:11
    vacu=(jj-1)./10;
    V_Im_Pop=vacu.*(v(1).*epsI(1)+v(2).*epsI(2));
    
    w1=1-vacu;
    w2=vacu.*v(1).*(1-epsI(1));
    w3=vacu.*v(2).*(1-epsI(2));
    w=w1+w2+w3;
    w1=w1./w;
    w2=w2./w;
    w3=w3./w;
    RIm(:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
end

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(RIm,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_ST,'Serial_Testing.xlsx','Sheet','RT-PCR_Delta_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

R_Immunity=zeros(14,3);
for V=1:3
    if(V==1)
        load(['VOC=' num2str(3) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=0.mat'],'RTotA','RTotS');
    else
        load(['VOC=' num2str(3) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=1.mat'],'RTotA','RTotS');
    end
    R_Immunity(:,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(14,11);
epsI=[0.38 0.82];
for jj=1:11
    vacu=(jj-1)./10;
    V_Im_Pop=vacu.*(v(1).*epsI(1)+v(2).*epsI(2));
    
    w1=1-vacu;
    w2=vacu.*v(1).*(1-epsI(1));
    w3=vacu.*v(2).*(1-epsI(2));
    w=w1+w2+w3;
    w1=w1./w;
    w2=w2./w;
    w3=w3./w;
    RIm(:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
end

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(RIm,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_ST,'Serial_Testing.xlsx','Sheet','RT-PCR_Omicron_BI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

R_Immunity=zeros(14,3);
for V=1:3
    if(V==1)
        load(['VOC=' num2str(3) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=0.mat'],'RTotA','RTotS');
    else
        load(['VOC=' num2str(3) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=1.mat'],'RTotA','RTotS');
    end
    R_Immunity(:,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
end

RIm=zeros(14,11);
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

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(RIm,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_ST,'Serial_Testing.xlsx','Sheet','RT-PCR_Omicron_BI3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% PanBio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0=[2.79 4.19 5.08 6.93];
pAv=[0.351 0.351 0.2204 0.275];
R_No_Vac=zeros(14,4);
for VOC=0:3
    
    load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(VOC) '_VAC= ' num2str(0) '.mat'],'RTotA','RTotS');
    pA=pAv(VOC+1);
    R_No_Vac(:,VOC+1)=pA.*RTotA+(1-pA).*RTotS;
end
Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(R_No_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];
writetable(T_ST,'Serial_Testing.xlsx','Sheet','PanBio_Non_Vac');

pAv=[0.9611 0.9546 0.9064 0.275];
R_Vac=zeros(14,4);
for VOC=0:3
    load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(VOC) '_VAC= ' num2str(1) '.mat'],'RTotA','RTotS');
    pA=pAv(VOC+1);
    R_Vac(:,VOC+1)=pA.*RTotA+(1-pA).*RTotS;
end
Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(R_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];
writetable(T_ST,'Serial_Testing.xlsx','Sheet','PanBio_Vac');


pAv=[0.9789 0.9754 0.9493 0.7629];
R_Third_Vac=zeros(14,4);
for VOC=0:3
    load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(VOC) '_VAC= ' num2str(1) '.mat'],'RTotA','RTotS');
    pA=pAv(VOC+1);
    R_Third_Vac(:,VOC+1)=pA.*RTotA+(1-pA).*RTotS;
end

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(R_Third_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];
writetable(T_ST,'Serial_Testing.xlsx','Sheet','PanBio_Third_Vac');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Delta  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.2204 0.9064 0.9493];

R_Immunity=zeros(14,3);
for V=1:3
    if(V==1)
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(2) '_VAC= ' num2str(0) '.mat'],'RTotA','RTotS');
    else
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(2) '_VAC= ' num2str(1) '.mat'],'RTotA','RTotS');
    end
    R_Immunity(:,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(14,11);
epsI=[0.76 0.94];
for jj=1:11
    vacu=(jj-1)./10;
    V_Im_Pop=vacu.*(v(1).*epsI(1)+v(2).*epsI(2));
    
    w1=1-vacu;
    w2=vacu.*v(1).*(1-epsI(1));
    w3=vacu.*v(2).*(1-epsI(2));
    w=w1+w2+w3;
    w1=w1./w;
    w2=w2./w;
    w3=w3./w;
    RIm(:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
end

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(RIm,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_ST,'Serial_Testing.xlsx','Sheet','PanBio_Delta_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

R_Immunity=zeros(14,3);
for V=1:3
    if(V==1)
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(3) '_VAC= ' num2str(0) '.mat'],'RTotA','RTotS');
    else
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(3) '_VAC= ' num2str(1) '.mat'],'RTotA','RTotS');
    end
    R_Immunity(:,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(14,11);
epsI=[0.38 0.82];
for jj=1:11
    vacu=(jj-1)./10;
    V_Im_Pop=vacu.*(v(1).*epsI(1)+v(2).*epsI(2));
    
    w1=1-vacu;
    w2=vacu.*v(1).*(1-epsI(1));
    w3=vacu.*v(2).*(1-epsI(2));
    w=w1+w2+w3;
    w1=w1./w;
    w2=w2./w;
    w3=w3./w;
    RIm(:,jj)=[w1.*R_Immunity(:,1)+w2.*R_Immunity(:,2)+w3.*R_Immunity(:,3)].*(1-V_Im_Pop);
end

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(RIm,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_ST,'Serial_Testing.xlsx','Sheet','PanBio_Omicron_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

R_Immunity=zeros(14,3);
for V=1:3
    if(V==1)
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(3) '_VAC= ' num2str(0) '.mat'],'RTotA','RTotS');
    else
        load(['Testing_Frequency_Abbot PanBio_VOC= ' num2str(3) '_VAC= ' num2str(1) '.mat'],'RTotA','RTotS');
    end
    R_Immunity(:,V)=pA(V).*RTotA+(1-pA(V)).*RTotS;
end

RIm=zeros(14,11);
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

Freq=[1:14]';
T_ST_q=table(Freq,'Variablenames',{'Freq.'});
T_ST=[T_ST_q array2table(RIm,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_ST,'Serial_Testing.xlsx','Sheet','PanBio_Omicron_BI3');

