clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT-PCR (Exit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TestingonExit_RTPCR_24hrDelay_Vaccinated=0.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.351 0.351 0.2204 0.275];
dist='Negative Binomial';
R_No_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_No_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Non_Vac= Probability_Onward(R_No_Vac,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Non_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_No_Vaccine');

load('TestingonExit_RTPCR_24hrDelay_Vaccinated=1.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.9611 0.9546 0.9064 0.275];
dist='Negative Binomial';
R_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Vac=Probability_Onward(R_Vac,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_Vaccine');


load('TestingonExit_RTPCR_24hrDelay_Vaccinated=1.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.9789 0.9754 0.9493 0.7629];
dist='Negative Binomial';
R_Third_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_Third_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Third_Vac=Probability_Onward(R_Third_Vac,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Third_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_Booster');



load('TestingonExit_RTPCR_24hrDelay_Vaccinated=0.mat','VOCv','qv','IDSLA','IDSLS');
pAv=[0.05 0.351 0.7];
R_pA=zeros(15,3);
for ii=1:3
    pA=pAv(ii);
    R_pA(:,ii)=pA.*IDSLA(VOCv==0)+(1-pA).*IDSLS(VOCv==0);
end

PQT_pA= Probability_Onward(R_pA,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_pA,'Variablenames',{['pA=' num2str(pAv(1)) ],['pA=' num2str(pAv(2)) ],['pA=' num2str(pAv(3)) ]})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_No_Vaccine_vary_pA');



R_RTPCR=zeros(15,2);
for VAC=0:1
    pA=0.351;
    load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=' num2str(VAC) '.mat'],'VOCv','qv','IDSLA','IDSLS');
    R_RTPCR(:,VAC+1)=pA.*IDSLA(VOCv==0)+(1-pA).*IDSLS(VOCv==0);
end

PQT_RTPCR= Probability_Onward(R_RTPCR,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_RTPCR,'Variablenames',{['Non-Vaccinated (pA=' num2str(pA) ],['Vaccinated (pA=' num2str(pA) ]})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_Vac_Status');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Original (Non-VOC)  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.351 0.9611 0.9789];

  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==0)+(1-pA(V)).*IDSLS(VOCv==0));
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
epsI=[0.92 0.98];
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
PQT_Immunity_Delta= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Delta,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_nonVOC_BI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Delta  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.2204 0.9064 0.9493];

  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==2)+(1-pA(V)).*IDSLS(VOCv==2));
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
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
PQT_Immunity_Delta= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Delta,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_Delta_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

R_Immunity=zeros(15,3);
for V=1:3
    if(V==1)
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
    R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
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
PQT_Immunity_Omicron= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Omicron,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_Omicron_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity three doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

R_Immunity=zeros(15,3);
for V=1:3
    if(V==1)
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_RTPCR_24hrDelay_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
    R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
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
PQT_Immunity_Omicron= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Omicron,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','RT-PCR_Exit_Omicron_BI_3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PanBio (Exit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TestingonExit_Abbot PanBio_Vaccinated=0.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.351 0.351 0.2204 0.275];
dist='Negative Binomial';
R_No_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_No_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Non_Vac= Probability_Onward(R_No_Vac,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Non_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_No_Vaccine');

load('TestingonExit_Abbot PanBio_Vaccinated=1.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.9611 0.9546 0.9064 0.275];
dist='Negative Binomial';
R_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Vac=Probability_Onward(R_Vac,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_Vaccine');


load('TestingonExit_Abbot PanBio_Vaccinated=1.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.9789 0.9754 0.9493 0.7629];
dist='Negative Binomial';
R_Third_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_Third_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Third_Vac=Probability_Onward(R_Third_Vac,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Third_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_Booster');



load('TestingonExit_Abbot PanBio_Vaccinated=0.mat','VOCv','qv','IDSLA','IDSLS');
pAv=[0.05 0.351 0.7];
R_pA=zeros(15,3);
for ii=1:3
    pA=pAv(ii);
    R_pA(:,ii)=pA.*IDSLA(VOCv==0)+(1-pA).*IDSLS(VOCv==0);
end

PQT_pA= Probability_Onward(R_pA,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_pA,'Variablenames',{['pA=' num2str(pAv(1)) ],['pA=' num2str(pAv(2)) ],['pA=' num2str(pAv(3)) ]})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_No_Vaccine_vary_pA');



R_RTPCR=zeros(15,2);
for VAC=0:1
    pA=0.351;
    load(['TestingonExit_Abbot PanBio_Vaccinated=' num2str(VAC) '.mat'],'VOCv','qv','IDSLA','IDSLS');
    R_RTPCR(:,VAC+1)=pA.*IDSLA(VOCv==0)+(1-pA).*IDSLS(VOCv==0);
end

PQT_RTPCR= Probability_Onward(R_RTPCR,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_RTPCR,'Variablenames',{['Non-Vaccinated (pA=' num2str(pA) ],['Vaccinated (pA=' num2str(pA) ]})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_Vac_Status');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Delta  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.2204 0.9064 0.9493];

  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==2)+(1-pA(V)).*IDSLS(VOCv==2));
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
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
PQT_Immunity_Delta= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Delta,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_Delta_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

% https://www.medpagetoday.com/special-reports/exclusives/96172
  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
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
PQT_Immunity_Omicron= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Omicron,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_Omicron_BI');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

% https://www.medpagetoday.com/special-reports/exclusives/96172
  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
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
PQT_Immunity_Omicron= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Omicron,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_Exit_Omicron_BI_3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PanBio (Entry_Exit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TestingonEntryExit_Abbot PanBio_Vaccinated=0.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.351 0.351 0.2204 0.275];
dist='Negative Binomial';
R_No_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_No_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Non_Vac= Probability_Onward(R_No_Vac,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Non_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_No_Vaccine');

load('TestingonEntryExit_Abbot PanBio_Vaccinated=1.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.9611 0.9546 0.9064 0.275];
dist='Negative Binomial';
R_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Vac=Probability_Onward(R_Vac,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_Vaccine');


load('TestingonEntryExit_Abbot PanBio_Vaccinated=1.mat','VOCv','qv','IDSLA','IDSLS');

pAv=[0.9789 0.9754 0.9493 0.7629];
dist='Negative Binomial';
R_Third_Vac=zeros(15,4);
for VOC=0:3
    pA=pAv(VOC+1);
    R_Third_Vac(:,VOC+1)=pA.*IDSLA(VOCv==VOC)+(1-pA).*IDSLS(VOCv==VOC);
end

PQT_Third_Vac=Probability_Onward(R_Third_Vac,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Third_Vac,'Variablenames',{'Non-VOC','Alpha VOC','Delta VOC','Omicron VOC'})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_Booster');



load('TestingonEntryExit_Abbot PanBio_Vaccinated=0.mat','VOCv','qv','IDSLA','IDSLS');
pAv=[0.05 0.351 0.7];
R_pA=zeros(15,3);
for ii=1:3
    pA=pAv(ii);
    R_pA(:,ii)=pA.*IDSLA(VOCv==0)+(1-pA).*IDSLS(VOCv==0);
end

PQT_pA= Probability_Onward(R_pA,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_pA,'Variablenames',{['pA=' num2str(pAv(1)) ],['pA=' num2str(pAv(2)) ],['pA=' num2str(pAv(3)) ]})];


writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_No_Vac_vary_pA');



R_RTPCR=zeros(15,2);
for VAC=0:1
    pA=0.351;
    load(['TestingonEntryExit_Abbot PanBio_Vaccinated=' num2str(VAC) '.mat'],'VOCv','qv','IDSLA','IDSLS');
    R_RTPCR(:,VAC+1)=pA.*IDSLA(VOCv==0)+(1-pA).*IDSLS(VOCv==0);
end

PQT_RTPCR= Probability_Onward(R_RTPCR,dist);


T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_RTPCR,'Variablenames',{['Non-Vaccinated (pA=' num2str(pA) ],['Vaccinated (pA=' num2str(pA) ]})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_Vac_Status');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Delta  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.2204 0.9064 0.9493];

  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonEntryExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonEntryExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==2)+(1-pA(V)).*IDSLS(VOCv==2));
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
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
PQT_Immunity_Delta= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Delta,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_Delta_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonEntryExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonEntryExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
    
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
    
end
v=[1-(88.48/211.36) (88.48/211.36)]; % Break down of 2-dose and booster

RIm=zeros(15,11);
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
PQT_Immunity_Omicron= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Omicron,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_Omicron_BI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Omicron  background immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
pA=[0.275 0.275 0.7629];

  R_Immunity=zeros(15,10);
for V=1:3
    if(V==1)
        load(['TestingonEntryExit_Abbot PanBio_Vaccinated=0.mat'],'VOCv','qv','IDSLA','IDSLS');
    else
        load(['TestingonEntryExit_Abbot PanBio_Vaccinated=1.mat'],'VOCv','qv','IDSLA','IDSLS');
    end
    
        R_Immunity(:,V)=(pA(V).*IDSLA(VOCv==3)+(1-pA(V)).*IDSLS(VOCv==3));
    
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
PQT_Immunity_Omicron= Probability_Onward(RIm,dist);

T_PQT_q=table(qv(VOCv==0),'Variablenames',{'Quarantine duration'});
T_PQT=[T_PQT_q array2table(PQT_Immunity_Omicron,'Variablenames',{'Vac_Up=0%','Vac_Up=10%','Vac_Up=20%','Vac_Up=30%','Vac_Up=40%','Vac_Up=50%','Vac_Up=60%','Vac_Up=70%','Vac_Up=80%','Vac_Up=90%','Vac_Up=100%'})];

writetable(T_PQT,'Post-Quarantine-Transmission.xlsx','Sheet','PanBio_EntryExit_Omicron_BI_3');

