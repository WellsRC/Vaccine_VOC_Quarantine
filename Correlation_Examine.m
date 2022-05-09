VOC_N={'NV','Alpha','Delta'};
VACB={'No','Yes'};
tsB=[6.3 5 4.3];

Data=readtable('ct_dat_refined.csv');
t_DATA=Data.TestDateIndex;
Ct_DATA=Data.CtT1;
VOC_DATA=Data.GreekLineage;
VB_DATA=Data.VaccineBreakthrough;

r_Kissler=zeros(3,2);
p_Kissler=zeros(3,2);
for VOC=0:2
    for VAC=0:1
        Org_Ct=Ct_DATA(strcmp(VOC_DATA,VOC_N{VOC+1}) & strcmp(VB_DATA,VACB{VAC+1}));
        Org_Time=t_DATA(strcmp(VOC_DATA,VOC_N{VOC+1}) & strcmp(VB_DATA,VACB{VAC+1}));
        ts=tsB(VOC+1);
        td=ts+20;
        load('Peak_Infection.mat','mmv','tsv');
        mm=pchip(tsv,mmv,ts);
        V_Org_NV = ViralShedding_Correlation(Org_Time+mm,td,ts,VAC,VOC);
        [r_Kissler(VOC+1,VAC+1),p_Kissler(VOC+1,VAC+1)]=corr(V_Org_NV,Org_Ct);
    end
end

VOC_N={'Delta';'Omicron';'Suspected Omicron'};
Data=readtable('ct_dat_refined-Omicron.csv');
t_DATA=Data.TestDateIndex;
Ct_DATA=Data.CtT1;
VOC_DATA=Data.GreekLineage;

r_Kissler_O=zeros(4,2);
p_Kissler_O=zeros(4,2);


tsB=[4.3 3.2 3.2 3.2];
for VOC_IND=1:4
    if(VOC_IND>=2)
        VOC=3;
    else
        VOV=2;
    end
    for VAC=0:1
        if(VOC_IND<4)
            Org_Ct=Ct_DATA(strcmp(VOC_DATA,VOC_N{VOC_IND}));
            Org_Time=t_DATA(strcmp(VOC_DATA,VOC_N{VOC_IND}));
        else            
            Org_Ct=Ct_DATA(strcmp(VOC_DATA,VOC_N{2}) | strcmp(VOC_DATA,VOC_N{3}));
            Org_Time=t_DATA(strcmp(VOC_DATA,VOC_N{2})| strcmp(VOC_DATA,VOC_N{3}));
        end
        ts=tsB(VOC_IND);
        td=ts+20;
        load('Peak_Infection.mat','mmv','tsv');
        mm=pchip(tsv,mmv,ts);
        V_Org_NV = ViralShedding_Correlation(Org_Time+mm,td,ts,VAC,VOC);
        [r_Kissler_O(VOC_IND,VAC+1),p_Kissler_O(VOC_IND,VAC+1)]=corr(V_Org_NV,Org_Ct);
    end
end

% 
% VOCN={'Pre-VOC','Delta','Omicron'};
% Data=readtable('41591_2022_1816_MOESM3_ESM.xlsx');
% Data=Data(strcmp(Data.IsolationSuccess,'yes'),:);
% VS=Data.VaccinationStatus;
% ND=Data.NumberOfDoses;
% DPOS=Data.DPOS;
% VL=Data.RNALoad_ml;
% FFUL=Data.FFU_ml;
% VARIANT=Data.Variant;
% 
% r_VL=zeros(3,2);
% p_VL=zeros(3,2);
% 
% 
% r_FFU=zeros(3,2);
% p_FFU=zeros(3,2);
% 
% tsB=[6.3 4.3 3];
% 
% for VOC=0:2
%     for VAC=0:1
%         figure('units','normalized','outerposition',[0. 0. 1 1]);
%         Org_VL=(VL(strcmp(VARIANT,VOCN{VOC+1}) & VS==VAC));
%         Org_FFU=(FFUL(strcmp(VARIANT,VOCN{VOC+1}) & VS==VAC));
%         Org_Time=DPOS(strcmp(VARIANT,VOCN{VOC+1}) & VS==VAC);
%         
%         scatter(Org_Time,log10(Org_VL))
%         
%         ts=tsB(VOC+1);
%         td=ts+20;
%         load('Peak_Infection.mat','mmv','tsv');
%         mm=pchip(tsv,mmv,ts);
%         if(VOC==0)
%             V_Org_NV = ViralShedding_Symptomatic(Org_Time+ts,td,ts,VAC,VOC);
%         else
%             V_Org_NV = ViralShedding_Symptomatic(Org_Time+ts,td,ts,VAC,VOC+1);
%         end
%         if(~isempty(V_Org_NV))
%             [r_VL(VOC+1,VAC+1),p_VL(VOC+1,VAC+1)]=corr(V_Org_NV,Org_VL);
%             [r_FFU(VOC+1,VAC+1),p_FFU(VOC+1,VAC+1)]=corr(V_Org_NV,Org_FFU);
%         end
%     end
% end
