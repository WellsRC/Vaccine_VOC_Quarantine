function S = TestSensitivity(t,ts,testtype,VAC,VOC)
%SensitivityvsViralLoad(V,asym) - Returns the sensitivity for a given viral
%load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ts- incubation period
% testtype - test type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S - Probability of that it is a true positive based on the viral load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Peak_Infection.mat','mmv','tsv');
% opts=optimset('TolX',10^(-16));
%  tsv=[0.1:0.1:28];
%  mmv=zeros(length(tsv),1);
%  for ii=1:280
%  mmv(ii)=fminbnd(@(x)-Relative_Infection_PCR(x,tsv(ii),0,0),0,tsv(ii),opts);
%  end
% save('Peak_Infection.mat');
mm=pchip(tsv,mmv,ts);


load('PCR_mapping.mat','VmI','SmI','VmS','SmS');
% load('MLE-Estimate-RTPCR.mat','beta')
% opts=optimset('TolX',10^(-16));
% ts=5.723;
% mm=fminbnd(@(x)-Relative_Infection_PCR(x,ts,0,0),0,ts,opts);
% t1=linspace(0,mm,1001);
% t2=linspace(mm,50,1001);
% [SmI]=PCRSens(t1,beta);
% [SmS]=PCRSens(t2,beta);
% VmI=Relative_Infection_PCR(t1,5.723,0,0);
% VmS=Relative_Infection_PCR(t2,5.723,0,0);
% save('PCR_mapping.mat');
 
 S(t<=mm)=pchip(VmI,SmI,Relative_Infection_PCR(t(t<=mm),ts,VAC,VOC));
 S(t>mm)=pchip(VmS,SmS,Relative_Infection_PCR(t(t>mm),ts,VAC,VOC));
if(~isempty(testtype))
    
    V=Relative_Infection_PCR(t,ts,VAC,VOC); % Use inf as we need to construct the mapping
    tt=[mm:0.1:90]; % made coarsered to improve the search tt=[ts:0.001:90]; 
    Vx=Relative_Infection_PCR(tt,ts,0,0); % Use inf as we need to construct the mapping
    PPA=LR(tt-ts,testtype);
    S=pchip(Vx,PPA,V).*S;
end

end

